#!/usr/bin/env bash
# Run Nextflow with its cache (.nextflow: the LevelDB task cache + run history)
# on pod-local disk, snapshotting it to the PVC atomically.
#
# Why: the cache did not survive losing the node mid-write. On 2026-09-28 a
# node dropped the Nextflow pod, the cache on the PVC was left corrupt
# ("Corruption: no meta-nextfile entry in descriptor") and every retry failed
# until the DB was repaired by hand. Part of the fix is below (leveldb.mmap).
#
# How: NXF_CACHE_DIR points at an emptyDir. Every CACHE_SYNC_SECONDS the
# cache is copied (retrying until no file changed during the copy, i.e. a
# point-in-time copy) into $RUN_DIR/.nextflow-snapshots/.tmp-<ts>, fsynced,
# then renamed to snap-<ts>. The rename is atomic, so a snap-* directory is
# always complete; a crash leaves at most a .tmp-* that is ignored and cleaned
# up. The newest CACHE_SNAPSHOTS_KEEP snapshots are kept; a final one is taken
# after Nextflow exits. On start the newest snapshot is restored (the first
# time, an existing $RUN_DIR/.nextflow is migrated), then `-resume` as before.
#
# If the newest snapshot will not open, delete it and relaunch to fall back to
# the previous one.
set -uo pipefail

LOCAL=/nxf-local
SNAPS="$RUN_DIR/.nextflow-snapshots"
INTERVAL="${CACHE_SYNC_SECONDS:-300}"
KEEP="${CACHE_SNAPSHOTS_KEEP:-3}"

log() { echo "[nf-run $(date -u +%H:%M:%S)] $*"; }

cd "$RUN_DIR" || exit 1
mkdir -p "$SNAPS"
rm -rf "$SNAPS"/.tmp-*

latest=$(ls -1d "$SNAPS"/snap-* 2>/dev/null | sort | tail -1)
if [ -n "$latest" ]; then
    log "restoring cache from $latest"
    cp -a "$latest/.nextflow" "$LOCAL/"
elif [ -d "$RUN_DIR/.nextflow" ]; then
    # Run started before the cache moved off the PVC: migrate it once, and
    # rename the PVC copy so nothing resumes from the stale one later.
    log "migrating $RUN_DIR/.nextflow to local disk"
    cp -a "$RUN_DIR/.nextflow" "$LOCAL/"
    mv "$RUN_DIR/.nextflow" "$RUN_DIR/.nextflow.migrated-$(date -u +%Y%m%dT%H%M%S)"
fi
mkdir -p "$LOCAL/.nextflow"
export NXF_CACHE_DIR="$LOCAL/.nextflow"
# Nextflow's LevelDB (iq80) memory-maps its log and MANIFEST into 1 MB
# zero-filled files by default. Mapped writes are invisible to a file copy and
# are lost outright if the node dies (the 2026-09-28 corruption: all-zero
# MANIFEST and log). Plain file writes are copyable and recover after a crash.
export NXF_OPTS="${NXF_OPTS:-} -Dleveldb.mmap=false"

state() { ls -lAR --time-style=full-iso "$LOCAL/.nextflow" 2>/dev/null | md5sum; }

snapshot() {
    local ts tmp before i
    ts=$(date -u +%Y%m%dT%H%M%S)
    tmp="$SNAPS/.tmp-$ts"
    for i in 1 2 3 4 5 6 7 8 9 10; do
        rm -rf "$LOCAL/.snap"
        before=$(state)
        cp -a "$LOCAL/.nextflow" "$LOCAL/.snap" || return 1
        [ "$before" = "$(state)" ] && break
        [ "$i" = 10 ] && { log "cache kept changing; snapshot skipped"; return 1; }
        sleep 1
    done
    mkdir -p "$tmp"
    cp -a "$LOCAL/.snap" "$tmp/.nextflow" && sync -f "$tmp" || { rm -rf "$tmp"; return 1; }
    mv "$tmp" "$SNAPS/snap-$ts" || return 1
    ls -1d "$SNAPS"/snap-* | sort | head -n -"$KEEP" | while read -r old; do rm -rf "$old"; done
}

[ -f nf_run.log ] && mv nf_run.log "nf_run.$(date -u +%Y%m%dT%H%M%S).log"
# shellcheck disable=SC2086  # NF_EXTRA_ARGS is a word list
nextflow run "$PIPELINE" -profile "$PROFILE" -c "$SITE_CONFIG" \
    -params-file "$PARAMS_FILE" -resume $NF_EXTRA_ARGS > >(tee nf_run.log) 2>&1 &
NF_PID=$!

# `kubectl delete` sends SIGTERM here: pass it on so Nextflow cleans up its
# task pods, then fall through to the final snapshot.
trap 'log "SIGTERM: stopping Nextflow"; kill -TERM "$NF_PID" 2>/dev/null' TERM

while kill -0 "$NF_PID" 2>/dev/null; do
    sleep "$INTERVAL" & wait $! 2>/dev/null
    kill -0 "$NF_PID" 2>/dev/null && snapshot
done
wait "$NF_PID"; rc=$?
# The trap interrupts the first wait; wait again for Nextflow's real exit status.
while kill -0 "$NF_PID" 2>/dev/null; do wait "$NF_PID"; rc=$?; done
snapshot && log "final cache snapshot saved" || log "WARNING: final cache snapshot failed"
exit "$rc"
