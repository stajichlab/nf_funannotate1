#!/usr/bin/env bash
# record_versions.sh — write a per-task tool-version record (YAML).
#
# Usage:
#   record_versions.sh DEST_DIR PROCESS SAMPLE CONTAINER PIPELINE_REV NAME=CMD [NAME=CMD ...]
#
# Runs each CMD (via bash -c, 60 s limit), keeps the first non-empty output
# line, and writes DEST_DIR/<PROCESS>.versions.yml. A failing CMD is recorded
# as "unavailable"; this script never fails the calling task.
# Env: RECORD_RETROACTIVE=1 marks the record as written after the fact (adds
# retroactive: "true"); RECORD_NOTE adds a free-text note: field.
# Written to be called from Nextflow process scripts, inside the task container.

dest="$1"; proc="$2"; sample="$3"; container="$4"; rev="$5"; shift 5 || exit 0

yq() { printf '%s' "$1" | tr -d '\r' | sed -e 's/\\/\\\\/g' -e 's/"/\\"/g'; }

mkdir -p "$dest" 2>/dev/null || exit 0
step="${proc##*:}"
out="$dest/${step}.versions.yml"
tmp="$out.tmp.$$"

{
    printf '"%s":\n' "$(yq "$proc")"
    printf '  sample: "%s"\n'    "$(yq "$sample")"
    printf '  container: "%s"\n' "$(yq "$container")"
    printf '  pipeline: "%s"\n'  "$(yq "$rev")"
    printf '  host: "%s"\n'      "$(yq "$(hostname 2>/dev/null)")"
    printf '  recorded: "%s"\n'  "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  retroactive: "%s"\n' "$([ -n "${RECORD_RETROACTIVE:-}" ] && echo true || echo false)"
    [ -n "${RECORD_NOTE:-}" ] && printf '  note: "%s"\n' "$(yq "$RECORD_NOTE")"
    printf '  tools:\n'
    for spec in "$@"; do
        name="${spec%%=*}"; cmd="${spec#*=}"
        raw=$(timeout 60 bash -c "$cmd" 2>&1); rc=$?
        v=$(printf '%s\n' "$raw" | tr -d '\r' | grep -v -e '^WARNING' -e 'setlocale' -e '^bash: warning' | sed -n '/[^[:space:]]/{p;q}')
        [ "$rc" -eq 127 ] && v=""
        [ -z "$v" ] && v="unavailable"
        printf '    %s: "%s"\n' "$name" "$(yq "$v")"
    done
} > "$tmp" 2>/dev/null && mv -f "$tmp" "$out"
exit 0
