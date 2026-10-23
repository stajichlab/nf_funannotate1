#!/usr/bin/env bash
# Regression test: every image that defaults to <sif_dir>/<name>.sif must follow a sif_dir set
# AFTER nextflow.config is read -- on the command line, in a params file, by a site profile
# (ucr_hpcc), by an environment variable, or in a `-c` site config file (read last). Before the fix
# these paths were strings built when nextflow.config / provision_singularity.config were parsed, so
# a later sif_dir was ignored (the funannotate image stayed in <launchDir>/singularity_cache under
# `-profile ...,ucr_hpcc,singularity`; found by the UCR HPCC end-to-end test).
#
# Runs tests/nf/container_resolution_test.nf (one process per label) under the real configs with the
# apptainer engine switched off, and reads each label's resolved container from the task command.
# Usage: bash tests/test_sif_dir_resolution.sh     (needs nextflow; ~1-2 min)
set -uo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HARNESS="$PROJECT_DIR/tests/nf/container_resolution_test.nf"
command -v nextflow >/dev/null || { echo "SKIP: nextflow not found"; exit 0; }
T="$(cd "$(mktemp -d)" && pwd -P)"; trap 'rm -rf "$T"' EXIT
pass=0; fail=0
check() { if eval "$2"; then pass=$((pass+1)); echo "ok   $1"; else fail=$((fail+1)); echo "FAIL $1"; fi; }
printf "apptainer.enabled = false\n" > "$T/off.config"
printf "apptainer.enabled = false\nparams.sif_dir = '/site/cfg'\n" > "$T/site.config"
printf "sif_dir: /site/pfile\n" > "$T/p.yaml"

n=0
resolve() {  # <name> <env assignments...> -- <nextflow run args...>: writes $T/<name>.txt (label=container)
    local name="$1"; shift; local envs=()
    while [ "$1" != "--" ]; do envs+=("$1"); shift; done; shift
    n=$((n+1)); mkdir -p "$T/run$n"
    ( cd "$T/run$n" && env -u NXF_APPTAINER_CACHE -u NXF_SINGULARITY_CACHE -u APPTAINER_CACHE ${envs[@]+"${envs[@]}"} \
        nextflow -c "$PROJECT_DIR/nextflow.config" run "$HARNESS" -ansi-log false -w "$T/run$n/w" "$@" 2>&1 ) \
        | grep -oE '(funannotate|genome_clean|antismash|signalp|deeptmhmm)=[^ ]+' | sort -u > "$T/$name.txt"
}
all_under() {  # all five labels resolve under <dir>/
    local f="$1" d="$2"; [ "$(cut -d= -f1 "$f" | sort -u | wc -l | tr -d ' ')" -eq 5 ] && ! grep -v "=$d/" "$f" | grep -q .
}
P="-profile annotate,local,singularity"

resolve default -- $P -c "$T/off.config"
check "default: all five images under <launchDir>/singularity_cache" 'grep -q "singularity_cache/" "$T/default.txt" && all_under "$T/default.txt" "$(sed -n "s/^funannotate=\(.*\)\/funannotate-1.9.0.sif$/\1/p" "$T/default.txt")"'

resolve cfg -- $P -c "$T/site.config"
check "sif_dir in a -c site config: all five images follow it"   'all_under "$T/cfg.txt" /site/cfg'

resolve cli -- $P -c "$T/off.config" --sif_dir /site/cli
check "--sif_dir: all five images follow it"                      'all_under "$T/cli.txt" /site/cli'

resolve pfile -- $P -c "$T/off.config" -params-file "$T/p.yaml"
check "sif_dir in a params file: all five images follow it"       'all_under "$T/pfile.txt" /site/pfile'

resolve prof -- -profile annotate,local,ucr_hpcc,singularity -c "$T/off.config"
check "site profile (ucr_hpcc) sets sif_dir: all five follow it"  'all_under "$T/prof.txt" /bigdata/stajichlab/shared/lib/singularity_cache'

resolve envv NXF_APPTAINER_CACHE=/site/env -- $P -c "$T/off.config"
check "NXF_APPTAINER_CACHE: all five images follow it"            'all_under "$T/envv.txt" /site/env'

resolve expl -- $P -c "$T/site.config" --container_funannotate /x/my.sif
check "an explicit --container_funannotate wins"                  'grep -q "^funannotate=/x/my.sif$" "$T/expl.txt" && grep -q "^antismash=/site/cfg/" "$T/expl.txt"'

echo "passed=$pass failed=$fail"
[ "$fail" -eq 0 ] || { for f in default cfg cli pfile prof envv expl; do echo "--- $f"; cat "$T/$f.txt"; done; exit 1; }
