#!/usr/bin/env bash
# Regression test: the funannotate image and the apptainer cache must follow a sif_dir set by a site
# profile, like the other <sif_dir>/<name>.sif images do. Before the fix container_funannotate was a
# "${params.sif_dir}/..." string built when nextflow.config was parsed, so under
# `-profile annotate,slurm,ucr_hpcc,singularity` it stayed at <launchDir>/singularity_cache while
# container_antismash etc. went to the site's shared cache (found by the HPCC end-to-end test: the
# funannotate image was "not found" in the launch directory). Reads the resolved config only.
# Usage: bash tests/test_sif_dir_resolution.sh   (needs nextflow)
set -uo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
command -v nextflow >/dev/null || { echo "SKIP: nextflow not found"; exit 0; }
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
pass=0; fail=0
check() { if eval "$2"; then pass=$((pass+1)); echo "ok   $1"; else fail=$((fail+1)); echo "FAIL $1"; fi; }
cfg() { ( cd "$T" && env -u NXF_APPTAINER_CACHE -u NXF_SINGULARITY_CACHE -u APPTAINER_CACHE "$@" \
        nextflow ${NF_C:+-c "$NF_C"} config "$PROJECT_DIR" -flat ${NF_PROFILE:+-profile "$NF_PROFILE"} 2>/dev/null ); }
val() { grep -E "^(params\.)?$1 = " | head -1 | sed "s/^[^=]*= *//; s/'//g"; }   # one resolved key

# 1. default: everything under <launchDir>/singularity_cache
out=$(NF_PROFILE=annotate,slurm,singularity cfg env)
check "default: funannotate image under singularity_cache" '[[ "$(echo "$out" | val container_funannotate)" == */singularity_cache/funannotate-1.9.0.sif ]]'
check "default: apptainer.cacheDir is sif_dir"             '[ "$(echo "$out" | val apptainer.cacheDir)" = "$(echo "$out" | val sif_dir)" ]'

# 2. a site profile sets sif_dir (ucr_hpcc): funannotate image, antismash image and cacheDir all follow
out=$(NF_PROFILE=annotate,slurm,ucr_hpcc,singularity cfg env)
D=$(echo "$out" | val sif_dir)
check "ucr_hpcc profile: sif_dir is the shared cache ($D)"        '[ -n "$D" ] && [[ "$D" != */nf_funannotate1/singularity_cache && "$D" != "$T"* ]]'
check "ucr_hpcc profile: funannotate image follows sif_dir"       '[ "$(echo "$out" | val container_funannotate)" = "$D/funannotate-1.9.0.sif" ]'
check "ucr_hpcc profile: antismash image follows sif_dir"          '[ "$(echo "$out" | val container_antismash)" = "$D/antismash-standalone-8.0.4-procps.sif" ]'
check "ucr_hpcc profile: apptainer.cacheDir follows sif_dir"       '[ "$(echo "$out" | val apptainer.cacheDir)" = "$D" ]'

# 3. NXF_APPTAINER_CACHE wins everywhere
out=$(NF_PROFILE=annotate,slurm,singularity cfg env NXF_APPTAINER_CACHE=/env/cache)
check "env var: funannotate image and cacheDir under it" '[ "$(echo "$out" | val container_funannotate)" = "/env/cache/funannotate-1.9.0.sif" ] && [ "$(echo "$out" | val apptainer.cacheDir)" = "/env/cache" ]'

# 4. an explicit container_funannotate (here from a -c file) is not overridden
printf "params.container_funannotate = '/x/my.sif'\n" > "$T/c.config"
out=$(NF_C="$T/c.config" NF_PROFILE=annotate,slurm,singularity cfg env)
check "explicit container_funannotate wins" '[ "$(echo "$out" | val container_funannotate)" = "/x/my.sif" ]'

# 5. other provisioning axes keep the unset sentinel; mariadbImage resolves it when a task runs
out=$(NF_PROFILE=annotate,slurm,conda cfg env)
check "conda axis: container_funannotate unset (resolved lazily by mariadbImage)" '[ -z "$(echo "$out" | val container_funannotate)" ]'

echo "passed=$pass failed=$fail"
[ "$fail" -eq 0 ] || exit 1
