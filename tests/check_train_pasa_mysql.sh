#!/usr/bin/env bash
# Check a launch directory made by tests/submit_train_pasa_mysql.sh.
# Reports, for each FUNANNOTATE_TRAIN task attempt: status, which MariaDB branch
# ran (in-image or sidecar), MariaDB readiness, the PASA backend, and whether
# training finished.
#
# Usage: bash tests/check_train_pasa_mysql.sh <launch dir> [<launch dir> ...]
# Exit status: 0 if every launch dir has a FUNANNOTATE_TRAIN task with exit 0,
# pasa_mysql true, a MariaDB branch reached, PASA finished, and no MariaDB
# errors; 1 otherwise.

set -uo pipefail

[ $# -ge 1 ] || { echo "usage: $0 <launch dir> [...]" >&2; exit 2; }

overall=0
for LAUNCH in "$@"; do
    echo "=== ${LAUNCH}"
    [ -f "${LAUNCH}/run_info.txt" ] && sed 's/^/    /' "${LAUNCH}/run_info.txt"

    # Newest trace file written by the pipeline for this launch dir.
    trace=$(ls -t "${LAUNCH}"/logs/nextflow/*trace*.txt "${LAUNCH}"/*trace*.txt 2>/dev/null | head -1)
    workdirs=()
    if [ -n "$trace" ]; then
        echo "    trace: $trace"
        # columns vary by config: locate them from the header.
        mapfile -t rows < <(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
            $h["name"] ~ /FUNANNOTATE_TRAIN/ {print $h["hash"] "\t" $h["status"] "\t" $h["name"]}' "$trace")
        for r in "${rows[@]}"; do
            hash=${r%%$'\t'*}
            # work dir is <launch>/work/annotate/xx/yyyy... (workDir set by the annotate profile)
            d=$(find "${LAUNCH}/work" -mindepth 2 -maxdepth 3 -type d -path "*/${hash}*" 2>/dev/null | head -1)
            [ -n "$d" ] && workdirs+=("$d|${r#*$'\t'}")
        done
    fi
    if [ ${#workdirs[@]} -eq 0 ]; then
        # No trace yet (run still starting) -> search work/ directly.
        while IFS= read -r f; do
            workdirs+=("$(dirname "$f")|unknown")
        done < <(find "${LAUNCH}/work" -maxdepth 4 -name .command.sh -exec grep -l 'Running funannotate train\|pasa_db_arg=' {} + 2>/dev/null)
    fi
    if [ ${#workdirs[@]} -eq 0 ]; then
        echo "    no FUNANNOTATE_TRAIN task found yet"
        overall=1
        continue
    fi

    ok_here=0
    for entry in "${workdirs[@]}"; do
        d=${entry%%|*}; info=${entry#*|}
        logs=("$d/.command.out" "$d/.command.err" "$d/.command.log")
        # first match across the logs (.command.out and .command.log repeat lines)
        g() { grep -h -m1 -E "$1" "${logs[@]}" 2>/dev/null | head -1; }
        echo "  - task: ${info//$'\t'/  }"
        echo "    workdir: $d"
        exitcode=$(cat "$d/.exitcode" 2>/dev/null || echo "running/none")
        echo "    exit code: ${exitcode}"
        if g 'Using in-image .* PASA MariaDB backend' >/dev/null; then branch="in-image"
        elif g "Initializing fresh MariaDB system tables via sidecar" >/dev/null; then branch="sidecar"
        else branch="not reached"; fi
        echo "    MariaDB branch: ${branch}"
        ready=$(g 'did not become ready' || true)
        [ -n "$ready" ] && echo "    MariaDB: NOT READY -> ${ready}"
        errs=$(grep -h -E 'ERROR: .*(mariadb|mysql|MariaDB|sidecar)|Can.t connect to (MySQL|server)|Read-only file system' "${logs[@]}" 2>/dev/null | sort -u | head -5)
        [ -n "$errs" ] && { echo "    MariaDB errors:"; echo "$errs" | sed 's/^/      /'; }
        # module prints "mysql is <params.pasa_mysql>"; funannotate prints
        # "PASA finished. PASAweb accessible via: ...?db=<db>" (a db name for
        # mysql, a file path for sqlite).
        mysql_flag=$(g '^mysql is ' | awk '{print $3}')
        echo "    pasa_mysql: ${mysql_flag:-not seen yet}"
        g 'Running funannotate train' | sed 's/^/    /'
        g 'PASA assigned' | sed 's/^/    /'
        pasa_done=$(g 'PASA finished' || true)
        [ -n "$pasa_done" ] && echo "    PASA db: ${pasa_done##*db=}"
        g 'Training cleanup complete' | sed 's/^/    /'
        if [ "$exitcode" = "0" ] && [ "$mysql_flag" = "true" ] && [ -n "$pasa_done" ] \
           && [ -z "$errs" ] && [ -z "$ready" ] && [ "$branch" != "not reached" ]; then
            ok_here=1
        fi
    done
    if [ $ok_here -eq 1 ]; then echo "  RESULT: PASS"; else echo "  RESULT: not passed (yet)"; overall=1; fi
done
exit $overall
