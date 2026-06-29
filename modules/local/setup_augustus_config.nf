// Seed a writable AUGUSTUS_CONFIG copy at params.augustus_config from the installed
// augustus config. Augustus (via funannotate train/predict) writes new species parameter
// sets into its config dir, so it cannot use the read-only config that ships with a
// module/conda/singularity install — every run needs its own writable copy. Runs under the
// 'funannotate' label so the source is the augustus that the active provisioning profile
// supplies; the install's config is located via AUGUSTUS_CONFIG_PATH (set by the module/
// conda env) or by resolving ../config from the augustus binary, or an explicit override
// (params.augustus_config_source). storeDir caches the populated dir at params.augustus_config,
// so this runs at most once across all pipeline runs; if the directory already exists the
// task is skipped entirely but still emits `config`.
process SETUP_AUGUSTUS_CONFIG {
    label 'funannotate'
    label 'process_single'

    // Closure defers evaluation to task runtime so a missing pipeline profile is
    // caught by the workflow guard (clear message) instead of throwing here.
    storeDir { file(params.augustus_config).parent }

    output:
    path "${cfg_dir}", emit: config

    script:
    cfg_dir = file(params.augustus_config).name
    def override = params.augustus_config_source ? params.augustus_config_source : ''
    """
    set -euo pipefail

    SRC="${override}"
    if [ -z "\$SRC" ] && [ -n "\${AUGUSTUS_CONFIG_PATH:-}" ] && [ -d "\${AUGUSTUS_CONFIG_PATH}" ]; then
        SRC="\${AUGUSTUS_CONFIG_PATH}"
    fi
    if [ -z "\$SRC" ] && command -v augustus >/dev/null 2>&1; then
        cand="\$(dirname "\$(command -v augustus)")/../config"
        [ -d "\$cand" ] && SRC="\$(readlink -f "\$cand")"
    fi
    if [ -z "\$SRC" ] || [ ! -d "\$SRC" ]; then
        echo "[ERROR] Could not locate an installed augustus config to copy." >&2
        echo "        Set AUGUSTUS_CONFIG_PATH in the provisioning environment, ensure 'augustus' is on PATH," >&2
        echo "        or pass --augustus_config_source /path/to/augustus/config." >&2
        exit 1
    fi

    echo "[INFO] Seeding writable augustus config at ${cfg_dir} from \$SRC"
    mkdir -p ${cfg_dir}
    cp -a "\$SRC/." ${cfg_dir}/
    echo "[INFO] augustus config ready at ${cfg_dir}"
    """

    stub:
    cfg_dir = file(params.augustus_config).name
    """
    mkdir -p ${cfg_dir}/species
    """
}
