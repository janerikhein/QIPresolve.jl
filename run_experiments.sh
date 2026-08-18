#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  ./run_experiments.sh INSTANCE_DIR STATS_DIR [presolve_lp_scip_stats options...]

Submits one SLURM job per top-level *.lp instance. Stats JSON files are written
to STATS_DIR as <instance>_stats.json.

Do not pass --output, --out-dir, or --time-limit. This wrapper sets --output
per instance, and SCIP time limits should come from --scip-config.

Cluster defaults can be overridden with:
  QIP_SBATCH_PARTITION       default: opt_int
  QIP_SBATCH_CONSTRAINT      default: Gold5122
  QIP_SBATCH_TIME            default: 02:05:00
  QIP_SBATCH_EXCLUSIVE       default: 1
  QIP_SBATCH_CPUS_PER_TASK   optional
USAGE
}

if [[ $# -lt 2 ]]; then
    usage >&2
    exit 2
fi

instance_dir_arg=$1
stats_dir_arg=$2
shift 2
script_args=("$@")

for arg in "${script_args[@]}"; do
    case "$arg" in
        --output|--output=*|--out-dir|--out-dir=*|--time-limit|--time-limit=*)
            echo "error: do not pass $arg; use --scip-config for SCIP time limits, and let run_experiments.sh set per-instance output paths" >&2
            exit 2
            ;;
    esac
done

if [[ ! -d "$instance_dir_arg" ]]; then
    echo "error: instance directory not found: $instance_dir_arg" >&2
    exit 1
fi

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
instance_dir=$(cd -- "$instance_dir_arg" && pwd)
mkdir -p -- "$stats_dir_arg"
stats_dir=$(cd -- "$stats_dir_arg" && pwd)
jobs_dir="$stats_dir/_jobs"
logs_dir="$stats_dir/_logs"
mkdir -p -- "$jobs_dir" "$logs_dir"

shopt -s nullglob
instances=("$instance_dir"/*.lp)
shopt -u nullglob

if [[ ${#instances[@]} -eq 0 ]]; then
    echo "No top-level .lp instances found in $instance_dir" >&2
    exit 0
fi

partition=${QIP_SBATCH_PARTITION:-opt_int}
constraint=${QIP_SBATCH_CONSTRAINT:-Gold5122}
time_limit=${QIP_SBATCH_TIME:-02:05:00}
exclusive=${QIP_SBATCH_EXCLUSIVE:-1}

submitted=0
skipped=0

for instance_path in "${instances[@]}"; do
    instance_file=$(basename -- "$instance_path")
    instance_name=${instance_file%.lp}
    stats_path="$stats_dir/${instance_name}_stats.json"
    job_script="$jobs_dir/${instance_name}.sh"
    log_path="$logs_dir/${instance_name}.out"

    if [[ -f "$stats_path" ]]; then
        echo "Skipping $instance_file: stats already exist at $stats_path"
        ((skipped += 1))
        continue
    fi

    {
        echo '#!/usr/bin/env bash'
        echo 'set -euo pipefail'
        printf 'exec '
        printf '%q ' \
            julia \
            "--project=$repo_root" \
            "$repo_root/scripts/presolve_lp_scip_stats.jl" \
            "$instance_path" \
            --output \
            "$stats_path" \
            "${script_args[@]}"
        printf '\n'
    } > "$job_script"
    chmod +x "$job_script"

    sbatch_args=(
        --partition "$partition"
        --constraint "$constraint"
        --time "$time_limit"
        --output "$log_path"
    )

    if [[ "$exclusive" != "0" ]]; then
        sbatch_args=(--exclusive "${sbatch_args[@]}")
    fi

    if [[ -n "${QIP_SBATCH_CPUS_PER_TASK:-}" ]]; then
        sbatch_args+=(--cpus-per-task "$QIP_SBATCH_CPUS_PER_TASK")
    fi

    echo "Submitting $instance_file -> $stats_path"
    sbatch "${sbatch_args[@]}" "$job_script"
    ((submitted += 1))
done

echo "Submitted $submitted job(s), skipped $skipped existing stats file(s)."
