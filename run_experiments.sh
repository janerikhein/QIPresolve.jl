#!/usr/bin/env bash

set -euo pipefail

instance_dir="$1"
results_dir="$2"
args=("${@:3}")

project_dir="/home/htc/jehein/QIPresolve.jl"

# Create output directories
mkdir -p "$results_dir/logs"
mkdir -p "$results_dir/stats"

instances=("$instance_dir"/*.lp)

for instance_path in "${instances[@]}"; do
    instance_name="$(basename "$instance_path")"

    echo "Instance: $instance_name"

    log_file="$results_dir/logs/$instance_name.log"

    if [ ! -f "$log_file" ]; then

        # Build the Julia command with shell-safe quoting
        printf -v command \
            'cd %q && julia --project=. scripts/presolve_lp_scip_stats.jl %q --out-dir %q' \
            "$project_dir" \
            "$instance_path" \
            "$results_dir/stats"

        # Append arbitrary arguments passed after INSTANCE_DIR RESULTS_DIR
        for arg in "${args[@]}"; do
            printf -v quoted_arg '%q' "$arg"
            command+=" $quoted_arg"
        done

        echo "Submitting: $instance_name"

        #--exclusive \
        sbatch \
            --partition=small \
            --constraint=Gold6338 \
            --time=02:05:00 \
            --nodes=1 \
            --cpu-freq=medium-medium:Performance \
            --threads-per-core=1 \
            --cpus-per-task=1 \
            --exclude=htc-cmp[101-102],htc-cmp104,htc-cmp126,htc-cmp[145-148] \
            --memory=12000 \
            --job-name= xx \
            --mail-user="hein@zib.de" \
            --mail-type=END,FAIL \
            --wrap="$command" \
            --output="$log_file" \
            --error="$error_file"

    else
        echo "Already submitted/solved: $instance_name"
    fi
done