#!/bin/bash
# Load the repo-root .env into the environment so the preprocessing scripts (and sbatch, which
# propagates the submitting shell's environment) can see it.  Usage, from anywhere:
#   source /path/to/repo/env.sh
_env_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -f "${_env_here}/.env" ]; then
    echo "No .env found at ${_env_here}/.env. Copy .env.example to .env and fill it in." >&2
    return 1 2>/dev/null || exit 1
fi
set -a
# shellcheck disable=SC1091
source "${_env_here}/.env"
set +a
unset _env_here
