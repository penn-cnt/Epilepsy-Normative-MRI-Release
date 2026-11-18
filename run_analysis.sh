#!/bin/bash
# Run the five analysis notebooks in order. Executed copies, tables and figures go to outputs/.
# supp_table_discrimination runs last because it reads tables the other notebooks write.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/analysis"

for nb in composite_figures metadata_summary mde_power_analysis abnormality_burden_summary supp_table_discrimination; do
    echo "Running ${nb}"
    python -m nbconvert --to notebook --execute --output "../outputs/${nb}_executed.ipynb" "${nb}.ipynb"
done
