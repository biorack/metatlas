#!/bin/bash --login
set -ef -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
NAME="metatlas-targeted-$(date --iso-8601)"
BASE_DIR="/global/common/software/m2650"
ENV_DIR="${BASE_DIR}/${NAME}"
ENV_FILE="${SCRIPT_DIR}/env.yaml"

echo "name: $NAME
channels:
- conda-forge
dependencies:
- python=3.8
- pip
- pip:" > "$ENV_FILE"
awk '{ print  "  - " $0 }' requirements.txt >> "$ENV_FILE"

conda env create \
  --prefix "$ENV_DIR" \
  --file "$ENV_FILE"
rm "$ENV_FILE"

cat >"${REPO_DIR}/notebooks/kernels/metatlas-targeted.kernel.json" <<EOL
{
 "argv": [
  "${ENV_DIR}/bin/python",
  "-m",
  "ipykernel_launcher",
  "-f",
  "{connection_file}"
 ],
 "display_name": "Metatlas Targeted",
 "language": "python"
}
EOL
