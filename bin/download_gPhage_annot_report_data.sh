#!/usr/bin/env bash
# Downloads the gPhage Annotation Report Data from Zenodo (10.5281/zenodo.20941568)
# into dataSet/gPhage_annot_report_data/, verifies MD5 checksums and unzips the files.
set -euo pipefail

ZENODO_RECORD_ID="20941568"
BASE_URL="https://zenodo.org/records/${ZENODO_RECORD_ID}/files"
DEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/dataSet/gPhage_annot_report_data"

declare -A FILES_MD5=(
    ["proteome_hits.tsv.zip"]="c2bda4d2d3c17b4cad746259b6d4edf8"
    ["iedb_tcruzi_epitopes_hits.tsv.zip"]="e0f710a7d02dfb1602e40b8ecfecad9f"
    ["iedb_tcruzi_epitopes_1747317719_15052025.csv.zip"]="ac09ac6e3e2cacb467ae05e46767ae97"
    ["iedb_human_epitopes_hits.tsv.zip"]="1134db0998baa6aea51d927e92f6bf9c"
    ["iedb_human_epitopes_1757095864.csv.zip"]="436cffd4c67c5144a8eb4fb61586adb0"
)

mkdir -p "${DEST_DIR}"
cd "${DEST_DIR}"

for filename in "${!FILES_MD5[@]}"; do
    expected_md5="${FILES_MD5[${filename}]}"

    echo "Downloading ${filename}..."
    curl -L -o "${filename}" "${BASE_URL}/${filename}?download=1"

    actual_md5="$(md5sum "${filename}" | awk '{print $1}')"
    if [[ "${actual_md5}" != "${expected_md5}" ]]; then
        echo "ERROR: MD5 mismatch for ${filename} (expected ${expected_md5}, got ${actual_md5})" >&2
        exit 1
    fi

    echo "Unzipping ${filename}..."
    unzip -o "${filename}"
done

echo "Done. Files are available in ${DEST_DIR}"
