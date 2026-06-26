#!/usr/bin/env bash
# Downloads the pre-processed gPhage dataset from Zenodo (10.5281/zenodo.20854186)
# into dataSet/gPhage_data/, verifies MD5 checksums and unzips the files.
set -euo pipefail

ZENODO_RECORD_ID="20854186"
BASE_URL="https://zenodo.org/records/${ZENODO_RECORD_ID}/files"
DEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/dataSet/gPhage_data"

declare -A FILES_MD5=(
    ["DNA_inserts_all_patients_groups.fastq.zip"]="9fcc12426e2ebf6c520edc5dd42e818f"
    ["DNA_insert_IDs_by_group.tsv.zip"]="785fb95877242e021d3b49772b5209aa"
    ["peptides_encoded_by_DNA_inserts.fasta.zip"]="a3ac336dc898e18a5571e5851e40f9bb"
    ["peptides_encoded_by_DNA_inserts.tsv.zip"]="29e64a654489eea63a33690d9fe067fe"
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
