#!/usr/bin/env python3
"""
extract_aabc_paths.py

Extracts Aspera source paths (.zip only, no .md5) from a BALSA
downloadToken API JSON response.

Usage:
  1. In Chrome DevTools (Network tab), find the "downloadToken" request
  2. Right-click on it > Copy > Copy Response
  3. Save the response to a file:
       pbpaste > /tmp/balsa_response.json    (macOS)
  4. Run this script:
       python3 extract_aabc_paths.py /tmp/balsa_response.json

  The output will be written to:
       $PROJECT_DIR/hcp_tractometry/code/aabc_filelist.txt

Alternative: If you saved the response on your Mac, you can scp it:
    scp /tmp/balsa_response.json <user>@<cluster>:/tmp/balsa_response.json
    Then run this script on the cluster.
"""

import json
import sys
import os

OUTPUT_PATH = os.environ["PROJECT_DIR"] + "/hcp_tractometry/code/aabc_filelist.txt"


def extract_paths(json_file):
    """Parse BALSA transfer spec JSON and extract .zip source paths."""
    with open(json_file, "r") as f:
        data = json.load(f)

    # The BALSA downloadToken response has a "transfer_specs" array.
    # Each element has a "transfer_spec" with a "paths" array.
    # Each path entry has a "source" key.
    paths = []

    # Handle different possible JSON structures
    if isinstance(data, dict):
        # Structure 1: {"transfer_specs": [{"transfer_spec": {"paths": [...]}}]}
        transfer_specs = data.get("transfer_specs", [])
        if not transfer_specs:
            # Structure 2: top-level {"paths": [...]}
            transfer_specs = [{"transfer_spec": data}]

        for ts_wrapper in transfer_specs:
            ts = ts_wrapper.get("transfer_spec", ts_wrapper)
            for p in ts.get("paths", []):
                src = p.get("source", "")
                if src:
                    paths.append(src)
    elif isinstance(data, list):
        # Structure 3: just a list of path strings
        paths = [p for p in data if isinstance(p, str)]

    # Filter: only .zip files, exclude .md5
    zip_paths = sorted(set(
        p for p in paths
        if p.endswith(".zip") and not p.endswith(".md5")
    ))

    return zip_paths


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("ERROR: Please provide the path to the BALSA JSON response file.")
        sys.exit(1)

    json_file = sys.argv[1]
    if not os.path.exists(json_file):
        print(f"ERROR: File not found: {json_file}")
        sys.exit(1)

    print(f"Reading BALSA response from: {json_file}")
    zip_paths = extract_paths(json_file)

    if not zip_paths:
        print("WARNING: No .zip paths found in the JSON response!")
        print("The JSON structure may be different than expected.")
        print("Dumping top-level keys for debugging...")
        with open(json_file, "r") as f:
            data = json.load(f)
        if isinstance(data, dict):
            print(f"  Top-level keys: {list(data.keys())}")
            for k, v in data.items():
                if isinstance(v, list):
                    print(f"  '{k}' is a list with {len(v)} items")
                    if v:
                        print(f"    First item type: {type(v[0])}")
                        if isinstance(v[0], dict):
                            print(f"    First item keys: {list(v[0].keys())}")
        sys.exit(1)

    # Write output
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        for p in zip_paths:
            f.write(p + "\n")

    print(f"Wrote {len(zip_paths)} .zip paths to: {OUTPUT_PATH}")
    print(f"First path: {zip_paths[0]}")
    print(f"Last path:  {zip_paths[-1]}")

    # Print some stats
    subjects = set()
    for p in zip_paths:
        parts = p.split("/")
        for part in parts:
            if part.startswith("HCA"):
                subjects.add(part)
                break
    print(f"Unique subjects: {len(subjects)}")


if __name__ == "__main__":
    main()
