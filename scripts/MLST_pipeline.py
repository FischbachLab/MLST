import argparse
import base64
import csv
import glob
import hashlib
import json
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path
from datetime import datetime


# ------------------------------------------------------------
# PubMLST configuration
# ------------------------------------------------------------

PUBMLST_DB = "pubmlst_sepidermidis_seqdef"
SCHEME_ID = "1"

PUBMLST_URL = (
    f"https://rest.pubmlst.org/db/{PUBMLST_DB}"
    f"/schemes/{SCHEME_ID}/sequence"
)

# Same resource without /sequence: returns scheme metadata rather than
# performing a query. Used to record which version of the database a run
# was made against.
PUBMLST_SCHEME_URL = (
    f"https://rest.pubmlst.org/db/{PUBMLST_DB}"
    f"/schemes/{SCHEME_ID}"
)

LOCI = ["arcC", "aroE", "gtr", "mutS", "pyrR", "tpiA", "yqiL"]


# ------------------------------------------------------------
# Submit genome to PubMLST
# ------------------------------------------------------------

REQUEST_TIMEOUT = 60       # seconds to wait for a single API response
MAX_RETRIES = 3            # retry attempts for transient failures
RETRY_BACKOFF = 5          # seconds, doubled after each retry
REQUEST_DELAY = 1          # seconds to wait between successive submissions

# HTTP codes worth retrying (rate limiting / transient server issues)
RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}


def submit_genome(fasta_path, write_log=print):
    """Submit one FASTA assembly to PubMLST and return JSON.

    Retries on rate-limiting (429) and transient server errors, with
    exponential backoff. Raises the last error if all attempts fail.
    """

    with open(fasta_path, "rb") as f:
        fasta_data = f.read()

    encoded_fasta = base64.b64encode(fasta_data).decode("ascii")

    payload = json.dumps({
        "base64": True,
        "details": True,
        "sequence": encoded_fasta
    }).encode("utf-8")

    last_error = None

    for attempt in range(1, MAX_RETRIES + 1):
        request = urllib.request.Request(
            PUBMLST_URL,
            data=payload,
            headers={"Content-Type": "application/json"},
            method="POST"
        )

        try:
            with urllib.request.urlopen(
                request, timeout=REQUEST_TIMEOUT
            ) as response:
                return json.loads(response.read().decode("utf-8"))

        except urllib.error.HTTPError as e:
            last_error = e

            if e.code in RETRYABLE_STATUS_CODES and attempt < MAX_RETRIES:
                wait = RETRY_BACKOFF * (2 ** (attempt - 1))
                write_log(
                    f"    HTTP {e.code} on attempt {attempt}/{MAX_RETRIES}, "
                    f"retrying in {wait}s..."
                )
                time.sleep(wait)
                continue

            raise

        except (urllib.error.URLError, TimeoutError) as e:
            last_error = e

            if attempt < MAX_RETRIES:
                wait = RETRY_BACKOFF * (2 ** (attempt - 1))
                write_log(
                    f"    Network error on attempt {attempt}/{MAX_RETRIES} "
                    f"({e}), retrying in {wait}s..."
                )
                time.sleep(wait)
                continue

            raise

    # Should not be reached, but just in case
    raise last_error


# ------------------------------------------------------------
# Database version
# ------------------------------------------------------------

# Keys worth recording, in the order they should appear. PubMLST does not
# publish a version number for a scheme, so 'last_updated' plus the profile
# count are what actually identify the database state a run was made against.
SCHEME_INFO_KEYS = [
    "description",
    "last_updated",
    "records",
    "profiles",
    "locus_count",
]


def fetch_scheme_info():
    """Return PubMLST scheme metadata as a dict, or None if unavailable.

    Never raises: failing to record the version must not abort a run that
    would otherwise succeed.
    """

    try:
        request = urllib.request.Request(
            PUBMLST_SCHEME_URL,
            headers={"Accept": "application/json"},
            method="GET"
        )

        with urllib.request.urlopen(
            request, timeout=REQUEST_TIMEOUT
        ) as response:
            return json.loads(response.read().decode("utf-8"))

    except Exception:
        return None


def log_scheme_info(info, write_log):
    """Write the database version block to the run log."""

    write_log(f"PubMLST database: {PUBMLST_DB}")
    write_log(f"PubMLST scheme: MLST (scheme {SCHEME_ID})")
    write_log(f"PubMLST endpoint: {PUBMLST_URL}")

    if info is None:
        write_log(
            "PubMLST scheme version: UNAVAILABLE "
            "(metadata request failed; results below are still valid, "
            "but the database state was not recorded)"
        )
        return

    write_log("PubMLST scheme version:")

    reported = False

    for key in SCHEME_INFO_KEYS:
        if key in info:
            write_log(f"    {key}: {info[key]}")
            reported = True

    # Field names on the API can change; fall back to dumping whatever
    # scalar metadata came back rather than silently recording nothing.
    if not reported:
        for key, value in sorted(info.items()):
            if not isinstance(value, (dict, list)):
                write_log(f"    {key}: {value}")


# ------------------------------------------------------------
# Local database (offline mode, --db-dir)
# ------------------------------------------------------------

# Length of the allele prefix used to find candidate hit positions. Alleles
# of one locus mostly share their first bases, so there are few distinct
# anchors to search for, and each hit is then checked against the full
# allele sequences.
ANCHOR_LENGTH = 20

COMPLEMENT = str.maketrans("ACGTRYKMBVDHN", "TGCAYRMKVBHDN")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    """Return a list of (name, uppercase sequence) records."""

    records = []
    name = None
    chunks = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                name = (line[1:].split() or [""])[0]
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        records.append((name, "".join(chunks).upper()))

    return records


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_local_db(db_dir):
    """Load a scheme snapshot written by download_pubmlst_db.py.

    Returns a dict with the manifest, per-locus allele indexes, and the
    profile table. Raises ValueError if the snapshot is partial, corrupted,
    or for a different scheme than LOCI.
    """

    db_dir = Path(db_dir)

    if not (db_dir / ".install_complete").exists():
        raise ValueError(
            f"{db_dir} has no .install_complete marker; the download is "
            f"missing or partial (re-run download_pubmlst_db.py)"
        )

    manifest = json.loads((db_dir / "manifest.json").read_text())

    if set(manifest.get("loci", [])) != set(LOCI):
        raise ValueError(
            f"Local database loci {manifest.get('loci')} do not match "
            f"this script's LOCI {LOCI}"
        )

    for rel_path, expected in manifest.get("sha256", {}).items():
        if sha256(db_dir / rel_path) != expected:
            raise ValueError(
                f"Checksum mismatch for {db_dir / rel_path}; the local "
                f"database has been modified or corrupted"
            )

    # Per locus: {allele sequence: allele id}, the set of allele lengths,
    # and the set of anchor prefixes
    loci = {}

    for locus in LOCI:
        alleles = {}

        for header, seq in read_fasta(db_dir / "alleles" / f"{locus}.fasta"):
            # Headers look like ">arcC_12"
            alleles[seq] = header.rsplit("_", 1)[-1]

        loci[locus] = {
            "alleles": alleles,
            "lengths": sorted({len(s) for s in alleles}),
            "anchors": {s[:ANCHOR_LENGTH] for s in alleles},
        }

    # profiles.tsv: ST, <loci...>, [clonal_complex, ...]
    profiles = {}

    with open(db_dir / "profiles.tsv", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            profiles[tuple(row[locus] for locus in LOCI)] = row["ST"]

    return {
        "dir": db_dir,
        "manifest": manifest,
        "loci": loci,
        "profiles": profiles,
    }


def type_genome_local(fasta_path, db):
    """Call exact allele matches and ST against a local database.

    Returns a dict shaped like the PubMLST /sequence response
    ("exact_matches" and "fields"), so extract_results handles both modes.
    """

    contigs = read_fasta(fasta_path)

    if not contigs:
        raise ValueError(f"No sequences found in {fasta_path}")

    exact_matches = {}

    for locus, index in db["loci"].items():
        hits = []

        for contig_name, forward in contigs:
            for strand, seq in (("+", forward),
                                ("-", reverse_complement(forward))):
                for anchor in index["anchors"]:
                    pos = seq.find(anchor)

                    while pos != -1:
                        for length in index["lengths"]:
                            allele_id = index["alleles"].get(
                                seq[pos:pos + length]
                            )

                            if allele_id is not None:
                                # 1-based forward-strand coordinates
                                start = (pos if strand == "+"
                                         else len(seq) - pos - length)
                                hits.append({
                                    "allele_id": allele_id,
                                    "contig": contig_name,
                                    "start": start + 1,
                                    "end": start + length,
                                    "orientation": strand,
                                })

                        pos = seq.find(anchor, pos + 1)

        if hits:
            exact_matches[locus] = hits

    fields = {}

    if all(locus in exact_matches for locus in LOCI):
        profile = tuple(
            exact_matches[locus][0]["allele_id"] for locus in LOCI
        )
        st = db["profiles"].get(profile)

        if st is not None:
            fields["ST"] = st

    return {"exact_matches": exact_matches, "fields": fields}


def log_local_db_info(db, write_log):
    """Write the local database version block to the run log."""

    manifest = db["manifest"]

    write_log(f"Local database: {db['dir']}")
    write_log(f"PubMLST database: {manifest.get('database')}")
    write_log(f"PubMLST scheme: {manifest.get('scheme_description')} "
              f"(scheme {manifest.get('scheme_id')})")
    write_log("PubMLST scheme version:")
    write_log(f"    last_updated: {manifest.get('scheme_last_updated')}")
    write_log(f"    records: {manifest.get('scheme_records')}")
    write_log(f"    profiles: {manifest.get('n_profiles')}")
    write_log(f"    alleles: {manifest.get('n_alleles')}")
    write_log(f"    downloaded_at_utc: {manifest.get('downloaded_at_utc')}")
    write_log("    checksums: verified")


# ------------------------------------------------------------
# Extract MLST results
# ------------------------------------------------------------

def extract_results(data):
    """Extract allele calls and ST from PubMLST response."""

    row = {}

    exact_matches = data.get("exact_matches", {})

    for locus in LOCI:
        matches = exact_matches.get(locus, [])

        if matches:
            row[locus] = matches[0].get("allele_id", "")
        else:
            row[locus] = ""

    fields = data.get("fields", {})
    row["ST"] = fields.get("ST", "")

    missing = [locus for locus in LOCI if not row[locus]]

    if not missing and row["ST"]:
        row["status"] = "complete"
        row["notes"] = ""
    elif missing:
        row["status"] = "incomplete"
        row["notes"] = "Missing: " + ", ".join(missing)
    else:
        row["status"] = "unknown"
        row["notes"] = "No ST returned"

    return row


# ------------------------------------------------------------
# Main program
# ------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Determine S. epidermidis MLST sequence types "
            "using the PubMLST API."
        )
    )

    parser.add_argument(
        "--input",
        nargs="+",
        default=["genomes"],
        help=(
            "One or more FASTA files, directories, and/or glob patterns "
            "(e.g. --input genomes/ or --input strainA.fasta strainB.fasta "
            "or --input 'batch1/*.fasta' 'batch2/*.fa'). Default: genomes"
        )
    )

    parser.add_argument(
        "--output",
        default="ST-outputs/MLST_results.csv",
        help=(
            "Output CSV filename "
            "(default: ST-outputs/MLST_results.csv)"
        )
    )

    parser.add_argument(
        "--details",
        action="store_true",
        help="Save raw PubMLST JSON responses (or local hit details)"
    )

    parser.add_argument(
        "--db-dir",
        help=(
            "Type against a local scheme snapshot made by "
            "download_pubmlst_db.py instead of the PubMLST API "
            "(no network access needed)"
        )
    )

    args = parser.parse_args()

    output_file = Path(args.output)

    # Create output directory
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Find FASTA files. Each --input entry can be a single file, a
    # directory (searched for *.fasta/*.fa/*.fna), or a glob pattern.
    fasta_extensions = ("*.fasta", "*.fa", "*.fna")
    fasta_files = []

    for entry in args.input:
        path = Path(entry)

        if path.is_dir():
            for pattern in fasta_extensions:
                fasta_files.extend(path.glob(pattern))
        elif path.is_file():
            fasta_files.append(path)
        else:
            # Treat as a glob pattern (e.g. "batch1/*.fasta")
            matches = [Path(p) for p in sorted(glob.glob(entry))]
            fasta_files.extend(matches)

    # De-duplicate while preserving order, then sort for a stable run order
    seen = set()
    unique_fasta_files = []

    for f in fasta_files:
        resolved = f.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique_fasta_files.append(f)

    fasta_files = sorted(unique_fasta_files, key=lambda p: p.name)

    # --------------------------------------------------------
    # Create run log
    # --------------------------------------------------------

    timestamp = datetime.now()
    timestamp_string = timestamp.strftime("%Y-%m-%d_%H-%M-%S")

    log_file = (
        output_file.parent /
        f"{output_file.stem}_run_{timestamp_string}.log"
    )

    with open(log_file, "w") as log:

        def write_log(message):
            print(message)
            log.write(message + "\n")

        write_log("=" * 70)
        write_log("PubMLST S. epidermidis MLST run")
        write_log("=" * 70)

        write_log(f"Run date/time: {timestamp.isoformat(timespec='seconds')}")
        write_log(f"Input source(s): {', '.join(args.input)}")
        write_log(f"Output CSV: {output_file}")

        # Record which version of the reference database this run used, so the
        # results stay interpretable after PubMLST adds new alleles/profiles.
        local_db = None

        if args.db_dir:
            write_log("Mode: local database")

            try:
                local_db = load_local_db(args.db_dir)
            except Exception as e:
                write_log(f"ERROR: Cannot use local database: {e}")
                return 1

            scheme_info = local_db["manifest"]
            log_local_db_info(local_db, write_log)
        else:
            write_log("Mode: PubMLST API")
            scheme_info = fetch_scheme_info()
            log_scheme_info(scheme_info, write_log)

        write_log(f"MLST loci: {', '.join(LOCI)}")
        write_log(f"Details requested: {args.details}")
        write_log("")

        # ----------------------------------------------------
        # Check for FASTA files
        # ----------------------------------------------------

        if not fasta_files:
            write_log(f"ERROR: No FASTA files found for input: {args.input}")
            return

        write_log(f"Found {len(fasta_files)} FASTA file(s).")
        write_log("")

        results = []

        # Directory for raw JSON files
        details_dir = output_file.parent / "details"

        if args.details:
            details_dir.mkdir(parents=True, exist_ok=True)

            # Machine-readable copy of the database version alongside the
            # per-isolate responses
            if scheme_info is not None:
                with open(details_dir / "_scheme_info.json", "w") as f:
                    json.dump(scheme_info, f, indent=2)

        # ----------------------------------------------------
        # Process genomes
        # ----------------------------------------------------

        for i, fasta_path in enumerate(fasta_files, start=1):

            isolate = fasta_path.stem

            write_log(
                f"[{i}/{len(fasta_files)}] Processing {isolate}..."
            )

            try:
                if local_db is not None:
                    data = type_genome_local(fasta_path, local_db)
                else:
                    data = submit_genome(fasta_path, write_log=write_log)

                # Save raw JSON
                if args.details:
                    json_file = details_dir / f"{isolate}.json"

                    with open(json_file, "w") as f:
                        json.dump(data, f, indent=2)

                result = extract_results(data)
                result["isolate"] = isolate

                results.append(result)

                if result["status"] == "complete":
                    write_log(
                        f"    ST: {result['ST']} "
                        f"(complete)"
                    )

                else:
                    write_log(
                        f"    Status: {result['status']}"
                    )
                    write_log(
                        f"    Notes: {result['notes']}"
                    )

            except urllib.error.HTTPError as e:

                message = f"HTTP {e.code}"

                write_log(
                    f"    ERROR: PubMLST returned {message}"
                )

                results.append({
                    "isolate": isolate,
                    "arcC": "",
                    "aroE": "",
                    "gtr": "",
                    "mutS": "",
                    "pyrR": "",
                    "tpiA": "",
                    "yqiL": "",
                    "ST": "",
                    "status": "API error",
                    "notes": message
                })

            except Exception as e:

                write_log(f"    ERROR: {e}")

                results.append({
                    "isolate": isolate,
                    "arcC": "",
                    "aroE": "",
                    "gtr": "",
                    "mutS": "",
                    "pyrR": "",
                    "tpiA": "",
                    "yqiL": "",
                    "ST": "",
                    "status": "error",
                    "notes": str(e)
                })

            # Be polite to the API between submissions
            if local_db is None and i < len(fasta_files):
                time.sleep(REQUEST_DELAY)

        # ----------------------------------------------------
        # Write CSV
        # ----------------------------------------------------

        fieldnames = [
            "isolate",
            "arcC",
            "aroE",
            "gtr",
            "mutS",
            "pyrR",
            "tpiA",
            "yqiL",
            "ST",
            "status",
            "notes"
        ]

        with open(output_file, "w", newline="") as f:

            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames
            )

            writer.writeheader()
            writer.writerows(results)

        # ----------------------------------------------------
        # Summary
        # ----------------------------------------------------

        complete = sum(
            r["status"] == "complete"
            for r in results
        )

        incomplete = sum(
            r["status"] == "incomplete"
            for r in results
        )

        errors = sum(
            r["status"] in ["error", "API error"]
            for r in results
        )

        write_log("")
        write_log("=" * 70)
        write_log("SUMMARY")
        write_log("=" * 70)

        write_log(f"Total genomes: {len(results)}")
        write_log(f"Complete MLST results: {complete}")
        write_log(f"Incomplete results: {incomplete}")
        write_log(f"Errors: {errors}")
        write_log("")
        write_log(f"Results CSV: {output_file}")
        write_log(f"Run log: {log_file}")
        write_log("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
