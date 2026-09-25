#!/usr/bin/env python3
"""Download a PubMLST scheme (profiles + allele FASTAs) for offline use.

PubMLST serves the underlying data as flat files, so a whole MLST scheme can
be snapshotted locally and then installed on EFS / shared storage. This keeps
ST calling reproducible: PubMLST is a live database, and an assembly that
returns no ST today may be assigned one after someone deposits the profile.

Downloads, per scheme:
    profiles.tsv                 ST -> allele profile table
    alleles/<locus>.fasta        every allele sequence for that locus
    manifest.json                scheme metadata, download date, checksums
    .install_complete            written last; its presence means "not partial"

Usage:
    python3 download_pubmlst_db.py --outdir pubmlst_sepidermidis
    python3 download_pubmlst_db.py --db pubmlst_saureus_seqdef --scheme 1 \
        --outdir pubmlst_saureus

Standard library only.
"""

import argparse
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REST_ROOT = "https://rest.pubmlst.org"

REQUEST_TIMEOUT = 60
MAX_RETRIES = 3
RETRY_BACKOFF = 5
REQUEST_DELAY = 1

RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}


def fetch(url, binary=False):
    """GET a URL with retry/backoff. Returns bytes or str."""

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            request = urllib.request.Request(url, method="GET")
            with urllib.request.urlopen(
                request, timeout=REQUEST_TIMEOUT
            ) as response:
                payload = response.read()
                return payload if binary else payload.decode("utf-8")

        except urllib.error.HTTPError as e:
            if e.code in RETRYABLE_STATUS_CODES and attempt < MAX_RETRIES:
                wait = RETRY_BACKOFF * (2 ** (attempt - 1))
                print(f"    HTTP {e.code}, retrying in {wait}s "
                      f"({attempt}/{MAX_RETRIES})", file=sys.stderr)
                time.sleep(wait)
                continue
            raise

        except (urllib.error.URLError, TimeoutError) as e:
            if attempt < MAX_RETRIES:
                wait = RETRY_BACKOFF * (2 ** (attempt - 1))
                print(f"    Network error ({e}), retrying in {wait}s "
                      f"({attempt}/{MAX_RETRIES})", file=sys.stderr)
                time.sleep(wait)
                continue
            raise


def locus_name(locus_url):
    """https://rest.pubmlst.org/db/<db>/loci/arcC -> arcC"""
    return locus_url.rstrip("/").rsplit("/", 1)[-1]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--db", default="pubmlst_sepidermidis_seqdef",
                   help="PubMLST seqdef database "
                        "(default: pubmlst_sepidermidis_seqdef)")
    p.add_argument("--scheme", default="1",
                   help="Scheme id; 1 is MLST for most organisms (default: 1)")
    p.add_argument("--outdir", required=True,
                   help="Destination directory (created if absent)")
    p.add_argument("--force", action="store_true",
                   help="Re-download even if .install_complete is present")
    args = p.parse_args(argv)

    outdir = Path(args.outdir)
    marker = outdir / ".install_complete"

    if marker.exists() and not args.force:
        print(f"Already downloaded: {outdir}")
        print(marker.read_text())
        return 0

    scheme_url = f"{REST_ROOT}/db/{args.db}/schemes/{args.scheme}"

    print(f"Scheme: {scheme_url}")
    scheme = json.loads(fetch(scheme_url))

    loci = [locus_name(u) for u in scheme.get("loci", [])]
    if not loci:
        print(f"ERROR: scheme {args.scheme} in {args.db} lists no loci.",
              file=sys.stderr)
        return 1

    print(f"  description : {scheme.get('description', '?')}")
    print(f"  last_updated: {scheme.get('last_updated', '?')}")
    print(f"  records     : {scheme.get('records', '?')}")
    print(f"  loci ({len(loci)}) : {', '.join(loci)}")
    print()

    alleles_dir = outdir / "alleles"
    alleles_dir.mkdir(parents=True, exist_ok=True)

    files = {}

    # ── Profiles ────────────────────────────────────────────────────────────
    # profiles_csv is actually TAB-delimited despite the endpoint name.
    profiles_path = outdir / "profiles.tsv"
    print("Downloading profiles ...")
    profiles_path.write_text(fetch(f"{scheme_url}/profiles_csv"))
    n_profiles = max(0, len(profiles_path.read_text().splitlines()) - 1)
    print(f"  profiles.tsv  ({n_profiles} profiles)")
    files["profiles.tsv"] = sha256(profiles_path)

    # ── Allele sequences ────────────────────────────────────────────────────
    print("\nDownloading alleles ...")
    total_alleles = 0

    for locus in loci:
        url = f"{REST_ROOT}/db/{args.db}/loci/{locus}/alleles_fasta"
        dest = alleles_dir / f"{locus}.fasta"
        dest.write_text(fetch(url))

        n = dest.read_text().count(">")
        total_alleles += n
        print(f"  {locus:<10} {n:>6} alleles")

        files[f"alleles/{locus}.fasta"] = sha256(dest)
        time.sleep(REQUEST_DELAY)

    # ── Manifest ────────────────────────────────────────────────────────────
    downloaded_at = datetime.now(timezone.utc).isoformat(timespec="seconds")

    manifest = {
        "database": args.db,
        "scheme_id": args.scheme,
        "scheme_description": scheme.get("description"),
        "scheme_last_updated": scheme.get("last_updated"),
        "scheme_records": scheme.get("records"),
        "loci": loci,
        "n_profiles": n_profiles,
        "n_alleles": total_alleles,
        "downloaded_at_utc": downloaded_at,
        "source": scheme_url,
        "sha256": files,
    }

    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    # Written last: its presence means the directory is whole, not half-copied
    marker.write_text(
        f"downloaded {downloaded_at} from {scheme_url}\n"
        f"scheme last_updated {scheme.get('last_updated')}\n"
    )

    print(f"\nDone: {outdir}")
    print(f"  {n_profiles} profiles, {total_alleles} alleles across "
          f"{len(loci)} loci")
    print(f"  manifest.json records scheme last_updated="
          f"{scheme.get('last_updated')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
