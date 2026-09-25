# MLST Analysis

Assigns *Staphylococcus epidermidis* MLST sequence types (STs) to assembled
genomes by submitting each assembly to the [PubMLST REST API](https://rest.pubmlst.org)
and collecting the allele calls into a single CSV.

No BLAST, no dependencies beyond the Python standard library. By default it
requires outbound internet access to `rest.pubmlst.org`; with `--db-dir` it
instead types against a local snapshot made by `download_pubmlst_db.py` and
needs no network at all (see [Local database mode](#local-database-mode)).

The local database path on efs using an ec2 instance:
```bash
 /mnt/efs/databases/pubmlst_sepidermidis/
 ```

## Quick start

```bash
# every FASTA in ./genomes/ -> ST-outputs/MLST_results.csv
python3 MLST_pipeline.py

# explicit inputs and output
python3 MLST_pipeline.py --input assemblies/ --output results/epi_mlst.csv

# keep the raw API responses
python3 MLST_pipeline.py --input assemblies/ --details

# offline: snapshot the scheme once, then type against it
python3 download_pubmlst_db.py --outdir pubmlst_sepidermidis
python3 MLST_pipeline.py --input assemblies/ --db-dir  /mnt/efs/databases/pubmlst_sepidermidis/
```

## What it does

For each FASTA file it:

1. base64-encodes the whole assembly and POSTs it to the PubMLST
   `sepidermidis_seqdef` scheme-1 sequence endpoint
2. reads back the exact allele match for each of the seven MLST loci
3. records the ST that PubMLST assigns to that allele profile
4. writes one row per isolate to the output CSV, plus a timestamped run log

Submissions are sequential with a 1-second pause between genomes.

## Requirements

- Python 3.6+ — standard library only (`argparse`, `urllib`, `json`, `csv`, `pathlib`)
- Network access to `https://rest.pubmlst.org` (no API key needed for this
  public read endpoint)

## Options

| Option | Default | Description |
|---|---|---|
| `--input` | `genomes` | One or more FASTA files, directories, and/or glob patterns |
| `--output` | `ST-outputs/MLST_results.csv` | Output CSV path; parent directories are created |
| `--details` | off | Also save each raw JSON response to `<output-dir>/details/<isolate>.json` |
| `--db-dir` | unset | Type against a local snapshot from `download_pubmlst_db.py` instead of the API |

### How `--input` is resolved

Each entry is handled independently:

- **a directory** → searched for `*.fasta`, `*.fa`, `*.fna` (**non-recursive** — subdirectories are not descended into)
- **an existing file** → used as-is, whatever its extension
- **anything else** → treated as a shell glob pattern (quote it so your shell doesn't expand it first)

Entries can be mixed, duplicates are removed by resolved path, and files are
processed in alphabetical order of filename.

```bash
python3 MLST_pipeline.py --input genomes/
python3 MLST_pipeline.py --input strainA.fasta strainB.fasta
python3 MLST_pipeline.py --input 'batch1/*.fasta' 'batch2/*.fa'
python3 MLST_pipeline.py --input assemblies/ extras/one_more.fna
```

## Output

### `MLST_results.csv`

```bash
isolate,arcC,aroE,gtr,mutS,pyrR,tpiA,yqiL,ST,status,notes
LM088,16,1,2,1,2,1,1,184,complete,
LM087,8,2,2,4,9,6,9,72,complete,
```

| Column | Description |
|---|---|
| `isolate` | FASTA filename without its extension |
| `arcC` `aroE` `gtr` `mutS` `pyrR` `tpiA` `yqiL` | allele ID per locus, blank if no exact match |
| `ST` | sequence type, blank if PubMLST assigns none |
| `status` | `complete`, `incomplete`, `unknown`, `error`, or `API error` |
| `notes` | which loci are missing, or the error message |

Status meanings:

| Status | Meaning |
|---|---|
| `complete` | all 7 loci called **and** an ST assigned |
| `incomplete` | one or more loci had no exact match — `notes` lists them |
| `unknown` | all 7 loci called but no ST returned; usually a **novel allele combination** not yet in PubMLST, worth following up |
| `API error` | PubMLST returned an HTTP error for this genome |
| `error` | anything else (unreadable file, malformed response, …) |

A failure on one genome never aborts the run — it becomes a row with empty
allele columns and the reason in `notes`.

### Run log

`<output-dir>/<output-stem>_run_<YYYY-MM-DD_HH-MM-SS>.log` — the same text
printed to the terminal, including the parameters used, the database version
(below), per-genome progress, and a summary of complete / incomplete / error
counts. A new log is written per run; the CSV is overwritten.

### Database version

PubMLST does not publish a version number for a scheme, so each run queries
`/db/<db>/schemes/<id>` and records the scheme metadata — `last_updated` and
the profile count being the fields that actually pin down the database state:

```
======================================================================
PubMLST S. epidermidis MLST run
======================================================================
Run date/time: 2026-09-24T22:25:02
Input source(s): assemblies/
Output CSV: ST-outputs/MLST_results.csv
Mode: PubMLST API
PubMLST database: pubmlst_sepidermidis_seqdef
PubMLST scheme: MLST (scheme 1)
PubMLST endpoint: https://rest.pubmlst.org/db/pubmlst_sepidermidis_seqdef/schemes/1/sequence
PubMLST scheme version:
    description: MLST
    last_updated: 2026-09-18
    records: 1366
    profiles: https://rest.pubmlst.org/db/pubmlst_sepidermidis_seqdef/schemes/1/profiles
    locus_count: 7
MLST loci: arcC, aroE, gtr, mutS, pyrR, tpiA, yqiL
```

This matters because PubMLST is a live database: an assembly that returns
`unknown` today may be assigned an ST once someone deposits the novel profile,
and ST assignments are only reproducible relative to a database state. Keep
the log with the CSV.

If the metadata request fails, the run continues and the log records
`PubMLST scheme version: UNAVAILABLE` rather than aborting. With `--details`,
the raw metadata is also written to `<output-dir>/details/_scheme_info.json`.

### `details/` (with `--details`)

Full JSON per isolate, including partial and inexact matches that the CSV
discards. Useful for investigating `incomplete` and `unknown` rows.

## Local database mode

`--db-dir DIR` replaces the API call with an exact-match search against a
snapshot written by `download_pubmlst_db.py` (`profiles.tsv`,
`alleles/<locus>.fasta`, `manifest.json`, `.install_complete`).

For each locus, every allele is searched for as an exact substring of each
contig, on both strands; the resulting allele profile is then looked up in
`profiles.tsv`. The CSV columns and status values are identical to API mode.
It takes well under a second per genome and there is no delay between genomes.

Before typing anything, the snapshot is checked:

- `.install_complete` must exist (a partial download is rejected)
- the manifest's loci must match `LOCI`
- every file's SHA-256 must match `manifest.json`

If any check fails the run logs `ERROR: Cannot use local database: …` and
**exits 1** without writing a CSV.

The log records the snapshot instead of querying PubMLST: the directory,
`last_updated`, profile and allele counts, and `downloaded_at_utc`. With
`--details`, `details/_scheme_info.json` is a copy of the manifest, and each
`details/<isolate>.json` lists every exact hit per locus with contig,
1-based start/end, and orientation.

Differences from API mode:

- Only exact, full-length matches are reported; there are no partial or
  inexact matches in `details/`.
- An allele split across two contigs is not found (it is reported missing).
- If a locus has several distinct exact hits (e.g. a mixed sample), the first
  one found is used for the ST, as in API mode; check `details/` for such cases.
- Results reflect the snapshot, not the live database — refresh it with
  `download_pubmlst_db.py --force` to pick up newly deposited STs.

## Scheme configuration

Hard-coded near the top of the script:

```python
PUBMLST_DB = "pubmlst_sepidermidis_seqdef"
SCHEME_ID  = "1"
LOCI       = ["arcC", "aroE", "gtr", "mutS", "pyrR", "tpiA", "yqiL"]
```

To target a different organism, change all three together — the loci list must
match the scheme, or every locus will come back blank.

## Network behavior

| Setting | Value |
|---|---|
| `REQUEST_TIMEOUT` | 60 s per request |
| `MAX_RETRIES` | 3 attempts |
| `RETRY_BACKOFF` | 5 s, doubling (5 s, 10 s) |
| Retried on | HTTP 429, 500, 502, 503, 504 |
| `REQUEST_DELAY` | 1 s between genomes |

Other HTTP errors (e.g. 400, 404) are not retried and are recorded against
that genome.

## Known limitations

- **Gzipped FASTA is not supported.** The file is read as raw bytes and
  base64-encoded, so a `.gz` would be submitted as binary and rejected.
  Decompress first.
- **Directory search is not recursive.** `--input genomes/` misses
  `genomes/sub/x.fasta`. Note also that `'genomes/**/*.fasta'` descends exactly
  one level, not arbitrarily deep; use `find` and pass explicit paths if you
  need a deep tree.
- **Duplicate basenames collide.** `a/S1.fasta` and `b/S1.fasta` are both
  submitted, but both produce `isolate = S1`, so the CSV gets two rows with the
  same name.
- **Exits 0 on nearly every failure** (the one exception is a rejected
  `--db-dir`), including when no FASTA files are found (it logs
  `ERROR: No FASTA files found` and stops). If you wrap this in a pipeline,
  check the CSV row count rather than the exit status.
- **Whole assemblies are uploaded** to a third-party server in API mode. Fine
  for public isolate genomes; for anything restricted, use `--db-dir`.
- **Serial, one genome at a time.** In API mode, roughly 1 s plus API
  round-trip each, so a few hundred genomes takes minutes, not seconds. The
  delay is deliberate politeness toward a free public API — raising it is not
  recommended.
- API mode requires internet egress, so it will not run on an isolated HPC
  compute node or an AWS Batch queue without a NAT route. Use `--db-dir` there,
  with the snapshot on shared storage.


## Contacts:

[Allison Weakley](),
[Xiandong Meng]()