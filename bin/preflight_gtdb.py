#!/usr/bin/env python3
"""Preflight GTDB version check for HiFiTaxa.

Before the pipeline runs, compare the locally-built GTDB BLCA database against
the LATEST GTDB release on the website. Behaviour:

  * local == latest                -> proceed.
  * newer release available        -> ASK whether to download/parse/build it;
                                       if built, ASK whether to delete the
                                       previous raw download.
  * GTDB website not reachable     -> just proceed (use whatever DB exists).
  * no local DB and reachable      -> offer to build the latest.
  * no local DB and unreachable    -> error (nothing to run against).

Exit code 0 = OK to proceed; non-zero = abort.
"""
import argparse, os, re, subprocess, sys
from pathlib import Path

LATEST_VERSION_URL = "https://data.gtdb.ecogenomic.org/releases/latest/VERSION.txt"
SSU_URL_TMPL = "https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_all/ssu_all_r{rel}.fna.gz"


def msg(s): print(f"[gtdb-preflight] {s}", flush=True)


def fetch_latest_release(timeout):
    """Return latest GTDB release number (int) or None if unreachable."""
    import urllib.request
    try:
        with urllib.request.urlopen(LATEST_VERSION_URL, timeout=timeout) as r:
            text = r.read().decode("utf-8", "replace")
        m = re.search(r"v?(\d+)", text.strip().splitlines()[0])
        return int(m.group(1)) if m else None
    except Exception as e:
        msg(f"could not reach GTDB website ({e.__class__.__name__}: {e}).")
        return None


def local_release(db_dir, blca_db):
    """Return locally-built release number (int) or None."""
    vf = Path(db_dir) / "GTDB_VERSION.txt"
    if vf.is_file():
        m = re.search(r"(\d+)", vf.read_text())
        if m:
            return int(m.group(1))
    return None


def db_present(blca_db, blca_tax):
    """A usable BLCA DB needs the taxonomy file and a blast index (.n*)."""
    if not Path(blca_tax).is_file():
        return False
    p = Path(blca_db)
    return any(p.parent.glob(p.name + ".n*"))


def ask(question, assume=None):
    """Yes/No prompt. `assume` (True/False) bypasses interaction."""
    if assume is not None:
        msg(f"{question} -> {'yes' if assume else 'no'} (non-interactive)")
        return assume
    try:
        return input(f"[gtdb-preflight] {question} [y/N] ").strip().lower() in ("y", "yes")
    except EOFError:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db-dir", default=os.environ.get("GTDB_DB_DIR", "db"))
    ap.add_argument("--blca-db", default=os.environ.get("BLCA_DB", "db/gtdb_ssu_BLCAparsed.fasta"))
    ap.add_argument("--blca-tax", default=os.environ.get("BLCA_TAX", "db/gtdb_ssu_BLCAparsed.taxonomy"))
    ap.add_argument("--build-script", default=str(Path(__file__).with_name("build_gtdb_blca_db.sh")))
    ap.add_argument("--timeout", type=float, default=15.0)
    ap.add_argument("--assume-yes", action="store_true", help="auto-confirm download/build")
    ap.add_argument("--assume-no", action="store_true", help="never download/build; just proceed")
    ap.add_argument("--min-ref-len", type=int, default=1000,
                    help="when building, remove reference seqs shorter than this (0 = keep all). "
                         "In interactive mode the user is asked; this is the non-interactive default.")
    ap.add_argument("--gtdb-release", default="232",
                    help="GTDB release to target/build (default 232, pinned for reproducibility). "
                         "Use 'latest' to query the GTDB website for the newest release.")
    args = ap.parse_args()

    assume = True if args.assume_yes else (False if args.assume_no else None)
    have_db = db_present(args.blca_db, args.blca_tax)
    local = local_release(args.db_dir, args.blca_db)
    msg(f"local GTDB DB: {'r'+str(local) if local else ('present (unknown release)' if have_db else 'none')}")

    # Target release: a pinned number (reproducible default) or 'latest' (queries GTDB).
    if str(args.gtdb_release).lower() == "latest":
        latest = fetch_latest_release(args.timeout)
        if latest is None:                       # website unreachable -> just proceed (per spec)
            if have_db:
                msg("GTDB website unreachable; proceeding with the existing local database.")
                return 0
            msg("ERROR: no local GTDB database and GTDB website unreachable - cannot proceed.")
            return 2
        msg(f"latest GTDB release on website: r{latest}")
    else:
        latest = int(args.gtdb_release)
        msg(f"target GTDB release (pinned): r{latest}")

    # --- up to date ---
    if have_db and local == latest:
        msg(f"database is up to date (r{latest}). Proceeding.")
        return 0

    # --- newer available (or nothing local) ---
    if have_db:
        msg(f"a newer GTDB release is available: local r{local} -> latest r{latest}.")
    else:
        msg(f"no local GTDB database found; latest available is r{latest}.")

    if not ask(f"Download, parse and build GTDB r{latest} now?", assume):
        if have_db:
            msg("keeping existing database. Proceeding.")
            return 0
        msg("ERROR: no database to run against. Aborting.")
        return 2

    # --- decide reference length filtering (the <1000bp removal) ---
    minlen = args.min_ref_len
    if assume is None:   # interactive: ask the user
        minlen = 1000 if ask("Remove reference sequences shorter than 1000 bp when building? (recommended)") else 0
    msg("reference length filter: " + (f"remove sequences < {minlen} bp" if minlen > 0 else "keep ALL sequences"))

    # --- build ---
    msg(f"building GTDB r{latest} BLCA database into {args.db_dir} (min-len {minlen}) ...")
    rc = subprocess.run(["bash", args.build_script, str(latest), args.db_dir, str(minlen)]).returncode
    if rc != 0:
        msg("ERROR: GTDB database build failed.")
        return rc
    msg(f"GTDB r{latest} database built.")

    # --- offer to delete previous raw download(s) ---
    olds = sorted(Path(args.db_dir).glob("ssu_all_r*.fna.gz"))
    olds = [p for p in olds if f"_r{latest}." not in p.name]
    if olds and ask(f"Delete previous raw download(s): {', '.join(p.name for p in olds)}?", assume):
        for p in olds:
            try:
                p.unlink(); msg(f"deleted {p.name}")
            except OSError as e:
                msg(f"could not delete {p.name}: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
