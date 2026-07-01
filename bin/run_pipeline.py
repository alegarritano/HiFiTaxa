#!/usr/bin/env python3
"""HiFiTaxa launcher.

- Runs the GTDB preflight (version check / optional build), then the Nextflow
  pipeline with a Snakemake-style tqdm progress bar.
- `--test` runs the bundled example (example/fastq/*.fastq.gz) end-to-end and
  compares the result to the bundled reference (example/reference), printing a
  GREEN/RED verdict.
- On the FIRST interactive run it offers to do a test run first.

Examples:
  run_pipeline.py --input s.tsv --metadata m.tsv --outdir out --qiime2_env <env> --blca_env <env>
  run_pipeline.py --asv_fasta ASVs.fasta --outdir tax_out          # taxonomy-only
  run_pipeline.py --test --qiime2_env <env> --blca_env <env>       # example smoke test
Unknown options are passed straight through to `nextflow run`.
"""
import argparse, glob, os, re, shutil, subprocess, sys, tempfile, threading, time
from pathlib import Path

DEFAULT_FWD = "AGRGTTYGATYMTGGCTCAG"   # 27F
DEFAULT_REV = "AAGTCGTAACAAGGTARCY"    # 1492R
# ITS (fungal) marker defaults — must match the contract / nextflow.config and
# HiFiTaxa_Fungi/bin/denoise_pipeline.sh (FWD / REVRC).
DEFAULT_ITS_FWD = "GTACACACCGCCCGTC"   # 1391F
DEFAULT_ITS_REV = "GCATATHANTAAGSGSAGG"  # revcomp(ITS4ngsUni)
READS_DIR = "00_Reads"                 # drop your .fastq.gz reads here

# Latest UNITE general release we know about. Bump this when UNITE ships a newer
# sh_general_release; the launcher prints it so the user can decide to rebuild.
LATEST_UNITE = "sh_general_release_dynamic_19.02.2025"
# UNITE releases are downloaded by hand from the UNITE repository page: pick the
# 'Current' general FASTA release for Fungi (e.g. version 10.0, 19.02.2025). The
# ITS preflight points the user here, then asks for the downloaded FASTA and
# builds the references from it (no auto-download / version API).
UNITE_REPO_URL = "https://unite.ut.ee/repository.php"

PROJDIR = Path(__file__).resolve().parent.parent
MARKER = PROJDIR / ".first_run_complete"

# ANSI colour for the resource-cost note in the welcome wizard.
ANSI_RED = "\033[1;31m"
ANSI_RST = "\033[0m"


def welcome_banner():
    """Just the static banner header — the interview is in welcome_wizard()."""
    print()
    print("💻  Welcome to HiFiTaxa v1.0  🧬")
    print()
    print("This appears to be a fresh install. We'll set it up in two phases:")
    print("  Phase 1 (now): Set-up on what is useful for you")
    print("  Phase 2 (unattended): install everything, then run a mock community")
    print("                       to validate the install. ~20-60 min total.")
    print()


def _ask_yn(prompt, default="y"):
    """Ask a yes/no question. `default` is 'y' or 'n' — used when the user
    just hits Enter. Returns True for yes."""
    suffix = "[Y/n]" if default == "y" else "[y/N]"
    try:
        ans = input(f"    {prompt} {suffix} ").strip().lower()
    except EOFError:
        ans = ""
    if ans == "":
        return default == "y"
    return ans in ("y", "yes")


def welcome_wizard(args):
    """First-run interview. Asks ALL questions upfront so the install phase
    can then run fully unattended (no surprise prompts in the middle of a
    30-min env build). Mutates `args` in-place based on answers.

    Returns a dict of decisions (used by main() to know whether to run the
    example validation run at the end).
    """
    welcome_banner()

    print("────────── Phase 1: Interview ──────────")
    print()

    # ---- Q1: GTDB ----------------------------------------------------------
    print("Q1. The GTDB SSU reference database is required (used by every classifier).")
    print(f"    Release: {args.gtdb_release}   (~3 GB download + ~10 min parsing, one-time)")
    if not _ask_yn("Continue with this GTDB release?", default="y"):
        print("[launcher] Aborted by user (re-run with --gtdb-release <N> to pin a release).")
        sys.exit(0)

    # ---- Q2: filter <1000 bp ----------------------------------------------
    print()
    print("Q2. Drop reference sequences shorter than 1000 bp from the GTDB DB?")
    print("    Recommended for full-length 16S — short reference seqs add noise.")
    keep_short = not _ask_yn("Apply the 1000 bp filter?", default="y")
    args.min_ref_len = 0 if keep_short else 1000

    # The Emu / two-step NB prompts only apply to 16S. For the fungal ITS
    # marker the read-level EM classifier is EMITS (not Emu) and the NB design is
    # single-step, so we skip these prompts and default to blca + nb + emits.
    is_its = (getattr(args, "marker", "16S") == "ITS")

    if is_its:
        # ---- ITS: fixed fungal classifier set --------------------------------
        print()
        print("Marker = ITS (fungal): classifiers default to BLCA + single-step")
        print("    Naive-Bayes (nb) + EMITS read-level profiler. Emu and the 16S")
        print("    two-step NB do not apply to ITS, so those prompts are skipped.")
        use_emu = False           # Emu is 16S-only
        use_nb = True             # ITS single-step NB
        args.classifier = "blca,nb,emits"
    else:
        # ---- Q3: Emu -----------------------------------------------------------
        print()
        print("Q3. Add the Emu classifier?")
        print("    EM-based species profiling on raw reads (Curry 2022 Nat. Methods).")
        print("    Adds: ~3-5 min env install + ~5 min Emu DB build.")
        use_emu = _ask_yn("Enable Emu?", default="n")

        # ---- Q4: DADA2 NB ------------------------------------------------------
        print()
        print("Q4. Add the DADA2 Naive-Bayes classifier (nb)?")
        print("    DADA2 produces ASVs; DADA2's assignTaxonomy names them to genus")
        print("    (Wang 2007 AEM, the RDP-style 16S Naive-Bayes), then exact-match")
        print("    addSpecies recovers species. DADA2 + QIIME2 are already in the")
        print("    amplicon image, so on -profile singularity/docker there's no extra")
        print("    env install. No training: the two GTDB references build in ~1-2 min.")
        print(f"    Adds: {ANSI_RED}~1-2 min reference build, ~11 GB RAM at run time{ANSI_RST}.")
        use_nb = _ask_yn("Enable nb?", default="n")

        # ---- Build the classifier list ----------------------------------------
        branches = ["blca"]
        if use_emu: branches.append("emu")
        if use_nb:  branches.append("nb")
        args.classifier = ",".join(branches)

    # ---- Q5: validation run -----------------------------------------------
    # The bundled reference-compared validation (example_test) is the 16S ATCC
    # mock; there's no ITS equivalent wired into example_test, so only offer it
    # for 16S. ITS users validate via `-profile test_its` (see nextflow.config).
    if is_its:
        run_example = False
        print()
        print("    (Skipping the bundled 16S mock validation — marker is ITS.")
        print("     Validate the ITS path with:  nextflow run . -profile test_its)")
    else:
        print()
        print("Q5. Run the bundled 8-sample ATCC mock community to validate the install?")
        print("    Recommended — confirms every classifier you chose works end-to-end.")
        run_example = _ask_yn("Run validation at the end?", default="y")

    # ---- Make subsequent preflights unattended ----------------------------
    # The user has now confirmed every install upfront, so the GTDB/Emu/NB
    # preflights should proceed without re-asking. assume_yes=True does that.
    args.assume_yes = True

    # ---- Summary -----------------------------------------------------------
    print()
    print("────────── Your choices ──────────")
    print(f"  • GTDB release         : {args.gtdb_release}")
    print(f"  • Reference filter     : "
          f"{'drop <1000 bp' if args.min_ref_len else 'keep all'}")
    print(f"  • Classifiers          : {args.classifier}")
    print(f"  • Validation run after : {'yes' if run_example else 'no'}")
    print()
    print(f"────────── Phase 2: Installing (unattended) ──────────")
    print(f"  Sit back — total ~20-30 min on a fresh Linux machine.")
    print()

    return {
        "run_example": run_example,
        "use_emu": use_emu,
        "use_nb": use_nb,
    }


def count_terminated(trace):
    try:
        with open(trace) as fh:
            return max(0, sum(1 for _ in fh) - 1)
    except OSError:
        return 0


def run_with_progress(cmd, trace):
    try:
        from tqdm import tqdm
    except ImportError:
        print("[launcher] tqdm not installed; running without progress bar.")
        return subprocess.run(cmd).returncode
    if Path(trace).exists():
        try: Path(trace).unlink()
        except OSError: pass
    bar = tqdm(total=0, unit="task", desc="submitting", dynamic_ncols=True)
    lock = threading.Lock(); stop = threading.Event()

    def poll():
        while not stop.is_set():
            with lock:
                bar.n = min(count_terminated(trace), bar.total or count_terminated(trace))
                bar.refresh()
            time.sleep(2)

    t = threading.Thread(target=poll, daemon=True); t.start()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    submitted = 0
    for line in proc.stdout:
        if re.search(r"(Submitted process|Cached process) >", line):
            submitted += 1
            m = re.search(r"process > (\S+)", line)
            with lock:
                bar.total = submitted
                if m: bar.set_description(m.group(1)[:38])
                bar.refresh()
        tqdm.write(line.rstrip())
    proc.wait(); stop.set(); t.join(timeout=3)
    with lock:
        done = count_terminated(trace)
        bar.total = max(bar.total, done); bar.n = done; bar.refresh()
    bar.close()
    return proc.returncode


def _auto_find_env(tool, candidate_names):
    """Search for a conda env containing the given tool.

    Two-strategy lookup so this works on every conda layout we've seen:
      1. If `conda` is on PATH, parse `conda env list` and check each env's
         bin/<tool>. This finds envs no matter where conda lives (HPC scratch,
         /apps/, /srv/scratch/$USER/, etc.) and ignores the candidate name list.
      2. Fall back to a hard-coded set of common prefix dirs + the candidate
         names — useful when `conda` isn't on PATH yet but the env happens to
         sit in a standard place.

    Returns the conda-env path (e.g. /scratch/x/miniconda3/envs/hifitax_emu)
    or None if no env has `tool` on PATH.
    """
    # Strategy 1: `conda env list` — most robust because it walks whatever
    # conda already knows about, including envs in scratch, /apps/, etc.
    try:
        result = subprocess.run(
            ["conda", "env", "list"], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            for line in result.stdout.splitlines():
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split()
                env_path = parts[-1]   # last column is always the prefix path
                if os.path.isdir(env_path) and \
                   os.path.isfile(os.path.join(env_path, "bin", tool)):
                    return env_path
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        pass

    # Strategy 2: hard-coded prefix dirs × candidate names (legacy fallback).
    prefixes = [
        os.path.expanduser("~/miniconda3/envs"),
        os.path.expanduser("~/anaconda3/envs"),
        os.path.expanduser("~/.conda/envs"),
        "/opt/conda/envs",
        "/srv/scratch/" + os.environ.get("USER", "") + "/miniconda3/envs",
        "/scratch/" + os.environ.get("USER", "") + "/miniconda3/envs",
    ]
    if os.environ.get("CONDA_ENVS_PATH"):
        prefixes = os.environ["CONDA_ENVS_PATH"].split(":") + prefixes
    for name in candidate_names:
        for prefix in prefixes:
            env_path = os.path.join(prefix, name)
            if os.path.isfile(os.path.join(env_path, "bin", tool)):
                return env_path
    return None


def _install_qiime2_env(target_name="qiime2-amplicon-2024.10"):
    """Install QIIME2 amplicon distro using THEIR own always-working install
    YAML (https://data.qiime2.org/distro/...). Adds seqkit+csvtk after.

    Beats relying on a static lock file in the repo because conda-forge churns
    every few weeks and pinned builds get GC'd. QIIME2 maintain their YAML
    to actually solve against current conda-forge state.

    Returns the env path on success, None on failure.
    """
    builder = _ensure_mamba() or shutil.which("conda")
    if not builder:
        print(f"[launcher] need mamba or conda first.")
        return None

    url = ("https://data.qiime2.org/distro/amplicon/"
           "qiime2-amplicon-2024.10-py310-linux-conda.yml")
    print(f"[launcher] downloading QIIME2 spec (the only spec that actually solves)")
    print(f"[launcher]   {url}")

    import tempfile, urllib.request
    fd, tmp_yml = tempfile.mkstemp(suffix=".yml", prefix="qiime2_")
    os.close(fd)
    try:
        urllib.request.urlretrieve(url, tmp_yml)
    except Exception as exc:
        print(f"[launcher] download failed: {exc}")
        try: os.unlink(tmp_yml)
        except OSError: pass
        return None

    print(f"[launcher] creating env '{target_name}' via "
          f"{os.path.basename(builder)} (~10-15 min on Linux native)")
    rc = subprocess.run(
        [builder, "env", "create", "-n", target_name, "-f", tmp_yml, "-y"]
    ).returncode
    try: os.unlink(tmp_yml)
    except OSError: pass

    if rc != 0:
        print(f"[launcher] ✗ QIIME2 env create failed (exit {rc})")
        return None
    print(f"[launcher] ✓ QIIME2 base env created")

    # Add the two extras the QC step needs (seqkit, csvtk).
    print(f"[launcher] adding seqkit + csvtk to {target_name} (~1 min)")
    rc2 = subprocess.run(
        [builder, "install", "-n", target_name, "-c", "bioconda",
         "--freeze-installed", "-y", "seqkit", "csvtk"]
    ).returncode
    if rc2 != 0:
        print(f"[launcher] WARNING: seqkit/csvtk add failed (exit {rc2}); "
              f"QC step may need: mamba install -n {target_name} -c bioconda seqkit csvtk")
    else:
        print(f"[launcher] ✓ seqkit + csvtk added")

    return _auto_find_env("qiime", [target_name])


def _install_emu_env(target_name="hifitax_emu"):
    """Install Emu env directly from bioconda packages (no YAML to break).
    Emu has few dependencies so conflicts are rare. Returns env path / None."""
    builder = _ensure_mamba() or shutil.which("conda")
    if not builder:
        print(f"[launcher] need mamba or conda first.")
        return None

    print(f"[launcher] creating env '{target_name}' with "
          f"emu + minimap2 + samtools (~3-5 min)")
    # --override-channels: use ONLY conda-forge + bioconda. Mixing in the user's
    # `defaults` channels under strict channel priority makes the emu/samtools/
    # pysam graph unsolvable on many setups.
    rc = subprocess.run(
        [builder, "create", "-n", target_name, "--override-channels",
         "-c", "conda-forge", "-c", "bioconda", "-y",
         "python>=3.10", "emu>=3.6.2", "minimap2>=2.24", "samtools>=1.17"],
    ).returncode
    if rc != 0:
        print(f"[launcher] ✗ Emu env create failed (exit {rc})")
        return None
    print(f"[launcher] ✓ Emu env created")
    return _auto_find_env("emu", [target_name])


def _pick_env_spec_for_platform(yml_path):
    """On Linux, prefer a sibling `*.linux-lock.yml` if one exists — it's the
    fully-pinned, validated spec that won't trip conda's solver into the
    'Found conflicts!' loop that the unpinned yml hits with QIIME2 etc.
    On other platforms or when no lock exists, use the original yml.
    """
    if sys.platform.startswith("linux") and yml_path.endswith(".yml"):
        lock = yml_path[:-4] + ".linux-lock.yml"
        if os.path.isfile(lock):
            return lock
    return yml_path


def _ensure_mamba():
    """Make sure `mamba` is on PATH. If only `conda` is available, install
    mamba into the conda base env first — this is a one-time ~3-5 min cost
    that makes every future env-create ~10× faster.

    Returns the path to mamba (str) on success, or None if mamba couldn't
    be made available (in which case the caller falls back to conda).
    """
    mamba = shutil.which("mamba")
    if mamba:
        return mamba
    conda = shutil.which("conda")
    if not conda:
        return None
    print(f"[launcher] mamba not found — bootstrapping via conda (one-time, ~3-5 min)")
    print(f"[launcher] {conda} install -n base mamba -c conda-forge -y")
    rc = subprocess.run(
        [conda, "install", "-n", "base", "mamba", "-c", "conda-forge", "-y"],
    ).returncode
    if rc != 0:
        print(f"[launcher] mamba bootstrap failed (exit {rc}); falling back to conda")
        return None
    # Mamba should now be available — search PATH AND conda base's bin/.
    mamba = shutil.which("mamba")
    if not mamba:
        conda_base = os.path.dirname(os.path.dirname(conda))
        candidate = os.path.join(conda_base, "bin", "mamba")
        if os.path.isfile(candidate):
            mamba = candidate
    if mamba:
        print(f"[launcher] ✓ mamba installed at {mamba}")
    return mamba


def _create_conda_env(yml_path, label):
    """Create a conda env from a YAML spec, FORCING mamba use (bootstraps
    mamba via conda if needed). Streams progress straight to the terminal
    so the user can see it work.

    Returns True on success, False otherwise. Used by preflight_emu /
    preflight_nb to bootstrap a missing env on a fresh install.
    """
    if not os.path.isfile(yml_path):
        print(f"[launcher] no env spec found at {yml_path}; skipping auto-create")
        return False
    builder = _ensure_mamba() or shutil.which("conda")
    if not builder:
        print(f"[launcher] neither 'mamba' nor 'conda' on PATH — "
              f"can't auto-create the {label} env.")
        print(f"[launcher] install miniforge (ships mamba) first, then re-run:")
        print(f"             https://github.com/conda-forge/miniforge")
        return False
    builder_name = os.path.basename(builder)
    print(f"[launcher] ── creating {label} env from {yml_path} via "
          f"{builder_name} ──")
    if builder_name == "mamba":
        print(f"[launcher] (mamba is fast — expect ~3-10 min depending on env size)")
    else:
        print(f"[launcher] (using conda fallback — expect ~10-30 min)")
    # No capture — let the builder's progress bar paint live on the terminal.
    rc = subprocess.run(
        [builder, "env", "create", "-f", yml_path, "-y"],
    ).returncode
    if rc == 0:
        print(f"[launcher] ✓ {label} env created.")
        return True
    print(f"[launcher] ✗ {label} env creation failed (exit {rc}).")
    return False


def _dir_size(*paths):
    """Total bytes of files under the given dirs (missing dirs count as 0)."""
    total = 0
    for path in paths:
        if not path:
            continue
        for root, _, files in os.walk(path):
            for f in files:
                try:
                    total += os.path.getsize(os.path.join(root, f))
                except OSError:
                    pass
    return total


def _watch_build(watch_dirs, stop_event, poll=8):
    """Report the steps `<engine> pull` shows no progress bar for:

      - the OCI *extraction* (these dirs grow as the image is unpacked), and
      - the post-build *cleanup* of the big uncompressed temp (these dirs shrink)
        -- the silent gap before the next image that can look like a crash.

    The blob download and the SIF creation have the engine's own bars, and both
    happen without shrinking the temp, so this won't fight those bars."""
    base = _dir_size(*watch_dirs)
    last = base
    start = time.time()
    while not stop_event.wait(poll):
        sz = _dir_size(*watch_dirs)
        if sz > last + 5_000_000:            # growing = unpacking the image
            print(f"[images]   extracting + building image… {(sz - base) / 1048576:,.0f} MB "
                  f"({int(time.time() - start)}s)", flush=True)
        elif sz < last - 50_000_000:         # shrinking = cleaning up the build temp
            print(f"[images]   finalising image (cleaning up build temp)… "
                  f"({int(time.time() - start)}s)", flush=True)
        last = sz


def _pipeline_images(passthrough):
    """The container images this pipeline uses, read from the actual source --
    the modules' `container "..."` directives plus the blca_container param --
    so the list can't drift when an image is retagged. (We deliberately do NOT
    keep a hardcoded list.)"""
    imgs = []
    for nf in sorted(glob.glob(os.path.join(str(PROJDIR), "modules", "*.nf"))):
        try:
            text = open(nf).read()
        except OSError:
            continue
        for m in re.finditer(r'container\s+"([^"$]+)"', text):   # literal images only
            imgs.append(m.group(1))
    # the BLCA image comes from a param: honour --blca_container, else the config default
    blca = None
    if "--blca_container" in passthrough:
        i = passthrough.index("--blca_container")
        if i + 1 < len(passthrough):
            blca = passthrough[i + 1]
    if not blca:
        try:
            m = re.search(r"blca_container\s*=\s*'([^']+)'", (PROJDIR / "nextflow.config").read_text())
            if m:
                blca = m.group(1)
        except OSError:
            pass
    if blca:
        imgs.append(blca)
    seen, out = set(), []
    for i in imgs:
        if i not in seen:
            seen.add(i); out.append(i)
    return out


def prepull_images(args, passthrough):
    """Pull the container images BEFORE Nextflow starts, so the user sees the
    engine's own colour / real-time download progress. Nextflow otherwise pulls
    them silently mid-run (layers land in APPTAINER_CACHEDIR, not the Nextflow
    cache, so there's nothing to watch -- it just looks frozen)."""
    if args.no_image_progress:
        return
    imgs = _pipeline_images(passthrough)
    if not imgs:
        return
    if args.profile == "singularity":
        engine = shutil.which("apptainer") or shutil.which("singularity")
        if not engine:
            print("[images] apptainer/singularity not on PATH; skipping pre-pull.")
            return
        cache_dir = os.environ.get("NXF_SINGULARITY_CACHEDIR", os.path.join(os.getcwd(), "singularity_cache"))
        print(f"[images] pre-pulling {len(imgs)} image(s) into {cache_dir} "
              f"(once; re-used on later runs)")
        print("[images] NOTE: the first build converts each image to a SIF, which is slow "
              "on HPC shared filesystems -- roughly QIIME2 ~45 min, pb-16s-nf-tools ~15 min, "
              "BLCA ~5 min (~1 h total). One-time: later runs reuse the cache.")
        for img in imgs:
            # Match Nextflow's Singularity cache filename so it re-uses our file
            # instead of pulling again: strip protocol, ':' and '/' -> '-', + .img
            dest = os.path.join(cache_dir, img.replace("/", "-").replace(":", "-") + ".img")
            if os.path.exists(dest):
                print(f"[images]   already cached: {img}")
                continue
            print(f"[images]   pulling {img} …")
            print( "[images]   (sequence: download bars → 'extracting + building image' (no bar) "
                   "→ Apptainer's SIF bar → a quiet temp cleanup before the next image)")
            # best-effort progress for the silent SIF-conversion step
            watch_dirs = [cache_dir, os.environ.get("APPTAINER_TMPDIR")]
            stop = threading.Event()
            t = threading.Thread(target=_watch_build, args=(watch_dirs, stop), daemon=True)
            t.start()
            try:
                rc = subprocess.run([engine, "pull", dest, f"docker://{img}"]).returncode
            finally:
                stop.set(); t.join(timeout=3)
            if rc != 0:
                print(f"[images]   pull failed for {img}; Nextflow will retry during the run.")
                if os.path.exists(dest):
                    try: os.remove(dest)
                    except OSError: pass
    elif args.profile == "docker":
        if not shutil.which("docker"):
            return
        print(f"[images] pre-pulling {len(imgs)} image(s) with docker")
        for img in imgs:
            print(f"[images]   pulling {img} …")
            subprocess.run(["docker", "pull", "--platform", "linux/amd64", img])


def preflight(args, blca_db, blca_tax):
    if args.skip_gtdb_check:
        return 0
    pf = [sys.executable, str(PROJDIR / "bin" / "preflight_gtdb.py"),
          "--db-dir", args.gtdb_db_dir, "--blca-db", blca_db, "--blca-tax", blca_tax,
          "--timeout", str(args.gtdb_timeout), "--gtdb-release", str(args.gtdb_release)]
    if args.assume_yes: pf.append("--assume-yes")
    if args.assume_no: pf.append("--assume-no")
    if args.min_ref_len is not None: pf += ["--min-ref-len", str(args.min_ref_len)]
    return subprocess.run(pf).returncode


def _classifier_set(s):
    """Parse the --classifier value into a set: respects 'all'/'both' shorthands
    and comma-lists. Mirrors main.nf's parsing.

    Canonical user-facing name for the Naive-Bayes branch is `nb`
    (DADA2 ASVs -> DADA2 assignTaxonomy genus + exact-match addSpecies).
    Backwards-compatible aliases: `dada2_nb`, `dada2-nb`, `qiime2_nb`,
    `qiime2-nb`. Internally everything reduces to `nb`.
    """
    raw = (s or "").lower().strip()
    if raw == "all":
        return {"blca", "emu", "nb"}
    if raw == "both":
        return {"blca", "emu"}
    pieces = {x.strip() for x in raw.split(",") if x.strip()}
    nb_aliases = {"nb", "dada2_nb", "dada2-nb", "qiime2_nb", "qiime2-nb"}
    return {("nb" if p in nb_aliases else p) for p in pieces}


def _detect_profile(args, passthrough):
    """Detect which Nextflow profile the user picked, so we can route the NB
    step through the QIIME2 container instead of a conda env.
    Returns 'singularity' | 'apptainer' | 'docker' | 'standard' | 'conda' | None.

    Resolution order:
      1. args.profile   — the launcher's own --profile flag (most common path;
         the launcher captures it before passthrough is even built).
      2. passthrough    — fallback for `-profile` / `-profile=` / `--profile=`
         smuggled past the launcher (e.g. by users who only know Nextflow flags).
    """
    # 1. Launcher's own --profile (this is where `--profile singularity` lands)
    if getattr(args, "profile", None):
        return args.profile.split(",")[0].strip().lower() or None
    # 2. Passthrough fallback
    for i, tok in enumerate(passthrough):
        if tok in ("-profile", "--profile") and i + 1 < len(passthrough):
            return passthrough[i + 1].split(",")[0].strip().lower() or None
        for prefix in ("-profile=", "--profile="):
            if tok.startswith(prefix):
                return tok.split("=", 1)[1].split(",")[0].strip().lower() or None
    return None


def _emu_container_ref():
    """The pinned Emu image, read from nextflow.config (single source of truth)."""
    try:
        m = re.search(r"emu_container\s*=\s*'([^']+)'",
                      (PROJDIR / "nextflow.config").read_text())
        if m:
            return m.group(1)
    except OSError:
        pass
    return None


def preflight_emu(args, blca_db, blca_tax, emu_db_dir):
    """If the classifier needs Emu and the Emu DB is missing, build it from the
    BLCA-parsed GTDB files via bin/build_gtdb_emu_db.sh. Keeps BLCA and Emu
    anchored on the same GTDB release.

    On -profile singularity/docker the build runs INSIDE the Emu image (which
    ships emu + python3), so no conda env is created — this avoids the fragile
    emu/samtools/pysam conda solve that fails on many HPCs. The conda path is
    used only for -profile standard/conda, or when no container engine is found."""
    if "emu" not in _classifier_set(args.classifier):
        return 0
    if args.skip_emu_db_check:
        return 0
    # treat the DB as present if the dir exists and is non-empty
    if os.path.isdir(emu_db_dir) and any(os.scandir(emu_db_dir)):
        return 0
    # missing -- decide whether to build
    build = args.assume_yes
    if not build and not args.assume_no and sys.stdin.isatty():
        try:
            ans = input(f"[launcher] Emu DB not found at {emu_db_dir}. Build it from "
                        f"the GTDB BLCA DB now? (~10-20 min) [y/N] ").strip().lower()
        except EOFError:
            ans = ""
        build = ans in ("y", "yes")
    if not build:
        print(f"[launcher] Emu DB missing at {emu_db_dir}; "
              f"pass --assume-yes to auto-build, or run bin/build_gtdb_emu_db.sh manually.")
        return 2
    build_script = str(PROJDIR / "bin" / "build_gtdb_emu_db.sh")
    emu_img = _emu_container_ref()
    profile = (getattr(args, "profile", "") or "").split(",")[0].strip().lower()

    # Containerized profiles: build the Emu DB inside the Emu image (it already
    # has emu + python3), sidestepping the conda solve. Honoured unless the user
    # explicitly pointed at a conda env with --emu_env.
    if not args.emu_env and emu_img and profile in ("singularity", "apptainer", "docker"):
        os.makedirs(emu_db_dir, exist_ok=True)
        binds = []
        for p in (str(PROJDIR), os.path.dirname(os.path.abspath(blca_db)),
                  os.path.dirname(os.path.abspath(blca_tax)),
                  os.path.abspath(emu_db_dir), os.getcwd()):
            if p and os.path.isdir(p) and p not in binds:
                binds.append(p)
        if profile == "docker" and shutil.which("docker"):
            vols = []
            for d in binds:
                vols += ["-v", f"{d}:{d}"]
            cmd = (["docker", "run", "--rm", "--platform", "linux/amd64",
                    "-u", f"{os.getuid()}:{os.getgid()}", "-w", os.getcwd()]
                   + vols + [emu_img, "bash", build_script, blca_db, blca_tax, emu_db_dir])
            print("[launcher] building the Emu DB inside the Emu container (docker; no conda needed)")
            print("[launcher] " + " ".join(cmd))
            return subprocess.run(cmd).returncode
        engine = shutil.which("apptainer") or shutil.which("singularity")
        if profile in ("singularity", "apptainer") and engine:
            # Apptainer needs its cache/tmp dirs to exist before it can convert
            # the OCI image to SIF. set_apptainer_cache.sh exports these, but the
            # dir may be absent (fresh clone, or removed by a cleanup) while the
            # env var still points at it -> ensure they exist here.
            for var in ("APPTAINER_TMPDIR", "SINGULARITY_TMPDIR",
                        "APPTAINER_CACHEDIR", "SINGULARITY_CACHEDIR",
                        "NXF_SINGULARITY_CACHEDIR"):
                d = os.environ.get(var)
                if d:
                    try:
                        os.makedirs(d, exist_ok=True)
                    except OSError:
                        pass
            bind_flags = []
            for d in binds:
                bind_flags += ["--bind", d]
            cmd = ([engine, "exec"] + bind_flags + [f"docker://{emu_img}",
                    "bash", build_script, blca_db, blca_tax, emu_db_dir])
            print("[launcher] building the Emu DB inside the Emu container (no conda needed)")
            print("[launcher] " + " ".join(cmd))
            return subprocess.run(cmd).returncode
        # no container engine on PATH -> fall through to the conda path below

    # --- conda path: -profile standard / conda, or when no container engine is found ---
    cmd = ["bash", build_script, blca_db, blca_tax, emu_db_dir]
    # Prefer --emu_env if user gave it; else auto-detect any conda env that has
    # `emu` on PATH. If still nothing, auto-CREATE one from envs/emu.yml using
    # mamba (preferred, much faster) or conda. The user doesn't have to know
    # any env names — fresh install Just Works.
    env = os.environ.copy()
    emu_env = args.emu_env or _auto_find_env(
        "emu", ["hifitax_emu", "hifi_emu", "emu_nf", "emu"])
    if not emu_env:
        print(f"[launcher] ── no Emu conda env found — auto-installing one ──")
        print(f"[launcher] (one-time setup; takes ~3-5 min on a fresh Linux machine)")
        # 1) Try direct bioconda install (no YAML to break). This is the
        #    reliable path for first-time users with no bioinformatics setup.
        emu_env = _install_emu_env()
        # 2) Fallback to local YAML if the direct install failed for some
        #    weird reason (e.g. offline mirror).
        if not emu_env:
            print(f"[launcher] direct install failed — trying local envs/emu.yml as fallback")
            emu_yml = _pick_env_spec_for_platform(str(PROJDIR / "envs" / "emu.yml"))
            if _create_conda_env(emu_yml, "Emu"):
                emu_env = _auto_find_env(
                    "emu", ["hifitax_emu", "hifi_emu", "emu_nf", "emu"])
    if emu_env:
        emu_bin = os.path.join(emu_env, "bin")
        env["PATH"] = emu_bin + os.pathsep + env.get("PATH", "")
        print(f"[launcher] using emu env at {emu_env}")
    else:
        print(f"[launcher] WARNING: emu env not available; build step may fail.")
    print("[launcher] " + " ".join(cmd))
    return subprocess.run(cmd, env=env).returncode


def preflight_nb(args, blca_db, blca_tax, genus_db, species_db, passthrough=None):
    """DADA2 two-step NB taxonomy, GTDB-only: ensure BOTH the genus-level
    (assignTaxonomy) and species-assignment (addSpecies) DADA2 references exist;
    build both from the BLCA-parsed GTDB if either is missing. No training step
    (DADA2 counts k-mers / exact-matches on the fly).

    `genus_db` / `species_db` are the gzipped DADA2 FASTAs (default dir:
    <projdir>/db_nb/, names gtdb_ssu_dada2_genus.fa.gz + _species.fa.gz).
    """
    if "nb" not in _classifier_set(args.classifier):
        return 0
    if args.skip_nb_classifier_check:
        return 0
    have_genus   = os.path.isfile(genus_db)   and os.path.getsize(genus_db)   > 0
    have_species = os.path.isfile(species_db) and os.path.getsize(species_db) > 0
    if have_genus and have_species:
        return 0

    print(f"[launcher] ── DADA2 NB references missing (genus={have_genus}, species={have_species}) ──")
    print(f"[launcher] building both from the BLCA-parsed GTDB (1-2 min reformat, no training)")

    os.makedirs(os.path.dirname(genus_db),   exist_ok=True)
    os.makedirs(os.path.dirname(species_db), exist_ok=True)
    cmd = ["bash", str(PROJDIR / "bin" / "build_gtdb_dada2_db.sh"),
           blca_db, blca_tax, genus_db, species_db]
    print(f"[launcher] {' '.join(cmd)}")
    return subprocess.run(cmd, env=os.environ.copy()).returncode


def _detect_unite_version(unite_db_dir):
    """Best-effort: report the installed UNITE reference version.

    Resolution order (matches what bin/build_unite_*.sh write):
      1. <unite_db_dir>/UNITE_VERSION.txt   (written by build_unite_blca_db.sh)
      2. the date in a staged FASTA filename / header, e.g.
         'sh_general_release_dynamic_19.02.2025'
    Returns the version string, or None if nothing UNITE-shaped is installed.
    """
    vtxt = os.path.join(unite_db_dir, "UNITE_VERSION.txt")
    if os.path.isfile(vtxt):
        try:
            tag = open(vtxt).read().strip()
            if tag:
                return tag.splitlines()[0].strip()
        except OSError:
            pass
    # Fall back to parsing a release tag (or a DD.MM.YYYY date) out of a staged
    # FASTA's name or its first header line.
    pat = re.compile(r"sh_general_release[\w.]*?\d{2}\.\d{2}\.\d{4}")
    datepat = re.compile(r"\b(\d{2}\.\d{2}\.\d{4})\b")
    for name in ("unite.fasta", "unite_BLCAparsed.fasta"):
        fa = os.path.join(unite_db_dir, name)
        if not os.path.isfile(fa):
            continue
        m = pat.search(os.path.basename(fa))
        if m:
            return m.group(0)
        try:
            with open(fa) as fh:
                head = fh.readline()
        except OSError:
            head = ""
        m = pat.search(head) or datepat.search(head)
        if m:
            return m.group(0)
    return None


def preflight_unite(args, unite_db_dir, blca_db, blca_tax,
                    emits_db, singlestep_db):
    """ITS-only: report the installed UNITE reference version, tell the user the
    latest release we know about, and offer to (re)build the three UNITE-derived
    references (BLCA, single-step DADA2 NB, EMITS) from a user-supplied UNITE
    FASTA. The user already has the FASTA — we just call the build scripts.

    Build scripts (bin/):
      build_unite_blca_db.sh   <unite_fasta> <db-dir> [min-len]  -> unite_BLCAparsed.{fasta,taxonomy}
      build_unite_dada2_db.sh  <blca_fasta> <blca_tax> <out_fa_gz> -> unite_full_singlestep_ref.fa.gz
      build_emits_db.sh        <unite_fasta> <db-dir>            -> unite.fasta (minimap2 target)

    Returns 0 on success / when the user keeps the existing DB, non-zero on a
    build failure. Skipped entirely with --skip-gtdb-check (the umbrella
    'use existing DBs as-is' flag).
    """
    if args.skip_gtdb_check:
        return 0

    installed = _detect_unite_version(unite_db_dir)
    have_blca = (os.path.isfile(blca_db) and os.path.isfile(blca_tax))
    have_all = have_blca and os.path.isfile(emits_db) and os.path.isfile(singlestep_db)

    print()
    print("────────── UNITE reference (ITS marker) ──────────")
    if installed:
        print(f"[unite] installed UNITE version : {installed}")
    else:
        print(f"[unite] no UNITE reference found under {unite_db_dir}")
    print(f"[unite] latest known release     : {LATEST_UNITE}")
    print(f"[unite] download releases from   : {UNITE_REPO_URL}")
    if have_all:
        print(f"[unite] all three UNITE references present "
              f"(BLCA + single-step NB + EMITS).")

    # Decide whether to (re)build.
    if have_all:
        # Present: only rebuild on explicit opt-in.
        if args.assume_yes or args.assume_no:
            return 0
        if not sys.stdin.isatty():
            return 0
        if not _ask_yn("(Re)build the UNITE references from a UNITE FASTA?", default="n"):
            print("[unite] keeping the installed UNITE references as-is.")
            return 0
    else:
        # Missing/incomplete: must build (unless told never to).
        if args.assume_no:
            print(f"[unite] UNITE references missing/incomplete and --assume-no given; "
                  f"build them with bin/build_unite_blca_db.sh / build_unite_dada2_db.sh "
                  f"/ build_emits_db.sh, then re-run.")
            return 2
        if not (args.assume_yes or sys.stdin.isatty()):
            print(f"[unite] UNITE references missing/incomplete; download a release from "
                  f"{UNITE_REPO_URL},")
            print(f"[unite] then pass --unite-fasta <path> (or build manually) and re-run.")
            return 2

    # Resolve the source UNITE FASTA: --unite-fasta wins; else prompt.
    unite_fasta = args.unite_fasta
    if not unite_fasta and sys.stdin.isatty():
        print(f"[unite] Download a UNITE general FASTA release for Fungi from:")
        print(f"[unite]   {UNITE_REPO_URL}")
        print(f"[unite] Pick the 'Current' row (whichever release you want), extract the")
        print(f"[unite] archive, then paste the path to the sh_general_release_dynamic_*.fasta below.")
        try:
            unite_fasta = input(
                "[unite] Path to your UNITE general-release FASTA "
                "(sh_general_release_dynamic_*.fasta[.gz]): ").strip()
        except EOFError:
            unite_fasta = ""
    if not unite_fasta:
        print(f"[unite] no UNITE FASTA supplied; cannot build. Download a release from "
              f"{UNITE_REPO_URL}")
        print(f"[unite] and pass --unite-fasta <path> (or paste the path when prompted).")
        return 2
    unite_fasta = os.path.abspath(os.path.expanduser(unite_fasta))
    if not os.path.isfile(unite_fasta):
        print(f"[unite] UNITE FASTA not found: {unite_fasta}")
        return 2

    os.makedirs(unite_db_dir, exist_ok=True)
    env = os.environ.copy()
    # Tag the build with the release parsed from the filename if it looks like a
    # UNITE release, so UNITE_VERSION.txt records the actual release rather than
    # the build script's hardcoded default.
    m = re.search(r"sh_general_release[\w.]*?\d{2}\.\d{2}\.\d{4}",
                  os.path.basename(unite_fasta))
    if m:
        env["UNITE_VERSION"] = m.group(0)

    print(f"[unite] building UNITE references from {unite_fasta}")
    # 1) BLCA reference (+ BLAST index) and the taxonomy that the DADA2 single-step
    #    build consumes.
    rc = subprocess.run(
        ["bash", str(PROJDIR / "bin" / "build_unite_blca_db.sh"),
         unite_fasta, unite_db_dir, str(args.min_ref_len or 0)],
        env=env).returncode
    if rc:
        print(f"[unite] ✗ build_unite_blca_db.sh failed (exit {rc})")
        return rc
    # 2) single-step 7-rank DADA2 NB reference, from the BLCA outputs.
    rc = subprocess.run(
        ["bash", str(PROJDIR / "bin" / "build_unite_dada2_db.sh"),
         blca_db, blca_tax, singlestep_db],
        env=env).returncode
    if rc:
        print(f"[unite] ✗ build_unite_dada2_db.sh failed (exit {rc})")
        return rc
    # 3) EMITS minimap2 target (raw UNITE FASTA staged verbatim).
    rc = subprocess.run(
        ["bash", str(PROJDIR / "bin" / "build_emits_db.sh"),
         unite_fasta, unite_db_dir],
        env=env).returncode
    if rc:
        print(f"[unite] ✗ build_emits_db.sh failed (exit {rc})")
        return rc
    print(f"[unite] ✓ UNITE references built under {unite_db_dir}")
    return 0


def _setup_container_caches(args):
    """Point Apptainer/Singularity + Nextflow image caches at repo-local dirs
    (off $HOME) and make sure they exist. Idempotent: respects anything already
    exported (e.g. by set_apptainer_cache.sh) but still creates the target dir,
    so a stale env var pointing at a removed dir cannot break image conversion.
    Only acts for the singularity profile. Run this BEFORE the preflights so the
    Emu DB container build uses the same caches — users don't have to remember to
    `source set_apptainer_cache.sh`."""
    if args.profile != "singularity":
        return
    cache_dir = os.environ.get("NXF_SINGULARITY_CACHEDIR") or os.path.join(os.getcwd(), "singularity_cache")
    os.environ["NXF_SINGULARITY_CACHEDIR"] = cache_dir   # match nextflow.config
    cache_root = os.path.dirname(cache_dir)
    targets = {cache_dir}
    for var, sub in (("APPTAINER_CACHEDIR", "apptainer_cache"),
                     ("SINGULARITY_CACHEDIR", "apptainer_cache"),
                     ("APPTAINER_TMPDIR", "apptainer_tmp"),
                     ("SINGULARITY_TMPDIR", "apptainer_tmp")):
        d = os.environ.get(var) or os.path.join(cache_root, sub)
        os.environ[var] = d
        targets.add(d)
    for d in targets:
        try:
            os.makedirs(d, exist_ok=True)
        except OSError:
            pass


def nf_run(args, blca_db, blca_tax, outdir, passthrough,
           input=None, metadata=None, asv_fasta=None, entry=None):
    cmd = ["nextflow", "run", str(PROJDIR / "main.nf"),
           "-profile", args.profile,
           "--outdir", outdir, "--gtdb_db_dir", args.gtdb_db_dir,
           "--blca_db", blca_db, "--blca_tax", blca_tax]
    if entry: cmd += ["-entry", entry]
    if asv_fasta: cmd += ["--asv_fasta", asv_fasta]
    if input: cmd += ["--input", input]
    if metadata: cmd += ["--metadata", metadata]
    cmd += passthrough
    os.environ.setdefault("NXF_DISABLE_CHECK_LATEST", "true")
    os.environ.setdefault("NXF_VER", "24.10.9")
    print("[launcher] " + " ".join(cmd))

    # Keep all container caches inside the clone (not $HOME) and PRE-PULL the
    # images up front, so the user sees the engine's own colour / real-time
    # download progress instead of Nextflow pulling them silently mid-run.
    _setup_container_caches(args)
    prepull_images(args, passthrough)

    if args.progress_bar:
        # tqdm bar needs line-based logs + a piped stdout (disables Nextflow's ANSI view)
        return run_with_progress(cmd + ["-ansi-log", "false"], os.path.join(outdir, "pipeline_trace.txt"))
    # default: Nextflow's native coloured, aggregated live display (inherits the terminal)
    return subprocess.run(cmd).returncode


def export_feature_table(outdir):
    """Best-effort: export the filtered feature table to TSV (needs qiime+biom)."""
    qza = Path(outdir) / "dada2" / "dada2-ccs_table_filtered.qza"
    out = Path(outdir) / "feature_table.tsv"
    if not qza.is_file() or not shutil.which("qiime") or not shutil.which("biom"):
        return
    tmp = tempfile.mkdtemp()
    try:
        subprocess.run(["qiime", "tools", "export", "--input-path", str(qza),
                        "--output-path", tmp], check=True, capture_output=True)
        subprocess.run(["biom", "convert", "-i", str(Path(tmp) / "feature-table.biom"),
                        "-o", str(out), "--to-tsv"], check=True, capture_output=True)
    except Exception as e:
        print(f"[launcher] (feature-table export skipped: {e})")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def example_test(args, blca_db, blca_tax, passthrough):
    """Run the bundled example end-to-end and compare to the reference."""
    man = PROJDIR / "example" / "samples.tsv"
    meta = PROJDIR / "example" / "metadata.tsv"
    if not man.is_file():
        print(f"[launcher] missing example manifest {man}"); return 2
    outdir = str(PROJDIR / "example_test_results")
    print(f"[launcher] running bundled example -> {outdir}")
    rc = preflight(args, blca_db, blca_tax)
    if rc: return rc
    emu_db_dir = args.emu_db_dir or str(PROJDIR / "db_emu")
    nb_dir = args.nb_db_dir or str(PROJDIR / "db_nb")
    nb_genus   = os.path.join(nb_dir, "gtdb_ssu_dada2_genus.fa.gz")
    nb_species = os.path.join(nb_dir, "gtdb_ssu_dada2_species.fa.gz")
    rc = preflight_emu(args, blca_db, blca_tax, emu_db_dir)
    if rc: return rc
    rc = preflight_nb(args, blca_db, blca_tax, nb_genus, nb_species, passthrough=passthrough)
    if rc: return rc
    # the example test uses its OWN input/metadata; drop any of the user's so
    # they don't get appended twice (which would override the example's).
    clean = _strip(passthrough, ("--input", "--metadata", "--asv_fasta", "--outdir"))
    # forward classifier choice + Emu/NB DB paths + env paths to Nextflow
    clean = clean + ["--classifier", args.classifier,
                     "--emu_db_dir", emu_db_dir,
                     "--gtdb_dada2_genus_db", nb_genus,
                     "--gtdb_dada2_species_db", nb_species]
    if args.qiime2_env: clean += ["--qiime2_env", args.qiime2_env]
    if args.emu_env:    clean += ["--emu_env",    args.emu_env]
    rc = nf_run(args, blca_db, blca_tax, outdir, clean, input=str(man), metadata=str(meta))
    if rc:
        print("[launcher] example test pipeline failed."); return rc
    export_feature_table(outdir)
    # reference comparison covers the BLCA branch; Emu output is just produced.
    if args.classifier == "emu":
        print("[launcher] (Emu-only run; no reference comparison performed.)")
        return 0
    return subprocess.run([sys.executable, str(PROJDIR / "bin" / "compare_to_reference.py"),
                           outdir, str(PROJDIR / "example" / "reference")]).returncode


def main():
    ap = argparse.ArgumentParser(
        allow_abbrev=False,
        description="HiFiTaxa launcher. Interactive by default; for batch / HPC "
                    "scheduler jobs pass --input/--metadata plus the flags below so it "
                    "runs with no prompts. Unknown options pass through to `nextflow run` "
                    "(e.g. -resume, --min_len, --pooling_method).")
    # --- inputs ---
    ap.add_argument("--input", default=None,
                    help="samples TSV (sample-id <tab> path); required for non-interactive runs")
    ap.add_argument("--metadata", default=None, help="metadata TSV (first column 'sample_name')")
    ap.add_argument("--asv_fasta", default=None, help="taxonomy-only: classify this ASV fasta (skips QC/denoise)")
    ap.add_argument("--reads-dir", default=None,
                    help="path to a directory containing your .fastq[.gz] / .fq[.gz] reads. "
                         "If omitted, the launcher looks for a '00_Reads' folder in the current "
                         "directory or the project directory; if not found in interactive mode, "
                         "you'll be prompted for the path.")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--profile", default="standard", help="standard | conda | docker | singularity")
    # --- marker (drives primers, length window, reference DBs, classifiers) ---
    ap.add_argument("--marker", default="16S", choices=["16S", "ITS"],
                    help="amplicon marker: '16S' (GTDB/Emu, default) or 'ITS' (fungal; UNITE/EMITS, "
                         "itsxrust extraction + single-step NB). Forwarded to Nextflow as --marker.")
    # --- primers / trimming (forwarded to Nextflow; defaults come from nextflow.config) ---
    ap.add_argument("--forward_primer", default=None,
                    help=f"forward primer (16S default {DEFAULT_FWD} = 27F; "
                         f"ITS default {DEFAULT_ITS_FWD} = 1391F)")
    ap.add_argument("--reverse_primer", default=None,
                    help=f"reverse primer (16S default {DEFAULT_REV} = 1492R; "
                         f"ITS default {DEFAULT_ITS_REV} = revcomp(ITS4ngsUni))")
    ap.add_argument("--skip_primer_trim", action="store_true", help="skip cutadapt primer trimming")
    # --- GTDB database ---
    ap.add_argument("--gtdb_db_dir", default=str(PROJDIR / "db"))
    ap.add_argument("--blca_db", default=None)
    ap.add_argument("--blca_tax", default=None)
    ap.add_argument("--skip-gtdb-check", action="store_true",
                    help="use the existing DB as-is (no version check, no build)")
    ap.add_argument("--assume-yes", action="store_true",
                    help="non-interactive: build/rebuild the DB if missing or outdated")
    ap.add_argument("--assume-no", action="store_true",
                    help="non-interactive: never build; proceed with the existing DB")
    ap.add_argument("--min-ref-len", type=int, default=None, help="GTDB build: drop refs < N bp (0 = keep all)")
    ap.add_argument("--gtdb-release", default="latest",
                    help="GTDB release to target. Default 'latest' = query GTDB at run time "
                         "for the newest release. Pass an explicit release (e.g. 232) to pin.")
    ap.add_argument("--gtdb-timeout", type=float, default=15.0)
    # --- UNITE database (ITS marker only) ---
    ap.add_argument("--unite_db_dir", default=None,
                    help="directory holding the UNITE-derived ITS references (unite_BLCAparsed.{fasta,"
                         "taxonomy}, unite_full_singlestep_ref.fa.gz, unite.fasta). Default: "
                         "<projdir>/db_unite. Only used when --marker ITS.")
    ap.add_argument("--unite-fasta", default=None,
                    help="path to your UNITE general-release FASTA "
                         "(sh_general_release_dynamic_*.fasta[.gz]), downloaded from "
                         "https://unite.ut.ee/repository.php; used to (re)build the ITS "
                         "references via bin/build_unite_blca_db.sh / build_unite_dada2_db.sh / "
                         "build_emits_db.sh. Only used when --marker ITS.")
    # --- run control ---
    ap.add_argument("--test", action="store_true", help="run the bundled example + compare to reference, then exit")
    ap.add_argument("--skip-test", action="store_true", help="never offer/run the first-run example test")
    ap.add_argument("--progress-bar", action="store_true",
                    help="show a tqdm progress bar instead of Nextflow's native coloured display")
    ap.add_argument("--no-image-progress", action="store_true",
                    help="don't pre-pull container images (Nextflow pulls them silently during the run)")
    # --- classifier choice ---
    ap.add_argument("--classifier", default="blca",
                    help="taxonomic classifier(s): single name (blca|emu|dada2_nb), a comma-list "
                         "(e.g. 'blca,emu,nb'), or shorthand 'all' / 'both' (== blca,emu). "
                         "Default: blca.")
    ap.add_argument("--emu_db_dir", default=None,
                    help="path to a prebuilt Emu DB dir; auto-built from the GTDB BLCA DB if missing "
                         "(default: <projdir>/db_emu)")
    ap.add_argument("--emu_env", default=None,
                    help="path to a conda env (or any dir with bin/emu) used when the launcher needs "
                         "to build the Emu DB; prepended to PATH for the build script. Also forwarded "
                         "to Nextflow for -profile standard.")
    ap.add_argument("--skip-emu-db-check", action="store_true",
                    help="use the existing Emu DB as-is; never build")
    ap.add_argument("--nb_db_dir", default=None,
                    help="directory holding the DADA2 NB references (genus + species gzipped "
                         "FASTAs: gtdb_ssu_dada2_genus.fa.gz, gtdb_ssu_dada2_species.fa.gz); both "
                         "auto-built from the GTDB BLCA DB if missing (default: <projdir>/db_nb). "
                         "DADA2 two-step genus assignTaxonomy + exact-match addSpecies (GTDB-only).")
    ap.add_argument("--qiime2_env", default=None,
                    help="path to a conda env with QIIME2 amplicon installed; used by the launcher "
                         "for the QIIME2 export/tabulate steps in the NB branch (which need `qiime` "
                         "on PATH). Also forwarded to Nextflow for -profile standard.")
    ap.add_argument("--skip-nb-classifier-check", action="store_true",
                    help="use the existing NB classifier as-is; never train")
    args, passthrough = ap.parse_known_args()

    # Did the user set these explicitly on the CLI? (Used so the interactive
    # marker prompt and the ITS classifier default don't clobber explicit flags.)
    marker_explicit = any(a == "--marker" or a.startswith("--marker=") for a in sys.argv[1:])
    classifier_explicit = any(
        a == "--classifier" or a.startswith("--classifier=") for a in sys.argv[1:])

    # ── Marker selection (interactive) ─────────────────────────────────────
    # Ask 16S vs ITS BEFORE resolving reference DBs / primers / classifiers,
    # since every one of those switches on the marker. Only prompt when we're
    # interactive and the user didn't pin --marker on the CLI; otherwise honour
    # whatever --marker holds (default '16S'). Skipped for taxonomy-only FASTA
    # runs handled the same as a normal marker default.
    _interactive_early = sys.stdin.isatty() and not (args.assume_yes or args.assume_no)
    if _interactive_early and not marker_explicit:
        print()
        print("Which amplicon marker are you classifying?")
        print("  16S = bacterial/archaeal full-length 16S rRNA (GTDB, Emu)")
        print("  ITS = fungal ITS (UNITE, EMITS; itsxrust extraction + single-step NB)")
        try:
            ans = input("    Marker [16S/ITS] (default 16S): ").strip().lower()
        except EOFError:
            ans = ""
        args.marker = "ITS" if ans in ("its", "fungal", "fungi") else "16S"
        print(f"[launcher] marker = {args.marker}")

    is_its = (args.marker == "ITS")
    # For ITS, default the classifier set to BLCA + single-step NB + EMITS and
    # skip the Emu / two-step NB prompts (Emu is 16S-only; ITS NB is single-step).
    # An explicit --classifier always wins.
    if is_its and not classifier_explicit:
        args.classifier = "blca,nb,emits"

    # marker-aware reference resolution. For ITS, blca_db/blca_tax point at the
    # UNITE references (matching the nextflow.config marker switch); 16S keeps GTDB.
    unite_db_dir = args.unite_db_dir or str(PROJDIR / "db_unite")
    if is_its:
        blca_db = args.blca_db or os.path.join(unite_db_dir, "unite_BLCAparsed.fasta")
        blca_tax = args.blca_tax or os.path.join(unite_db_dir, "unite_BLCAparsed.taxonomy")
    else:
        blca_db = args.blca_db or os.path.join(args.gtdb_db_dir, "gtdb_ssu_BLCAparsed.fasta")
        blca_tax = args.blca_tax or os.path.join(args.gtdb_db_dir, "gtdb_ssu_BLCAparsed.taxonomy")
    emits_db = os.path.join(unite_db_dir, "unite.fasta")
    unite_singlestep_db = os.path.join(unite_db_dir, "unite_full_singlestep_ref.fa.gz")
    emu_db_dir = args.emu_db_dir or str(PROJDIR / "db_emu")
    nb_dir = args.nb_db_dir or str(PROJDIR / "db_nb")
    nb_genus   = os.path.join(nb_dir, "gtdb_ssu_dada2_genus.fa.gz")
    nb_species = os.path.join(nb_dir, "gtdb_ssu_dada2_species.fa.gz")

    # Set up repo-local container caches up front (before any preflight that
    # might exec a container, e.g. the Emu DB build), so users don't have to
    # source set_apptainer_cache.sh first.
    _setup_container_caches(args)

    # explicit test mode
    if args.test:
        return example_test(args, blca_db, blca_tax, passthrough)

    interactive = sys.stdin.isatty() and not (args.assume_yes or args.assume_no)
    # batch jobs (no TTY) with no DB decision -> build if needed, so the preflight
    # never blocks waiting for an answer
    if not interactive and not (args.assume_yes or args.assume_no or args.skip_gtdb_check):
        args.assume_yes = True

    input_tsv, metadata_tsv = args.input, args.metadata
    fwd, rev = args.forward_primer, args.reverse_primer

    # Welcome wizard: first-run + interactive only. Asks ALL setup questions
    # upfront (classifiers, GTDB filter, validation) so the install phase
    # below runs fully unattended — no surprise prompts halfway through a
    # 30-min env build. The validation example run is queued for AFTER the
    # installs complete (see end of main()).
    wizard_decisions = None
    if not MARKER.exists():
        if interactive:
            wizard_decisions = welcome_wizard(args)
        else:
            # Non-interactive fresh install (e.g. PBS job, CI): silently
            # mark the welcome as done and proceed with whatever flags the
            # user passed.
            welcome_banner()
            MARKER.write_text("done\n")

    # resolve inputs: explicit --input/--asv_fasta, else interactive sample-sheet build
    if not (args.asv_fasta or input_tsv):
        if interactive:
            input_tsv, metadata_tsv, fwd, rev = interactive_setup(
                reads_dir=args.reads_dir, marker=args.marker)
        else:
            print("[launcher] no --input or --asv_fasta given.")
            print("[launcher] for a batch job pass, e.g.:")
            print("[launcher]   --input samples.tsv --metadata metadata.tsv \\")
            print("[launcher]   --forward_primer <FWD> --reverse_primer <REV> --profile singularity --assume-no")
            return 0

    # forward primer / trim choices to Nextflow (unset -> nextflow.config defaults)
    extra = list(passthrough)
    extra += ["--marker", args.marker]
    if fwd: extra += ["--forward_primer", fwd]
    if rev: extra += ["--reverse_primer", rev]
    if args.skip_primer_trim: extra += ["--skip_primer_trim", "true"]
    # classifier choice; forward the env paths the launcher consumed
    # (--qiime2_env, --emu_env) so -profile standard sees them too.
    extra += ["--classifier", args.classifier]
    if is_its:
        # ITS reference paths (UNITE): EMITS minimap2 target + single-step NB ref.
        extra += ["--unite_db_dir", unite_db_dir,
                  "--emits_db", emits_db,
                  "--unite_dada2_singlestep_db", unite_singlestep_db]
    else:
        # 16S reference paths (GTDB): Emu DB + two-step DADA2 NB references.
        extra += ["--emu_db_dir", emu_db_dir,
                  "--gtdb_dada2_genus_db", nb_genus,
                  "--gtdb_dada2_species_db", nb_species]
    if args.qiime2_env: extra += ["--qiime2_env", args.qiime2_env]
    if args.emu_env:    extra += ["--emu_env",    args.emu_env]

    if is_its:
        # ITS: report the installed UNITE version, offer a (re)build from a
        # user-supplied UNITE FASTA, then ensure the three UNITE references
        # exist. Emu / two-step GTDB NB preflights do not apply.
        rc = preflight_unite(args, unite_db_dir, blca_db, blca_tax,
                             emits_db, unite_singlestep_db)
        if rc: return rc
    else:
        # 16S: GTDB version check + optional Emu / two-step NB reference builds.
        rc = preflight(args, blca_db, blca_tax)
        if rc: return rc
        rc = preflight_emu(args, blca_db, blca_tax, emu_db_dir)
        if rc: return rc
        rc = preflight_nb(args, blca_db, blca_tax, nb_genus, nb_species, passthrough=passthrough)
        if rc: return rc

    # ── Phase 3: validation run (only on first install, only if user opted in) ──
    # The wizard collected this decision upfront so we can run unattended.
    # We run validation AFTER every install succeeds and BEFORE the user's real
    # data, so any setup problem surfaces against the known-good 8-sample ATCC
    # mock instead of disguising itself as a problem with the user's data.
    if (wizard_decisions and wizard_decisions["run_example"]
            and not args.skip_test and not MARKER.exists()):
        print()
        print("────────── Phase 3: Validation run (bundled 8-sample ATCC mock) ──────────")
        rc = example_test(args, blca_db, blca_tax, passthrough)
        if rc:
            print(f"[launcher] validation failed (exit {rc}). Not marking first-run "
                  f"complete; fix the issue and re-run.")
            return rc
        MARKER.write_text("done\n")
        print()
        print("[launcher] ✓ Install validated. Proceeding to your data run.")
        print()
    elif not MARKER.exists():
        # Fresh install but user declined / skipped validation — still mark
        # first-run done so we don't ask again.
        MARKER.write_text("done\n")

    entry = "taxonomy_only" if args.asv_fasta else None
    # taxonomy_only supports all three classifiers (BLCA, Emu, NB) on a FASTA.
    return nf_run(args, blca_db, blca_tax, args.outdir, extra,
                  input=input_tsv, metadata=metadata_tsv,
                  asv_fasta=args.asv_fasta, entry=entry)


def _sample_id(fname):
    """Derive a clean sample id from a fastq filename."""
    b = fname
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if b.endswith(ext):
            b = b[:-len(ext)]
            break
    b = b.replace(".hifi_reads", "")
    m = re.search(r"For_(bc\d+)--16S_Rev_(bc\d+)", fname)   # PacBio demux pattern
    if m:
        return f"{m.group(1)}_{m.group(2)}"
    return re.sub(r"[^A-Za-z0-9._-]", "_", b)


def _ask(q):
    try:
        return input(q).strip().lower() in ("y", "yes")
    except EOFError:
        return False


def _find_reads_dir(explicit=None):
    """Locate the reads folder.

    Search order:
      1. `explicit` if given (e.g. --reads-dir) — error if it doesn't exist.
      2. `./00_Reads/` (current dir).
      3. `<projdir>/00_Reads/` (repo dir).
    Returns the absolute path or None if not found.
    """
    if explicit:
        d = os.path.abspath(os.path.expanduser(explicit))
        return d if os.path.isdir(d) else None
    for base in (os.getcwd(), str(PROJDIR)):
        d = os.path.join(base, READS_DIR)
        if os.path.isdir(d):
            return d
    return None


def _find_reads(d):
    """Read files in a folder: gzipped or plain. (QC_fastq gzip-normalises both,
    so plain .fastq works end-to-end.)"""
    fqs = []
    for pat in ("*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"):
        fqs += glob.glob(os.path.join(d, pat))
    return sorted(set(fqs))


def _metadata_blanks(path):
    """Return the sample ids whose 'condition' (2nd column) is still empty."""
    blanks = []
    with open(path) as fh:
        next(fh, None)  # skip header
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            cond = parts[1].strip() if len(parts) > 1 else ""
            if not cond:
                blanks.append(parts[0] if parts else "?")
    return blanks


def interactive_setup(reads_dir=None, marker="16S"):
    """Generate samples.tsv + metadata.tsv from the reads folder, have the
    user fill the metadata 'condition' column (left blank on purpose), confirm
    it's filled, then confirm primers. Returns (samples, metadata, fwd, rev).

    `marker` ('16S' or 'ITS') selects which primer defaults are offered: 16S
    27F/1492R, or the fungal ITS 1391F / revcomp(ITS4ngsUni) pair.

    If `reads_dir` is given (i.e. user passed --reads-dir), use that.
    Otherwise look for ./00_Reads and <projdir>/00_Reads; if neither exists,
    in interactive mode ask the user for the path."""
    d = _find_reads_dir(explicit=reads_dir)
    if not d and reads_dir:
        # explicit path was given but doesn't exist
        print(f"\n[setup] --reads-dir not found: "
              f"{os.path.abspath(os.path.expanduser(reads_dir))}")
        sys.exit(2)
    if not d and sys.stdin.isatty():
        # No 00_Reads folder anywhere — let the user point us at one.
        print(f"\n[setup] No '{READS_DIR}' folder found in the current directory "
              f"or the project directory.")
        try:
            ans = input(f"          Enter the path to your reads directory "
                        f"(or press Enter to create ./{READS_DIR} and re-run): "
                        ).strip()
        except EOFError:
            ans = ""
        if not ans:
            print(f"          Create the folder and add reads, then re-run:")
            print(f"            mkdir {READS_DIR} && cp /path/to/*.fastq.gz {READS_DIR}/")
            sys.exit(2)
        d = os.path.abspath(os.path.expanduser(ans))
        if not os.path.isdir(d):
            print(f"          Not found: {d}")
            sys.exit(2)
    if not d:
        print(f"\n[setup] No '{READS_DIR}' folder found here and no --reads-dir given. "
              f"Create the folder and add reads, then re-run:")
        print(f"          mkdir {READS_DIR} && cp /path/to/*.fastq.gz {READS_DIR}/")
        sys.exit(2)
    fqs = _find_reads(d)
    if not fqs:
        print(f"\n[setup] No .fastq[.gz] / .fq[.gz] files found in: {d}")
        sys.exit(2)

    print(f"\n[setup] Found {len(fqs)} read file(s) in {d}")
    samples = os.path.abspath("samples.tsv")
    metadata = os.path.abspath("metadata.tsv")

    # If either file already exists, ask before clobbering. Default = keep,
    # because the user has often already curated the metadata 'condition'
    # column and would lose it on overwrite.
    pre_existing = [p for p in (samples, metadata) if os.path.isfile(p)]
    write_files = True
    if pre_existing and sys.stdin.isatty():
        for p in pre_existing:
            print(f"[setup] Found existing  {p}")
        try:
            ans = input(
                "[setup] Overwrite these with fresh ones from the reads folder? "
                "[y/N] ").strip().lower()
        except EOFError:
            ans = ""
        write_files = ans in ("y", "yes")
        if not write_files:
            print(f"[setup] Keeping existing samples.tsv + metadata.tsv as-is.")

    if write_files:
        with open(samples, "w") as s, open(metadata, "w") as m:
            s.write("sample-id\tabsolute-filepath\n")
            m.write("sample_name\tcondition\n")
            for fq in fqs:
                sid = _sample_id(os.path.basename(fq))
                s.write(f"{sid}\t{os.path.abspath(fq)}\n")
                m.write(f"{sid}\t\n")    # condition left blank for the user to fill
        print(f"[setup] Wrote sample sheet : {samples}")
        print(f"[setup] Wrote metadata     : {metadata}  (sample_name filled in)")
        print(f"\n[setup] >>> Now open {metadata} and fill in the 'condition' column —")
        print(f"[setup]     one group label per sample — then save the file. <<<\n")

    while True:
        if not _ask("[setup] Have you filled in the condition column and saved metadata.tsv? [y/N] "):
            print(f"[setup] Please edit {metadata} first, then answer y.")
            continue
        blanks = _metadata_blanks(metadata)
        if not blanks:
            break
        print(f"[setup] These samples still have an empty condition: {', '.join(blanks)}")
        print(f"[setup] Fill them in, save, then answer y.")
    print("[setup] metadata.tsv looks complete.")

    # Offer the marker-appropriate primer defaults.
    if marker == "ITS":
        def_fwd, def_rev = DEFAULT_ITS_FWD, DEFAULT_ITS_REV
        prompt = (f"[setup] Are your primers the fungal ITS 1391F ({def_fwd}) / "
                  f"ITS4ngsUni-revcomp ({def_rev})? [y/N] ")
    else:
        def_fwd, def_rev = DEFAULT_FWD, DEFAULT_REV
        prompt = (f"[setup] Are your primers the full-length 16S 27F ({def_fwd}) / "
                  f"1492R ({def_rev})? [y/N] ")
    if _ask(prompt):
        fwd, rev = def_fwd, def_rev
    else:
        fwd = input("[setup] Forward primer sequence (5'->3'): ").strip()
        rev = input("[setup] Reverse primer sequence (5'->3'): ").strip()
    return samples, metadata, fwd, rev


def _strip(args_list, flags):
    """Drop the given '--flag value' pairs from a passthrough arg list."""
    out, skip = [], False
    for a in args_list:
        if skip:
            skip = False
            continue
        if a in flags:
            skip = True
            continue
        out.append(a)
    return out


if __name__ == "__main__":
    sys.exit(main())
