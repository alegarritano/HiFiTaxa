# Marker-aware HiFiTaxa (16S and fungal ITS)

HiFiTaxa runs on two markers from one pipeline, selected with `--marker`:

| | `16S` (default) | `ITS` |
|---|---|---|
| reference | GTDB SSU (`db/`) | UNITE (`db_unite/`) |
| primers (cutadapt) | 27F `AGRGTTYGATYMTGGCTCAG` / 1492R `AAGTCGTAACAAGGTARCY` | **1391F `GTACACACCGCCCGTC` / ITS4ngsUni-rc `GCATATHANTAAGSGSAGG`** |
| read prep | cutadapt → DADA2 | cutadapt → **itsxrust ITS extraction** → DADA2 |
| DADA2 window | 1000–1600 bp | 200–1300 bp |
| Naive Bayes | two-step (genus + addSpecies) | **single-step** (7-rank `assignTaxonomy`) |
| read-level EM | Emu | **EMITS** (fungal Emu analogue) |

The 16S behaviour is unchanged. See `docs/pipeline_scheme.svg` for the flow.

## Why the ITS choices differ
- **itsxrust before DADA2**: ITS amplicons carry variable SSU/LSU flanks; itsxrust (HMM-anchored)
  extracts the ITS region robustly and primer-agnostically, so BLCA's coverage filter and DADA2
  both see clean ITS. itsxrust does NOT denoise — DADA2 still runs after it.
- **single-step NB**: the 16S two-step adds species only on an exact match; ITS varies within a
  genus so addSpecies fires for ~0.1 % of queries. Single-step `assignTaxonomy` recovers ~40 %.
- **EMITS** is competitive on clean HiFi ITS (the read-level EM ≈ a good best-hit mapper).

## Running ITS
The default-primer values resolve from `params.marker` in a trailing `params{}` block (after
`profiles{}`), so the sanctioned ITS entry points are a profile or the launcher — a bare
`--marker ITS` on the CLI does not re-trigger the primer/length/DB switches.

```bash
# bundled ITS smoke test (2,000-read fungal mock fixture)
nextflow run main.nf -profile test_its,standard

# interactive (asks marker, confirms primers, checks UNITE version, offers metadata)
python bin/run_pipeline.py --profile standard

# non-interactive ITS run
python bin/run_pipeline.py --profile standard --marker ITS \
  --input samples.tsv --metadata metadata.tsv --classifier blca,nb,emits
```

## UNITE database
The launcher checks the installed UNITE version (`db_unite/UNITE_VERSION.txt`) against the latest
known release and offers to (re)build from a user-supplied UNITE general-release FASTA:

```bash
bin/build_unite_blca_db.sh   <unite.fasta> db_unite           # BLCA reference + BLAST index
bin/build_unite_dada2_db.sh  db_unite/unite_BLCAparsed.fasta db_unite/unite_BLCAparsed.taxonomy \
                             db_unite/unite_full_singlestep_ref.fa.gz   # single-step NB ref
bin/build_emits_db.sh        <unite.fasta> db_unite           # stages unite.fasta for minimap2/EMITS
```

`itsxrust` and `EMITS` build from source (`github.com/ayobi/itsxrust`, `github.com/ayobi/emits`);
see `envs/itsxrust.yml` and `envs/emits.yml`.

## Verified
Smoke-tested on real reads: 16S example reads classify against GTDB (BLAST, 99.8–99.9 % identity);
the 2,000-read ITS fixture runs itsxrust → minimap2 → EMITS (recovers the mock's *Rhodocollybia*,
*Cortinarius*, *Russula*, *Inocybe* …) and the single-step NB script classifies to fungal lineages.
