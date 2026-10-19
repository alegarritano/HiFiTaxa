# Third party notices

The MIT licence in `LICENSE` covers the original code of HiFiTaxa. The components
listed here were adapted from other projects and keep their own licences, reproduced
below as published by their authors.

## PacBio HiFi-16S-workflow (BSD 3-Clause Clear)

The quality filtering, primer removal and DADA2 denoising modules
(`modules/qc.nf`, `modules/dada2.nf`, `scripts/run_dada_2023.2.R`,
`scripts/learnError.R`, and the related process definitions in `main.nf`) were
adapted from the PacBio HiFi-16S-workflow
(https://github.com/PacificBiosciences/HiFi-16S-workflow), whose licence reads:

```
Copyright (c) 2020, Pacific Biosciences of California, Inc.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted (subject to the limitations in the
disclaimer below) provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

 * Neither the name of Pacific Biosciences nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.
```

## BLCA (GNU licence, as stated by its authors)

`bin/blca_main.py` is the BLCA classifier script `2.blca_main.py` of Gao, Lin,
Revanna and Dong (https://github.com/qunfengdong/BLCA; Gao X, Lin H, Revanna K,
Dong Q. A Bayesian taxonomic classification method for 16S rRNA gene sequences
with improved species-level accuracy. BMC Bioinformatics 18, 247 (2017),
doi:10.1186/s12859-017-1670-4). The upstream repository carries no licence file;
its README states the licence as "GNU". The copyright remains with the BLCA
authors. HiFiTaxa's copy differs from the upstream script in two places only:

1. `run_blastn` takes the identity threshold (`--iset`) and passes it to BLAST
   as `-perc_identity`, where the upstream script fixes `-perc_identity` at 90
   whatever `--iset` is requested (commit 49aceef, 17 July 2026).
2. A `--seed` option (default 42) seeds Python's random number generator before
   the bootstrap resampling, so that repeated runs give the same confidence
   values for the same input and chunking (commit b768817).

Everything else, including the BLAST, Clustal Omega and bootstrap logic, is the
upstream code.
