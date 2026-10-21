# Third-party data and software

## Research data

Gene counts and metadata: Paculova H and colleagues, **IKAROS Gene Regulatory Network Reveal ERG as a Vulnerability in B-cell Acute Lymphoblastic Leukemia**, bioRxiv version 1, 2025-12-29, DOI [10.64898/2025.12.29.696402](https://doi.org/10.64898/2025.12.29.696402); GEO [GSE310958](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310958), including reused samples from [GSE90659](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90659). The preprint is offered under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The site credits the authors, identifies transformations and does not imply their endorsement. Published count values are retained; computed selection and normalization are the demo's reanalysis.

## Gene annotation

GENCODE v36 annotation via the UCSC Genome Browser `hg38.wgEncodeGencodeAttrsV36` table. Credit the GENCODE consortium and UCSC Genome Browser team. [GENCODE data access](https://www.gencodegenes.org/pages/data_access.html) states that project data are open access; [UCSC conditions](https://genome.ucsc.edu/conditions.html) require appropriate contributor credit, and its [data licensing explanation](https://genome.ucsc.edu/license/) allows use of raw tables while retaining source restrictions. Only gene IDs/names/types are redistributed here. No browser binaries, figures or proprietary clinical annotation tracks are used.

## ssvQC-derived Pearson primitive

Repository: [FrietzeLabUVM/ssvQC](https://github.com/FrietzeLabUVM/ssvQC), commit `90a38bf72a97db43cdcdbb4e68bea9a4ceb389c5`, `R/ssvQC_correlation.R`.

Package DESCRIPTION declares `MIT + file LICENSE`. Its LICENSE identifies:

Copyright (c) 2021 Joseph R Boyd.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

The demo adapts the source-defined Pearson correlation primitive to gene counts; it does not ship or execute the entire R package. Cite Boyd J, Rodriguez P, Schjerven H, Frietze S. **ssvQC: an integrated CUT&RUN quality control workflow for histone modifications and transcription factors.** BMC Research Notes 14, 366 (2021), [doi:10.1186/s13104-021-05781-8](https://doi.org/10.1186/s13104-021-05781-8).
