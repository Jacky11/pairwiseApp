# pairwiseApp — differential analysis of label-free proteomics data between two conditions

Shiny application (and command-line script) to compare two biological conditions from a
MaxQuant `proteinGroups.txt` file or a Spectronaut protein-group report (tab-, semicolon-
or comma-separated). It returns an Excel workbook
with all the result tables and a PowerPoint report with the figures. Optionally, it runs
SAINTexpress (AP-MS experiments) and analyses of the overlap between samples and of the
robustness of the sample clustering (Jaccard and overlap indices, PERMANOVA, betadisper,
pvclust AU/BP p-values).

## Analysis pipeline

1. **Quantified values** – a quantified value is a non-zero intensity measured with at least *u* unique peptides (default 1); intensities (Raw or LFQ) are log2-transformed.
2. Proteins flagged `+` in the *Reverse*, *Potential contaminant* and *Only identified by site* columns are removed.
3. **Filter 1 (number of quantified values)** – for each condition *k* with *n_k* samples, a protein is *sufficiently quantified* if it has at least max(*m*, ⌈*p*·*n_k*⌉) quantified values (`min_observed` *m*, default 2; `min_percent_observed` *p*, default 0%). A protein is **tested** if it is sufficiently quantified in at least one condition (`quantified_in: one`, default) or in each condition (`quantified_in: each`), and has more than `absent_max` quantified values in each condition.
4. **Normalisation** – within each condition, the samples are centred on the mean of their medians (medians computed on the proteins kept by the filters that have no missing value in the condition).
5. **Filter 2 (present/absent proteins)** – a protein with at most `absent_max` quantified values (default 0) in one condition and at least `present_min` values (default: the number required in filter 1) in the other one is **present/absent**: it is not tested and is reported as more abundant in the condition where it is present (with its average iBAQ). All other proteins are **removed**; the sheet `without_enough_values` gives the reason for each of them, and all the result tables give the number of quantified values of each protein in each condition.
6. **Imputation** of the missing values of the tested proteins with `impute.mle` (or `impute.slsa`) of the R package [imp4p](https://CRAN.R-project.org/package=imp4p).
7. **Test** – moderated t-test (limma, or paired limma); the p-values of the proteins whose |log2 fold-change| exceeds the threshold (default log2(2) = 1) are adjusted with the adaptive Benjamini–Hochberg procedure, the proportion of non-differential proteins being estimated by the *pounds* (or *slim*) method of [cp4p](https://CRAN.R-project.org/package=cp4p). A protein is significant if its adjusted p-value is below the FDR level (default 1%).
   *Option* `conditional_pvalues`: before the adjustment, the p-value of each selected protein is made conditional on the selection |log2FC| > threshold, by dividing it by the p-value the protein would have with |log2FC| equal to the threshold (given its moderated standard error). Under the null hypothesis, these conditional p-values are uniformly distributed among the selected proteins. The conditional p-values are reported in the `ConditionalPvalues` sheet.

The PowerPoint report also contains Venn/UpSet plots, missing-value counts, intensity
distributions, correlation- and Jaccard-based sample clustering, PCA, volcano plots and
heatmaps of the proteins of interest.

## Repository content

```
ui.R, server.R, global.R     Shiny application
load_app.R                   loads the packages and the functions of R/
R/pipeline.R                 run_pairwise(): the whole analysis
R/read_data.R                reading of MaxQuant / Spectronaut files
R/params.R                   analysis parameters (defaults, YAML files, checks)
R/overlap_clustering.R       sample overlap indices and pvclust (optional)
R/saint.R                    SAINTexpress step (optional)
R/limma-test.R, R/adjustp.R, R/venn.R, R/utils.R
                             helper functions
scripts/run_analysis.R       runs an analysis without the interface
scripts/install_packages.R   installs the R packages
params/example.yml           parameter file with all parameters and their default values
bin/                         place for the SAINTexpress executable (optional)
```

## Installation

R ≥ 4.3 is required. Install the packages with

```r
source("scripts/install_packages.R")
```

## Usage

### Shiny application

```r
shiny::runApp("path/to/pairwiseApp")
```

1. Tick *Spectronaut input file* if needed, choose the input file and the kind of intensities (Raw or LFQ).
2. Give a name to each condition and select its samples.
3. Set the parameters (filters, normalisation, imputation, test, fold-change and FDR thresholds) and the optional analyses.
4. Click **Run the analysis**, then **Download Analysis Files**.

The zip file contains:

| File | Content |
|------|---------|
| `Final_<cond2>_vs_<cond1>.xlsx` | result tables: proteins without enough values, present/absent proteins, differentially abundant proteins, non-significant proteins, imputed data, correlations between samples, parameters |
| `Differential_Analysis_<cond2>_vs_<cond1>.pptx` | report with the description of the pipeline and all the figures |
| `parameters.yml` | all the parameters of the analysis (can be used as input of `scripts/run_analysis.R`) |
| `sessionInfo.txt` | R and package versions, size and md5 of the input file |
| `<prefix>_*` | sample overlap and pvclust results (if selected) |
| `saint/` | SAINTexpress input and output files (if selected) |

### Without the interface

```bash
# list the sample columns of a file
Rscript scripts/run_analysis.R --list-samples data/proteinGroups.txt Raw

# run an analysis described by a parameter file
Rscript scripts/run_analysis.R params/example.yml results/my_analysis
```

or from R:

```r
source("load_app.R")
p <- read_params("params/example.yml")
res <- run_pairwise(p, out_dir = "results/my_analysis")
res$counts
```

All the parameters are described in `params/example.yml`. Relative paths in a parameter
file are relative to the location of this file.

### Reproducibility

The MLE imputation, the permutation tests (seed 123) and pvclust (seed 123) are
deterministic; the SLSA imputation uses the `seed` parameter. The pvclust AU/BP values
also depend on the number of parallel workers (`n_cores`, default 5), which must be
kept identical to obtain the same values.

## SAINTexpress (optional)

SAINTexpress (Teo *et al.*, J Proteomics 2014) is not distributed with this repository.
Download it from <https://saint-apms.sourceforge.net> and put `SAINTexpress-int.exe`
(Windows) or `SAINTexpress-int` (Linux/macOS, made executable) in `bin/`, or give its path
with the `saint_exe` parameter. Condition 1 must contain the control samples and
condition 2 the bait samples; the name of condition 2 is used as bait name.

* MaxQuant files: preys are identified by their majority protein IDs; protein lengths are read in the `Sequence length` column.
* Spectronaut files: preys are identified by the first accession of their protein group; the
  protein lengths required by the prey file are estimated from `PG.MolecularWeight` (mean
  residue mass of 110 Da).

## License

GPL-3 (see [LICENSE](LICENSE)).
