## Robustness and applicability of transcription factor and pathway analysis tools on single-cell RNA-seq data

### Abstract
**Background**

Many functional analysis tools have been developed to extract functional and mechanistic insight from bulk transcriptome data. With the advent of single-cell RNA sequencing (scRNA-seq), it is in principle possible to do such an analysis for single cells. However, scRNA-seq data has characteristics such as drop-out events and low library sizes. It is thus not clear if functional analysis tools established for bulk sequencing can be applied to scRNA-seq in a meaningful way.

**Results**

To address this question, we performed benchmark studies on simulated and real scRNA-seq data. We included the bulk-RNA tools PROGENy, GO enrichment and DoRothEA that estimate pathway and transcription factor (TF) activities, respectively, and compared them against the tools SCENIC/AUCell and metaVIPER, designed for scRNA-seq. For the in silico study, we simulated single cells from TF/pathway perturbation bulk RNA-seq experiments. We complemented the simulated data with real scRNA-seq data upon CRISPR-mediated knock-out. Our benchmarks on both simulated and real data revealed comparable performance to the original bulk data. Additionally, we showed that the TF and pathway activities preserve cell-type-specific variability by analyzing a mixture sample sequenced with 13 scRNA-seq protocols. We also provide the benchmark data for further use by the community.

**Conclusions**

Our analyses suggest that bulk-based functional analysis tools that use manually curated footprint gene sets can be applied to scRNA-seq data, partially outperforming dedicated single-cell tools. Furthermore, we found that the performance of functional analysis tools is more sensitive to the gene sets than to the statistic used.


***

### Availabilty of data
The datasets supporting the conclusions of this publication are available at Zenodo: [![DOI](https://doi.org/10.5281/zenodo.3564179.svg)](https://doi.org/10.5281/zenodo.3564179).

From Zenodo you can download these (zipped) folders: 

 * `data` - Contains raw data for the analyses
 * `output` - Contains intermediate and final results
 
 Please deposit the unzipped folders in the root directory of this R-project.
 
 **Exceptions:**
 
 * The raw Human Cell Atlas data (`"data/hca_data/expression_data/sce.all.technologies.RData"`) are not available on Zenodo. Instead the users can work with the normalized data stored in `"output/hca_data/expression/norm.rds"`. The raw data are accessible on GEO: [GSE133549](GSE133549)
 * Only a small fraction of the raw and normalized expression data of the simulated single cells are  avaiable. All data together exceed the limitation of Zenodo's maximal upload size.

***

### How to cite
> Holland CH, Tanevski J, Gleixner J, Kumar MP, Mereu E, Perales-Patón J, Joughin BA, Stegle O, Lauffenburger DA, Heyn H, Szalai B, Saez-Rodriguez, J. "Robustness and applicability of functional genomics tools on scRNA-seq data." _bioRxiv._ 2019. DOI: [10.1101/753319](https://doi.org/10.1101/753319).

***

### Analyses & Scripts
#### Testing robustness with respect to low gene coverage (Fig. 1)
Analysis script available [here](https://github.com/saezlab/FootprintMethods_on_scRNAseq/blob/master/analyses/general_robustness.Rmd).

#### In-silico benchmark (Fig. 2)
Analysis script available [here](https://github.com/saezlab/FootprintMethods_on_scRNAseq/blob/master/analyses/in_silico_benchmark.Rmd).

#### In-vitro benchmark (Fig. 3)
Analysis script available [here](https://github.com/saezlab/FootprintMethods_on_scRNAseq/blob/master/analyses/in_vitro_benchmark.Rmd).

#### Analysis of HCA-data (Fig. 4)
Analysis script available [here](https://github.com/saezlab/FootprintMethods_on_scRNAseq/blob/master/analyses/hca_data_analysis.Rmd).

#### Plotting
We provide also scripts for plotting [individual](https://github.com/saezlab/FootprintMethods_on_scRNAseq/blob/master/analyses/plot_figures.Rmd) figures and [collages](https://github.com/saezlab/FootprintMethods_on_scRNAseq/blob/master/analyses/figure_arrangement.Rmd).
