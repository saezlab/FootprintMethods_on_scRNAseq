## Robustness and applicability of functional tools on scRNA-seq data

### Abstract
Many tools have been developed to extract functional and mechanistic insight from bulk transcriptome profiling data. With the advent of single-cell RNA sequencing (scRNA-seq), it is in principle possible to do such an analysis for single cells. However, scRNA-seq data has specific characteristics such as drop-out events, low library sizes and a comparatively large number of samples/cells. It is thus not clear if functional genomic tools established for bulk sequencing can be applied on scRNA-seq in a meaningful way. To address this question, we performed benchmark studies on in silico and in vitro single-cell RNA-seq data. We focused on the tools PROGENy and VIPER that estimate pathway and transcription factor (TF) activities, respectively. For the in silico study we simulated single cells from TF/pathway perturbation bulk RNA-seq experiments. Our simulation strategy guaranties that the information of the original perturbation is preserved while resembling the characteristics of scRNA-seq data. We complemented the in silico data with in vitro scRNA-seq data upon CRISPR-mediated knock-out. Our benchmarks on both the simulated and in vitro data revealed comparable performance to the original bulk data. Additionally, we showed that the TF and pathways activities preserve cell-type specific variability by analysing a mixture sample sequenced with 13 scRNA-seq different protocols. Our analyses suggest that functional genomics tools can be used on scRNA-seq data, and provide a benchmark for further methods development by the community. 

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
