# mpolymorpha_meristem_activity_transcriptomics

This repository contains the code and metadata files that are required to reproduce the transcriptomic analysis performed in [Spencer V., Casey C., Dolan L. et al. 2025]() to identify genes that are differentially expressed between inactive and active meristems in *Marchantia polymorpha*. 


## Analysis workflow

1. `cluster_scripts/01_trimmomatic.sbatch`  
Trim Illumina adaptors and low-quality reads from raw reads.

2. `cluster_scripts/02_salmon_index.sbatch`  
Generate salmon index for MpTak_v6.1r1 transcriptome. A k-mer size of `-k 21` was used due to the short read lengths.

3. `cluster_scripts/03_salmon_quant.sbatch`  
Perform transcript-level read quantification in quasi-mapping based mode with salmon.  

4. `deseq2_analysis/04_mpolymorpha_meristem_activity_transcriptomic_analysis.R`  
Extract gene-level TPMs for the Mp1g14150 (Mp*CYP78E1*) gene.  
Perform differential gene expression analysis to identify genes that are significantly differentially expressed (padj < 0.05, |log2FC| > 2) between inactive and active *M. polymorpha* meristems at each stage of meristem inactivation.   


## Environment

The analysis was performed using the [CLIP](https://clip.science/) high-performance computing cluster at the Vienna BioCenter (steps 1-3) and in R version 4.2.0 within RStudio 2022.02.2+485 (step 4). The differential gene expression analysis (step 4) uses `renv` in order to ensure reproducibility.  


## Reproducing the differential gene expression analysis

1. Clone the git repository.  

```
git clone https://github.com/calcasey/m_polymorpha_meristem_activity_transcriptomics.git
```

2. Download the transcript-level read quantification `quant.sf` files from the Gene Expression Omnibus (GEO) data repository using the accession ID: [GSE295341](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295341). Download the files into `deseq2_analysis/data/`.  

3. Edit the `quant_file` column in `deseq2_analysis/data/sampleinfo_meristems.csv` to associate each `quant.sf` file with the relevant stage, meristem activity, and sample pair.  

4. Ensure that `deseq2_analysis/data/` contains:  
- `quant.sf` files for each sample  
- `sampleinfo_meristems.csv` with the correct `quant.sf` filenames  
- `tx2gene.csv` to map each transcript ID from the MpTak_v6.1 reference transcriptome to the corresponding gene ID 

5. Install R version 4.2.0. This can be found in the [CRAN archives](https://cran.r-project.org/bin/).  

6. Run `deseq2_analysis/04_mpolymorpha_meristem_activity_transcriptomic_analysis.R` either in RStudio or using the terminal. The R environment needed to run the analysis is automatically reconstructed upon execution using `renv`.

### Running from the terminal

```
cd deseq2_analysis
Rscript 04_mpolymorpha_meristem_activity_transcriptomic_analysis.R
```

### Running from RStudio

In RStudio, open the project by navigating to `deseq2_analysis/` and clicking `deseq2_analysis.Rproj`. Then, run `deseq2_analysis/04_mpolymorpha_meristem_activity_transcriptomic_analysis.R` in RStudio.  

The expected outputs are gene-level TPMs for the Mp1g14150 (Mp*CYP78E1*) gene and lists of significantly differentially expressed genes between inactive and active *M. polymorpha* meristems at each stage of meristem inactivation.  
