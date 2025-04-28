# R version 4.2.0

# Package versions:
# - readr 2.1.5
# - jsonlite 2.0.0
# - tximport 1.26.1
# - DESeq2 1.38.3

setwd(".")

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

renv::activate(project = ".")
renv::restore(prompt = FALSE)

library(readr)
library(jsonlite)
library(tximport)
library(DESeq2)

mpolymorpha_meristem_activity_transcriptomic_analysis <- function(stage, in_dir, out_dir, sample_info_filename, tpm = FALSE, deg = FALSE) {
    print(paste0("Processing stage: ", stage))

    sampleinfo_all <- read.csv(sample_info_filename, header = TRUE, sep = ",")
    
    if (stage == "all") {
      sampleinfo <- sampleinfo_all
    } else {
      print(paste0("Filtering for stage: ", stage))
      sampleinfo <- sampleinfo_all[sampleinfo_all$stage == stage, ]
    }

    files <- file.path(in_dir, sampleinfo$quant_file)
    names(files) <- sampleinfo$quant_file
    if (all(file.exists(files))) {
        print("All files exist. Continuing with analysis")
    } else {
        stop("Some files are missing... Check the input folder")
    }
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    tx2gene <- read_csv(paste0(in_dir, "tx2gene.csv"))
    txi <- tximport(files, type = "salmon", tx2gene=tx2gene)

    if (tpm == TRUE) {
        gene_tpms <- txi$abundance
        single <- gene_tpms["Mp1g14150", ]
        write.csv(single, file = paste0(out_dir, "gene_level_TPMs_Mp1g14150.csv"))
    }                       

    if (deg == TRUE) {
        sampleinfo$activity <- factor(sampleinfo$activity, levels = c("inactive", "active"))
        sampleinfo$pair <- factor(sampleinfo$pair)

        ddsTxi <- DESeqDataSetFromTximport(txi,
                                        colData = sampleinfo,
                                        design = ~ pair + activity)

        ddsTxi$activity <- relevel(ddsTxi$activity, "inactive")
        ddsDESeq <- DESeq(ddsTxi)

        res <- results(ddsDESeq)
        res <- res[order(res$padj),]
        results_filename <- paste0(out_dir, stage, "_DEG.csv")
        write.csv(as.data.frame(res), file=results_filename)

        res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
        res_sig <- res[abs(res$log2FoldChange) > 2 & res$padj < 0.05, ]
        res_sig <- res_sig[order(res_sig$log2FoldChange, decreasing=TRUE), ]
        sig_results_filename <- paste0(out_dir, stage, "_DEG_sig.csv")
        write.csv(as.data.frame(res_sig), file=sig_results_filename)
    
    }

    print(paste0("Finished processing stage: ", stage))
}

# Defining inputs and outputs
in_dir <- 'data/'
sample_info_filename <- paste0(in_dir, "sampleinfo_meristems.csv")
out_dir <- 'output/'

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Extracting gene-level TPMs for Mp1g14150
mpolymorpha_meristem_activity_transcriptomic_analysis("all", 
    in_dir = in_dir, 
    out_dir = out_dir,
    sample_info_filename = sample_info_filename, 
    tpm = TRUE)

# Performing differential gene expression analysis for each meristem stage
for (stage in list('S1', 'S2', 'S3')) {
    mpolymorpha_meristem_activity_transcriptomic_analysis(stage, 
        in_dir = in_dir, 
        out_dir = out_dir,
        sample_info_filename = sample_info_filename,  
        deg = TRUE)  
}
