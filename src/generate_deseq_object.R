#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(data.table)
library(DESeq2)
library(tximport)

tx2gene_file <- snakemake@input[["tx2gene"]]
file_list <- snakemake@input[["quant_files"]]

dds_file <- snakemake@output[["dds"]]

# dev
# tx2gene_file <- "output/030_trinity/trinity.raw/Trinity.fasta.gene_trans_map"
# file_list <- list.files("output/040_trinity-abundance/raw",
#                         full.names = TRUE,
#                         recursive = TRUE,
#                         pattern = "quant.sf$")


# tx2gene
tx2gene <- fread(tx2gene_file, header = FALSE)[, .(V2, V1)]

txi <- tximport(file_list, type = "salmon", tx2gene = tx2gene)


# generate col_data
col_data <- data.table(rn = gsub(".*/", "", dirname(file_list)))

# datasets
col_data[grepl("^n[[:digit:]]+", rn), dataset := "angry"]
col_data[!grepl("^n[[:digit:]]+", rn), dataset := "castes"]

# castes
col_data[dataset == "castes", c("caste", "indiv") := tstrsplit(rn, "_")]

# angrywasps
col_data[dataset == "angry",
         splitname := gsub("^(n[[:digit:]]+)([d|f])([[:digit:]]+)",
                           "\\1_\\2_\\3",
                           rn)]
col_data[dataset == "angry",
         c("nest", "caste", "indiv") := tstrsplit(splitname, "_")]
col_data[, splitname := NULL]

# generate deseq object
dds <- DESeqDataSetFromTximport(
  txi,
  colData = data.frame(col_data, row.names = "rn"),
  design = ~ caste )

# write output
saveRDS(dds, dds_file)

# log
sessionInfo()
