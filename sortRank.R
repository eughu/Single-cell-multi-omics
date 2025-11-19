args<-commandArgs(T)


library(plyr)
library(reshape2)
library(pheatmap)
library(HMM)
library(ggplot2)
library(grid)
library(gridExtra)

if(length(args)<1)
{
    cat("Usage:\n")
    cat("\tRscript sortRank.R train\tUsing the training set to get the probability matrix.\n")
    cat("\tRscript sortRank.R sort\t\tUsing the probability matrix to rank the input matrix based on a naive bayes method.\n")
    q()
}

source('/data/Analysis/lilin/bin/MethGC/bin/motif_NDR/nucleo_rank/nucleoRank.R')
if(args[1] == "train")
{
    if(length(args)<2){
        cat("First step, training a dataset.\n")
        cat("\t/usr/local/bin/Rscript sortRank.R train /data/Analysis/huboqiang/tmp/lilin_20151026/motif_based/nucleo_rank/test_hPGC_Week13_rep2_1.NDR_100bin.rmNA_half.withGroup.xls\n")
        q()
    }
    train_set(args[2])
}

if(args[1] == "sort")
{
    if(length(args)<4){
        cat("Second step, using the training result.\n")
        cat("\t/usr/local/bin/Rscript sortRank.R sort /data/Analysis/huboqiang/tmp/lilin_20151026/motif_based/merge_mES_150913_col_1234.CTCF.ChIP_motif_based.mat probMat.0.xls probMat.1.xls\n")
        q()
    }
    sortMat(args[2], file_probMat0=args[3], file_probMat1=args[4], mean_cnt=10)
}



