source("/datd/lilin/Project/NOM-Seq/mouse_NOME_seq/StatInfo_steps/01.plot_geneRegion/03.heatmap/01.motif_within_Chip_peak/motif_based/RenLab_CTCF/try_blue_yellow/tmp20151111/bin/NormSeqHmmTrain.R")

########################
##### Training Set #####
########################
train_set <- function(train_file, mean_cnt = 10){
    file_Info<-get_args_infoV1(train_file)
    setwd(file_Info$prefix)
    data<-read.table(file_Info$infile)
    
    color_palette <- colorRampPalette(c("blue", "yellow"))(100)
    df<-read.table(train_file)
    # LOADING DATA MATRIX
    row.names(df)<-df$V1
    df$V1<-NULL
    rowGroup <- df$V2
    df$V2 <- NULL

    # MATRIX TO VECTOR
    mean_col <- as.vector(unlist(t(df)))
    mean_surround <- avg_n_surrounding(mean_col, mean_cnt)
    x_lab <- get_x_lab(mean_col)
    mean_col <- mean_col[x_lab]
    mean_surround <- mean_surround[x_lab]

    mean_delta <- mean_col - mean_surround
    df_delta <- nucleo_vec2df(is_nucleotide = mean_delta, x_lab = x_lab, out_df_dim = dim(df))

    #plot(colMeans(df_delta[rowGroup==1,], na.rm = 1))
    #lines(colMeans(df_delta[rowGroup==0,], na.rm = 1))


    df_delta$rowGroup <- rowGroup
    df_delta_melt <- melt(df_delta, id.vars = "rowGroup")

    theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(5,1,1,1),"line"))

    p <- ggplot(df_delta_melt, aes(factor(variable), value))
    p<-p+geom_boxplot()+theme + facet_wrap(~rowGroup, nrow=2, ncol=1)
    out_pdf<-paste(file_Info$outprefix,"PCA13","pdf",sep=".")
    pdf(  out_pdf,width = 7,height = 7)
    yy <- grid.arrange(p,nrow=1)
    op <- par(no.readonly=TRUE)
    par(op)
    dev.off()
    l_dfTag <- lapply( unique(df_delta_melt$variable), function(i, bw=5){
      df_delta_melt_test0 <- df_delta_melt[rowGroup == "0",]
      df_delta_melt_test0 <- df_delta_melt_test0[!is.na(df_delta_melt_test0$value),]
      dens_tag0 <- density(df_delta_melt_test0[df_delta_melt_test0$variable==i, ]$value, bw=bw, from=-100, to=100, n = 201)
  
      df_delta_melt_test1 <- df_delta_melt[rowGroup == "1",]
      df_delta_melt_test1 <- df_delta_melt_test1[!is.na(df_delta_melt_test1$value),]
      dens_tag1 <- density(df_delta_melt_test1[df_delta_melt_test1$variable==i, ]$value, bw=bw, from=-100, to=100, n = 201)
  
      df_tag_tmp<-data.frame("x" = round(dens_tag0$x, 0), "y0" = dens_tag0$y/sum(dens_tag0$y), "y1" = dens_tag1$y/sum(dens_tag1$y))
      df_tag_tmp
    } )

    df_pos_prob_0 <- sapply(l_dfTag, "[[", "y0")
    row.names(df_pos_prob_0) <- as.character(l_dfTag[[1]]$x)

    df_pos_prob_1 <- sapply(l_dfTag, "[[", "y1")
    row.names(df_pos_prob_1) <- as.character(l_dfTag[[1]]$x)

    out_pdf1<-paste(file_Info$outprefix,"probMat.0.pdf",sep=".")
    pdf(out_pdf1,width = 7,height = 7)
    pheatmap(df_pos_prob_0, cluster_rows = F, cluster_cols = F, show_rownames = F, color = color_palette)
    dev.off()
    
    out_xls0 <- paste(file_Info$outprefix, "probMat.0.xls", sep = "/")
    write.table(df_pos_prob_0, file = out_xls0, sep = "\t", quote = FALSE)
    out_pdf2<-paste(file_Info$outprefix,"probMat.1.pdf",sep=".")
    pdf(out_pdf2,width = 7,height = 7)
    pheatmap(df_pos_prob_1, cluster_rows = F, cluster_cols = F, show_rownames = F, color = color_palette)
    dev.off()
    
    out_xls1 <- paste(file_Info$outprefix, "probMat.1.xls", sep = "/")
    write.table(df_pos_prob_1, file = out_xls1, sep = "\t", quote = FALSE)
    
    df_out <- data.frame("probMat0" = out_xls0, "probMat1" = out_xls1)
    return(df_out)
}


########################
##### Used for all #####
########################
sortMat <- function(matFile, file_probMat0, file_probMat1, pval = 0.35, pval_peak = 0.1, mean_cnt = 10){
    color_palette <- colorRampPalette(c( "blue", "yellow"))(100)
    df_probMat0 <- read.table(file_probMat0, sep = "\t", header = TRUE)
    df_probMat1 <- read.table(file_probMat1, sep = "\t", header = TRUE)
    
    file_Info<-get_args_infoV1(matFile)
    setwd(file_Info$prefix)
    
    df <- read.table(file_Info$infile)
    row.names(df) <- df$V1
    df$V1 <- NULL
    
    
    df_idx<-(as.numeric(apply(df, 1, cnt_NA)) < dim(df)[2]*0.5)
    df_idx[is.na(df_idx)] = FALSE
    df<-df[df_idx,]
    
    df <- df*100
    
    mean_col <- as.vector(unlist(t(df)))
    mean_surround <- avg_n_surrounding(mean_col, mean_cnt)
    x_lab <- get_x_lab(mean_col)
    mean_col <- mean_col[x_lab]
    mean_surround <- mean_surround[x_lab]
    
    mean_delta <- mean_col - mean_surround
    df_delta <- nucleo_vec2df(is_nucleotide = mean_delta, x_lab = x_lab, out_df_dim = dim(df))

    prob_line <- function(i){
      i <- as.numeric(i)
      prob_0 <- 0
      prob_1 <- 0
      for(pos_i in c(1:100)){
        if (!is.na(i[pos_i])){
          prob_0 = prob_0 + log(df_probMat0[[pos_i]][ row.names(df_probMat0) == as.character(round(i[pos_i], 0))])
          prob_1 = prob_1 + log(df_probMat1[[pos_i]][ row.names(df_probMat1) == as.character(round(i[pos_i], 0))])
        }
      }
      df_out <- data.frame("prob_0" = prob_0, "prob_1" = prob_1)
      return(df_out)
    }
    
    
    rowGroup2 <- apply(df_delta, 1, prob_line)
    prob1 <- sapply(rowGroup2, "[[", "prob_1")
    prob0 <- sapply(rowGroup2, "[[", "prob_0")
    rowGroup2_prob <- ((prob1-prob0)>0)
    delta_prob <- (prob1-prob0)
    order_delta <- as.integer(order(delta_prob, decreasing = T))
    is_nucleotide = avg2nucleosome(delta = mean_delta, pval = pval, pval_peak = pval_peak)
    df_isNucleotide <- nucleo_vec2df(is_nucleotide = is_nucleotide, x_lab = x_lab, out_df_dim = dim(df))
    
    df_state_emit <- data.frame( "state" = is_nucleotide, "value" = round(mean_delta, 0) )
    df_state_emit$state<-factor(df_state_emit$state, level=unique(df_state_emit$state))
    
    emit_mat <- get_emit(df_state_emit, bw=10)
    tran_mat <- as.matrix(get_tran(is_nucleotide, order_list = as.character(unique(df_state_emit$state))))

    out_pdf3<-paste(file_Info$outprefix,"EmitMatrix.pdf", sep=".")    
    pdf( out_pdf3,width = 7,height = 7)
    pheatmap(emit_mat, cluster_rows = F, cluster_cols = F, show_rownames = T, color = color_palette)
    dev.off()

    out_pdf4<-paste(file_Info$outprefix,"TranMatrix.pdf", sep=".")
    pdf( out_pdf4 ,width = 7,height = 7)
    pheatmap(tran_mat, cluster_rows = F, cluster_cols = F, show_rownames = T, color = color_palette)
    dev.off()
    
    hmm = initHMM(rownames(emit_mat), colnames(emit_mat), transProbs = as.matrix(tran_mat), emissionProbs = as.matrix(emit_mat))
    
    viterbi_multiline<-function(in_vec, mean_cnt){
        in_vec <- as.numeric(in_vec[! is.na(in_vec)])
        in_vec_surround <- avg_n_surrounding(in_vec, mean_cnt)
        viterbi = HMM::viterbi(hmm, as.character(round(in_vec-in_vec_surround, 0)))
        return(viterbi)
    }
    
    viterbi <- as.character(unlist(apply(df, 1, viterbi_multiline, mean_cnt = mean_cnt  )))
    df_HMM_viterbi <- nucleo_vec2df(is_nucleotide = viterbi, x_lab = x_lab, out_df_dim = dim(df))
   
    out_xls2<- paste(file_Info$outprefix,"TranMatrix.xls",sep=".") 
    write.table(as.data.frame(hmm$transProbs), file=out_xls2,sep="\t",quote = FALSE)
    out_xls3<- paste(file_Info$outprefix,"EmitMatrix.xls",sep=".")
    write.table(as.data.frame(hmm$emissionProbs), file=out_xls3,sep="\t",quote = FALSE)
    
    df_raw_out <- df[order_delta,]
    df_raw_out$pos <- row.names(df)[order_delta]
    df_raw_out$group <- rowGroup2_prob[order_delta]
    df_raw_out <- df_raw_out[ ,c( dim(df_raw_out)[2]-1, dim(df_raw_out)[2],c(1:(dim(df_raw_out)[2]-2))) ]
    write.table(df_raw_out, file=paste(file_Info$outprefix, "rawDataSort.xls", sep="."), sep="\t", quote = FALSE, row.names=F)
    
    
    df_HMM_viterbi_out <- df_HMM_viterbi[order_delta,]
    df_HMM_viterbi_out$pos <- row.names(df)[order_delta]
    df_HMM_viterbi_out$group <- rowGroup2_prob[order_delta]
    df_HMM_viterbi_out <- df_HMM_viterbi_out[ ,c( dim(df_HMM_viterbi_out)[2]-1, dim(df_HMM_viterbi_out)[2],c(1:(dim(df_HMM_viterbi_out)[2]-2))) ]
    write.table(df_HMM_viterbi_out, file=paste(file_Info$outprefix, "train.xls", sep="."), sep="\t", quote = FALSE, row.names=F)
    
    df_HMM_viterbi_out$pos <- NULL
    group_idx <- df_HMM_viterbi_out$group
    df_HMM_viterbi_out$group <- NULL
    
    
    pdf(paste(file_Info$outprefix, "train.all", "pdf", sep="."))
    pheatmap(  df_HMM_viterbi_out, cluster_rows=F,cluster_cols=F,show_rownames = F,color=color_palette)
    dev.off()
    
    pdf(paste(file_Info$outprefix, "train.Positioned", "pdf", sep="."))
    pheatmap(  df_HMM_viterbi_out[group_idx==1,], cluster_rows=F,cluster_cols=F,show_rownames = F,color=color_palette)
    dev.off()
    
    pdf(paste(file_Info$outprefix, "train.nonPositioned", "pdf", sep="."))
    pheatmap(  df_HMM_viterbi_out[group_idx==0,], cluster_rows=F,cluster_cols=F,show_rownames = F,color=color_palette)
    dev.off()
    
}
