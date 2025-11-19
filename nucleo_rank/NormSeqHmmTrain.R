join_str<-function(in_vec,join_str){
  out_str<-in_vec[1]
  if (length(in_vec) > 1){
    for (i in 2:length(in_vec))
    {
      out_str = paste(out_str,in_vec[i],sep=join_str  )
    }
  }
  return(out_str)
}

get_args_infoV1<-function(args_in){
  split.arg<-strsplit(args_in[1],'/',perl=TRUE)
  prefix_strs<-split.arg[[1]][1:length(split.arg[[1]])]
  file_Info<-NULL
  file_Info$prefix<-join_str( prefix_strs[1:length(split.arg[[1]])-1],"/" )
  if ( is.na(file_Info$prefix) ){
    file_Info$prefix<-c("./")
  }
  file_Info$infile<-tail(prefix_strs,n=1)
  split.out<-strsplit(file_Info$infile,'\\.',perl=TRUE)
  prefix_out_strs<-split.out[[1]][1:length(split.out[[1]])]
  file_Info$outprefix<-join_str( prefix_out_strs[1:length(split.out[[1]])-1],"." )
  return( file_Info )
}

get_x_lab<-function(mean_col){
  x_lab <- seq(1, length(mean_col), 1)
  idx_naF <- ! is.na(mean_col)
  return(x_lab[idx_naF])
}

avg_n_surrounding<-function(mean_col, mean_cnt){
  BigMat <- matrix(NA, mean_cnt, length(mean_col)+mean_cnt-1)
  len_col<-length(mean_col)
  for (i in 1:mean_cnt){
    BigMat[i, c(i:(len_col+i-1))] <- mean_col
  }
  x_avg <- colMeans(BigMat, na.rm=TRUE)
  beg <- as.integer(mean_cnt/2)
  x_avg <- x_avg[c(beg:(beg+len_col-1))]
  return(x_avg)
}

avg2nucleosome<-function(delta, pval = 0.25, pval_peak = 0.1){
  
  cutoff_upper = quantile(delta, 1 - pval)
  cutoff_lower = quantile(delta, pval)
  
#  peak_cutoff_upper = quantile(delta, 1 - pval_peak)
#  peak_cutoff_lower = quantile(delta, pval_peak)
  
  is_nucleotide <- array("N", dim = c(1, length(delta)))
  is_nucleotide[c(delta > cutoff_upper)] <- "H"
  is_nucleotide[c(delta < cutoff_lower)] <- "L"
#  is_nucleotide[c(delta > peak_cutoff_upper)] <- "PH"
#  is_nucleotide[c(delta < peak_cutoff_lower)] <- "PL"
  
  is_nucleotide <- as.vector(is_nucleotide)
  return(is_nucleotide)
}

nucleo_vec2df<-function(is_nucleotide, x_lab, out_df_dim){
  vec_val <- array(dim = c(1, out_df_dim[1]*out_df_dim[2]) )
  vec_val[x_lab] <- is_nucleotide
#  vec_val[ vec_val == "PL"] <- 0
  vec_val[ vec_val == "L"] <- 2
  vec_val[ vec_val == "N"] <- 3
  vec_val[ vec_val == "H"] <- 4
#  vec_val[ vec_val == "PH"] <- 4
  vec_val <- as.numeric(vec_val)
  mat_isNucleotide <- matrix(vec_val, nrow = out_df_dim[1], ncol = out_df_dim[2], byrow = TRUE)
  df_isNucleotide <- as.data.frame(mat_isNucleotide)
  return(df_isNucleotide)
}

get_emit <- function(df_state_emit, bw=10){
  l_emit <- lapply( unique(df_state_emit$state), function(i){
    df_state_emit_local <- df_state_emit[df_state_emit$state == i,]
    df_state_emit_local <- df_state_emit_local[!is.na(df_state_emit_local$value),]
    dens_tag <- density(df_state_emit_local$value, bw=bw, from=-100, to=100, n = 201)
    df_tag_tmp <- data.frame("x" = round(dens_tag$x, 0), "y" = dens_tag$y/sum(dens_tag$y))
    df_tag_tmp
  } )
  emit_mat <- sapply( l_emit, "[[", "y" )
  row.names(emit_mat) <- l_emit[[1]]$x
  colnames(emit_mat) <- unique(df_state_emit$state)
  emit_mat <- t(emit_mat)
  return(emit_mat)
}

get_tran<-function(is_nucleotide, order_list = order_list){
    df_pre_next <-data.frame( "next_char" = c(is_nucleotide, "NA"), "pres_char" = c("NA", is_nucleotide) )
    df_pre_next <- df_pre_next[c(2:(dim(df_pre_next)[1]-1)), ]
    df_pre_next_ply <- ddply(df_pre_next, .(pres_char, next_char), summarise, sum=length(pres_char))
    df_pre_next_ply_dcasted<-dcast( df_pre_next_ply,pres_char~next_char,value.var="sum")
    row.names(df_pre_next_ply_dcasted) <- df_pre_next_ply_dcasted$pres_char
    df_pre_next_ply_dcasted$pres_char <- NULL
    df_pre_next_ply_dcasted / rowSums(df_pre_next_ply_dcasted)
    df_pre_next_ply_dcasted[order_list, order_list]
}

cnt_NA<-function(x){
    x<-as.numeric(x)
    return( sum(is.nan(x)) )
}

