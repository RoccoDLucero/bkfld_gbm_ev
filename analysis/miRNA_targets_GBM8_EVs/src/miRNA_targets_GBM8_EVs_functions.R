#miRNA_targets_GBM8_EVs_functions.R
#Rocco Lucero 7.13.2017


convert.to.numeric <- function(rna_seq_rpm_df){
    
    rpm <- rna_seq_rpm_df
    
    cnms <- colnames(rpm)
    rnms <- rownames(rpm)
    
    rpm <- apply(X = rpm, MARGIN = 2 , FUN = as.numeric)
    
    dimnames(rpm) <- list(rnms,cnms)
    
    return(rpm)
    
}



filter.by.abundance <- function(sample_expression, min_exp = 0, remove_zeros = T){
    
    sample_expression <- sample_expression[sample_expression >= min_exp]
    
    if(remove_zeros){sample_expression <- sample_expression[sample_expression > 0]}

    return(sample_expression)
}
    
get.quantiles.of.rna.expression <- function(sample_expression_vals, n_quantiles = 100){
    
    qs <- quantile(sample_expression_vals, probs = seq(0, 1, (1/n_quantiles)), names = T )
    
    exp_qntls <- sapply(sample_expression_vals, findInterval, vec = qs)
    
    q_e <- as.data.frame(rbind(sample_expression_vals, exp_qntls))
    
    return(q_e[,order(q_e[1,])])
        
}




match.rna.seq.to.array <- function(array, rnaseq){
    
    get.matched_probes <- function(matched_probes){
        
        matched_probes_len <- sapply(X = matched_probes, length)
        matched_probes <- matched_probes[matched_probes_len == 1]
        matched_probes_arry <- names(matched_probes)
        matched_probes_rseq <- unlist(matched_probes)
        
        cbind(matched_probes_arry, matched_probes_rseq)
        
    }
    
    seq_nms <- colnames(rnaseq)
    arr_nms <- colnames(array)
    
    mtch_5p <- sapply(X = arr_nms,
                      FUN = function(x){grep(paste(x, '-5p', sep = ''),
                                             seq_nms,
                                             ignore.case = T,
                                             value = T)})
    
    mtch_3p <- sapply(X = arr_nms,
                      FUN = function(x){grep(paste(x, '-3p', sep = ''),
                                             seq_nms,
                                             ignore.case = T,
                                             value = T)})
    
    mtch_root <- sapply(X = arr_nms, 
                        FUN = function(x){grep(paste(x, '$', sep = ''),
                                               seq_nms,
                                               ignore.case = T,
                                               value = T)})
    
    mtchd <- lapply(X = list(mtch_root, mtch_5p, mtch_3p), FUN = get.matched_probes)
    
    mtchd <- Reduce(f = rbind, x = mtchd)
    
    mtchd <- mtchd[order(mtchd[,1]),]
    
}


sum.rna.seq.expr.by.array.probe <- function(array_probes, seq_samp){
    
    base_probes <- unique(array_probes)
    
    trg_probes <- names(seq_samp)
    
    sum_rpm <-  sapply(X = base_probes,
                       
                       FUN = function(prb){
                           
                           prbs <- grep(prb, trg_probes, value = T)
                           
                           sum(seq_samp[prbs])
                           
                       }
                       
    )
    
    sum_rpm <- sum_rpm[sum_rpm > 0]
    
    
    return(sum_rpm) 
    
}
