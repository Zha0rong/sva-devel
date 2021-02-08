#' Adjust for batch effects using an empirical Bayes framework in RNA-seq raw counts
#' 
#' ComBat_seq is an improved model from ComBat using negative binomial regression, 
#' which specifically targets RNA-Seq count data.
#' 
#' @param counts Raw count matrix from genomic studies (dimensions gene x sample) 
#' @param batch Vector / factor for batch
#' @param group Vector / factor for biological condition of interest 
#' @param covar_mod Model matrix for multiple covariates to include in linear model (signals from these variables are kept in data after adjustment) 
#' @param full_mod Boolean, if TRUE include condition of interest in model
#' @param shrink Boolean, whether to apply shrinkage on parameter estimation
#' @param shrink.disp Boolean, whether to apply shrinkage on dispersion
#' @param gene.subset.n Number of genes to use in empirical Bayes estimation, only useful when shrink = TRUE
#' 
#' @return data A gene x sample count matrix, adjusted for batch effects.
#' 
#' @importFrom edgeR DGEList estimateGLMCommonDisp estimateGLMTagwiseDisp glmFit glmFit.default getOffset
#' @importFrom stats dnbinom lm pnbinom qnbinom
#' @importFrom utils capture.output
#' 
#' @examples 
#' 
#' count_matrix <- matrix(rnbinom(400, size=10, prob=0.1), nrow=50, ncol=8)
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0,1), 4)
#' 
#' # include condition (group variable)
#' adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group, full_mod=TRUE)
#' 
#' # do not include condition
#' adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=NULL, full_mod=FALSE)
#' 
#' @export
#' 

ComBat_seq <- function(counts, batch, group=NULL, covar_mod=NULL, full_mod=TRUE, 
                       shrink=FALSE, shrink.disp=FALSE, gene.subset.n=NULL){  
  ########  Preparation  ########  
  ## Does not support 1 sample per batch yet
  #batch <- as.factor(batch)
  #if(any(table(batch)<=1)){
  #  stop("ComBat-seq doesn't support 1 sample per batch yet")
  #}
  
  consensus_gene_vector=list()
  for (i in 1:ncol(batch)) {
    individual_batch_name=colnames(batch)[i]
    individual_batch <- as.factor(batch[,i])
    individual_keep_lst <- lapply(c(levels(individual_batch)), function(b){
      which(apply(counts[, individual_batch==b], 1, function(x){!all(x==0)}))
    })
    individual_keep_lst=Reduce(intersect,individual_keep_lst)
    consensus_gene_vector[[individual_batch_name]]=individual_keep_lst
  }
  consensus_gene_vector=Reduce(intersect,consensus_gene_vector)
  counts <- counts[consensus_gene_vector, ]
  
  multiple_batch_list=list()
  for (i in 1:ncol(batch)) {
    individual_batch_name=colnames(batch)[i]
    individual_batch <- as.factor(batch[,i])


    individual_counts <- counts
    individual_dge_obj <- DGEList(counts=individual_counts)
    
    individual_n_batch <- nlevels(individual_batch)  # number of batches
    individual_batches_ind <- lapply(1:individual_n_batch, function(i){which(individual_batch==levels(individual_batch)[i])}) # list of samples in each batch  
    individual_n_batches <- sapply(individual_batches_ind, length)
    individual_n_sample <- sum(individual_n_batches)
  
    
    
    
    individual_batchmod <- model.matrix(~-1+individual_batch)  # colnames: levels(batch)
    individual_group <- as.factor(group)
    if(full_mod & nlevels(individual_group)>1){
      cat("Using full model in ComBat-seq.\n")
      mod <- model.matrix(~individual_group)
    }
    else{
      cat("Using null model in ComBat-seq.\n")
      mod <- model.matrix(~1, data=as.data.frame(t(counts)))
    }
    if(!is.null(covar_mod)){
      if(is.data.frame(covar_mod)){
        covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), function(i){model.matrix(~covar_mod[,i])}))
      }
      covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x){all(x==1)})]
    }
    individual_mod <- cbind(mod, covar_mod)
    individual_design <- cbind(individual_batchmod, individual_mod)
    
    check <- apply(individual_design, 2, function(x) all(x == 1))
    #if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
    individual_design <- as.matrix(individual_design[,!check])
    
    individual_disp_common <- sapply(1:individual_n_batch, function(i){
      if((individual_n_batches[i] <= ncol(individual_design)-ncol(individual_batchmod)+1) | qr(individual_mod[individual_batches_ind[[i]], ])$rank < ncol(individual_mod)){ 
        # not enough residual degree of freedom
        return(estimateGLMCommonDisp(individual_counts[, individual_batches_ind[[i]]], design=NULL, subset=nrow(individual_counts)))
      }else{
        return(estimateGLMCommonDisp(individual_counts[, individual_batches_ind[[i]]], design=individual_mod[individual_batches_ind[[i]], ], subset=nrow(individual_counts)))
      }
    })
    individual_genewise_disp_lst <- lapply(1:individual_n_batch, function(j){
      if((individual_n_batches[j] <= ncol(individual_design)-ncol(individual_batchmod)+1) | qr(individual_mod[individual_batches_ind[[j]], ])$rank < ncol(individual_mod)){
        # not enough residual degrees of freedom - use the common dispersion
        return(rep(individual_disp_common[j], nrow(individual_counts)))
      }
      else{
        return(estimateGLMTagwiseDisp(individual_counts[, individual_batches_ind[[j]]], design=individual_mod[individual_batches_ind[[j]], ], 
                                      dispersion=individual_disp_common[j], prior.df=0))
      }
    })
    names(individual_genewise_disp_lst) <- paste0('batch', levels(individual_batch))
    individual_phi_matrix <- matrix(NA, nrow=nrow(individual_counts), ncol=ncol(individual_counts))
    for(k in 1:individual_n_batch){
      individual_phi_matrix[, individual_batches_ind[[k]]] <- vec2mat(individual_genewise_disp_lst[[k]], individual_n_batches[k]) 
    }
    
    
    
    
    cat("Fitting the GLM model\n")
    glm_f <- glmFit(individual_dge_obj, design=individual_design, 
                    dispersion=individual_phi_matrix, prior.count=1e-4) #no intercept - nonEstimable; compute offset (library sizes) within function
    alpha_g <- glm_f$coefficients[, 1:individual_n_batch] %*% as.matrix(individual_n_batches/individual_n_sample) #compute intercept as batch-size-weighted average from batches
    new_offset <- t(vec2mat(getOffset(individual_dge_obj),
                            nrow(individual_counts))) +   # original offset - sample (library) size
      vec2mat(alpha_g, ncol(individual_counts))  # new offset - gene background expression # getOffset(dge_obj) is the same as log(dge_obj$samples$lib.size)
    glm_f2 <- glmFit.default(individual_dge_obj$counts, 
                             design=individual_design, 
                             dispersion=individual_phi_matrix, 
                             offset=new_offset, prior.count=1e-4) 
    
    gamma_hat <- glm_f2$coefficients[, 1:individual_n_batch]
    mu_hat <- glm_f2$fitted.values
    phi_hat <- do.call(cbind, individual_genewise_disp_lst)
    if(shrink){
      cat("Apply shrinkage - computing posterior estimates for parameters\n")
      mcint_fun <- monte_carlo_int_NB
      monte_carlo_res <- lapply(1:individual_n_batch, function(ii){
        if(ii==1){
          mcres <- mcint_fun(dat=counts[, individual_batches_ind[[ii]]], mu=mu_hat[, individual_batches_ind[[ii]]], 
                             gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)
        }else{
          invisible(capture.output(mcres <- mcint_fun(dat=counts[, individual_batches_ind[[ii]]], mu=mu_hat[, individual_batches_ind[[ii]]], 
                                                      gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)))
        }
        return(mcres)
      })
      names(monte_carlo_res) <- paste0('batch', levels(individual_batch))
      
      gamma_star_mat <- lapply(monte_carlo_res, function(res){res$gamma_star})
      gamma_star_mat <- do.call(cbind, gamma_star_mat)
      phi_star_mat <- lapply(monte_carlo_res, function(res){res$phi_star})
      phi_star_mat <- do.call(cbind, phi_star_mat)
      
      if(!shrink.disp){
        cat("Apply shrinkage to mean only\n")
        phi_star_mat <- phi_hat
      }
    }else{
      cat("Shrinkage off - using GLM estimates for parameters\n")
      gamma_star_mat <- gamma_hat
      phi_star_mat <- phi_hat
    }
    multiple_batch_list[[individual_batch_name]]=list(common_dispersion=individual_disp_common,
                                                      gene_dispersion=individual_genewise_disp_lst,
                                                      phi_matrix=individual_phi_matrix,
                                                      gamma_hat=gamma_hat,
                                                      mu_hat=mu_hat,
                                                      phi_hat=phi_hat,
                                                      gamma_star_mat=gamma_star_mat,
                                                      phi_star_mat=phi_star_mat)
    mu_star <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
    for(jj in 1:individual_n_batch){
      mu_star[, individual_batches_ind[[jj]]] <- exp(log(mu_hat[, individual_batches_ind[[jj]]])-vec2mat(gamma_star_mat[, jj], individual_n_batches[jj]))
    }
    phi_star <- (phi_star_mat[,2])
    multiple_batch_list[[individual_batch_name]]=list(common_dispersion=individual_disp_common,
                                                      gene_dispersion=individual_genewise_disp_lst,
                                                      phi_matrix=individual_phi_matrix,
                                                      gamma_hat=gamma_hat,
                                                      mu_hat=mu_hat,
                                                      phi_hat=phi_hat,
                                                      gamma_star_mat=gamma_star_mat,
                                                      phi_star_mat=phi_star_mat,
                                                      mu_star=mu_star,
                                                      phi_star=phi_star)
    
    
    
    
    
  }
  
  #The multiple_batch_list is a list that follows the following structures:
    # Batch name:
      ## common_dispersion vector
      ## gene dispersion matrix
      ## phi_matrix
  
  Concatenated_batches=(apply(batch[,colnames(batch)],1,paste,collapse="-"))
  
  Concatenated_batches <- as.factor(Concatenated_batches)
  
  Concatenated_n_batch <- nlevels(Concatenated_batches)  # number of batches
  Concatenated_batches_ind <- lapply(1:Concatenated_n_batch, function(i){which(Concatenated_batches==levels(Concatenated_batches)[i])}) # list of samples in each batch  

  ########  Adjust the data  ########  
  cat("Adjusting the data\n")
  adjust_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(i in 1:Concatenated_n_batch){
    counts_sub <- counts[, Concatenated_batches_ind[[i]]]

    old_mu=multiple_batch_list[[colnames(batch)[1]]][['mu_hat']][, Concatenated_batches_ind[[i]]]
    for ( j in 2:ncol(batch)) {
      old_mu=old_mu+multiple_batch_list[[colnames(batch)[j]]][['mu_hat']][, Concatenated_batches_ind[[i]]]
    }
    
    

    batch_components=as.character(unlist(strsplit(levels(Concatenated_batches)[i],split = '-')))
    old_phi=rep(0,nrow(counts))
    for ( j in 1:ncol(batch)) {
      old_phi=old_phi+multiple_batch_list[[colnames(batch)[j]]][['phi_hat']][, paste('batch',batch_components[j],sep = '')]
    }
    

    new_mu=multiple_batch_list[[colnames(batch)[1]]][['mu_star']][, Concatenated_batches_ind[[i]]]
    for ( j in 2:ncol(batch)) {
      new_mu=new_mu+multiple_batch_list[[colnames(batch)[j]]][['mu_star']][, Concatenated_batches_ind[[i]]]
   }
    
    new_phi <- multiple_batch_list[[colnames(batch)[1]]][['phi_star']]
   for ( j in 2:ncol(batch)) {
      new_phi=new_phi+multiple_batch_list[[colnames(batch)[j]]][['phi_star']]
    }
    
    adjust_counts[, Concatenated_batches_ind[[i]]] <- match_quantiles(counts_sub=counts_sub, 
                                                          old_mu=old_mu, old_phi=old_phi, 
                                                          new_mu=new_mu, new_phi=new_phi)
  }
  
  dimnames(adjust_counts) <- dimnames(counts)
  #return(adjust_counts)
  
  ## Add back genes with only 0 counts in any batch (so that dimensions won't change)
  #adjust_counts_whole <- matrix(NA, nrow=nrow(countsOri), ncol=ncol(countsOri))
  #dimnames(adjust_counts_whole) <- dimnames(countsOri)
  #adjust_counts_whole[keep, ] <- adjust_counts
  #adjust_counts_whole[rm, ] <- countsOri[rm, ]
  return(adjust_counts)
}
