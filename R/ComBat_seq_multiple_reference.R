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
#' @param reference_set A vector of batch that will be used as reference in the batch effect correction.

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










ComBat_seq_multiple_reference <- function(counts, batch, group=NULL, covar_mod=NULL, full_mod=TRUE, 
                       shrink=FALSE, shrink.disp=FALSE, gene.subset.n=NULL,reference_set=NULL) {
if (isFalse(all(reference_set%in%batch))) {
stop('There are batch in reference set that are not in the overall batch information.')
}
#### Subset out the referencebatch. ####
referencebatch_batch=batch[which(batch%in%reference_set)]
referencebatch_count=counts[,which(batch%in%reference_set)]
if (is.null(group)) {referencebatch_group=NULL}
else{referencebatch_group=group[which(batch%in%reference_set)]
}
if (is.null(covar_mod)) {referencebatch_covar_mod=NULL}
else{referencebatch_covar_mod=covar_mod[which(batch%in%reference_set),]}

corrected_reference_batch=ComBat_seq(referencebatch_count,referencebatch_batch,referencebatch_group,referencebatch_covar_mod,full_mod=full_mod,
shrink=shrink,shrink.disp=shrink.disp,gene.subset.n=gene.subset.n)

#### Next step is to batch correct all other dataset individually against this one reference batch. ####
results=corrected_reference_batch
for (i in unique(batch[batch%in%reference_set])) {
query_batch=batch[which(batch%in%i)]
query_count=counts[,which(batch%in%i)]
if (is.null(group)) {query_group=NULL}
else{query_group=group[which(batch%in%i)]
query_group=c(referencebatch_group,query_group)
}
if (is.null(covar_mod)) {query_covar_mod=NULL}
else{query_covar_mod=covar_mod[which(batch%in%i),]
query_covar_mod=rbind(referencebatch_covar_mod,query_covar_mod)
}
query_batch=c(referencebatch_batch,query_batch)
query_count=cbind(corrected_reference_batch,query_count)

query=ComBat_seq(query_count,query_batch,query_group,query_covar_mod,full_mod=full_mod,
shrink=shrink,shrink.disp=shrink.disp,gene.subset.n=gene.subset.n)
query=query[,!colnames(query)%in%]

}







}