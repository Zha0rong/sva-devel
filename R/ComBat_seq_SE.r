Combat_Seq_SE <- function(Summarize_experiment_object,assay_name=NULL, batch, group=NULL, covar_mod=NULL, full_mod=TRUE, 
                          shrink=FALSE, shrink.disp=FALSE, gene.subset.n=NULL) {
  Results=Summarize_experiment_object
  if (all(c(batch,group,covar)%in%colnames(Results@colData))) { # Check whether all variables are in the coldata, exit if not.
    
    if (assay_name%in%names(assays(Results))){ # Check whether the designated assay is in the assay, use if exist.
      adjusted=ComBat_seq(assays(Results)[[assay_name]],Results@colData[,batch],Results@colData[,group],
                          Results@colData[,covar_mod],full_mod=full_mod,shrink = shrink,shrink.disp = shrink.disp,gene.subset.n = gene.subset.n)
      
      assays(Results)$adjusted=adjusted
    }
    else { # Check whether the designated assay is in the assay, use counts assay if not exist.
      cat(paste('Assay ',assay_name,' not presented, using the counts assay now.'))
      adjusted=ComBat_seq(assays(Results)[['counts']],Results@colData[,batch],Results@colData[,group],
                          Results@colData[,covar_mod],full_mod=full_mod,shrink = shrink,shrink.disp = shrink.disp,gene.subset.n = gene.subset.n)
      
      assays(Results)$adjusted=adjusted
    }
    return(Results)
  }
  else {
    if (!batch%in%colnames(Results@colData)) {
      cat('Batch Variable not included in the summarizeexperiment coldata, check the batch variable name.')
    }
    if (!group%in%colnames(Results@colData)) {
      cat('Group Variable not included in the summarizeexperiment coldata, check the group variable name.')
    }
    if (!all(group%in%colnames(Results@colData))) {
      Missing_covariate = group[!group%in%colnames(Results@colData)]
      for (missing in Missing_covariate) {
        cat(paste('Covariate ',missing, ' not included in the summarizeexperiment coldata, check the Covariate variable name.',sep = ''))
      }
    }
  }
}