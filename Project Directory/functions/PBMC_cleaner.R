pbmc_QC <- function(data,RNA.feat,RNA.count,exclusions) {

  subset(data, subset = nFeature_RNA < RNA.feat &
    nCount_RNA < RNA.count &
    !(data$Sample_Tag %in% exclusions))
    
}