#' Compare AMA tables by strain only
#'
#' Compares all concentrations to other strains at each time point
#' @description Takes in the output of AMA table then compares different strains at each timepoint to the same concentration using a T-Test.
#' @param mydata The output of amaTable. Required that headers are included for each well seperated by an underscore and in this format: Concentration_BacterialStrain
#' @param outGroup A regex expression that matches the non-wildtype strain. Default is PS392.
#' @param alpha  The p-value below which results are no longer considered significant.
#' @param pval Raw p-values returned. By default no
#' @return A matrix of each concentration and each p-value.
#' @export

compare2strain = function(mydata,outGroup ="^PS392.*", alpha =0.05,pval = F){
  colum = colnames(mydata)
  colum = strsplit(colum,split = "_")
  df =data.frame()
  for (i in 1:length(mydata)){
    df =rbind(df, c(colum[[i]][1],colum[[i]][2]), stringsAsFactors =F)
    #help = rbind(help,df,stringsAsFactors =F)
  }
  df[,1] = factor(df[,1])
  rm(colum)
  mx = matrix(nrow = nrow(mydata), ncol = nlevels(df[,1]))
  colnames(mx)= levels(df[,1])
  rownames(mx) = rownames(mydata)
  for (i in 1:nlevels(df[,1])){
    wt = grep(pattern = levels(df[,1])[i],df[,1])
    mut = wt[grep(pattern = outGroup,df[wt,2])]
    wt = wt[wt!=mut]
    for (j in 1:nrow(mydata)){
      mx[j,i] = t.test(mydata[j,wt],mydata[j,mut],alternative = "t")$p.value
    }
  }
  if (pval==F){
    mx = mx<alpha
  }
  mx
}
