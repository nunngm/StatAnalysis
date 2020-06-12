#' Compare bacterial assays
#'
#' Compares all concentrations to other strains and to the controls at each time point
#' @description Takes in the output of assay table then compares different strains at each timepoint to the same concentration using T-Tests.
#' @param mydata The output of amaTable. Required that headers are included for each well seperated by an underscore and in this format: Concentration_BacterialStrain
#' @param outGroup A regex expression that matches the non-wildtype strain. Default is PS392.
#' @param alpha  The p-value below which results are no longer considered significant.
#' @param pval Raw p-values returned. By default is off.
#' @return 3 data frames that correspond to: Strain to strain, Wt to EtOH, Outgroup to EtOH. In that order.
#' @export

compare2all = function(mydata,outGroup = "^PS392", alpha = 0.05,pval = F){
  colum = colnames(mydata)
  colum = strsplit(colum,split = "_")
  df =data.frame()
  for (i in 1:length(mydata)){
    df =rbind(df, c(colum[[i]][1],strsplit(colum[[i]][2], split = ".", fixed = T)[[1]][1]), stringsAsFactors =F)
    #df =rbind(df, c(colum[[i]][1],colum[[i]][2]), stringsAsFactors =F)
  }
  df[,1] = factor(df[,1])
  df[,2] = factor(df[,2])
  rm(colum)
  mx = matrix(nrow = nrow(mydata), ncol = nlevels(df[,1]))
  colnames(mx)= levels(df[,1])
  rownames(mx) = rownames(mydata)
  wtTT = mx
  mutTT = mx
  for (i in 1:nlevels(df[,1])){
    wt = grep(pattern = paste0("^",levels(df[,1])[i]),df[,1])
    mut = wt[grep(pattern = outGroup,df[wt,2])]

    for (j in 1:nrow(mydata)){
      mx[j,i] = t.test(mydata[j,wt],mydata[j,mut],alternative = "t")$p.value
    }
    if (i>=2){
      EtOH_wt = grep(pattern = paste0("^",levels(df[,1])[1]),df[,1])
      EtOH_mut = EtOH_wt[grep(pattern = outGroup,df[EtOH_wt,2])]
      EtOH_wt =EtOH_wt[EtOH_wt!=EtOH_mut]
      for (j in 1:nrow(mydata)){
        wtTT[j,i] = t.test(mydata[j,EtOH_wt],mydata[j,wt])$p.value
        mutTT[j,i] = t.test(mydata[j,EtOH_mut],mydata[j,wt])$p.value
      }
    }
  }
  if (pval ==F){
    df2 = list(mx<alpha,wtTT<alpha,mutTT<alpha)
  }else{
    df2 = list(mx,wtTT,mutTT)
  }
  df2
}
