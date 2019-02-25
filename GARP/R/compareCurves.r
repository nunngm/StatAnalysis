#' Compare all growth curves to each other
#'
#' Wrapper for statmod's compareGrowthCurves
#' @importFrom statmod compareGrowthCurves
#' @import metaSEM
#' @import multcompView
#' @param mydata The output file from baTable or amaTable.
#' @param runs The number of permutations performed by compareGrowthCurves default is 1000.
#' @param order The final reordering of your rows so that it goes from highest to lowest conc. By default set-up for basic antimicrbial growth assay.
#' @param raw If true, returns the raw letters and groups to be checked against the final table.
#' @return Returns a table that shows the letter codes of each of the corresponding groups. Writes the table to an excel file titled "letters_test.csv" by default.
#' @export
#'
compareCurves = function(mydata,filename = "test.csv", runs = 1000,raw =F, order = c(1,3,6,9,2,5,8,4,7,10)){
  mydata = t(mydata) #transpose data to make compatible with compareGrowthCurves
  g= vector()
  for (i in 1:80){
    g = c(g,strsplit(rownames(mydata), split =".", fixed =T)[[i]][1])
  } #to make all of the levels to be analyzed by compareGrowthCurves
  lev = levels(as.factor(g))
  res = compareGrowthCurves(g,mydata, nsim = runs)
  res = vec2symMat(res$P.Value, diag =F) #Turns a matrix of pvals into a symmetrical matrix which can be used by multcompLetters to generate letter codes
  rownames(res) = lev
  x = multcompLetters(res)
  if (raw ==T){
    return(x)
  }
  x = x$Letters
  # A whole bunch of formatting stuff to turn
  y = matrix(data = NA, nrow = 2, ncol = length(x)/2)
  for (i in 1:length(x)) {
    y[((i%%2 -1)*-1)+1,(i+1)%/%2] = x[i]
  }
  lev = levels(as.factor(unlist(strsplit(lev, split ="_"))))
  colnames(y) = grep(pattern = "^[0-9]",lev,value =T)
  rownames(y) = lev[grep(pattern = "^[0-9]", lev, invert = T)]
  y = y[,order]
  colnames(y)= colnames(y)[order]
  write.csv(y,paste0("letters_",filename))
}
