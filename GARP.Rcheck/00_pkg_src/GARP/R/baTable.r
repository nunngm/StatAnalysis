#' Blank and summarize biofilm formation
#'
#' Takes in the results from the machince reading, narrows down the time points and blanks all values.
#' @param filename A string that corresponds to the name fo the raw data file. Timepoints should be entered manually before analysis.
#' @param blankName The regex expression used to find the blank wells
#' @param empty The regex expression used to find/remove empty wells
#' @note Omits and NA data by removing that entire timepoint.
#' @return A table of summarized and blanked readings. Writes this table to a csv file titled "analyzed_[your filename]"
#' @export
baTable = function(filename, blankName="^blank.*", empty ="empty"){
  mydata <- read.table(filename,sep = ",", header=TRUE,row.names = 1, check.names = F)
  #rownames(mydata)= (0:(length(mydata[,1])-1))*15/60 #Make the row names representative of hours
  mydata = na.omit(mydata)
  #Fix column names
  mydata = mydata[colnames(mydata)!=empty] #remove columns labelled empty
  #Mean subtract the mean of the blanks from each row
  blankAvg = rowMeans(mydata[,grepl(pattern = blankName,colnames(mydata))])
  mydata = mydata[,grep(pattern = blankName,colnames(mydata),invert =T)] - blankAvg
  #mydata = mydata[,order(colnames(mydata))] #Sort by column name
  write.csv(mydata,paste0("analyzed_",filename))
  mydata
}
