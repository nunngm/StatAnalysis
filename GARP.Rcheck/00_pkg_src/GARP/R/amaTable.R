#' Blank and summarize growth data
#'
#' Takes in the results from the machince reading, narrows down the time points and blanks all values.
#' @param filename A string that corresponds to the name fo the raw data file
#' @param everyHour How often should a data point be used, default is every other hour
#' @param blankName The regex expression used to find the blank wells
#' @param empty The regex expression used to find/remove empty wells
#' @return A table of summarized and blanked readings. Writes this table to a csv file titled "analyzed_[your filename]"
#' @export
#' @note Omits any NA data by simply removing that timepoint.

amaTable = function(filename, everyHour=2,blankName="^blank.*",empty= "empty"){
  mydata <- read.table(filename,sep = ",", header=TRUE,row.names = 1, check.names = F)
  rownames(mydata)= (0:(length(mydata[,1])-1))*15/60 #Make the row names representative of hours
  mydata = na.omit(mydata)
  #Fix column names
  mydata = mydata[colnames(mydata)!=empty] #remove columns labelled empty
  #Mean subtract the mean of the blanks from each row
  blankAvg = rowMeans(mydata[,grepl(pattern = blankName,colnames(mydata))])
  mydata = mydata[,grep(pattern = blankName,colnames(mydata),invert =T)] - blankAvg
  #Select only desired rows

  ###this shit doesn't work yet
  mydata = mydata[as.numeric(rownames(mydata))%%everyHour==0,]
  #mydata = mydata[,order(colnames(mydata))] #Sort by column name
  write.csv(mydata,paste0("analyzed_",filename))
  mydata
}
