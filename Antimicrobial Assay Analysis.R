#Garrett's Antimicrobial Assay R Script

#Required packages, only needs to be done once per computer
install.packages("car")
install.packages("agricolae") 

#Data needs to be in a table of 96 columns (one for each sample) and all the recorded times below

#Reset and read the data (row,col)
rm(list=ls())
makeAMATable = function(filename, everyHour=2,blankName="Blank"){
  mydata <- read.table(filename,sep = ",", header=TRUE,row.names =1, check.names = F)
  rownames(mydata)= (0:192*15)/60 #Make the row names representative of hours
  #Mean subtract the mean of the blanks from each row
  blankAvg = rowMeans(mydata[,colnames(mydata)==blankName])
  newData = mydata[,colnames(mydata)!=blankName] - blankAvg
  #Select only desired rows
  newData = newData[as.numeric(rownames(newData))%%everyHour==0,]
  newData
}
filename = "AMA-18-2v2.csv"
data = makeAMATable("AMA-18-2v2.csv")
df_new <- as.data.frame(lapply(data, FUN=function(x) x-x[1]))
mydata <- read.table("AMA-18-2v2.csv",sep = ",", header=TRUE,row.names =1, check.names = F)
rownames(mydata)= (0:192*15)/60 #Make the row names representative of hours
#Mean subtract the mean of the blanks from each row
blankAvg = rowMeans(mydata[,colnames(mydata)=="Blank"])
mydata = mydata[9:193,colnames(mydata)!="Blank"] - blankAvg
blankAvg = as.numeric(mydata[1,])
df = data.frame()
for (i in 1:185){
  df = rbind(df,blankAvg)
}
df = mydata-df
df = df+round(mean(as.numeric(mydata[1,])),digits=5)
df = df[as.numeric(rownames(df))%%2==0,]
write.csv(df, "analyzed.csv")


=as.double(newData[1,])+blankAvg
newData2=newData
for (hell in 1:185){
  newData2[hell,] =newData[hell,]-blankAvg
}
newData2= newData2+0.04
#Select only desired rows
newData = newData2[as.numeric(rownames(newData))%%2==0,]

write.csv(newData,"analyzed.csv")
