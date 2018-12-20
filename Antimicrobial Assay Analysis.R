#Garrett's Antimicrobial Assay R Script

#Required packages, only needs to be done once per computer
install.packages("car")
install.packages("agricolae") 

#Data needs to be in a table of 96 columns (one for each sample) and all the recorded times below

#Reset and read the data (row,col)
rm(list=ls())
#before running this program label blank wells "blank" and empty wells "empty"
makeAMATable = function(filename, everyHour=2,blankName="^blank.*"){
  mydata <- read.table(filename,sep = ",", header=TRUE,row.names = 1, check.names = F)
  rownames(mydata)= (0:(length(mydata[,1])-1))*15/60 #Make the row names representative of hours
  mydata = na.omit(mydata)
  #Fix column names
  mydata = mydata[colnames(mydata)!="empty"] #remove columns labelled empty
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
filename = "AMA-18-6.csv"
setwd("C:/Users/garre/OneDrive/Documents/Undergrad and Masters/Data/R WD/RawData")
mydata = makeAMATable("AMA-18-5.csv")

AMATTest = function(mydata){
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
    mut = wt[grep(pattern = "^PS392.*",df[wt,2])]
    wt = wt[wt!=mut]
    for (j in 1:nrow(mydata)){
      mx[j,i] = t.test(mydata[j,wt],mydata[j,mut])$p.value
    }
  }
  mx = mx<0.05
  
}



best = mydata[2,]
best = best+0.05
write.csv(best,paste0("messed_",filename))


x = c("blank","blank.1","blank.2","blank.3")
y = strsplit(x,"\\.")
y = as.data.frame(strsplit(x,"\\."))[1,]

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
