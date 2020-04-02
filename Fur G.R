reg1= read.table(file= "clipboard",sep= "\t",header =F) #read from clipboard
reg2 = read.table(file= "clipboard",sep= "\t",header =F)

reg1 = as.character(reg1[,1])
reg2= as.character(reg2[,1])
x=0
output = vector()
for (i in reg2){
  if (sum(i==reg1)==1){output=c(output,i)  }
  x=x+1
}

givename = read.table(file= "clipboard",sep= "\t",header =F)
names = vector()
x=1
for (i in givename[,1]){
  if (sum(i==output)==1){names=c(names,as.character(givename[x,2]))}
  x=x+1
}
output = cbind(output,names)
write.csv(output, file="output2.csv")

reg1[,2]= as.character(reg1[,2])
reg2 = reg1[grepl("fur",reg1[,2],ignore.case = T),]


fur= read.table(file= "clipboard",sep= "\t",header =T)
output = cbind(as.character(fur[,1]),0)
for (i in reg1){
  output[grepl(i,output[,1],ignore.case = T),2]= 1
}
sum(as.integer(output[,2]))
reg2 = cbind(reg2,0)
for (i in 1:nrow(reg2)){
  reg2[i,3]= trunc((reg2[i,1]+reg2[i,2])/2)
}
write.csv(reg2, file="output2.csv")
