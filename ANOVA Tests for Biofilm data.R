##AMA-Wrapper
library(GARP)
library(growthcurver)
setwd("C:\\Users\\garre\\OneDrive\\Documents\\Undergrad and Masters\\Data\\Data-PTIGE and Bacterial assays for PTI\\Antimicrobial Growth Assays\\rawdata")
mydata =amaTable("AMA-quad-19-2.csv", everyHour = 4)



splitme= function(x, sep = "|",fix = T, select = 1){
  x = strsplit(y, split = sep, fixed = fix)
  res = as.character()
  for (i in 1:length(x)){
    res = c(res, x[[i]][select])
    print(i)
  }
}
mydata = read.csv("analyzed_AMA-19-2.csv", row.names = 1)
data_time = as.data.frame(lapply(mydata[,1:4],function(x) as.integer( rownames(mydata))))
mat = matrix(nrow=4,ncol= length(mydata)/4)
mat = as.data.frame(mat)
for (i in 1:(length(mydata)/4)){
  res = growthRates(mydata,i)
  name = as.character(names(res)[1])
  mat[,i] = res
  colnames(mat)[i] = name
}
colnames(mat) = lapply(colnames(mat), function(x) gsub('X','',x))

growthRates= function(data,group){
  name = colnames(data)[(group*4)-3]
  name = rep(name,4)
  subset = data[,((group*4)-3):(group*4)]
  res = list()
  for (i in 1:4){
    res = c(res,SummarizeGrowth(as.integer( rownames(data)),subset[,i])$vals$r)
  }
  res = unlist(res)
  names(res) = name
  return(res)
}
lapply()

data_time = as.data.frame(lapply(mydata[,1:4],function(x) as.integer( rownames(mydata))))
SummarizeGrowth(data_time,mydata[,2:89])




mydata= read.table(file= "clipboard",sep= "\t",header =T, check.names = F)

od = c(lapply(unique(colnames(mydata)),function(x){tmp = mydata[1,colnames(mydata)==x]}), lapply(unique(colnames(mydata)),function(x){tmp = mydata[2,colnames(mydata)==x]}))
repnum =4
treatnum = (length(od)/2)-1
df = matrix(nrow = (length(od)-2)*repnum, ncol = 3) # nrow's = length of od starting after defining the strains (Number of concentrations) * number of reps
colnames(df) = c("Density","Strain","Dose")
for (i in 1:treatnum){
  startpos = (i-1)*repnum*2
  for (j in 1:repnum){
    startpos = startpos +1
    df[startpos,1] = as.double(od[[i+1]][j])
    df[startpos,2] = as.character(od[[1]])
    df[startpos,3]=as.integer(colnames(od[[i+1]])[1])
    
    df[startpos+repnum,1] = as.double(od[[treatnum+i+2]][j])
    df[startpos+repnum,2] = as.character(od[[treatnum+2]])
    df[startpos+repnum,3]=as.integer(colnames(od[[treatnum+i+2]])[1])
  }
}
df = as.data.frame(df)
df[,1] =as.double(as.character( df[,1] ))

mixed <- with(df, interaction(Strain, Dose))

anovaModel = aov(Density ~ mixed, data = df)
summary(anovaModel)
#test for normality of the residuals
df.lm = lm(Density ~ mixed, data = df)
par(mfrow = c(2,2))
plot(df.lm)
par(mfrow = c(1,1))

df.res = resid(df.lm)
names(df.res)= df$Strain
shapiro.test(df.res)

shapiro.test(resid(anovaModel)) 
#pass if p>0.05

#tests for homegeneity of variance, null hypothesis: homogeneity of variance 
library(car) 
leveneTest(Density ~ mixed, data = df) #passes if p>0.05 thus variance if homogenous
bartlett.test(Density ~ mixed, data = df) #if you can run this run this but is more fickle

#POST-HOC TESTING 
#only if your data satisfy the assumptions and the ANOVA returned a significant result
library(agricolae)
TukeyHSD(anovaModel) #gives p values for pair-wise comparisons 
print(HSD.test(anovaModel, alpha = 0.05,trt ="mixed")) #gives letter codes but no p values (default alpha = 0.05)
  anovaModel$model$mixed = make.unique(as.character(anovaModel$model$mixed))
