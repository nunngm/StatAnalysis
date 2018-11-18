#ANOVA, testing assumptions, and Tukey's post-hoc test 

#install packages, only needs to be done once per computer
install.packages("car")
install.packages("agricolae") 

#clear the global environment
rm(list=ls()) 

#import data file, an excel spreadsheet saved in csv format 
setwd ("C:/Users/garre/Desktop") #set the working directory to wherever your file is located 
mydata <- read.csv("mydata.csv", header=TRUE) #tell R to read the data file and add it to the global environment 
#label variable
mylabel =data.frame(groups=c("mut_m","mut_y","wt_m","wt_y"))

#set the variables as factors or numeric
mydata$treatment <- factor(mydata$treatment)
mydata$cfu <- as.numeric(mydata$cfu)
mean_cfu= tapply(mydata$cfu,mydata$treatment,mean)
sd_cfu=tapply(mydata$cfu,mydata$treatment,sd)
data_plot=rbind(mean_cfu,sd_cfu)
data_plot=t(data_plot)
data_plot= cbind(mylabel,data_plot)

#plotting the graph
library(ggplot2)
ggplot(data=data_plot,aes(x= groups,y=mean_cfu))+geom_bar(stat="identity")+scale_x_discrete(limits=c("wt_m","mut_m","wt_y","mut_y"))

#obtain some basic outputs  
print(mydata) 
plot(mydata$treatment,mydata$cfu) 
tapply(mydata$cfu,mydata$treatment,mean)
tapply(mydata$cfu,mydata$treatment,sd)
tapply(mydata$cfu,mydata$treatment,length) 

#ANOVA
my_anova_model = aov(mydata$cfu ~ mydata$treatment, data=mydata) #the numeric variable must be first in this expression 
summary(my_anova_model) #show results, NB, before interpreting the results you need to check that your data meet the assumptions of ANOVA, however, to do this you must first run the ANOVA

#TESTING ASSUMPTIONS

#test for (distribution) normality of the data (numeric variable), null hypothesis: data are normally distributed 
shapiro.test(mydata$cfu) 

#test for normality of the residuals
shapiro.test(resid(my_anova_model)) 

#tests for homegeneity of variance, null hypothesis: homogeneity of variance 
library(car) 
leveneTest(mydata$cfu ~ mydata$treatment, data=mydata) 
bartlett.test(mydata$cfu ~ mydata$treatment, data=mydata)

#POST-HOC TESTING 
#only if your data satisfy the assumptions and the ANOVA returned a significant result
library(agricolae)
TukeyHSD(my_anova_model) #gives p values for pair-wise comparisons 
print(HSD.test(my_anova_model, alpha=0.05, "mydata$treatment")) #gives letter codes but no p values


--

#if your data fail to meet any of the assumptions try a log transformation and run an ANOVA on the log-transformed data  
mydata$logcfu = log10(mydata$cfu)
mydata$logcfu <- as.numeric(mydata$logcfu)

#ANOVA (log-tranformed data)
my_anova_model_2 = aov(mydata$logcfu ~ mydata$treatment, data=mydata) 
summary(my_anova_model_2)

#TESTING ASSUMPTIONS (log-transformed data)
shapiro.test(mydata$logcfu)  
shapiro.test(resid(my_anova_model_2))
levene.test(mydata$logcfu ~ mydata$treatment, data=mydata)
bartlett.test(mydata$logcfu ~ mydata$treatment, data=mydata)

#POST-HOC TESTING (log-transformed data)
#only if your data satisfy the assumptions and the ANOVA returned a significant result
library(agricolae)
TukeyHSD(my_anova_model_2) #gives p values for pair-wise comparisons 
print(HSD.test(my_anova_model_2, alpha=0.05, "mydata$treatment")) #gives letter codes but no p values

