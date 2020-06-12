#Note to self: [Rows,Column]

#Reset and read the data
rm(list=ls())
mydata= read.table(file= "clipboard",sep= "\t",header =T) #header row included with data
newdata=mydata
mydata=newdata

##Transformations of the data to make the reisduals better
mydata[,1] = log(mydata[,1], base = 10) #ln of the growth data to make it usable by a linear model

bf = lm(formula = Biofilm ~ Growth + Source, data = mydata)
bf = lm(formula = Biofilm ~ Source + Growth, data = mydata)
bf2 = lm(formula = Biofilm ~ Growth, data = mydata)
bf3 = lm(formula = Biofilm ~ Source, data = mydata)
summary(bf)
summary(bf2)
summary(bf3)
so = lm(formula = Growth ~ Biofilm + Source, data = mydata)
so = lm(formula = Growth ~ Biofilm, data = mydata)
summary(so)

plot(bf$residuals, pch = 16, col = "red")
plot(mydata$Growth, newdata$Biofilm, pch = 16, cex = 1, col = "blue")

library(ggplot2)
p = ggplot(data = mydata, mapping = aes(x = Growth, y = Biofilm))+geom_point(aes(colour = factor(Source)),cex=2)
p+theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size=1),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank()),
        axis.ticks=element_line(colour = "black", size =1),
        axis.ticks.length = unit(5,"points") ,
        axis.title.y = element_blank(),
        legend.position = "none"
) +geom_smooth(method=glm, color='#2C3E50',linetype="dashed",se=T)

##Drawing a line baesed on pregenerated linear model
p+geom_abline(slope = bf2$coefficients[2], intercept = bf2$coefficients[1],col= "black",cex=1,linetype="dashed")
