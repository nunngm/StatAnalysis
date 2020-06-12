setwd("C:\\Users\\garre\\OneDrive\\Documents\\Cameron Lab- McMaster University\\Data\\Data-ARR RNA-seq\\Exp-R workshop")
files <- file.path("counts", list.files("counts"))


samples <- c("m_mg_0h_s1", "m_mg_0h_s2","m_mg_0h_s3","m_mg_12h_s1","m_mg_12h_s2","m_mg_12h_s3","m_mg_24h_s1","m_mg_24h_s2","m_mg_24h_s3","m_pst_0h_s1","m_pst_0h_s2","m_pst_0h_s3","m_pst_12h_s1","m_pst_12h_s2","m_pst_12h_s3","m_pst_24h_s1","m_pst_24h_s2","m_pst_24h_s3","y_mg_0h_s1","y_mg_0h_s2","y_mg_0h_s3","y_mg_12h_s1","y_mg_12h_s2","y_mg_12h_s3","y_mg_24h_s1","y_mg_24h_s2","y_mg_24h_s3","y_pst_0h_s1","y_pst_0h_s2","y_pst_0h_s3","y_pst_12h_s1","y_pst_12h_s2","y_pst_12h_s3","y_pst_24h_s1","y_pst_24h_s2","y_pst_24h_s3")
names(files) <- samples

sub.0 = c(1:3,10:12,19:21,28:30)
sub.12 = c(4:6,13:15,22:24,31:33)
sub.24 = c(7:9,16:18,25:27,34:36)

##Full data set

infection = c(rep("mg",9),rep("pst",9),rep("mg",9),rep("pst",9))
infection <- as.factor(infection)
length(infection)

hpi = c(rep(0,3),rep(12,3),rep(24,3))
hpi = c(hpi,hpi,hpi,hpi)
hpi = as.factor(hpi)
length(hpi)

age = c(rep("m",18),rep("y",18))
age = as.factor(age)
length(age)


design.full <- data.frame(sample=names(files),
                          file=files,
                          age=age,
                          infection =infection,
                          hpi=hpi
                )

design.full

##very general model that just has all the sample information
model.full <- formula(~age+infection+hpi)

all.data <- DESeqDataSetFromHTSeqCount(design.full,design=model.full)
all.data$infection = factor(all.data$infection, levels =  c("mg","pst"))
all.data$hpi = factor(all.data$hpi, levels = c("0","12","24"))
all.data$age = factor(all.data$age, levels = c("y","m"))

#Make group fator and trim any row that has less than 10 counts total
## This allows us to pick our treatment group (ie ymg0h rather than all of each type (all of 0 hpi or all of youbg or all of mock))
all.data$group <- factor(paste0(all.data$age,all.data$infection,all.data$hpi))
all.data$group = factor(all.data$group,levels=c("ymg0","ymg12","ymg24","ypst0","ypst12","ypst24","mmg0","mmg12","mmg24","mpst0","mpst12","mpst24"))
keep <- rowSums(counts(all.data)) >= 360
all.data <- all.data[keep,]

#Lets redi the model design to encompass every sample individually
all.data@design = ~group
all.data2 <- DESeq(all.data)

#Set-up the raw data comparison
#contrast = c(intgroup,Var1, Var2) -> contrasts Var 1 to reference (Var2)
#pAdjustMethod ="BH" -> is the false discovery rate method
## lfc Shrinkage doesn't matter here as p-values are calculated based on unshrunken values therefore shrinkage should only be used for data visualisation and for ranking of genes
res.y= results(all.data2,contrast = c("group", "ypst24", "ymg24"), alpha = 0.05, pAdjustMethod="BH")
temp = as.data.frame(res.y@listData)
rownames(temp) = res.y@rownames
res.y = temp

res.m = results(all.data2,contrast = c("group","mpst24","mmg24"),alpha =0.05, pAdjustMethod = "BH")
temp = as.data.frame(res.m@listData)
rownames(temp) = res.m@rownames
res.m = temp

##Collect up and down genes
y.up = res.y[res.y$log2FoldChange > 0,]$padj #collect genes where fold change is positive (up-regulated)
names(y.up) <- rownames(res.y[res.y$log2FoldChange > 0,]) #collect the genenames
y.up <- y.up[complete.cases(y.up)] #remove any NAs

y.down = res.y[res.y$log2FoldChange < 0,]$padj #see above but for down-reg genes
names(y.down) <- rownames(res.y[res.y$log2FoldChange < 0,]) 
y.down <- y.down[complete.cases(y.down)]

m.up = res.m[res.m$log2FoldChange > 0,]$padj 
names(m.up) <- rownames(res.m[res.m$log2FoldChange > 0,]) 
m.up <- m.up[complete.cases(m.up)]

m.down = res.m[res.m$log2FoldChange < 0,]$padj 
names(m.down) <- rownames(res.m[res.m$log2FoldChange < 0,]) 
m.down <- m.down[complete.cases(m.down)]