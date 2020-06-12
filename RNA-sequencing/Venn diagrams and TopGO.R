goi.y = y.down #sets the genes of interest for young plants
goi.m = m.down #sets the genes of interested for mature plants
goi.mm = mm.up
goi.yy = yy.up

#Filter out genes that are not significant
sg.y = names(goi.y[goi.y<0.05])
sg.m = names(goi.m[goi.m<0.05])
sg.mm = names(goi.mm[goi.mm<0.05])
sg.yy =names(goi.yy[goi.yy<0.05])

#All genes sets all the unique genes in the environment and can be used to compare back too
allGenes = c(sg.y,sg.m)
allGenes = unique(allGenes)

#Set the colours for the venn diagram
colors = c("#c3db0f","#2cff21","#6b7fff","#ff4059","#de4dff")
venn.diagram(x = list(sg.y,sg.m,sg.mm),
             category.names = c("Y.Ps>Y.Mo","M.Ps>M.Mo","M.Mo>Y.Mo"),
             filename = "3venn12h.tiff",
             output = T,
             imagetype = "tiff",
             scaled = F,
             col = "black",
             fill = colors[1:3],
             cat.col = "black",
             cat.cex = 2,
             cat.dist = c(0.12, 0.12,0.04),
             margin =0.15)
#Two-way venn diagram - More useful
venn.diagram(x = list(sg.y,sg.m),
             category.names = c("Y.Mo>Y.Ps","M.Mo>M.Ps"),
             filename = "temp.tiff",
             output = T,
             inverted = T,
             rotation.degree = 180, #In this program the labels are stationary while the venn diagram is not. The circle with more genes will be put on the left therefore we have to rotate it so that the labels match the circle
             imagetype = "tiff",
             scaled = F,
             col = "black",
             fill = colors[1:2],
             cat.col = "black",
             cat.cex = 2,
             cat.dist = c(0.16, 0.16),
             margin =0.15)

#Specifies which part of the Venn diagram you would like to complete GO enrichment on
geneList = as.integer(allGenes %in% sg.y)
names(geneList) = allGenes
geneList[geneList==1] = as.integer( names(geneList[geneList==1])%in% sg.m)
geneList[geneList==1] = as.integer( names(geneList[geneList==1])%in% sg.mm)
sum(geneList)

geneList = as.factor(geneList) #This has to be a factor for TopGO

#Instead of just comparing to an environment of DEGs this allow comparison to all genes (which have an average expression >10 reads).
totalGenes = unique(c(names(goi.y),names(goi.m)))
geneList2 = as.integer(totalGenes %in% names(geneList[geneList==1]))
names(geneList2) = totalGenes
sum(geneList2)
geneList2 = as.factor(geneList2)

#Does the GO enrichment (BP, MF, or CC)
GOdata <- new("topGOdata",
              ontology = "BP", 
              allGenes = geneList2,  
              annotationFun = annFUN.gene2GO, 
              gene2GO = gene_GO) 

fisher_test <- new("classicCount", testStatistic = GOFisherTest, name = "fisher_test")
fisher.GO <- getSigGroups(GOdata, fisher_test)
fisher.GO

#Lets see the most significant ones
table.GO <- GenTable(GOdata, Fisher = fisher.GO, topNodes = 50)
table.GO #A table of significant terms

#A tree of the significant terms
par(cex = 0.2)
showSigOfNodes(GOdata, score(fisher.GO), firstSigNodes = 15, useInfo = 'all')

#GO term -> significant genes in goi list
allGO = genesInTerm(GOdata)
sigGenes = lapply(allGO,function(x) x[x %in% names(geneList[geneList==1])] )
objectSymbol[sigGenes[["GO:0009725"]]]

#This outputs a file with all of the genes of interest from each part of a two-way venn diagram as a list with gene names and descriptions
venn.y = as.integer(allGenes %in% sg.y)
names(venn.y) = allGenes
venn.y[venn.y==1] = as.integer(! names(venn.y[venn.y==1])%in% sg.m)
venn.m = as.integer(allGenes %in% sg.m)
names(venn.m) = allGenes
venn.m[venn.m==1] = as.integer(! names(venn.m[venn.m==1])%in% sg.y)
venn.both = as.integer(allGenes %in% sg.y)
names(venn.both) = allGenes
venn.both[venn.both==1] = as.integer( names(venn.both[venn.both==1])%in% sg.m)
sum(venn.y)
sum(venn.both)
sum(venn.m)
venn.y = names(venn.y[venn.y==1])
venn.m = names(venn.m[venn.m==1])
venn.both = names(venn.both[venn.both==1])

venn = list(venn.y,venn.both,venn.m)
venn =data.frame(lapply(venn, "length<-", max(lengths(venn))))
venn = cbind(as.character(venn[,1]),as.character(venn[,1]),as.character(venn[,2]),as.character(venn[,2]),as.character(venn[,3]),as.character(venn[,3]))
venn = cbind(venn[,1:2],desvec[venn[,1]],venn[,3:4],desvec[venn[,3]],venn[,5:6],desvec[venn[,5]])
venn[,2] = objectSymbol[venn[,2]]
venn[,5] = objectSymbol[venn[,5]]
venn[,8] = objectSymbol[venn[,8]]
colnames(venn) = c("Young only","Gene name","Description","Both","Gene name","Description","Mature only","Gene name","Description")
write.csv(venn,"temp.csv")