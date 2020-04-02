#' Take leaf disc information in long format and analyze
#'
#' Takes in the results from the machince reading, narrows down the time points and blanks all values.
#' @import car
#' @import agricolae
#' @param db Database of data in wide format
#' @param check False returns letter codes and assumption tests. True returns long format data, grouped data w std devs, the generated linear model, and raw statistical analysis.
#' @param logIt Log (base 10) the organized data before stats
#' @param twoVar Is there two independent variables? Default is FALSE.
#' @note If ANOVA assumptions are failed, set logIt = T, then reset logIt and check the residuals of the linear model using check=T
#' @return Returns letter code reults of ANOVA and also results of the statistical assumption tests
#' @export
ldStats = function(db,check=F, logIt=F,twoVar=F){
  if (twoVar==F){

    ##First function takes imputted wide format data and reorganizes it to long format
    numcol= ncol(db)
    mydp=data.frame()
    orgData= data.frame()
    orgData.names = NULL
    mean_cfu=NULL
    numrep=NULL
    ##Complicated mess I can parse, but this is the workhorse of this function
    #If something goes wrong double check that the format of the imported data is correct
    #This function also counts the number of replicates that are in each sample
    for (i in 1:numcol){
      db[,i]= as.numeric(db[,i])
      mydp = rbind(mydp,sd(db[,i],na.rm=T))
      mean_cfu = rbind(mean_cfu, mean(db[,i],na.rm=T))
      x=0
      for (j in 1:nrow(db)){
        orgData = rbind(orgData, db[j,i])
        orgData.names = rbind(orgData.names, names(db)[i])
        x=x+(1-as.integer(is.na(db[j,i])))
      }
      numrep = c(numrep,rep(x,each = j))
    }
    newDB = NULL
    orgData = cbind(orgData.names,orgData,numrep)
    names(orgData)= c("type","cfu","reps")
    #Removal of bad values
    orgData = subset(orgData, reps >1 ) #Removes data that only has one replicate
    orgData = delete.na(orgData,"greater") #'greater' specifies to remove any rows that contain an NA value
    rownames(orgData) = 1:nrow(orgData)

    #Makes an easy to read (and graph) version of the data
    mydp = cbind(names(db),mean_cfu,mydp)
    names(mydp)= c("type","cfu","stdev")

    #log the data before beginning
    if (logIt ==T){
      orgData$cfu = log10(orgData$cfu)
    }


    ##start of second function

    #Generates the output matrix
    results = matrix(NA,ncol =4,nrow =2)
    colnames(results) = c("ANOVA","Shapiro-Wilk","Levene","Bartlett")
    rownames(results) = c("P-value","Result")
    #Generating the ANOVA model
    anovaModel = aov(cfu ~ type, data=orgData) #the numeric variable must be first in this expression
    r = summary(anovaModel)
    results[1,1] = r[[1]]$`Pr(>F)`[1]
    orgData.lm = lm(cfu ~ type, data = orgData) #Using same model generates a linear model
    orgData.res = resid(orgData.lm) #rather residuals of linear model to see if their are any outliers
    results[1,2]=shapiro.test(orgData.res)[2]$p.value #Can be run with two replicates but runs under the assumption of 3. If you only have 2 replicates be careful
    results[1,3]=leveneTest(cfu ~ type, data=orgData)[3][1,1]
    results[1,4] = bartlett.test(cfu ~ type, data=orgData)[3]$p.value
    par(mfrow = c(2,2))
    plot(orgData.lm)
    par(mfrow = c(1,1))
    results[2,] = rep("FAIL",each = ncol(results))
    if (as.numeric(results[1,1])<=0.05){
      results[2,1]= "PASS"
    }
    for (i in 2:ncol(results)){
      if (as.numeric(results[1,i])>0.05){
        results[2,i]="PASS"
      }
    }


    if (check ==T){
      newDB$long = orgData
      newDB$lm = orgData.lm
      newDB$grouped = mydp
      newDB$rawstats = results
      return(newDB)
    }

    y=NULL #this was designed to be an output variable if you need to output both stats table and tukey grouping
    y$groups= HSD.test(anovaModel, alpha=0.05, "type", console=T)$groups#gives letter codes but no p values
    y$stats=results
    return(y)
  } else{
    ##This uses a different input format than above
    #Makes the assumption that each treatment type has at least 2 replicates
    #colVar = column based variable
    #reorganizing database to be usable
    colVar=t(db[1,2:ncol(db)]) #takes top row of database as a copy of which samples are var2- and which are mature
    rownames(db)=db[,1]
    db = t(db[2:nrow(db),2:ncol(db)])
    db <- as.data.frame(apply(db,FUN= as.numeric,MARGIN=c(1,2)))

    mydp= NULL
    orgData= data.frame()
    mean_cfu=NULL
    mydp.colVar = NULL
    newDB = NULL


    ##This block makes both the grouped and long format data
    for (i in 1:ncol(db)){
      mean_cfu = rbind(mean_cfu,mean(db[1:3,i],na.rm=T))
      mean_cfu = rbind(mean_cfu,mean(db[4:6,i],na.rm=T))
      mydp = rbind(mydp,sd(db[1:3,i],na.rm=T))
      mydp = rbind(mydp,sd(db[4:6,i],na.rm =T))
      mydp.colVar =rbind(mydp.colVar,colVar[1])
      mydp.colVar = rbind(mydp.colVar,colVar[4])
    }
    mydp = cbind(rep(colnames(db[,1:ncol(db)]),each =2),mean_cfu,mydp,mydp.colVar)
    colnames(mydp) = c("type","cfu","stdev","colVar")
    orgData.names = as.data.frame(rep(colnames(db),each = 6))
    orgData.colVar = as.data.frame(rep(colVar, times = ncol(db)))
    orgData = unlist(db)
    names(orgData)= rep_len(NA, length(orgData))
    orgData = cbind(orgData.names,orgData,orgData.colVar)
    colnames(orgData) = c("type","cfu","colVar")

    rm(orgData.names,colVar,orgData.colVar,mydp.colVar,i,mean_cfu) #clean-up on aisle 404

    if (logIt == T){
      orgData$cfu = log10(orgData$cfu)
    }

    #Removes any rows that contain any NA values
    orgData = delete.na(orgData,"greater") #'greater' specifies to remove any rows that contain an NA value

    ##ANOVA and assumption testing for young and old
    results = matrix(NA,ncol =4,nrow =2)
    colnames(results) = c("ANOVA","Shapiro-Wilk","Levene","Bartlett")
    rownames(results) = c("P-value","Result")


    #model set up
    mixed <- with(orgData, interaction(type, colVar))
    anovaModel = aov(cfu ~ mixed, data = orgData)
    summary(anovaModel)
    #test for normality of the residuals
    orgData.lm = lm(cfu ~ mixed, data = orgData)
    orgData.res = resid(orgData.lm)
    names(orgData.res)= orgData[,1]
    results[1,2] =shapiro.test(orgData.res)
    results[2,2]= results[1,2]>0.05 #pass if p>0.05


    #tests for homegeneity of variance, null hypothesis: homogeneity of variance, only need to pass one
    results[1,3] = leveneTest(cfu ~ mixed, data=orgData) #passes if p>0.05 thus variance if homogenous
    results[2,3] = results[1,3] >0.05
    results[1,4] = bartlett.test(cfu ~ mixed, data=orgData) #if you can run this run this but is more fickle
    results[2,4] = results[1,4]>0.05
    #POST-HOC TESTING
    #only if your data satisfy the assumptions and the ANOVA returned a significant result
    library(agricolae)
    TukeyHSD(anovaModel) #gives p values for pair-wise comparisons
    print(HSD.test(anovaModel, alpha = 0.05,trt ="mixed")) #gives letter codes but no p values
    anovaModel$model$mixed = make.unique(as.character(anovaModel$model$mixed))

    if (check ==T){
      newDB$long = orgData
      newDB$grouped = mydp
      newDB$lm = orgData.lm
      newDB$rawStats = results
      return(newDB)
    }
    y = NULL
    y$groups = HSD.test(anovaModel, alpha = 0.05,trt ="mixed")
    y$stats = results
    return(y)
  }
}
