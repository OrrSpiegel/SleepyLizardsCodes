
#this code reads the behavioral scores for 2015 agression trails, 
#and looks for the best combination of behaviors for repeatabilty

### STEP1.1: some parameters to choose ####
Randomizations=5000;#how many randomizations to do on that repeatability test
PCweightMinVal=0.25;#if using David's approach - only some predicotrs to calcscores in PC1, and only thier direction
IndxBehaveToPCA=17:32; #all==17:32 indices of all behave in the long form (17:32), and in the agreggate dataframe 3:18
#IndxBehaveToPCA=c(1,6,11,14)+16; #to work only on: "NumStrikes" "NoRxn" "ShelterLatency" "FleeLatency": was 0.2 for Repeatability
WorkOnBehavsWithUnivarRepeatOnly=2;#change this to 0- use all behaviors
                                                  #1 to work only on those behaiovrs that i found to be with a sig repeatability using my randomizaitob test
                                                  #2 work only on Number of strikes,Backed Off (corrected),No Reaction (corrected),Flee Rating

CorrectForNumOfStrikes=1 #since many of the behaviors are a rate, we need to correct for how many strikes were actually done
#11 behaviors showing high repeatbility:"NumStrikes" "AlertPosture"  "BackedOff" "TopHead" "NoRxn" "TongueFlick" "Gaping" "ShelterLatency" "FleeRating" "FleeLatency"
# marginal repeatability : "Wince"
# not repeatable: Biting , Hissing (becomes repeatable after correction) Flee_Y.N FleeLatency EndDistance FurthestDistance 
if (CorrectForNumOfStrikes==0 & WorkOnBehavsWithUnivarRepeatOnly==1) IndxBehaveToPCA=IndxBehaveToPCA [! IndxBehaveToPCA %in% (16+c(3,9,10,12,14,15,16))]# here i drop out those behaviors that were not repeatable at all  in the randomization. looks good, it give PC1var of 0.398 and repeatalby 0.31
if (CorrectForNumOfStrikes==1 & WorkOnBehavsWithUnivarRepeatOnly==1) IndxBehaveToPCA=IndxBehaveToPCA [! IndxBehaveToPCA %in% (16+c(3,9,   12,14,15,16))]# here i drop out those behaviors that were not repeatable at all  in the randomization. looks good, it give PC1var of 0.364 and repeatalby 0.31
if (                            WorkOnBehavsWithUnivarRepeatOnly==2) IndxBehaveToPCA=c(1,   4,6,    13)+16# here i drop out everything but a behaviors listed for option 2, all corrected
if (                            WorkOnBehavsWithUnivarRepeatOnly==3) IndxBehaveToPCA=c(1,3, 4,6, 11,13)+16# here i drop out everything but a behaviors listed for option 3, all corrected


#loading packages stats lme4
require(stats)
require(lme4)


### STEP1.2: helper functions for the randomizaiton test ####
#this function gets a dataframe with two cols- IDs (or pseudo IDs) and a trait to estimate.called by the operator function
STDcalculator <- function(Dat) {
  UniqNames=unique(Dat$LizardID);#the lizards names included in this randomize dataset
  STDs= rep(NA, length(UniqNames))#a parameter to store the STD of each "lizard"
  for(CerNum in 1:length(UniqNames)){#loop on all lizard names to calc their STD
    STDs[CerNum]= sd(Dat$TrialScore[which(Dat$LizardID==UniqNames[CerNum])])
  }#loop  
  #plot(density(STDs))
  return(mean(STDs, na.rm=TRUE))
}#function

#this function is gets the ids and a trait, calls the STDcalculator for the real data and then for randomized data
RndmztnTestOperator<- function(DataFrRndmztn){
  CurrTraitName=names(DataFrRndmztn)[2];names(DataFrRndmztn)[2]='TrialScore' #just logging the name and converting to a generic name for column
  ObsSTDPhenotypScore=STDcalculator(DataFrRndmztn) #sending real data to helper function
  
  #finding the null distribution- in a loop with N=Randomizations 
  rndSTDPhenotypScore=rep(NA, Randomizations) #here i store the result for each randomization
  for(RndNum in 1:Randomizations){#loop on Randomizations
    Dat=as.data.frame(cbind(sample(DataFrRndmztn$LizardID),DataFrRndmztn$TrialScore ))
    colnames(Dat)=c('LizardID','TrialScore')
    rndSTDPhenotypScore[RndNum]=STDcalculator(Dat)
  }#loop on Randomizations
  
  #looking on results
  PV=sum(rndSTDPhenotypScore<=ObsSTDPhenotypScore)/Randomizations
  plot(density(rndSTDPhenotypScore), xlab=('mean STD'),main=(paste(CurBehv-16, CurrTraitName,PV))  );
  abline(v=ObsSTDPhenotypScore, col=2); #redline for observed value
  #main=paste(CurrTraitName,': randomization test: null distribution. Observed value=',round(ObsSTDPhenotypScore,digit=2))
  print(paste(CurrTraitName,':the observed mean STD was larger than a prop of ',
              PV ,'of the randomizations, out of ',Randomizations))
  rm(CurrTraitName,PV)
  
}#end of RndmztnTestOperator
  


### STEP1.3: read csv, set some relevant columns as factor ####
#setwd('C:\\Users\\ors\\Dropbox\\R codes\\BehavioralAssayAnalyses')#lab pc:
#setwd('D:\\Dropbox\\R codes\\BehavioralAssayAnalyses')#home pc
setwd('D:\\OrrS2\\Box Sync\\R codes\\Lizards\\BehavioralAssayAnalyses')
AggrsScors2015= read.csv('aggressionScoresLongform2015.csv', header = TRUE) 

#AggrsScors2015$LizardID=as.factor(AggrsScors2015$LizardID)
AggrsScors2015$TestedWith=as.factor(AggrsScors2015$TestedWith)
AggrsScors2015$CodedBy=as.factor(AggrsScors2015$CodedBy)
AggrsScors2015$TestArena=as.factor(AggrsScors2015$TestArena)
#AggrsScors2015$Flee_Y.N=as.factor(AggrsScors2015$Flee_Y.N)


### STEP1.4: CorrectForNumOfStrikes and removing a test that failed  (what about lizards with one trial who has no GPS?) ####
AggrsScors2015 <- AggrsScors2015[-c(which( AggrsScors2015$LizardID== 41221 & AggrsScors2015$Week== 5)), ]
#liz id: 4247 (no gps)  7618 (29 days) 10167 (13 days) 14111 (39 days with GPS) 41227 (25 days of GPS data)
#AggrsScors2015 <- AggrsScors2015[-which( AggrsScors2015$LizardID==c( 4247) ), ];
str(AggrsScors2015)

##CorrectForNumOfStrikes- behaviors corrected for number of strikes
if (CorrectForNumOfStrikes==1){
  AggrsScors2015$AlertPosture=AggrsScors2015$AlertPosture/AggrsScors2015$NumStrikes
  AggrsScors2015$Wince       =AggrsScors2015$Wince/       AggrsScors2015$NumStrikes
  AggrsScors2015$BackedOff   =AggrsScors2015$BackedOff/   AggrsScors2015$NumStrikes
  AggrsScors2015$TopHead     =AggrsScors2015$TopHead/     AggrsScors2015$NumStrikes
  AggrsScors2015$NoRxn       =AggrsScors2015$NoRxn/       AggrsScors2015$NumStrikes
  AggrsScors2015$TongueFlick =AggrsScors2015$TongueFlick/ AggrsScors2015$NumStrikes
  AggrsScors2015$Gaping      =AggrsScors2015$Gaping/      AggrsScors2015$NumStrikes
  AggrsScors2015$Biting      =AggrsScors2015$Biting/      AggrsScors2015$NumStrikes
  AggrsScors2015$Hissing     =AggrsScors2015$Hissing/     AggrsScors2015$NumStrikes
}


### STEP2: for each behavior- is its repeatble across trials of the same indivdual? #### 
#sending the different traits to the operator function
DataFrRndmztn=as.data.frame(AggrsScors2015$LizardID);names(DataFrRndmztn)='LizardID'
par(mfrow=c(3,3))

for (CurBehv in 17:32 ){
  names(AggrsScors2015)[CurBehv]
  DataFrRndmztn[,2]=AggrsScors2015[,CurBehv]; 
  names(DataFrRndmztn)[2]=names(AggrsScors2015)[CurBehv];
  RndmztnTestOperator(DataFrRndmztn)
  DataFrRndmztn[,2]=NULL
}
par(mfrow=c(1,1));
#DataFrRndmztn$NumStrikes=AggrsScors2015$NumStrikes;


### STEP3: building the short form with average for each lizard ####
Indiv=unique(AggrsScors2015$LizardID)
AggrsScors2015Aggrgt=data.frame(Indiv);colnames(AggrsScors2015Aggrgt)='LizardID'
AggrsScors2015Aggrgt$NumTrials=sapply( Indiv, function(i){length(AggrsScors2015$LizardID[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$NumStrikes=sapply( Indiv, function(i){mean(AggrsScors2015$NumStrikes[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$AlertPosture=sapply( Indiv, function(i){mean(AggrsScors2015$AlertPosture[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$Wince=sapply( Indiv, function(i){mean(AggrsScors2015$Wince[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$BackedOff=sapply( Indiv, function(i){mean(AggrsScors2015$BackedOff[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$TopHead=sapply( Indiv, function(i){mean(AggrsScors2015$TopHead[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$NoRxn=sapply( Indiv, function(i){mean(AggrsScors2015$NoRxn[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$TongueFlick=sapply( Indiv, function(i){mean(AggrsScors2015$TongueFlick[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$Gaping=sapply( Indiv, function(i){mean(AggrsScors2015$Gaping[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$Biting=sapply( Indiv, function(i){mean(AggrsScors2015$Biting[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$Hissing=sapply( Indiv, function(i){mean(AggrsScors2015$Hissing[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$ShelterLatency=sapply( Indiv, function(i){mean(AggrsScors2015$ShelterLatency[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$Flee_Y.N=sapply( Indiv, function(i){mean(AggrsScors2015$Flee_Y.N[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$FleeRating=sapply( Indiv, function(i){mean(AggrsScors2015$FleeRating[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$FleeLatency=sapply( Indiv, function(i){mean(AggrsScors2015$FleeLatency[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$EndDistance=sapply( Indiv, function(i){mean(AggrsScors2015$EndDistance[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$FurthestDistance=sapply( Indiv, function(i){mean(AggrsScors2015$FurthestDistance[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$NewAGGRating=sapply( Indiv, function(i){mean(AggrsScors2015$NewAGGRating[c(which(AggrsScors2015$LizardID==i))])})
str(AggrsScors2015Aggrgt)

print(paste('lizards with 1 trials are', AggrsScors2015Aggrgt$LizardID[which(AggrsScors2015Aggrgt$NumTrials==1)]))   
AggrsScors2015Aggrgt2=as.data.frame(scale(AggrsScors2015Aggrgt[3:18]))


### STEP4.1: preps for PCA ####
#http://www.r-bloggers.com/computing-and-visualizing-pca-in-r/
#https://cran.r-project.org/web/packages/HSAUR/vignettes/Ch_principal_components_analysis.pdf

IndxBehaveToPCA2=IndxBehaveToPCA-14; #indices in the second dataframe
#just making sure i have all behaviors:
rbind(names(AggrsScors2015)[IndxBehaveToPCA],names(AggrsScors2015Aggrgt)[IndxBehaveToPCA2])


### STEP4.2: PCA analysis and plotting ####
#on Agreggated Data (1 val per lizard), includes also scaling and centering the data. log transform (had no effect i could see on the first 15 traits pca)?
PCA_Aggregat=prcomp(AggrsScors2015Aggrgt[,IndxBehaveToPCA2],center = TRUE, scale. = TRUE) ;

# plotting & printing resutls of the PCA
print(PCA_Aggregat)[1]
#the Standard deviations are the loading of the PCs, the rotation is the weights. 
plot(PCA_Aggregat,main='variance explained by PC axes');#plot(PCA_Aggregat, type = "l")
biplot(PCA_Aggregat)
summary(PCA_Aggregat)
ProVarExplPC1=summary(PCA_Aggregat)$importance[2,1]
print(paste('the first PC explains',ProVarExplPC1, 'of the variance'))


### STEP4.3: calculating PC1 the loading (wieghts) for the two methods ####
#for each behavior for the PC1 scores with David's approach- just + -1 or 0, only larger values
PC1weigthAll   <- PCA_Aggregat$rotation[,1]  #the weights for PCA 1
PC1weigthDavid <- PCA_Aggregat$rotation[,1] #the weights for PCA 1
PC1weigth2=PC1weigthDavid;
PC1weigth2[abs(PC1weigth2) < PCweightMinVal]=0;#working only on PCweight bigger than threshold
PC1weigthDavid=PC1weigthDavid/abs(PC1weigthDavid);#ignoring weights, just direction (plus minus 1)
PC1weigthDavid[PC1weigth2==0]=0;
rm(PC1weigth2)
hist(PC1weigthAll,breaks=10,main=(paste('loadings of PC1, David Ss method will use those > |', PCweightMinVal,'|')))

PCA_AggregatD=PCA_Aggregat;PCA_AggregatD$rotation[,1]=PC1weigthDavid; #replacing loadings with David's loading

AggrsScors2015Aggrgt$PC1scoreAll=predict(PCA_Aggregat)[,1] #All loading. this is what we need as the Aggresivness of the lizard!!!
AggrsScors2015Aggrgt$PC1scoreDavids=predict(PCA_AggregatD,newdata=AggrsScors2015Aggrgt[,IndxBehaveToPCA2])[,1] #Davids Loading. this is what we need a the Aggresivness of the lizard!!!


### STEP4.4: Scaling the longform data, and calculating the score for the trial from PC1 ####
#now scaling the longform data (i.e., the relevant cols in AggrsScors2015)
#what scaling param the PCA used?
centerUsed <- PCA_Aggregat$center;scaleUsed <- PCA_Aggregat$scale
#scaling appropriately also the long form data"
ScaledAggrsScors2015=as.data.frame(scale(as.matrix(AggrsScors2015[IndxBehaveToPCA]), center = centerUsed, scale = scaleUsed))

summary(ScaledAggrsScors2015) #the scaling for the Aggregated data are not a good fit for long form
print('the summary shows that the scaling is not perfect since it was done for Aggregt scores and now applied to trial data')

#adding the scores from the PCA to the trials in long form dataset- so we can calc repeatability
AggrsScors2015$PCA1stScoreByAllLoading=   drop(as.matrix(ScaledAggrsScors2015)  %*% PC1weigthAll)
AggrsScors2015$PCA1stScoreByDavidsLoading=drop(as.matrix(ScaledAggrsScors2015)  %*% PC1weigthDavid)
#now the data frame has a colum for each lizard&trial for the PC1 scores. one for each method


### STEP4.5: just validation, checking correlation between methods and ####
#just for curiosity are the two methods correlated? yes. 
CorrelationInScoresBtwnApproaches=cor(AggrsScors2015$PCA1stScoreByDavidsLoading,AggrsScors2015$PCA1stScoreByAllLoading)
print(paste('now, each trail has a score, based on PC1, using some loading (Davids method)',
'and the alternative of full loadings. the corelation between methods:', round(CorrelationInScoresBtwnApproaches,digit=3)))
cor.test(AggrsScors2015$PCA1stScoreByDavidsLoading,AggrsScors2015$PCA1stScoreByAllLoading)

par(mfrow=c(1,2))
plot(x=AggrsScors2015$PCA1stScoreByDavidsLoading,y=AggrsScors2015$PCA1stScoreByAllLoading,
     typ='p',main=('cor between menthods by trial'))
plot(x=AggrsScors2015Aggrgt$PC1scoreDavids,y=AggrsScors2015Aggrgt$PC1scoreAll,
     typ='p',main=('cor between menthods for lizard'))
par(mfrow=c(1,1))

#trying with another method to get the scores- gives the same result- remeber NOT to scale before- it scales automaticly
AggrsScors2015$PCA1stScoreByAllLoading2   =predict(PCA_Aggregat,  newdata=AggrsScors2015[IndxBehaveToPCA])[,1]
AggrsScors2015$PCA1stScoreByDavidsLoading2=predict(PCA_AggregatD, newdata=AggrsScors2015[IndxBehaveToPCA])[,1] 

plot(x=AggrsScors2015$PCA1stScoreByDavidsLoading2   ,y=AggrsScors2015$PCA1stScoreByAllLoading2,
     typ='p',main=('cor between menthods for trials- a differenc calc same result'))
     
#CorrelationInScoresBtwnApproaches=cor(AggrsScors2015$PCA1stScoreByDavidsLoading2,AggrsScors2015$PCA1stScoreByAllLoading2)
#cor.test(AggrsScors2015$PCA1stScoreByDavidsLoading2,AggrsScors2015$PCA1stScoreByAllLoading2)
#cor(AggrsScors2015$PCA1stScoreByAllLoading2,AggrsScors2015$PCA1stScoreByAllLoading)
#cor(AggrsScors2015$PCA1stScoreByDavidsLoading2  ,AggrsScors2015$PCA1stScoreByDavidsLoading2)





#just for validation will i get the same reults with another calculation? YES:
#PedictedPCAScores=drop(scale(as.matrix(AggrsScors2015Aggrgt[IndxBehaveToPCA2]), center = centerUsed, scale = scaleUsed)
#                       %*% PCA_Aggregat$rotation[,1])
#almost the same as for the original data: PCA_Aggregat$x[,1];cor(PCA_Aggregat$x[,1],PedictedPCAScores)
#almost the same as for predicting resuts for the same data: predict(PCA_Aggregat)[,1] 
# PedictedPCAScores=drop(scale(as.matrix(AggrsScors2015[17:32]), center = centerUsed, scale = scaleUsed)
#                             %*% PCA_Aggregat$rotation[,1])
AggrsScors2015Aggrgt$PCA1stScoreByAllLoading   =sapply( Indiv, function(i){mean(AggrsScors2015$PCA1stScoreByAllLoading[c(which(AggrsScors2015$LizardID==i))])})
AggrsScors2015Aggrgt$PCA1stScoreByDavidsLoading=sapply( Indiv, function(i){mean(AggrsScors2015$PCA1stScoreByDavidsLoading[c(which(AggrsScors2015$LizardID==i))])})
#hist(AggrsScors2015Aggrgt$PCA1stScoreByDavidsLoading)
#hist(AggrsScors2015Aggrgt$PCA1stScoreByAllLoading)

par(mfrow=c(1,2))
plot(x=AggrsScors2015Aggrgt$PC1scoreDavids,y=AggrsScors2015Aggrgt$PCA1stScoreByDavidsLoading,
     typ='p',main=('cor between calculation within Davids loading menthod for lizard'))
plot(x=AggrsScors2015Aggrgt$PC1scoreAll,y=AggrsScors2015Aggrgt$PCA1stScoreByAllLoading ,
     typ='p',main=('cor between calculations within the All loading menthod for lizard'))
par(mfrow=c(1,1))


### STEP5.1: are scores repeatable? using the a randomization test ####
### aplying randomization for the two scores- are repatable?
par(mfrow=c(1,2))
DataFrRndmztn=AggrsScors2015[c("LizardID","PCA1stScoreByAllLoading")]
RndmztnTestOperator(DataFrRndmztn)
DataFrRndmztn=AggrsScors2015[c("LizardID","PCA1stScoreByDavidsLoading")]
RndmztnTestOperator(DataFrRndmztn)
par(mfrow=c(1,1))
rm(RndmztnTestOperator)


### STEP5.2: GLMM for trial scores to see the repeatbility for both methods ####
#modeling the PC1 scores for each trials, is it repeatble for individuals?
Model1AllLoading= lmer(PCA1stScoreByAllLoading ~ 0 +   (1|LizardID),data=AggrsScors2015,REML = T)
summary(Model1AllLoading);

Model1DavidsLoading= lmer(PCA1stScoreByDavidsLoading ~ 0 +   (1|LizardID),data=AggrsScors2015,REML = T)
summary(Model1DavidsLoading);

#getting repeatability estimate and number of observations and individuals for the best model and the second Best model
vars=as.data.frame(VarCorr(Model1AllLoading))$vcov;  
LizRptbltModel1AllLoading=round(vars[1]/(vars[1]+vars[2]),digit=2)#

vars=as.data.frame(VarCorr(Model1DavidsLoading))$vcov;  
LizRptbltModel1DavidsLoading=round(vars[1]/(vars[1]+vars[2]),digit=2)#

Nobs=nobs(Model1DavidsLoading)
Nindiv=dim(coef(Model1DavidsLoading)$LizardID)[1]                        

print(paste('Repeatability estimate from GLMMs are ',LizRptbltModel1AllLoading,' (All loading) and ',
            LizRptbltModel1DavidsLoading,'(davids method). for ', Nindiv,'lizards and ',Nobs,'trials'))

print('the Aggressivenes score for the lizards, based on PCA1 with the two methods are stored in ')
print('AggrsScors2015Aggrgt$:  ')
print('(All loading) $PCA1stScoreByAllLoading')
print('(Davids approach): $PCA1stScoreByDavidsLoading')


### Step 6 agreement between the scoring and JK video rating scores ####
print('by each test independently'); cor.test(AggrsScors2015$NewAGGRating,AggrsScors2015$PCA1stScoreByAllLoading)
print('for each lizard aggregated'); cor.test(AggrsScors2015Aggrgt$NewAGGRating,AggrsScors2015Aggrgt$PCA1stScoreByAllLoading)
par(mfrow=c(1,2))
plot(x=AggrsScors2015$NewAGGRating,y=AggrsScors2015$PCA1stScoreByAllLoading,
     typ='p',main=('cor between methods by trial'))
plot(x=AggrsScors2015Aggrgt$NewAGGRating,y=AggrsScors2015Aggrgt$PCA1stScoreByAllLoading,
     typ='p',main=('cor between menthods for lizard'), xlab='JK video rating',ylab='Agg scores from agg PCA')
abline(lm(AggrsScors2015Aggrgt$PCA1stScoreByAllLoading ~ AggrsScors2015Aggrgt$NewAGGRating),col='red')

par(mfrow=c(1,1))

Model1JKratingVideoAGG= lmer(NewAGGRating ~ 0 +   (1|LizardID),data=AggrsScors2015,REML = T)
summary(Model1JKratingVideoAGG);

#getting repeatability estimate and number of observations and individuals for the best model 
vars=as.data.frame(VarCorr(Model1JKratingVideoAGG))$vcov;  
LizRptbltModel1JKratingVideoAGG=round(vars[1]/(vars[1]+vars[2]),digit=2)#it doesnt make sense, whats wrong with my formula?!!!

AggrsScors2015$LizardAsCer=NA
for (i in 1:length(AggrsScors2015$LizardID)){
  AggrsScors2015$LizardAsCer[i]=which(levels(as.factor(AggrsScors2015$LizardID))==AggrsScors2015$LizardID[i])
}
plot(x=AggrsScors2015$LizardAsCer ,y=AggrsScors2015$NewAGGRating,
     typ='p',main=('JKs video rating by lizards'),col='blue',pch=16,xlab='Lizards as cerial id')
abline(v=c(1:73))


#### Step 7 intraclass correlation coefficient for JKs score (repeatability using 1-way ANOVA methods) #####
#which lizards had 3 trials, creating a wide form, each trial is a column, each row i a lizard
Liz=unique(AggrsScors2015$LizardID)
WideForm=NULL;Rowname=NULL
for(L in 1:length(Liz)){
   JKsRating=AggrsScors2015$NewAGGRating[AggrsScors2015$LizardID==Liz[L]]
   #print(JKsRating)
   if (length(JKsRating)==3) {WideForm=rbind(WideForm,JKsRating);Rowname=rbind(Rowname,Liz[L])}
   if (length(JKsRating)==2) {WideForm=rbind(WideForm,c(JKsRating,NA));Rowname=rbind(Rowname,Liz[L])} #comment this line to work only on liz with three trials
   if (length(JKsRating)==1) {WideForm=rbind(WideForm,c(JKsRating,NA,NA));Rowname=rbind(Rowname,Liz[L])}  #comment this line to work only on liz with three trials
}#loop on lizards
rownames(WideForm)=Rowname
library(psych)
ICC(WideForm, missing=FALSE, alpha=.05)
#Shrout and Fleiss (1979) consider six cases of reliability of ratings done by k raters on n targets. 
#ICC1: Each target is rated by a different judge and the judges are selected at random. (This is a one-way ANOVA fixed effects model and is found by (MSB- MSW)/(MSB+ (nr-1)*MSW))
#ICC2: A random sample of k judges rate each target. The measure is one of absolute agreement in the ratings. Found as (MSB- MSE)/(MSB + (nr-1)*MSE + nr*(MSJ-MSE)/nc)
#ICC3: A fixed set of k judges rate each target. There is no generalization to a larger population of judges. (MSB - MSE)/(MSB+ (nr-1)*MSE)
#Then, for each of these cases, is reliability to be estimated for a single rating or for the average of k ratings? (The 1 rating case is equivalent to the average intercorrelation, the k rating case to the Spearman Brown adjusted reliability.)
#ICC1 is sensitive to differences in means between raters and is a measure of absolute agreement.
#ICC2 and ICC3 remove mean differences between judges, but are sensitive to interactions of raters by judges. The difference between ICC2 and ICC3 is whether raters are seen as fixed or random effects.
#ICC1k, ICC2k, ICC3K reflect the means of k raters.

