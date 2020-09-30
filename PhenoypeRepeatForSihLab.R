## Code for SIHLABBERS / Orr SPiegel May 2016###

###### IS MY PEHNOTYPE REPEATABLE? A Randomization test  #########
# i've also asked Dave Harris on this test. 
# this is a simple test to see if having individual identity improves predictability (i.e., reduces STD between samples)
Randomizations=5000;#how many randomizations to do on that

# a helper function #####

# this function gets a dataframe with IDs and scores for each lizard, and calculate the STD for each lizard, 
#and the mean STD for across all lizards. can be applied to randomly assinged lizards names to generate random expectation  
#if individuals are repeatable, than STD of real lizards should be lower than the values for randomized data
STDcalculator <- function(Dat) {
  UniqNames=unique(Dat$Lizard_ID);#the lizards names included in this randomize dataset
  STDs= rep(NA, length(UniqNames))#a parameter to store the STD of each "lizard"
  for(CerNum in 1:length(UniqNames)){#loop on all lizard names to calc their STD
    STDs[CerNum]= sd(Dat$TrialScore[which(Dat$Lizard_ID==UniqNames[CerNum])])
  }#loop  
  #plot(density(STDs))
  return(mean(STDs))
}#function

### loading data #####
setwd('D:\\Dropbox\\R codes'); setwd('C:\\Users\\ors\\Dropbox\\R codes')
load ('dataforsihlab.rdata')
View(DatPhenotypScore); #a dataframe with lizard ID and phenotype score in a given assay, so each lizard appear 3 times

##  TrialsPhenotypScore
#finding the observed STD (mean over all lizards)
#DatPhenotypScore=as.data.frame(cbind(Data2$Lizard_ID,Data2$TrialsPhenotypScore ));colnames(DatPhenotypScore)=c('Lizard_ID','TrialScore')


#getting the observed value to compare with
ObsSTDPhenotypScore=STDcalculator(DatPhenotypScore) #sending real data to helper function

##### 
#finding the null distribution- in a loop with N=Randomizations 
rndSTDPhenotypScore=rep(NA, Randomizations) #here i store the result for each randomization
for(RndNum in 1:Randomizations){#loop on Randomizations
  Dat=as.data.frame(cbind(sample(DatPhenotypScore$Lizard_ID),DatPhenotypScore$TrialScore ))
  colnames(Dat)=c('Lizard_ID','TrialScore')
  rndSTDPhenotypScore[RndNum]=STDcalculator(Dat)
}#loop on Randomizations

  #looking on results
  plot(density(rndSTDPhenotypScore), xlab=('mean STD'),
       main=paste('randomization test: null distribution. Observed value=',round(ObsSTDPhenotypScore,digit=2)));
  abline(v=ObsSTDPhenotypScore, col=2); #redline for observed value

print(paste('the observed mean STD was larger than a prop of ',
            sum(rndSTDPhenotypScore<=ObsSTDPhenotypScore)/Randomizations ,'of the randomizations, out of ',Randomizations))
##this is the one sided p value for our test

##  another validation : #####
  # different means for each individual: DatPhenotypScore$TrialScore =DatPhenotypScore$TrialScore+ rep(1:60,each=3)
  # center individuals around zero so it is the same mean for all individuals: 
#   UniqNames=unique(DatPhenotypScore$Lizard_ID);#the lizards names included in this randomize dataset
#   for(CerNum in 1:length(UniqNames)){#loop on all lizard names to calc their STD
#     indx=which(DatPhenotypScore$Lizard_ID==UniqNames[CerNum])
#     DatPhenotypScore$TrialScore[indx]= DatPhenotypScore$TrialScore[indx]-mean(DatPhenotypScore$TrialScore[indx])
#   }#loop
 
  
## simulating data with normal distribution ####
  #generating pseudo data
  Nindiv=60;
  Mtrials=3;
  
  SpreadWithinIndiv=10;#how spread the results for each individual
  DiffAmongIndiv=3; #how different the individuals are
  DatPhenotypScore$Lizard_ID=rep(1:Nindiv,each=Mtrials);
  
  indviScore=rep(rnorm(Nindiv,mean= 0,sd=DiffAmongIndiv),each=Mtrials) #different value for each individual same for all his trials
  DatPhenotypScore$TrialScore=indviScore +rnorm(Nindiv*Mtrials,mean= 0,sd=SpreadWithinIndiv);#different values for the trials
    #plotting the data
  plot(x=DatPhenotypScore$Lizard_ID,
       y=DatPhenotypScore$TrialScore,typ='p',
       xlab='individual #',ylab='trial score', 
       main=paste('SpreadWithinIndiv=',SpreadWithinIndiv,'DiffAmongIndiv=',DiffAmongIndiv, 'repeatb=?',DiffAmongIndiv^2/(DiffAmongIndiv^2+SpreadWithinIndiv^2) ))

  #now run the code above again! lines 33-47  ###
 
  
## simulating data with uniform distributrion ####
  #generating pseudo data
  SpreadWithinIndiv=3;#how spread the results for each individual
  DiffAmongIndiv=0.5; #how different the individuals are
  DatPhenotypScore$Lizard_ID=rep(1:Nindiv,each=Mtrials);
  #DatPhenotypScore$TrialScore=10+rep(1:Nindiv,each=Mtrials)*DiffAmongIndiv+runif(Nindiv*Mtrials,-1,1)*Spread
  DatPhenotypScore$TrialScore=rep(runif(Nindiv,-0.5,0.5)*DiffAmongIndiv,each=Mtrials)+ runif(Nindiv*Mtrials,-0.5,0.5)*SpreadWithinIndiv
  #plotting the data
  plot(x=DatPhenotypScore$Lizard_ID,
       y=DatPhenotypScore$TrialScore,typ='p',
       xlab='individual #',ylab='trial score', 
       main=paste('SpreadWithinIndiv=',SpreadWithinIndiv,'DiffAmongIndiv=',DiffAmongIndiv))
  
  #now run the code above again! lines 33-47  ###
  
  
  #### calculate repeatability #####
  
  