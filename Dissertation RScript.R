#poisson model for the number of goals scored 
library('coda')
library('rjags')
library('dplyr')
library('reshape2')
library('knitr')
library('ggplot2')
library('writexl')
library('ggrepel')

#ENTERING THE DATA
#attacking data including goals, attempts on target and corners
#make sure data is sorted by StatsName, then by HomeTeamName, then by AwayTeamName
data<-read.csv(file.choose(),header=TRUE)

#GOALS MODELS  
goalsdat<-data[1063:1596,]  
str(goalsdat)
View(goalsdat)

#Player and team identification 
playersgoals<-unique(goalsdat$Player.Name)
teamsgoals<-unique(c(goalsdat$HomeTeamName,goalsdat$AwayTeamName))

#POISSON MODEL FOR NUMBER OF GOALS 
#BUILDING THE MODEL 
model="
model{
#DEFINING THE LIKELIHOOD AND RANDOM EFFECT MODEL FOR THE SCORING INTENSITY 
for (j in 1:ngames){
X[j]~dpois(eta[j]*tau[j])

#Predictive distribution for the number of goals scored 
X.pred[j]~dpois(eta[j]*tau[j])

#eta[j]
log(eta[j])<-delta[j]+tau[j]*(lambda1[Team[j]]-lambda2[Opp[j]])+home*kronecker[j]
}

#BASIC MODEL FOR THE HYPERPARAMETERS 
#prior on the home effect 
home~dnorm(0,0.0001)

#Trick to code the sum-to-zero constraint 
for (i in 1:nplayers){
  delta.star[i]~dnorm(m,s) #s is the precision i.e. 1/variance  
  delta[i]<-delta.star[i]-mean(delta.star[])
}

for (j in 1:nteams){
  lambda1.star[j]~dnorm(mu.lambda1,tau.lambda1)  #consider this the attack effect 
  lambda2.star[j]~dnorm(mu.lambda2,tau.lambda2)
  lambda1[j]<-lambda1.star[j]-mean(lambda1.star[])
  lambda2[j]<-lambda2.star[j]-mean(lambda2.star[])
}

#priors in the random effects 
m~dnorm(0,0.0001)
mu.lambda1~dnorm(0,0.0001)
mu.lambda2~dnorm(0,0.0001)

s~dgamma(0.1,0.01)
tau.lambda1~dgamma(0.01,0.01)
tau.lambda2~dgamma(0.01,0.01)
}
"
attach(goalsdat)
nplayers<-n_distinct(goalsdat[,4])
ngames<-n_distinct(goalsdat[,1])
nteams<-n_distinct(c(goalsdat$HomeTeamName,goalsdat$AwayTeamName))
datalist<-list('nplayers'=nplayers, 'ngames'=ngames, 'nteams'=nteams, 'X'=goalsdat[,10], 'tau'=FractionTime, 'kronecker'=Home.Away, 'Team'=Team.Index, 'Opp'=Opp.Index)

#FIT THE MODEL IN JAGS 

inits1<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=2000)
inits2<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=6000)
inits3<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=10000)
initial<-list(inits1,inits2,inits3)

output<-jags.model(textConnection(model),data=datalist,inits=initial,n.chains=3,n.adapt=5000)
update(output,n.iter=5000)
codaSamples2<-coda.samples(output,variable.names=c("delta","home","X.pred","lambda1","lambda2"),n.iter=10000)
results<-summary(codaSamples2)
#removing burnin 
results1<-summary(window(codaSamples2,15001))
results_matrix<-as.matrix(window(codaSamples2,15001)) #3*(10000-5000)=15000 rows 275 coloumns
df<-as.data.frame(results_matrix)

#EXTRACTING RESULTS TO CSV
resultsdf1<-data.frame(results1[1]) #converting the results to a data frame 
resultsdf2<-data.frame(results1[2])
write.csv(resultsdf1,file="goals_results.csv")
write.csv(resultsdf2,file="goals_results_2.csv")
write.csv(df,file="goals_noburnin.csv")

#TRACEPLOTS
plot(window(codaSamples2[,"delta[45]"],15001),main="delta[45]")
plot(window(codaSamples2[,"delta[82]"],15001),main="delta[82]")
#delta[45] - Ciro Immobile, delta[82] - Jonas Wind

#CHECK FOR CONVERGENCE 
geweke.diag(window(codaSamples2[,"delta[45]"],15001)) 
geweke.diag(window(codaSamples2[,"delta[82]"],15001)) 

gelman.diag(window(codaSamples2[,"delta[45]"],15001))
gelman.diag(window(codaSamples2[,"delta[82]"],15001))

#finding an average for the 3 chains 
mat<-mat.or.vec(5000,275)
for (i in 1:5000){
  for (j in 1:275){
    mat[i,j]=mean(df[i,j],df[i+5000,j],df[i+10000,j])
  }
}

#RMSE
newdf<-mat[,1:50] #50 coloumns, 5000 rows
matrix<-rbind(newdf,goalsdat$Value) #Select X.pred only --> exclude deltas, att, def and home
diff.square<-mat.or.vec(5000,50) 
for(i in 1:5000){
  for(j in 1:50){
    diff.square[i,j]=(matrix[i,j]-matrix[5001,j])^2
  }
} #Calculating the difference squared between each elt of coloumn j with the j^th entry of dat$Value
rmse<-numeric(5000)
for (i in 1:5000){
  rmse[i]=sqrt(sum(diff.square[i,])/50)
} #Calculating RMSE for each coloumn 
finalrmse<-mean(rmse)

#NEGATIVE BINOMIAL MODEL FOR NUMBER OF GOALS 
#BUILDING THE MODEL 
neg_bin_model="
model{
#DEFINING THE LIKELIHOOD AND RANDOM EFFECT MODEL FOR THE SCORING INTENSITY 
for (j in 1:ngames){
X[j]~dnegbin(p[j],p[j]*tau[j]*eta[j]/(1-p[j]))

#Predictive distribution for the number of goals scored 
X.pred[j]~dnegbin(p[j],p[j]*tau[j]*eta[j]/(1-p[j]))

#eta[j]
eta[j]<-exp(delta[j]+tau[j]*(lambda1[Team[j]]-lambda2[Opp[j]])+home*kronecker[j])
p[j]~dunif(0,1)
}

#BASIC MODEL FOR THE HYPERPARAMETERS 
#prior on the home effect 
home~dnorm(0,0.0001)

#Trick to code the sum-to-zero constraint 
for (i in 1:nplayers){
  delta.star[i]~dnorm(m,s) #s is the precision i.e. 1/variance  
  delta[i]<-delta.star[i]-mean(delta.star[])
}

for (j in 1:nteams){
  lambda1.star[j]~dnorm(mu.lambda1,tau.lambda1)  #consider this the attack effect 
  lambda2.star[j]~dnorm(mu.lambda2,tau.lambda2)
  lambda1[j]<-lambda1.star[j]-mean(lambda1.star[])
  lambda2[j]<-lambda2.star[j]-mean(lambda2.star[])
}

#priors in the random effects 
m~dnorm(0,0.0001)
mu.lambda1~dnorm(0,0.0001)
mu.lambda2~dnorm(0,0.0001)

s~dgamma(0.1,0.01)
tau.lambda1~dgamma(0.01,0.01)
tau.lambda2~dgamma(0.01,0.01)
}
"
attach(goalsdat)
nplayers<-n_distinct(goalsdat[,4])
ngames<-n_distinct(goalsdat[,1])
nteams<-n_distinct(c(goalsdat$HomeTeamName,goalsdat$AwayTeamName))
datalist<-list('nplayers'=nplayers, 'ngames'=ngames, 'nteams'=nteams, 'X'=goalsdat[,10], 'tau'=FractionTime, 'kronecker'=Home.Away, 'Team'=Team.Index, 'Opp'=Opp.Index)

#FIT THE MODEL IN JAGS 
inits1<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=2000)
inits2<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=6000)
inits3<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=10000)
initial<-list(inits1,inits2,inits3)

neg_bin_output<-jags.model(textConnection(neg_bin_model),inits=initial, data=datalist,n.chains=3,n.adapt=5000)
update(neg_bin_output,n.iter=5000)
neg_bin_codaSamples2<-coda.samples(neg_bin_output,variable.names=c("delta","home","X.pred","lambda1","lambda2"),n.iter=10000)
neg_bin_results<-summary(neg_bin_codaSamples2)
#removing burnin 
neg_bin_results1<-summary(window(neg_bin_codaSamples2,15001))
neg_bin_results_matrix<-as.matrix(window(neg_bin_codaSamples2,15001))
neg_bin_df<-as.data.frame(neg_bin_results_matrix)

#EXTRACTING TO CSV FORMAT IN EXCEL
neg_bin_resultsdf1<-data.frame(neg_bin_results1[1]) #converting the results to a data frame 
neg_bin_resultsdf2<-data.frame(neg_bin_results1[2])
write.csv(neg_bin_resultsdf1,file="neg_bin_goals_results.csv")
write.csv(neg_bin_resultsdf2,file="neg_bin_goals_quantiles.csv")
write.csv(neg_bin_df,file="neg_bin_Goals_noburnin.csv")

#TRACEPLOT
plot(window(neg_bin_codaSamples2[,"delta[45]"],15001),main="delta[45]")
plot(window(neg_bin_codaSamples2[,"delta[82]"],15001),main="delta[82]")

#CHECK FOR CONVERGENCE 
geweke.diag(window(neg_bin_codaSamples2[,"delta[45]"],15001))
geweke.diag(window(neg_bin_codaSamples2[,"delta[82]"],15001))

gelman.diag(window(neg_bin_codaSamples2[,"delta[45]"],15001))
gelman.diag(window(neg_bin_codaSamples2[,"delta[82]"],15001))

#finding an average for the 3 chains 
neg_bin_mat<-mat.or.vec(5000,275)
for (i in 1:5000){
  for (j in 1:275){
    neg_bin_mat[i,j]=mean(neg_bin_df[i,j],neg_bin_df[i+5000,j],neg_bin_df[i+10000,j])
  }
}

#RMSE
neg_bin_newdf<-neg_bin_df[,1:50] #50 coloumns, 3000 rows
neg_bin_matrix<-rbind(neg_bin_newdf,goalsdat$Value) #Select X.pred only --> exclude deltas, att, def and home
neg_bin_diff.square<-mat.or.vec(5000,50) 
for(i in 1:5000){
  for(j in 1:50){
    neg_bin_diff.square[i,j]=(neg_bin_matrix[i,j]-neg_bin_matrix[5001,j])^2
  }
}

#Calculating the difference squared between each elt of coloumn j with the j^th entry of dat$Value
neg_bin_rmse<-numeric(5000)
for (i in 1:5000){
  neg_bin_rmse[i]=sqrt(sum(neg_bin_diff.square[i,])/50)
} #Calculating RMSE for each coloumn 
neg_bin_finalrmse<-mean(neg_bin_rmse)

#ATTEMPTS ON TARGET MODELS 
attdat<-data[1:531,] 
str(attdat)
View(attdat)

#Player and team identification 
attplayers<-unique(attdat$Player.Name)
attteams<-unique(c(attdat$HomeTeamName,attdat$AwayTeamName))

#POISSON MODEL FOR NUMBER OF ATTEMPTS ON TARGET 
#BUILDING THE MODEL 
att_pois_model="
model{
#DEFINING THE LIKELIHOOD AND RANDOM EFFECT MODEL FOR THE SCORING INTENSITY 
for (j in 1:ngames){
X[j]~dpois(eta[j]*tau[j])

#Predictive distribution for the number of goals scored 
X.pred[j]~dpois(eta[j]*tau[j])

#eta[j]
log(eta[j])<-delta[j]+tau[j]*(lambda1[Team[j]]-lambda2[Opp[j]])+home*kronecker[j]
p[j]~dunif(0,1)
}

#BASIC MODEL FOR THE HYPERPARAMETERS 
#prior on the home effect 
home~dnorm(0,0.0001)

#Trick to code the sum-to-zero constraint 
for (i in 1:nplayers){
  delta.star[i]~dnorm(m,s) #s is the precision i.e. 1/variance  
  delta[i]<-delta.star[i]-mean(delta.star[])
}

for (j in 1:nteams){
  lambda1.star[j]~dnorm(mu.lambda1,tau.lambda1)  #consider this the attack effect 
  lambda2.star[j]~dnorm(mu.lambda2,tau.lambda2)
  lambda1[j]<-lambda1.star[j]-mean(lambda1.star[])
  lambda2[j]<-lambda2.star[j]-mean(lambda2.star[])
}

#priors in the random effects 
m~dnorm(0,0.0001)
mu.lambda1~dnorm(0,0.0001)
mu.lambda2~dnorm(0,0.0001)

s~dgamma(0.1,0.01)
tau.lambda1~dgamma(0.01,0.01)
tau.lambda2~dgamma(0.01,0.01)
}
"
attach(attdat)
nplayers<-n_distinct(attdat[,4])
ngames<-n_distinct(attdat[,1])
nteams<-n_distinct(c(attdat$HomeTeamName,attdat$AwayTeamName))
datalist<-list('nplayers'=nplayers, 'ngames'=ngames, 'nteams'=nteams, 'X'=attdat[,10], 'tau'=FractionTime, 'kronecker'=Home.Away, 'Team'=Team.Index, 'Opp'=Opp.Index)

#FIT THE MODEL IN JAGS 
inits1<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=2000)
inits2<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=6000)
inits3<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=10000)
initial<-list(inits1,inits2,inits3)

att_pois_output<-jags.model(textConnection(att_pois_model),data=datalist,inits=initial,n.chains=3,n.adapt=5000)
update(att_pois_output,n.iter=5000)
att_pois_codaSamples2<-coda.samples(att_pois_output,variable.names=c("delta","home","X.pred","lambda1","lambda2"),n.iter=10000)
att_pois_results<-summary(att_pois_codaSamples2)

#removing burnin
att_pois_results1<-summary(window(att_pois_codaSamples2,15001))
att_pois_results_matrix<-as.matrix(window(att_pois_codaSamples2,15001))
att_pois_df<-as.data.frame(att_pois_results_matrix)

#EXTRACTING RESULTS TO CSV
att_pois_resultsdf1<-data.frame(att_pois_results1[1]) #converting the results to a data frame 
att_pois_resultsdf2<-data.frame(att_pois_results1[2])
write.csv(att_pois_resultsdf1,file="attempts_pois_results.csv")
write.csv(att_pois_resultsdf2,file="attempts_pois_results2.csv")
write.csv(att_pois_df,file="attempts_noburnin.csv")

#TRACEPLOTS
plot(window(att_pois_codaSamples2[,"delta[45]"],15001),main="delta[45]")
plot(window(att_pois_codaSamples2[,"delta[82]"],15001),main="delta[82]")
#delta[45] - Ciro Immobile, delta[82] - Jonas Wind

#CHECK FOR CONVERGENCE 
geweke.diag(window(att_pois_codaSamples2[,"delta[45]"],15001))
geweke.diag(window(att_pois_codaSamples2[,"delta[82]"],15001))

gelman.diag(window(att_pois_codaSamples2[,"delta[45]"],15001))
gelman.diag(window(att_pois_codaSamples2[,"delta[82]"],15001))

#RMSE
#finding an average for the 3 chains 
att_pois_mat<-mat.or.vec(5000,275)
for (i in 1:5000){
  for (j in 1:275){
    att_pois_mat[i,j]=mean(att_pois_df[i,j],att_pois_df[i+5000,j],att_pois_df[i+10000,j])
  }
}

att_pois_newdf<-att_pois_mat[,1:50] #50 coloumns, 4500 rows
att_pois_matrix<-rbind(att_pois_newdf,attdat$Value) #Select X.pred only --> exclude deltas, att, def and home
att_pois_diff.square<-mat.or.vec(5000,50) 
for(i in 1:5000){
  for(j in 1:50){
    att_pois_diff.square[i,j]=(att_pois_matrix[i,j]-att_pois_matrix[5001,j])^2
  }
} #Calculating the difference squared between each elt of coloumn j with the j^th entry of dat$Value
att_pois_rmse<-numeric(5000)
for (i in 1:5000){
  att_pois_rmse[i]=sqrt(sum(att_pois_diff.square[i,])/50)
} #Calculating RMSE for each coloumn 
att_pois_finalrmse<-mean(att_pois_rmse)

#NEGATIVE BINOMIAL MODEL FOR NUMBER OF ATTEMPTS ON TARGET
#BUILDING THE MODEL 
att_neg_bin_model="
model{
#DEFINING THE LIKELIHOOD AND RANDOM EFFECT MODEL FOR THE SCORING INTENSITY 
for (j in 1:ngames){
X[j]~dnegbin(p[j],p[j]*tau[j]*eta[j]/(1-p[j]))

#Predictive distribution for the number of goals scored 
X.pred[j]~dnegbin(p[j],p[j]*tau[j]*eta[j]/(1-p[j]))

#eta[j]
eta[j]<-exp(delta[j]+tau[j]*(lambda1[Team[j]]-lambda2[Opp[j]])+home*kronecker[j])
p[j]~dunif(0,1)
}

#BASIC MODEL FOR THE HYPERPARAMETERS 
#prior on the home effect 
home~dnorm(0,0.0001)

#Trick to code the sum-to-zero constraint 
for (i in 1:nplayers){
  delta.star[i]~dnorm(m,s) #s is the precision i.e. 1/variance  
  delta[i]<-delta.star[i]-mean(delta.star[])
}

for (j in 1:nteams){
  lambda1.star[j]~dnorm(mu.lambda1,tau.lambda1)  #consider this the attack effect 
  lambda2.star[j]~dnorm(mu.lambda2,tau.lambda2)
  lambda1[j]<-lambda1.star[j]-mean(lambda1.star[])
  lambda2[j]<-lambda2.star[j]-mean(lambda2.star[])
}

#priors in the random effects 
m~dnorm(0,0.0001)
mu.lambda1~dnorm(0,0.0001)
mu.lambda2~dnorm(0,0.0001)

s~dgamma(0.1,0.01)
tau.lambda1~dgamma(0.01,0.01)
tau.lambda2~dgamma(0.01,0.01)
}
"
attach(dat)
nplayers<-n_distinct(attdat[,4])
ngames<-n_distinct(attdat[,1])
nteams<-n_distinct(c(attdat$HomeTeamName,attdat$AwayTeamName))
datalist<-list('nplayers'=nplayers, 'ngames'=ngames, 'nteams'=nteams, 'X'=attdat[,10], 'tau'=FractionTime, 'kronecker'=Home.Away, 'Team'=Team.Index, 'Opp'=Opp.Index)

#FIT THE MODEL IN JAGS 
inits1<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=2000)
inits2<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=6000)
inits3<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=10000)
initial<-list(inits1,inits2,inits3)

att_neg_bin_output<-jags.model(textConnection(att_neg_bin_model),inits=initial, data=datalist,n.chains=3,n.adapt=5000)
update(att_neg_bin_output,n.iter=5000)
att_neg_bin_codaSamples2<-coda.samples(att_neg_bin_output,variable.names=c("delta","home","X.pred","lambda1","lambda2"),n.iter=10000) 
att_neg_bin_results<-summary(att_neg_bin_codaSamples2)
#removing burnin
att_neg_bin_results1<-summary(window(att_neg_bin_codaSamples2,15001))
att_neg_bin_results_matrix<-as.matrix(window(att_neg_bin_codaSamples2,15001))
att_neg_bin_df<-as.data.frame(att_neg_bin_results_matrix) 

#EXTRACTING TO CSV FORMAT IN EXCEL
att_neg_bin_resultsdf1<-data.frame(att_neg_bin_results1[1]) #converting the results to a data frame 
att_neg_bin_resultsdf2<-data.frame(att_neg_bin_results1[2])
write.csv(att_neg_bin_resultsdf1,file="neg_bin_attempts_results.csv")
write.csv(att_neg_bin_resultsdf2,file="neg_bin_attempts_quantiles.csv")
write.csv(att_neg_bin_df,file="noburnin_neg_bin_attempts.csv")

#TRACEPLOTS
plot(att_neg_bin_codaSamples2[,"delta[45]"])
plot(att_neg_bin_codaSamples2[,"delta[82]"])

#CONVERGENCE
geweke.diag(window(att_neg_bin_codaSamples2[,"delta[45]"],15001))
geweke.diag(window(att_neg_bin_codaSamples2[,"delta[82]"],15001))

gelman.diag(window(att_neg_bin_codaSamples2[,"delta[45]"],15001))
gelman.diag(window(att_neg_bin_codaSamples2[,"delta[82]"],15001))

#finding an average for the 3 chains 
att_neg_bin_mat<-mat.or.vec(5000,275)
for (i in 1:5000){
  for (j in 1:275){
    att_neg_bin_mat[i,j]=mean(att_neg_bin_df[i,j],att_neg_bin_df[i+5000,j],att_neg_bin_df[i+10000,j])
  }
}

#GOODNESS OF FIT
#RMSE in R
att_neg_bin_newdf<-att_neg_bin_df[,1:50] #50 coloumns, 3000 rows
att_neg_bin_matrix<-rbind(att_neg_bin_newdf,attdat$Value) #Select X.pred only --> exclude deltas, att, def and home
att_neg_bin_diff.square<-mat.or.vec(5000,50) 
for(i in 1:5000){
  for(j in 1:50){
    att_neg_bin_diff.square[i,j]=(att_neg_bin_matrix[i,j]-att_neg_bin_matrix[5001,j])^2
  }
}

#Calculating the difference squared between each elt of coloumn j with the j^th entry of dat$Value
att_neg_bin_rmse<-numeric(5000)
for (i in 1:5000){
  att_neg_bin_rmse[i]=sqrt(sum(att_neg_bin_diff.square[i,])/50)
} #Calculating RMSE for each coloumn 
att_neg_bin_finalrmse<-mean(att_neg_bin_rmse)

#CORNERS MODELS 
cordat<-data[532:1062,] 
str(cordat)
View(cordat)

#Player and team identification 
playerscor<-unique(cordat$Player.Name)
teamscor<-unique(c(cordat$HomeTeamName,cordat$AwayTeamName))

#BUILDING THE MODEL 
cor_pois_model="
model{
#DEFINING THE LIKELIHOOD AND RANDOM EFFECT MODEL FOR THE SCORING INTENSITY 
for (j in 1:ngames){
X[j]~dpois(eta[j]*tau[j])

#Predictive distribution for the number of goals scored 
X.pred[j]~dpois(eta[j]*tau[j])

#eta[j]
log(eta[j])<-delta[j]+tau[j]*(lambda1[Team[j]]-lambda2[Opp[j]])+home*kronecker[j]
}

#BASIC MODEL FOR THE HYPERPARAMETERS 
#prior on the home effect 
home~dnorm(0,0.0001)

#Trick to code the sum-to-zero constraint 
for (i in 1:nplayers){
  delta.star[i]~dnorm(m,s) #s is the precision i.e. 1/variance  
  delta[i]<-delta.star[i]-mean(delta.star[])
}

for (j in 1:nteams){
  lambda1.star[j]~dnorm(mu.lambda1,tau.lambda1)  #consider this the attack effect 
  lambda2.star[j]~dnorm(mu.lambda2,tau.lambda2)
  lambda1[j]<-lambda1.star[j]-mean(lambda1.star[])
  lambda2[j]<-lambda2.star[j]-mean(lambda2.star[])
}

#priors in the random effects 
m~dnorm(0,0.0001)
mu.lambda1~dnorm(0,0.0001)
mu.lambda2~dnorm(0,0.0001)

s~dgamma(0.1,0.01)
tau.lambda1~dgamma(0.01,0.01)
tau.lambda2~dgamma(0.01,0.01)
}
"
attach(cordat)
nplayers<-n_distinct(cordat[,4])
ngames<-n_distinct(cordat[,1])
nteams<-n_distinct(c(cordat$HomeTeamName,cordat$AwayTeamName))
datalist<-list('nplayers'=nplayers, 'ngames'=ngames, 'nteams'=nteams, 'X'=cordat[,10], 'tau'=FractionTime, 'kronecker'=Home.Away, 'Team'=Team.Index, 'Opp'=Opp.Index)

#FIT THE MODEL IN JAGS 

inits1<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=2000)
inits2<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=6000)
inits3<-list(.RNG.name="base::Marsaglia-Multicarry",.RNG.seed=10000)
initial<-list(inits1,inits2,inits3)

cor_pois_output<-jags.model(textConnection(cor_pois_model),data=datalist,inits=initial,n.chains=3,n.adapt=5000)
update(cor_pois_output,n.iter=5000)
cor_pois_codaSamples2<-coda.samples(cor_pois_output,variable.names=c("delta","home","X.pred","lambda1","lambda2"),n.iter=10000) #set seed inside coda.samples
#thinning allows the user to specify how much the MCMC chain should be thinned out before storing them e.g. thin=10 would keep every
#10th value and discard the rest --> thin=1 keeps all the values. 
cor_pois_results<-summary(cor_pois_codaSamples2)
cor_pois_results1<-summary(window(cor_pois_codaSamples2,15001))
cor_pois_matrix<-as.matrix(window(cor_pois_codaSamples2,15001))
cor_pois_df<-as.data.frame(cor_pois_matrix)

#EXTRACTING TO CSV
cor_pois_resultsdf1<-data.frame(cor_pois_results1[1]) #converting the results to a data frame 
cor_pois_resultsdf2<-data.frame(cor_pois_results1[2])
write.csv(cor_pois_resultsdf1,file="corners_poisson_results.csv")
write.csv(cor_pois_resultsdf2,file="corners_poisson_quantiles.csv")
write.csv(cor_pois_df,file="corners_noburnin.csv")

#TRACEPLOTS
plot(window(cor_pois_codaSamples2[,"delta[45]"],15001),main="delta[45]")
plot(window(cor_pois_codaSamples2[,"delta[82]"],15001),main="delta[82]")
#delta[45] - Ciro Immobile, delta[82] - Jonas Wind

#CHECK FOR CONVERGENCE 
geweke.diag(window(cor_pois_codaSamples2[,"delta[45]"],15001))
geweke.diag(window(cor_pois_codaSamples2[,"delta[82]"],15001))

gelman.diag(window(cor_pois_codaSamples2[,"delta[45]"],15001))
gelman.diag(window(cor_pois_codaSamples2[,"delta[82]"],15001))

#finding an average for the 3 chains 
cor_pois_mat<-mat.or.vec(5000,275)
for (i in 1:5000){
  for (j in 1:275){
    cor_pois_mat[i,j]=mean(cor_pois_df[i,j],cor_pois_df[i+5000,j],cor_pois_df[i+10000,j])
  }
}

#RMSE
cor_pois_newdf<-cor_pois_mat[,1:50] #50 coloumns, 4500 rows
cor_pois_matrix<-rbind(cor_pois_newdf,cordat$Value) #Select X.pred only --> exclude deltas, att, def and home
cor_pois_diff.square<-mat.or.vec(5000,50) 
for(i in 1:5000){
  for(j in 1:50){
    cor_pois_diff.square[i,j]=(cor_pois_matrix[i,j]-cor_pois_matrix[5001,j])^2
  }
} #Calculating the difference squared between each elt of coloumn j with the j^th entry of dat$Value
cor_pois_rmse<-numeric(5000)
for (i in 1:5000){
  cor_pois_rmse[i]=sqrt(sum(cor_pois_diff.square[i,])/50)
} #Calculating RMSE for each coloumn 
cor_pois_finalrmse<-mean(cor_pois_rmse)

#CLUSTER ANALYSIS 
#packages
library('scatterplot3d')
library('factoextra')
library('NbClust')
set.seed(2023)

#first check that list of players is the same for all 3 events
data<-read.csv(file.choose(),header=TRUE)
goalsdat<-data[1063:1596,]
attdat<-data[1:531,]
cordat<-data[532:1062,]
playersgoals<-unique(goalsdat$Player.Name)
playersattempts<-unique(attdat$Player.Name)
playerscorners<-unique(cordat$Player.Name)
playernumber<-c(1:1:176)
players<-cbind(playernumber,playersgoals,playersattempts,playerscorners)

#average of all the chains minus burn in 
#goals
df1<-read.csv(file.choose(),header=TRUE)
df2<-df1[,-1]
goalsvec<-numeric(176)
goals_dataframe<-df2[,51:226]
for (i in 1:176){
  goalsvec[i]=mean(goals_dataframe[,i])
}

#attempts
attvec<-numeric(176)
att_pois_df1<-read.csv(file.choose(),header=TRUE)
att_pois_df2<-att_pois_df1[,-1]
attempts_dataframe<-att_pois_df2[,51:226]
for (i in 1:176){
  attvec[i]=mean(attempts_dataframe[,i])
}

#corners
corvec<-numeric(176)
cor_pois_df1<-read.csv(file.choose(),header=TRUE)
cor_pois_df2<-cor_pois_df1[,-1]
corners_dataframe<-cor_pois_df2[,51:226]
for (i in 1:176){
  corvec[i]=mean(corners_dataframe[,i])
}

#plot the 3d graph
s3d<-scatterplot3d(attvec,goalsvec,corvec)
s3d.coords<-s3d$xyz.convert(goalsvec,attvec,corvec)
text(s3d.coords$x,s3d$coords$z,labels=players[,1])

#find total number of goals, attempts and corners for each player in the entire tournament 
totalgoals<-numeric(176)
for (i in 1:176){
  totalgoals[i]=with(goalsdat,sum(goalsdat$Value[goalsdat$Player.Name==playersgoals[i]]))
}

totalattempts<-numeric(176)
for (i in 1:176){
  totalattempts[i]=with(attdat,sum(attdat$Value[attdat$Player.Name==playersattempts[i]]))
}

totalcorners<-numeric(176)
for (i in 1:176){
  totalcorners[i]=with(cordat,sum(cordat$Value[cordat$Player.Name==playerscorners[i]]))
}

#turn goalsplayers vector to data frame 
playersgoals_dataframe<-as.data.frame(playersgoals)
#add labels to points
colnames(goals_dataframe)<-playersgoals_dataframe[,1]

#Cluster Analysis 
#Create Data Frame
clusterdata1<-cbind(playersgoals,goalsvec,attvec,corvec) #the list of playersgoals is equivalent to that of att and cor
write.csv(clusterdata1,file="clusterdata.csv") #saved to csv and re-entered to be able to scale

clusterdata<-read.csv(file.choose(),header=TRUE) 
clusterdata2<-clusterdata[,-1] #removed first coloumn (not needed)
clusterdata3<-data.frame(clusterdata2[,1:4]) 

#standardize data to make it comparable
clusterdata4<-scale(clusterdata3[,2:4]) 

#K-MEANS CLUSTERING
#determining the optimal number of clusters
#Elbow method 
fviz_nbclust(clusterdata4,kmeans,method="wss") #5 clusters
#Silhouette method 
fviz_nbclust(clusterdata4,kmeans,method="silhouette") #2 clusters 
#Gap statistic method 
fviz_nbclust(clusterdata4,kmeans,nstart=100,method="gap_stat") #2 clusters 


set.seed(2023)
clusterdatakm<-kmeans(clusterdata4,5)
fviz_cluster(clusterdatakm,data=clusterdata4)
o=order(clusterdatakm$cluster)
dataframe10<-data.frame(clusterdata3[,1][o],clusterdatakm$cluster)
dataframe11<-dataframe10[order(dataframe10$clusterdatakm.cluster),]
write.csv(dataframe11,file='clustered.csv')

#extract the names of players in each cluster
dataframe11[dataframe11$clusterdatakm.cluster=='1',][,1]
dataframe11[dataframe11$clusterdatakm.cluster=='2',][,1]
dataframe11[dataframe11$clusterdatakm.cluster=='3',][,1]
dataframe11[dataframe11$clusterdatakm.cluster=='4',][,1]
dataframe11[dataframe11$clusterdatakm.cluster=='5',][,1]