# Function to estimate power for linear regresison with binary predictor

#beta: your coefficient from regression
#Y.sd: standard deviation of your outcome
# proportion.cases: ratio between number of cases and total number of samples
# alpha: significance level
# power: power you want. 
#k: number of predictors,i.e. number of covariates plus one. Note that for multilevel factors it is Nlevels - 1,
# so for the model: outcome ~ phenotype_binary + DNA_concentration + BATCH , for 4 batches, this "k" equals to 5

install.packages("pwrss")

which_sampleSize = function(beta, Y.sd, proportion.cases, power=0.9,alpha=0.05,k=5){
  library(pwrss)
  X.sd = sqrt(proportion.cases * (1 - proportion.cases))
  pwrss.t.reg(beta1 = beta,
              k = k,
              power=power,
              alpha = alpha,
              sdx = X.sd,
              sdy = Y.sd,
              alternative = "not equal")
}



##### ALL COHORTS ################

metadata$antibiotics=as.factor(metadata$antibiotics)
summary (metadata$antibiotics)
#Proportion of cases 37/79

# For AR resistance gene load (all cohorts)
# First compute standard deviation 
sd (AR_load$total_AB_load, na.rm = T) #6847.063
which_sampleSize(4111.88,6847.063, 0.47, power=0.90,alpha=0.05,k=5)
#n = 109

# For Alpha diversity (all cohorts)
sd (my_alpha_diversity$Shannon)
which_sampleSize(0.13,0.599, 0.47, power=0.90,alpha=0.05,k=5 )
# n = 887 



