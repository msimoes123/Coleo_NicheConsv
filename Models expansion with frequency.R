### Running regression models for invasive species
## Author: M. Krull
## last update: 10/11/2022

##importing the data and loading the packages
require(brms)
require(ggplot2)
require(dplyr)


#importing the data 
ndat<-read.csv("data.csv")



## Function to calculate LOO weights
LOO.weights<-function(LOO)
{  
  min.LOO<-min(LOO) 
  diff<-LOO-min.LOO
  sum.d<-sum(exp(-diff/2))
  LOO.w<-numeric(length(LOO))
  for (i in 1:length(LOO)) LOO.w[i]<-round(exp(-diff[i]/2)/sum(exp(-diff/2)),4)
  return(cbind(LOO,LOO.w))
}



### Subsetting the dataset and log transforming skewed variables
### and creating ratios
dat<-data.frame(Family=ndat$Family,Status_Biome=ndat$Status_Biome,E_Biome=ndat$E_Biome,U_Biome=ndat$U_Biome,Feeding.habits=ndat$Feeding.habits2,Size_adult..mm.=ndat$Size_adult..mm.,size.log=ndat$size.log,nb=ndat$nb,ib=ndat$ib)
dat$nhii<-log(ndat$nhii/ndat$ihii)
dat$npop<-log(ndat$npop/ndat$ipop)
dat$nFreq<-log(ndat$nFreq/ndat$iFreq)


###########Fitting models
##one variable
### GLM assuming a beta distribution


### null model
# set prior
p <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
       set_prior("gamma(0.01, 0.01)",class="phi"))


null <- brm(E_Biome ~ 1+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(null,"models/expansion with freq/null.RDS")

### adding factors
# setting priors
p <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
       set_prior("normal(0, 2)", class = "b"),
       set_prior("gamma(0.01, 0.01)",class="phi"))

mod1 <- brm(E_Biome  ~ size.log+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod1,"models/expansion with freq/mod1.RDS")

mod2 <- brm(E_Biome  ~ Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod2,"models/expansion with freq/mod2.RDS")

mod3 <- brm(E_Biome ~ Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod3,"models/expansion with freq/mod3.RDS")

mod4 <- brm(E_Biome ~ Status_Biome+Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod4,"models/expansion with freq/mod4.RDS")

mod5 <- brm(E_Biome ~ Status_Biome+size.log+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod5,"models/expansion with freq/mod5.RDS")

mod6 <- brm(E_Biome ~ size.log+Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod6,"models/expansion with freq/mod6.RDS")

mod7 <- brm(E_Biome ~ Status_Biome+size.log+Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod7,"models/expansion with freq/mod7.RDS")

mod8 <- brm(E_Biome ~ Status_Biome*Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod8,"models/expansion with freq/mod8.RDS")

mod9 <- brm(E_Biome ~ Status_Biome*size.log+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod9,"models/expansion with freq/mod9.RDS")

mod10 <- brm(E_Biome  ~ size.log*Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod10,"models/expansion with freq/mod10.RDS")

mod11 <- brm(E_Biome ~ Status_Biome*size.log*Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod11,"models/expansion with freq/mod11.RDS")

mod12 <- brm(E_Biome ~ nhii+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod12,"models/expansion with freq/mod12.RDS")

mod13 <- brm(E_Biome ~ nb+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod13,"models/expansion with freq/mod13.RDS")

mod14 <- brm(E_Biome ~ npop+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod14,"models/expansion with freq/mod14.RDS")

mod15 <- brm(E_Biome ~ nhii+nb+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod15,"models/expansion with freq/mod15.RDS")

mod16 <- brm(E_Biome ~ nb+npop+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod16,"models/expansion with freq/mod16.RDS")

mod17 <- brm(E_Biome ~ npop+nb+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod17,"models/expansion with freq/mod17.RDS")

mod18 <- brm(E_Biome ~ npop+nb+nhii+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod18,"models/expansion with freq/mod18.RDS")

mod19 <- brm(E_Biome ~ npop+nb+nhii+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod19,"models/expansion with freq/mod19.RDS")

mod20 <- brm(E_Biome ~ npop+nb+nhii+Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod20,"models/expansion with freq/mod20.RDS")

mod21 <- brm(E_Biome ~ npop+nb+nhii+size.log+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod21,"models/expansion with freq/mod21.RDS")

mod22 <- brm(E_Biome ~ ib+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod22,"models/expansion with freq/mod22.RDS")

mod23 <- brm(E_Biome ~ nhii+ib+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod23,"models/expansion with freq/mod23.RDS")

mod24 <- brm(E_Biome ~ nb+ib+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod24,"models/expansion with freq/mod24.RDS")

mod25 <- brm(E_Biome ~ ib+npop+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod25,"models/expansion with freq/mod25.RDS")

mod26 <- brm(E_Biome ~ npop+nb+nhii+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod26,"models/expansion with freq/mod26.RDS")

mod27 <- brm(E_Biome ~ ib+nb+nhii+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod27,"models/expansion with freq/mod27.RDS")

mod28 <- brm(E_Biome ~ npop+ib+nhii+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod28,"models/expansion with freq/mod28.RDS")

mod29 <- brm(E_Biome ~ npop+nb+ib+Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod29,"models/expansion with freq/mod29.RDS")

mod30 <- brm(E_Biome ~ npop+ib+nhii+size.log+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod30,"models/expansion with freq/mod30.RDS")

mod31 <- brm(E_Biome ~ npop+ib+nhii+Status_Biome+Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod31,"models/expansion with freq/mod31.RDS")

mod32 <- brm(E_Biome ~ npop+nb+ib+Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod32,"models/expansion with freq/mod32.RDS")

mod33 <- brm(E_Biome ~ npop+ib+nhii+size.log*Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod33,"models/expansion with freq/mod33.RDS")

mod34 <- brm(E_Biome ~ nb+ib+size.log*Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod34,"models/expansion with freq/mod34.RDS")

mod35 <- brm(E_Biome ~ ib+size.log*Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod35,"models/expansion with freq/mod35.RDS")

mod36 <- brm(E_Biome ~ nb+size.log*Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod36,"models/expansion with freq/mod36.RDS")

mod37 <- brm(E_Biome ~ nb+ib+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod37,"models/expansion with freq/mod37.RDS")

mod38 <- brm(E_Biome ~ nb+ib+size.log*Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod38,"models/expansion with freq/mod38.RDS")

mod39 <- brm(E_Biome ~ ib+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod39,"models/expansion with freq/mod39.RDS")

mod40 <- brm(E_Biome ~ nb+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod40,"models/expansion with freq/mod40.RDS")

mod41 <- brm(E_Biome ~ nb+Status_Biome+nhii+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod41,"models/expansion with freq/mod41.RDS")

mod42 <- brm(E_Biome ~ ib+Status_Biome+size.log+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod42,"models/expansion with freq/mod42.RDS")


###############################################################################################
################################## Random intercept only models ###############################
###############################################################################################

### null model
# set prior
p <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
       set_prior("gamma(0.01, 0.01)",class="phi"))


null_random <- brm(E_Biome ~ 1+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(null_random,"models/expansion with freq/null_random.RDS")

### adding factors
# set prior
p <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
       set_prior("normal(0, 2)", class = "b"),
       set_prior("gamma(0.01, 0.01)",class="phi"))



mod1_random <- brm(E_Biome  ~ size.log+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod1_random,"models/expansion with freq/mod1_random.RDS")

mod2_random <- brm(E_Biome  ~ Feeding.habits+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod2_random,"models/expansion with freq/mod2_random.RDS")

mod3_random <- brm(E_Biome ~ Status_Biome+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod3_random,"models/expansion with freq/mod3_random.RDS")

mod4_random <- brm(E_Biome ~ Status_Biome+Feeding.habits+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod4_random,"models/expansion with freq/mod4_random.RDS")

mod5_random <- brm(E_Biome ~ Status_Biome+size.log+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod5_random,"models/expansion with freq/mod5_random.RDS")

mod6_random <- brm(E_Biome ~ size.log+Feeding.habits+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod6_random,"models/expansion with freq/mod6_random.RDS")

mod7_random <- brm(E_Biome ~ Status_Biome+size.log+Feeding.habits+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod7_random,"models/expansion with freq/mod7_random.RDS")

mod8_random <- brm(E_Biome ~ Status_Biome*Feeding.habits+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod8_random,"models/expansion with freq/mod8_random.RDS")

mod9_random <- brm(E_Biome ~ Status_Biome*size.log+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod9_random,"models/expansion with freq/mod9_random.RDS")

mod10_random <- brm(E_Biome  ~ size.log*Feeding.habits+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod10_random,"models/expansion with freq/mod10_random.RDS")

mod11_random <- brm(E_Biome ~ Status_Biome*size.log*Feeding.habits+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod11_random,"models/expansion with freq/mod11_random.RDS")

mod12_random <- brm(E_Biome ~ nhii+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod12_random,"models/expansion with freq/mod12_random.RDS")

mod13_random <- brm(E_Biome ~ nb+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod13_random,"models/expansion with freq/mod13_random.RDS")

mod14_random <- brm(E_Biome ~ npop+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod14_random,"models/expansion with freq/mod14_random.RDS")

mod15_random <- brm(E_Biome ~ nhii+nb+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod15_random,"models/expansion with freq/mod15_random.RDS")

mod16_random <- brm(E_Biome ~ nb+npop+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod16_random,"models/expansion with freq/mod16_random.RDS")

mod17_random <- brm(E_Biome ~ npop+nb+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod17_random,"models/expansion with freq/mod17_random.RDS")

mod18_random <- brm(E_Biome ~ npop+nb+nhii+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod18_random,"models/expansion with freq/mod18_random.RDS")

mod19_random <- brm(E_Biome ~ npop+nb+nhii+Status_Biome+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod19_random,"models/expansion with freq/mod19_random.RDS")

mod20_random <- brm(E_Biome ~ npop+nb+nhii+Feeding.habits+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod20_random,"models/expansion with freq/mod20_random.RDS")

mod21_random <- brm(E_Biome ~ npop+nb+nhii+size.log+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod21_random,"models/expansion with freq/mod21_random.RDS")

mod22_random <- brm(E_Biome ~ ib+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod22_random,"models/expansion with freq/mod22_random.RDS")

mod23_random <- brm(E_Biome ~ nhii+ib+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod23_random,"models/expansion with freq/mod23_random.RDS")

mod24_random <- brm(E_Biome ~ nb+ib+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod24_random,"models/expansion with freq/mod24_random.RDS")

mod25_random <- brm(E_Biome ~ ib+npop+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod25_random,"models/expansion with freq/mod25_random.RDS")

mod26_random <- brm(E_Biome ~ npop+nb+nhii+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod26_random,"models/expansion with freq/mod26_random.RDS")

mod27_random <- brm(E_Biome ~ ib+nb+nhii+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod27_random,"models/expansion with freq/mod27_random.RDS")

mod28_random <- brm(E_Biome ~ npop+ib+nhii+Status_Biome+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod28_random,"models/expansion with freq/mod28_random.RDS")

mod29_random <- brm(E_Biome ~ npop+nb+ib+Feeding.habits+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod29_random,"models/expansion with freq/mod29_random.RDS")

mod30_random <- brm(E_Biome ~ npop+ib+nhii+size.log+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod30_random,"models/expansion with freq/mod30_random.RDS")

mod31_random <- brm(E_Biome ~ npop+ib+nhii+Status_Biome+Feeding.habits+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod31_random,"models/expansion with freq/mod31_random.RDS")

mod32_random <- brm(E_Biome ~ npop+nb+ib+Feeding.habits+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod32_random,"models/expansion with freq/mod32_random.RDS")

mod33_random <- brm(E_Biome ~ npop+ib+nhii+size.log*Status_Biome+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod33_random,"models/expansion with freq/mod33_random.RDS")

mod34_random <- brm(E_Biome ~ nb+ib+size.log*Status_Biome+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod34_random,"models/expansion with freq/mod34_random.RDS")

mod35_random <- brm(E_Biome ~ ib+size.log*Status_Biome+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod35_random,"models/expansion with freq/mod35_random.RDS")

mod36_random <- brm(E_Biome ~ nb+size.log*Status_Biome+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod36_random,"models/expansion with freq/mod36_random.RDS")

mod37_random <- brm(E_Biome ~ nb+ib+Status_Biome+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod37_random,"models/expansion with freq/mod37_random.RDS")

mod38_random <- brm(E_Biome ~ nb+ib+size.log*Status_Biome+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod38_random,"models/expansion with freq/mod38_random.RDS")

mod39_random <- brm(E_Biome ~ ib+Status_Biome+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod39_random,"models/expansion with freq/mod39_random.RDS")

mod40_random <- brm(E_Biome ~ nb+Status_Biome+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod40_random,"models/expansion with freq/mod40_random.RDS")

mod41_random <- brm(E_Biome ~ nb+Status_Biome+nhii+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod41_random,"models/expansion with freq/mod41_random.RDS")

mod42_random <- brm(E_Biome ~ ib+Status_Biome+size.log+nFreq+(1|Family), data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod42_random,"models/expansion with freq/mod42_random.RDS")


###############################################################################################
############################ Random slopes and intercept models ###############################
###############################################################################################

p <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
       set_prior("normal(0, 2)", class = "b"))

mod1_random_slope <- brm(E_Biome ~ size.log+ (1+size.log|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod1_random_slope,"models/expansion with freq/mod1_random_slope.RDS")

mod2_random_slope <- brm(E_Biome ~ Feeding.habits+ (1+Feeding.habits|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod2_random_slope,"models/expansion with freq/mod2_random_slope.RDS")

mod3_random_slope <- brm(E_Biome ~ Status_Biome+ (1+Status_Biome|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod3_random_slope,"models/expansion with freq/mod3_random_slope.RDS")

mod4_random_slope <- brm(E_Biome ~ Status_Biome+Feeding.habits+ (1+Status_Biome+Feeding.habits|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod4_random_slope,"models/expansion with freq/mod4_random_slope.RDS")

mod5_random_slope <- brm(E_Biome ~ Status_Biome+size.log+ (1+size.log|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod5_random_slope,"models/expansion with freq/mod5_random_slope.RDS")

mod6_random_slope <- brm(E_Biome ~ size.log+Feeding.habits+ (1+size.log|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod6_random_slope,"models/expansion with freq/mod6_random_slope.RDS")

mod7_random_slope <- brm(E_Biome ~ Status_Biome+size.log+Feeding.habits+ (1+size.log|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod7_random_slope,"models/expansion with freq/mod7_random_slope.RDS")

mod8_random_slope <- brm(E_Biome ~ Status_Biome*Feeding.habits+ (1+Feeding.habits|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod8_random_slope,"models/expansion with freq/mod8_random_slope.RDS")

mod9_random_slope <- brm(E_Biome ~ Status_Biome*size.log+ (1+size.log|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod9_random_slope,"models/expansion with freq/mod9_random_slope.RDS")

mod10_random_slope <- brm(E_Biome ~ size.log*Feeding.habits+ (1+size.log|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod10_random_slope,"models/expansion with freq/mod10_random_slope.RDS")

mod11_random_slope <- brm(E_Biome ~ Status_Biome*size.log*Feeding.habits+ (1+size.log|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod11_random_slope,"models/expansion with freq/mod11_random_slope.RDS")

mod12_random_slope <- brm(E_Biome ~ nhii+ (1+nhii|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod12_random_slope,"models/expansion with freq/mod12_random_slope.RDS")

mod13_random_slope <- brm(E_Biome ~ nb+ (1+nb|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod13_random_slope,"models/expansion with freq/mod13_random_slope.RDS")

mod14_random_slope <- brm(E_Biome ~ npop+ (1+npop|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod14_random_slope,"models/expansion with freq/mod14_random_slope.RDS")

mod15_random_slope <- brm(E_Biome ~ nhii+nb+ (1+nhii+nb|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod15_random_slope,"models/expansion with freq/mod15_random_slope.RDS")

mod16_random_slope <- brm(E_Biome ~ nb+npop+ (1+nb+npop|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod16_random_slope,"models/expansion with freq/mod16_random_slope.RDS")

mod17_random_slope <- brm(E_Biome ~ npop+nb+ (1+npop+nb|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod17_random_slope,"models/expansion with freq/mod17_random_slope.RDS")

mod18_random_slope <- brm(E_Biome ~ npop+nb+nhii+ (1+npop+nb+nhii|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod18_random_slope,"models/expansion with freq/mod18_random_slope.RDS")

mod19_random_slope <- brm(E_Biome ~ npop+nb+nhii+ (1+npop+nb+nhii|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod19_random_slope,"models/expansion with freq/mod19_random_slope.RDS")

mod20_random_slope <- brm(E_Biome ~ npop+nb+nhii+Status_Biome+ (1+npop+nb+nhii|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod20_random_slope,"models/expansion with freq/mod20_random_slope.RDS")

mod21_random_slope <- brm(E_Biome ~ npop+nb+nhii+Feeding.habits+ (1+npop+nb+nhii|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod21_random_slope,"models/expansion with freq/mod21_random_slope.RDS")

mod22_random_slope <- brm(E_Biome ~ npop+nb+nhii+size.log+ (1+npop+nb+nhii+size.log|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod22_random_slope,"models/expansion with freq/mod22_random_slope.RDS")

mod23_random_slope <- brm(E_Biome ~ npop+nb+nhii+Status_Biome+ (1+npop+nb+nhii|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod23_random_slope,"models/expansion with freq/mod23_random_slope.RDS")

mod24_random_slope <- brm(E_Biome ~ npop+nb+nhii+Feeding.habits+ (1+npop+nb+nhii|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod24_random_slope,"models/expansion with freq/mod24_random_slope.RDS")

mod25_random_slope <- brm(E_Biome ~ npop+nb+nhii+size.log+ (1+npop+nb+nhii+size.log|Family)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod25_random_slope,"models/expansion with freq/mod25_random_slope.RDS")


###############################################################################################
######################################## GAM models ###########################################
###############################################################################################

mod1_gam <- brm(E_Biome  ~ s(size.log,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod1_gam,"models/expansion with freq/mod1_gam.RDS")

mod2_gam <- brm(E_Biome ~ Status_Biome+s(size.log,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod2_gam,"models/expansion with freq/mod2_gam.RDS")

mod3_gam <- brm(E_Biome ~ s(size.log,k=3)+Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod3_gam,"models/expansion with freq/mod3_gam.RDS")

mod4_gam <- brm(E_Biome ~ Status_Biome+s(size.log,k=3)+Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod4_gam,"models/expansion with freq/mod4_gam.RDS")

mod5_gam <- brm(E_Biome ~ s(nhii,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod5_gam,"models/expansion with freq/mod5_gam.RDS")

mod6_gam <- brm(E_Biome ~ s(nb,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod6_gam,"models/expansion with freq/mod6_gam.RDS")

mod7_gam <- brm(E_Biome ~ s(npop,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod7_gam,"models/expansion with freq/mod7_gam.RDS")

mod8_gam <- brm(E_Biome ~ s(nhii,k=3)+s(nb,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod8_gam,"models/expansion with freq/mod8_gam.RDS")

mod9_gam <- brm(E_Biome ~ s(nb,k=3)+s(npop,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod9_gam,"models/expansion with freq/mod9_gam.RDS")

mod10_gam <- brm(E_Biome ~ s(npop,k=3)+s(ib,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod10_gam,"models/expansion with freq/mod10_gam.RDS")

mod11_gam <- brm(E_Biome ~ s(npop,k=3)+s(nb,k=3)+s(nhii,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod11_gam,"models/expansion with freq/mod11_gam.RDS")

mod12_gam <- brm(E_Biome ~ s(npop,k=3)+s(nb,k=3)+s(nhii,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod12_gam,"models/expansion with freq/mod12_gam.RDS")

mod13_gam <- brm(E_Biome ~ s(npop,k=3)+s(nb,k=3)+s(nhii,k=3)+Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod13_gam,"models/expansion with freq/mod13_gam.RDS")

mod14_gam <- brm(E_Biome ~ s(npop,k=3)+s(nb,k=3)+s(nhii,k=3)+s(size.log,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod14_gam,"models/expansion with freq/mod14_gam.RDS")

mod15_gam <- brm(E_Biome ~ s(npop,k=3)+s(nb,k=3)+s(nhii,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod15_gam,"models/expansion with freq/mod15_gam.RDS")

mod16_gam <- brm(E_Biome ~ s(npop,k=3)+s(nb,k=3)+s(nhii,k=3)+Feeding.habits+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod16_gam,"models/expansion with freq/mod16_gam.RDS")

mod17_gam <- brm(E_Biome ~ s(npop,k=3)+s(nb,k=3)+s(nhii,k=3)+s(size.log,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod17_gam,"models/expansion with freq/mod17_gam.RDS")

mod18_gam <- brm(E_Biome ~ s(ib,k=3)+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod18_gam,"models/expansion with freq/mod18_gam.RDS")

mod19_gam <- brm(E_Biome ~ s(ib,k=3)+s(nb,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod19_gam,"models/expansion with freq/mod19_gam.RDS")

mod20_gam <- brm(E_Biome ~ s(ib,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod20_gam,"models/expansion with freq/mod20_gam.RDS")

mod21_gam <- brm(E_Biome ~ s(nb,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod21_gam,"models/expansion with freq/mod21_gam.RDS")

mod22_gam <- brm(E_Biome ~ s(size.log,k=3)+s(ib,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod22_gam,"models/expansion with freq/mod22_gam.RDS")

mod23_gam <- brm(E_Biome ~ s(size.log,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod23_gam,"models/expansion with freq/mod23_gam.RDS")

mod24_gam <- brm(E_Biome ~ s(size.log,by=Status_Biome,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod24_gam,"models/expansion with freq/mod24_gam.RDS")

mod25_gam <- brm(E_Biome ~ s(size.log,k=3)+s(ib,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod25_gam,"models/expansion with freq/mod25_gam.RDS")

mod26_gam <- brm(E_Biome ~ s(nb,k=3)+s(nhii,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod26_gam,"models/expansion with freq/mod26_gam.RDS")

mod27_gam <- brm(E_Biome ~ s(size.log,k=3)+s(ib,k=3)+s(nb,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod27_gam,"models/expansion with freq/mod27_gam.RDS")

mod28_gam <- brm(E_Biome ~ s(size.log,k=3,by=Status_Biome)+s(ib,k=3)+Status_Biome+nFreq, data=dat, family="beta",prior = p,iter = 5000, control=list(adapt_delta=0.99))
saveRDS(mod28_gam,"models/expansion with freq/mod28_gam.RDS")

##model selection
## calculating the LOO information criteria
null<-add_criterion(null, "loo",moment_match = F)
mod1<-add_criterion(mod1, "loo",moment_match = F)
mod2<-add_criterion(mod2, "loo",moment_match = F)
mod3<-add_criterion(mod3, "loo",moment_match = F)
mod4<-add_criterion(mod4, "loo",moment_match = F)
mod5<-add_criterion(mod5, "loo",moment_match = F)
mod6<-add_criterion(mod6, "loo",moment_match = F)
mod7<-add_criterion(mod7, "loo",moment_match = F)
mod8<-add_criterion(mod8, "loo",moment_match = F)
mod9<-add_criterion(mod9, "loo",moment_match = F)
mod10<-add_criterion(mod10, "loo",moment_match = F)
mod11<-add_criterion(mod11, "loo",moment_match = F)
mod12<-add_criterion(mod12, "loo",moment_match = F)
mod13<-add_criterion(mod13, "loo",moment_match = F)
mod14<-add_criterion(mod14, "loo",moment_match = F)
mod15<-add_criterion(mod15, "loo",moment_match = F)
mod16<-add_criterion(mod16, "loo",moment_match = F)
mod17<-add_criterion(mod17, "loo",moment_match = F)
mod18<-add_criterion(mod18, "loo",moment_match = F)
mod19<-add_criterion(mod19, "loo",moment_match = F)
mod20<-add_criterion(mod20, "loo",moment_match = F)
mod21<-add_criterion(mod21, "loo",moment_match = F)
mod22<-add_criterion(mod22, "loo",moment_match = F)
mod23<-add_criterion(mod23, "loo",moment_match = F)
mod24<-add_criterion(mod24, "loo",moment_match = F)
mod25<-add_criterion(mod25, "loo",moment_match = F)
mod26<-add_criterion(mod26, "loo",moment_match = F)
mod27<-add_criterion(mod27, "loo",moment_match = F)
mod28<-add_criterion(mod28, "loo",moment_match = F)
mod29<-add_criterion(mod29, "loo",moment_match = F)
mod30<-add_criterion(mod30, "loo",moment_match = F)
mod31<-add_criterion(mod31, "loo",moment_match = F)
mod32<-add_criterion(mod32, "loo",moment_match = F)
mod33<-add_criterion(mod33, "loo",moment_match = F)
mod34<-add_criterion(mod34, "loo",moment_match = F)
mod35<-add_criterion(mod35, "loo",moment_match = F)
mod36<-add_criterion(mod36, "loo",moment_match = F)
mod37<-add_criterion(mod37, "loo",moment_match = F)
mod38<-add_criterion(mod38, "loo",moment_match = F)
mod39<-add_criterion(mod39, "loo",moment_match = F)
mod40<-add_criterion(mod40, "loo",moment_match = F)
mod41<-add_criterion(mod41, "loo",moment_match = F)
mod42<-add_criterion(mod42, "loo",moment_match = F)

null_random<-add_criterion(null_random, "loo",moment_match = F)
mod1_random<-add_criterion(mod1_random, "loo",moment_match = F)
mod2_random<-add_criterion(mod2_random, "loo",moment_match = F)
mod3_random<-add_criterion(mod3_random, "loo",moment_match = F)
mod4_random<-add_criterion(mod4_random, "loo",moment_match = F)
mod5_random<-add_criterion(mod5_random, "loo",moment_match = F)
mod6_random<-add_criterion(mod6_random, "loo",moment_match = F)
mod7_random<-add_criterion(mod7_random, "loo",moment_match = F)
mod8_random<-add_criterion(mod8_random, "loo",moment_match = F)
mod9_random<-add_criterion(mod9_random, "loo",moment_match = F)
mod10_random<-add_criterion(mod10_random, "loo",moment_match = F)
mod11_random<-add_criterion(mod11_random, "loo",moment_match = F)
mod12_random<-add_criterion(mod12_random, "loo",moment_match = F)
mod13_random<-add_criterion(mod13_random, "loo",moment_match = F)
mod14_random<-add_criterion(mod14_random, "loo",moment_match = F)
mod15_random<-add_criterion(mod15_random, "loo",moment_match = F)
mod16_random<-add_criterion(mod16_random, "loo",moment_match = F)
mod17_random<-add_criterion(mod17_random, "loo",moment_match = F)
mod18_random<-add_criterion(mod18_random, "loo",moment_match = F)
mod19_random<-add_criterion(mod19_random, "loo",moment_match = F)
mod20_random<-add_criterion(mod20_random, "loo",moment_match = F)
mod21_random<-add_criterion(mod21_random, "loo",moment_match = F)
mod22_random<-add_criterion(mod22_random, "loo",moment_match = F)
mod23_random<-add_criterion(mod23_random, "loo",moment_match = F)
mod24_random<-add_criterion(mod24_random, "loo",moment_match = F)
mod25_random<-add_criterion(mod25_random, "loo",moment_match = F)
mod26_random<-add_criterion(mod26_random, "loo",moment_match = F)
mod27_random<-add_criterion(mod27_random, "loo",moment_match = F)
mod28_random<-add_criterion(mod28_random, "loo",moment_match = F)
mod29_random<-add_criterion(mod29_random, "loo",moment_match = F)
Mod30_random <-add_criterion(Mod30_random, "loo",moment_match = F)
mod31_random <-add_criterion(mod31_random, "loo",moment_match = F)
mod32_random <-add_criterion(mod32_random, "loo",moment_match = F)
mod33_random <-add_criterion(mod33_random, "loo",moment_match = F)
mod34_random <-add_criterion(mod34_random, "loo",moment_match = F)
mod35_random <-add_criterion(mod35_random, "loo",moment_match = F)
mod36_random <-add_criterion(mod36_random, "loo",moment_match = F)
mod37_random <-add_criterion(mod37_random, "loo",moment_match = F)
mod38_random <-add_criterion(mod38_random, "loo",moment_match = F)
mod39_random <-add_criterion(mod39_random, "loo",moment_match = F)
mod40_random <-add_criterion(mod40_random, "loo",moment_match = F)
mod41_random <-add_criterion(mod41_random, "loo",moment_match = F)
mod42_random <-add_criterion(mod42_random, "loo",moment_match = F)
mod43_random <-add_criterion(mod43_random, "loo",moment_match = F)
mod44_random <-add_criterion(mod44_random, "loo",moment_match = F)
mod45_random <-add_criterion(mod45_random, "loo",moment_match = F)
mod46_random <-add_criterion(mod46_random, "loo",moment_match = F)
mod47_random <-add_criterion(mod47_random, "loo",moment_match = F)
mod48_random <-add_criterion(mod48_random, "loo",moment_match = F)
mod49_random <-add_criterion(mod49_random, "loo",moment_match = F)


null_random_slope<-add_criterion(null_random_slope, "loo",moment_match = F)
mod1_random_slope<-add_criterion(mod1_random_slope, "loo",moment_match = F)
mod2_random_slope<-add_criterion(mod2_random_slope, "loo",moment_match = F)
mod3_random_slope<-add_criterion(mod3_random_slope, "loo",moment_match = F)
mod4_random_slope<-add_criterion(mod4_random_slope, "loo",moment_match = F)
mod5_random_slope<-add_criterion(mod5_random_slope, "loo",moment_match = F)
mod6_random_slope<-add_criterion(mod6_random_slope, "loo",moment_match = F)
mod7_random_slope<-add_criterion(mod7_random_slope, "loo",moment_match = F)
mod8_random_slope<-add_criterion(mod8_random_slope, "loo",moment_match = F)
mod9_random_slope<-add_criterion(mod9_random_slope, "loo",moment_match = F)
mod10_random_slope<-add_criterion(mod10_random_slope, "loo",moment_match = F)
mod11_random_slope<-add_criterion(mod11_random_slope, "loo",moment_match = F)
mod12_random_slope<-add_criterion(mod12_random_slope, "loo",moment_match = F)
mod13_random_slope<-add_criterion(mod13_random_slope, "loo",moment_match = F)
mod14_random_slope<-add_criterion(mod14_random_slope, "loo",moment_match = F)
mod15_random_slope<-add_criterion(mod15_random_slope, "loo",moment_match = F)
mod16_random_slope<-add_criterion(mod16_random_slope, "loo",moment_match = F)
mod17_random_slope<-add_criterion(mod17_random_slope, "loo",moment_match = F)
mod18_random_slope<-add_criterion(mod18_random_slope, "loo",moment_match = F)
mod19_random_slope<-add_criterion(mod19_random_slope, "loo",moment_match = F)
mod20_random_slope<-add_criterion(mod20_random_slope, "loo",moment_match = F)
mod21_random_slope<-add_criterion(mod21_random_slope, "loo",moment_match = F)
mod22_random_slope<-add_criterion(mod22_random_slope, "loo",moment_match = F)
mod23_random_slope<-add_criterion(mod23_random_slope, "loo",moment_match = F)
mod24_random_slope<-add_criterion(mod24_random_slope, "loo",moment_match = F)
mod25_random_slope<-add_criterion(mod25_random_slope, "loo",moment_match = F)
mod26_random_slope<-add_criterion(mod26_random_slope, "loo",moment_match = F)
mod27_random_slope<-add_criterion(mod27_random_slope, "loo",moment_match = F)
mod28_random_slope<-add_criterion(mod28_random_slope, "loo",moment_match = F)
mod29_random_slope<-add_criterion(mod29_random_slope, "loo",moment_match = F)


mod1_gam<-add_criterion(mod1_gam, "loo",moment_match = F)
mod2_gam<-add_criterion(mod2_gam, "loo",moment_match = F)
mod3_gam<-add_criterion(mod3_gam, "loo",moment_match = F)
mod4_gam<-add_criterion(mod4_gam, "loo",moment_match = F)
mod5_gam<-add_criterion(mod5_gam, "loo",moment_match = F)
mod6_gam<-add_criterion(mod6_gam, "loo",moment_match = F)
mod7_gam<-add_criterion(mod7_gam, "loo",moment_match = F)
mod8_gam<-add_criterion(mod8_gam, "loo",moment_match = F)
mod9_gam<-add_criterion(mod9_gam, "loo",moment_match = F)
mod10_gam<-add_criterion(mod10_gam, "loo",moment_match = F)
mod11_gam<-add_criterion(mod11_gam, "loo",moment_match = F)
mod12_gam<-add_criterion(mod12_gam, "loo",moment_match = F)
mod13_gam<-add_criterion(mod13_gam, "loo",moment_match = F)
mod14_gam<-add_criterion(mod14_gam, "loo",moment_match = F)
mod15_gam<-add_criterion(mod15_gam, "loo",moment_match = F)
mod16_gam<-add_criterion(mod16_gam, "loo",moment_match = F)
mod17_gam<-add_criterion(mod17_gam, "loo",moment_match = F)
mod18_gam<-add_criterion(mod18_gam, "loo",moment_match = F)
mod19_gam<-add_criterion(mod19_gam, "loo",moment_match = F)
mod20_gam<-add_criterion(mod20_gam, "loo",moment_match = F)
mod21_gam<-add_criterion(mod21_gam, "loo",moment_match = F)
mod22_gam<-add_criterion(mod22_gam, "loo",moment_match = F)
mod23_gam<-add_criterion(mod23_gam, "loo",moment_match = F)
mod24_gam<-add_criterion(mod24_gam, "loo",moment_match = F)
mod25_gam<-add_criterion(mod25_gam, "loo",moment_match = F)
mod26_gam<-add_criterion(mod26_gam, "loo",moment_match = F)
mod27_gam<-add_criterion(mod27_gam, "loo",moment_match = F)
mod28_gam<-add_criterion(mod28_gam, "loo",moment_match = F)

### comparing all models
comp <- loo_compare(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,
                    mod12, mod13, mod14, mod15, mod16,mod17, mod18, mod19, mod20,mod21,   
                    mod22, mod23,mod24, mod25, mod26,mod27,mod28,  mod29, mod30, mod31,
                    mod32, mod33, mod34, mod35, mod36, mod37, mod38, mod39, mod40, 
                    mod41,mod42,mod43, mod44, mod45, mod46, mod47, mod48, mod49,mod50,mod51,
                    
                    ## random intercept
                    null_random,mod1_random, mod2_random, mod3_random, mod4_random,mod5_random,
                    mod6_random, mod7_random, mod8_random, mod9_random, mod10_random,
                    mod11_random, mod12_random, mod13_random,mod14_random,mod15_random,mod16_random,
                    mod17_random, mod18_random,  mod19_random, mod20_random, mod21_random,   
                    mod22_random, mod23_random, mod24_random, mod25_random, mod26_random, 
                    mod27_random, mod28_random, mod29_random, Mod30_random, mod31_random,
                    mod32_random, mod33_random, mod34_random, mod35_random, mod36_random, 
                    mod37_random, mod38_random, mod39_random, mod40_random, mod41_random,
                    mod42_random, mod43_random, mod44_random, mod45_random, mod46_random, 
                    mod47_random, mod48_random, mod49_random,
                    
                    #random intercept and slope
                    null_random_slope,mod1_random_slope, mod2_random_slope, mod3_random_slope, mod4_random_slope,
                    mod5_random_slope, mod6_random_slope, mod7_random_slope, mod8_random_slope,
                    mod9_random_slope, mod10_random_slope, mod11_random_slope, mod12_random_slope,
                    mod13_random_slope,mod14_random_slope,mod15_random_slope,mod16_random_slope,
                    mod17_random_slope,mod18_random_slope, mod19_random_slope,mod20_random_slope,
                    mod21_random_slope,mod22_random_slope,mod23_random_slope,mod24_random_slope,
                    mod25_random_slope,mod26_random_slope,mod27_random_slope,mod28_random_slope,mod28_random_slope,mod29_random_slope,
                    
                    # gamm
                    mod1_gam, mod2_gam, mod3_gam, mod4_gam, mod5_gam, mod6_gam, mod7_gam,
                    mod8_gam, mod9_gam, mod10_gam, mod11_gam, mod12_gam, mod13_gam,
                    mod14_gam, mod15_gam, mod16_gam, mod17_gam, mod18_gam, mod19_gam,
                    mod20_gam, mod21_gam, mod22_gam, mod23_gam, mod24_gam, mod25_gam,
                    mod26_gam, mod27_gam, mod28_gam, mod29_gam, mod30_gam)

### comparing the models                     
result<-data.frame(print(comp,simplify = FALSE, digits = 3))
### calculating LOOic weights
result<-data.frame(result,weights=LOO.weights(result$looic)[,2])
head(results)