#=========================== LDA1 - MALARIA MODELLING ============
#====== CAMPBELL MCDULING
#====== 16 MAY 2023
rm(list = ls())
setwd("~/OneDrive - University of Cape Town/2023/MSc Biostatistics/Coursework Semester 1/LDA/A1")
library(lme4)
require(tidyverse)
library(dplyr)
library(psych)
library(ggplot2)
library(nlme)
library(lattice)
library(mice)
library(ggeffects)

#==================== DATA MANAGEMENT =============
dat0 = read.csv("malariadata.csv")
head(dat0)
summary(dat0)

describe(dat0)
str(dat0)

#define variable types accordingly
dat0$site = as.factor(dat0$site)
dat0$arm = as.factor(dat0$arm)
dat0$pid = as.factor(dat0$pid)
dat0$gender = as.factor(dat0$gender)
dat0$country = as.factor(dat0$country)
dat0$PIoutcome = as.factor(dat0$PIoutcome)
summary(dat0)

#create weight and age quartiles
quart <- quantile(dat0$weight, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = T)
dat0$weight_quart <- cut(dat0$weight, breaks = quart, labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = T)
quart <- quantile(dat0$age, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = T)
dat0$age_quart <- cut(dat0$age, breaks = quart, labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = T)



#--------- sample 75% of patients
n_distinct(dat0$pid)

unique_participants <- unique(dat0$pid)

set.seed(19980723)
sampled_participants <- sample(unique_participants, size=0.75*length(unique_participants), replace = FALSE)

# Filter the original dataframe to retain observations for the sampled participants
dat <- dat0 %>%
  filter(pid %in% sampled_participants)

# Check - should have 306 participants now
n_distinct(dat$pid)


#--------- create wide df
dat.wide = reshape(dat, timevar = "pday",  idvar="pid",direction="wide", v.names = c("Hb"), 
                          drop=c("pardens", "gamedens", "Pyrconcentration", "Sulfconcentration", "Studyyear", "PIoutcome", "PIoutcomeday"))

#==================== DATA VALIDATION =============
#--------- MISSING VALUES
colSums(is.na(dat)) #check missing patterns across variables -> many missing outcomes, especially Hb
na.row = rowSums(is.na(dat)) #idenfity observations with many missings
na.row.many = which(na.row>4)
View(dat[na.row.many,]) # 52 observations with no outcomes recorded
n_distinct(dat[na.row.many,]$pid) #of these, 45 unique pids have no outcomes recorded at at least 1 time point

md.pattern(dat, plot = F)
md.pairs(dat)

#missing Hb values
sum(is.na(dat$Hb)) #1100 missing Hb observations
n_distinct(dat[is.na(dat$Hb),]$pid) #over 301 subjects (only 5 have no NA Hb vals, probably because they didnt have the missing visit (1 and 2) recorded)

#missing Hb values per time point - there are no observations at day1, 2, and 21
sum(is.na(dat.wide$Hb.0))
sum(is.na(dat.wide$Hb.1))
sum(is.na(dat.wide$Hb.2))
sum(is.na(dat.wide$Hb.3))
sum(is.na(dat.wide$Hb.7))
sum(is.na(dat.wide$Hb.14))
sum(is.na(dat.wide$Hb.21))
sum(is.na(dat.wide$Hb.28))
sum(is.na(dat.wide$Hb.42))




### REMOVE all observations at days 1, 2, and 21
dat <- filter(dat, !pday %in% c(1, 2, 21))
n_distinct(dat$pid) #still 306 pids
sum(is.na(dat$Hb)) #now only 295 missing Hb vals

#create wide df again
#--------- create wide df
dat.wide = reshape(dat, timevar = "pday",  idvar="pid",direction="wide", v.names = c("Hb"), 
                   drop=c("pardens", "gamedens", "Pyrconcentration", "Sulfconcentration", "Studyyear", "PIoutcome", "PIoutcomeday"))

#missing Hb per time point - which participants have many missings
sum(rowSums(is.na(dat.wide[,10:15]))==5) #16 are missing 5/6 Hb 
index.onlybase = which(rowSums(is.na(dat.wide[,10:15]))==5)
dat.wide[index.onlybase,] #these only have baseline
sum(rowSums(is.na(dat.wide[,10:15]))==4) #26 are missing 4/6 Hb measurements
which(rowSums(is.na(dat.wide[,10:15]))==4)

# remove all 16 participants who only have baseline Hb 
pid.onlybase = dat.wide[index.onlybase,]$pid
index.onlybase.long = which(dat$pid %in% pid.onlybase)
View(dat[index.onlybase.long,]) #confirmed these only have baseline Hb
dat = dat[-index.onlybase.long,]
n_distinct(dat$pid) #sample size reduced to 290

### REMOVE THOSE OBS MISSING Hb from other timepoints -- just for rest of exploratory analysis
na.row.Hb = which(is.na(dat$Hb)==T)
n_distinct(dat$pid) #290 pids
colSums(is.na(dat))
dat.fullHB = dat[-na.row.Hb,] 
colSums(is.na(dat.fullHB))
n_distinct(dat.fullHB$pid) #290 pids still

#--------- create wide df from df with no missing Hb
dat.wide.fullHB = reshape(dat.fullHB, timevar = "pday",  idvar="pid",direction="wide", v.names = c("Hb"), 
                          drop=c("pardens", "gamedens", "Pyrconcentration", "Sulfconcentration", "Studyyear", "PIoutcome", "PIoutcomeday"))





#====================  DATA EXPLORATION =============
#------- DESCRIPTIVE STATS (gender, age, weight, site, arm, study year)

describe(dat)

#baseline statistics by study arm
dat.baseline = dat %>% filter(pday==0)

#continuous variables
library(purrr)
dat.baseline[,c(2, 7, 11, 12)] %>% split(.$arm) %>% map(summary)

#discrete variables
dat.baseline[,c(1,2, 10)] %>% split(.$arm) %>% map(summary)

#descriptive plots
ggplot(dat, aes(x=age, col=arm)) + geom_density()
ggplot(dat, aes(x=weight, col=arm)) + geom_density()
ggplot(dat, aes(x=gender)) + geom_bar()


#------- OUTCOMES
### UNIVARIATE
ggplot(dat, aes(x=Hb)) + geom_density() #this is highly skewed, must be an outlying small value
ggplot(dat, aes(y=Hb, x=pid)) + geom_point() #there is one zero value which is outlying

dat = scrub(dat, where = "Hb", isvalue=0, newvalue = NA) #replace outlying value with NA
dat.fullHB = scrub(dat.fullHB, where = "Hb", isvalue=0, newvalue = NA) #replace outlying value with NA
na.row.Hb = which(is.na(dat.fullHB$Hb))
dat.fullHB = dat.fullHB[-na.row.Hb,] #remove this new NA from the complete data


ggplot(dat, aes(x=Hb)) + geom_density() #this looks more realistic
ggplot(dat.baseline, aes(x=Hb)) + geom_density() #this looks more realistic
summary(dat.baseline$Hb)


### BIVARIATE
dat.baseline[,c(2,7)] %>% split(.$arm) %>% map(summary)
ggplot(dat.baseline, aes(x=Hb, col=arm)) + geom_density()
ggplot(dat.baseline, aes(x=Hb, col=gender)) + geom_density()
ggplot(dat.baseline, aes(x=Hb, col=site)) + geom_density()
dat.baseline[,c(10,7)] %>% split(.$gender) %>% map(summary)
dat.baseline[,c(1,7)] %>% split(.$site) %>% map(summary)

ggplot(dat.baseline, aes(x=age, y=Hb)) + geom_point() + geom_smooth(method="lm")
ggplot(dat.baseline, aes(x=weight, y=Hb)) + geom_point() + geom_smooth(method="lm")
tmp=lm(Hb~weight, data=dat.baseline)
summary(tmp)
tmp=lm(Hb~age, data=dat.baseline)
summary(tmp)



#------- LONGITUDINAL
#--------- create grouped data object
dat.grp = groupedData(data=dat.fullHB, Hb ~ pday|arm/site/pid, outer = ~gender * weight  )
head(dat.grp)

#individual profiles
xyplot(Hb ~ pday, dat.fullHB, type="l", groups = pid) 

#by study arm
xyplot(Hb ~ pday|arm, dat.fullHB, type="l", groups = pid) 

#by study site
xyplot(Hb ~ pday|site, dat.fullHB, type="l", groups = pid) 


#mean profiles
plot(dat.grp, collapse=1) #across study arm
plot(dat.grp, collapse=2) #across site
plot(dat.grp, inner=~gender, collapse=1)
plot(dat.grp, inner=~weight_quart, collapse=1)
plot(dat.grp, inner=~age_quart, collapse=1)

par(mfrow=c(1,1))
interaction.plot(dat.fullHB$pday, dat.fullHB$arm, dat.fullHB$Hb, las=1)

#same but in ggplot

tmp = dat.grp[!c(is.na(dat.grp$weight_quart)),] #data subset with no missing weight measurements
#study arm
mean_Hbs <- dat.fullHB %>%
  group_by(pday, arm) %>%
  summarize(mean_value = mean(Hb))
ggplot(mean_Hbs, aes(x = pday, y = mean_value, color = arm)) +
  geom_line() +
  labs(x = "pday", y = "Hb", color = "Arm") +
  scale_color_discrete(name = "Arm") + theme_light()
#weight_quart
mean_Hbs <- tmp %>%
  group_by(pday, weight_quart) %>%
  summarize(mean_value = mean(Hb))
ggplot(mean_Hbs, aes(x = pday, y = mean_value, color = weight_quart)) +
  geom_line() +
  labs(x = "pday", y = "Hb", color = "Weight quartile") +
  scale_color_discrete(name = "Weight quartile") + theme_light()
#age_quart
mean_Hbs <- tmp %>%
  group_by(pday, age_quart) %>%
  summarize(mean_value = mean(Hb))
ggplot(mean_Hbs, aes(x = pday, y = mean_value, color = age_quart)) +
  geom_line() +
  labs(x = "pday", y = "Hb", color = "Age quartile") +
  scale_color_discrete(name = "Age quartile") + theme_light()

#assess correlation between age and weight
cor(tmp$age, tmp$weight)
plot(tmp$age, tmp$weight)


#variance profiles
plot(dat.grp, collapse=1, FUN=function(x) sqrt(var(x)), ylab="sd(Hb)") #by arm
#plot(dat.grp, inner=~gender, FUN=function(x) sqrt(var(x)), ylab="sd(Hb)") 


#==================== MODELLING   ==============
dat.baseline$age2 = dat.baseline$age^2
#create age-relative weight
#regress weight against age and age^2 (baseline)
lm.weightadjust = lm(weight~age+age2, dat.baseline)
summary(lm.weightadjust)
#plot(lm.weightadjust) #residuals are not great, double check implications of this
#extract residuals
adjusted.weight = lm.weightadjust$residuals

#bind with pids (MOM2004_103 is missing a weight)
tst = names(adjusted.weight)
integer_vector <- as.numeric(tst)
# Find the missing integer
missing_integer <- which(diff(integer_vector) != 1)
#insert the missing NA
adjusted.weight = c(adjusted.weight[1:missing_integer], NA, adjusted.weight[(missing_integer + 1):length(integer_vector)])
names(adjusted.weight)[170] = '170'
adjusted.weight = data.frame(pid=dat.baseline$pid, adjusted.weight)

#add to dataframe
dat = merge(dat, adjusted.weight, by="pid")
dat.baseline = merge(dat.baseline, adjusted.weight, by="pid")

#Visualise decorrelation of age and weight
ggplot(dat.baseline, aes(x=age, y=weight)) + geom_point() + geom_smooth(method="lm")
ggplot(dat.baseline, aes(x=age, y=adjusted.weight)) + geom_point() + geom_smooth(method="lm")


#remember to center AGE for interpretable intercept
dat$age.centered = scale(dat$age, scale=F)
#add squared day term
dat$pday2 = dat$pday^2

#REMOVE participant with no weight measurement
index.tmp = which(is.na(dat$weight))
dat = dat[-index.tmp,]
n_distinct(dat$pid) #sample size reduced by 1

#write.csv(dat, file="dat.csv")
#==================== MODEL DATA
dat = read.csv("dat.csv")
dat$arm = as.factor(dat$arm)
dat$gender = as.factor(dat$gender)

#==================== MODEL 2 - incorporate splines, fixed = study arm and day, random = intercept|pid
#remove all missing Hb
dat.mod2 = dat
index.na=which(is.na(dat.mod2$Hb))
dat.mod2 = dat.mod2[-index.na,]
dat.mod2 = groupedData(data=dat.mod2, Hb ~ pday|site/pid  )

#using b-splines ------- only random effect on pid
library(splines)
b = bs(dat.mod2$pday, knots=7, degree = 2, intercept = FALSE)
summary(b)

mod2 = lme(data=dat.mod2, Hb ~ arm+ bs(pday, knots=7, degree=2, intercept=FALSE), random=~1|pid,
             na.action = na.omit , method = "ML")
mod2.reml = update(mod2, method="REML")

summary(mod2.reml)
intervals(mod2.reml)
plot(mod2.reml)
densityplot(as.numeric(mod2.reml$residuals))

#check augPred for generating fitted values - two pids from different arms, with different sexes
tst=augPred(mod2.reml, level=0:1)
tmp=tst %>% filter(.$.groups=="MOB2004_002" | .$.groups=="MOB2004_017")
plot(tmp)

ggemmeans(mod2, terms=c("pday [all]"))
ggeffect(mod2, terms=c("pday [all]")) 
ggpredict(mod2, terms=c("pday [all]")) 

ggemmeans(mod2, terms="pday [all]") %>% plot()
ggemmeans(mod2, terms=c("pday [all]", "arm"), type = "fe") %>% plot()

ggemmeans(mod2, terms="pday [0:42 by=1]") %>% plot() #smooth function



#==================== MODEL 3 - model 2 but with arm-day interaction
mod3 = lme(data=dat.mod2, Hb ~ arm*bs(pday, knots=7, degree=2), 
           random=~1|pid, na.action = na.omit , method = "ML")
mod3.reml = update(mod3, method="REML")

summary(mod3)
intervals(mod3)
plot(mod3)
densityplot(as.numeric(mod3$residuals))

tst=augPred(mod3, level=0:1)
tmp=tst %>% filter(.$.groups=="MOB2004_002" | .$.groups=="MOB2004_017")
plot(tmp)

ggemmeans(mod3, terms=c("pday [all]"))
ggeffect(mod3, terms=c("pday [all]")) 
ggpredict(mod3, terms=c("pday [all]")) 

ggemmeans(mod3, terms="pday [all]") %>% plot()
ggemmeans(mod3, terms=c("pday [all]", "arm"), type = "fe") %>% plot()

#compare models 2 and 3
anova(mod2, mod3)  #interaction terms not necessary

#==================== MODEL 4 - adjusted for covariates
mod4 = update(mod2, fixed=Hb~ arm + bs(pday, knots=7, degree=2) +age.centered+adjusted.weight+gender)
mod4.reml = update(mod4, method="REML")
mod4b = update(mod2, fixed=Hb~ arm * bs(pday, knots=7, degree=2) +age.centered+adjusted.weight+gender)

summary(mod4.reml)
intervals(mod4.reml)
plot(mod4)
densityplot(as.numeric(mod4$residuals))

#check augPred for generating fitted values
tst=augPred(mod4.reml, level=0:1)
tmp=tst %>% filter(.$.groups=="MOB2004_002" | .$.groups=="MOB2004_017")
plot(tmp)

ggemmeans(mod4, terms="pday [all]") %>% plot()
ggemmeans(mod4, terms=c("pday [all]", "gender", "arm"), type = "fe") %>% plot()
ggemmeans(mod4b, terms="pday [all]") %>% plot()
ggemmeans(mod4b, terms=c("pday [all]", "gender", "arm"), type = "fe") %>% plot()

#compare adjusted models to mod2
anova(mod2, mod4, mod4b) #model is significantly improved when adjusted for covariates

#==================== MODEL 5 - model4 but with day-sex-arm interaction
mod5 = update(mod4, fixed=Hb~ arm * bs(pday, knots=7, degree=2)*gender +age.centered+adjusted.weight)
mod5.reml = update(mod5, method="REML")
summary(mod5.reml)
intervals(mod5.reml)
plot(mod5)
densityplot(as.numeric(mod5$residuals))

tst=augPred(mod5, level=0:1)
tmp=tst %>% filter(.$.groups=="MOB2004_002" | .$.groups=="MOB2004_017")
plot(tmp)

ggemmeans(mod5, terms="pday [all]") %>% plot()
ggemmeans(mod5, terms=c("pday [all]", "arm"), type = "fe") %>% plot()
ggemmeans(mod5, terms = c("pday [all]", "gender", "arm")) %>% plot()


anova(mod4, mod5) #model fit not improved


#==================== MODEL 6 - RANDOM EFFECT ON SITE
mod4 = update(mod4, method="REML") #now fixed effects have been decided on, use REML
mod6 = update(mod4, random=~1|site/pid)
summary(mod6)

intervals(mod6)
anova(mod4, mod6) #model fit not improved

#visualise random effects for pids from different sites
tst=augPred(mod6, level=0:2)
tmp=tst %>% filter(.$.groups=="Boane/MOB2004_002" | .$.groups=="Namaacha/MON2003_089")
plot(tmp)

#==================== MODEL 7 - PID RE ON SPLINE
mod7 = update(mod4.reml, random=~bs(pday, knots=7, degree=2)|pid)
summary(mod7)
intervals(mod7)


anova(mod4, mod6, mod7) #model fit improved by RE on slope terms

#visualise the effect of adding random effects to slope and quadratic terms
tst=augPred(mod7, level=0:1)
tmp=tst %>% filter(.$.groups=="MOB2004_005" | .$.groups=="MOM2004_194")
plot(tmp)
tmp=tst %>% filter(.$.groups=="MOB2004_005" | .$.groups=="MOB2004_017")
plot(tmp)
tmp=tst %>% filter(.$.groups=="MOB2004_005" | .$.groups=="MOB2004_074")
plot(tmp)
tmp=tst %>% filter(.$.groups=="MOB2004_005" | .$.groups=="MOB2004_064" | .$.groups=="MOM2004_194"
                   | .$.groups=="MOB2004_074")
plot(tmp)



ggeffect(mod7, terms=c("pday [all]", "arm")) 
ggemmeans(mod7, terms=c("pday [all]", "arm")) 
ggpredict(mod7, terms=c("pday [all]", "arm")) 

#marginal means (adjusted for mean age.centered and mean adjusted.weight; average over sexes)
ggemmeans(mod7, terms = c("pday [all]", "arm")) %>% plot()

ggemmeans(mod7, terms = c("pday [all]", "arm", "gender")) %>% plot()



#==================== MODEL 7 DIAGNOSTICS
plot(mod7) #fitted vs residuals
densityplot(as.numeric(mod7$residuals))
plot(fitdist(as.numeric(mod7$residuals), distr="norm"))
qqnorm(mod7)  
plot(mod7, pid~resid(.), abline=0)

#fitted vs resids
plot(mod7, resid(., type="p")~fitted(.)|arm)
plot(mod7, resid(., type="p")~fitted(.)|pday)
plot(mod7, resid(., type="p")~fitted(.)|gender)
plot(mod7, resid(., type="p")~fitted(.)|weight_quart)
plot(mod7, resid(., type="p")~fitted(.)|age_quart)

#resids vs covariates
plot(mod7, resid(.)~pday)
plot(mod7, resid(.)~age)
plot(mod7, resid(.)~adjusted.weight)
plot(mod7, resid(.)~as.numeric(arm), xlab="arm")
plot(mod7, resid(.)~as.numeric(gender), xlab="sex")

plot(mod7, Hb~fitted(.), ylab="Observed") #high association between fitted and observed


#------- random effects
rf.m7 = ranef(mod7)
boxplot(ranef(mod7))

densityplot(rf.m7$`(Intercept)`, xlab="Intercept")
densityplot(rf.m7$`bs(pday, knots = 7, degree = 2)1`, xlab="Spline term 1")
densityplot(rf.m7$`bs(pday, knots = 7, degree = 2)2`, xlab="Spline term 2")
densityplot(rf.m7$`bs(pday, knots = 7, degree = 2)3`, xlab="Spline term 3")



#==================== MODEL 8 - model variance
#different variance for arms
vf1 = varIdent(form=~1|arm)
mod8a = update(mod7, weights=vf1)
summary(mod8a)
plot(mod8a)
plot(mod8a, resid(., type="p")~fitted(.)|arm)
anova(mod7, mod8a) #not improved

#as power of fitted values, conditional on day
vf2 = varPower(form=~fitted(.)|pday)
mod8b = update(mod7, weights=vf2)
summary(mod8b)
plot(mod8b)
plot(mod7, resid(., type="p")~fitted(.)|pday) #doesnt seem to have removed the association
plot(mod8b, resid(., type="p")~fitted(.)|pday) #doesnt seem to have removed the association


plot(mod7, resid(.)~pday)
plot(mod8b, resid(.)~pday)

anova(mod7, mod8b) #not improved



#==================== MODEL 9 - model correlation
mod9a = update(mod7, correlation=corCompSymm())
summary(mod9a)
mod9b = update(mod7, correlation=corAR1())
summary(mod9b)
anova(mod7, mod8b,mod9a, mod9b)

plot(mod9a)
plot(mod9b)



# compare model 7 to the simpler cubic polynomial model
mod7.ml = update(mod7, method="ML")
mod.tst = update(mod7.ml, fixed=Hb~ arm + pday + I(pday^2) + I(pday^3) +age.centered+adjusted.weight+gender,
                 random=~(pday+pday2)|pid)
anova(mod7.ml, mod.tst) #IC are worsened, LR is improved - mixed results
summary(mod.tst)
tst=augPred(mod.tst, level=0:1)
tmp=tst %>% filter(.$.groups=="MOB2004_005" | .$.groups=="MOM2004_194")
plot(tmp) #trend

ggemmeans(mod.tst, terms = c("pday [all]", "gender", "arm")) %>% plot()
ggemmeans(mod7.ml, terms = c("pday [all]", "gender", "arm")) %>% plot()
