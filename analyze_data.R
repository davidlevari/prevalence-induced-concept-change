# Prevalence-induced concept change in human judgment: Code for analysis and figures
# David Levari, david.levari at gmail dot com
# last updated: 2 May 2018

# How to use this file ----------------------------------------------------

# Start a fresh R session by clicking on the ".Rproj" file 
# located in the same directory as this file. Then, open this file, 
# and the code below should work.

# Please note that the GLMMs (with the glmer function) can take a long time to 
# run on some machines. I have commented out some secondary diagnostic checks in 
# the script which take even longer.

# TO DO
# Add legends back in
# fix x-axis?

# Setup -------------------------------------------------------------------

# Packages needed to run code
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
} # function to load the packages or install them if not present
packages <- c("foreign", "car", "ggplot2", "reshape", "lattice",
              "plyr", "lme4", "boot","minpack.lm","piecewiseSEM")
ipak(packages) # install packages

options(digits=3)
set.seed(101) 

# cbPalette <- c("#000000", "#839192") # set a black and gray palette for plots
cbPalette <- c("#3b92c5", "#f6bb0b") # set a blue and yellow palette for plots

theme_science <- function () { 
  theme_bw(base_family="Helvetica") 
  theme(
    axis.ticks.x = element_blank(),
    legend.key = element_blank(),
    legend.background = element_blank(),
    panel.grid = element_blank()
  )
}

# define complementary error function for CDF curve-fitting in plots
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE) 

# Study 1: Decreasing prevalence of colors --------------------------------

# IMPORT DATA
study1 <- read.csv("study1.csv", head=TRUE, na.strings="NA") #load data file
study1 <- study1[!(study1$trial %in% c(0:10)),] #remove practice trials
study1$color <- (study1$blue-154) #create zeroed color variable
study1$trial <- (study1$trial-10) #zero trial var after practice trials removed
study1$response <- as.factor(study1$response-1) #zero response var
study1$responsenum <- as.numeric(study1$response)-1 #save response var as numeric

# RECODE VARS AS FACTORS FOR GLMM
study1$participant <- as.factor(study1$participant) #treat subject as factor
study1$condition <- as.factor(study1$condition) #treat condition as factor
study1$condition <- revalue(study1$condition, 
                            c("0"="Stable Prevalence Condition", 
                              "1"="Decreasing Prevalence Condition"))

# CREATE BINNED TIME VARIABLE FOR GRAPHS
study1$timebin5 <- as.factor(apply(as.matrix(study1$trial), 2, 
                                   cut, c(seq(0,1000,200)), labels=FALSE)) #5 time bins

# BASIC DESCRIPTIVE STATS
length(unique(study1$participant)) # number of participants
length(unique(study1$participant[study1$gender=="male"])) # gender breakdown
length(unique(study1$participant[study1$gender=="female"])) # gender breakdown
mean(study1[study1$trial==1,]$age) # mean age
sd(study1[study1$trial==1,]$age) # SD of age
length(unique(study1$trial)) # how many trials per participant

# EXCLUSIONS
study1exclude <- c(116) # participants to exclude
study1 <- study1[!(study1$participant %in% study1exclude),] # exclude participants

# STUDY 1 FIGURE

# create means table for plotting
study1.means <- ddply(study1, .(condition, timebin5, color), summarise, 
                    meanr=mean(responsenum,na.rm=TRUE),.drop=FALSE) 

plotStudy1 <- ggplot(study1.means[(study1.means$timebin5==1)|
                                     (study1.means$timebin5==5),], 
                      aes(x=color, y = meanr, colour = timebin5)) + 
  facet_wrap(~condition) +
  geom_point(alpha=.7) +
  geom_line(stat="smooth",method = "nls", method.args = list(
    formula = y  ~ 0.5 * erfc( - (x - mu) / (sig * sqrt(2))),
    start=list(mu=40, sig=.5)), se=FALSE,size=1,alpha=.7) +
  xlab("Objective color") + 
  ylab("% Dots identified as blue") +
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), 
                      name=NULL,values=cbPalette) +
  scale_x_continuous(breaks = c(10,90),
                     label = c("Very purple","Very blue")) + 
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() + theme_science()
plotStudy1

### FIT A GENERALIZED LINEAR MIXED MODEL (GLMM)

# Rescale predictors between 0 and 1 for GLMM
study1zeroed <- study1
study1zeroed$color <- study1zeroed$color/100
study1zeroed$trial <- study1zeroed$trial/1000

# MODEL 1: full model with random intercepts for participant and random slopes for trial
study1.model1 <- glmer(response ~ (condition + trial + color)^3 + (trial | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study1zeroed)
summary(study1.model1)

# MODEL 1A: without 3-way interaction
study1.model1a <- glmer(response ~ condition + trial + color +
                          condition*trial +
                          condition*color +
                          trial*color + (trial | participant),
                        family = binomial,
                        control=glmerControl(optimizer="bobyqa"),
                        data = study1zeroed)
# LRT compared to model with 3-way interaction
anova(study1.model1a,study1.model1) # maximal - 1

# MODEL 2: without random slopes
study1.model2 <- glmer(response ~ (condition + trial + color)^3 + (1 | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study1zeroed)
# LRT compared to model with random slopes
anova(study1.model2, study1.model1)

# MODEL 3: without random intercepts
study1.model3 <- glmer(response ~ (condition + trial + color)^3 + (trial - 1| participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study1zeroed)
# LRT compared to model with random intercepts
anova(study1.model3, study1.model1) 

# Study 2: Decreasing prevalence of colors with forewarning ---------------

# IMPORT DATA
study2 <- read.csv("study2.csv", head=TRUE, na.strings="NA") #load data file
study2 <- study2[!(study2$trial %in% c(0:10)),] #remove practice trials
study2$color <- (study2$blue-154) #create zeroed color variable
study2$trial <- (study2$trial-10) #zero trial var after practice trials removed
study2$response <- as.factor(study2$response-1) #zero response var
study2$responsenum <- as.numeric(study2$response)-1 #save response var as numeric

# RECODE VARS AS FACTORS FOR GLMM
study2$participant <- as.factor(study2$participant) #treat subject as factor
study2$condition <- as.factor(study2$condition) #treat condition as factor
study2$condition <- revalue(study2$condition, 
                            c("0"="Stable Prevalence with Warning Condition", 
                              "1"="Decreasing Prevalence with Warning Condition"))

# CREATE BINNED TIME VARIABLE FOR GRAPHS
study2$timebin5 <- as.factor(apply(as.matrix(study2$trial), 2, 
                                   cut, c(seq(0,1000,200)), labels=FALSE)) #5 time bins

# BASIC DESCRIPTIVE STATS
length(unique(study2$participant)) # number of participants
length(unique(study2$participant[study2$gender=="male"])) # gender breakdown
length(unique(study2$participant[study2$gender=="female"])) # gender breakdown
mean(study2[study2$trial==1,]$age) # mean age
sd(study2[study2$trial==1,]$age) # SD of age
length(unique(study2$trial)) # how many trials per participant

# EXCLUSIONS
study2exclude <- c(101,104,113,131) # participants to exclude
study2 <- study2[!(study2$participant %in% study2exclude),] # exclude participants

# STUDY 2 FIGURE

# create means table for plotting
study2.means <- ddply(study2, .(condition, timebin5, color), summarise, 
                    meanr=mean(responsenum,na.rm=TRUE),.drop=FALSE) 

plotStudy2 <- ggplot(study2.means[(study2.means$timebin5==1)|
                                    (study2.means$timebin5==5),], 
                     aes(x=color, y = meanr, colour = timebin5)) + 
  facet_wrap(~condition) +
  geom_point(alpha=.7) +
  geom_line(stat="smooth",method = "nls", method.args = list(
    formula = y  ~ 0.5 * erfc( - (x - mu) / (sig * sqrt(2))),
    start=list(mu=40, sig=.5)), se=FALSE,size=1,alpha=.7) +
  xlab("Objective color") + 
  ylab("% Dots identified as blue") +
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), 
                      name=NULL,values=cbPalette) +
  scale_x_continuous(breaks = c(10,90),
                     label = c("Very purple","Very blue")) + 
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() + theme_science() 
plotStudy2

### FIT A GENERALIZED LINEAR MIXED MODEL (GLMM)

# Rescale predictors between 0 and 1 for GLMM
study2zeroed <- study2
study2zeroed$color <- study2zeroed$color/100
study2zeroed$trial <- study2zeroed$trial/1000

# MODEL 1: full model with random intercepts for participant and random slopes for trial
study2.model1 <- glmer(response ~ (condition + trial + color)^3 + (trial | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study2zeroed)
summary(study2.model1)

# MODEL 1A: Model 1 without 3-way interaction
study2.model1a <- glmer(response ~ condition + trial + color +
                          condition*trial +
                          condition*color +
                          trial*color + (trial | participant),
                        family = binomial,
                        control=glmerControl(optimizer="bobyqa"),
                        data = study2zeroed)
# LRT compared to model without 3-way interaction
anova(study2.model1a,study2.model1) # maximal - 1

# MODEL 2: Model without random slopes
study2.model2 <- glmer(response ~ (condition + trial + color)^3 + (1 | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study2zeroed)
# LRT compared to model without random slopes
anova(study2.model2, study2.model1) 

# MODEL 3: Model without random intercepts
study2.model3 <- glmer(response ~ (condition + trial + color)^3 + (trial - 1| participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study2zeroed)
# LRT compared to model without random intercepts
anova(study2.model3, study2.model1) 

# Study 3: Decreasing prevalence of colors with incentives  --------

### IMPORT DATA
study3 <- read.csv("study3.csv", head=TRUE, na.strings="NA") #load data file
study3$color <- (study3$color-154) #create zeroed color variable
study3$response <- as.factor(study3$response) #zero response var
study3$responsenum <- as.numeric(study3$response)-1 #save response var as numeric
study3 <- study3[!(study3$trial %in% c(1001:1010)),] #remove practice trials

# RECODE VARS AS FACTORS FOR GLMM
study3$participant <- as.factor(study3$participant) #treat subject as factor
study3$condition <- as.factor(study3$condition) #treat condition as factor
study3$condition <- revalue(study3$condition, 
                            c("0"="Stable Prevalence Condition", 
                              "1"="Decreasing Prevalence Condition",
                              "2"="Decreasing Prevalence + Instruction Condition",
                              "3"="Decreasing Prevalence + Instruction + Incentive Condition"
                            ))

# CREATE BINNED TIME VARIABLE FOR GRAPHS
study3$timebin4 <- as.factor(apply(as.matrix(study3$trial), 2,
                                   cut, c(seq(0,800,200)), labels=FALSE)) #4 time bins

# BASIC DESCRIPTIVE STATS
length(unique(study3$participant)) # number of participants
length(unique(study3$participant[study3$gender=="Male"])) # gender breakdown
length(unique(study3$participant[study3$gender=="Female"])) # gender breakdown
mean(study3[study3$trial==1,]$age,na.rm=TRUE) # mean age
sd(study3[study3$trial==1,]$age,na.rm=TRUE) # SD of age
length(unique(study3$trial)) # how many trials per participant

# EXCLUSIONS
study3exclude <- c(164) # participants to exclude
study3 <- study3[!(study3$participant %in% study3exclude),] # exclude participants

# STUDY 3 FIGURE
study3.means <- ddply(study3, .(condition, timebin4, color), summarise, 
                    meanr=mean(responsenum,na.rm=TRUE),.drop=FALSE) # create means table for plotting

plotStudy3 <- ggplot(study3.means[(study3.means$timebin4==1)|
                                    (study3.means$timebin4==4),], 
                     aes(x=color, y = meanr, colour = timebin4)) + 
  facet_wrap(~condition) +
  geom_point(alpha=.7) +
  geom_line(stat="smooth",method = "nls", method.args = list(
    formula = y  ~ 0.5 * erfc( - (x - mu) / (sig * sqrt(2))),
    start=list(mu=40, sig=.5)), se=FALSE,size=1,alpha=.7) +
  xlab("Objective color") + 
  ylab("% Dots identified as blue") +
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), 
                      name=NULL,values=cbPalette) +
  scale_x_continuous(breaks = c(10,90),
                     label = c("Very purple","Very blue")) + 
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() + theme_science()
plotStudy3

### FIT A GENERALIZED LINEAR MIXED MODEL (GLMM)

# RESCALE DATASET BETWEEN 0 AND 1 FOR GLMM
study3zeroed <- study3
study3zeroed$color <- study3zeroed$color/100
study3zeroed$trial <- study3zeroed$trial/800

# set reference level for condition variable
study3zeroed <- within(study3zeroed, condition <- relevel(condition,                                                     
                        ref = "Stable Frequency"))

# MODEL 1: full model with random intercepts for participant and random slopes for trial
study3.model1 <- glmer(response ~ (condition + trial + color)^3 + (trial | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study3zeroed)
summary(study3.model1)

# change reference level for condition variable and rerun model
study3zeroed <- within(study3zeroed, condition <- relevel(condition, 
                        ref = "Decreasing Prevalence + Instruction + Incentive Condition"))
study3.model1.releveled1 <- glmer(response ~ (condition + trial + color)^3 + 
                                    (trial | participant),
                                  family = binomial,
                                  control=glmerControl(optimizer="bobyqa"),
                                  data = study3zeroed)
summary(study3.model1.releveled1)

# change reference level for condition variable and rerun model
study3zeroed <- within(study3zeroed, condition <- relevel(condition, 
                        ref = "Decreasing Prevalence + Instruction Condition"))
study3.model1.releveled2 <- glmer(response ~ (condition + trial + color)^3 + 
                                    (trial | participant),
                                  family = binomial,
                                  control=glmerControl(optimizer="bobyqa"),
                                  data = study3zeroed)
summary(study3.model1.releveled2)

# collect releveled model p-values and correct for multiple comparisons
study3.model1.pvalues <- c(coef(summary(study3.model1))[14:16,4],
                      coef(summary(study3.model1.releveled1))[14:16,4],
                      coef(summary(study3.model1.releveled2))[14:16,4])
p.adjust(study3.model1.pvalues,method="holm")

# MODEL 1A: without 3-way interaction
study3.model1a <- glmer(response ~ condition + trial + color +
                          condition*trial +
                          condition*color +
                          trial*color + (trial | participant),
                        family = binomial,
                        control=glmerControl(optimizer="bobyqa"),
                        data = study3zeroed)
# LRT compared to model without 3-way interaction
anova(study3.model1a,study3.model1) # maximal - 1

# MODEL 2: without random slopes
study3.model2 <- glmer(response ~ (condition + trial + color)^3 + 
                         (1 | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study3zeroed)
# LRT compared to model without random slopes
anova(study3.model2, study3.model1) 

# MODEL 3: without random intercepts
study3.model3 <- glmer(response ~ (condition + trial + color)^3 + 
                         (trial - 1| participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study3zeroed)
# LRT compared to model with random intercepts
anova(study3.model3, study3.model1) 


# Study 4: Abrupt or gradual decreasing prevalence of colors --------------

# IMPORT DATA
study4 <- read.csv("study4.csv", head=TRUE, na.strings="NA") #load data file
study4$color <- (study4$blue-154) #create zeroed color variable
study4$response <- as.factor(study4$response-1) #zero response var
study4$responsenum <- as.numeric(study4$response)-1 #save response var as numeric

# RECODE VARS AS FACTORS FOR GLMM
study4$participant <- as.factor(study4$participant) #treat subject as factor
study4$condition <- as.factor(study4$condition) #treat condition as factor
study4$condition <- revalue(study4$condition, 
                            c("0"="Stable Prevalence Condition", 
                              "1"="Gradually Decreasing Prevalence Condition",
                              "2" = "Abruptly Decreasing Prevalence Condition"))

# CREATE BINNED TIME VARIABLE FOR GRAPHS
study4$timebin4 <- as.factor(apply(as.matrix(study4$trial), 2, 
                                   cut, c(seq(0,800,200)), labels=FALSE)) #5 time bins

# BASIC DESCRIPTIVE STATS
length(unique(study4$participant)) # number of participants
length(unique(study4$participant[study4$gender=="male"])) # gender breakdown
length(unique(study4$participant[study4$gender=="female"])) # gender breakdown
mean(study4[study4$trial==1,]$age) # mean age
sd(study4[study4$trial==1,]$age) # SD of age
length(unique(study4$trial)) # how many trials per participant

# EXCLUSIONS
# study4exclude <- c() # no exclusions
# study4 <- study4[!(study4$participant %in% study4exclude),] 

# STUDY 4 FIGURE

# create means table for plotting
study4.means <- ddply(study4, .(condition, timebin4, color), summarise, 
                      meanr=mean(responsenum,na.rm=TRUE),.drop=FALSE) 
plotStudy4 <- ggplot(study4.means[(study4.means$timebin4==1)|
                                    (study4.means$timebin4==4),], 
                     aes(x=color, y = meanr, colour = timebin4)) + 
  facet_wrap(~condition) +
  geom_point(alpha=.7) +
  geom_line(stat="smooth",method = "nls", method.args = list(
    formula = y  ~ 0.5 * erfc( - (x - mu) / (sig * sqrt(2))),
    start=list(mu=40, sig=.5)), se=FALSE,size=1,alpha=.7) +
  xlab("Objective color") + 
  ylab("% Dots identified as blue") +
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), 
                      name=NULL,values=cbPalette) +
  scale_x_continuous(breaks = c(10,90),
                     label = c("Very purple","Very blue")) + 
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() + theme_science()
plotStudy4

# FIT A GENERALIZED LINEAR MIXED MODEL (GLMM)

# Rescale predictors between 0 and 1 for GLMM
study4zeroed <- study4
study4zeroed$color <- study4zeroed$color/100
study4zeroed$trial <- study4zeroed$trial/800

# MODEL 1: full model with random intercepts for participant and random slopes for trial
study4.model1 <- glmer(response ~ (condition + trial + color)^3 + (trial | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study4zeroed)
summary(study4.model1)

# change reference level for condition variable and rerun model
study4zeroed <- within(study4zeroed, condition <- relevel(condition, ref = 3))
study4.model1.releveled <- glmer(response ~ (condition + trial + color)^3 + 
                                   (trial | participant),
                                 family = binomial,
                                 control=glmerControl(optimizer="bobyqa"),
                                 data = study4zeroed)
summary(study4.model1.releveled)

# collect releveled model p-values and correct for multiple comparisons
study4.model1.pvalues <- c(coef(summary(study4.model1))[11:12,4],
                      coef(summary(study4.model1.releveled))[11:12,4])
p.adjust(study4.model1.pvalues,method="holm")

# MODEL 1A: without 3-way interaction
study4.model1a <- glmer(response ~ condition + trial + color +
                          condition*trial +
                          condition*color +
                          trial*color + (trial | participant),
                        family = binomial,
                        control=glmerControl(optimizer="bobyqa"),
                        data = study4zeroed)
# LRT compared to model without 3-way interaction
anova(study4.model1a,study4.model1) # maximal - 1

# MODEL 2: without random slopes
study4.model2 <- glmer(response ~ (condition + trial + color)^3 + (1 | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study4zeroed)
# LRT compared to model without random slopes
anova(study4.model2, study4.model1) 

# MODEL 3: without random intercepts
study4.model3 <- glmer(response ~ (condition + trial + color)^3 + 
                         (trial - 1| participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study4zeroed)
# LRT compared to model without random intercepts
anova(study4.model3, study4.model1) 

# Study 5: Increasing prevalence of colors --------------------------------

# IMPORT DATA
study5 <- read.csv("study5.csv", head=TRUE, na.strings="NA") #load data file
study5 <- study5[!(study5$trial %in% c(0:10)),] #remove practice trials
study5$color <- (study5$blue-154) #create zeroed color variable
study5$trial <- (study5$trial-10) #zero trial var after practice trials removed
study5$response <- as.factor(study5$response-1) #zero response var
study5$responsenum <- as.numeric(study5$response)-1 #save response var as numeric

# RECODE VARS AS FACTORS FOR GLMM
study5$participant <- as.factor(study5$participant) #treat subject as factor
study5$condition <- as.factor(study5$condition) #treat condition as factor
study5$condition <- revalue(study5$condition, 
                            c("0"="Stable Prevalence Condition", 
                              "1"="Increasing Prevalence Condition"))

# CREATE BINNED TIME VARIABLE FOR GRAPHS
study5$timebin5 <- as.factor(apply(as.matrix(study5$trial), 2, 
                                cut, c(seq(0,1000,200)), labels=FALSE)) # 5 time bins

# BASIC DESCRIPTIVE STATS
length(unique(study5$participant)) # number of participants
length(unique(study5$participant[study5$gender=="male"])) # gender breakdown
length(unique(study5$participant[study5$gender=="female"])) # gender breakdown
mean(study5[study5$trial==1,]$age) # mean age
sd(study5[study5$trial==1,]$age) # SD of age
length(unique(study5$trial)) # how many trials per participant

# EXCLUSIONS
study5exclude <- c(111) # participants to exclude
study5 <- study5[!(study5$participant %in% study5exclude),] # exclude participants

# STUDY 5 FIGURE

# create means table for plotting
study5.means <- ddply(study5, .(condition, timebin5, color), summarise, 
                      meanr=mean(responsenum,na.rm=TRUE),.drop=FALSE) 

plotStudy5 <- ggplot(study5.means[(study5.means$timebin5==1)|
                                    (study5.means$timebin5==5),], 
                     aes(x=color, y = meanr, colour = timebin5)) + 
  facet_wrap(~condition) +
  geom_point(alpha=.7) +
  geom_line(stat="smooth",method = "nls", method.args = list(
    formula = y  ~ 0.5 * erfc( - (x - mu) / (sig * sqrt(2))),
    start=list(mu=40, sig=.5)), se=FALSE,size=1,alpha=.7) +
  xlab("Objective color") + 
  ylab("% Dots identified as blue") +
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), 
                      name=NULL,values=cbPalette) +
  scale_x_continuous(breaks = c(10,90),
                     label = c("Very purple","Very blue")) + 
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() + theme_science()
plotStudy5

# FIT A GENERALIZED LINEAR MIXED MODEL (GLMM)

# RESCALE DATASET BETWEEN 0 AND 1 FOR GLMM
study5zeroed <- study5
study5zeroed$color <- study5zeroed$color/100
study5zeroed$trial <- study5zeroed$trial/1000

# MODEL 1: full model with random intercepts for participant and random slopes for trial
study5.model1 <- glmer(response ~ (condition + trial + color)^3 + (trial | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study5zeroed)
summary(study5.model1)

# MODEL 1A: without 3-way interaction
study5.model1a <- glmer(response ~ condition + trial + color +
                          condition*trial +
                          condition*color +
                          trial*color + (trial | participant),
                        family = binomial,
                        control=glmerControl(optimizer="bobyqa"),
                        data = study5zeroed)
# LRT compared to model without 3-way interaction
anova(study5.model1a,study5.model1) # maximal - 1

# MODEL 2: without random slopes
study5.model2 <- glmer(response ~ (condition + trial + color)^3 + (1 | participant),
                        family = binomial,
                        control=glmerControl(optimizer="bobyqa"),
                        data = study5zeroed)
# LRT compared to model without random slopes
anova(study5.model2, study5.model1) 

# MODEL 3: without random intercepts
study5.model3 <- glmer(response ~ (condition + trial + color)^3 + (trial - 1| participant),
                        family = binomial,
                        control=glmerControl(optimizer="bobyqa"),
                        data = study5zeroed)
# LRT compared to model without random intercepts
anova(study5.model3, study5.model1) 

# Study 6: Decreasing prevalence of threatening faces ---------------------

# IMPORT DATA
study6 <- read.csv("study6.csv", head=TRUE, na.strings="NA") #load data file
study6 <- study6[!(study6$trial %in% c(0:10)),] #remove practice trials
study6$trial <- (study6$trial-10) #zero trial var after practice trials removed
study6$response <- as.factor(study6$response-1) #zero response var
study6$responsenum <- as.numeric(study6$response)-1 #save response var as numeric

# RECODE VARS AS FACTORS FOR GLMM
study6$participant <- as.factor(study6$participant) #treat subject as factor
study6$condition <- as.factor(study6$condition) #treat condition as factor
study6$condition <- revalue(study6$condition, 
                            c("0"="Stable Prevalence Condition", 
                              "1"="Decreasing Prevalence Condition"))

# CREATE BINNED TIME VARIABLE FOR GRAPHS
study6$timebin4 <- as.factor(apply(as.matrix(study6$trial), 2, 
                                   cut, c(seq(0,800,200)), labels=FALSE)) #5 time bins

# BASIC DESCRIPTIVE STATS
length(unique(study6$participant)) # number of participants
length(unique(study6$participant[study6$gender=="male"])) # gender breakdown
length(unique(study6$participant[study6$gender=="female"])) # gender breakdown
mean(study6[study6$trial==1,]$age) # mean age
sd(study6[study6$trial==1,]$age) # SD of age
length(unique(study6$trial)) # how many trials per participant

# EXCLUSIONS
study6exclude <- c(125) # participants to exclude
study6 <- study6[!(study6$participant %in% study6exclude),] # exclude participants

# STUDY 6 FIGURE 
study6.means <- ddply(study6, .(condition, timebin4, face), summarise, 
                    meanr=mean(responsenum,na.rm=TRUE),.drop=FALSE) # create means table for plotting

plotStudy6 <- ggplot(study6.means[(study6.means$timebin4==1)|
                                    (study6.means$timebin4==4),], 
                     aes(x=face, y = meanr, colour = timebin4)) + 
  facet_wrap(~condition) +
  geom_point(alpha=.7) +
  geom_line(stat="smooth",method = "nls", method.args = list(
    formula = y  ~ 0.5 * erfc( - (x - mu) / (sig * sqrt(2))),
    start=list(mu=40, sig=.5)), se=FALSE,size=1,alpha=.7) +
  xlab("Objective threateningness") + 
  ylab("% Faces identified as threatening") +
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), 
                      name=NULL,values=cbPalette) +
  scale_x_continuous(breaks = c(10,50),
                     label = c("Very nonthreatening","Very threatening")) + 
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() + theme_science()
plotStudy6

# FIT A GENERALIZED LINEAR MIXED MODEL (GLMM)

# RESCALE DATASET BETWEEN 0 AND 1 FOR GLMM
study6zeroed <- study6
study6zeroed$face <- study6zeroed$face/60
study6zeroed$trial <- study6zeroed$trial/800

# full model with random intercepts for participant and random slopes for trial
study6.model1 <- glmer(response ~ (condition + trial + face)^3 + (trial | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study6zeroed)
summary(study6.model1)

# LRT compared to model without 3-way interaction
study6.model1a <- glmer(response ~ condition + trial + face +
                          condition*trial +
                          condition*face +
                          trial*face + (trial | participant),
                        family = binomial,
                        control=glmerControl(optimizer="bobyqa"),
                        data = study6zeroed)

anova(study6.model1a,study6.model1) # maximal - 1

# LRT compared to model without random slopes
  
study6.model2 <- glmer(response ~ (condition + trial + face)^3 + (1 | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study6zeroed)
anova(study6.model2, study6.model1) 

# LRT compared to model without random intercepts
  
study6.model3 <- glmer(response ~ (condition + trial + face)^3 + (trial - 1| participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study6zeroed)
anova(study6.model3, study6.model1) 


# Study 7: Decreasing prevalence of unethical studies ---------------------

### IMPORT DATA
study7 <- read.csv("study7.csv", head=TRUE, na.strings="NA") #load data file
study7$responsenum <- study7$response #save response var as numeric
study7$response <- as.factor(study7$response) #zero response var

# RECODE VARS AS FACTORS FOR GLMM
study7$participant <- as.factor(study7$participant) #treat subject as factor
study7$condition <- as.factor(study7$prevalence) #treat condition as factor
study7$condition <- factor(study7$condition, levels(study7$condition)[c(2,1)])
study7$condition <- revalue(study7$condition, 
                            c("Stable Condition"="Stable Prevalence Condition", 
                              "Decreasing Condition"="Decreasing Prevalence Condition"))
study7$norm_mean <- -study7$norm_mean + 8 # reverse score pretested mean ethicality

# CREATE BINNED TIME AND INTENSITY VARIABLE FOR GRAPHS
study7$timebin5 <- factor(study7$timebin5)
study7$normbin24 <- apply(as.matrix(study7$norm_mean), 2, 
                          cut, 24, labels=FALSE)

# BASIC DESCRIPTIVE STATS
length(unique(study7$participant)) # number of participants
length(unique(study7$participant[study7$gender=="male"])) # gender breakdown
length(unique(study7$participant[study7$gender=="female"])) # gender breakdown
mean(study7[study7$trial==1,]$age,na.rm=TRUE) # mean age
sd(study7[study7$trial==1,]$age,na.rm=TRUE) # SD of age
length(unique(study7$trial)) # how many trials per participant

# EXCLUSIONS
study7exclude <- c(101,119,133,139,143,160,177,183,188) # participants to exclude
study7 <- study7[!(study7$participant %in% study7exclude),] # exclude participants

# Figure 5: Results of Study 7

# create means table for plotting
study7.means <- ddply(study7, .(condition, timebin5, normbin24), summarise, 
                    meanr=mean(responsenum,na.rm=TRUE),.drop=FALSE) 

plotStudy7 <- ggplot(study7.means[(study7.means$timebin5==1)|
                                    (study7.means$timebin5==5),], 
                     aes(x=normbin24, y = meanr, colour = timebin5)) + 
  facet_wrap(~condition) +
  geom_point(alpha=.7) +
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1,alpha=.7) +
  xlab("Objective ethicality") + 
  ylab("% Proposals identified as unethical") +
  scale_colour_manual(labels=c("Initial 48 trials", "Final 48 trials"), 
                      name=NULL,values=cbPalette) +
  scale_x_continuous(breaks = c(4,21),
                     label = c("Very ethical","Very unethical")) + 
  scale_y_continuous(labels = scales::percent) + 
  theme_bw() + theme_science()
plotStudy7

# FIT A GENERALIZED LINEAR MIXED MODEL (GLMM)

# RESCALE DATASET BETWEEN 0 AND 1 FOR GLMM
study7zeroed <- study7
study7zeroed$trial <- study7zeroed$trial/240
study7zeroed$norm_mean <- study7zeroed$norm_mean/max(study7zeroed$norm_mean)

# MODEL 1: full model with random intercepts for participant and random slopes for trial
study7.model1 <- glmer(response ~ (condition + trial + norm_mean)^3 + (trial | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study7zeroed)
summary(study7.model1)

# MODEL 1A: without 3-way interaction
study7.model1a <- glmer(response ~ condition + trial + norm_mean +
                          condition*trial +
                          condition*norm_mean +
                          trial*norm_mean + (trial | participant),
                        family = binomial,
                        control=glmerControl(optimizer="bobyqa"),
                        data = study7zeroed)
# LRT compared to model without 3-way interaction
anova(study7.model1a,study7.model1) # maximal - 1

# MODEL 2: without random slopes
study7.model2 <- glmer(response ~ (condition + trial + norm_mean)^3 + (1 | participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study7zeroed)
# LRT compared to model without random slopes
anova(study7.model2, study7.model1) 

# MODEL 3: without random intercepts
study7.model3 <- glmer(response ~ (condition + trial + norm_mean)^3 + (trial - 1| participant),
                       family = binomial,
                       control=glmerControl(optimizer="bobyqa"),
                       data = study7zeroed)
# LRT compared to model without random intercepts
anova(study7.model3, study7.model1) 

# GLMM Conditional R-Squared, all models ----------------------------------

allstudies.modelist <- list(study1.model1,
                            study2.model1,
                            study3.model1,
                            study4.model1,
                            study5.model1,
                            study6.model1,
                            study7.model1)
sem.model.fits(allstudies.modelist)

# Checks against model convergence failures in lme4 -----------------------

# source(system.file("utils", "allFit.R", package="lme4"))
# 
# # function to check models with multiple optimizers 
# # to assess severity of convergence failures
# allFits <- function(model){
#   modelname <- paste(deparse(substitute(model)),"allFit",sep=".")
#   summaryname <- paste("ss",deparse(substitute(model)),sep=".")
#   model.allFit <- allFit(model)
#   ss <- summary(model.allFit)
#   assign(modelname, model.allFit, env=.GlobalEnv)
#   assign(summaryname, summary(model.allFit), env=.GlobalEnv)
#   print(ss$fixef)
#   print(ss$which.OK)
# } 
# 
# Use multiple optimizers on models with convergence failures (all look OK)
# allFits(study2.model2)
# allFits(study3.model1)
# allFits(study3.model1.releveled1)
# allFits(study3.model1.releveled2)
# allFits(study3.model2)
# allFits(study3.model3)
# allFits(study4.model1)
# allFits(study4.model2)
# allFits(study4.model3)
# allFits(study5.model1)
# allFits(study5.model2)
# allFits(study5.model3)

# Save graphs -------------------------------------------------------------

ggsave(file="study1fig.pdf", plotStudy1,height = 2, width = 4.75,dpi=300)
ggsave(file="study2fig.pdf", plotStudy2,height = 1.58, width = 4.75,dpi=300)
ggsave(file="study3fig.pdf", plotStudy3,height = 3.17, width = 4.75,dpi=300)
ggsave(file="study4fig.pdf", plotStudy4,height = 1.58, width = 4.75,dpi=300)
ggsave(file="study5fig.pdf", plotStudy5,height = 1.58, width = 4.75,dpi=300)
ggsave(file="study6fig.pdf", plotStudy6,height = 1.58, width = 4.75,dpi=300)
ggsave(file="study7fig.pdf", plotStudy7,height = 1.58, width = 4.75,dpi=300)