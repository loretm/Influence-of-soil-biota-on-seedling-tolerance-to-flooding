##### BALD CYPRESS SEEDLING FLOODING EXPERIMENT DATA ANALYSES#######

#* Setting your working directory
setwd("~/Desktop/GrowthAnalysis_BCFII/GrowthAnalyses") # This code sets your working directory

############################################################################################# 
## PACKAGES AND LIBRARIES NEEDED 
#############################################################################################

install.packages("Rmisc")
install.packages("ggplot2", dependencies=TRUE)
install.packages ("grid")
install.packages ("gridExtra")
install.packages ("grid")
install.packages("MASS")
install.packages("lme4")
install.packages("agridat")
install.packages("car")
install.packages("lsmeans")
install.packages("lmerTest")
install.packages("lme4")
install.packages("agricolae")
install.packages( "multcomp")
install.packages ("emmeans")
install.packages("pls")
install.packages("multcompView")
install.packages("car")  # For testing homogeneity of variances

library (car)
library(ggplot2)
library(gridExtra)
library(grid)
install.packages("ggpubr")
library(ggpubr)
install.packages("Rmisc")
library(Rmisc)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library (car) ## Needed for testing homogeneity of variances
library(agricolae) # Needed to obtain Tukey groupings
library(MASS)
library(lme4)
library(lsmeans)
library(lmerTest)
library(foreign)
library(multcomp)
library(multcompView)
library(pls) ## Needed for partial least-square analyses
library(tidyr)

############################################################################################# 
## DATA SETS TO BE USED:
#############################################################################################

#**** Plant Physiological Measurements:
P <- read.table("PlantPhysiologicalMesuarements.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 
# Creating additional variables within this dataset:
P$DryBiomassRatio <- P$ADW/P$BDW  ### Calculating Shoot (Above) to root (belowground) ratios
P$TDB <- P$ADW + P$BDW  #### Calculating Total Dry Biomass (TDB)
P$RADW <- P$ADW/P$IH    #### Calcutlating Relative Aboveground dry biomass with respect to the initial height of plant before experiment (H)   
P$RBDW <- P$BDW/P$IH    #### Calculating Relative Belowground dry biomass
P$RatioRAB <- P$RADW/P$RBDW  ### Calculating Relative Above-Belowground biomass ratio
P$TCSH <- ((P$FH - P$IH)/P$IH)
#**** Beneficial Symbiont Abundances (AMF & DSE) #########
E <- read.table("Fungal.Colonization.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 

#***** Soil Nutrient Analyses: 
Nu <- read.table("SoilNutrients.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 

#### Setting up contrasts:
options("contrasts")  # To veryfy which type of contrast is R using for the ANOVAs tables
options(contrasts=c("contr.sum","contr.poly")) ### To set up correct contrasts for Type III ANOVA's
#options(contrasts=c("contr.treatment","contr.poly"))
#To comeback to default options: backup_options <- options() and then type options(backup_options). Do this before changing any options when you start R


################################################################################################################################################
## TESTING ASSUMPTIONS OF NORMALITY & HETEROSCEDASCITY
###############################################################################################################################################


#****************************************************************************
#             GROWTH VARIABLES: CHANGE IN HEIGHT, SLA
#****************************************************************************
#---------------------------- Relative Total Change in Height (TCSH):
hist(P$TCSH, probability = T, main= "Histogram of normal data - Relative Total Change in Height", xlab="Approximately normally distributed data")
lines(density(P$TCSH), col=2)
qqnorm(P$TCSH,main="QQ plot of normal data- Relative Total Change in Height",pch=19) 
qqline(P$TCSH) 
 #* Testing Normality assuming Mixed-Model:
lmmTCSH <- lmer(TCSH ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P, REML = TRUE) #Linear Mixed Model
Anova(lmmTCSH, type = "III")  
hist(resid(lmmTCSH,type="deviance"))
qqnorm(resid(lmmTCSH,type="deviance"))
qqline(resid(lmmTCSH,type="deviance"))   # Distribution of residuals is not normal, we need to transform data
shapiro.test(resid(lmmTCSH))
#* Testing Normality assuming Mixed-Model & log-transformed data:
P$logTCSH <- log(P$TCSH + 1)
lmmTCSH <- lmer(logTCSH ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P, REML = TRUE) #Linear Mixed Model
hist(resid(lmmTCSH,type="deviance"))
qqnorm(resid(lmmTCSH,type="deviance"))
qqline(resid(lmmTCSH,type="deviance"))   # Distribution of residuals is normal, proceed with analyses
shapiro.test(resid(lmmTCSH))
Anova(lmmTCSH, type = "III")  #testing significance of fixed factors

boxplot(logTCSH ~ Inundation + Soil + Soil*Inundation, data = P)


#TCSH= Relative Total Change in Height
lmmTCSH <- lmer(logTCSH ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P, REML=TRUE) #Linear Mixed Model
Anova(lmmTCSH, type = "III")  # testing significance of fixed factors

#testing significance of random factor
lmTCSH <- lm(logTCSH ~ Inundation + Soil + Soil*Inundation, data = P)
anova(lmmTCSH,lmTCSH)

#*Posthoc tests & LSmeans
library(multcomp)
library(multcompView)
lsmeans(lmmTCSH, pairwise ~ Inundation*Soil, adjust ="") #all pairwise comparisons
lsmeans(lmmTCSH, pairwise ~ Inundation|Soil, adjust ="tukey") #within soil treatments among hydrological regimes

#*Obtaining Tukey-groupings:
#-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
TCSH.lsm <- lsmeans(lmmTCSH, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
cld(TCSH.lsm, alpha = .50)



#---------------------------- SLAFL - specific leave area of firts leave
hist(P$SLAFL, probability = T, main= "Histogram of normal data - SLAFL", xlab="Approximately normally distributed data")
lines(density(P$SLAFL), col=2)
qqnorm(P$SLAFL,main="QQ plot of normal data- SLAFL",pch=19) 
qqline(P$SLAFL) 
shapiro.test(P$SLAFL) # Test of Normality - Not normality in this data
leveneTest(SLAFL ~ Soil*Inundation, data=P) ## There is homogeneity of variance
#* Testing Normality assuming Mixed-Model:
lmmSLA <- lmer(SLAFL ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P, REML = TRUE) #Linear Mixed Model
Anova(lmmSLA, type = "III")  
hist(resid(lmmSLA,type="deviance"))
qqnorm(resid(lmmSLA,type="deviance"))
qqline(resid(lmmSLA,type="deviance"))   # Distribution of residuals is not normal, we need to transform data
shapiro.test(resid(lmmSLA))
#* Testing Normality assuming Mixed-Model & log-transformed data:            
P$LogSLA <- log(P$SLAFL)
lmmSLA <- lmer(LogSLA ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P, REML = TRUE)
hist(resid(lmmSLA,type="deviance"))
qqnorm(resid(lmmSLA,type="deviance"))
qqline(resid(lmmSLA,type="deviance"))   # Distribution of residuals is normal, proceed with analyses
shapiro.test(resid(lmmSLA))          
    
#SLA = Specific Leaf Area Mixed-model
P$LogSLA <- log(P$SLAFL)
lmmSLA <- lmer(LogSLA ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
               REML = FALSE)
Anova(lmmSLA, type = "III") #testing significance of fixed factors
plot(lmmSLA)
#*Testing random factor effect
lmmSLA1 <- lmer(LogSLA ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                REML = FALSE)
lmmSLA0 <- lm(LogSLA ~ Inundation + Soil + Soil*Inundation, data = P)
anova(lmmSLA1, lmmSLA0)


#*Posthoc tests & LSmeans
library(multcomp)
library(multcompView)
lsmeans(lmmSLA, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
lsmeans(lmmSLA, pairwise ~ Inundation|Soil, adjust ="tukey") #within soil treatments among hydrological regimes

#*Obtaining Tukey-groupings:
#-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
SLA.lsm <- lsmeans(lmmSLA, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
cld(SLA.lsm, alpha = .50)



#---------------------------- RHOAF = Change in Height after flooding

# RHOAF = Change in Height after flooding
hist(P$RHOAF, probability = T, main= "Histogram of normal data - Relative Total Change in Height", xlab="Approximately normally distributed data")
lines(density(GAF$RHOAF), col=2)
qqnorm(GAF$RHOAF,main="QQ plot of normal data- Relative Total Change in Height",pch=19) 
qqline(GAF$RHOAF) 
shapiro.test(GAF$RHOAF)
leveneTest(RHOAF ~ Soil*Inundation, data=GAF) # Testing heteroscesdacity

          #* Transforming data & testing assumptions using Linear Mixed Model:
          GAF$LogRHOAF <- log(GAF$RHOAF + 1)
          lmmRHOAF <- lmer(LogRHOAF ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = GAF,
                REML = FALSE)
          hist(resid(lmmRHOAF,type="deviance"))
          qqnorm(resid(lmmRHOAF,type="deviance"))
          qqline(resid(lmmRHOAF,type="deviance"))
          
          #* Conclusion: residuals are not normal. Proceeded with GLMMS!


# RHIN = Change in height after inoculation
hist(GAF$RHIN, probability = T, main= "Histogram of normal data - Relative Total Change in Height", xlab="Approximately normally distributed data")
lines(density(GAF$RHIN), col=2)
qqnorm(GAF$RHIN,main="QQ plot of normal data- Relative Total Change in Height",pch=19) 
qqline(GAF$RHIN) 
shapiro.test(GAF$RHIN)
leveneTest(RHIN ~ Soil*Inundation, data=GAF)

          #* Transforming data & testing assumptions using Linear Mixed Model:
          GAF$LogRHIN <- log(GAF$RHIN + 1)
          lmmRHIN <- lmer(LogRHIN ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = GAF,
                 REML = FALSE)
          hist(resid(lmmRHIN,type="deviance"))
          qqnorm(resid(lmmRHIN,type="deviance"))
          qqline(resid(lmmRHIN,type="deviance"))
          #* Conclusion: residuals are normal. Proceeded with LM!


#IH = Initial Height of plants
hist(P$IH, probability = T, main= "Histogram of normal data - Leave Dry Biomass", xlab="Approximately normally distributed data")
lines(density(P$IH), col=2)
qqnorm(P$IH,main="QQ plot of normal data- Leave Dry Biomass",pch=19) 
qqline(P$IH) 
shapiro.test(P$IH) # Test of Normality: met assumptions
leveneTest(IH ~ Soil*Inundation, data=P)

#FH = Final Height of plants
hist(P$IH, probability = T, main= "Histogram of normal data - Leave Dry Biomass", xlab="Approximately normally distributed data")
lines(density(P$FH), col=2)
qqnorm(P$FH,main="QQ plot of normal data- Leave Dry Biomass",pch=19) 
qqline(P$FH) 
shapiro.test(P$FH) # Test of Normality : met assumptions
leveneTest(FH ~ Soil*Inundation, data=P) # 

#****************************************************************************
#            PERFORMANCE VARIABLES: DRY BIOMASS MEASUREMENTS
#****************************************************************************

# Aboveground biomass (AB)
hist(P$ADW, probability = T, main= "Histogram of normal data - Aboveground Dry Biomass", xlab="Approximately normally distributed data")
lines(density(P$SDW), col=2)
qqnorm(P$ADW,main="QQ plot of normal data- Aboveground Dry Biomass",pch=19) 
qqline(P$ADW) 
shapiro.test(P$ADW) # Test of Normality
leveneTest(ADW ~ Soil*Inundation, data=P)

          #* Transforming data & testing assumptions using Linear Mixed Model:
          P$LogADW <- log(P$ADW)
          lmmADW <- lmer(LogADW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                REML = FALSE)
          hist(resid(lmmADW,type="deviance"))
          qqnorm(resid(lmmADW,type="deviance"))
          qqline(resid(lmmADW,type="deviance"))
          #* Conclusion: residuals are normal. Proceeded with LM!


# Belowground biomass (BB)
hist(P$BDW, probability = T, main= "Histogram of normal data - Belowground Dry Biomass", xlab="Approximately normally distributed data")
lines(density(P$BDW), col=2)
qqnorm(P$BDW,main="QQ plot of normal data- Belowground Dry Biomass",pch=19) 
qqline(P$BDW) 
shapiro.test(P$BDW) # Test of Normality
leveneTest(BDW ~ Soil*Inundation, data=P)

          #* Transforming data & testing assumptions using Linear Mixed Model:
          P$LogBDW <- log(P$BDW)
          lmmBDW <- lmer(LogBDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
               REML = TRUE)
          hist(resid(lmmBDW,type="deviance"))
          qqnorm(resid(lmmBDW,type="deviance"))
          qqline(resid(lmmBDW,type="deviance"))
          #* Conclusion: residuals are normal. Proceeded with LM!
          
          Anova(lmmBDW)



# Total Dry Biomass (TDB)
hist(P$TDB, probability = T, main= "Histogram of normal data - Shoot to root ratio", xlab="Approximately normally distributed data")
lines(density(P$TDB), col=2)
qqnorm(P$TDB,main="QQ plot of normal data- Shoot to root ratio",pch=19) 
qqline(P$TDB) 
shapiro.test(P$TDB) # Test of Normality
leveneTest(TDB ~ Soil*Inundation, data=P)

        #* Transforming data & testing assumptions using Linear Mixed Model:
        P$LogTDB <- log(P$TDB)
        lmmTDB <- lmer(LogTDB ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                   REML = TRUE)
        hist(resid(lmmTDB,type="deviance"))
        qqnorm(resid(lmmTDB,type="deviance"))
        qqline(resid(lmmTDB,type="deviance"))
        #* Conclusion: residuals are normal. Proceeded with LM!
        plot(lmmTDB)
        
        #*testing significance of fixed factors:
        Anova(lmmTDB)
        
        #*testing significance of random factor:
        lmTDB <- lm(LogTDB ~ Inundation + Soil + Soil*Inundation, data = P)
        anova(lmmTDB,lmTDB)                      
        
        

# DryBiomassRatio (Shoot:Root ratio)
hist(P$DryBiomassRatio, probability = T, main= "Histogram of normal data - Shoot to root ratio", xlab="Approximately normally distributed data")
lines(density(P$BDW), col=2)
qqnorm(P$DryBiomassRatio,main="QQ plot of normal data- Shoot to root ratio",pch=19) 
qqline(P$DryBiomassRatio) 
shapiro.test(P$DryBiomassRatio) # Test of Normality
leveneTest(DryBiomassRatio ~ Soil*Inundation, data=P)
        
        #* Transforming data & testing assumptions using Linear Mixed Model:
        P$LogABratio <- log(P$DryBiomassRatio)
        lmmABratio <- lmer(LogABratio ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = TRUE)
        hist(resid(lmmABratio,type="deviance"))
        qqnorm(resid(lmmABratio,type="deviance"))
        qqline(resid(lmmABratio,type="deviance"))
        #* Conclusion: residuals are normal. Proceeded with LM!
        plot(lmmABratio)
        
        #*testing significance of fixed factor
        Anova(lmmABratio)
        
        #testing significance of random factor
        lmABratio <- lm(LogABratio ~ Inundation + Soil + Soil*Inundation, data = P)
        anova(lmmABratio,lmABratio)
        
        
        
# Leaf biomass
hist(P$LDW, probability = T, main= "Histogram of normal data - Shoot to root ratio", xlab="Approximately normally distributed data")
lines(density(P$LDW), col=2)
qqnorm(P$LDW,main="QQ plot of normal data- Shoot to root ratio",pch=19) 
qqline(P$LDW) 
shapiro.test(P$LDW) # Test of Normality
leveneTest(LDW ~ Soil*Inundation, data=P)
        
        #* Transforming data & testing assumptions using Linear Mixed Model: 
        P$LogLDW <- log(P$LDW)
        lmmLDW <- lmer(LogLDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = TRUE)
        hist(resid(lmmLDW,type="deviance"))
        qqnorm(resid(lmmLDW,type="deviance"))
        qqline(resid(lmmLDW,type="deviance"))
        #* Conclusion: residuals are normal. Proceeded with LM!
        
        
# Stem biomass
hist(P$SDW, probability = T, main= "Histogram of normal data - Shoot to root ratio", xlab="Approximately normally distributed data")
lines(density(P$SDW), col=2)
qqnorm(P$SDW,main="QQ plot of normal data- Shoot to root ratio",pch=19) 
qqline(P$SDW) 
shapiro.test(P$SDW) # Test of Normality
leveneTest(SDW ~ Soil*Inundation, data=P)
        
        #* Transforming data & testing assumptions using Linear Mixed Model: 
        P$LogSDW <- log(P$SDW)
        lmmSDW <- lmer(LogSDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                       REML = TRUE)
        hist(resid(lmmSDW,type="deviance"))
        qqnorm(resid(lmmSDW,type="deviance"))
        qqline(resid(lmmSDW,type="deviance"))
        #* Conclusion: residuals are normal. Proceeded with LM!
        

        
############################################################################################################
# Running Linear Mixed Models & extracting LSmeans
#############################################################################################################
#TDB = Total Dry Biomass
P$LogTDB <- log(P$TDB)

lmmTDB <- lmer(LogTDB ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P, REML = TRUE)   
Anova(lmmTDB) #testing significance of fixed effects
            #*Posthoc tests & LSmeans
            library(multcomp)
            library(multcompView)
            lsmeans(lmmTDB, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            lsmeans(lmmTDB, pairwise ~ Inundation|Soil, adjust ="tukey") #within soil treatments among hydrological regimes
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            TDB.lsm <- lsmeans(lmmTDB, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            cld(TDB.lsm, alpha = .50)
            
            #*Testing random factor effect
            lmmTDB1 <- lmer(LogTDB ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                            REML = FALSE)
            lmmTDB0 <- lm(LogTDB ~ Inundation + Soil + Soil*Inundation, data = P)
            anova(lmmTDB1, lmmTDB0)
            

#Shoot:Root ratio
            
P$LogABratio <- log(P$DryBiomassRatio)
lmmABratio <- lmer(LogABratio ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = FALSE)   
            
            #*Posthoc tests & LSmeans
            library(multcomp)
            library(multcompView)
            lsmeans(lmmABratio, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            lsmeans(lmmABratio, pairwise ~ Inundation|Soil, adjust ="tukey") #within soil treatments among hydrological regimes
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            ABratio.lsm <- lsmeans(lmmABratio, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            cld(ABratio.lsm, alpha = .50)
            
            #*Testing random factor effect
            lmmABratio1 <- lmer(LogABratio ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                            REML = FALSE)
            lmmABratio0 <- lm(LogABratio ~ Inundation + Soil + Soil*Inundation, data = P)
            anova(lmmABratio1, lmmABratio0)
            

#########################################################################################
            ##  ANALYSES OF FUNGAL COLONIZATION
##########################################################################################
            

#******************************************************************************************************************************************************* 
            #*        Analyses in total % fungal colonization for non-sterile treatment only
#*******************************************************************************************************************************************************
            
            #Dataset:re-structuring our original data to only include the total AMF and DSE colonization in the non-sterile treatments
            E <- read.table("Fungal.Colonization.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 

            totals <- E[c(1:16),c("ID", "Soil", "Inundation", "Seed","TAMF","TotalFungi", "TDSE", "TDB")]      
            
            short.fungi.totals <- gather(totals, Fungi, Percentage, -ID, -Soil, -Inundation,-Seed, -TDB) #This reshapes the matrix for Strain
            
            totals2<- E[c(1:16),c("ID", "Soil", "Inundation", "Seed","Oomycetes", "TDB", "RTHO", "FH")]   
            short.oomycetes <- gather(totals2, Fungi, Percentage, -ID, -Soil, -Inundation,-Seed, -TDB, -RTHO, -FH)
            
#******************************************************************************************************************           
# Objective 2: Is the growth benefit provided by mutualistic interactions with soil microbes, such as AMF and DSE? 
#******************************************************************************************************************
# Evaluating relationship between total dry biomass (plant performance) and total fungal colonization:
                        
   lmTDB.FP <- lm(TDB ~ TotalFungi, data= totals)                     
    plot( lmTDB.FP) 
    shapiro.test(resid(lmTDB.FP,type="deviance"))
    anova( lmTDB.FP)
    summary(lmTDB.FP)
    
    #Let's try transforming TDB because of lack of normality of residuals
    lmTDB.FP2 <- lm(log(TDB + 1) ~ TotalFungi, data= totals)                     
    plot( lmTDB.FP2) 
    shapiro.test(resid(lmTDB.FP2,type="deviance"))
    anova(lmTDB.FP2)
    summary(  lmTDB.FP2) 
  
    library(ggplot2)
    install.packages("ggpmisc")
    library(ggpmisc)
    formula <- y ~ x
    FungiTDB.regression.plot <-ggplot(totals, aes(x=TotalFungi, y=TDB)) +
      geom_point(shape=1) +
      scale_colour_hue(l=0) + 
      #geom_smooth (alpha=0.3, size=0) +
      #stat_smooth (geom="line", alpha=0.3, size=3) +
      geom_smooth(method=lm,se=TRUE, alpha=0.3,colour="black") + 
      stat_poly_eq(formula = formula, 
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                   parse = TRUE) +  
      labs(x="Colonization (%)", y="Plant total dry biomass (mg)") +
      theme(panel.background = element_rect(fill = "white", colour = "black"), 
            axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
            axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
            axis.text.x = element_text(size=16, colour ="black"),
            axis.text.y = element_text(size=16, colour ="black"),
            legend.position="top",
            panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
    
    
    library("ggpubr")
    ggscatter(totals, x="TotalFungi", y="TDB", 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "%Colonization", ylab = "Plant dry biomass (mg)")
    
#Evaluating if the colonization percentage of each fungal group influenced seedling growth.    
  #filtering by fungal group:
    amf <- subset(short.fungi.totals, Fungi=="TAMF")
    dse <- subset(short.fungi.totals, Fungi=="TDSE")
    
            lmAMF.TDB <- lm(TDB ~ Percentage, data=amf)
            plot( lmAMF.TDB)
            anova(lmAMF.TDB)
            summary(lmAMF.TDB)
            
            ggscatter(amf, x="Percentage", y="TDB",
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "%Colonization", ylab = "Plant dry biomass (mg)")
            
            
                 
                  
            ggplot(amf, aes(Percentage,TDB)) + geom_smooth(method = "lm") + geom_point() + 
                    facet_grid(~ Inundation, margins=TRUE)
            
            
            lmAMF.DSE <- lm(TDB ~ Percentage, data=dse)
            plot(lmAMF.DSE)
            anova(lmAMF.DSE)
            summary(lmAMF.DSE)
            
            ggscatter(dse, x="Percentage", y="TDB",
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "%Colonization", ylab = "Plant dry biomass (mg)")
            
            amf.dse <- subset(short.fungi.totals, Fungi != "TotalFungi")
            fit <- lm(TDB~Percentage, amf.dse)
            
            Fungi.regression.plot <-ggplot(amf.dse, aes(x=Percentage, y=TDB, color=Fungi, shape=Fungi)) +
              geom_point() +
              #scale_colour_hue(l=50) + 
              scale_color_manual (values = c("black","grey")) +
              scale_shape_manual(values=c(21,22)) +
              geom_smooth(method=lm,se=TRUE, alpha=0.3) + 
              stat_poly_eq(formula = formula, 
                          aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                          parse = TRUE) +  
              labs(x="Colonization (%)", y="Plant total dry biomass (mg)") +
              theme(panel.background = element_rect(fill = "white", colour = "black"), 
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=16, colour ="black"),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
    
                             
                        
#What about the presence of pathogens such as Oomycetes?:      
            lmOomycetes <- lm(TDB ~ Percentage, data=short.oomycetes)
            plot(lmOomycetes)
            anova(lmOomycetes)
            summary(lmOomycetes) 
                
            
            Oomycetes.regression.plot <-ggplot(short.oomycetes, aes(x=Percentage, y=TDB, color=Fungi, shape=Fungi)) +
              geom_point() +
              #scale_colour_hue(l=50) + 
              scale_color_manual (values = c("grey")) +
              scale_shape_manual(values=c(23)) +
              geom_smooth(method=lm,se=TRUE) + 
              stat_poly_eq(formula = formula, 
                           aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                           parse = TRUE) + 
              labs(x="Colonization (%)", y="Plant total dry biomass (mg)") +
              theme(panel.background = element_rect(fill = "white", colour = "black"), 
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=16, colour ="black"),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            oomycetes.tdb <- ggscatter(short.oomycetes, x="Percentage", y="TDB",
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "%Colonization", ylab = "Plant dry biomass (mg)")
            
            lmOomycetes.RTHO <- lm(RTHO ~ Percentage, data=short.oomycetes)
            plot(lmOomycetes)
            anova(lmOomycetes.RTHO)
            summary(lmOomycetes.RTHO) 
            
            
            
            Oomycetes.RTHO.plot <-ggplot(short.oomycetes, aes(x=Percentage, y=RTHO, color=Fungi, shape=Fungi)) +
              geom_point() +
              #scale_colour_hue(l=50) + 
              scale_color_manual (values = c("grey")) +
              scale_shape_manual(values=c(23)) +
              geom_smooth(method=lm,se=TRUE) + 
              stat_poly_eq(formula = formula, 
                           aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                           parse = TRUE) + 
              labs(x="Colonization (%)", y="Relative change in height (mm)") +
              theme(panel.background = element_rect(fill = "white", colour = "black"), 
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=16, colour ="black"),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            
            
            oomycetes.RTHO <- ggscatter(short.oomycetes, x="Percentage", y="RTHO",
                                       add = "reg.line", conf.int = TRUE, 
                                       cor.coef = TRUE, cor.method = "pearson",
                                       xlab = "%Colonization", ylab = "Plant Change in Height (mg)")
            
            
            lmOomycetes.FH <- lm(sqrt(FH) ~ sqrt(Percentage), data=short.oomycetes)
            plot(lmOomycetes)
            anova(lmOomycetes)
            summary(lmOomycetes) 
            
            #It might be best to model the percentage of oomycetes with a poisson distribution:
            str(short.oomycetes)
            short.oomycetes$Percentage <- as.numeric(short.oomycetes$Percentage)
            hist(short.oomycetes$Percentage) #yep nice zero- inclined distribution
            mean(short.oomycetes$Percentage)
            var(short.oomycetes$Percentage)  #variance is twice the mean, we might have over-dispersion
            
            poisson.model.FH <- glm(Percentage ~ FH, family = poisson(link = "log"), data=short.oomycetes)
            summary(poisson.model.FH)
            
            poisson.model.RTHO <- glm(Percentage ~ RTHO, family = poisson(link = "log"), data=short.oomycetes)
            summary(poisson.model.RTHO)
            plot(poisson.model.RTHO)
            
            oomycetes.FH <- ggscatter(short.oomycetes, x="Percentage", y="FH",
                                        add = "reg.line", conf.int = TRUE, 
                                        cor.coef = TRUE, cor.method = "pearson",
                                        xlab = "%Colonization", ylab = "Plant Change in Height (mg)")
            
            
           Oomycetes.reg <- ggarrange(Oomycetes.regression.plot, Oomycetes.RTHO.plot, nrow=1, ncol=2,common.legend = TRUE, legend="top",labels=c("(b)", "(c)"))
            ggsave("oomycetes.reg.pdf", plot=Oomycetes.reg, device="pdf", scale = 1, width = 20, height = 10, units = "cm", dpi = 300)
                     
            FIGURE4 <- ggarrange(FungiTDB.regression.plot, Fungi.regression.plot, Oomycetes.regression.plot, Oomycetes.RTHO.plot, ncol=4, nrow=1, common.legend=FALSE, legend="top",labels=c("(a)", "(b)", "(c)", "(d)"))           
            ggsave("FIGURE4.pdf", plot= FIGURE4, device="pdf", scale = 1, width = 40, height = 10, units = "cm", dpi = 300)
            
                     
#******************************************************************************************************************           
# Objective 3: Is the colonization of beneficial microbes in host seedlings affected by flooding?
#******************************************************************************************************************            

            lmerfungitotal <- lmer(Percentage ~ Inundation + (1|Seed), data= short.fungi.totals)
            lmfungitotal <- lm(Percentage ~ Inundation, data= short.fungi.totals)
            anova(lmerfungitotal,lmfungitotal) # testing significance of random factor
            
            Anova(lmfungitotal, Type="III")
            
            #Testing Normality:
            hist(resid(lmfungitotal,type="deviance"))
            qqnorm(resid(lmfungitotal,type="deviance"))
            qqline(resid(lmfungitotal,type="deviance"))
            shapiro.test(resid(lmfungitotal,type="deviance")) #it is normal, but let's see if transformation helps
            plot(lmfungitotal)
            
            short.fungi.totals$sqrt.percentage <- sqrt(short.fungi.totals$Percentage) #using square-root transformation since it is a percentage
            lmfungitotal.1 <- lm(sqrt.percentage ~ Inundation, data= short.fungi.totals)
            hist(resid(lmfungitotal.1,type="deviance"))
            qqnorm(resid(lmfungitotal.1,type="deviance"))
            qqline(resid(lmfungitotal.1,type="deviance"))
            shapiro.test(resid(lmfungitotal.1,type="deviance")) #it is normal, but let's see if transformation helps
            plot(lmfungitotal.1)
            
            #let's use log to reduce effect of outliers:
            short.fungi.totals$log.percentage <- log(short.fungi.totals$Percentage + 1) #using square-root transformation since it is a percentage
            lmfungitotal.2 <- lm(log.percentage ~ Inundation, data= short.fungi.totals)
            hist(resid(lmfungitotal.1,type="deviance"))
            qqnorm(resid(lmfungitotal.1,type="deviance"))
            qqline(resid(lmfungitotal.1,type="deviance"))
            shapiro.test(resid(lmfungitotal.1,type="deviance")) #it is normal, but let's see if transformation helps
            plot(lmfungitotal.1)
            
            #*Posthoc tests
            library(lsmeans)
            lsmeans(lmfungitotal, pairwise ~ Inundation, adjust ="Tukey") #all pairwise comparisons
                      #*Conclusion no differences in fungal colonization among flooding treatments
            
            
#******************************************************************************************************************           
# Objective 4: Does flooding differentially affect the colonization of AMF and DSE?. 
#******************************************************************************************************************            
           #Let's visualize the data and potential outliers
            short.fungi.totals$Fungi = as.factor(short.fungi.totals$Fungi)
            with(short.fungi.totals, plot(Percentage ~ Inundation))
            with(short.fungi.totals, plot(Percentage ~ Fungi))
            
            #checking outlier
            b =boxplot(Percentage ~ Inundation*Fungi, data=short.fungi.totals)$out  #This tells the specific values that are outlier
            o <- short.fungi.totals[ short.fungi.totals$Percentage %in% b,] #This let me see which values are considered outliers
            return(o)
            
            #removing outlier to test effect further
            short.fungi.totals.noout <- short.fungi.totals[!short.fungi.totals$ID %in% o$ID,]
            
      ###Models to test: does the abundance between fungal groups differs among the different treatments?
            short.fungi.totals$Fungi = as.factor(short.fungi.totals$Fungi)
            lmerfungi <- lmer(Percentage ~ Inundation*Fungi + (1|Seed), data= short.fungi.totals)
            lmfungi <- lm(Percentage ~ Inundation*Fungi, data= short.fungi.totals)
            anova(lmerfungi,lmfungi)
            anova(lmerfungi)
            
            #Testing Normality:
            hist(resid(lmfungi,type="deviance"))
            qqnorm(resid(lmfungi,type="deviance"))
            qqline(resid(lmfungi,type="deviance"))
            shapiro.test(resid(lmfungi,type="deviance"))
            plot(lmfungi)
            
            #Testing same model removing the outlier:
            
            lmerfungi.noout <- lm(Percentage ~ Inundation*Fungi, data= short.fungi.totals.noout)
            anova(lmerfungi.noout)
            Anova(lmerfungi.noout, type="III")
            plot(lmerfungi.noout)
            
            
            #testing with square-root transformation:
            lmfungi.1 <- lm(sqrt(Percentage) ~ Inundation*Fungi, data= short.fungi.totals)
            plot( lmfungi.1)
            hist(resid(lmfungi.1,type="deviance"))
            qqnorm(resid(lmfungi.1,type="deviance"))
            qqline(resid(lmfungi.1,type="deviance"))
            shapiro.test(resid(lmfungi.1,type="deviance"))
            Anova(lmfungi.1, type="III")
            
            lmfungi.1.noout <- lm(sqrt(Percentage) ~ Inundation*Fungi, data= short.fungi.totals.noout)
            plot(lmfungi.1.noout)
            
            shapiro.test(resid(lmfungi.1.noout,type="deviance"))
            Anova(lmfungi.1.noout, type="III")
            
            
            #We are staying with the un-transformed data with the removal of the one outlier because this model fits better our data
            # and provides better assumptions of normality
            
           
            
            #FINAL MODEL:
            lmerfungi.noout <- lm(Percentage ~ Inundation*Fungi, data= short.fungi.totals.noout)
            anova(lmerfungi.noout)
            Anova(lmerfungi.noout, type="III")
            plot(lmerfungi.noout)
            
            #*Posthoc tests
            library(lsmeans)
            lsmeans(lmerfungi.noout, pairwise ~ Inundation*Fungi, adjust ="Tukey") #all pairwise comparisons
            lsmeans(lmerfungi.noout, pairwise ~ Inundation|Fungi, adjust ="Tukey") #within soil treatments among hydrological regimes
            lsmeans(lmerfungi.noout, pairwise ~ Fungi|Inundation, adjust ="Tukey") #within soil treatments among hydrological regimes
            lsmeans(lmerfungi.noout, pairwise ~ Fungi, adjust ="Tukey") #within soil treatments among hydrological regimes
            
            
            library(multcomp)
            library(multcompView)
            lsmean.fungi = lsmeans(lmerfungi.noout,pairwise ~ Fungi*Inundation, adjust="tukey")
            fungi.cld <- cld(lsmean.fungi$lsmeans,alpha=0.05, Letters=letters) 
            fungi.cld$Hydrology <- factor(fungi.cld$Inundation, levels=c("None","Typical","Extreme")) ### Order the levels for printing
            fungi.cld$.group=gsub(" ", "", fungi.cld$.group)  ###  Remove spaces in .group 
            
            #Plotting lsmeans tukey groupings for model excluding outlier
            library(ggplot2)
            fungifitness.plot <- ggplot(data = fungi.cld, mapping = aes(x=Hydrology, y=lsmean, fill=Fungi)) + 
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              scale_fill_manual (values = c("#3a3336ff", "grey")) +
              labs(x="Flooding Conditions", y= "Colonization(%)") +
              theme(panel.background = element_rect(fill = "white", colour = "black"), 
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=16, colour ="black"),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
            geom_text(aes(label = .group),
            position = position_dodge(0.9),
            vjust = -5)
            
            
            #Plotting lsmeans tukey groupings for model including outlier
            lmerfungi <- lm(Percentage ~ Inundation*Fungi, data= short.fungi.totals)
            anova(lmerfungi)
            Anova(lmerfungi, type="III")
            plot(lmerfungi)
            
            #*Posthoc tests
            library(lsmeans)
            lsmeans(lmerfungi, pairwise ~ Inundation*Fungi, adjust ="Tukey") #all pairwise comparisons
            lsmeans(lmerfungi, pairwise ~ Inundation|Fungi, adjust ="Tukey") #within soil treatments among hydrological regimes
            lsmeans(lmerfungi, pairwise ~ Fungi|Inundation, adjust ="Tukey") #within soil treatments among hydrological regimes
            lsmeans(lmerfungi, pairwise ~ Fungi, adjust ="Tukey") #within soil treatments among hydrological regimes
            library(multcomp)
            library(multcompView)
            lsmean.fungi = lsmeans(lmerfungi,pairwise ~ Fungi*Inundation, adjust="tukey")
            fungi.cld <- cld(lsmean.fungi$lsmeans,alpha=0.05, Letters=letters) 
            fungi.cld$Hydrology <- factor(fungi.cld$Inundation, levels=c("None","Typical","Extreme")) ### Order the levels for printing
            fungi.cld$.group=gsub(" ", "", fungi.cld$.group)  ###  Remove spaces in .group 
            
            
            library(ggplot2)
            fungifitness.plot2 <- ggplot(data = fungi.cld, mapping = aes(x=Hydrology, y=lsmean, fill=Fungi)) + 
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              scale_fill_manual (values = c("#3a3336ff", "grey")) +
              labs(x="Flooding Conditions", y= "Colonization(%)") +
              theme(panel.background = element_rect(fill = "white", colour = "black"), 
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=16, colour ="black"),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),
                    legend.position = "top") +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -5)
            
            
        ggsave("fungifitness.pdf", plot= fungifitness.plot2, device="pdf", scale = 1, width = 10, height = 10, units = "cm", dpi = 300)
            

        
#******************************************************************************************************************           
        # Objective 4.1: Does the abundance of fungal structures differed among flooding treatments?. 
#******************************************************************************************************************         
        
        

#Summary of all types of fungal structures

        #Firts we need to wraggle our data matrix a bit:
        E1 <- E[,c(1,2,3,4,10,11,12,13,14)] #this matrix contains data for both sterile and non-sterile treatment
        
        library(tidyr)
        short.fungi <- gather(E1, Fungi, Percentage, -ID, -Soil, -Inundation,-Seed) #This reshapes the matrix for Strain
        #obtaining averages of the percentages for each fungal group observed:
        sum.fungi <-summarySE(data=short.fungi, measurevar= "Percentage", groupvars = c("Soil", "Inundation", "Fungi"), na.rm = TRUE)
        sum.fungi$Inundation <- factor(sum.fungi$Inundation, levels=c("None","Typical","Extreme"))
        # The palette with grey:
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        fungi.diver.plot <- ggplot(data = sum.fungi, aes(x = Inundation, y = Percentage, fill = Fungi)) + 
          geom_bar(colour= "black", stat= "identity", size=0.5) +
          #geom_errorbar(aes(ymin=RGR.2-se, ymax=RGR.2+se), width=0.3, colour="black", position=position_dodge(0.9)) +
          facet_grid(~Soil, switch = "x", scales = "free_x", space = "free_x") +
          scale_fill_manual(values =cbPalette) +
          scale_fill_grey(start=0.2, end = 0.8, name="Fungal structures",
                         breaks=c("ARB_percentage", "HYP_percentage", "VES_percentage", "MS_percentage", "SHYP_percentage"),
                        labels=c("ARB", "HYP", "VES", "MS", "SHYP")) +
          labs(x="Treatments", y= "Colonization(%)") +
          theme(panel.spacing = unit(0.5, "lines"), 
                strip.background = element_blank(),
                strip.placement = "outside",
                strip.text = element_text(face="bold", vjust=1.0, size = 14),
                axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                axis.text.y = element_text(size=16, colour ="black"),
                panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
        
        #Let' perform some accumulation curves to see if we capture the diversity with just five samples per treatment combination
        library(vegan)
        install.packages("labdsv")
        library(labdsv)
        short.fungi$plant_treatment <- paste0(short.fungi$ID, short.fungi$Inundation)
        
        E2 <- subset(short.fungi, Soil == "Non-sterile") #This matrix contains only data for
        #for typical treatment
        E2_typical <- subset(E2, Inundation == "Typical")
        E2_typical <- E2_typical[,c(7,5,6)]
        typical_matrix <- matrify(as.data.frame(E2_typical))
        accum.TYP.0 <- specaccum(typical_matrix)
        plot(accum.TYP.0, ci.type="polygon", ci.col="yellow")
        H_typical <- diversity(typical_matrix, "simpson") #estimating Simpson Index within each individual plant
        
        
        #for none flooding
        E2_none <- subset(E2, Inundation == "None")
        E2_none <- E2_none[,c(7,5,6)]
        none_matrix <- matrify(as.data.frame(E2_none))
        accum.NONE <- specaccum(none_matrix)
        plot(accum.NONE, ci.type="polygon", ci.col="yellow")
        H_typical <- diversity(typical_matrix, "simpson") #estimating Simpson Index within each individual plant
        
        
        #for extreme flooding
        E2_extreme <- subset(E2, Inundation == "Extreme")
        E2_extreme <- E2_extreme[,c(7,5,6)]
        extreme_matrix <- matrify(as.data.frame(E2_extreme))
        accum.extreme <- specaccum(extreme_matrix)
        plot(accum.extreme, ci.type="polygon", ci.col="yellow")
        H_typical <- diversity(typical_matrix, "simpson") #estimating Simpson Index within each individual plant
        
        #wrangle into a single data frame
        fungi.accum <- rbind(data.frame(plants=accum.TYP.0$sites, richness=accum.TYP.0$richness, SD=accum.TYP.0$sd, Flooding="Typical"), 
                             data.frame(plants=accum.NONE$sites, richness=accum.NONE$richness, SD=accum.NONE$sd, Flooding="None"), 
                             data.frame(plants=accum.extreme$sites, richness=accum.extreme$richness, SD=accum.extreme$sd, Flooding="Extreme"))
        
        #plotting accumulation curve for each individual flooding treatment:
        flooding.accum.curve <- ggplot(data=fungi.accum, aes(x=plants, y=richness, group=Flooding)) +
          geom_errorbar(aes(x=plants, ymin=richness-SD, ymax=richness+SD), color="black", alpha=0.5, width=0.1)+
          geom_point(aes(color=Flooding)) +
          geom_line(aes(color=Flooding))+
          #scale_colour_manual(values=c("blue","red","green")) +
          scale_colour_manual(values = c("#999999ff","grey","#3a3336ff")) +
          #scale_color_brewer(palette="Dark2") +
          xlab("Plants sampled (no.)")+ylab("Fungal structures (No)")+
          scale_x_continuous(limits=c(0,6))+
          theme(panel.background = element_rect(fill = "white", colour = "black"), 
                axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                axis.text.x = element_text(size=14, colour ="black"),
                axis.text.y = element_text(size=14, colour ="black"))
        
        Fig3Sb <- ggarrange(flooding.accum.curve, fungi.diver.plot, ncol=2,nrow=1, common.legend = FALSE, legend="top", labels=c("(a)", "(b)"))
        ggsave("Fig3Sb.pdf", plot=Fig3Sb , device="pdf", scale = 1, width = 30, height = 10, units = "cm", dpi = 300)
        
        
        
        
#Analyzing differences in the fungal structures
        E2 <- subset(short.fungi, Soil == "Non-sterile") #This matrix contains only data for
        E2 <- subset(E2, ID != "32")  #removing outlier from Typical 
        E2.ARB <- subset(E2, Fungi == "ARB")
        E2.VES <- subset(E2, Fungi == "VES")
        E2.HYP <- subset(E2, Fungi == "HYP")
        E2.SHYP <- subset(E2, Fungi == "SHYP")
        E2.MS <- subset(E2, Fungi == "MS") 
           
            lmARB <- lm(sqrt(Percentage) ~ Inundation, data=E2.ARB)
            summary(lmARB)
            plot(lmARB)
            bartlett.test(Percentage ~ Inundation, data=E2.ARB)
            shapiro.test(resid(lmARB,type="deviance"))
            Anova(lmARB, type="III") #Significant
            
            hist(resid(lmARB,type="deviance"))
            qqnorm(resid(lmARB,type="deviance"))
            qqline(resid(lmARB,type="deviance"))
            
            hist(E2.ARB$Percentage)
            qqnorm(E2.ARB$Percentage)
            qqline(E2.ARB$Percentage)
            shapiro.test(E2.ARB$Percentage)
            
            lsmean.ARB = lsmeans(lmARB,pairwise ~ Inundation, adjust="tukey")
            ARB.cld <- cld(lsmean.ARB$lsmeans,alpha=0.05, Letters=letters) 
            
            
            
            
            lmVES <- lm(sqrt(Percentage) ~ Inundation, data=E2.VES)
            summary(lmVES)
            plot(lmVES)
            Anova(lmVES, type="III") #No-significant
            hist(resid(lmVES,type="deviance"))
            qqnorm(resid(lmVES,type="deviance"))
            qqline(resid(lmVES,type="deviance"))
            shapiro.test(resid(lmVES,type="deviance"))
          
            
            lsmean.VES = lsmeans(lmVES,pairwise ~ Inundation, adjust="tukey")
            VES.cld <- cld(lsmean.VES$lsmeans,alpha=0.05, Letters=letters) 
            
            
            lmHYP <- lm(sqrt(Percentage) ~ Inundation, data=E2.HYP)
            summary(lmHYP)
            plot(lmHYP)
            Anova(lmHYP, type="III")  #significant
            hist(resid(lmHYP,type="deviance"))
            qqnorm(resid(lmHYP,type="deviance"))
            qqline(resid(lmHYP,type="deviance"))
            shapiro.test(resid(lmHYP,type="deviance"))
            
            
            lsmean.HYP = lsmeans(lmHYP,pairwise ~ Inundation, adjust="tukey")
            HYP.cld <- cld(lsmean.HYP$lsmeans,alpha=0.05, Letters=letters) 
            
            
            
            lmSHYP <- lm(sqrt(Percentage) ~ Inundation, data=E2.SHYP)
            summary(lmSHYP)
            anova(lmSHYP)
            Anova(lmSHYP, type="III")  #No-significant
            hist(resid(lmSHYP,type="deviance"))
            qqnorm(resid(lmSHYP,type="deviance"))
            qqline(resid(lmSHYP,type="deviance"))
            shapiro.test(resid(lmSHYP,type="deviance"))
            
            
            lmMS <- lm(sqrt(Percentage) ~ Inundation, data=E2.MS)
            summary(lmMS)
            Anova(lmMS, type="III")  #No-significant
            hist(resid(lmMS,type="deviance"))
            qqnorm(resid(lmMS,type="deviance"))
            qqline(resid(lmMS,type="deviance"))
            shapiro.test(resid(lmMS,type="deviance"))
            
            
            #Let's try a Poisson distribution for some of the structures with zeros
            glmMS <- glm(Percentage ~ Inundation,family=poisson(log), data=E2.MS)
            summary(glmMS)
            plot(glmMS)
            residuals(glmMS, type = "pearson")
            plot(glmMS$y ~ fitted(glmMS)) #fitted vs predicted values
            abline(a=0,b=1)
            residualPlots(glmMS)
            library(effects)
            plot(allEffects(glmMS))
            Anova(glmMS, type="III") #Significant
            
            
            glmVES <- glm(Percentage ~ Inundation,family=poisson(log), data=E2.VES)
            summary(glmVES)
            plot(glmVES)
            residuals( glmVES, type = "pearson")
            plot(glmVES$y ~ fitted( glmVES)) #fitted vs predicted values
            abline(a=0,b=1)
            residualPlots(glmVES)
            library(effects)
            plot(allEffects(glmVES))
            Anova(glmVES, type="III") #Significant
            
            glmSHYP <- glm(Percentage ~ Inundation,family=poisson(log), data=E2.SHYP)
            summary(glmSHYP)
            library(effects)
            plot(allEffects(glmSHYP))
            Anova(glmSHYP, type="III") #Significant
            
            
            
            
            sumfungistructure = summarySE(data=E2, measurevar= "Percentage", groupvars = c("Fungi", "Inundation"))
            sumfungistructure$Inundation <- factor(sumfungistructure$Inundation, levels=c("None", "Typical", "Extreme")) ### Order the levels for printing
            sumfungistructure$Fungi <- factor(sumfungistructure$Fungi,levels=c("ARB", "HYP", "VES", "SHYP", "MS"))
            sumfungistructure$Fungalgroup <- c("AMF", "AMF", "AMF", "AMF", "AMF", "AMF", "DSE", "DSE", "DSE", "DSE", "DSE", "DSE", "AMF", "AMF", "AMF")
            
            fungi.diver.plot.2 <- ggplot(data =  sumfungistructure, aes(x = Fungi, y = Percentage, fill = Inundation)) + 
              geom_bar(stat="identity", position=position_dodge(), colour="black") +
              geom_errorbar(aes(ymin=Percentage-se, ymax=Percentage+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~Fungalgroup, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey84","grey72","grey30")) +
              #scale_fill_manual(values =cbPalette) +
              labs(x="Fungal structures", y= "Colonization(%)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    #panel.background = element_rect(fill = "white", colour = "black"),
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
              
            
            
 
 
 #******************************************************************************************************************
 # Did all fungi, including Oomycetes colonization differed among flooding treatments?
 #*********************************************************************************************************************
 #Dataset:re-structuring our original data to only include the total AMF and DSE colonization in the non-sterile treatments
 E <- read.table("Fungal.Colonization.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 
 totals.all <- E[c(1:16),c("ID", "Soil", "Inundation", "Seed","TAMF", "Oomycetes", "TDSE", "TDB", "RHIN", "RTHO")]      
 short.fungi.totals.all <- gather(totals.all, Fungi, Percentage, -ID, -Soil, -Inundation,-Seed, -TDB, -RHIN, -RTHO) #This reshapes the matrix for Strain
 
 lmerfungi.all <- lm(Percentage ~ Inundation*Fungi, data= short.fungi.totals.all)
 plot(lmerfungi.all)
 anova(lmerfungi.all)
 Anova(lmerfungi.all, type="III")

 #let's remove the one outlier:
 b1 =boxplot(Percentage ~ Inundation*Fungi, data=short.fungi.totals.all)$out  #This tells the specific values that are outlier
 o1 <- short.fungi.totals.all[ short.fungi.totals.all$Percentage %in% b1,] #This let me see which values are considered outliers
 return(o1)
 
 #removing outlier to test effect further
 short.fungi.totals.all.noout <- short.fungi.totals.all[!short.fungi.totals.all$ID %in% o1$ID,]
 #short.fungi.totals.all.noout <- subset(short.fungi.totals.all, ID != 5)
 
 lmerfungi.all.noout <- lm(log(Percentage + 1) ~ Inundation*Fungi, data=  short.fungi.totals.all.noout)
 plot(lmerfungi.all.noout)
 anova(lmerfungi.all.noout)
 Anova(lmerfungi.all.noout, type="III")
 
 #*Posthoc tests
 library(lsmeans)
 lsmeans(lmerfungi.all.noout, pairwise ~ Inundation*Fungi, adjust ="Tukey") #all pairwise comparisons
 lsmeans(lmerfungi.all.noout, pairwise ~ Inundation|Fungi, adjust ="Tukey") #within soil treatments among hydrological regimes
 lsmeans(lmerfungi.all.noout, pairwise ~ Fungi|Inundation, adjust ="Tukey") #within soil treatments among hydrological regimes
 lsmeans(lmerfungi.all.noout, pairwise ~ Fungi, adjust ="Tukey") #within soil treatments among hydrological regimes
 
 
 library(multcomp)
 library(multcompView)
 lsmean.fungi = lsmeans(lmerfungi.all.noout,pairwise ~ Fungi*Inundation, adjust="tukey", type="response")
 fungi.cld <- cld(lsmean.fungi$lsmeans,alpha=0.05, Letters=letters) 
 fungi.cld$Hydrology <- factor(fungi.cld$Inundation, levels=c("None","Typical","Extreme")) ### Order the levels for printing
 fungi.cld$Microbe <- factor(fungi.cld$Fungi, levels=c("TAMF","TDSE","Oomycetes")) ### Order the levels for printing
 fungi.cld$.group=gsub(" ", "", fungi.cld$.group)  ###  Remove spaces in .group 
 
 #Plotting lsmeans tukey groupings for model excluding outlier
 library(ggplot2)
 fungifitness.plot.all <- ggplot(data = fungi.cld, mapping = aes(x=Hydrology, y=response, fill=Microbe)) + 
   geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
   geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
   scale_fill_manual (values = c("grey18", "grey30", "grey72")) +
   labs(x="Flooding Conditions", y= "Colonization(%)") +
   theme(panel.background = element_rect(fill = "white", colour = "black"), 
         axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
         axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
         axis.text.x = element_text(size=16, colour ="black"),
         axis.text.y = element_text(size=16, colour ="black"),
         panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
         panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
   geom_text(aes(label = .group),
             position = position_dodge(0.9),
             vjust = -5)
 
 FIGURE5 <- ggarrange(fungifitness.plot.all, fungi.diver.plot.2, ncol=2, nrow=1, common.legend=FALSE, legend="top",labels=c("(a)", "(b)"))           
 ggsave("FIGURE5.pdf", plot= FIGURE5, device="pdf", scale = 1, width = 25, height = 15, units = "cm", dpi = 300)
 
 
 #Finding out differences among structures:
 #Firts we need to wraggle our data matrix a bit:
 E1 <- E[,c(1,2,3,4,10,11,12,13,14,15)] #this matrix contains data for both sterile and non-sterile treatment
 
 library(tidyr)
 short.fungi <- gather(E1, Fungi, Percentage, -ID, -Soil, -Inundation,-Seed) #This reshapes the matrix for Strain
 E2 <- subset(short.fungi, Soil == "Non-sterile") #This matrix contains only data for
 E2 <- subset(E2, ID != "32" & ID != "5")  #removing outlier from Typical 
 E2.ARB <- subset(E2, Fungi == "ARB")
 E2.VES <- subset(E2, Fungi == "VES")
 E2.HYP <- subset(E2, Fungi == "HYP")
 E2.SHYP <- subset(E2, Fungi == "SHYP")
 E2.MS <- subset(E2, Fungi == "MS") 
 
 lmARB <- lm(sqrt(Percentage) ~ Inundation, data=E2.ARB)
 summary(lmARB)
 plot(lmARB)
 bartlett.test(Percentage ~ Inundation, data=E2.ARB)
 shapiro.test(resid(lmARB,type="deviance"))
 Anova(lmARB, type="III") #Significant
 
 hist(resid(lmARB,type="deviance"))
 qqnorm(resid(lmARB,type="deviance"))
 qqline(resid(lmARB,type="deviance"))
 
 hist(E2.ARB$Percentage)
 qqnorm(E2.ARB$Percentage)
 qqline(E2.ARB$Percentage)
 shapiro.test(E2.ARB$Percentage)
 
 lsmean.ARB = lsmeans(lmARB,pairwise ~ Inundation, adjust="tukey")
 ARB.cld <- cld(lsmean.ARB$lsmeans,alpha=0.05, Letters=letters) 
 
 
 
 
 lmVES <- lm(sqrt(Percentage) ~ Inundation, data=E2.VES)
 summary(lmVES)
 plot(lmVES)
 Anova(lmVES, type="III") #No-significant
 hist(resid(lmVES,type="deviance"))
 qqnorm(resid(lmVES,type="deviance"))
 qqline(resid(lmVES,type="deviance"))
 shapiro.test(resid(lmVES,type="deviance"))
 
 
 lsmean.VES = lsmeans(lmVES,pairwise ~ Inundation, adjust="tukey")
 VES.cld <- cld(lsmean.VES$lsmeans,alpha=0.05, Letters=letters) 
 
 
 lmHYP <- lm(sqrt(Percentage) ~ Inundation, data=E2.HYP)
 summary(lmHYP)
 plot(lmHYP)
 Anova(lmHYP, type="III")  #significant
 hist(resid(lmHYP,type="deviance"))
 qqnorm(resid(lmHYP,type="deviance"))
 qqline(resid(lmHYP,type="deviance"))
 shapiro.test(resid(lmHYP,type="deviance"))
 
 
 lsmean.HYP = lsmeans(lmHYP,pairwise ~ Inundation, adjust="tukey")
 HYP.cld <- cld(lsmean.HYP$lsmeans,alpha=0.05, Letters=letters) 
 
 
 lmSHYP <- lm(sqrt(Percentage) ~ Inundation, data=E2.SHYP)
 summary(lmSHYP)
 anova(lmSHYP)
 Anova(lmSHYP, type="III")  #No-significant
 hist(resid(lmSHYP,type="deviance"))
 qqnorm(resid(lmSHYP,type="deviance"))
 qqline(resid(lmSHYP,type="deviance"))
 shapiro.test(resid(lmSHYP,type="deviance"))
 
 
 lmMS <- lm(sqrt(Percentage) ~ Inundation, data=E2.MS)
 summary(lmMS)
 Anova(lmMS, type="III")  #No-significant
 hist(resid(lmMS,type="deviance"))
 qqnorm(resid(lmMS,type="deviance"))
 qqline(resid(lmMS,type="deviance"))
 shapiro.test(resid(lmMS,type="deviance"))
 
 
 #Let's try a Poisson distribution for some of the structures with zeros
 glmMS <- glm(Percentage ~ Inundation,family=poisson(log), data=E2.MS)
 summary(glmMS)
 plot(glmMS)
 residuals(glmMS, type = "pearson")
 plot(glmMS$y ~ fitted(glmMS)) #fitted vs predicted values
 abline(a=0,b=1)
 residualPlots(glmMS)
 library(effects)
 plot(allEffects(glmMS))
 Anova(glmMS, type="III") #Significant
 
 
 glmVES <- glm(Percentage ~ Inundation,family=poisson(log), data=E2.VES)
 summary(glmVES)
 plot(glmVES)
 residuals( glmVES, type = "pearson")
 plot(glmVES$y ~ fitted( glmVES)) #fitted vs predicted values
 abline(a=0,b=1)
 residualPlots(glmVES)
 library(effects)
 plot(allEffects(glmVES))
 Anova(glmVES, type="III") #Significant
 
 glmSHYP <- glm(Percentage ~ Inundation,family=poisson(log), data=E2.SHYP)
 summary(glmSHYP)
 library(effects)
 plot(allEffects(glmSHYP))
 Anova(glmSHYP, type="III") #Significant
 
 #plotting structures:
 sumfungistructure = summarySE(data=lmerfungi.all.noout, measurevar= "Percentage", groupvars = c("Fungi", "Inundation"))
 sumfungistructure$Inundation <- factor(sumfungistructure$Inundation, levels=c("None", "Typical", "Extreme")) ### Order the levels for printing
 sumfungistructure$Fungi <- factor(sumfungistructure$Fungi,levels=c("ARB", "HYP", "VES", "SHYP", "MS"))
 sumfungistructure$Fungalgroup <- c("AMF", "AMF", "AMF", "AMF", "AMF", "AMF", "DSE", "DSE", "DSE", "DSE", "DSE", "DSE", "AMF", "AMF", "AMF")
 
 fungi.diver.plot.2 <- ggplot(data =  sumfungistructure, aes(x = Fungi, y = Percentage, fill = Inundation)) + 
   geom_bar(stat="identity", position=position_dodge(), colour="black") +
   geom_errorbar(aes(ymin=Percentage-se, ymax=Percentage+se), width=0.3, colour="black", position=position_dodge(0.9)) +
   facet_grid(~Fungalgroup, switch = "x", scales = "free_x", space = "free_x") +
   scale_fill_manual(values = c("#999999ff","grey","#3a3336ff")) +
   #scale_fill_manual(values =cbPalette) +
   labs(x="Fungal structures", y= "Colonization(%)") +
   theme(panel.spacing = unit(0.5, "lines"), 
         #panel.background = element_rect(fill = "white", colour = "black"),
         strip.background = element_blank(),
         strip.placement = "outside",
         strip.text = element_text(face="bold", vjust=1.0, size = 14),
         axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
         axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
         axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=1),
         axis.text.y = element_text(size=16, colour ="black"),
         panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
         panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
 
 
 
 #Correlation between Oomycetes colonization and plant biomass:
 oomycetes <- subset(short.fungi.totals.all.noout, Fungi == "Oomycetes")
 
 lmerOomycetes <- lm(Percentage ~ Inundation, data=  oomycetes)
 anova( lmerOomycetes)
 
 lmOomycetes <- lm(TDB ~ Percentage, data= oomycetes)
 plot(lmOomycetes)
 anova(lmOomycetes)
 summary(lmOomycetes) 
 plot(TDB ~ Percentage, data= oomycetes)

  #HOW ABOUT WITH PLANT HEIGHT -since oomycetes can stunt plant growth:
 lmOomycetes.RTHO <- lm(RTHO ~ Percentage, data= oomycetes)
 plot(lmOomycetes.RTHO)
 anova(lmOomycetes.RTHO)
 summary(lmOomycetes.RTHO) 
 plot(RTHO ~ Percentage, data= oomycetes)
 
 lmOomycetes.RHIN <- lm(RHIN ~ Percentage, data= oomycetes)
 plot(lmOomycetes.RHIN)
 anova(lmOomycetes.RHIN)
 summary(lmOomycetes.RHIN) 
 plot(RHIN ~ Percentage, data= oomycetes)
 
 
 library(ggplot2)
 Oomycetes.regression.plot <-ggplot(oomycetes, aes(x=Percentage, y=TDB, color=Fungi, shape=Fungi)) +
   geom_point() +
   #scale_colour_hue(l=50) + 
   scale_color_manual (values = c("grey")) +
   scale_shape_manual(values=c(23)) +
   geom_smooth(method=lm,se=FALSE) + 
   labs(x="Colonization (%)", y="Plant total dry biomass (mg)") +
   theme(panel.background = element_rect(fill = "white", colour = "black"), 
         axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
         axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
         axis.text.x = element_text(size=16, colour ="black"),
         axis.text.y = element_text(size=16, colour ="black"),
         panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
         panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
 
 
 #Plotting correlations with plant biomass:
 Fungi.regression.plot.all <-ggplot(short.fungi.totals.all.noout, aes(x=Percentage, y=TDB, color=Fungi, shape=Fungi)) +
   geom_point() +
   #scale_colour_hue(l=50) + 
   scale_color_manual (values = c("grey18", "grey30", "grey72")) +
   scale_shape_manual(values=c(21,22,23)) +
   geom_smooth(method=lm,se=FALSE) + 
   labs(x="Colonization (%)", y="Plant total dry biomass (mg)") +
   theme(panel.background = element_rect(fill = "white", colour = "black"), 
         axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
         axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
         axis.text.x = element_text(size=16, colour ="black"),
         axis.text.y = element_text(size=16, colour ="black"),
         panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
         panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
 
 
 
 
 
 #Checking correlations for all fungal groups X treatments with plant biomass:
 facet.bytreatment <- ggplot(short.fungi.totals.all.noout, aes(Percentage,TDB)) + geom_smooth(method = "lm") + geom_point() + 
   facet_grid(Fungi ~ Inundation, margins=TRUE)
 
ggsave("SI.4.pdf", plot=facet.bytreatment, device="pdf", scale = 1, width = 35, height = 15, units = "cm", dpi = 300)
 
 
##############################################################################################################
#NUTRIENT PLOTS & ANALYSES 
##############################################################################################################
            
            Nu <- read.table("SoilNutrient.edited_12.3.19.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY")           
            
            
            #* This is the data set I end up using:
            #Nu2 <- read.table("Nutrient.Analyses.12.14.18.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY")           
            #Nu2$InocTime <-factor(Nu2$Time,levels=levels(Nu2$Time)[c(2,1)])                        
            
            #Obtaining averages for each treatment
            library(Rmisc)
            library(ggplot2)
            Averages = summarySE(data=Nu, measurevar = "Value", groupvars = c("Nutrient","Inundation", "Soil", "CollectionTime"))
            
            #a general linear model testing effects of inundation and soil, removing the soil before flooding:
            after.flooding.nutrients <- subset(Nu, CollectionTime == "After flooding")
            after.flooding.nutrients$CollectionTime <- factor(after.flooding.nutrients$CollectionTime)
            
            lmNuAF <- lm(Value ~ Nutrient*Inundation*Soil, data=after.flooding.nutrients)
            Anova(lmNuAF, type ="III")
            
            lsmeansNu <- lsmeans(lmNuAF, pairwise ~ Nutrient*Inundation*Soil)
            cld.lsmeansNu <- cld(lsmeansNu$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            
            library(RColorBrewer)
            #Plotting all nutrients at once:
            cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
            Nutrients.plot <- ggplot(data = Averages, aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~Nutrient, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = cbPalette) +
              labs(x="Flooding Treatment", y= "Value") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            
            
            #Let's subset some of the nutrients with higher scales to be able to see the other nutrients
            Averages.Calcium <- subset(Averages, Nutrient == "Calcium")
            Averages.Calcium$Inundation <- factor(Averages.Calcium$Inundation, levels= c("Before", "None", "Typical", "Extreme"))
            Averages.Calcium$CollectionTime <- factor(Averages.Calcium$CollectionTime, levels= c("Before flooding", "After flooding"))
            
            Averages.Carbon <- subset(Averages, Nutrient == "Carbon")
            Averages.Carbon$Inundation <- factor(Averages.Carbon$Inundation, levels= c("Before", "None", "Typical", "Extreme"))
            Averages.Carbon$CollectionTime <- factor(Averages.Carbon$CollectionTime, levels= c("Before flooding", "After flooding"))
            
            
            Averages.Copper <- subset(Averages, Nutrient == "Copper")
            Averages.Copper$Inundation <- factor(Averages.Copper$Inundation, levels= c("Before", "None", "Typical", "Extreme"))
            Averages.Copper$CollectionTime <- factor(Averages.Copper$CollectionTime, levels= c("Before flooding", "After flooding"))
            
            
            Averages.Magnesium <- subset(Averages, Nutrient == "Magnesium")
            Averages.Magnesium$Inundation <- factor(Averages.Magnesium$Inundation, levels= c("Before", "None", "Typical", "Extreme"))
            Averages.Magnesium$CollectionTime <- factor(Averages.Magnesium$CollectionTime, levels= c("Before flooding", "After flooding"))
            
            
            Averages.Nitrogen <- subset(Averages, Nutrient == "Nitrogen")
            Averages.Nitrogen$Inundation <- factor(Averages.Nitrogen$Inundation, levels= c("Before", "None", "Typical", "Extreme"))
            Averages.Nitrogen$CollectionTime <- factor(Averages.Nitrogen$CollectionTime, levels= c("Before flooding", "After flooding"))
            
            Averages.pH <- subset(Averages, Nutrient == "pH")
            Averages.pH$Inundation <- factor(Averages.pH$Inundation, levels= c("Before", "None", "Typical", "Extreme"))
            Averages.pH$CollectionTime <- factor(Averages.pH$CollectionTime, levels= c("Before flooding", "After flooding"))
            
            
            Averages.Phosphorus <- subset(Averages, Nutrient == "Phosphorus")
            Averages.Phosphorus$Inundation <- factor(Averages.Phosphorus$Inundation, levels= c("Before", "None", "Typical", "Extreme"))
            Averages.Phosphorus$CollectionTime <- factor(Averages.Phosphorus$CollectionTime, levels= c("Before flooding", "After flooding"))
            
            Averages.Potassium <- subset(Averages, Nutrient == "Potassium")
            Averages.Potassium$Inundation <- factor(Averages.Potassium$Inundation, levels= c("Before", "None", "Typical", "Extreme"))
            Averages.Potassium$CollectionTime <- factor(Averages.Potassium$CollectionTime, levels= c("Before flooding", "After flooding"))
            
            Averages.Sodium <- subset(Averages, Nutrient == "Sodium")
            Averages.Sodium$Inundation <- factor(Averages.Sodium$Inundation, levels= c("Before", "None", "Typical", "Extreme"))
            Averages.Sodium$CollectionTime <- factor(Averages.Sodium$CollectionTime, levels= c("Before flooding", "After flooding"))
            
            Averages.Sulfur <- subset(Averages, Nutrient == "Sulfur")
            Averages.Sulfur$Inundation <- factor(Averages.Sulfur$Inundation, levels= c("Before", "None", "Typical", "Extreme"))
            Averages.Sulfur$CollectionTime <- factor(Averages.Sulfur$CollectionTime, levels= c("Before flooding", "After flooding"))
            
            Averages.Zinc <- subset(Averages, Nutrient == "Zinc")
            Averages.Zinc$Inundation <- factor(Averages.Zinc$Inundation, levels= c("Before", "None", "Typical", "Extreme"))
            Averages.Zinc$CollectionTime <- factor(Averages.Zinc$CollectionTime, levels= c("Before flooding", "After flooding"))
            
            #Let's make different plots:
            Calcium.plot <- ggplot(data = Averages.Copper, aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "grey60")) +
              labs(x="Treatment", y= "Calcium (ppm)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            
            Carbon.plot <- ggplot(data = Averages.Carbon, aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "grey60")) +
              labs(x="Treatment", y= "Carbon (%)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Copper.plot <- ggplot(data = Averages.Copper, aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "grey60")) +
              labs(x="Treatment", y= "Copper (ppm)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Magnesium.plot <- ggplot(data = Averages.Magnesium, aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "grey60")) +
              labs(x="Treatment", y= "Magnesium (ppm)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Nitrogen.plot <- ggplot(data = Averages.Nitrogen, aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "grey60")) +
              labs(x="Treatment", y= "Nitrogen (%)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            pH.plot <- ggplot(data = Averages.pH, aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "grey60")) +
              labs(x="Treatment", y= "pH (1:1 Water)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Phosphorus.plot <- ggplot(data = Averages.Phosphorus , aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "grey60")) +
              labs(x="Treatment", y= "Phosphorus (ppm)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            
            Potassium.plot <- ggplot(data = Averages.Potassium , aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "grey60")) +
              labs(x="Treatment", y= "Potassium (ppm)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            
            Sulfur.plot <- ggplot(data = Averages.Sulfur, aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "grey60")) +
              labs(x="Treatment", y= "Sulfur (ppm)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            
            Zinc.plot <- ggplot(data = Averages.Zinc, aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "grey60")) +
              labs(x="Treatment", y= "Zinc (ppm)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            
            Sodium.plot <- ggplot(data = Averages.Sodium, aes(x = Inundation, y = Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "grey60")) +
              labs(x="Treatment", y= "Sodium (ppm)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=45,hjust=1),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            
            # **Plotting all graphs together with a shared legend  
            Nutrients.Macro <- ggarrange(Carbon.plot,Nitrogen.plot,Phosphorus.plot,Potassium.plot,Magnesium.plot,Calcium.plot,Sulfur.plot,pH.plot, ncol=2,nrow=4, common.legend = TRUE, legend="top")
            ggsave("Nutrients.Macro.pdf", plot=Nutrients.Macro, device="pdf", scale = 1, width = 40, height = 40, units = "cm", dpi = 300)
            
            Nutrients.Micro <- ggarrange(Sodium.plot, Zinc.plot,Copper.plot, ncol=1,nrow=3, common.legend = TRUE, legend="top")
            ggsave("Nutrients.Micro.pdf", plot=Nutrients.Micro , device="pdf", scale = 1, width = 20, height = 25, units = "cm", dpi = 300)
            
            
            #Let's run the same general model but separately for each nutrient:
            
            #Model for carbon:
            Carbon <- subset(Nu, Nutrient == "Carbon" & Inundation != "Before")
            Carbon$Nutrient <- factor(Carbon$Nutrient)
            
            lmCarbon <- lm(Value ~ Inundation*Soil, data=Carbon)
            anova(lmCarbon)
            plot(lmCarbon)
            qqPlot(Carbon$Value)
            hist(resid(lmCarbon,type="deviance"))
            qqnorm(resid(lmCarbon,type="deviance"))
            qqline(resid(lmCarbon,type="deviance"))   # Distribution of residuals is normal, proceed with analyses
            shapiro.test(resid(lmCarbon))
            
            lsmeansCarbon <- lsmeans(lmCarbon, pairwise ~ Inundation*Soil)
            cld.Carbon <- cld(lsmeansCarbon$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Carbon$Inundation <- factor(cld.Carbon$Inundation, levels=c("None", "Typical", "Extreme")) ### Order the levels for printing
            cld.Carbon$.group=gsub(" ", "", cld.Carbon$.group)  ###  Remove spaces in .group 
            
            Averages.Carbon.after <-  Averages.Carbon[c(3:8),]  #subseting only for those with averages after flooding
            Averages.Carbon.before <- Averages.Carbon[c(1:2),]  #subsetting only for those before flooding because they were not tested
            #Averages.Carbon.after$Inundation <- factor(Averages.Carbon.after$Inundation, levels= c("None", "Typical", "Extreme"))
            newdata <-cld.Carbon [order(cld.Carbon$Inundation,cld.Carbon$Soil),] #Sorting so we can add the tukey-groupings to the raw means to be plotted
            Averages.Carbon.after.sorted <- Averages.Carbon.after[order(Averages.Carbon.after$Inundation, Averages.Carbon.after$Soil),]
            Averages.Carbon.after.sorted$.group <- newdata$.group
            Averages.Carbon.before$.group <- c("NA", "NA")
            
            Averages.Carbon.after.sorted.letters <- rbind(Averages.Carbon.after.sorted,Averages.Carbon.before)
            
            Carbon.plot2 <- ggplot(data =Averages.Carbon.after.sorted.letters, aes(x = Inundation, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Flooding Treatment", y= "Carbon (%)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -4,hjust=0.5,
                        color="black",
                        size=4)
            
            #testing main effects
            lmCarbon.simple <- lm(Value ~ Inundation + Soil, data=Carbon)
            anova(lmCarbon.simple)
            lsmeansCarbon.soil <- lsmeans(lmCarbon.simple, pairwise ~ Soil)
            cld.Carbon.soil <- cld(lsmeansCarbon.soil$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Carbon.soil$Nutrient <- c("Carbon", "Carbon")
            
            
            lsmeansCarbon.flooding <- lsmeans(lmCarbon.simple, pairwise ~ Inundation)
            cld.Carbon.flooding <- cld(lsmeansCarbon.flooding$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Carbon.flooding$Nutrient <- c("Carbon", "Carbon", "Carbon")
            cld.Carbon.flooding$Inundation <- factor(cld.Carbon.flooding$Inundation, levels = c("None", "Typical", "Extreme"))
            
            
            Carbon.plot.soil <- ggplot(data =cld.Carbon.soil, aes(x = Soil, y=lsmean, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Carbon (%)") +
              ylim(0,11.0) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Carbon.plot.flooding <- ggplot(data =cld.Carbon.flooding, aes(x = Inundation, y=lsmean, fill = Inundation)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey68", "grey92", "gray19")) +
              labs(x="Flooding Treatment", y= "Carbon (%)") +
              ylim(0,11.0) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Averages.Carbon.before <- Averages.Carbon.before[,c(1:7)]
            
            Carbon.plot.before <- ggplot(data = Averages.Carbon.before, aes(x = Soil, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              #geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Carbon (%)") +
              ylim(0,11.0) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Carbon.all <- ggarrange(Carbon.plot.before, Carbon.plot.soil,Carbon.plot.flooding, ncol=3,nrow=1, common.legend = TRUE, legend="top")
            
            
            
            
            
            
            #Model for Calcium
            Calcium <- subset(Nu, Nutrient == "Calcium" & Inundation != "Before")
            Calcium$Nutrient <- factor(Calcium$Nutrient)
            
            lmCalcium <- lm(Value ~ Inundation*Soil, data=Calcium)
            anova(lmCalcium)
            plot(lmCalcium)
            qqPlot(Calcium$Value)
            hist(resid(lmCalcium,type="deviance"))
            qqnorm(resid(lmCalcium,type="deviance"))
            qqline(resid(lmCalcium,type="deviance"))   # Distribution of residuals is normal, proceed with analyses
            shapiro.test(resid(lmCalcium))
            
            lsmeansCalcium <- lsmeans(lmCalcium, pairwise ~ Inundation*Soil)
            cld.Calcium <- cld(lsmeansCalcium$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Calcium$Inundation <- factor(cld.Calcium$Inundation, levels=c("None", "Typical", "Extreme")) ### Order the levels for printing
            cld.Calcium$.group=gsub(" ", "", cld.Calcium$.group)  ###  Remove spaces in .group 
            
            Averages.Calcium.after <-  Averages.Calcium[c(3:8),]  #subseting only for those with averages after flooding
            Averages.Calcium.before <- Averages.Calcium[c(1:2),]  #subsetting only for those before flooding because they were not tested
            #Averages.Calcium.after$Inundation <- factor(Averages.Calcium.after$Inundation, levels= c("None", "Typical", "Extreme"))
            newdata <-cld.Calcium [order(cld.Calcium$Inundation,cld.Calcium$Soil),] #Sorting so we can add the tukey-groupings to the raw means to be plotted
            Averages.Calcium.after.sorted <- Averages.Calcium.after[order(Averages.Calcium.after$Inundation, Averages.Calcium.after$Soil),]
            Averages.Calcium.after.sorted$.group <- newdata$.group
            Averages.Calcium.before$.group <- c("NA", "NA")
            
            Averages.Calcium.after.sorted.letters <- rbind(Averages.Calcium.after.sorted,Averages.Calcium.before)
            
            
            Calcium.plot2 <- ggplot(data =Averages.Calcium.after.sorted.letters, aes(x = Inundation, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Flooding Treatment", y= "Calcium (ppm)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -4,hjust=0.5,
                        color="black",
                        size=4)
            
            #testing main effects only
            lmCalcium.simple <- lm(Value ~ Inundation + Soil, data=Calcium)
            anova(lmCalcium.simple)
            lsmeansCalcium.soil <- lsmeans(lmCalcium.simple, pairwise ~ Soil)
            cld.Calcium.soil <- cld(lsmeansCalcium.soil$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Calcium.soil$Nutrient <- c("Calcium", "Calcium")
            
            
            lsmeansCalcium.flooding <- lsmeans(lmCalcium.simple, pairwise ~ Inundation)
            cld.Calcium.flooding <- cld(lsmeansCalcium.flooding$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Calcium.flooding$Nutrient <- c("Calcium", "Calcium", "Calcium")
            cld.Calcium.flooding$Inundation <- factor(cld.Calcium.flooding$Inundation, levels = c("None", "Typical", "Extreme"))
            
            
            Calcium.plot.soil <- ggplot(data =cld.Calcium.soil, aes(x = Soil, y=lsmean, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Calcium (ppm)") +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            palette <- c()
            Calcium.plot.flooding <- ggplot(data =cld.Calcium.flooding, aes(x = Inundation, y=lsmean, fill = Inundation)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey68", "grey92", "gray19")) +
              labs(x="Flooding Treatment", y= "Calcium (ppm)") +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Averages.Calcium.before <- Averages.Calcium.before[,c(1:7)]
            
            Calcium.plot.before <- ggplot(data = Averages.Calcium.before, aes(x = Soil, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              #geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Calcium (ppm)") +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Calcium.all <- ggarrange(Calcium.plot.before, Calcium.plot.soil,Calcium.plot.flooding, ncol=3,nrow=1, common.legend = TRUE, legend="top")
            
            
            
            
            #Model for Phosphorus:
            Phosphorus <- subset(Nu, Nutrient == "Phosphorus" & Inundation != "Before")
            Phosphorus$Nutrient <- factor(Phosphorus$Nutrient)
            
            lmPhosphorus <- lm(Value ~ Inundation*Soil, data=Phosphorus)
            anova(lmPhosphorus)
            plot(lmPhosphorus)
            qqPlot(Phosphorus$Value)
            hist(resid(lmPhosphorus,type="deviance"))
            qqnorm(resid(lmPhosphorus,type="deviance"))
            qqline(resid(lmPhosphorus,type="deviance"))   # Distribution of residuals is normal, proceed with analyses
            shapiro.test(resid(lmPhosphorus))
            
            lsmeansPhosphorus <- lsmeans(lmPhosphorus, pairwise ~ Inundation*Soil)
            #lsmeansPhosphorus.2 <- lsmeans(lmPhosphorus, pairwise ~ Soil|Inundation)
            cld.Phosphorus <- cld(lsmeansPhosphorus$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Phosphorus$Inundation <- factor(cld.Phosphorus$Inundation, levels=c("None", "Typical", "Extreme")) ### Order the levels for printing
            cld.Phosphorus$.group=gsub(" ", "", cld.Phosphorus$.group)  ###  Remove spaces in .group 
            
            Averages.Phosphorus.after <-  Averages.Phosphorus[c(3:8),]  #subseting only for those with averages after flooding
            Averages.Phosphorus.before <- Averages.Phosphorus[c(1:2),]  #subsetting only for those before flooding because they were not tested
            #Averages.Phosphorus.after$Inundation <- factor(Averages.Phosphorus.after$Inundation, levels= c("None", "Typical", "Extreme"))
            newdata <-cld.Phosphorus [order(cld.Phosphorus$Inundation,cld.Phosphorus$Soil),] #Sorting so we can add the tukey-groupings to the raw means to be plotted
            Averages.Phosphorus.after.sorted <- Averages.Phosphorus.after[order(Averages.Phosphorus.after$Inundation, Averages.Phosphorus.after$Soil),]
            Averages.Phosphorus.after.sorted$.group <- newdata$.group
            Averages.Phosphorus.before$.group <- c("NA", "NA")
            
            Averages.Phosphorus.after.sorted.letters <- rbind(Averages.Phosphorus.after.sorted,Averages.Phosphorus.before)
            
            Phosphorus.plot2 <- ggplot(data =Averages.Phosphorus.after.sorted.letters, aes(x = Inundation, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Flooding Treatment", y= "Phosphorus (ppm)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            #testing main effects:
            lmPhosphorus.simple <- lm(Value ~ Inundation + Soil, data=Phosphorus)
            anova(lmPhosphorus.simple)
            lsmeansPhosphorus.soil <- lsmeans(lmPhosphorus.simple, pairwise ~ Soil)
            cld.Phosphorus.soil <- cld(lsmeansPhosphorus.soil$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Phosphorus.soil$Nutrient <- c("Phosphorus", "Phosphorus")
            
            
            lsmeansPhosphorus.flooding <- lsmeans(lmPhosphorus.simple, pairwise ~ Inundation)
            cld.Phosphorus.flooding <- cld(lsmeansPhosphorus.flooding$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Phosphorus.flooding$Nutrient <- c("Phosphorus", "Phosphorus", "Phosphorus")
            cld.Phosphorus.flooding$Inundation <- factor(cld.Phosphorus.flooding$Inundation, levels = c("None", "Typical", "Extreme"))
            
            
            Phosphorus.plot.soil <- ggplot(data =cld.Phosphorus.soil, aes(x = Soil, y=lsmean, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Phosphorus (ppm)") +
              ylim(0,50) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Phosphorus.plot.flooding <- ggplot(data =cld.Phosphorus.flooding, aes(x = Inundation, y=lsmean, fill = Inundation)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey68", "grey92", "gray19")) +
              labs(x="Flooding Treatment", y= "Phosphorus (ppm)") +
              ylim(0,50) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Averages.Phosphorus.before <- Averages.Phosphorus.before[,c(1:7)]
            
            Phosphorus.plot.before <- ggplot(data = Averages.Phosphorus.before, aes(x = Soil, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              #geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Phosphorus (ppm)") +
              ylim(0,50) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Phosphorus.all <- ggarrange(Phosphorus.plot.before, Phosphorus.plot.soil,Phosphorus.plot.flooding, ncol=3,nrow=1, common.legend = TRUE, legend="top")
            
            
            
            
            
            #Model for Nitrogen:
            Nitrogen <- subset(Nu, Nutrient == "Nitrogen" & Inundation != "Before")
            Nitrogen$Nutrient <- factor(Nitrogen$Nutrient)
            
            lmNitrogen <- lm(Value ~ Inundation*Soil, data=Nitrogen)
            anova(lmNitrogen)
            plot(lmNitrogen)
            qqPlot(Nitrogen$Value)
            hist(resid(lmNitrogen,type="deviance"))
            qqnorm(resid(lmNitrogen,type="deviance"))
            qqline(resid(lmNitrogen,type="deviance"))   # Distribution of residuals is normal, proceed with analyses
            shapiro.test(resid(lmNitrogen))
            
            lsmeansNitrogen <- lsmeans(lmNitrogen, pairwise ~ Inundation*Soil)
            #lsmeansNitrogen.2 <- lsmeans(lmNitrogen, pairwise ~ Soil|Inundation)
            cld.Nitrogen <- cld(lsmeansNitrogen$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Nitrogen$Inundation <- factor(cld.Nitrogen$Inundation, levels=c("None", "Typical", "Extreme")) ### Order the levels for printing
            cld.Nitrogen$.group=gsub(" ", "", cld.Nitrogen$.group)  ###  Remove spaces in .group 
            
            Averages.Nitrogen.after <-  Averages.Nitrogen[c(3:8),]  #subseting only for those with averages after flooding
            Averages.Nitrogen.before <- Averages.Nitrogen[c(1:2),]  #subsetting only for those before flooding because they were not tested
            #Averages.Nitrogen.after$Inundation <- factor(Averages.Nitrogen.after$Inundation, levels= c("None", "Typical", "Extreme"))
            newdata <-cld.Nitrogen [order(cld.Nitrogen$Inundation,cld.Nitrogen$Soil),] #Sorting so we can add the tukey-groupings to the raw means to be plotted
            Averages.Nitrogen.after.sorted <- Averages.Nitrogen.after[order(Averages.Nitrogen.after$Inundation, Averages.Nitrogen.after$Soil),]
            Averages.Nitrogen.after.sorted$.group <- newdata$.group
            Averages.Nitrogen.before$.group <- c("NA", "NA")
            
            Averages.Nitrogen.after.sorted.letters <- rbind(Averages.Nitrogen.after.sorted,Averages.Nitrogen.before)
            
            Nitrogen.plot2 <- ggplot(data =Averages.Nitrogen.after.sorted.letters, aes(x = Inundation, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Flooding Treatment", y= "Nitrogen (%)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            #testing main effects:
            lmNitrogen.simple <- lm(Value ~ Inundation + Soil, data=Nitrogen)
            anova(lmNitrogen.simple)
            lsmeansNitrogen.soil <- lsmeans(lmNitrogen.simple, pairwise ~ Soil)
            cld.Nitrogen.soil <- cld(lsmeansNitrogen.soil$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Nitrogen.soil$Nutrient <- c("Nitrogen", "Nitrogen")
            
            
            lsmeansNitrogen.flooding <- lsmeans(lmNitrogen.simple, pairwise ~ Inundation)
            cld.Nitrogen.flooding <- cld(lsmeansNitrogen.flooding$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Nitrogen.flooding$Nutrient <- c("Nitrogen", "Nitrogen", "Nitrogen")
            cld.Nitrogen.flooding$Inundation <- factor(cld.Nitrogen.flooding$Inundation, levels = c("None", "Typical", "Extreme"))
            
            
            Nitrogen.plot.soil <- ggplot(data =cld.Nitrogen.soil, aes(x = Soil, y=lsmean, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Nitrogen (%)") +
              ylim(0, 0.4) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Nitrogen.plot.flooding <- ggplot(data =cld.Nitrogen.flooding, aes(x = Inundation, y=lsmean, fill = Inundation)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey68", "grey92", "gray19")) +
              labs(x="Flooding Treatment", y= "Nitrogen (%)") +
              ylim(0, 0.4) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Averages.Nitrogen.before <- Averages.Nitrogen.before[,c(1:7)]
            
            Nitrogen.plot.before <- ggplot(data = Averages.Nitrogen.before, aes(x = Soil, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              #geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              ylim(0,0.4) +
              labs(x="Soil Treatment", y= "Nitrogen (%)") +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Nitrogen.all <- ggarrange(Nitrogen.plot.before, Nitrogen.plot.soil,Nitrogen.plot.flooding, ncol=3,nrow=1, common.legend = TRUE, legend="top")
            
            
            #Model for Sulfur:
            Sulfur <- subset(Nu, Nutrient == "Sulfur" & Inundation != "Before")
            Sulfur$Nutrient <- factor(Sulfur$Nutrient)
            
            lmSulfur <- lm(Value ~ Inundation*Soil, data=Sulfur)
            anova(lmSulfur)
            plot(lmSulfur)
            qqPlot(Sulfur$Value)
            hist(resid(lmSulfur,type="deviance"))
            qqnorm(resid(lmSulfur,type="deviance"))
            qqline(resid(lmSulfur,type="deviance"))   # Distribution of residuals is normal, proceed with analyses
            shapiro.test(resid(lmSulfur))
            
            lsmeansSulfur <- lsmeans(lmSulfur, pairwise ~ Inundation*Soil)
            lsmeansSulfur.2 <- lsmeans(lmSulfur, pairwise ~ Inundation|Soil)
            lsmeansSulfur.3 <- lsmeans(lmSulfur, pairwise ~ Soil|Inundation)
            lsmeansSulfur.4 <- lsmeans(lmSulfur, pairwise ~ Soil)
            
            #lsmeansSulfur.2 <- lsmeans(lmSulfur, pairwise ~ Soil|Inundation)
            cld.Sulfur <- cld(lsmeansSulfur$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Sulfur$Inundation <- factor(cld.Sulfur$Inundation, levels=c("None", "Typical", "Extreme")) ### Order the levels for printing
            cld.Sulfur$.group=gsub(" ", "", cld.Sulfur$.group)  ###  Remove spaces in .group 
            
            Averages.Sulfur.after <-  Averages.Sulfur[c(3:8),]  #subseting only for those with averages after flooding
            Averages.Sulfur.before <- Averages.Sulfur[c(1:2),]  #subsetting only for those before flooding because they were not tested
            #Averages.Sulfur.after$Inundation <- factor(Averages.Sulfur.after$Inundation, levels= c("None", "Typical", "Extreme"))
            newdata <-cld.Sulfur [order(cld.Sulfur$Inundation,cld.Sulfur$Soil),] #Sorting so we can add the tukey-groupings to the raw means to be plotted
            Averages.Sulfur.after.sorted <- Averages.Sulfur.after[order(Averages.Sulfur.after$Inundation, Averages.Sulfur.after$Soil),]
            Averages.Sulfur.after.sorted$.group <- newdata$.group
            Averages.Sulfur.before$.group <- c("NA", "NA")
            
            Averages.Sulfur.after.sorted.letters <- rbind(Averages.Sulfur.after.sorted,Averages.Sulfur.before)
            
            Sulfur.plot2 <- ggplot(data =Averages.Sulfur.after.sorted.letters, aes(x = Inundation, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Flooding Treatment", y= "Sulfur (ppm)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            lmSulfur.simple <- lm(Value ~ Inundation + Soil, data=Sulfur)
            anova(lmSulfur.simple)
            lsmeansSulfur.soil <- lsmeans(lmSulfur.simple, pairwise ~ Soil)
            cld.Sulfur.soil <- cld(lsmeansSulfur.soil$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Sulfur.soil$Nutrient <- c("Sulfur", "Sulfur")
            
            
            lsmeansSulfur.flooding <- lsmeans(lmSulfur.simple, pairwise ~ Inundation)
            cld.Sulfur.flooding <- cld(lsmeansSulfur.flooding$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Sulfur.flooding$Nutrient <- c("Sulfur", "Sulfur", "Sulfur")
            cld.Sulfur.flooding$Inundation <- factor(cld.Sulfur.flooding$Inundation, levels = c("None", "Typical", "Extreme"))
            
            
            Sulfur.plot.soil <- ggplot(data =cld.Sulfur.soil, aes(x = Soil, y=lsmean, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Sulfur (ppm)") +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            palette <- c()
            Sulfur.plot.flooding <- ggplot(data =cld.Sulfur.flooding, aes(x = Inundation, y=lsmean, fill = Inundation)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey68", "grey92", "gray19")) +
              labs(x="Flooding Treatment", y= "Sulfur (ppm)") +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Averages.Sulfur.before <- Averages.Sulfur.before[,c(1:7)]
            
            Sulfur.plot.before <- ggplot(data = Averages.Sulfur.before, aes(x = Soil, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              #geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Sulfur (ppm)") +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Sulfur.all <- ggarrange(Sulfur.plot.before, Sulfur.plot.soil,Sulfur.plot.flooding, ncol=3,nrow=1, common.legend = TRUE, legend="top")
            
            #Model for Copper:
            Copper <- subset(Nu, Nutrient == "Copper" & Inundation != "Before")
            Copper$Nutrient <- factor(Copper$Nutrient)
            
            lmCopper <- lm(Value ~ Inundation*Soil, data=Copper)
            anova(lmCopper)
            plot(lmCopper)
            qqPlot(Copper$Value)
            hist(resid(lmCopper,type="deviance"))
            qqnorm(resid(lmCopper,type="deviance"))
            qqline(resid(lmCopper,type="deviance"))   # Distribution of residuals is normal, proceed with analyses
            shapiro.test(resid(lmCopper))
            
            lsmeansCopper <- lsmeans(lmCopper, pairwise ~ Inundation*Soil)
            cld.Copper <- cld(lsmeansCopper$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Copper$Inundation <- factor(cld.Copper$Inundation, levels=c("None", "Typical", "Extreme")) ### Order the levels for printing
            cld.Copper$.group=gsub(" ", "", cld.Copper$.group)  ###  Remove spaces in .group 
            
            Averages.Copper.after <-  Averages.Copper[c(3:8),]  #subseting only for those with averages after flooding
            Averages.Copper.before <- Averages.Copper[c(1:2),]  #subsetting only for those before flooding because they were not tested
            #Averages.Copper.after$Inundation <- factor(Averages.Copper.after$Inundation, levels= c("None", "Typical", "Extreme"))
            newdata <-cld.Copper [order(cld.Copper$Inundation,cld.Copper$Soil),] #Sorting so we can add the tukey-groupings to the raw means to be plotted
            Averages.Copper.after.sorted <- Averages.Copper.after[order(Averages.Copper.after$Inundation, Averages.Copper.after$Soil),]
            Averages.Copper.after.sorted$.group <- newdata$.group
            Averages.Copper.before$.group <- c("NA", "NA")
            
            Averages.Copper.after.sorted.letters <- rbind(Averages.Copper.after.sorted,Averages.Copper.before)
            
            Copper.plot2 <- ggplot(data =Averages.Copper.after.sorted.letters, aes(x = Inundation, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Flooding Treatment", y= "Copper (%)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -4,hjust=0.5,
                        color="black",
                        size=4)
            
            #testing main effects
            lmCopper.simple <- lm(Value ~ Inundation + Soil, data=Copper)
            anova(lmCopper.simple)
            lsmeansCopper.soil <- lsmeans(lmCopper.simple, pairwise ~ Soil)
            cld.Copper.soil <- cld(lsmeansCopper.soil$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Copper.soil$Nutrient <- c("Copper", "Copper")
            
            
            lsmeansCopper.flooding <- lsmeans(lmCopper.simple, pairwise ~ Inundation)
            cld.Copper.flooding <- cld(lsmeansCopper.flooding$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Copper.flooding$Nutrient <- c("Copper", "Copper", "Copper")
            cld.Copper.flooding$Inundation <- factor(cld.Copper.flooding$Inundation, levels = c("None", "Typical", "Extreme"))
            
            
            Copper.plot.soil <- ggplot(data =cld.Copper.soil, aes(x = Soil, y=lsmean, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Copper (ppm)") +
              ylim(0, 1.8) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Copper.plot.flooding <- ggplot(data =cld.Copper.flooding, aes(x = Inundation, y=lsmean, fill = Inundation)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey68", "grey92", "gray19")) +
              labs(x="Flooding Treatment", y= "Copper (ppm)") +
              ylim(0, 1.8) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Averages.Copper.before <- Averages.Copper.before[,c(1:7)]
            
            Copper.plot.before <- ggplot(data = Averages.Copper.before, aes(x = Soil, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              #geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Copper (ppm)") +
              ylim(0, 1.8) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Copper.all <- ggarrange(Copper.plot.before, Copper.plot.soil,Copper.plot.flooding, ncol=3,nrow=1, common.legend = TRUE, legend="top")
            
            
            #Model for Zinc:
            Zinc <- subset(Nu, Nutrient == "Zinc" & Inundation != "Before")
            Zinc$Nutrient <- factor(Zinc$Nutrient)
            
            lmZinc <- lm(Value ~ Inundation*Soil, data=Zinc)
            anova(lmZinc)
            plot(lmZinc)
            qqPlot(Zinc$Value)
            hist(resid(lmZinc,type="deviance"))
            qqnorm(resid(lmZinc,type="deviance"))
            qqline(resid(lmZinc,type="deviance"))   # Distribution of residuals is normal, proceed with analyses
            shapiro.test(resid(lmZinc))
            
            lsmeansZinc <- lsmeans(lmZinc, pairwise ~ Inundation*Soil)
            cld.Zinc <- cld(lsmeansZinc$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Zinc$Inundation <- factor(cld.Zinc$Inundation, levels=c("None", "Typical", "Extreme")) ### Order the levels for printing
            cld.Zinc$.group=gsub(" ", "", cld.Zinc$.group)  ###  Remove spaces in .group 
            
            Averages.Zinc.after <-  Averages.Zinc[c(3:8),]  #subseting only for those with averages after flooding
            Averages.Zinc.before <- Averages.Zinc[c(1:2),]  #subsetting only for those before flooding because they were not tested
            #Averages.Zinc.after$Inundation <- factor(Averages.Zinc.after$Inundation, levels= c("None", "Typical", "Extreme"))
            newdata <-cld.Zinc [order(cld.Zinc$Inundation,cld.Zinc$Soil),] #Sorting so we can add the tukey-groupings to the raw means to be plotted
            Averages.Zinc.after.sorted <- Averages.Zinc.after[order(Averages.Zinc.after$Inundation, Averages.Zinc.after$Soil),]
            Averages.Zinc.after.sorted$.group <- newdata$.group
            Averages.Zinc.before$.group <- c("NA", "NA")
            
            Averages.Zinc.after.sorted.letters <- rbind(Averages.Zinc.after.sorted,Averages.Zinc.before)
            
            Zinc.plot2 <- ggplot(data =Averages.Zinc.after.sorted.letters, aes(x = Inundation, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Flooding Treatment", y= "Zinc (%)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -4,hjust=0.5,
                        color="black",
                        size=4)
            
            #testing main effects
            lmZinc.simple <- lm(Value ~ Inundation + Soil, data=Zinc)
            anova(lmZinc.simple)
            lsmeansZinc.soil <- lsmeans(lmZinc.simple, pairwise ~ Soil)
            cld.Zinc.soil <- cld(lsmeansZinc.soil$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Zinc.soil$Nutrient <- c("Zinc", "Zinc")
            
            
            lsmeansZinc.flooding <- lsmeans(lmZinc.simple, pairwise ~ Inundation)
            cld.Zinc.flooding <- cld(lsmeansZinc.flooding$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Zinc.flooding$Nutrient <- c("Zinc", "Zinc", "Zinc")
            cld.Zinc.flooding$Inundation <- factor(cld.Zinc.flooding$Inundation, levels = c("None", "Typical", "Extreme"))
            
            
            Zinc.plot.soil <- ggplot(data =cld.Zinc.soil, aes(x = Soil, y=lsmean, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Zinc (ppm)") +
              ylim(0, 16) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Zinc.plot.flooding <- ggplot(data =cld.Zinc.flooding, aes(x = Inundation, y=lsmean, fill = Inundation)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey68", "grey92", "gray19")) +
              labs(x="Flooding Treatment", y= "Zinc (ppm)") +
              ylim(0, 16) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Averages.Zinc.before <- Averages.Zinc.before[,c(1:7)]
            
            Zinc.plot.before <- ggplot(data = Averages.Zinc.before, aes(x = Soil, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              #geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Zinc (ppm)") +
              ylim(0, 16) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Zinc.all <- ggarrange(Zinc.plot.before, Zinc.plot.soil,Zinc.plot.flooding, ncol=3,nrow=1, common.legend = TRUE, legend="top")
            
            
            
            
            #Model for Magnesium:
            Magnesium <- subset(Nu, Nutrient == "Magnesium" & Inundation != "Before")
            Magnesium$Nutrient <- factor(Magnesium$Nutrient)
            
            lmMagnesium <- lm(Value ~ Inundation*Soil, data=Magnesium)
            anova(lmMagnesium)
            plot(lmMagnesium)
            qqPlot(Magnesium$Value)
            hist(resid(lmMagnesium,type="deviance"))
            qqnorm(resid(lmMagnesium,type="deviance"))
            qqline(resid(lmMagnesium,type="deviance"))   # Distribution of residuals is normal, proceed with analyses
            shapiro.test(resid(lmMagnesium))
            
            lsmeansMagnesium <- lsmeans(lmMagnesium, pairwise ~ Inundation*Soil)
            cld.Magnesium <- cld(lsmeansMagnesium$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Magnesium$Inundation <- factor(cld.Magnesium$Inundation, levels=c("None", "Typical", "Extreme")) ### Order the levels for printing
            cld.Magnesium$.group=gsub(" ", "", cld.Magnesium$.group)  ###  Remove spaces in .group 
            
            Averages.Magnesium.after <-  Averages.Magnesium[c(3:8),]  #subseting only for those with averages after flooding
            Averages.Magnesium.before <- Averages.Magnesium[c(1:2),]  #subsetting only for those before flooding because they were not tested
            #Averages.Magnesium.after$Inundation <- factor(Averages.Magnesium.after$Inundation, levels= c("None", "Typical", "Extreme"))
            newdata <-cld.Magnesium [order(cld.Magnesium$Inundation,cld.Magnesium$Soil),] #Sorting so we can add the tukey-groupings to the raw means to be plotted
            Averages.Magnesium.after.sorted <- Averages.Magnesium.after[order(Averages.Magnesium.after$Inundation, Averages.Magnesium.after$Soil),]
            Averages.Magnesium.after.sorted$.group <- newdata$.group
            Averages.Magnesium.before$.group <- c("NA", "NA")
            
            Averages.Magnesium.after.sorted.letters <- rbind(Averages.Magnesium.after.sorted,Averages.Magnesium.before)
            
            Magnesium.plot2 <- ggplot(data =Averages.Magnesium.after.sorted.letters, aes(x = Inundation, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Flooding Treatment", y= "Magnesium (%)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -4,hjust=0.5,
                        color="black",
                        size=4)
            
            #testing main effects
            lmMagnesium.simple <- lm(Value ~ Inundation + Soil, data=Magnesium)
            anova(lmMagnesium.simple)
            lsmeansMagnesium.soil <- lsmeans(lmMagnesium.simple, pairwise ~ Soil)
            cld.Magnesium.soil <- cld(lsmeansMagnesium.soil$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Magnesium.soil$Nutrient <- c("Magnesium", "Magnesium")
            
            
            lsmeansMagnesium.flooding <- lsmeans(lmMagnesium.simple, pairwise ~ Inundation)
            cld.Magnesium.flooding <- cld(lsmeansMagnesium.flooding$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.Magnesium.flooding$Nutrient <- c("Magnesium", "Magnesium", "Magnesium")
            cld.Magnesium.flooding$Inundation <- factor(cld.Magnesium.flooding$Inundation, levels = c("None", "Typical", "Extreme"))
            
            
            Magnesium.plot.soil <- ggplot(data =cld.Magnesium.soil, aes(x = Soil, y=lsmean, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Magnesium (ppm)") +
              ylim(0, 800) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Magnesium.plot.flooding <- ggplot(data =cld.Magnesium.flooding, aes(x = Inundation, y=lsmean, fill = Inundation)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey68", "grey92", "gray19")) +
              labs(x="Flooding Treatment", y= "Magnesium (ppm)") +
              ylim(0, 800) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Averages.Magnesium.before <- Averages.Magnesium.before[,c(1:7)]
            
            Magnesium.plot.before <- ggplot(data = Averages.Magnesium.before, aes(x = Soil, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              #geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "Magnesium (ppm)") +
              ylim(0, 800) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            Magnesium.all <- ggarrange(Magnesium.plot.before, Magnesium.plot.soil,Magnesium.plot.flooding, ncol=3,nrow=1, common.legend = TRUE, legend="top")
            
            
            #Model for pH:
            pH <- subset(Nu, Nutrient == "pH" & Inundation != "Before")
            pH$Nutrient <- factor(pH$Nutrient)
            
            lmpH <- lm(Value ~ Inundation*Soil, data=pH)
            anova(lmpH)
            plot(lmpH)
            qqPlot(pH$Value)
            hist(resid(lmpH,type="deviance"))
            qqnorm(resid(lmpH,type="deviance"))
            qqline(resid(lmpH,type="deviance"))   # Distribution of residuals is normal, proceed with analyses
            shapiro.test(resid(lmpH))
            
            lsmeanspH <- lsmeans(lmpH, pairwise ~ Inundation*Soil)
            cld.pH <- cld(lsmeanspH$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.pH$Inundation <- factor(cld.pH$Inundation, levels=c("None", "Typical", "Extreme")) ### Order the levels for printing
            cld.pH$.group=gsub(" ", "", cld.pH$.group)  ###  Remove spaces in .group 
            
            Averages.pH.after <-  Averages.pH[c(3:8),]  #subseting only for those with averages after flooding
            Averages.pH.before <- Averages.pH[c(1:2),]  #subsetting only for those before flooding because they were not tested
            #Averages.pH.after$Inundation <- factor(Averages.pH.after$Inundation, levels= c("None", "Typical", "Extreme"))
            newdata <-cld.pH [order(cld.pH$Inundation,cld.pH$Soil),] #Sorting so we can add the tukey-groupings to the raw means to be plotted
            Averages.pH.after.sorted <- Averages.pH.after[order(Averages.pH.after$Inundation, Averages.pH.after$Soil),]
            Averages.pH.after.sorted$.group <- newdata$.group
            Averages.pH.before$.group <- c("NA", "NA")
            
            Averages.pH.after.sorted.letters <- rbind(Averages.pH.after.sorted,Averages.pH.before)
            
            pH.plot2 <- ggplot(data =Averages.pH.after.sorted.letters, aes(x = Inundation, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Flooding Treatment", y= "pH (%)") +
              theme(panel.spacing = unit(0.5, "lines"), 
                    strip.background = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face="bold", vjust=1.0, size = 14),
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -4,hjust=0.5,
                        color="black",
                        size=4)
            
            #testing main effects
            lmpH.simple <- lm(Value ~ Inundation + Soil, data=pH)
            anova(lmpH.simple)
            lsmeanspH.soil <- lsmeans(lmpH.simple, pairwise ~ Soil)
            cld.pH.soil <- cld(lsmeanspH.soil$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.pH.soil$Nutrient <- c("pH", "pH")
            
            
            lsmeanspH.flooding <- lsmeans(lmpH.simple, pairwise ~ Inundation)
            cld.pH.flooding <- cld(lsmeanspH.flooding$lsmeans, alpha=0.05, Letters=letters, adjust="tukey")
            cld.pH.flooding$Nutrient <- c("pH", "pH", "pH")
            cld.pH.flooding$Inundation <- factor(cld.pH.flooding$Inundation, levels = c("None", "Typical", "Extreme"))
            
            
            pH.plot.soil <- ggplot(data =cld.pH.soil, aes(x = Soil, y=lsmean, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "pH") +
              ylim(0, 6) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            pH.plot.flooding <- ggplot(data =cld.pH.flooding, aes(x = Inundation, y=lsmean, fill = Inundation)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey68", "grey92", "gray19")) +
              labs(x="Flooding Treatment", y= "pH") +
              ylim(0, 6) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
              geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -2.5,hjust=0.5,
                        color="black",
                        size=4)
            
            Averages.pH.before <- Averages.pH.before[,c(1:7)]
            
            pH.plot.before <- ggplot(data = Averages.pH.before, aes(x = Soil, y=Value, fill = Soil)) + 
              geom_bar(colour= "black", stat= "identity", size=0.5, position=position_dodge()) +
              #geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9), na.rm = TRUE) +
              #facet_grid(~CollectionTime, switch = "x", scales = "free_x", space = "free_x") +
              scale_fill_manual(values = c("grey39", "white")) +
              labs(x="Soil Treatment", y= "pH") +
              ylim(0, 6) +
              theme(axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=12, colour ="black", angle=360,hjust=0.5),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            pH.all <- ggarrange(pH.plot.before, pH.plot.soil,pH.plot.flooding, ncol=3,nrow=1, common.legend = TRUE, legend="top")
            
            #Plotting all figures:
            All.nutrients.single.effects <- ggarrange(Carbon.all,Phosphorus.all,Nitrogen.all,Sulfur.all,Magnesium.all,Copper.all,Zinc.all, pH.all, ncol=1,nrow=8, common.legend = TRUE, legend="top", labels= c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)","(g)", "(h)"))
            ggsave("All.nutrients.single.effects.pdf", plot=All.nutrients.single.effects , device="pdf", scale = 1, width = 20, height = 60, units = "cm", dpi = 300)
            
            
            
            
            
            
            
            
            
            
            
            
            library(gridExtra)
            library(grid)
            install.packages("ggpubr")
            library(ggpubr)
            
            Nutrients <- ggarrange(Cplot,Nplot,Phosphorusplot,Splot, Znplot,Coplot,Mgplot, pHsplot,ncol=2,nrow=4, common.legend = TRUE, legend="top")
            ggsave("Nutrients.pdf", plot=Nutrients, device="pdf", scale = 1, width = 30, height = 60, units = "cm", dpi = 300)
            
            
            Nutrients.Micro <- ggarrange(Sodium.plot, Zinc.plot,Copper.plot, ncol=1,nrow=3, common.legend = TRUE, legend="top")
            ggsave("Nutrients.Micro.pdf", plot=Nutrients.Micro , device="pdf", scale = 1, width = 20, height = 25, units = "cm", dpi = 300)
            
            
            # Subsetting data for after soil to test for differences between sterile and non-sterile soil:
            after <- subset(Nu2, InocTime == "after")
            
            #paired-t-test between sterile and non-sterile treatments after inoculation:
            install.packages("PairedData") 
            library(PairedData)
            
            Ct.test <- t.test(Calcium..  ~ Inoculum.treatment, data = after, paired = TRUE)
            Ct.test 
            Nt.test <- t.test(Nitrogen..  ~ Inoculum.treatment, data = after, paired = TRUE)
            Nt.test 
            pH.test <- t.test(pH..1.1.Water.  ~ Inoculum.treatment, data = after, paired = TRUE)
            pH.test 
            Co.test <- t.test(Copper..ppm  ~ Inoculum.treatment, data = after, paired = TRUE)
            Co.test 
            Phos.test <- t.test(Phosphorus..ppm  ~ Inoculum.treatment, data = after, paired = TRUE)
            Phos.test
            Mg.test <- t.test(Magnesium..ppm  ~ Inoculum.treatment, data = after, paired = TRUE)
            Mg.test
            Zn.test <- t.test(Zinc..ppm  ~ Inoculum.treatment, data = after, paired = TRUE)
            Zn.test 
            S.test <- t.test(Sulfur..ppm  ~ Inoculum.treatment, data = after, paired = TRUE)
            S.test 
            
            
                     
            
###############################################################################
            ###### MAIN MANUSCRIPT FIGURES
################################################################################
            
            # To change order of Inundation levels
            GAF$Hydrology <-factor(GAF$Inundation,levels=levels(GAF$Inundation)[c(2,1,3)])
            P$Hydrology <-factor(P$Inundation,levels=levels(P$Inundation)[c(2,1,3)])
            
            
            ## TCSH PLOT
            
            #Need to recalculate lsmeans
            glmm_TCSH <- glmer(TCSH.t ~ Hydrology + Soil + Hydrology*Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            library(lsmeans)
            lmmTCSH.rg1 <- ref.grid(glmm_TCSH, type="response") # here we type = "response" to obtain back-transformed lsmeans
            lsmeansTCSH <- summary(lmmTCSH.rg1)
            
            # Making plot in ggplot2
            sumTCSHplot <- ggplot(lsmeansTCSH, aes(x= Soil, y= response, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "TCSH (cm)") 
            p5 <- sumTCSHplot + theme_bw()
            #p5.1 <- p5 + theme(legend.position="none")
            p5.3 <- p5 + ylim(0,5)
            p5.2 <- p5.3 + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                 axis.text.x = element_blank(),
                                 axis.text.y = element_text(size=14, colour ="black"))
            pTCSH <- p5.2 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                  panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            ggsave("pTCSH.pdf", plot=pTCSH, device="pdf", scale = 1, width = 20, height = 20, units = "cm", dpi = 300)
            
            
            
            ########## SLA plot
            
            #Need to recalculate lsmeans
            lmmSLAFL <- lmer(SLAFL ~ Hydrology + Soil + Soil*Hydrology +  (1 | Seed), data = P,
                             REML = FALSE) # Full model for SLAFL
            library(lsmeans)
            lmmSLA.rg1 <- ref.grid(lmmSLAFL) 
            lsmeansSLA <- summary(lmmSLA.rg1)
            
            # plot
            SLAFLplot <- ggplot(lsmeansSLA, aes(x= Soil, y= prediction, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=prediction-SE, ymax=prediction+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y="SLA (m2 Kg-1)") 
            p11 <- SLAFLplot + theme_bw() 
            p12 <- p11 + ylim(0,2.5)
            p13 <- p12 + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                               axis.title.x = element_blank(),
                               axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                               axis.text.x = element_blank(),
                               axis.text.y = element_text(size=14, colour ="black"))
            pSLAFLplot <- p13 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                      panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            ggsave("pSLAFLplot.pdf", plot=pSLAFLplot, device="pdf", scale = 1, width = 20, height = 20, units = "cm", dpi = 300)
            
            
            #### TOTAL DRY BIOMASS PLOT (TDB)
            
            #Need to recalculate lsmeans
            glmm_TDB2 <- glmer(TDB ~  Soil + Hydrology + Soil*Hydrology + (1 | Seed), data = P, family = gaussian(link = "log")) 
            lmmTDB.rg1 <- ref.grid(glmm_TDB2, type="response") # here we type = "response" to obtain back-transformed lsmeans
            lsmeansTDB <- summary(lmmTDB.rg1)
            
            # Plot TDB
            TDBplot <- ggplot(lsmeansTDB, aes(x= Soil, y= response, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "Total Dry Biomass (gr)") 
            TDBplot1 <- TDBplot + ylim(0,3)
            pTDB <- TDBplot1 + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                     axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                     axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                     axis.text.x = element_text(size=14, colour ="black"),
                                     axis.text.y = element_text(size=14, colour ="black"))
            p6 <- pTDB + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pTDW <- p6 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                               panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            ggsave("pTDW.pdf", plot=pTDW, device="pdf", scale = 1, width = 20, height = 20, units = "cm", dpi = 300)
            
            ################ PLOT A:B ratio 
            
            #Need to recalculate lsmeans
            library(lsmeans)
            glmm_ABratio <- glmer(DryBiomassRatio ~ Hydrology + Soil + Hydrology*Soil + (1 | Seed), data = P, family = gaussian(link = "log"))  
            lmmABratio.rg1 <- ref.grid(glmm_ABratio, type="response") # here we type = "response" to obtain back-transformed lsmeans
            lsmeansABratio <- summary(lmmABratio.rg1)
            
            #Plotting A:B ratio
            
            DryBiomassRatioplot <- ggplot(lsmeansABratio, aes(x= Soil, y= response, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "A:B ratio (gr)") 
            p8 <- DryBiomassRatioplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                              axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                              axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                              axis.text.x = element_text(size=14, colour ="black"),
                                              axis.text.y = element_text(size=14, colour ="black"))
            p9 <- p8 + ylim(0,3)
            p10 <- p9 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pRatio <- p10 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                  panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            # **Plotting all graphs together with a shared legend  
            Growth <- ggarrange(pIH, pFH, pTCSH, ncol=3, common.legend = TRUE, legend="top")
            #grid.arrange(pIH, pFH, pTCSH, pSLAFLplot, ncol=4)
            ggsave("Growth.pdf", plot=Growth, device="pdf", scale = 1, width = 20, height = 10, units = "cm", dpi = 300)
            
            Figure4 <- ggarrange(pTCSH,pSLAFLplot,pTDW, pRatio, ncol=2,nrow=2, common.legend = TRUE, legend="top")
            ggsave("Figure3.pdf", plot=Figure3, device="pdf", scale = 1, width = 20, height = 20, units = "cm", dpi = 300)
            
            
            
##### TDSE & TAMF for Inoculant only
            
            RootFungi <- read.csv("RootFungi.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 
            RootFungi$Hydrology <- factor(RootFungi$Inundation,levels=levels(RootFungi$Inundation)[c(2,1,3)]) # To arrange inundation at specific order
            sumRootFungi = summarySE(data=RootFungi, measurevar= "Colonization", groupvars = c("Fungi", "Hydrology"))
            sumRootFungiType = summarySE(data=RootFungi, measurevar= "Colonization", groupvars = c("Fungi"))
            
   dev.off()         
            
            Fungiplot <- ggplot(sumRootFungi, aes(x= Hydrology, y= Colonization, fill= Fungi)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=Colonization-se, ymax=Colonization+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              #scale_fill_manual(values = c("#999999ff","#f2f2f7ff","#3a3336ff")) +
              scale_fill_manual(values = c("#3a3336ff","gray")) +
              labs(x="Flooding treatments", y= "Root Colonization (%)") +
              ylim(0,0.45) +
              scale_x_discrete(labels=c("None", "Typical", "Extreme")) +
           theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                    axis.text.x = element_text(size=14, colour ="black"),
                                    axis.text.y = element_text(size=14, colour ="black"),
                                    legend.position = "top") 
            
            non.sterile2 <- non.sterile[,c(3,15,18,24)]
            fungi.hostperformance <- gather(non.sterile2, Fungi, Colonization, -TDB, -Inundation) #This reshapes the matrix for Strain
            Fungi.regression.plot <-ggplot(fungi.hostperformance, aes(x=Colonization, y=TDB, color=Fungi, shape=Fungi)) +
              geom_point() +
              scale_shape_manual(values=c(21,17), breaks=c("TAMF", "TDSE"), labels=c("AMF", "DSE")) +
              scale_colour_manual(values=c("#3a3336ff","gray"), breaks=c("TAMF", "TDSE"), labels=c("AMF", "DSE")) + 
              geom_smooth(method=lm,se=FALSE) + 
              labs(x="Colonization (%)", y="Plant total dry biomass (g)") +
              theme(panel.background = element_rect(fill = "white", colour = "black"), 
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=16, colour ="black"),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 
            
            lsmean.fungi.str = lsmeans(lmfungi.structures,pairwise ~ Inundation*Fungi, adjust="tukey")
            fungi.cld <- cld(lsmean.fungi.str$lsmeans,alpha=0.05, Letters=letters) 
            fungi.cld$Inundation <- factor(fungi.cld$Inundation, levels=c("None","Typical","Extreme")) ### Order the levels for printing
            fungi.cld$Fungi <- factor(fungi.cld$Fungi, levels=c("ARB_percentage", "HYP_percentage", "VES_percentage", "MS_percentage", "SHYP_percentage"))
            fungi.cld$.group=gsub(" ", "", fungi.cld$.group)  ###  Remove spaces in .group 
            
            lsmean.fungi.str.1 = lsmeans(lmfungi.structures,pairwise ~ Inundation|Fungi, adjust="tukey")
            fungi.cld.1 <- cld(lsmean.fungi.str.1$lsmeans,alpha=0.05, Letters=letters) 
            
            lsmean.fungi.str.2 = lsmeans(lmfungi.structures,pairwise ~ Fungi|Inundation, adjust="tukey")
            fungi.cld.1 <- cld(lsmean.fungi.str.1$lsmeans,alpha=0.05, Letters=letters) 
            
            fungi.diver.plot.2 <- ggplot(data = fungi.cld, aes(x = Inundation, y = lsmean, fill = Fungi)) + 
              geom_bar(stat="identity", position=position_dodge(), colour="black") +
              geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              #facet_grid(~Soil, switch = "x", scales = "free_x", space = "free_x") +
              #scale_fill_manual(values =cbPalette) +
              scale_fill_grey(start=0.2, end = 0.8, name="Fungal structures",
                              breaks=c("ARB_percentage", "HYP_percentage", "VES_percentage", "MS_percentage", "SHYP_percentage"),
                              labels=c("ARB", "HYP", "VES", "MS", "SHYP")) +
              labs(x="Flooding Treatments", y= "Colonization(%)") +
              theme(panel.background = element_rect(fill = "white", colour = "black"), 
                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                    axis.text.x = element_text(size=16, colour ="black"),
                    axis.text.y = element_text(size=16, colour ="black"),
                    panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) + 
            geom_text(aes(label = .group),
                        position = position_dodge(0.9),
                        vjust = -4.5,hjust=0.5,
                        color="black",
                        size=4)
            
            
            
           
            Figure4 <- ggarrange(fungi.diver.plot.2,Fungiplot, Fungi.regression.plot, ncol=3, common.legend = FALSE, legend="top", labels=c("(a)", "(b)", "(c)"))
            ggsave("Figure4_RootFungi_allplots.pdf", plot=Figure4, device="pdf", scale = 1, width = 40, height = 15, units = "cm", dpi = 300)
            
            
            
            
            
            
                     
###############################################################################
            ###### SUPPLEMENTAL MATERIALS ANALYSES
################################################################################          
           
            
            head(GAF)
#IH - models
            lmmIH <- lmer(IH ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = GAF,
                          REML = FALSE) # Simple model for IH
            summary(lmmIH)
            Anova(lmmIH)
            
            #* testing the effect adding inundation to the model
            lmmIH.null <- lmer(IH ~ Soil +  (1 | Seed), data = GAF,
                               REML = FALSE) # Simple model for IH
            
            lmmIH.inu <- lmer(IH ~ Soil + Inundation +  (1 | Seed), data = GAF,
                              REML = FALSE) # Adding inundation
            
            anova (lmmIH.null,lmmIH.inu) # Likelihood ratio test to demonstrate no effect of inundation
            
            #* testing the effect adding soil to the model
            lmmIH.0 <- lmer(IH ~ 1 + (1 | Seed), data = GAF,
                            REML = FALSE) # Simple model for IH
            
            lmmIH.soil <- lmer(IH ~ Soil + (1 | Seed), data = GAF,
                               REML = FALSE) # Adding inundation
            
            anova (lmmIH.0,lmmIH.soil)
            
            #* testing the effect of random factors in the model
            
            lmmIH.0 <- lmer(IH ~ Soil + Hydrology + (1 | Seed), data = GAF,
                            REML = FALSE) # Simple model for IH
            
            lmmIH.1 <- lmer(IH ~ Soil + Hydrology + , data = GAF,
                            REML = FALSE)
            
            anova(lmmIH.1,lmmIH.0)
            
            
            # lsmeans IH
            library(lsmeans)
            lmmIH.rg1 <- ref.grid(lmmIH)
            summary(lmmIH.rg1)
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            IH.lsm <- lsmeans(lmmIH, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            cld(IH.lsm, alpha = .50) # To obtain groupings! 
            
            
#FH -model
            lmmFH <- lmer(FH ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = GAF,
                          REML = FALSE) # Simple model FH
            summary(lmmFH)
            Anova(lmmFH)
            
            #* testing the effect adding inundation to the model
            lmmFH.null <- lmer(FH ~ Soil +  (1 | Seed), data = GAF,
                               REML = FALSE) # Simple model for IH
            
            lmmFH.inu <- lmer(FH ~ Soil + Inundation +  (1 | Seed), data = GAF,
                              REML = FALSE) # Adding inundation
            
            anova (lmmFH.null,lmmFH.inu) # Likelihood ratio test to demonstrate no effect of inundation
            
            #* testing the effect adding soil to the model
            lmmFH.0 <- lmer(FH ~ 1 + (1 | Seed), data = GAF,
                            REML = FALSE) # Simple model for IH
            
            lmmFH.soil <- lmer(FH ~ Soil + (1 | Seed), data = GAF,
                               REML = FALSE) # Adding inundation
            
            anova (lmmFH.0,lmmFH.soil)
            summary(lmmFH.soil)
            
            
            # easy backward selection
            
            stimm.step1 <- step(lmmFH) #backward step model selection -random factor dropped
            stimm.step1
            
            #lsmeans FH
            library(lsmeans)
            lmmFH.rg1 <- ref.grid(lmmFH)
            summary(lmmFH.rg1)
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            FH.lsm <- lsmeans(lmmFH, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            cld(FH.lsm, alpha = .50) # To obtain groupings! 
            
            
            
#######    RTIHN --- GLMM ###################

head(GAF)
            
GAF$RHIN.t <- GAF$RHIN + 1
            
            
            #full model including random effect but using Laplace approximation instead 
            glmm_RHIN.t <- glmer(RHIN.t ~ Inundation + Soil + Inundation*Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            summary(glmm_RHIN.t)
            Anova(glmm_RHIN.t)
            
            #*testing relevance of inundation and interaction term
            glmm_RHIN.t0 <- glmer(RHIN.t ~  Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            glmm_RHIN.t1 <- glmer(RHIN.t ~  Soil + Inundation + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            glmm_RHIN.t2 <- glmer(RHIN.t ~  Soil + Inundation + Soil*Inundation + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            anova(glmm_RHIN.t0,glmm_RHIN.t1)
            anova(glmm_RHIN.t0,glmm_RHIN.t1,glmm_RHIN.t2, test="Chisq") 
            
            #*testing importance of random effect
            RHIN.nrandom <- glm(RHIN.t ~ Inundation + Soil + Inundation*Soil, family = gaussian(link = "log"),data = GAF) 
            RHIN.random <- glmer(RHIN.t ~ Inundation + Soil + Inundation*Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            anova(RHIN.random, RHIN.nrandom,test="Chisq") #* X2=0, P=1
            #*Posthoc tests
            lsmeans(glmm_RHIN.t2, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            lsmeans(glmm_RHIN.t2, pairwise ~ Inundation|Soil, adjust ="tukey") #within soil treatments among hydrological regimes
            
            
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            RIHN.lsm <- lsmeans(glmm_RHIN.t, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            cld(RIHN.lsm, alpha = .50) # To obtain groupings!

            
#######    RHOAF --- GLMM ###################
            
            GAF$RHOAF.t
            
            GAF$RHOAF.t <- GAF$RHOAF + 1
            
            #full model including random effect but using Laplace approximation instead 
            glmm_RHOAF.t <- glmer(RHOAF.t ~ Inundation + Soil + Inundation*Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            summary(glmm_RHOAF.t)
            Anova(glmm_RHOAF.t)
            
            #*testing relevance of inundation and interaction term
            glmm_RHOAF.t0 <- glmer(RHOAF.t ~  Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            glmm_RHOAF.t1 <- glmer(RHOAF.t ~  Soil + Inundation + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            glmm_RHOAF.t2 <- glmer(RHOAF.t ~  Soil + Inundation + Soil*Inundation + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            anova(glmm_RHOAF.t0,glmm_RHOAF.t1)
            anova(glmm_RHOAF.t0,glmm_RHOAF.t1,glmm_RHOAF.t2, test="Chisq") 
            
            #*testing importance of random effect
            RHOAF.nrandom <- glm(RHOAF.t ~ Inundation + Soil + Inundation*Soil, family = gaussian(link = "log"),data = GAF) 
            RHOAF.random <- glmer(RHOAF.t ~ Inundation + Soil + Inundation*Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            anova(RHIN.random, RHIN.nrandom,test="Chisq") #* X2=0, P=1
            #*Posthoc tests
            lsmeans(glmm_RHOAF.t2, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            lsmeans(glmm_RHOAF.t2, pairwise ~ Inundation|Soil, adjust ="tukey") #within soil treatments among hydrological regimes
           
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            RHOAF.lsm <- lsmeans(glmm_RHOAF.t, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            cld(RHOAF.lsm, alpha = .50) # To obtain groupings! 
            
           
####### Abroveground biomass (AB) 
            
            P$LogADW <- log(P$ADW)
            lmmADW <- lmer(LogADW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = TRUE)
            Anova(lmmADW)
            
            #*LSmeans
            lsmeans(lmmADW, pairwise ~ Inundation*Soil, adjust ="tukey")
            lsmeans(lmmADW, pairwise ~ Inundation|Soil, adjust ="tukey")
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            ADW.lsm <- lsmeans(lmmADW, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            cld(ADW.lsm, alpha = .50) # To obtain groupings! 

#### Belowground Biomass (BB) 
            
            P$LogBDW <- log(P$BDW)
            lmmBDW <- lmer(LogBDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = TRUE)
            Anova(lmmBDW, type ="III")
            
            #*LSmeans
            lsmeans(lmmBDW, pairwise ~ Inundation*Soil, adjust ="tukey")
            lsmeans(lmmBDW, pairwise ~ Inundation|Soil, adjust ="tukey")
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            BDW.lsm <- lsmeans(lmmBDW, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
            cld(BDW.lsm, alpha = .50) # To obtain groupings! 
            
###### Leaf biomass
            
P$LogLDW <- log(P$LDW)
lmmLDW <- lmer(LogLDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = TRUE)
Anova(lmmLDW, type ="III")

#*LSmeans
lsmeans(lmmLDW, pairwise ~ Inundation*Soil, adjust ="tukey")
lsmeans(lmmLDW, pairwise ~ Inundation|Soil, adjust ="tukey")

#*Obtaining Tukey-groupings:
#-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
LDW.lsm <- lsmeans(lmmLDW, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
cld(LDW.lsm, alpha = .50) # To obtain groupings! 


##### Stem biomass
P$LogSDW <- log(P$SDW)
lmmSDW <- lmer(LogSDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
               REML = TRUE)
Anova(lmmSDW, type = "III")

          #*LSmeans
          lsmeans(lmmSDW, pairwise ~ Inundation*Soil, adjust ="tukey")
          lsmeans(lmmSDW, pairwise ~ Inundation|Soil, adjust ="tukey")

          #*Obtaining Tukey-groupings:
          #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
          SDW.lsm <- lsmeans(lmmSDW, pairwise ~ Inundation*Soil, adjust ="tukey") #all pairwise comparisons
          cld(SDW.lsm, alpha = .50) # To obtain groupings! 

            
###############################################################################
            ###### SUPPLEMENTAL MATERIALS FIGURES
################################################################################

# To change order of Inundation levels
GAF$Hydrology <-factor(GAF$Inundation,levels=levels(GAF$Inundation)[c(2,1,3)])
P$Hydrology <-factor(P$Inundation,levels=levels(P$Inundation)[c(2,1,3)])

            
            
            ## INITIAL HEIGHT PLOT
            
            # Need to re-calculate lsmeans from new re-organized data set
            lmmIH <- lmer(IH ~ Hydrology + Soil + Soil*Hydrology +  (1 | Seed), data = GAF,REML = FALSE) # Simple model for IH
            lmmIH.rg1 <- ref.grid(lmmIH)
            lsmeansHI <- summary(lmmIH.rg1) #This way ggplot2 can recognize the output as a dataframe
            
            # Making plot in ggplot2
            sumIHplot <- ggplot(lsmeansHI, aes(x= Soil, y= prediction, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=prediction-SE, ymax=prediction+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "Initial Height (cm)") 
            IHplot <- sumIHplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                        axis.text.x = element_text(size=14, colour ="black"),
                                        axis.text.y = element_text(size=14, colour ="black"),
                                        legend.position = "top")
            p7 <-  IHplot + scale_x_discrete(name = "Soil type", labels = c("IN" = "Non sterile", "NO"= "Sterile"))
            #IHplot1<- IHplot + theme(legend.position="none")
            IHplot2 <- p7 + ylim(0,30) 
            pX <- IHplot2 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pIH <- pX + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                              panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            ## FINAL HEIGHT PLOT
            
            # Need to re-calculate lsmeans from new re-organized data set
            lmmFH <- lmer(FH ~ Hydrology + Soil + Soil*Hydrology +  (1 | Seed), data = GAF,REML = FALSE) # Simple model for IH
            lmmFH.rg1 <- ref.grid(lmmFH)
            lsmeansFI <- summary(lmmFH.rg1) #This way ggplot2 can recognize the output as a dataframe
            
            ### FINAL HEIGHT PLOT
            sumFHplot <- ggplot(lsmeansFI, aes(x= Soil, y= prediction, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=prediction-SE, ymax=prediction+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "Final Height (cm)") 
            FHplot <- sumFHplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                        axis.text.x = element_text(size=14, colour ="black"),
                                        axis.text.y = element_text(size=14, colour ="black"),
                                        legend.position = "top")
            p7 <-  FHplot + scale_x_discrete(name = "Soil type", labels = c("IN" = "Non sterile", "NO"= "Sterile"))
            #FHplot1<- FHplot + theme(legend.position="none")
            FHplot2 <- p7 + ylim(0,30) 
            pX <- FHplot2 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pFH <- pX + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                              panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            # RHOAF PLOT
            
            # Need to re-calculate lsmeans from new re-organized data set
            library(lsmeans)
            PQL_TCSH <- glmmPQL(TCSH.t ~ Hydrology + Soil + Hydrology*Soil, ~1 | Seed, family = gaussian(link = "log"),
                                data = GAF, verbose = FALSE) #Note that I am using TCSH.t due to zeros
            
            lmmTCSH.rg1 <- ref.grid(glmm_TCSH2, type="response") # here we type = "response" to obtain back-transformed lsmeans
            summary(lmmTCSH.rg2)
            
            sumRHOAFplot <- ggplot(lsmeansTCSH, aes(x= Soil , y= response, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Inoculation Treatments", y= "RHOAF (mm)") 
            p6 <- sumRHOAFplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                       axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                       axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                       axis.text.x = element_text(size=14, colour ="black"),
                                       axis.text.y = element_text(size=14, colour ="black"),
                                       legend.position = "top")
            p7 <- p6 + scale_x_discrete(name = "Soil type", labels = c("IN" = "Non sterile", "NO"= "Sterile"))
            p12 <- p7 + ylim(0,5.0)
            p13 <- p12 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pRHOAF <- p13 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                  panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            # RHIN PLOT
            
            # Need to re-calculate lsmeans from new re-organized data set
            ????
              
              
              sumRHINplot <- ggplot(sumRHIN, aes(x= Soil , y= RHIN, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=RHIN-se, ymax=RHIN+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("white","light gray","dark gray")) +
              labs(x="Inoculation Treatments", y= "RHIN (mm)") 
            p5 <- sumRHINplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                      axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                      axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                      axis.text.x = element_text(size=14, colour ="black"),
                                      axis.text.y = element_text(size=14, colour ="black"),
                                      legend.position = "top")
            p7 <- p5 + scale_x_discrete(name = "Soil type", labels = c("IN" = "Non sterile", "NO"= "Sterile"))
            p12 <- p7 + ylim(0,3.0)
            p13 <- p12 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pRHIN <- p13 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                 panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            
#################################### FUNGI ###########################################################################

E$Hydrology <-factor(E$Inundation,levels=levels(E$Inundation)[c(2,1,3)])           
#*Re-calculating raw means not categories- basic percentages
            sumTAMF = summarySE(data=E, measurevar = "TAMF", groupvars = c("Soil", "Hydrology"))
            sumTDSE = summarySE(data=E, measurevar = "TDSE", groupvars = c("Soil", "Hydrology"))
            sumTEndo = summarySE(data=E, measurevar = "Tendo", groupvars = c("Soil", "Hydrology"))
            
#### TAMF ######
            
            #* Plot
            TAMFplot <- ggplot(sumTAMF, aes(x= Soil, y= TAMF, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=TAMF-se, ymax=TAMF+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "AMF colonization (%)") 
            p8 <- TAMFplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                   axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                   axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                   axis.text.x = element_text(size=14, colour ="black"),
                                   axis.text.y = element_text(size=14, colour ="black"))
            p9 <- p8 + ylim(0,0.45)
            p10 <- p9 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pTAMF <- p10 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                 panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
#### TDSE
            
            TDSEplot <- ggplot(sumTDSE, aes(x= Soil, y= TDSE, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=TDSE-se, ymax=TDSE+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "DSE colonization (%)") 
            p11 <- TDSEplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                    axis.text.x = element_text(size=14, colour ="black"),
                                    axis.text.y = element_text(size=14, colour ="black"))
            p12 <- p11 + ylim(0,0.45)
            p13 <- p12 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pTDSE <- p13 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                 panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            
    # **Plotting all graphs together with a shared legend  
            FIGS3 <- ggarrange(pTAMF, pTDSE, ncol=2, common.legend = TRUE, legend="top")
            #grid.arrange(pIH, pFH, pTCSH, ncol=3)
            ggsave("FIGS3.pdf", plot=FIGS3, device="pdf", scale = 1, width = 20, height = 10, units = "cm", dpi = 300)
            
