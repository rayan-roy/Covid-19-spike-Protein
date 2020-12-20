#Rcode for Protein Train data

options(max.print = 20000000)
set.seed(12345678)
library(MASS)
library(car)
library(dplyr)

trainprotein <- read.csv("Stat331/Final\ Project/protein-train.csv")

#####################################################

# ************ Exploratory Data Analysis ************

# Box Plot and histogram of Angles
boxplot(trainprotein$angles,)
hist(trainprotein$angles, xlab="Angles", main = "Histogram of Angles")

# Finding the distribution of aromatic and non aromatic atom types
aromatic_count <- ncol(select(trainprotein,contains("aromatic")))
non_count <-684 - aromatic_count
slices <- c(aromatic_count, non_count)
lbls <- c("aromatic", "non-aromatic")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Percentage of Aromatic")

# Distribution of various types of distances types in predictors
vshort <- ncol(select(trainprotein,contains("vshort")))
short <- ncol(select(trainprotein,contains("short")))
medshort <- ncol(select(trainprotein,contains("medshort")))
medlong <- ncol(select(trainprotein,contains("medlong")))
long <- ncol(select(trainprotein,contains("long")))
vlong <- ncol(select(trainprotein,contains("vlong")))

# Pie Chart
slices <- c(vshort, short, medshort, medlong, long,vlong)
lbls <- c("vshort", "short", "medshort", "medlong", "long", "vlong")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Distribution of distances")


#####################################################

# Methods
(summary(trainprotein)) #Indicates no singularity present (no perfect multicollinearity)
min(trainprotein$angles)
table(sum(is.na(trainprotein))) #No empty values present
full <- lm(accuracy ~., data = trainprotein)
summary(full)
alias(full)
trainprotein$scArgN_bbC_medshort <- NULL
trainprotein$scArgN_bbO_short <- NULL


# VIF Calculations
# Segment dataset into half 
s1 <- lm(trainprotein[,1] ~ ., trainprotein[, 2:301])
s2 <- lm(trainprotein[,1] ~ ., trainprotein[, 302:684])

which.max(vif(s1)) #carbonylC_aliph3HC_medshort 

s2 <- lm(trainprotein[,1] ~ . -carbonylC_aliph3HC_medshort, trainprotein[, 2:301])
vif(s2)
which.max(vif(s2)) #aliph2HC_bbC_medlong

s3 <- lm(trainprotein[,1] ~ . -carbonylC_aliph3HC_medshort-aliph2HC_bbC_medlong, trainprotein[, 2:301])
which.max(vif(s3)) #carboxylC_aliph2HC_medshort

# Continue this process for s1 and s2 and then combine for both s1 and s2 (call it s3) 
# and repeat the VIF process until all predictors have VIF that are less than 10


# Final Predictors with VIF < 10 present in the model. A total of 115 predictors removed
vif_model <- lm(accuracy ~ . -carbonylC_aliph3HC_medshort -aliph2HC_bbC_medlong
               -carboxylC_aliph2HC_medshort -aliph1HC_bbCA_medlong
               -aliph2HC_bbCA_long -carbonylC_aliph1HC_medlong -aliph1HC_bbC_medlong
               -carbonylC_bbCA_medlong -aliph1HC_bbProN_long
               -carbonylC_aliph2HC_medshort -aliph1HC_aliph3HC_medshort
               -carbonylC_carbonylC_medshort -aliph3HC_scArgN_medlong
               -carboxylC_bbCA_medshort -carboxylC_scArgN_vlong
               -carbonylC_bbC_medshort -carboxylC_bbC_long
               -scArgN_bbCA_medlong -scArgN_bbC_medlong
               -aromaticC_scArgN_medlong -scArgN_hydroxylO_medlong
               -scArgN_bbN_medlong -aromaticC_scArgN_medshort
               -scArgN_bbC_short -bbCA_bbC_medlong -aromaticC_bbC_medlong
               -bbN_bbO_short -bbCA_bbCA_medlong -bbC_bbC_medlong
               -sulfur_bbC_medshort -bbC_bbO_medshort -aromaticC_scArgN_vlong
               -scArgN_bbO_medshort -bbN_bbC_medlong -scArgN_bbO_vlong
               -bbCA_bbO_medlong -aromaticC_bbCA_medlong -aliph3HC_bbCA_medlong
               -scArgN_bbCA_long -hydroxylO_bbC_medlong -bbProN_bbCA_medlong
               -bbN_bbN_medlong -scArgN_bbCA_medshort -aliph3HC_bbC_long
               -sulfur_bbCA_medlong -scArgN_bbC_vlong -carboxylO_bbCA_short
               -scAGN_bbCA_medshort -scArgN_bbN_medshort -aromaticC_scArgN_long
               -scLysN_bbC_long -bbN_bbC_short -scArgN_bbO_vshort -carboxylO_bbC_medlong
               -carbonylO_bbC_medlong -hydroxylO_bbCA_medshort -bbC_bbO_long
               -scAGN_bbC_medlong -bbN_bbCA_long -aromaticC_aromaticC_medlong
               -aromaticC_bbC_long -aromaticC_bbC_medshort -sulfur_bbC_long
               -scLysN_bbCA_medlong -carboxylC_hydroxylO_medshort -carbonylC_scAGN_medlong
               -carboxylC_aromaticC_medlong -carbonylC_sulfur_medshort -scAGN_carbonylO_medlong
               -aliph2HC_scArgN_medshort -carboxylC_carboxylO_long -carbonylC_hydroxylO_short
               -scArgN_bbC_long -carbonylC_bbN_long -scArgN_bbN_vlong -scArgN_carboxylO_vlong
               -carboxylO_bbN_long -aliph2HC_scArgN_medlong -carbonylC_aromaticC_medshort
               -carbonylO_bbCA_medlong -aliph3HC_aromaticC_short -carboxylC_bbC_medshort
               -carbonylC_bbO_medlong -hydroxylO_carbonylO_medlong -carbonylC_carbonylO_medlong
               -aliph2HC_bbN_long -aromaticC_aromaticC_medshort -aliph1HC_scAGN_medshort
               -aliph3HC_bbC_medshort -carboxylC_bbO_medlong -scArgN_bbCA_vlong
               -aliph1HC_bbN_long -aliph2HC_scArgN_long -aliph2HC_carboxylO_long
               -bbProN_bbC_long -carboxylC_scLysN_long -carbonylC_carboxylO_short
               -bbN_bbCA_medshort -aromaticC_aromaticC_long -aliph1HC_aliph2HC_medlong
               -carboxylC_carboxylO_medlong -scArgN_bbN_long -scLysN_bbN_medshort
               -aliph2HC_aromaticC_medlong -aliph1HC_aromaticC_short
               -bbO_bbO_medlong -carboxylO_bbCA_long -bbCA_bbC_long
               -aliph1HC_sulfur_medlong -aromaticC_bbCA_medshort -bbProN_bbO_medshort
               -aliph3HC_carbonylO_medlong -aliph3HC_scArgN_vlong -aliph3HC_aromaticC_medshort
               -carbonylC_bbC_long, data = trainprotein)

which.max(vif(vif_model)) #aliph2HC_scLysN_long 
max(vif(vif_model)) #9.880524

# Residual vs Fitted Plot
plot(vif_model$fitted.values, vif_model$residuals, xlab="Fitted Values", ylab="Residuals")

# Residual vs Index Plot
plot(1:nrow(trainprotein), vif_model$residuals, xlab = "Index", ylab = "Residuals")

#QQ Plot of Residuals
qqnorm(vif_model$residuals)
qqline(vif_model$residuals, col="blue", lwd=2)



# Hence we will remove the 115 explanatory variables 
old_dataset <- trainprotein
new_dataset <- select(old_dataset, -c(carbonylC_aliph3HC_medshort ,aliph2HC_bbC_medlong
                              ,carboxylC_aliph2HC_medshort ,aliph1HC_bbCA_medlong
                              ,aliph2HC_bbCA_long ,carbonylC_aliph1HC_medlong ,aliph1HC_bbC_medlong
                              ,carbonylC_bbCA_medlong ,aliph1HC_bbProN_long
                              ,carbonylC_aliph2HC_medshort ,aliph1HC_aliph3HC_medshort
                              ,carbonylC_carbonylC_medshort ,aliph3HC_scArgN_medlong
                              ,carboxylC_bbCA_medshort ,carboxylC_scArgN_vlong
                              ,carbonylC_bbC_medshort ,carboxylC_bbC_long
                              ,scArgN_bbCA_medlong ,scArgN_bbC_medlong
                              ,aromaticC_scArgN_medlong ,scArgN_hydroxylO_medlong
                              ,scArgN_bbN_medlong ,aromaticC_scArgN_medshort
                              ,scArgN_bbC_short ,bbCA_bbC_medlong ,aromaticC_bbC_medlong
                              ,bbN_bbO_short ,bbCA_bbCA_medlong ,bbC_bbC_medlong
                              ,sulfur_bbC_medshort ,bbC_bbO_medshort ,aromaticC_scArgN_vlong
                              ,scArgN_bbO_medshort ,bbN_bbC_medlong ,scArgN_bbO_vlong
                              ,bbCA_bbO_medlong ,aromaticC_bbCA_medlong ,aliph3HC_bbCA_medlong
                              ,scArgN_bbCA_long ,hydroxylO_bbC_medlong ,bbProN_bbCA_medlong
                              ,bbN_bbN_medlong ,scArgN_bbCA_medshort ,aliph3HC_bbC_long
                              ,sulfur_bbCA_medlong ,scArgN_bbC_vlong ,carboxylO_bbCA_short
                              ,scAGN_bbCA_medshort ,scArgN_bbN_medshort ,aromaticC_scArgN_long
                              ,scLysN_bbC_long ,bbN_bbC_short ,scArgN_bbO_vshort ,carboxylO_bbC_medlong
                              ,carbonylO_bbC_medlong ,hydroxylO_bbCA_medshort ,bbC_bbO_long
                              ,scAGN_bbC_medlong ,bbN_bbCA_long ,aromaticC_aromaticC_medlong
                              ,aromaticC_bbC_long ,aromaticC_bbC_medshort ,sulfur_bbC_long
                              ,scLysN_bbCA_medlong ,carboxylC_hydroxylO_medshort ,carbonylC_scAGN_medlong
                              ,carboxylC_aromaticC_medlong ,carbonylC_sulfur_medshort ,scAGN_carbonylO_medlong
                              ,aliph2HC_scArgN_medshort ,carboxylC_carboxylO_long ,carbonylC_hydroxylO_short
                              ,scArgN_bbC_long ,carbonylC_bbN_long ,scArgN_bbN_vlong ,scArgN_carboxylO_vlong
                              ,carboxylO_bbN_long ,aliph2HC_scArgN_medlong ,carbonylC_aromaticC_medshort
                              ,carbonylO_bbCA_medlong ,aliph3HC_aromaticC_short ,carboxylC_bbC_medshort
                              ,carbonylC_bbO_medlong ,hydroxylO_carbonylO_medlong ,carbonylC_carbonylO_medlong
                              ,aliph2HC_bbN_long ,aromaticC_aromaticC_medshort ,aliph1HC_scAGN_medshort
                              ,aliph3HC_bbC_medshort ,carboxylC_bbO_medlong ,scArgN_bbCA_vlong
                              ,aliph1HC_bbN_long ,aliph2HC_scArgN_long ,aliph2HC_carboxylO_long
                              ,bbProN_bbC_long ,carboxylC_scLysN_long ,carbonylC_carboxylO_short
                              ,bbN_bbCA_medshort ,aromaticC_aromaticC_long ,aliph1HC_aliph2HC_medlong
                              ,carboxylC_carboxylO_medlong ,scArgN_bbN_long ,scLysN_bbN_medshort
                              ,aliph2HC_aromaticC_medlong ,aliph1HC_aromaticC_short
                              ,bbO_bbO_medlong ,carboxylO_bbCA_long ,bbCA_bbC_long
                              ,aliph1HC_sulfur_medlong ,aromaticC_bbCA_medshort ,bbProN_bbO_medshort
                              ,aliph3HC_carbonylO_medlong ,aliph3HC_scArgN_vlong ,aliph3HC_aromaticC_medshort
                              ,carbonylC_bbC_long))
dim(new_dataset)
vif(vif_model)

# We will split the new_dataset in 80/20 with 80% allocated to training and rest to validation
N <- nrow(new_dataset)
trainInd <- sample(1:N , round(N *0.8), replace = F)
trainSet <- new_dataset[trainInd ,]
validSet <- new_dataset[ - trainInd ,]

#####################################################

#******* APPENDIX 1 *********************************

# Model number 1 : We will run stepwise forward with BIC with the penalty log(N)

full <- lm(accuracy ~., data = trainSet)
empty <- lm(accuracy ~ 1, data = trainSet)

m_aic <- stepAIC(object= empty, scope = list(upper = full, lower = empty), direction = "forward", k = log(nrow(trainSet)))
summary(m_aic)

model_1 <- lm(formula = accuracy ~ aliph1HC_aliph2HC_long + scLysN_bbC_vlong + 
                aliph2HC_bbN_medshort + aliph1HC_aromaticC_medshort + carbonylC_aromaticC_short + 
                aromaticC_hydroxylO_medlong + aliph1HC_aromaticC_vlong + 
                carboxylC_bbCA_medlong + bbC_bbC_medshort + aromaticC_sulfur_long + 
                bbN_bbCA_medlong + aliph1HC_aromaticC_long + aromaticC_sulfur_short + 
                aliph1HC_aliph1HC_vlong + scAGN_bbN_long + carboxylC_bbN_vlong + 
                aliph1HC_bbC_medshort + sulfur_bbC_medlong + aromaticC_hydroxylO_long + 
                aliph2HC_scLysN_vlong + aliph3HC_bbC_vlong + aliph1HC_bbN_vlong + 
                aliph2HC_aliph2HC_medshort + sulfur_bbC_vlong + scAGN_carbonylO_medshort + 
                aliph1HC_aromaticC_medlong + aliph3HC_bbN_short + aliph1HC_bbO_long + 
                aliph2HC_aromaticC_vlong + aliph1HC_sulfur_short + carboxylC_bbN_long + 
                aliph1HC_bbO_medlong + aliph3HC_aliph3HC_short + carboxylC_aromaticC_long + 
                carboxylC_scLysN_vlong + aromaticC_carbonylO_medlong + carbonylC_bbN_vlong + 
                aliph3HC_bbN_medlong + carbonylC_bbC_medlong + scAGN_bbN_medlong + 
                bbProN_bbCA_medshort + aliph2HC_hydroxylO_long + carboxylO_sulfur_medshort + 
                bbO_bbO_short + aliph2HC_bbN_medlong + bbN_bbC_vlong + bbN_bbO_vlong + 
                bbN_bbN_medshort + carboxylO_carboxylO_vlong + carbonylC_aliph2HC_short + 
                bbO_bbO_long + aliph1HC_bbProN_medlong + aromaticC_bbO_vshort + 
                aromaticC_bbC_short + aromaticC_aromaticC_vlong + aliph1HC_scArgN_long + 
                scArgN_bbO_medlong + carbonylC_sulfur_short + aliph2HC_aromaticC_medshort + 
                carboxylC_hydroxylO_short + aliph2HC_scLysN_medshort + aliph2HC_bbN_vlong + 
                scAGN_bbProN_medshort + hydroxylO_sulfur_short + hydroxylO_sulfur_medshort + 
                scLysN_carboxylO_long + aromaticC_hydroxylO_medshort + aliph3HC_bbN_long + 
                aliph1HC_aliph3HC_long + aliph3HC_hydroxylO_short + sulfur_bbO_vlong + 
                sulfur_bbCA_short + aliph2HC_scLysN_medlong + bbCA_bbO_vshort + 
                carbonylC_hydroxylO_vlong + bbProN_carboxylO_vlong + aliph2HC_hydroxylO_vlong + 
                aromaticC_bbO_vlong + aliph3HC_aromaticC_vlong + aliph1HC_hydroxylO_vlong + 
                aliph1HC_hydroxylO_long + aliph1HC_bbN_medshort + hydroxylO_sulfur_long + 
                hydroxylO_sulfur_vlong + aliph1HC_hydroxylO_medlong + aliph1HC_hydroxylO_medshort + 
                aliph1HC_bbO_medshort + aliph1HC_bbO_short + aliph3HC_bbN_medshort + 
                aliph2HC_scArgN_vlong + carbonylC_scLysN_long + bbC_bbO_medlong + 
                carbonylC_carboxylC_long + bbCA_bbCA_vlong + scLysN_bbN_vlong + 
                aliph1HC_bbProN_vlong + carboxylO_sulfur_medlong + bbO_bbO_vshort + 
                bbC_bbO_short + scAGN_scLysN_medlong, data = trainSet)


AIC(model_1, k = log(nrow(trainSet)))  #3039.501
BIC(model_1) #3039.501 Both models matched    
summary(model_1)$r.squared #0.8918248
summary(model_1)$adj.r.squared #0.8843952
pred9 <- predict(model_1, newdata = validSet)
sqrt(mean((validSet$accuracy - pred9)^2)) # RMSE on validation we get is 0.6171842
sqrt(mean(model_1$residuals^2)) # RMSE on train we get is 0.4849922


#####################################################################

#### ******************* APPENDIX 2 *********************************
# Model number 2: We will run stepwise forward with BIC with the penalty 2log(N)

stepAIC(object = empty, scope = list(upper = full, lower = empty), direction = "forward", k = 2*log(nrow(trainSet)))

model_2 <- lm(formula = accuracy ~ aliph1HC_aliph2HC_long + scLysN_bbC_vlong + 
               aliph2HC_bbN_medshort + aliph1HC_aromaticC_medshort + carbonylC_aromaticC_short + 
               aromaticC_hydroxylO_medlong + aliph1HC_aromaticC_vlong + 
               carboxylC_bbCA_medlong + bbC_bbC_medshort + aromaticC_sulfur_long + 
               bbN_bbCA_medlong + aliph1HC_aromaticC_long + aromaticC_sulfur_short + 
               aliph1HC_aliph1HC_vlong + scAGN_bbN_long + carboxylC_bbN_vlong + 
               aliph1HC_bbC_medshort + sulfur_bbC_medlong + aromaticC_hydroxylO_long + 
               aliph2HC_scLysN_vlong + aliph3HC_bbC_vlong + aliph1HC_bbN_vlong + 
               aliph2HC_aliph2HC_medshort + sulfur_bbC_vlong + scAGN_carbonylO_medshort + 
               aliph1HC_aromaticC_medlong + aliph3HC_bbN_short + aliph1HC_bbO_long + 
               aliph2HC_aromaticC_vlong + aliph1HC_sulfur_short + carboxylC_bbN_long + 
               aliph1HC_bbO_medlong + aliph3HC_aliph3HC_short + carboxylC_aromaticC_long + 
               carboxylC_scLysN_vlong + aromaticC_carbonylO_medlong + carbonylC_bbN_vlong + 
               aliph3HC_bbN_medlong + carbonylC_bbC_medlong + scAGN_bbN_medlong + 
               bbProN_bbCA_medshort + aliph2HC_hydroxylO_long + carboxylO_sulfur_medshort, 
             data = trainSet)

AIC(model_2, k = 2*log(nrow(trainSet)))  #3505.675
BIC(model_2) #3505.675 Both models matched
summary(model_2)$r.squared #0.8372157
summary(model_2)$adj.r.squared #0.8325893
# Calculating the RMSE scores for model 2
pred2 <- predict(model_2, newdata = validSet)
sqrt(mean((validSet$accuracy - pred2)^2)) #RMSE on validation we get is 0.6772338
sqrt(mean(model_2$residuals^2)) #RMSE on train we get is 0.5949453

#########################################################################

#### ******************* APPENDIX 3 *********************************

# Model number 3: We will run stepwise forward with AIC criterion with the penalty 4
stepAIC(object= empty, scope = list(upper = full, lower = empty), direction = "forward", k = 4)

model_3 <- lm(formula = accuracy ~ aliph1HC_aliph2HC_long + scLysN_bbC_vlong + 
     aliph2HC_bbN_medshort + aliph1HC_aromaticC_medshort + aromaticC_hydroxylO_medlong + 
     carbonylC_aromaticC_short + sulfur_bbC_vlong + aromaticC_hydroxylO_long + 
     aliph2HC_scLysN_vlong + bbC_bbC_medshort + aliph3HC_bbC_vlong + 
     aliph1HC_aromaticC_long + carboxylC_bbN_vlong + aliph1HC_aromaticC_vlong + 
     aromaticC_sulfur_long + bbN_bbCA_medlong + aliph1HC_bbC_medshort + 
     sulfur_bbC_medlong + carboxylC_bbN_long + aromaticC_bbCA_vlong + 
     aliph2HC_bbCA_medshort + aliph1HC_bbO_long + aliph1HC_aromaticC_medlong + 
     carboxylC_scLysN_vlong + carboxylO_carboxylO_vlong + aliph2HC_aromaticC_vlong + 
     bbN_bbC_vlong + aliph3HC_bbN_short + aliph1HC_aliph1HC_vlong + 
     scAGN_bbN_long + aliph1HC_aliph3HC_short + aliph1HC_bbO_medlong + 
     aliph3HC_aromaticC_long + bbCA_bbO_vshort + aliph2HC_bbN_vlong + 
     aliph1HC_sulfur_short + aliph1HC_aliph2HC_vlong + carboxylC_aromaticC_long + 
     aliph3HC_hydroxylO_long + carbonylC_bbC_medlong + aliph1HC_aliph3HC_medlong + 
     aliph1HC_hydroxylO_medlong + aliph1HC_bbProN_medlong + aliph1HC_bbN_medshort + 
     aliph1HC_aliph1HC_long + sulfur_bbCA_short + scAGN_carbonylO_medshort + 
     aromaticC_carbonylO_medlong + aliph2HC_aliph3HC_vlong + aliph1HC_scArgN_long + 
     scLysN_bbN_vlong + bbCA_bbCA_vlong + aromaticC_bbO_vlong + 
     aliph2HC_hydroxylO_long + scArgN_hydroxylO_short + aromaticC_scLysN_vlong + 
     aliph3HC_hydroxylO_short + bbO_bbO_long + aliph1HC_hydroxylO_vlong + 
     bbN_bbN_medshort + bbO_bbO_short + carboxylO_bbO_long + aliph3HC_aromaticC_vlong + 
     carbonylC_bbN_vlong + carbonylC_sulfur_short + aliph2HC_aromaticC_medshort + 
     scLysN_bbC_medlong + carbonylC_hydroxylO_vlong + carbonylC_bbProN_medlong + 
     carbonylO_carbonylO_long + carbonylC_bbProN_long + aliph2HC_bbCA_medlong + 
     aliph3HC_aromaticC_medlong + aliph2HC_aliph3HC_short + scAGN_bbC_long + 
     sulfur_bbCA_medshort + carbonylO_bbCA_short + carboxylO_sulfur_medshort + 
     bbN_bbO_vlong + aliph2HC_scArgN_vlong + aliph1HC_hydroxylO_medshort + 
     scLysN_hydroxylO_long + carbonylC_bbProN_vlong + scArgN_bbO_medlong + 
     scLysN_bbN_medlong + aliph2HC_aromaticC_long + scLysN_bbN_long + 
     carbonylC_aliph3HC_medlong + aliph1HC_hydroxylO_long + aliph1HC_bbO_medshort + 
     carbonylO_carbonylO_vlong + aromaticC_bbO_vshort + bbC_bbO_medlong + 
     aliph2HC_bbCA_short + hydroxylO_sulfur_short + carbonylO_bbO_long + 
     aliph3HC_hydroxylO_medshort + scLysN_carbonylO_vlong + hydroxylO_sulfur_long + 
     aliph1HC_sulfur_medshort + hydroxylO_sulfur_medshort + scArgN_carboxylO_long + 
     aliph3HC_bbCA_short + aromaticC_bbO_medlong + carbonylC_scLysN_vlong + 
     scLysN_carboxylO_long + hydroxylO_sulfur_vlong + sulfur_bbO_vlong + 
     aromaticC_bbC_short + aromaticC_aromaticC_vlong + bbN_bbC_vshort + 
     carboxylC_carboxylC_vlong + aliph2HC_scLysN_long + aliph3HC_aliph3HC_short + 
     aliph2HC_sulfur_medshort + carboxylC_bbCA_medlong + bbN_bbO_long + 
     aliph2HC_bbC_long + aliph2HC_bbO_medlong + bbProN_bbO_medlong + 
     scLysN_bbCA_long + carboxylO_sulfur_vlong, data = trainSet)

AIC(model_3, k = 4) #3032.292
summary(model_3)$r.squared #0.9000784
# Calculating the RMSE scores for model 3
pred_3 <- predict(model_3, newdata = validSet)
sqrt(mean((validSet$accuracy - pred_3)^2)) #RMSE on validation we get is 0.6170718
sqrt(mean(model_3$residuals^2)) #RMSE on train we get is 0.4658136

#########################################################################
# Running ICM model:

ICM <- function(trainSet,penalty) {
  pen <- penalty            
  varlist = c()
  varnames = names(trainSet)
  n = nrow(trainSet)
  varorder <- sample(1:ncol(trainSet)) # random order of variables
  minCrit = Inf
  noChange = F
  while (!noChange) {
   noChange = T
    for (i in varorder) { 
     if (i == 1)  #because col 1 is accuracy
       next
     if (i %in% varlist & length(varlist) > 1) {
       index = c(1, varlist[varlist != i]) 
       trainVars = trainSet[, index]
      
       fit = lm(accuracy ~ ., data = trainVars)
      
       if (AIC(fit, k = pen) < minCrit) {
         minCrit = AIC(fit, k = pen)
         varlist = varlist[varlist != i]
         print(paste0("Criterion: ", round(minCrit, 1), ", variables: ", paste0(varnames[varlist], collapse = " ")))
         best.model = fit
         noChange = F
        }
      } else if (!i %in% varlist) {
        index = c(1, varlist, i) 
        trainVars = trainSet[, index]
      
        fit = lm(accuracy ~ ., data = trainVars)
        
        if (AIC(fit, k = pen) < minCrit) {
         minCrit = AIC(fit, k = pen)
         varlist = c(varlist, i)
         print(paste0("Criterion: ", round(minCrit, 1), ", variables: ", paste0(varnames[varlist], collapse = " ")))
         best.model = fit
         noChange = F
        }      
      }
    }
   }
  return (best.model)
  }

#### ******************* APPENDIX 4 *********************************
#ICM for model 1
model1 <- ICM(trainSet,log(nrow(trainSet)))
summary(model1)
AIC(model1, k = log(nrow(trainSet))) #3132.292
BIC(model1) #3132.292 Both models matched
summary(model1)$r.squared #0.9000784
pred_1_ICM <- predict(model1, newdata = validSet)
sqrt(mean((validSet$accuracy - pred_1_ICM)^2)) # RMSE on validation 0.6431905
sqrt(mean(model1$residuals^2)) # RMSE on train 0.498459

#### ******************* APPENDIX 5 *********************************
#ICM for model 2
model2_icm <- ICM(trainSet,2*log(nrow(trainSet)))
summary(model2_icm)
AIC(model2_icm, k = 2*log(nrow(trainSet))) #3132.292
BIC(model2_icm) #3132.292 Both models matched
summary(model2_icm)$r.squared #0.9000784
pred_2_ICM <- predict(model2_icm, newdata = validSet)
sqrt(mean((validSet$accuracy - pred_2_ICM)^2)) # RMSE on validation 0.6878186
sqrt(mean(model2_icm$residuals^2)) # RMSE on train 0.5819979

#### ******************* APPENDIX 6 *********************************
#ICM for model 3
model3_icm <- ICM(trainSet,4)
summary(model3_icm)
AIC(model3_icm, k = 4) #3132.292
summary(model3_icm)$r.squared #0.9000784
pred_3_ICM <- predict(model3_icm, newdata = validSet)
sqrt(mean((validSet$accuracy - pred_3_ICM)^2)) # RMSE on validation 0.6200904
sqrt(mean(model3_icm$residuals^2)) # RMSE on train 0.4438868


##########################################################
# K-fold cross validation

N <- nrow(trainprotein)
K <- 5
validSetSplits <- sample((1:N)%%K + 1)
RMSE1 <- c() #Model_1 stepwise
RMSE2 <- c() #Model_2 stepwise
RMSE3 <- c() #Model_3 stepwise
RMSE4 <- c() #Model_1 ICM
RMSE3 <- c() #Model_2 ICM
RMSE3 <- c() #Model_3 ICM

for (k in 1:K) {
  validSet <- new_dataset[validSetSplits==k,] #Applying the model on the data with removed predictors of VIF
  trainSet <- new_dataset[validSetSplits!=k,]  
  
  full <- lm(accuracy ~ ., data = trainSet)
  empty <- lm(accuracy ~ 1, data = trainSet)
  
  m1 <- stepAIC(object = empty, scope = list(upper = full, lower = empty),
                direction = "forward", k = log(nrow(trainSet)))
  predk1 <- predict(m1, newdata = validSet)
  RMSE1[k] <- sqrt(mean((validSet$accuracy - predk1)^2))  
  
  m2 <- stepAIC(object = empty, scope = list(upper = full, lower = empty),
               direction = "forward", k = 2*log(nrow(trainSet)))
  predk2 <- predict(m2, newdata = validSet)
  RMSE2[k] <- sqrt(mean((validSet$accuracy - predk2)^2))  
  
  m3 <- stepAIC(object = empty, scope = list(upper = full, lower = empty),
                direction = "forward", k = 4)
  predk3 <- predict(m3, newdata = validSet)
  RMSE3[k] <- sqrt(mean((validSet$accuracy - predk3)^2))
  
  m4 <- ICM(trainSet,log(nrow(trainSet)))
  predk4 <- predict(m4, newdata = validSet)
  RMSE4[k] <- sqrt(mean((validSet$accuracy - predk4)^2))
  
  m5 <- ICM(trainSet,2*log(nrow(trainSet)))
  predk5 <- predict(m5, newdata = validSet)
  RMSE5[k] <- sqrt(mean((validSet$accuracy - predk5)^2))
  
  m6 <- ICM(trainSet,4)
  predk6 <- predict(m6, newdata = validSet)
  RMSE5[k] <- sqrt(mean((validSet$accuracy - predk6)^2))
  
}

RMSE1  # 0.5883980 0.5804227 0.6432637 0.6167076 0.6356070
mean(RMSE1)  #0.6128798
RMSE2 # 0.6757525 0.6384283 0.6923043 0.6767486 0.6753434
mean(RMSE2)  #0.6717154
RMSE3 # 0.5839171 0.6130055 0.6234234 0.5907642 0.5867903
mean(RMSE3) #0.60952342
RMSE4 # 0.654191 0.6145002 0.645601 0.589045 0.613904
mean(RMSE4) #0.6213211
RMSE5 # 0.780323 0.7459023 0.760291 0.77398 0.784210
mean(RMSE5) #0.761401
RMSE6 # 0.60792 0.624500 0.68703 0.57345 0.612147
mean(RMSE6) #0.6141432

# It turns out that the model_1 with forward stepwise regression is the best procedure among all 
# the models based on CV prediction error. 

###############################################################

newfull <- lm(accuracy ~., data = trainprotein)
newempty <- lm(accuracy ~ 1, data = trainprotein)

final_model <- stepAIC(object= newempty, scope = list(upper = newfull, lower = newempty), direction = "forward", k = log(nrow(trainSet)))
summary(final_model)

lm_final <- lm(formula = accuracy ~ aliph1HC_aliph2HC_long + scLysN_bbC_vlong + 
                 scArgN_bbN_medlong + aliph2HC_bbN_medshort + aliph1HC_aromaticC_medshort + 
                 sulfur_bbC_vlong + bbC_bbC_medshort + aliph1HC_aromaticC_vlong + 
                 aliph3HC_scArgN_medlong + aromaticC_hydroxylO_medlong + sulfur_bbC_medlong + 
                 aliph2HC_aromaticC_vlong + carboxylC_bbN_vlong + aliph1HC_aromaticC_medlong + 
                 carbonylC_aromaticC_medshort + carboxylC_aromaticC_long + 
                 aliph1HC_aromaticC_long + bbN_bbCA_medlong + aliph1HC_aliph1HC_vlong + 
                 scAGN_bbN_long + aromaticC_hydroxylO_long + aliph1HC_bbO_long + 
                 carboxylC_bbN_long + aliph3HC_aliph3HC_short + aliph3HC_bbN_short + 
                 aliph1HC_bbO_medlong + scLysN_carboxylO_long + aliph1HC_bbProN_medlong + 
                 aromaticC_sulfur_long + scLysN_carboxylO_vlong + scAGN_bbN_medlong + 
                 carbonylC_bbC_medlong + aliph2HC_bbN_medlong + aliph3HC_bbC_vlong + 
                 bbCA_bbO_vshort + bbN_bbC_vlong + aliph3HC_bbN_medlong + 
                 aliph2HC_hydroxylO_long + aromaticC_bbO_vlong + scAGN_carbonylO_medlong + 
                 aliph1HC_aliph3HC_long + aliph1HC_sulfur_short + aliph2HC_aliph3HC_vlong + 
                 carbonylC_aromaticC_short + carbonylC_bbProN_medlong + aliph1HC_hydroxylO_vlong + 
                 aliph2HC_aliph2HC_short + aliph1HC_aromaticC_short + scArgN_bbO_vshort + 
                 carbonylC_sulfur_medshort + carbonylC_bbN_vlong + aliph2HC_bbN_vlong + 
                 bbO_bbO_long + aromaticC_hydroxylO_medshort + aliph2HC_aromaticC_medshort + 
                 bbN_bbO_vlong + aliph2HC_scArgN_vlong + carbonylC_hydroxylO_long + 
                 aliph2HC_scArgN_medlong + scLysN_bbN_medlong + carbonylC_carbonylO_vlong + 
                 carboxylO_bbC_medshort + carboxylC_hydroxylO_short + carboxylO_carboxylO_vlong + 
                 scArgN_bbO_medlong + bbN_bbCA_medshort + aliph3HC_sulfur_long + 
                 aliph3HC_scArgN_long + aliph1HC_scAGN_medshort + carbonylC_bbCA_medlong + 
                 aliph3HC_aromaticC_short + bbN_bbCA_long + aromaticC_bbN_vlong + 
                 hydroxylO_bbC_medlong + bbProN_bbCA_medshort + aliph1HC_bbProN_long + 
                 carboxylC_carboxylC_vlong + aliph3HC_hydroxylO_short + aliph1HC_bbProN_vlong + 
                 sulfur_sulfur_vlong + scAGN_bbN_medshort + scLysN_bbC_medlong + 
                 aliph3HC_hydroxylO_medshort + aromaticC_aromaticC_vlong + 
                 scAGN_scAGN_medshort + scArgN_bbO_medshort + aromaticC_scLysN_vlong + 
                 aliph1HC_hydroxylO_medlong + scArgN_bbCA_medshort + aromaticC_bbO_vshort + 
                 aromaticC_bbC_short + carbonylC_hydroxylO_vlong + bbO_bbO_short + 
                 aliph1HC_aliph3HC_medlong + aliph1HC_aliph1HC_long + aliph1HC_hydroxylO_medshort + 
                 carbonylC_aliph3HC_medshort + scLysN_bbN_vlong + aliph2HC_bbN_long + 
                 bbCA_bbCA_vlong + carbonylC_aliph3HC_medlong + aliph2HC_bbC_medlong + 
                 aliph1HC_bbCA_medlong + aliph3HC_bbCA_medlong + aliph2HC_aliph3HC_long + 
                 aliph1HC_aliph3HC_short + scArgN_bbC_short + aliph1HC_hydroxylO_long + 
                 scArgN_bbCA_long + bbProN_carbonylO_long + scLysN_hydroxylO_long + 
                 carboxylO_bbO_vlong + hydroxylO_sulfur_short + aromaticC_scArgN_vlong + 
                 carboxylC_bbC_medlong + aromaticC_scArgN_medshort, data = trainprotein)

summary(lm_final)
summary(lm_final)$r.squared
length(coef(lm_final))-1

#QQ-Plot of Residual
qqnorm(lm_final$residuals)
qqline(lm_final$residuals, col =" blue ", lwd = 2)

# Residual vs Fitted Plot
plot(lm_final$fitted.values, lm_final$residuals, main= "Final Model Residual vs Fitted Plot", xlab="Fitted Values", ylab="Residuals")

# Residual vs Index Plot
plot(1:nrow(trainprotein), lm_final$residuals, main= "Final Model Residual vs Index Plot",xlab = "Index", ylab = "Residuals")

# Histogram of Residuals
hist(lm_final$residuals, main = "Histogram of Residuals")

# Studentized Residual
plot(fitted(lm_final), rstudent(lm_final), main = "Studentized residuals vs Fitted Values", xlab =" Fitted values ", ylab ="Studentized residuals" )
abline( h = c (3 , -3) , col = "red ", lty =2)
which(abs(studres(lm_final)) > 3)

testprotein <- read.csv("Stat331/Final\ Project/protein-test.csv")

prediction <- predict(lm_final, newdata = testprotein)
writeLines(as.character(prediction), "mypreds.txt")
lm_final$coefficients

final_pred <- predict(lm_final, newdata = testprotein)
df_pred <- data.frame(final_pred)
write.table(df_pred, file = "prediction.txt")
range(df_pred)

library(broom)
library(knitr)
values_produce <- tidy(lm_final)
kable(values_produce)

