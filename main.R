# Author: Chan Min Min
# Date Created: 31 Dec 2021

library(tableone)
library(MatchIt)
library(cobalt)
library(ggplot2)
library(dplyr)
library(survival)
library(survminer)
library(mice)

# Flags
workdir = "C:/Users/idscmm/Downloads/RWDDrugEfficacy/"
fmice = 1
fall = 1


# Load Event Data
Event_duration <- read.csv(paste(workdir,"/Biostats - RWD drug efficacy/Event_duration.csv",sep=""))
colnames(Event_duration) <- c('patient_id', 'label', 'time', 'treatment_variable')
Event_duration$time <- as.integer(Event_duration$time * 12)

# Plot Kaplan Meier of original data before imputation
Surv_unmatched <- Surv(time = Event_duration$time, event = Event_duration$label)
KM_unmatched <- survfit(Surv_unmatched ~ treatment_variable, data = Event_duration, type="kaplan-meier")
ggsurvplot(KM_unmatched, data = Event_duration, pval = TRUE, xlab = "Time to Bleeding (Months)",palette = c("red", "blue"),
           legend.labs = c("Drug_A", "Drug_B"), legend.title="")
Event_duration <-dplyr::select(Event_duration, -c('treatment_variable'))

# Load Processed Data processed in python
if (fmice == 1) {
  data_full <- read.csv(paste(workdir,"/Processed_Data/data_full_unscaled.csv",sep=""))
} else {
  if (fall == 1) {
    data_full <- read.csv(paste(workdir,"/Processed_Data/data_knn_allvariables.csv",sep=""))
  } else {
    data_full <- read.csv(paste(workdir,"/Processed_Data/data_knn_4variables.csv",sep=""))
  }

}

data_full<- merge(x = data_full, y = Event_duration, by = "patient_id")

#Rearrange Columns
data_full <- data_full[, c("patient_id","label","time","treatment_variable","sex","age","other_drugs_1","other_drugs_2",
                           "other_drugs_3", "other_drugs_4", "other_drugs_5", "other_drugs_6", "other_drugs_7", "other_drugs_8",
                           "diagnosis_1", "diagnosis_2", "diagnosis_3", "diagnosis_4", "diagnosis_5", "diagnosis_6",
                           "diagnosis_7", "diagnosis_8", "diagnosis_9", "diagnosis_10", "diagnosis_11", "diagnosis_12",
                           "diagnosis_13", "diagnosis_14", "diagnosis_15",
                           "lab_1", "lab_2","lab_3","lab_4","lab_5","lab_6", "lab_7", "lab_8",
                           "Diag_Score_1_0", "Diag_Score_1_1", "Diag_Score_1_2", "Diag_Score_1_3", "Diag_Score_1_4","Diag_Score_1_5",
                           "Diag_Score_2_0", "Diag_Score_2_1", "Diag_Score_2_2", "Diag_Score_2_3", "Diag_Score_2_4","Diag_Score_2_5", "Diag_Score_2_6", "Diag_Score_2_7")] # leave the row index blank to keep all rows

colmissing <- apply(data_full, 2,
                    function(x){ sum(is.na(x)) })
colmissing


if (fmice == 1) {
  if (fall == 1) {
    data.complete <- data_full
    imputed_Data <- mice(data.complete, m=5, maxit = 50, method = 'pmm', seed = 500)
  } else { # Drop Variables with more than 60% missingness before doing multiple imputation
    data.complete <-dplyr::select(data_full, -c('lab_2', 'lab_4', 'lab_3'))
    imputed_Data <- mice(data.complete, m=5, maxit = 50, method = 'pmm', seed = 500)
    saveRDS(imputed_Data, file = "C:/Users/idscmm/Downloads/RWDDrugEfficacy/Processed_Data/imputed_Data_mice_4variables.rds")    
  }
  data.complete <- complete(imputed_Data,5)
} else {
  if (fall == 1) {
    data.complete <- data_full
  } else {
    data.complete <-dplyr::select(data_full, -c('lab_2', 'lab_4', 'lab_3'))
  }

}

data.complete <- na.omit(data.complete)

# Prediction Performance after imputation
smp_size <- floor(0.8 * nrow(data.complete)) # 80% of the sample size
set.seed(123) # set the seed to make your partition reproducible
train_ind <- sample(seq_len(nrow(data.complete)), size = smp_size)
train <- data.complete[train_ind, ]
test <- data.complete[-train_ind, ]

# COX PH
fit_concordance<- coxph(Surv(time,label) ~ ., data=train[,-c(1)])
concordance(fit_concordance, timewt="n",newdata=test[,-c(1)]) 
# logistic regression
fit_auc <- glm(label ~ ., binomial, data= train[,-c(1,3)])
concordance(fit_auc, newdata=test[,-c(1,3)])  # equal to the AUC

# Check Common Support
psmodel <- glm(treatment_variable~., family=binomial(), data=data.complete[,-c(1,2,3)])
summary(psmodel)
pscore <- psmodel$fitted.values
data.complete$fitted_pscore <- pscore

# Plot Histogram of PS before trimming
data.complete %>%
  ggplot( aes(x=fitted_pscore, fill=treatment_variable,group=treatment_variable)) +
  geom_histogram(aes(y=-1*..density..),alpha=0.6,colour='black',
               data = ~ subset(., treatment_variable %in% c(1)))+
  geom_histogram(aes(y=..density..),alpha=0.6,colour='black',
               data = ~ subset(., treatment_variable %in% c(0)))+
  ylab('density') +
  xlab('propensity score')
quantile_99 <- quantile(pscore, c(.1, .9)) 


data.complete <- data.complete[which(data.complete$fitted_pscore > quantile_99[1] & data.complete$fitted_pscore < quantile_99[2]),]

# Plot Histogram of PS after trimming
data.complete %>%
  ggplot( aes(x=fitted_pscore, fill=treatment_variable,group=treatment_variable)) +
  geom_histogram(aes(y=-1*..density..),alpha=0.6,colour='black',
               data = ~ subset(., treatment_variable %in% c(1)))+
  geom_histogram(aes(y=..density..),alpha=0.6, colour = 'black',
               data = ~ subset(., treatment_variable %in% c(0)))+
  ylab('density') +
  xlab('propensity score')


m.out_final <- matchit(treatment_variable~.,
                 data=data.complete[,-c(1,2,3,dim(data.complete)[2])],
                 distance=log(data.complete$fitted_pscore/(1-data.complete$fitted_pscore)),
                 method='nearest',
                 replace=FALSE,
                 caliper = 0.1,
                 ratio=1)

bal.tab(m.out_final, m.threshold = 0.1, un = TRUE)
bal.tab(m.out_final, v.threshold = 2)
love.plot(bal.tab(m.out_final, m.threshold=0.1),
         stat = "mean.diffs", 
         grid=TRUE,
         stars="raw",
         abs = F)


# ==== Kaplan Meir Plots (Before and After Trimming) ===== #
m.out <- m.out_final
DtMatchedFinal <- data.complete[row.names(match.data(m.out)),]
Surv_unmatched <- Surv(time = data.complete$time, event = data.complete$label)
KM_unmatched <- survfit(Surv_unmatched ~ treatment_variable, data = data.complete, type="kaplan-meier")
ggsurvplot(KM_unmatched, data =data.complete, pval = TRUE, xlab="Time to Bleeding (Months)", legend.title="",
           legend.labs = c("Drug_B", "Drug_A"), palette = c("blue","red"))
 
Surv_matched <- Surv(time = DtMatchedFinal$time, event = DtMatchedFinal$label)
KM_matched <- survfit(Surv_matched ~ treatment_variable, data = DtMatchedFinal, type="kaplan-meier")
ggsurvplot(KM_matched, data = DtMatchedFinal, pval = TRUE, xlab="Time to Bleeding (Months)", legend.title="",
           legend.labs = c("Drug_B", "Drug_A"), palette = c("blue","red"))

# Absolute Average Survival Time Improvement
abs_surv = mean(KM_matched[1]$surv) - mean(KM_matched[2]$surv)
abs_surv 

#===== Sensitivity Analysis ===== #
library(rbounds)
m.pairs <- cbind(data.complete[row.names(m.out$match.matrix), 'label'], data.complete[m.out$match.matrix, 'label'])
m.pairs <- m.pairs[complete.cases(m.pairs), ]
x <- sum((m.pairs[,1]==FALSE) & (m.pairs[,2]==TRUE))
y <- sum((m.pairs[,1]==TRUE) & (m.pairs[,2]==FALSE))
psens(x=m.pairs[,1], y=m.pairs[,2], Gamma = 1.5, GammaInc 
      = 0.1)

#===== Hazard Ratio from CoxPH ===== #
fit_matched <- coxph(Surv(time, label) ~ treatment_variable ,data= DtMatchedFinal[,-c(1,dim(DtMatchedFinal)[2])])
summary(fit_matched)