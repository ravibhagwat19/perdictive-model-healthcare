
## Step 1  - Read in Data
data=read.csv("casestudydata.csv")
names(data)


#separating the data by CKD values
class(data)
summary(data)
out_sample=which(is.na(data$CKD)==1)
data_out=data[out_sample,]   ## the ones without a disease status
data_in=data[-out_sample,]   ## the ones with a disease status

#finding out how much data is missing for each variable
sum(is.na(data_in))
sum(complete.cases(data_in))
install.packages("mice")
library(mice)
misspattern = md.pattern(data_in)
write.csv(misspattern, "missingdata_pattern.csv")

## Step 2  - Missing Data
summary(data_in)
dim(data_in)
?na.omit
data_in=na.omit(data_in)
dim(data_in)

#separating female and ont female
out_sample1 = which((data$Female)==1)
data_female=data[-out_sample1,] #observations are not female
data_nfem = data[out_sample1,] # observations are female
summary(data_female)
summary(data_nfem)

## Step 3 and 4 - Correlation
cor(data_in)
summary(data_in)
data_new=model.matrix(~-1+Racegrp+CareSource,data=data_in)
summary(data_new)
data_in=data_in[,-c(4,8)]
data_in=cbind(data_in,data_new)
cor(data_in)
names(data_in)
data_in=data_in[,-33]
cor(data_in)

#removing the first column 'ID'
rownames(data_in)=data_in[,1]
data_in=data_in[,-1]
cor(data_in)

#Imputing data using mice
library("mice")
m_data = mice(data_imp, method = "pmm", m=5)
data_imp = complete(m_data)
write.csv(data_imp, file = "important variables.csv")


#Seeding data using mice
install.packages("missForest")
library("missForest")
data.data_in <- prodNA(data_in, noNA = 0.1)
summary(data.data_in)

#remove categorical variables
data.data_in <- subset(data.data_in, select = -c(CKD))
summary(data.data_in)
md.pattern(data.data_in)

## any highly correlated - large relation to PCA


## Step 5 - Run a Regression
model=lm(CKD~.,data=data_in)
summary(model)

summary(predict(model))


## Step 6 - Screening tool

#  GET CREATIVE :))))  :)) ha ha 

## Hint:  Do you have Diabetes?  (Yes = c points or No = d points)

## PCA intro

#How do we plot all 33 variables?

dim(data_in)
install.packages('pca3d')
library(pca3d)

pca <- prcomp( data_in[,-30], scale.= TRUE )
pca3d( pca, group= data_in[,30] )

