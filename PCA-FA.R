#PCA & FA

## Step 1  - Read in Data
data=read.csv("casestudydata_pca_fa.csv")
names(data)

rownames(data)=data[,1]
data=data[,-1]

#separating the data by CKD values
class(data)
summary(data)
out_sample=which(is.na(data$CKD)==1)
data_out=data[out_sample,]   ## the ones without a disease status
data_in=data[-out_sample,]   ## the ones with a disease status

#Imputing data using mice
library("mice")
m_data = mice(data_in, method = "pmm", m=5)
data_in = complete(m_data)

data_in=scale(data_in)
summary(data_in)

var1=var(data_in$Age)
var2=var(data_in$Female)
var3=var(data_in$Educ)
var4=var(data_in$Unmarried)
var5=var(data_in$Income)
var6=var(data_in$Insured)
var7=var(data_in$Weight)
var8=var(data_in$Height)
var9=var(data_in$BMI)
var10=var(data_in$Obese)
var11=var(data_in$Waist)
var12=var(data_in$SBP)
var13=var(data_in$DBP)
var14=var(data_in$HDL)
var15=var(data_in$LDL)
var16=var(data_in$Total.Chol)
var17=var(data_in$Dyslipidemia)
var18=var(data_in$PVD)
var19=var(data_in$Activity)
var20=var(data_in$PoorVision)
var21=var(data_in$Smoker)
var22=var(data_in$Hypertension)
var23=var(data_in$Fam.Hypertension)
var24=var(data_in$Diabetes)
var25=var(data_in$Fam.Diabetes)
var26=var(data_in$Stroke)
var27=var(data_in$CVD)
var28=var(data_in$Fam.CVD)
var29=var(data_in$CHF)
var30=var(data_in$Anemia)
var31 = var(data_in$CKD)


Total_var = var1 + var2 + var3 + var4 + var5 + var6 + var7 + var8 + var9 + var10 + var11 + var12+ var13 + var14 + var15 + var16 + var17 + var18 + var19 + var20 + var21+ var22 + var23 + var24 +var25 + var26 + var27 + var28 + var29 + var30
Total_var

summary(data_in)
rho=cor(data_in)
rho

sigma_ckd=var(data_in)
vars_ckd=diag(sigma_ckd)
round(vars_ckd,2)
percentvars_ckd=vars_ckd/sum(vars_ckd)
percentvars_ckd

## Step 3 - Compute all the eigenvalues and eigenvectors in R

eigenvalues=eigen(sigma_ckd)$values
eigenvectors=eigen(sigma_ckd)$vectors

round(eigenvalues,2)
round(eigenvectors,2)

# define principal componenets
y1=as.matrix(data_in)%*%(eigenvectors[,1])
y2=as.matrix(data_in)%*%(eigenvectors[,2])
y3=as.matrix(data_in)%*%(eigenvectors[,3])
y4=as.matrix(data_in)%*%(eigenvectors[,4])
y5=as.matrix(data_in)%*%(eigenvectors[,5])
y6=as.matrix(data_in)%*%(eigenvectors[,6])
y7=as.matrix(data_in)%*%(eigenvectors[,7])
y8=as.matrix(data_in)%*%(eigenvectors[,8])


y_aud=as.matrix(data_in)%*%eigenvectors
y_aud


## Step 4 - Check variance estimates of the pcs and all other properties

var1 + var2 + var3 + var4 + var5 + var6 + var7 + var8 + var9 + var10 + var11 + var12+ var13 + var14 + var15 + var16 + var17 + var18 + var19 + var20 + var21+ var22 + var23 + var24 +var25 + var26 + var27 + var28 + var29 + var30 + var31
percentvars_aud

percentvars_pc=eigenvalues/sum(eigenvalues)
round(percentvars_pc,6)


## C: sum of variances  
var1 + var2 + var3 + var4 + var5 + var6 + var7 + var8 + var9 + var10 + var11 + var12+ var13 + var14 + var15 + var16 + var17 + var18 + var19 + var20 + var21+ var22 + var23 + var24 +var25 + var26 + var27 + var28 + var29 + var30 + var31
var(y1)+var(y2)+var(y3)

## D: Magnitude of eigenvectors are importance of kth variable in the ith PC
eigenvectors
eigenvalues

## E:  correlation between Yi and Xk
eigenvectors[,1]
cor(y1,data_in)
cor(y2,data_in)
cor(y3,data_in)
eigenvectors[,1]*sqrt(eigenvalues[1])/sqrt(diag(vars_aud))


# E: are they uncorrelated?
sigma_aud=var(y_aud)
rho_aud=cor(y_aud)
round(rho_aud, 4)

## Step 5 - Nice plots and Interpret PCs

ts.plot(cbind(percentvars_aud,percentvars_pc),col=c("blue","red"),xlab="ith vector",ylab="percent variance")

plot(y1,y2)
text(y1,y2, labels = rownames(data_in), pos = 4)

#prcomp(data)
#autoplot(prcomp(data))

## Step 6 - regression, use as an input
set.seed(1000)
dv=rowSums(data_in)+rnorm(100,mean=0,sd=10)
summary(lm(dv~as.matrix(data_in)))

summary(lm(dv~as.matrix(y_aud)))

## let's pick the best two
cor(dv,data_in)
summary(lm(dv~data_in$LDL+data_in$Diabetes+data_in$CVD))
summary(lm(dv~y1+y2))

## Step 7 - standardize variables, see how PC changes.  "each is equally important to diet"
data_aud=scale(data_aud)
data_aud=as.data.frame(data_aud)

################# PCA done ##########################

################# FA Starts ########################

## Step 0 - Read in Data

#data_in=scale(data_in)

vartotal=var(data_in)
var_diagonal=diag(vartotal)
var_diagonal

percentmast=var_diagonal/sum(var_diagonal)
percentmast

rho_mast= cor(data_in)
rho_mast


## Step 2 - Compute the eigenvalues and eigenvactors of the correlation matrix
eigenvalues=eigen(rho_mast)$values
round(eigenvalues,2)
eigenvalues>2
m=4 ## can change this to include more factors


eigenvectors=eigen(rho_mast)$vectors
eigenvectors

## Step 3 - Compute Estimated Factor Loadings
L=matrix(nrow=31,ncol=m)
for (j in 1:m){
  L[,j]=sqrt(eigenvalues[j])*eigenvectors[,j]  
}

L  # first column is Factor 1, 2nd column is factor 2


## Step 4 - Compute common variance and unique variance

common=rowSums(L^2)
unique=1-common  ## this diagonal of error matrix

common
unique


## Step 5 - Check the model to reproduce correlation

phi=diag(31)*unique

recreate=L%*%t(L)+phi
recreate

rho_mast

## Step 6 - Create Residual Matrix

residual=rho_mast - recreate
residual  ## check to see if off-diagonal elements are "small"

sum(residual[lower.tri(residual)]^2)  ## sum of square of off-diagonal elements

#sum(sqrt(residual[lower.tri(residual)]^2))/(136-2)  ## sum of square of off-diagonal elements
sum(eigenvalues[18:length(eigenvalues)]^2) 

sum(eigenvalues[3:5]^2)  ## sum of square of non-used eigenvalues



## Step 7  - Plot pairs of loadings to interpret factor loadings
## if we can't tell, we may need to do a varimax rotation

plot(L[,1],L[,2],col=1:5,xlab="Loading 1",ylab="Loading 2")
text(L[,1],L[,2],colnames(data_in))


## Step 8

install.packages('psych')
library(psych)

## should reproduce our results
fit2 <- principal(data_in, nfactors=4, rotate="none")

fit2

fit <- principal(data_in, nfactors=31, rotate="varimax")
fit
plot(fit$loadings[,1],fit$loadings[,2],col=1:31)
text(fit$loadings[,1],fit$loadings[,2],colnames(data_in))

