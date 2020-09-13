## Mendelian Randomization 
	- **Author(s):** Sara Bandres-Ciga
	- **Date Last Updated:** 9.05.2019
	
### R code
```
library(TwoSampleMR)
install.packages("glmnet")
library(glmnet)
devtools::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
```

### Extract exposures
```
token <- "ya29.Glt2BSKxkLWJCzEV6_L_QL9fEfQnMdHDnAGeOy-4wl-6sSFz2HuyV2azTIlI73Wa9Qiehd2b7UhjyFDi1BcP-Uwd2xj3OnJAQKp-xZnOHvySYSjisnLdhOOVYaye"
exposure_dat <- mv_extract_exposures(c(300,301,783), clump_r2 = 0.001, clump_kb = 10000, harmonise_strictness = 2, access_token = token)
exposure_dat <- mv_extract_exposures(c(299,300,302), clump_r2 = 0.001, clump_kb = 10000, harmonise_strictness = 2, access_token = token)
#exposure_dat <- mv_extract_exposures(c(300,301,843,928,927,880,882,926,885,884,895,883,783,881,931,870,933,908,904,867,907,924,905,868,869,906,953,955,925,640,909,929,930,871), clump_r2 = 0.001, clump_kb = 10000, harmonise_strictness = 2, access_token = token) 
exposure_dat <- mv_extract_exposures(c(843,928,927,880,882,926,885,884,895,883,881,931,870,933,908,904,867,907,924,905,868,869,906,953,955,925,640,909,929,930,871), clump_r2 = 0.001, clump_kb = 10000, harmonise_strictness = 2, access_token = token) 
mvdat <- mv_harmonise_data(exposure_dat, Out_data)
```

### Do Lasso regression
```
a <- glmnet(x=mvdat$exposure_beta, y=mvdat$outcome_beta, weight=1/mvdat$outcome_se^2, intercept=0)
```

### PLot the coefficients over range of shrinkage parameters
```
plot(a)
```

### Run cross validation and choose lowest MSE
```
b <- cv.glmnet(x=mvdat$exposure_beta, y=mvdat$outcome_beta, weight=1/mvdat$outcome_se^2, intercept=0)
plot(b)
```
### What traits remain?
```
coef(b, s = "lambda.min"
```
