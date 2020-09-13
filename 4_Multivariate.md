## Mendelian Randomization 
	- **Author(s):** Sara Bandres-Ciga
	- **Date Last Updated:** 9.05.2019

### R code ###
```
library(data.table)
library(TwoSampleMR)
```

### Format PD sumStats ###
```
ALStemp <- fread("~/Desktop/mendelRandoAuto/SummaryStats/ALS/alsMetaSummaryStats_march21st2018.tab", header = T)
HRCtemp <- fread("~/Desktop/mendelRandoAuto/SummaryStats/PD/chrPosRs.tab", header = F)
names(HRCtemp) <- c("index","SNP")
HRC <- subset(HRCtemp, SNP != ".")
ALStemp$index <- paste(ALStemp$CHR,ALStemp$BP,sep = "_")
ALS <- merge(ALStemp, HRC, by = "index")
ALS$effect_allele <- toupper(as.character(ALS$Allele1))
ALS$other_allele <- toupper(as.character(ALS$Allele2))
Out_data <- format_data(ALS, type="outcome", beta_col = "Effect", se_col = "StdErr", eaf_col = "effectAlleleFreq", pval_col = "P", snp_col = "SNP.y")
#token <- get_mrbase_access_token()
token <- "ya29.GluLBbJT6HOco3fE_EHB9Gm4YIv1oQcdAITj3aEI8UGNBF2s0J61VY6PkHJGOrPNXJ__c-ixt_I4fcz4fLC_ikOm2E8ZfWipF7cjFdwbLO8JA-XdilAsM8-4y3a4"
possibleInstruments <- available_outcomes(access_token = token)
```
### Run Multivariate analysis
```
id_exposure <- c (300,"UKB-a:108")
exposure_dat <- mv_extract_exposures(id_exposure, clump_r2= 0.001, clump_kb = 10000, harmonise_strictness = 2, access_token= token)
mvdat <- mv_harmonise_data(exposure_dat, Out_data)
#write.table(mvdat, file = "Multivariable.table", na = "NA", quote = F, row.names = F, sep = "\t")
res <- mv_multiple(mvdat)
write.table(res, file ="Results.Multivariable.table", na = "NA", quote = F, row.names = F, sep = "\t")
#q(no)
```
