### Mendelian Randomization 

	- **Author(s):** Sara Bandres-Ciga
	- **Date Last Updated:** 9.05.2019

### R code
```
library(data.table)
library(TwoSampleMR)
install.packages("MRinstruments")
library(MRInstruments)
install.packages("WSpiller/RadialMR")
library("RadialMR")
```

### Format SumStats

```
ALStemp <- fread("~/Desktop/mendelRandoAuto/SummaryStats/ALS/alsMetaSummaryStats_march21st2018.tab", header = T)
HRCtemp <- fread("~/Desktop/mendelRandoAuto/SummaryStats/PD/chrPosRs.tab", header = F)
names(HRCtemp) <- c("index","SNP")
HRC <- subset(HRCtemp, SNP != ".")
ALStemp$index <- paste(ALStemp$CHR,ALStemp$BP,sep = "_")
ALS <- merge(ALStemp, HRC, by = "index")
ALS$effect_allele <- toupper(as.character(ALS$Allele1))
ALS$other_allele <- toupper(as.character(ALS$Allele2))
```

### Run Leave-one-out analyses
```
Out_data <- format_data(ALS, type="outcome", beta_col = "Effect", se_col = "StdErr", eaf_col = "effectAlleleFreq", pval_col = "P", snp_col = "SNP.y")
#token <- get_mrbase_access_token()
token <- "ya29.Glt2BSKxkLWJCzEV6_L_QL9fEfQnMdHDnAGeOy-4wl-6sSFz2HuyV2azTIlI73Wa9Qiehd2b7UhjyFDi1BcP-Uwd2xj3OnJAQKp-xZnOHvySYSjisnLdhOOVYaye"
possibleInstruments <- available_outcomes(access_token = token)
#### note, possible values of i range from 1 to 1688 ####
listOfGwasIds <- read.table("~/Desktop/listofexposures.txt", header = T)
for(i in 1)
{
  instrumentId <- as.character(listOfGwasIds$id[i])
  tag <- paste("INSTRUMENT IS ",instrumentId," AT i = ",i, sep = "")
  print(tag)
  flagged <- "nope"
  Exp_data <- extract_instruments(outcomes=instrumentId, p1 = 5e-08, clump = TRUE, p2 = 5e-08,
                                  r2 = 0.001, kb = 10000, access_token = token,
                                  force_server = TRUE)
  skip <- ifelse(length(Exp_data$beta.exposure) < 1, 1, 0)
  if(skip == 0)
  {
    dat <- harmonise_data(exposure_dat=Exp_data, outcome_dat=Out_data, action=2)
    res_loo <- mr_leaveoneout(dat)
    res_single <- mr_singlesnp(dat)
    p2 <- mr_forest_plot(res_loo)
    p2
  }
  else
  {
    print("FAIL")
  }
}
#q(no)
```



