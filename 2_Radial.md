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
DLBtemp< -fread("/data/LNG/MR/DLB.tbl", header = T)
HRCtemp <- fread("/data/LNG/MR/chrBpIndexRs.tab", header = F)
names(HRCtemp) <- c("index","SNP","b")
HRC <- subset(HRCtemp, SNP != ".")
#DLBtemp$index <- paste(DLBtemp$chr,DLBtemp$bp,sep = "_")
DLBtemp$index <- DLBtemp$MarkerName
DLB <- merge(DLBtemp, HRC, by = "index")
DLB$effect_allele <- toupper(as.character(DLB$Allele1))
DLB$other_allele <- toupper(as.character(DLB$Allele2))
Out_data <- format_data(DLB, type="outcome", beta_col = "Effect", se_col = "StdErr", eaf_col = "Freq1", pval_col = "P-value")
```

### Run MR analysis
```
#token <- get_mrbase_access_token()
token <- "ya29.Glt2BSKxkLWJCzEV6_L_QL9fEfQnMdHDnAGeOy-4wl-6sSFz2HuyV2azTIlI73Wa9Qiehd2b7UhjyFDi1BcP-Uwd2xj3OnJAQKp-xZnOHvySYSjisnLdhOOVYaye"
possibleInstruments <- available_outcomes(access_token = token)
listOfGwasIds <- read.table("/data/LNG/RESOLUTION/MR/TRAITS.txt", header = T)
for(i in 1:401)
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
   r_input <- format_radial(BXG = dat$beta.exposure, BYG = dat$beta.outcome, seBXG = dat$se.exposure,seBYG =dat$se.outcome, RSID = dat$SNP)
   ivw <- ivw_radial(r_input,0.05,2)
   egger <- egger_radial(r_input,0.05,2)
   plot_radial(c(ivw,egger),TRUE,TRUE,FALSE)
   funnel_radial(c(ivw,egger),TRUE)
   }
else
{
print("FAIL")
}
#q(no)
