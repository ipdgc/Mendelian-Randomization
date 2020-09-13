## Expression and methylation QTL Mendelian Randomization
```

The SMR software tool implements the SMR & HEIDI methods to test for pleiotropic association 
between the expression level of a gene and a complex trait of interest using summary-level data 
from GWAS and expression quantitative trait loci (eQTL) studies (Zhu et al. 2016 Nature Genetics). 
The methodology can be interpreted as an analysis to test if the effect size of a SNP on the 
phenotype is mediated by gene expression. This tool can therefore be used to prioritize genes 
underlying GWAS hits for follow-up functional studies.
```

## Convert gene names to ENSEMBL ID format.
```
library(data.table)
geneList <- fread("MR_list.txt", header = F)
names(geneList) <- "GENE"
ensg <- fread("genIdsFromHugo.txt", header = T, sep = "\t")
data <- merge(geneList, ensg, by.x = "GENE", by.y = "Approved Symbol")
names(data)[12] <- "ensgOut"
write.table(data$ensgOut, "ensembl_MR_list.txt", quote = F, sep = "\t", row.names = F, col.names = F)
```


## Run Summary-based Mendelian Randomization analyses
``` 
smr_Linux --bfile workingData  --smr-multi --gwas-summary toSMR_summarystats.txt --beqtl-summary eQTLGen/cis-eQTLsFDR0.05-ProbeLevel.txt_besd --out eQTLGen_multi_PD_PATHWAYS --genes ensembl_MR_list.txt --thread-num 12
smr_Linux --bfile workingData  --smr-multi --gwas-summary toSMR_summarystats.txt --beqtl-summary Brain-eMeta/Brain-eMeta --out Brain-eMeta_multi_PD_PATHWAYS --genes MR_list.txt --thread-num 12
smr_Linux --bfile workingData  --smr-multi --gwas-summary toSMR_summarystats.txt --beqtl-summary /Brain-mMeta/Brain-mMeta --out Brain-mMeta_multi_PD_PATHWAYS --genes MR_list.txt --thread-num 12
smr_Linux --bfile workingData  --smr-multi --gwas-summary toSMR_summarystats.txt --beqtl-summary GTEx_V7_cis_eqtl_summary_lite/Brain_Substantia_nigra_1e-05 --out GTEx-Brain_Substantia_nigra_1e-05_PD_PATHWAYS --genes MR_list.txt --thread-num 12
```

## NOTES
```
# MR_list.txt contains gene names
# genIdsFromHugo.txt is used to convert them into ESN IDs. The list genIdsFromHugo.txt can be found as part of the resources.
# --bfile reads individual-level SNP genotype data (workingData PLINK binary format) from a reference sample for LD estimation, i.e. .bed, .bim, and .fam files.
# --gwas-summary reads summary-level data from GWAS. The input format follows that for GCTA-COJO analysis (http://cnsgenomics.com/software/gcta/#COJO).
# --beqtl-summary reads summary-level data from a eQTL study in binary format. eQTL summary data is stored in three separate files .esi (SNP information, in the same format as the PLINK .bim file), .epi (probe information) and .besd (eQTL summary statistics in binary format). See Data Management from https://cnsgenomics.com/software/smr/#SMR&HEIDIanalysis for more information. Data from the Westra study (Westra et al. 2013 Nat Genet), the Genotype‐Tissue Expression (GTEx) Consortium v6, the CommonMind Consortium (CMC; dorsolateral prefrontal cortex), the Religious Orders Study and Memory and Aging Project (ROSMAP) and the Brain eQTL Almanac project (Braineac; 10 brain regions)), are available to download in https://cnsgenomics.com/software/smr/#SMR&HEIDIanalysis
# --thread-num specifies the number of OpenMP threads for parallel computing. The default value is 1.
# --smr-multi turns on set-based SMR test in the cis-region.
# --genes specifies your list of genes
# --out refers to your output name
```
