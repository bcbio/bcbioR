# Overview

This is an example of complex DE analysis with multiple covariates with multiple levels.

We have the SEX variable (2 levels) and the GENOTYPE VARIABLE (4 levels)

# Intercept Analysis

```
# Model design and creating dds object from the dataset
design = ~sex + genotype + sex:genotype
dds <- DESeqDataSet(se_Striatum, design)
```

## Filtering lowly expressed genes
We are filtering out genes with fewer than 10 raw counts in total and are present in fewer than 3 samples.

```
keep <- rowSums(counts(dds)>=10) >=4
dds <- dds[keep, ]
#dds # comment out this line to print the dds object and compare the dimension of the dataset before and after filtering is applied.
```

setting up WT as reference genotype and Male and reference sex. Otherwise DESeq2 will use the conditions in their alphabetical order.

```
dds$genotype <- relevel(dds$genotype, ref = "WT")
dds$sex <- relevel(dds$sex, ref = "Male")

#Checking model design and reference condition comment out the three lines below to print the design and order of genotype and sex
design(dds)
levels(dds$genotype)
levels(dds$sex)

#estimating size factors for normalization and fitting our model with DESeq model
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
resultsNames(dds) #This will print out the name of coefficients being compared, comment it out to view

# get coefficient matrix
mod_mat <- model.matrix(design(dds), data = colData(dds))
mod_mat

(Intercept) sexFemale   genotypeCR3KO   genotypeQ175    genotypeQ175_CR3KO  sexFemale:genotypeCR3KO sexFemale:genotypeQ175  sexFemale:genotypeQ175_CR3KO
a10_st_q175_m_r1    1   0   0   1   0   0   0   0
a12_st_q175_f_r1    1   1   0   1   0   0   1   0
a14_st_wt_m_r1  1   0   0   0   0   0   0   0
a16_st_wt_f_r1  1   1   0   0   0   0   0   0
```


coefficient weights extracted from the mod_mat above

```
WT_M <- c(1, 0, 0, 0, 0, 0, 0, 0)
WT_F <- c(1, 1, 0, 0, 0, 0, 0, 0)
WTCR3ko_M <- c(1, 0, 1, 0, 0, 0, 0, 0)
WTCR3ko_F <- c( 1, 1, 1, 0, 0, 1, 0, 0)
Q175_M <- c(1, 0, 0, 1, 0, 0, 0, 0)
Q175_F <- c(1, 1, 0, 1, 0, 0, 1, 0)
Q175CR3ko_M <- c(1, 0, 0, 0, 1, 0, 0, 0)
Q175CR3ko_F <- c(1, 1, 0, 0, 1, 0, 0, 1)
```

# Differential gene expression analysis

## Comp_2: Female vs Male : WTCR3ko
```
comp2_F.v.M_WTCR3ko <- results(dds, contrast = c(WTCR3ko_F - WTCR3ko_M))
comp2_F.v.M_WTCR3ko_shrink <- lfcShrink(dds, contrast = c(WTCR3ko_F - WTCR3ko_M), type = "ashr")
summary(comp2_F.v.M_WTCR3ko)
```

## Comp_3: Female vs Male : Q175
```
comp3_F.v.M_Q175 <- results(dds, contrast = c(Q175_F - Q175_M))
comp3_F.v.M_Q175_shrink <- lfcShrink(dds, contrast = c(Q175_F - Q175_M), type = "ashr")
summary(comp3_F.v.M_Q175)
```

## Comp_5: WTCR3ko vs WT : Male
```{r}
comp5_WTCR3ko.v.WT_Male <- results(dds, contrast = c(WTCR3ko_M - WT_M))
comp5_WTCR3ko.v.WT_Male_shrink <- lfcShrink(dds, contrast = c(WTCR3ko_M - WT_M), type = "ashr")
summary(comp5_WTCR3ko.v.WT_Male)
```

## Comp_6: WTCR3ko vs WT : Female
```{r}
comp6_WTCR3ko.v.WT_Female <- results(dds, contrast = c(WTCR3ko_F - WT_F))
comp6_WTCR3ko.v.WT_Female_shrink <- lfcShrink(dds, contrast = c(WTCR3ko_F - WT_F), type = "ashr")
summary(comp6_WTCR3ko.v.WT_Female)
```

## Comp_11: Q175CR3ko vs Q175 : Male
```
comp11_Q175CR3ko.v.Q175_Male <- results(dds, contrast = c(Q175CR3ko_M - Q175_M))
comp11_Q175CR3ko.v.Q175_Male_shrink <- lfcShrink(dds, contrast = c(Q175CR3ko_M - Q175_M), type = "ashr")
summary(comp11_Q175CR3ko.v.Q175_Male)
```

## Comp_12: Q175CR3ko vs Q175 : Female
```
comp12_Q175CR3ko.v.Q175_Female <- results(dds, contrast = c(Q175CR3ko_F - Q175_F))
comp12_Q175CR3ko.v.Q175_Female_shrink <- lfcShrink(dds, contrast = c(Q175CR3ko_F - Q175_F), type = "ashr")
summary(comp12_Q175CR3ko.v.Q175_Female)
```

## Comp_15: (Q175CR3ko-Q176) - (WTCR3ko - WT) : Male

Does the CR3 knockout in Q175 differ from CR3 knockout in WT for Males?

```
comp15_CR3koinQ175.v.CR3koinWT_Male <- results(dds,
                                               contrast = c(Q175CR3ko_M - Q175_M) - (WTCR3ko_M - WT_M))
comp15_CR3koinQ175.v.CR3koinWT_Male_shrink <- lfcShrink(dds,
                                                        contrast = c(Q175CR3ko_M - Q175_M) - (WTCR3ko_M - WT_M), type = "ashr")
summary(comp15_CR3koinQ175.v.CR3koinWT_Male)
```

## Comp_17: (WTCR3koall) - (WTall)

Does the average of the samples in WTCR3KO differ from average of the samples in WT

```
comp17_WTCR3all.v.WTall <- results(dds, contrast = c(WTCR3ko_M + WTCR3ko_F)/2 - (WT_M + WT_F)/2)
comp17_WTCR3all.v.WTall_shrink <- lfcShrink(dds, contrast = c((WTCR3ko_M + WTCR3ko_F)/2 - (WT_M + WT_F)/2), type = "ashr")
summary(comp17_WTCR3all.v.WTall)
```
