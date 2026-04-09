# **Hemoplasma epidemiological survey : R command lines and script**

We analyzed data from 626 individuals across 44 wild mammal species captured in French Guiana. The epidemiological dataset includes the following variables for each sampled individual:
- `species` : Species identity (one of the 44 wild mammal species included in this study)
- `order` : Species taxonomic order
- `hemoplasma` : Infection status with hemotropic mycoplasmas (0: Uninfected; 1: Infected)
- `anaplasmataceae` : Infection status with bacteria of the Anaplasmataceae family (*Anaplasma*, *Ehrlichia* and *Allocryptoplasma*) (0: Uninfected; 1: Infected)
- `apicomplexa` : Infection status with blood parasites, including piroplasmids (*Babesia* and *Theileria*) and haemogregarines (*Hepatozoon* and *Hemolivia*) (0: Uninfected; 1: Infected)
- `trypanosoma` : Infection status with trypanosomes (0: Uninfected; 1: Infected)
- `filaria` : Infection status with microfilariae (0: Uninfected; 1: Infected)
  
Details about all the experimental methods are available in the related manuscript.

## Table of contents 
- [Step 1. Retrieving the data](#step-1-retrieving-the-data)
- [Step 2. Prepare the data for analysis](#step-2-prepare-the-data-for-analysis)
- [Step 3. Calculate *Anaplasma* infection prevalence](#step-3-calculate-anaplasma-infection-prevalence)
- [Step 4. Test whether _Anaplasma_ infection prevalence in _Bradypus tridactylus_ (Bt) is influenced by sex, age, season, ticks and blood parasites (GLM model 1)](#step-4-test-whether-anaplasma-infection-prevalence-in-bradypus-tridactylus-bt-is-influenced-by-sex-age-season-ticks-and-blood-parasites-glm-model-1)
- [Step 5. Test whether _Anaplasma_ infection prevalence in _Choloepus didactylus_ (Cd) is influenced by sex, age, season, ticks and blood parasites (GLM model 2)](#step-5-test-whether-anaplasma-infection-prevalence-in-choloepus-didactylus-cd-is-influenced-by-sex-age-season-ticks-and-blood-parasites-glm-model-2)
- [Step 6. Test whether the proportion of sloths carrying ticks and blood parasites vary between seasons](#step-6-test-whether-the-proportion-of-sloths-carrying-ticks-and-blood-parasites-vary-between-seasons)
- [Step 7. Impact of _Anaplasma_ infections on Scale Mass Index (SMI) (GLM models 3 and 4)](#step-7-impact-of-anaplasma-infections-on-scale-mass-index-smi-glm-models-3-and-4)
- [Step 8. Impact of _Anaplasma_ infections on neck circumference (GLM models 5 and 6)](#step-8-impact-of-anaplasma-infections-on-neck-circumference-glm-models-5-and-6)
- [Step 9. Impact of _Anaplasma_ infections on hematocrit levels (GLM models 7, 8 and 9)](#step-9-impact-of-anaplasma-infections-on-hematocrit-levels-glm-models-7-8-and-9)
- [Step 10. Impact of _Anaplasma_ infections on body temperature (CLRM models 10 and 11)](#step-10-impact-of-anaplasma-infections-on-body-temperature-clrm-models-10-and-11)
- [Step 11. Impact of _Anaplasma_ infections on general health condition](#step-11-impact-of-anaplasma-infections-on-general-health-condition)
- [Step 12. Impact of _Anaplasma_ infections on female reproductive status](#step-12-impact-of-anaplasma-infections-on-female-reproductive-status)

## Step 1. Retrieving the data

All veterinary clinical data for the two sloth species are available here: (https://github.com/olivierduron/Hemoplasma_sloth_infections/blob/main/data_hemoplasma_sloth.csv)

This database will be referred to as `data_hemoplasma` throughout the R command lines and scripts provided below. It corresponds to the dataset provided in Table S1 of the related manuscript.

Load the dataset directly from the GitHub repository to R:
```
data_hemoplasma <- read.csv ("https://raw.githubusercontent.com/olivierduron/Hemoplasma_sloth_infections/main/data_hemoplasma_sloth.csv", sep = "\t")
```


## Step 2. Prepare the data for analysis

Convert categorical variables into factors:
```
data_hemoplasma$hemoplasma      <- as.factor(data_hemoplasma$hemoplasma)
data_hemoplasma$anaplasma      <- as.factor(data_hemoplasma$anaplasma)
data_hemoplasma$species        <- as.factor(data_hemoplasma$species)
data_hemoplasma$season         <- as.factor(data_hemoplasma$season)
data_hemoplasma$sex            <- as.factor(data_hemoplasma$sex)
data_hemoplasma$age            <- as.factor(data_hemoplasma$age)
data_hemoplasma$tick           <- as.factor(data_hemoplasma$tick)
data_hemoplasma$microfilaria   <- as.factor(data_hemoplasma$microfilaria)
data_hemoplasma$trypanosome    <- as.factor(data_hemoplasma$trypanosome)
data_hemoplasma$babesia        <- as.factor(data_hemoplasma$babesia)
data_hemoplasma$bloodparasite  <- as.factor(data_hemoplasma$bloodparasite)
```

Load libraries for analysis: 
```
library(binom)
library(dplyr)
library(MASS)
library(ggplot2)
library(patchwork)
library(smatr)
library(lmtest)
library(akima)
library(pwr)
library(survival)
library(RColorBrewer)
```

## Step 3. Calculate hemoplasma infection prevalence
Calculate hemoplasma infection prevalence and 95% confidence interval for _Bradypus tridactylus_ (Bt) and _Choloepus didactylus_ (Cd):

```
prevalence_results <- data_hemoplasma %>% group_by(species) %>% summarise(n = n(), positives = sum(hemoplasma == 1), prevalence = positives / n, conf_low = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$lower, conf_high = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$upper)
print(prevalence_results)
```

Results are:
```
A tibble: 2 × 6
species     n positives prevalence conf_low conf_high
  <fct>   <int>     <int>      <dbl>    <dbl>     <dbl>
1 Bt         92         4     0.0435   0.0120     0.108
2 Cd         83        68     0.819    0.720      0.895
```

