# **Hemoplasma epidemiological survey: R scripts and analysis pipeline**

We analyzed data from 626 individuals belonging to 44 species of wild mammals sampled in French Guiana. Details of sampling, examined traits and laboratory procedures are provided in the associated manuscript.

## Dataset description
-> For the epidemiological survey dataset ([here](data_hemoplasma_stat.csv)), each row corresponds to one individual and includes the following variables : 
- `species` : Species identity (44 wild mammal species)
- `order` : Taxonomic order
- `hemoplasma` : Infection status with hemotropic mycoplasmas (0 = Uninfected; 1 = Infected)
- `sex` : Sex of the individual (M = Male; F = Female)
- `anaplasmataceae` : Infection status with bacteria of the Anaplasmataceae family (*Anaplasma*, *Ehrlichia* and *Allocryptoplasma*) (0 = Uninfected; 1 = Infected)
- `apicomplexa` : Infection status with piroplasmids (*Babesia* and *Theileria*) and haemogregarines (*Hepatozoon* and *Hemolivia*) (0 = Uninfected; 1 = Infected)

-> For the life trait dataset ([here](https://github.com/olivierduron/Hemoplasma_infections/blob/main/data_mammal_traits.csv)), each row corresponds to one species and includes the following variables :
-  `species` : Species identity (44 wild mammal species)
-  `dietinv` : Percentage of the diet consisting of invertebrates (from 0 to 100%)
-  `dietvet` : Percentage of the diet consisting of vertebrates (from 0 to 100%)
-  `dietplant` : Percentage of the diet consisting of plants (fruits, leaves, seeds, nectar, etc; from 0 to 100%)
-  `strata` : Foraging stratum category (G = Ground level, including aquatic foraging; S = Scansorial; Ar = Arboreal)
-  `strataG` : Foraging activity at ground level (0 = No; 1 = Yes); species with scansorial foraging were classified as 1
-  `strataAr` : Foraging activity at arboreal level (0 = No; 1 = Yes); species with scansorial foraging were classified as 1
-  `activitynocturnal` : Foraging activity at night (0 = No; 1 = Yes)
-  `activitycrepuscular` : Foraging activity at at twilight (0 = No; 1 = Yes)
-  `activitydiurnal` : Foraging activity at day (0 = No; 1 = Yes)
-  `bodymass` : Mean adult body mass (grammes)
-  `longevity` : Mean longevity (years)
-  `femalematurity` : Mean age at maturity for females (days)
-  `littersize` : Mean litter size (number of offspring per litter)
   
## Table of contents A UPDATER !!!
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

## Step 1. Data retrieval

The epidemiological dataset is available in the GitHub repository [here](data_hemoplasma_stat.csv)
```
data_hemoplasma_stat <- read.csv2("https://raw.githubusercontent.com/olivierduron/Hemoplasma_infections/main/data_hemoplasma_stat.csv")
data_hemoplasma_stat
str(data_hemoplasma_stat)
get_modalities <- function(x) {sort(table(x), decreasing = TRUE)}
lapply(data_hemoplasma_stat, get_modalities)
```

## Step 2. Data preparation

Convert categorical variables :
```
data_hemoplasma_stat$species        <- as.factor(data_hemoplasma_stat$species)
data_hemoplasma_stat$order           <- as.factor(data_hemoplasma_stat$order)
data_hemoplasma_stat$hemoplasma      <- as.factor(data_hemoplasma_stat$hemoplasma)
data_hemoplasma_stat$sex   <- as.factor(data_hemoplasma_stat$sex)
data_hemoplasma_stat$strata   <- as.factor(data_hemoplasma_stat$strata)
data_hemoplasma_stat$activity   <- as.factor(data_hemoplasma_stat$activity)
data_hemoplasma_stat$diet   <- as.factor(data_hemoplasma_stat$diet)
data_hemoplasma_stat$sociality   <- as.factor(data_hemoplasma_stat$sociality)
data_hemoplasma_stat$anaplasmataceae      <- as.factor(data_hemoplasma_stat$anaplasmataceae)
data_hemoplasma_stat$apicomplexa         <- as.factor(data_hemoplasma_stat$apicomplexa)
```

Load required libraries : 
```
library(binom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bayestestR)
library(posterior)
library(rotl)
library(ape)
library(glmmTMB)
library(picante)
library(phytools)
library(MCMCglmm)
library(scales)
library(ggthemes)
library(ggtree)
library(lme4)
library(car)
library(emmeans)
library(brms)
library(patchwork)
```

## Step 3. Species-level summary

```
species_summary <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(as.numeric(as.character(hemoplasma)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_negative = n_sampled - n_positive,
    prevalence = n_positive / n_sampled
  ) %>%
  rowwise() %>%
  mutate(
    ci_low = binom::binom.confint(
      n_positive, n_sampled, method = "wilson"
    )$lower,
    ci_high = binom::binom.confint(
      n_positive, n_sampled, method = "wilson"
    )$upper
  ) %>%
  ungroup() %>%
  mutate(
    ci_low = ifelse(n_positive == 0, 0, ci_low)
  )
print(species_summary, n = Inf)
```

Results : 
```
# A tibble: 44 × 7
   species                   n_sampled n_positive n_negative prevalence  ci_low ci_high
   <fct>                         <int>      <dbl>      <dbl>      <dbl>   <dbl>   <dbl>
 1 Alouatta_macconnelli             22         20          2     0.909  0.722    0.975 
 2 Bradypus_tridactylus            108          4        104     0.0370 0.0145   0.0914
 3 Cabassous_unicinctus              2          0          2     0      0        0.658 
 4 Caluromys_philander               5          0          5     0      0        0.434 
 5 Sapajus_apella                    1          0          1     0      0        0.793 
 6 Choloepus_didactylus             90         72         18     0.8    0.706    0.870 
 7 Coendou_melanurus                 1          0          1     0      0        0.793 
 8 Coendou_prehensilis               3          1          2     0.333  0.0615   0.792 
 9 Cyclopes_didactylus               1          0          1     0      0        0.793 
10 Dasypus_novemcinctus             15          5         10     0.333  0.152    0.583 
11 Didelphis_marsupialis            51         24         27     0.471  0.341    0.605 
12 Eira_barbara                      4          0          4     0      0        0.490 
13 Leopardus_wiedii                  1          0          1     0      0        0.793 
14 Galictis_vittata                  4          3          1     0.75   0.301    0.954 
15 Holochilus_sciureus               5          1          4     0.2    0.0362   0.624 
16 Hydrochoerus_hydrochaeris         2          0          2     0      0        0.658 
17 Hylaeamys_megacephalus           15          0         15     0      0        0.204 
18 Hylaeamys_yunganus               10          0         10     0      0        0.278 
19 Lontra_longicaudis                1          1          0     1      0.207    1     
20 Makalata_didelphoides             8          0          8     0      0        0.324 
21 Marmosa_lepida                    1          1          0     1      0.207    1     
22 Marmosa_murina                   20          3         17     0.15   0.0524   0.360 
23 Marmosops_parvidens               5          0          5     0      0        0.434 
24 Mesomys_hispidus                 13          0         13     0      0        0.228 
25 Metachirus_nudicaudatus           5          0          5     0      0        0.434 
26 Micoureus_demerarae              16          2         14     0.125  0.0350   0.360 
27 Mus_musculus                     34          0         34     0      0        0.102 
28 Neacomys_dubosti                  1          0          1     0      0        0.793 
29 Neacomys_paracou                  8          0          8     0      0        0.324 
30 Nectomys_rattus                   4          2          2     0.5    0.150    0.850 
31 Oecomys_auyantepui               16          1         15     0.0625 0.0111   0.283 
32 Oecomys_bicolor                  16          0         16     0      0        0.194 
33 Oligoryzomys_fulvescens           7          1          6     0.143  0.0257   0.513 
34 Philander_opossum                20          9         11     0.45   0.258    0.658 
35 Pithecia_pithecia                 1          0          1     0      0        0.793 
36 Potos_flavus                      2          1          1     0.5    0.0945   0.905 
37 Proechimys_cuvieri               18          2         16     0.111  0.0310   0.328 
38 Proechimys_guyannensis           20          1         19     0.05   0.00888  0.236 
39 Puma_yagouaroundi                 5          0          5     0      0        0.434 
40 Rattus_rattus                    19          2         17     0.105  0.0294   0.314 
41 Saguinus_midas                   41         41          0     1      0.914    1     
42 Saimiri_sciureus                  1          0          1     0      0        0.793 
43 Sciurus_aestuans                  1          0          1     0      0        0.793 
44 Tamandua_tetradactyla             3          0          3     0      0        0.561
```

### Visualization
```
plot_data <- species_summary %>%
  mutate(
    species_label = gsub("_", " ", as.character(species)),
    CI_width = ci_high - ci_low
  ) %>%
  arrange(
    desc(prevalence),
    CI_width,
    species_label
  ) %>%
  mutate(
    species_label = factor(
      species_label,
      levels = rev(species_label)
    )
  )

ggplot(
  plot_data,
  aes(
    x = prevalence * 100,
    y = species_label
  )
) +
  geom_errorbarh(
    aes(
      xmin = ci_low * 100,
      xmax = ci_high * 100
    ),
    height = 0.15
  ) +
  geom_point(size = 3) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = "Hemoplasma prevalence",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    panel.grid.major.x = element_line(
      linewidth = 0.3,
      colour = "grey85"
    )
  )
```

### Spearman correlation test
```
cor.test(
  species_summary$n_sampled,
  species_summary$prevalence,
  method = "spearman",
  exact = FALSE
)
```

Results : 
```
Spearman's rank correlation rho
data:  species_summary$n_sampled and species_summary$prevalence
S = 9444.8, p-value = 0.02651
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho 
0.3344048 
```

### Interpretation
Hemoplasma prevalence increased significantly with sample size per species, suggesting that prevalence in mammals is likely underestimated in less sampled species.

### Visualization
```
plot_data <- species_summary %>%
  filter(
    is.finite(n_sampled),
    is.finite(prevalence),
    n_sampled > 0
  ) %>%
  mutate(
    log_n = log10(n_sampled)
  )

p <- ggplot(
  plot_data,
  aes(x = log_n, y = prevalence)
) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "blue",
    fill = "grey70"
  ) +
  scale_x_continuous(
    breaks = log10(c(1, 2, 5, 10, 20, 50, 100)),
    labels = c(1, 2, 5, 10, 20, 50, 100)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  theme_classic() +
  labs(
    x = "Sample size per species",
    y = "Hemoplasma prevalence"
  )
print(p)
```

## Step 4. Variation in hemoplasma infection status according to the host’s `sex` and the presence of other blood-borne pathogens (`anaplasmataceae` and `apicomplexa`)

### Data preparation: Create the `pathogens` variable by merging `anaplasmataceae` and `apicomplexa`
```
data_hemoplasma_stat <- data_hemoplasma_stat %>%
  mutate(
    pathogens = ifelse(
      anaplasmataceae == 1 | apicomplexa == 1,
      1, 0
    ),
    sex = factor(sex),
    species = factor(species)
  )
```

### Restrict all candidate models to the same complete-case dataset to ensure that likelihood-ratio tests and AIC comparisons are based on identical observations
```
model_data <- data_hemoplasma_stat %>%
  filter(
    !is.na(hemoplasma),
    !is.na(sex),
    !is.na(pathogens),
    !is.na(species),
    !is.na(anaplasmataceae),
    !is.na(apicomplexa),
  )
```

### Fit the full GLMM : model 1
Model 1 tests whether `hemoplasma` infection probability differs between sexes (`sex`) and infections by other blood-borne pathogens (`pathogens`) while accounting for species-level random effects (`1 | species`).
```
model1_a <- glmer(
  hemoplasma ~ sex * pathogens + (1 | species),
  data = model_data,
  family = binomial
)
summary(model1_a)
```

Results are:
```
Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: hemoplasma ~ sex * pathogens + (1 | species)
   Data: model_data

      AIC       BIC    logLik -2*log(L)  df.resid 
    311.1     331.0    -150.5     301.1       394 

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-3.2444 -0.2722 -0.1376  0.3538  8.0537 

Random effects:
 Groups  Name        Variance Std.Dev.
 species (Intercept) 8.045    2.836   
Number of obs: 399, groups:  species, 33

Fixed effects:
               Estimate Std. Error z value Pr(>|z|)    
(Intercept)     -2.6708     0.8035  -3.324 0.000888 ***
sexM             0.6537     0.4102   1.594 0.111022    
pathogens        1.5699     0.6340   2.476 0.013278 *  
sexM:pathogens  -0.9294     0.7785  -1.194 0.232553    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) sexM   pthgns
sexM        -0.372              
pathogens   -0.163  0.271       
sexM:pthgns  0.197 -0.513 -0.713
```

### Test and remove the `sex` × `pathogen` interaction

```
model1_b <- glmer(
  hemoplasma ~ sex + pathogens + (1 | species),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

anova(
  model1_b,
  model1_a,
  test = "Chisq"
)

AIC(model1_a, model1_b)
```

Results are:
```
Data: model_data
Models:
model1_b: hemoplasma ~ sex + pathogens + (1 | species)
model1_a: hemoplasma ~ sex * pathogens + (1 | species)
         npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
model1_b    4 310.53 326.48 -151.26    302.53                     
model1_a    5 311.07 331.02 -150.54    301.07 1.4547  1     0.2278
---
> AIC(model1_a, model1_b)
         df      AIC
model1_a  5 311.0715
model1_b  4 310.5262
```

Interpretation: The model including the `sex` × `pathogens` interaction did not provide a better fit than the additive model (likelihood-ratio test: χ² = 1.45, df = 1, p = 0.228). The additive model also had a lower AIC (310.53 vs. 311.07; ΔAIC = 0.55). The interaction was therefore removed, and the additive model was retained for subsequent model simplification.

### Test the `sex` effect
```
model1_c <- glmer(
  hemoplasma ~ pathogens + (1 | species),
  data = model_data,
  family = binomial,
)

anova(
  model1_b,
  model1_c,
  test = "Chisq"
)

AIC(model1_b, model1_c)
```

Results are: 
```
Data: model_data
Models:
model1_c: hemoplasma ~ pathogens + (1 | species)
model1_b: hemoplasma ~ sex + pathogens + (1 | species)
         npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
model1_c    3 309.90 321.86 -151.95    303.90                     
model1_b    4 310.53 326.48 -151.26    302.53 1.3717  1     0.2415
---
> AIC(model1_b, model1_c)
         df      AIC
model1_b  4 310.5262
model1_c  3 309.8978
```

Interpretation: Including `sex` did not significantly improve model fit (likelihood-ratio test: χ² = 1.37, df = 1, p = 0.242). The model without `sex` also had a lower AIC (309.90 vs. 310.53; ΔAIC = 0.63). `sex` was therefore removed, and the model including `pathogen` occurrence as the sole fixed effect, with species as a random intercept, was retained for subsequent analysis.

### Test the `pathogens` effect
```
model1_d <- glmer(
  hemoplasma ~ (1 | species),
  data = model_data,
  family = binomial,
)

anova(
  model1_c,
  model1_d,
  test = "Chisq"
)

AIC(model1_c, model1_d)
```

Results are: 
```
Data: model_data
Models:
model1_d: hemoplasma ~ (1 | species)
model1_c: hemoplasma ~ pathogens + (1 | species)
         npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)   
model1_d    2 315.06 323.04 -155.53    311.06                        
model1_c    3 309.90 321.86 -151.95    303.90 7.1671  1   0.007425 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
---
> AIC(model1_c, model1_d)
         df      AIC
model1_c  3 309.8978
model1_d  2 315.0649
```

Interpretation: Removing `pathogen` occurrence from the model significantly reduced model fit (likelihood-ratio test: χ² = 7.17, df = 1, p = 0.007). The model including `pathogen` occurrence also had a lower AIC (309.90 vs. 315.06; ΔAIC = 5.17), providing substantially stronger support for the model including `pathogen` occurrence. `pathogen` occurrence was therefore retained as a significant predictor in the minimal adequate model.

### Fit the full GLMM : model 2
Model 2 tests whether `hemoplasma` infection probability differs between `anaplasmataceae` and `apicomplexa` while accounting for species-level random effects (`1 | species`).
```
model2_a <- glmer(
  hemoplasma ~ anaplasmataceae * apicomplexa + (1 | species),
  data = model_data,
  family = binomial
)
summary(model2_a)
```

Results are:
```
Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: hemoplasma ~ anaplasmataceae * apicomplexa + (1 | species)
   Data: model_data

      AIC       BIC    logLik -2*log(L)  df.resid 
    307.2     327.2    -148.6     297.2       394 

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-2.8827 -0.2545 -0.1437  0.3469  6.9577 

Random effects:
 Groups  Name        Variance Std.Dev.
 species (Intercept) 7.581    2.753   
Number of obs: 399, groups:  species, 33

Fixed effects:
                              Estimate Std. Error z value Pr(>|z|)   
(Intercept)                    -2.2380     0.7242  -3.090   0.0020 **
anaplasmataceae1                1.1429     0.5426   2.106   0.0352 * 
apicomplexa1                    0.3152     0.7122   0.443   0.6580   
anaplasmataceae1:apicomplexa1  17.7390   128.0002   0.139   0.8898   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) anpls1 apcmp1
anaplasmtc1 -0.089              
apicomplex1 -0.038  0.055       
anplsmtc1:1 -0.001 -0.001 -0.001
```

### Test and remove the `anaplasmataceae` × `apicomplexa` interaction

```
model2_b <- glmer(
  hemoplasma ~ anaplasmataceae + apicomplexa + (1 | species),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

anova(
  model2_b,
  model2_a,
  test = "Chisq"
)

AIC(model2_a, model2_b)
```

Results are:
```
Data: model_data
Models:
model2_b: hemoplasma ~ anaplasmataceae + apicomplexa + (1 | species)
model2_a: hemoplasma ~ anaplasmataceae * apicomplexa + (1 | species)
         npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
model2_b    4 308.59 324.55 -150.30    300.59                       
model2_a    5 307.24 327.18 -148.62    297.24 3.3571  1    0.06692 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
---
> AIC(model2_a, model2_b)
         df      AIC
model2_a  5 307.2372
model2_b  4 308.5943
```

Interpretation: The model including the `anaplasmataceae` × `apicomplexa` interaction did not significantly improve model fit compared with the additive model (likelihood-ratio test: χ² = 3.36, df = 1, p = 0.067). The interaction model had a lower AIC (307.24 vs. 308.59; ΔAIC = 1.36), but not sufficient to support retaining the interaction based on the LRT criterion. The interaction was therefore removed, and the additive model was retained for subsequent model simplification.

### Test the `apicomplexa` effect
```
model2_c <- glmer(
  hemoplasma ~ anaplasmataceae + (1 | species),
  data = model_data,
  family = binomial,
)

anova(
  model2_b,
  model2_c,
  test = "Chisq"
)

AIC(model2_b, model2_c)
```

Results are: 
```
Data: model_data
Models:
model2_c: hemoplasma ~ anaplasmataceae + (1 | species)
model2_b: hemoplasma ~ anaplasmataceae + apicomplexa + (1 | species)
         npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
model2_c    3 308.12 320.09 -151.06    302.12                     
model2_b    4 308.59 324.55 -150.30    300.59 1.5256  1     0.2168
---
> AIC(model2_b, model2_c)
         df      AIC
model2_b  4 308.5943
model2_c  3 308.1199
```

Interpretation: Including `apicomplexa` did not significantly improve model fit compared with the model containing `anaplasmataceae` alone (likelihood-ratio test: χ² = 1.53, df = 1, p = 0.217). The model without `apicomplexa` also had a slightly lower AIC (308.12 vs. 308.59; ΔAIC = 0.47). `apicomplexa` infection was therefore removed, and the model including `anaplasmataceae` as the sole fixed-effect predictor, with species as a random intercept, was retained for subsequent model simplification.

### Test the `anaplasmataceae` effect
```
model2_d <- glmer(
  hemoplasma ~ (1 | species),
  data = model_data,
  family = binomial,
)

anova(
  model2_c,
  model2_d,
  test = "Chisq"
)

AIC(model2_c, model2_d)
```

Results are: 
```
Data: model_data
Models:
model2_d: hemoplasma ~ (1 | species)
model2_c: hemoplasma ~ anaplasmataceae + (1 | species)
         npar    AIC    BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)   
model2_d    2 315.06 323.04 -155.53    311.06                       
model2_c    3 308.12 320.09 -151.06    302.12 8.945  1   0.002782 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
---
> AIC(model2_c, model2_d)
         df      AIC
model2_c  3 308.1199
model2_d  2 315.0649
```

Interpretation: Removing `anaplasmataceae` significantly reduced model fit (likelihood-ratio test: χ² = 8.95, df = 1, p = 0.003). The model including `anaplasmataceae` also had a substantially lower AIC (308.12 vs. 315.06; ΔAIC = 6.94), providing strong support for retaining this predictor. `anaplasmataceae` infection was therefore retained as a significant predictor in the minimal adequate model.

### Odds ratios and 95% HDIs for `sex`, `pathogens`, `apicomplexa` and `anaplasmataceae`
```
model1_b_bayes <- brm(
  hemoplasma ~ sex + pathogens + (1 | species),
  data = model_data,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 1234
)
model2_b_bayes <- brm(
  hemoplasma ~ anaplasmataceae + apicomplexa + (1 | species),
  data = model_data,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 1234
)
posterior_model1 <- as_draws_df(model1_b_bayes)
posterior_model2 <- as_draws_df(model2_b_bayes)
or_sex <- exp(posterior_model1$b_sexM)
or_pathogens <- exp(posterior_model1$b_pathogens)
or_anaplasmataceae <- exp(
  posterior_model2$b_anaplasmataceae1
)
or_apicomplexa <- exp(
  posterior_model2$b_apicomplexa1
)
or_results <- data.frame(
  
  variable = c(
    "Sex (M vs F)",
    "Pathogens (1 vs 0)",
    "Anaplasmataceae (1 vs 0)",
    "Apicomplexa (1 vs 0)"
  ),
  OR = c(
    median(or_sex),
    median(or_pathogens),
    median(or_anaplasmataceae),
    median(or_apicomplexa)
  ),
  HDI_low = c(
    hdi(or_sex, ci = 0.95)$CI_low,
    hdi(or_pathogens, ci = 0.95)$CI_low,
    hdi(or_anaplasmataceae, ci = 0.95)$CI_low,
    hdi(or_apicomplexa, ci = 0.95)$CI_low
  ),
  HDI_high = c(
    hdi(or_sex, ci = 0.95)$CI_high,
    hdi(or_pathogens, ci = 0.95)$CI_high,
    hdi(or_anaplasmataceae, ci = 0.95)$CI_high,
    hdi(or_apicomplexa, ci = 0.95)$CI_high
  )
)
or_results
```

Results are: 
```
                  variable       OR   HDI_low  HDI_high
1             Sex (M vs F) 1.508622 0.6341916  2.826458
2       Pathogens (1 vs 0) 2.977485 1.0463991  6.634585
3 Anaplasmataceae (1 vs 0) 4.323103 1.0301907 11.702642
4     Apicomplexa (1 vs 0) 2.302248 0.2783500  7.422254
```

Interpretation: `sex` was not strongly associated with `hemoplasma` infection, with males showing higher odds of infection than females (OR = 1.51, 95% HDI: 0.63–2.83), but the HDI included 1. In contrast, individuals infected with other `pathogens` had approximately threefold higher odds of `hemoplasma` infection (OR = 2.98, 95% HDI: 1.05–6.63). `anaplasmataceae` infection was associated with an approximately fourfold increase in the odds of `hemoplasma` infection (OR = 4.32, 95% HDI: 1.03–11.70), whereas there was no association with `apicomplexa` infection (OR = 2.30, 95% HDI: 0.28–7.42), with the HDI broadly overlapping 1.

### Visualization of odds ratios and 95% HDIs for `sex`, `pathogens`, `apicomplexa` and `anaplasmataceae`
```
plot_or <- or_results %>%
  mutate(
    variable = factor(
      variable,
      levels = rev(c(
        "Sex (M vs F)",
        "Pathogens (1 vs 0)",
        "Anaplasmataceae (1 vs 0)",
        "Apicomplexa (1 vs 0)"
      ))
    )
  )

ggplot(
  plot_or,
  aes(
    x = OR,
    y = variable
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_segment(
    aes(
      x = HDI_low,
      xend = HDI_high,
      y = variable,
      yend = variable
    ),
    linewidth = 1
  ) +
  geom_point(
    size = 5
  ) +
  scale_x_log10() +
  labs(
    x = "Odds ratio (95% HDI)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11)
  )
```
### Fit the full GLMM : model 3
Model 3 tests whether `hemoplasma` infection is associated `anaplasmataceae` without the Linnaeus’s two-toed sloths (Choloepus didactylus) while accounting for species-level random effects (`1 | species`).
```
model_data_no_choloepus <- model_data %>%
  filter(species != "Choloepus_didactylus")

model3_a <- glmer(
  hemoplasma ~ anaplasmataceae + (1 | species),
  data = model_data_no_choloepus,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
summary(model3_a)
model3_b <- glmer(
  hemoplasma ~ (1 | species),
  data = model_data_no_choloepus,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
anova(
  model3_a,
  model3_b,
  test = "Chisq"
)

AIC(model3_a, model3_b)
```

Results are: 
```
Data: model_data_no_choloepus
Models:
model3_b: hemoplasma ~ (1 | species)
model3_a: hemoplasma ~ anaplasmataceae + (1 | species)
         npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
model3_b    2 230.26 237.77 -113.13    226.26                     
model3_a    3 229.92 241.19 -111.96    223.92 2.3393  1     0.1261
---
> AIC(model3_a, model3_b)
         df      AIC
model3_a  3 229.9228
model3_b  2 230.2621
```

Interpretation: After excluding Choloepus didactylus, Anaplasmataceae infection no longer significantly improved model fit (likelihood-ratio test: χ²₁ = 2.34, p = 0.126). The model including Anaplasmataceae had only a marginally lower AIC than the model without this predictor (229.92 vs. 230.26; ΔAIC = 0.34), providing little support for an independent association between Anaplasmataceae and hemoplasma infection.

### Odds ratio and 95% HDIs for `anaplasmataceae` after excluding Choloepus didactylus
```
model3_a_bayes <- brm(
  hemoplasma ~ anaplasmataceae + (1 | species),
  data = model_data_no_choloepus,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  seed = 1234
)
posterior_model3 <- as_draws_df(model3_a_bayes)
or_anaplasmataceae_model3 <- exp(
  posterior_model3$b_anaplasmataceae1
)
or_model3_results <- data.frame(
  variable = "Anaplasmataceae (1 vs 0)",
  
  OR = median(or_anaplasmataceae_model3),
  
  HDI_low = hdi(
    or_anaplasmataceae_model3,
    ci = 0.95
  )$CI_low,
  
  HDI_high = hdi(
    or_anaplasmataceae_model3,
    ci = 0.95
  )$CI_high
)
or_model3_results
```

Results are:
```
                  variable       OR   HDI_low HDI_high
1 Anaplasmataceae (1 vs 0) 3.958514 0.1981571 20.30309

```

Interpretation: After excluding Choloepus didactylus, Anaplasmataceae-positive individuals still showed higher estimated odds of hemoplasma infection (OR = 3.96), but the 95% HDI was very wide and included 1 (0.20–20.30), indicating substantial uncertainty and no clear evidence for an association.


### Fisher's exact test: hemoplasma × Anaplasmataceae in Choloepus didactylus
```
choloepus_data <- model_data %>%
  filter(species == "Choloepus_didactylus")
table_choloepus <- table(
  choloepus_data$hemoplasma,
  choloepus_data$anaplasmataceae
)
table_choloepus
fisher_choloepus <- fisher.test(
  table_choloepus
)
fisher_choloepus
```

Results are:
```
data:  table_choloepus
p-value = 0.02172
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  1.126489 28.191553
sample estimates:
odds ratio 
  4.689964 
```

Interpretation: In Choloepus didactylus, hemoplasma infection was significantly associated with Anaplasmataceae infection (Fisher’s exact test, p = 0.022).

### Visualization of odds ratios and 95% HDIs for `sex`, `pathogens`, `apicomplexa` and `anaplasmataceae`
```
or_results_extended <- bind_rows(
  or_results,
  or_model3_results %>%
    mutate(
      variable = "Anaplasmataceae (excluding Choloepus)"
    )
)
plot_or <- or_results_extended %>%
  mutate(
    variable = factor(
      variable,
      levels = rev(c(
        "Sex (M vs F)",
        "Pathogens (1 vs 0)",
        "Apicomplexa (1 vs 0)",
        "Anaplasmataceae (1 vs 0)",
        "Anaplasmataceae (excluding Choloepus)"
      ))
    )
  )
ggplot(
  plot_or,
  aes(
    x = OR,
    y = variable
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_segment(
    aes(
      x = HDI_low,
      xend = HDI_high,
      y = variable,
      yend = variable
    ),
    linewidth = 1
  ) +
  geom_point(
    size = 5
  ) +
  scale_x_log10() +
  labs(
    x = "Odds ratio (95% HDI)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11)
  )
```

### Visualization of `hemoplasma` prevalence et 95% CI (Wilson) by `anaplasmataceae` infection status in Choloepus didactylus
```
choloepus_data <- model_data %>%
  filter(species == "Choloepus_didactylus")
choloepus_prev <- choloepus_data %>%
  group_by(anaplasmataceae) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1),
    .groups = "drop"
  ) %>%
  mutate(
    prevalence = n_positive / n_sampled
  ) %>%
  rowwise() %>%
  mutate(
    ci = list(
      binom.confint(
        x = n_positive,
        n = n_sampled,
        method = "wilson"
      )
    ),
    ci_low = ci$lower,
    ci_high = ci$upper
  ) %>%
  ungroup() %>%
  mutate(
    Anaplasmataceae = ifelse(
      anaplasmataceae == 1,
      "Anaplasmataceae +",
      "Anaplasmataceae −"
    )
  ) %>%
  select(
    Anaplasmataceae,
    n_sampled,
    n_positive,
    prevalence,
    ci_low,
    ci_high
  )
choloepus_prev
ggplot(
  choloepus_prev,
  aes(
    x = prevalence * 100,
    y = Anaplasmataceae
  )
) +
  geom_segment(
    aes(
      x = ci_low * 100,
      xend = ci_high * 100,
      y = Anaplasmataceae,
      yend = Anaplasmataceae
    ),
    linewidth = 1
  ) +
  geom_point(
    size = 7
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = "Hemoplasma prevalence",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11)
  )
```

## Step 5. Phylogeny of the 44 mammalian species (Open Tree of Life & Grafen branch lengths) and other evolutionary metrics
### List of mammalian species
```
mammal_species <- c(
  "Alouatta macconnelli",
  "Saguinus midas",
  "Sapajus apella",
  "Saimiri sciureus",
  "Pithecia pithecia",
  "Bradypus tridactylus",
  "Choloepus didactylus",
  "Cyclopes didactylus",
  "Tamandua tetradactyla",
  "Cabassous unicinctus",
  "Dasypus novemcinctus",
  "Hydrochoerus hydrochaeris",
  "Holochilus sciureus",
  "Hylaeamys megacephalus",
  "Hylaeamys yunganus",
  "Neacomys dubosti",
  "Neacomys paracou",
  "Nectomys rattus",
  "Oecomys auyantepui",
  "Oecomys bicolor",
  "Oligoryzomys fulvescens",
  "Makalata didelphoides",
  "Mesomys hispidus",
  "Proechimys cuvieri",
  "Proechimys guyannensis",
  "Coendou melanurus",
  "Coendou prehensilis",
  "Mus musculus",
  "Rattus rattus",
  "Sciurus aestuans",
  "Leopardus wiedii",
  "Puma yagouaroundi",
  "Eira barbara",
  "Galictis vittata",
  "Lontra longicaudis",
  "Potos flavus",
  "Caluromys philander",
  "Didelphis marsupialis",
  "Marmosa lepida",
  "Marmosa murina",
  "Marmosops parvidens",
  "Metachirus nudicaudatus",
  "Micoureus demerarae",
  "Philander opossum"
)
```
### Match species names to the Open Tree Taxonomy, extract OTT IDs and check tree coverage
```
taxon_matches <- tnrs_match_names(
  names = mammal_species,
  context_name = "Mammals",
  do_approximate_matching = TRUE
)
taxon_matches[, c(
  "search_string",
  "unique_name",
  "score",
  "ott_id",
  "is_synonym",
  "number_matches"
)]
ott_ids <- ott_id(taxon_matches)

in_tree <- is_in_tree(ott_ids)

taxa_not_in_tree <- taxon_matches[
  !in_tree,
  c(
    "search_string",
    "unique_name",
    "score",
    "ott_id"
  )
]
```
### Extract the induced subtree and add branch lengths using Grafen's method
```
mammal_tree <- tol_induced_subtree(
  ott_ids = ott_ids[in_tree],
  label_format = "name"
)
mammal_tree$tip.label <- gsub(
  " ",
  "_",
  mammal_tree$tip.label
)
mammal_tree_grafen <- compute.brlen(
  mammal_tree,
  method = "Grafen"
)
Ntip(mammal_tree_grafen)

expected_species <- gsub(
  " ",
  "_",
  mammal_species
)
missing_species <- setdiff(
  expected_species,
  mammal_tree_grafen$tip.label
)
missing_species
```
### Plot and save the phylogeny
```
plot(
  mammal_tree_grafen,
  cex = 0.7,
  no.margin = TRUE
)
write.tree(
  mammal_tree_grafen,
  file = "mammal_phylogeny_Grafen.tre"
)
write.csv(
  taxon_matches,
  file = "OpenTree_taxon_matching.csv",
  row.names = FALSE
)
```

### Pairwise phylogenetic distances and heatmap
```
phylo_dist <- cophenetic.phylo(mammal_tree_grafen)
dim(phylo_dist)
write.csv(
  phylo_dist,
  file = "pairwise_phylogenetic_distances_Grafen.csv",
  row.names = TRUE
)
phylo_dist_long <- as.data.frame(phylo_dist) %>%
  mutate(
    species_1 = rownames(.)
  ) %>%
  pivot_longer(
    cols = -species_1,
    names_to = "species_2",
    values_to = "phylogenetic_distance"
  )
tree_order <- mammal_tree_grafen$tip.label
phylo_dist_long <- phylo_dist_long %>%
  mutate(
    species_1 = factor(
      species_1,
      levels = tree_order
    ),
    species_2 = factor(
      species_2,
      levels = rev(tree_order)
    )
  )
ggplot(
  phylo_dist_long,
  aes(
    x = species_1,
    y = species_2,
    fill = phylogenetic_distance
  )
) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "darkgreen",
    name = "Phylogenetic\ndistance"
  ) +
  coord_fixed() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7
    ),
    axis.text.y = element_text(
      size = 7
    ),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )
  "pairwise_phylogenetic_distances_heatmap.pdf",
  width = 10,
  height = 10
)
ggsave(
  "pairwise_phylogenetic_distances_heatmap.png",
  width = 10,
  height = 10,
  dpi = 300
)
```

### Evolutionary distinctiveness of the 44 mammalian species (Equal-splits method)
```
evolutionary_distinctiveness <- evol.distinct(
  mammal_tree_grafen,
  type = "equal.splits",
  scale = FALSE,
  use.branch.lengths = TRUE
)
evolutionary_distinctiveness
colnames(evolutionary_distinctiveness) <- c(
  "species",
  "evolutionary_distinctiveness"
)
evolutionary_distinctiveness <- evolutionary_distinctiveness %>%
  mutate(
    species = gsub(" ", "_", species)
  )
print(
  evolutionary_distinctiveness,
  n = 44
)
nrow(evolutionary_distinctiveness)
setdiff(
  mammal_species %>% gsub(" ", "_", .),
  evolutionary_distinctiveness$species
)
write.csv(
  evolutionary_distinctiveness,
  file = "evolutionary_distinctiveness_44_mammals.csv",
  row.names = FALSE
)
plot_ed <- evolutionary_distinctiveness %>%
  mutate(
    species = reorder(
      species,
      evolutionary_distinctiveness
    )
  )
ggplot(
  plot_ed,
  aes(
    x = evolutionary_distinctiveness,
    y = species
  )
) +
  geom_point(
    size = 4
  ) +
  labs(
    x = "Evolutionary distinctiveness",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11)
  )
ggsave(
  "evolutionary_distinctiveness_44_mammals.pdf",
  width = 8,
  height = 10
)
ggsave(
  "evolutionary_distinctiveness_44_mammals.png",
  width = 8,
  height = 10,
  dpi = 300
)
```

### Step 6. `hemoplasma` prevalence by mammalian `order`
## Observed prevalence and 95% Wilson CI by `order`
```
order_prevalence <- data_hemoplasma_stat %>%
  group_by(order) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    prevalence = n_positive / n_sampled,
    CI_low = binom.confint(
      n_positive,
      n_sampled,
      method = "wilson"
    )$lower,
    CI_high = binom.confint(
      n_positive,
      n_sampled,
      method = "wilson"
    )$upper
  ) %>%
  ungroup() %>%
  mutate(
    prevalence_percent = prevalence * 100,
    CI_low_percent = CI_low * 100,
    CI_high_percent = CI_high * 100
  )
order_prevalence
```

Results are:
```
# A tibble: 6 × 9
  order           n_sampled n_positive prevalence CI_low CI_high prevalence_percent CI_low_percent CI_high_percent
  <fct>               <int>      <int>      <dbl>  <dbl>   <dbl>              <dbl>          <dbl>           <dbl>
1 Carnivora              17          5     0.294  0.133   0.531               29.4           13.3            53.1 
2 Cingulata              17          5     0.294  0.133   0.531               29.4           13.3            53.1 
3 Didelphimorphia       123         39     0.317  0.241   0.404               31.7           24.1            40.4 
4 Pilosa                202         76     0.376  0.312   0.445               37.6           31.2            44.5 
5 Primates               66         61     0.924  0.835   0.967               92.4           83.5            96.7 
6 Rodentia              201         11     0.0547 0.0308  0.0953               5.47           3.08            9.53
```

## GLMM testing the effect of mammalian `order`
```
model_order <- glmer(
  hemoplasma ~ order + (1 | species),
  data = data_hemoplasma_stat,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa"
  )
)
summary(model_order)
model_order_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  data = data_hemoplasma_stat,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa"
  )
)
order_test <- anova(
  model_order_null,
  model_order,
  test = "Chisq"
)
order_test
AIC(
  model_order_null,
  model_order
)
```

Results are:
```
Data: data_hemoplasma_stat
Models:
model_order_null: hemoplasma ~ 1 + (1 | species)
model_order: hemoplasma ~ order + (1 | species)
                 npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
model_order_null    2 453.95 462.83 -224.97    449.95                       
model_order         7 452.52 483.60 -219.26    438.52 11.426  5    0.04355 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
---
> AIC(
+   model_order_null,
+   model_order
+ )
                 df      AIC
model_order_null  2 453.9505
model_order       7 452.5242
```

Interpretation: Hemoplasma infection probability differed among mammalian orders (GLMM, likelihood-ratio test: χ²₅ = 11.43, p = 0.044; ΔAIC = 1.43), although the improvement in model fit was modest.


## Post-hoc pairwise comparisons (odds ratios)
```
order_emmeans <- emmeans(
  model_order,
  ~ order
)
order_OR <- pairs(
  order_emmeans,
  adjust = "tukey"
)
order_OR
order_OR_results <- summary(
  order_OR,
  infer = TRUE
) %>%
  as.data.frame() %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(asymp.LCL),
    CI_high = exp(asymp.UCL)
  ) %>%
  select(
    contrast,
    OR,
    CI_low,
    CI_high,
    p.value
  )
order_OR_results
```

Results are: 
```
                      contrast          OR       CI_low     CI_high     p.value
1        Carnivora / Cingulata  1.64568553 0.0097803566  276.910237 0.999781286
2  Carnivora / Didelphimorphia  1.68480227 0.0498315840   56.963043 0.998297197
3           Carnivora / Pilosa  1.73031002 0.0277156411  108.024662 0.999002034
4         Carnivora / Primates  0.09942907 0.0015067896    6.561061 0.618377988
5         Carnivora / Rodentia  9.64033829 0.3582351226  259.427724 0.364879022
6  Cingulata / Didelphimorphia  1.02376927 0.0091349909  114.735035 1.000000000
7           Cingulata / Pilosa  1.05142203 0.0059761555  184.983187 0.999999998
8         Cingulata / Primates  0.06041802 0.0003319356   10.997126 0.640349199
9         Cingulata / Rodentia  5.85794680 0.0629637391  545.004811 0.876917680
10    Didelphimorphia / Pilosa  1.02701074 0.0275811172   38.241781 0.999999999
11  Didelphimorphia / Primates  0.05901527 0.0014915853    2.334967 0.241032395
12  Didelphimorphia / Rodentia  5.72194048 0.4180120223   78.324548 0.402198511
13           Pilosa / Primates  0.05746315 0.0008410547    3.926039 0.385397146
14           Pilosa / Rodentia  5.57145146 0.1921761003  161.524099 0.693822202
15         Primates / Rodentia 96.95694286 3.1685032997 2966.905154 0.001927086

```

Interpretation: After Tukey correction, only Primates and Rodentia differed significantly, with `hemoplasma` infection showing substantially higher odds in Rodentia than in Primates (OR = 96.96, 95% CI: 3.17–2966.91, p = 0.0019). All other pairwise comparisons were non-significant after correction.

## Step 7. Phylogenetic signal of hemoplasma prevalence (Pagel's lambda — 44 mammal species)
```
species_prev <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n_sampled,
    .groups = "drop"
  ) %>%
  mutate(
    species = as.character(species)
  )
species_prev <- species_prev %>%
  mutate(
    species = case_when(
      species == "Alouatta_macconnelli" ~
        "Alouatta_seniculus_macconnelli",
      TRUE ~ species
    )
  )
cat(
  "Number of species in data:",
  nrow(species_prev),
  "\n"
)
cat(
  "Number of species in phylogeny:",
  length(mammal_tree_grafen$tip.label),
  "\n"
)
cat("\nSpecies in data but not in phylogeny:\n")
print(
  setdiff(
    species_prev$species,
    mammal_tree_grafen$tip.label
  )
)
cat("\nSpecies in phylogeny but not in data:\n")
print(
  setdiff(
    mammal_tree_grafen$tip.label,
    species_prev$species
  )
)
stopifnot(nrow(species_prev) == 44)
stopifnot(length(mammal_tree_grafen$tip.label) == 44)
prevalence <- species_prev$prevalence
names(prevalence) <- species_prev$species
prevalence <- prevalence[
  mammal_tree_grafen$tip.label
]
stopifnot(
  identical(
    names(prevalence),
    mammal_tree_grafen$tip.label
  )
)
pagel_lambda <- phylosig(
  tree = mammal_tree_grafen,
  x = prevalence,
  method = "lambda",
  test = TRUE,
  nsim = 1000
)
pagel_lambda
```

Results are:
```
Phylogenetic signal lambda : 7.50593e-05 
logL(lambda) : -12.0333 
LR(lambda=0) : -0.000669289 
P-value (based on LR test) : 1 
```

Interpretation: Hemoplasma prevalence showed no detectable phylogenetic signal across the 44 mammalian species sampled (Pagel’s λ ≈ 0; likelihood-ratio test, p = 1), indicating that prevalence did not systematically covary with host evolutionary relatedness (closely related mammalian species did not exhibit more similar hemoplasma prevalence than expected under a model with no phylogenetic structure).

## Step 8. Exhaustive phylogenetic clade screening `hemoplasma` prevalence across 44 mammalian species

### Prepare species-level data
```
species_data <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n_sampled,
    .groups = "drop"
  ) %>%
  mutate(
    species = as.character(species)
  )
species_data <- species_data %>%
  mutate(
    species = case_when(
      species == "Alouatta_macconnelli" ~
        "Alouatta_seniculus_macconnelli",
      species == "Cebus_apella" ~
        "Sapajus_apella",
      species == "Coendou_sp" ~
        "Coendou_prehensilis",
      species == "Felis_wiedii" ~
        "Leopardus_wiedii",
      TRUE ~ species
    )
  )
stopifnot(
  all(species_data$species %in% mammal_tree_grafen$tip.label)
)
cat(
  "Number of species:",
  nrow(species_data),
  "\n"
)
```
### Function testing one internal node
```
test_clade <- function(node, tree, species_data) {
  clade_species <- extract.clade(
    tree,
    node = node
  )$tip.label
  test_data <- species_data %>%
    mutate(
      clade = ifelse(
        species %in% clade_species,
        "Clade",
        "Outside"
      ),
      clade = factor(
        clade,
        levels = c("Outside", "Clade")
      )
    )
  if (
    n_distinct(test_data$clade) < 2 ||
    sum(test_data$clade == "Clade") < 2 ||
    sum(test_data$clade == "Outside") < 2
  ) {
    return(NULL)
  }
  model_full <- tryCatch(
    glmer(
      cbind(n_positive, n_sampled - n_positive) ~
        clade + (1 | species),
      data = test_data,
      family = binomial,
      control = glmerControl(
        optimizer = "bobyqa"
      )
    ),
    error = function(e) NULL
  )

  if (is.null(model_full)) {
    return(NULL)
  }
  model_null <- tryCatch(
    glmer(
      cbind(n_positive, n_sampled - n_positive) ~
        1 + (1 | species),
      data = test_data,
      family = binomial,
      control = glmerControl(
        optimizer = "bobyqa"
      )
    ),
    error = function(e) NULL
  )

  if (is.null(model_null)) {
    return(NULL)
  }
  LRT <- anova(
    model_null,
    model_full,
    test = "Chisq"
  )

  chi2 <- LRT$Chisq[2]
  p_value <- LRT$`Pr(>Chisq)`[2]
  coef_table <- summary(model_full)$coefficients

  clade_row <- grep(
    "^clade",
    rownames(coef_table)
  )

  if (length(clade_row) != 1) {
    return(NULL)
  }

  estimate <- coef_table[
    clade_row,
    "Estimate"
  ]
  OR <- exp(estimate)
prevalence_clade <- test_data %>%
    filter(clade == "Clade") %>%
    summarise(
      prev = sum(n_positive) /
        sum(n_sampled)
    ) %>%
    pull(prev)
  prevalence_outside <- test_data %>%
    filter(clade == "Outside") %>%
    summarise(
      prev = sum(n_positive) /
        sum(n_sampled)
    ) %>%
    pull(prev)
  tibble(
    node = node,
    n_species_clade =
      sum(test_data$clade == "Clade"),
    n_species_outside =
      sum(test_data$clade == "Outside"),
    prevalence_clade =
      prevalence_clade,
    prevalence_outside =
      prevalence_outside,
    OR = OR,
    chi2 = chi2,
    p_value = p_value,
    species_clade =
      paste(
        clade_species,
        collapse = "; "
      )
  )
}
```
### Test all internal nodes
```
internal_nodes <- (
  Ntip(mammal_tree_grafen) + 1
):
(
  Ntip(mammal_tree_grafen) +
    Nnode(mammal_tree_grafen)
)
clade_results_list <- lapply(
  internal_nodes,
  test_clade,
  tree = mammal_tree_grafen,
  species_data = species_data
)
clade_results <- bind_rows(
  clade_results_list
)
```
### Multiple-testing correction
```
clade_results <- clade_results %>%
  mutate(
    p_adjusted = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  arrange(
    p_adjusted,
    p_value
  )
clade_results %>%
  select(
    node,
    n_species_clade,
    n_species_outside,
    prevalence_clade,
    prevalence_outside,
    OR,
    chi2,
    p_value,
    p_adjusted,
    species_clade
  ) %>%
  print(n = 30)
significant_clades <- clade_results %>%
  filter(
    p_adjusted < 0.05
  )
significant_clades
suggestive_clades <- clade_results %>%
  filter(
    p_adjusted < 0.10
  )
suggestive_clades
```
### Most interesting clades
```
clade_results %>%
  arrange(p_value) %>%
  select(
    node,
    n_species_clade,
    prevalence_clade,
    prevalence_outside,
    OR,
    chi2,
    p_value,
    p_adjusted,
    species_clade
  ) %>%
  print(n = 20)
```

Results are (most interesting clades):
```
    node n_species_clade prevalence_clade prevalence_outside      OR  chi2 p_value p_adjusted species_clade                                                                                                                          
   <int>           <int>            <dbl>              <dbl>   <dbl> <dbl>   <dbl>      <dbl> <chr>                                                                                                                                  
 1    68               4           0.938               0.242 6.55e+1  7.60 0.00583      0.164 Saguinus_midas; Saimiri_sciureus; Sapajus_apella; Alouatta_seniculus_macconnelli                                                       
 2    49              19           0.0547              0.438 9.23e-2  6.98 0.00823      0.164 Mus_musculus; Rattus_rattus; Oecomys_auyantepui; Oecomys_bicolor; Hylaeamys_megacephalus; Hylaeamys_yunganus; Oligoryzomys_fulvescens;…
 3    50              18           0.055               0.437 1.01e-1  6.36 0.0117       0.164 Mus_musculus; Rattus_rattus; Oecomys_auyantepui; Oecomys_bicolor; Hylaeamys_megacephalus; Hylaeamys_yunganus; Oligoryzomys_fulvescens;…
 4    67               5           0.924               0.243 3.21e+1  5.51 0.0189       0.199 Saguinus_midas; Saimiri_sciureus; Sapajus_apella; Alouatta_seniculus_macconnelli; Pithecia_pithecia                                    
 5    74               2           0.8                 0.311 1.79e+2  4.81 0.0282       0.237 Lontra_longicaudis; Galictis_vittata                                                                                                   
 6    54               4           0.0175              0.344 3.02e-2  4.07 0.0436       0.305 Oecomys_auyantepui; Oecomys_bicolor; Hylaeamys_megacephalus; Hylaeamys_yunganus                                                        
 7    56               2           0                   0.328 2.24e-9  3.52 0.0608       0.365 Hylaeamys_megacephalus; Hylaeamys_yunganus                                                                                             
 8    64               2           0                   0.326 2.52e-9  3.24 0.0717       0.376 Makalata_didelphoides; Mesomys_hispidus                                                                                                
 9    51              11           0.0519              0.387 1.53e-1  3.03 0.0817       0.381 Mus_musculus; Rattus_rattus; Oecomys_auyantepui; Oecomys_bicolor; Hylaeamys_megacephalus; Hylaeamys_yunganus; Oligoryzomys_fulvescens;…
10    69               3           0.953               0.268 2.99e+1  2.74 0.0977       0.401 Saguinus_midas; Saimiri_sciureus; Sapajus_apella                                                                                       
11    72               4           0.455               0.312 1.56e+1  2.63 0.105        0.401 Eira_barbara; Lontra_longicaudis; Galictis_vittata; Potos_flavus                                                                       
12    59               2           0                   0.319 6.21e-9  1.86 0.173        0.504 Neacomys_dubosti; Neacomys_paracou                                                                                                     
13    73               3           0.444               0.313 1.46e+1  1.84 0.175        0.504 Eira_barbara; Lontra_longicaudis; Galictis_vittata                                                                                     
14    53               9           0.0610              0.353 1.99e-1  1.83 0.176        0.504 Oecomys_auyantepui; Oecomys_bicolor; Hylaeamys_megacephalus; Hylaeamys_yunganus; Oligoryzomys_fulvescens; Neacomys_dubosti; Neacomys_p…
15    87               2           0.465               0.295 1.21e+1  1.78 0.182        0.504 Philander_opossum; Didelphis_marsupialis                                                                                               
16    75               2           0                   0.318 9.84e-9  1.57 0.211        0.504 Puma_yagouaroundi; Leopardus_wiedii                                                                                                    
17    84               2           0.190               0.319 1.67e+1  1.47 0.226        0.504 Marmosa_lepida; Marmosa_murina                                                                                                         
18    62               4           0.0508              0.342 1.42e-1  1.45 0.228        0.504 Proechimys_cuvieri; Proechimys_guyannensis; Makalata_didelphoides; Mesomys_hispidus                                                    
19    61               7           0.0615              0.344 2.03e-1  1.41 0.235        0.504 Proechimys_cuvieri; Proechimys_guyannensis; Makalata_didelphoides; Mesomys_hispidus; Coendou_melanurus; Coendou_prehensilis; Hydrochoe…
20    48              24           0.270               0.348 3.22e-1  1.38 0.240        0.504 Mus_musculus; Rattus_rattus; Oecomys_auyantepui; Oecomys_bicolor; Hylaeamys_megacephalus; Hylaeamys_yunganus; Oligoryzomys_fulvescens;…
```

Interpretation: Across the 44-species phylogeny, several clades showed nominal differences in hemoplasma prevalence, with the strongest contrasts observed for a four-species primate clade (93.8% vs. 24.2% outside; OR = 65.5, p = 0.0058) and a 19-species rodent-rich clade (5.5% vs. 43.8%; OR = 0.092, p = 0.0082). However, none remained significant after Benjamini–Hochberg correction (all adjusted p ≥ 0.164), providing no robust evidence for a specific phylogenetic clade associated with hemoplasma prevalence.

### Visualization:
```
species_prevalence <- data_hemoplasma_stat %>%
  group_by(species, order) %>%
  summarise(
    n_sampled = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n_sampled,
    .groups = "drop"
  ) %>%
  mutate(
    n_class = case_when(
      n_sampled <= 10 ~ "1–10",
      n_sampled <= 20 ~ "11–20",
      n_sampled <= 50 ~ "21–50",
      n_sampled <= 100 ~ "51–100",
      n_sampled > 100 ~ ">100"
    ),
    n_class = factor(
      n_class,
      levels = c(
        "1–10",
        "11–20",
        "21–50",
        "51–100",
        ">100"
      )
    )
  )
setdiff(
  species_prevalence$species,
  mammal_tree_grafen$tip.label
)

setdiff(
  mammal_tree_grafen$tip.label,
  species_prevalence$species
)
tree_44 <- drop.tip(
  mammal_tree_grafen,
  setdiff(
    mammal_tree_grafen$tip.label,
    species_prevalence$species
  )
)

Ntip(tree_44)

p_tree <- ggtree(
  tree_44,
  layout = "rectangular"
)

p_tree <- p_tree %<+% species_prevalence
x_max <- max(p_tree$data$x, na.rm = TRUE)
p <- p_tree +

  # Species names
  geom_tiplab(
    size = 2.7,
    hjust = 0,
    offset = 0.05
  ) +

  # Prevalence points
  geom_tippoint(
    aes(
      x = x_max + 1.5,
      colour = prevalence,
      size = n_class
    ),
    shape = 16
  ) +

  # Prevalence scale
  scale_colour_viridis_c(
    name = "Hemoplasma\nprevalence",
    limits = c(0, 1),
    labels = percent_format(accuracy = 1)
  ) +

  # Sample size classes
  scale_size_manual(
    name = "n sampled",
    values = c(
      "1–10" = 2.5,
      "11–20" = 3.5,
      "21–50" = 4.5,
      "51–100" = 5.5,
      ">100" = 7
    ),
    drop = FALSE
  ) +

  # Give space for labels and points
  xlim(
    0,
    x_max + 3
  ) +

  theme_tree2() +

  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
p
```

### Step 9. `hemoplasma` infection prevalence and mammal trait-based analyses

## Data retrieval and convert categorical variables. 

The life trait dataset is available in the GitHub repository [here](https://github.com/olivierduron/Hemoplasma_infections/blob/main/data_mammal_traits.csv).
```
data_mammal_traits <- read.csv2("https://raw.githubusercontent.com/olivierduron/Hemoplasma_infections/main/data_mammal_traits.csv")
data_mammal_traits
data_mammal_traits$species        <- as.factor(data_mammal_traits$species)
data_mammal_traits$strata        <- as.factor(data_mammal_traits$strata)
data_mammal_traits$strataG        <- as.factor(data_mammal_traits$strataG)
data_mammal_traits$strataAr        <- as.factor(data_mammal_traits$strataAr)
data_mammal_traits$activitynocturnal        <- as.factor(data_mammal_traits$activitynocturnal)
data_mammal_traits$activitycrepuscular        <- as.factor(data_mammal_traits$activitydiurnal)
data_mammal_traits$activitydiurnal        <- as.factor(data_mammal_traits$activitydiurnal)
```

# ============================================================
# SPECIES-LEVEL HAEMOPLASMA PREVALENCE ~ DIET
# BETA-BINOMIAL GLMM
# LRT + AIC + FDR + THREE-PANEL FIGURE
# ============================================================

library(dplyr)
library(ggplot2)
library(glmmTMB)
library(patchwork)


# ============================================================
# 1. Species-level dataset
# ============================================================

species_data <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_positive = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  left_join(
    data_mammal_traits %>%
      select(species, dietinv, dietvet, dietplant),
    by = "species"
  )


# ============================================================
# 2. Remove missing values
# ============================================================

data_dietinv <- species_data %>%
  filter(!is.na(dietinv))

data_dietvet <- species_data %>%
  filter(!is.na(dietvet))

data_dietplant <- species_data %>%
  filter(!is.na(dietplant))


# ============================================================
# 3. BETA-BINOMIAL MODELS
# ============================================================

# ------------------------------------------------------------
# DIETINV
# ------------------------------------------------------------

model_dietinv_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_dietinv,
  family = betabinomial(link = "logit")
)

model_dietinv <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ dietinv,
  data = data_dietinv,
  family = betabinomial(link = "logit")
)


# ------------------------------------------------------------
# DIETVET
# ------------------------------------------------------------

model_dietvet_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_dietvet,
  family = betabinomial(link = "logit")
)

model_dietvet <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ dietvet,
  data = data_dietvet,
  family = betabinomial(link = "logit")
)


# ------------------------------------------------------------
# DIETPLANT
# ------------------------------------------------------------

model_dietplant_null <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ 1,
  data = data_dietplant,
  family = betabinomial(link = "logit")
)

model_dietplant <- glmmTMB(
  cbind(n_positive, n - n_positive) ~ dietplant,
  data = data_dietplant,
  family = betabinomial(link = "logit")
)


# ============================================================
# 4. LIKELIHOOD-RATIO TESTS
# ============================================================

lrt_dietinv <- anova(
  model_dietinv_null,
  model_dietinv
)

lrt_dietvet <- anova(
  model_dietvet_null,
  model_dietvet
)

lrt_dietplant <- anova(
  model_dietplant_null,
  model_dietplant
)

print(lrt_dietinv)
print(lrt_dietvet)
print(lrt_dietplant)


# ============================================================
# 5. AIC
# ============================================================

AIC(model_dietinv_null, model_dietinv)

AIC(model_dietvet_null, model_dietvet)

AIC(model_dietplant_null, model_dietplant)


# ============================================================
# 6. Model summaries
# ============================================================

summary(model_dietinv)

summary(model_dietvet)

summary(model_dietplant)


# ============================================================
# 7. Extract statistical results
# ============================================================

# Function to extract coefficient information

extract_beta_results <- function(model, null_model, variable, lrt) {

  coef_table <- summary(model)$coefficients$cond

  estimate <- coef_table[variable, "Estimate"]

  SE <- coef_table[variable, "Std. Error"]

  z <- coef_table[variable, "z value"]

  p <- coef_table[variable, "Pr(>|z|)"]

  # LRT p-value
  LRT_p <- lrt$`Pr(>Chisq)`[2]

  # LRT statistic
  LRT_chisq <- lrt$Chisq[2]

  # AIC
  AIC_null <- AIC(null_model)
  AIC_model <- AIC(model)

  data.frame(
    variable = variable,
    n_species = nobs(model),
    estimate = estimate,
    SE = SE,
    z = z,
    p_coefficient = p,
    LRT_chisq = LRT_chisq,
    LRT_p = LRT_p,
    AIC_null = AIC_null,
    AIC_model = AIC_model,
    delta_AIC = AIC_null - AIC_model,
    CI_low = estimate - 1.96 * SE,
    CI_high = estimate + 1.96 * SE,
    OR = exp(estimate),
    OR_low = exp(estimate - 1.96 * SE),
    OR_high = exp(estimate + 1.96 * SE)
  )
}


# Extract results

results_diet <- bind_rows(

  extract_beta_results(
    model_dietinv,
    model_dietinv_null,
    "dietinv",
    lrt_dietinv
  ),

  extract_beta_results(
    model_dietvet,
    model_dietvet_null,
    "dietvet",
    lrt_dietvet
  ),

  extract_beta_results(
    model_dietplant,
    model_dietplant_null,
    "dietplant",
    lrt_dietplant
  )
)


# ============================================================
# 8. FDR correction
# ============================================================

results_diet$LRT_p_FDR <-
  p.adjust(
    results_diet$LRT_p,
    method = "BH"
  )

results_diet$p_coefficient_FDR <-
  p.adjust(
    results_diet$p_coefficient,
    method = "BH"
  )


# ============================================================
# 9. Display final results
# ============================================================

print(results_diet)
```

Results : 
Species-level haemoplasma prevalence was tested against three dietary variables using beta-binomial models, accounting for overdispersion in the number of infected individuals among species.
- Invertebrate consumption (dietinv): no association with haemoplasma prevalence (LRT: χ² = 0.36, p = 0.550; β = 0.0051, p = 0.543; ΔAIC = −1.64). The dietary model was not supported over the null model.
- Vertebrate consumption (dietvet): no significant association (LRT: χ² = 2.03, p = 0.155; β = 0.0138, p = 0.147; ΔAIC = +0.03). The slight positive effect was not supported after FDR correction (p_FDR = 0.232).
- Plant consumption (dietplant): no significant association (LRT: χ² = 2.19, p = 0.139; β = −0.0102, p = 0.135; ΔAIC = +0.19). The negative trend was not significant after FDR correction (p_FDR = 0.232).
After Benjamini–Hochberg correction for the three dietary tests, none of the dietary variables was significantly associated with haemoplasma prevalence (all LRT p_FDR ≥ 0.232).

Interpretation : 
These results provide no evidence that species-level haemoplasma prevalence is associated with dietary composition, whether considering the proportion of invertebrates, vertebrates, or plants in the diet. Although vertebrate consumption showed a weak positive trend and plant consumption a weak negative trend, neither was statistically supported.

Visualisation
```
# ============================================================
# 10. PREPARE DATA FOR PLOTTING
# ============================================================

# Function for predictions and 95% CI

get_predictions <- function(model, data, variable, label) {

  x <- seq(
    min(data[[variable]], na.rm = TRUE),
    max(data[[variable]], na.rm = TRUE),
    length.out = 100
  )

  newdata <- data.frame(x)
  names(newdata) <- variable

  # Prediction on link scale
  pred <- predict(
    model,
    newdata = newdata,
    type = "link",
    se.fit = TRUE
  )

  newdata$fit <- plogis(pred$fit)

  newdata$lower <- plogis(
    pred$fit - 1.96 * pred$se.fit
  )

  newdata$upper <- plogis(
    pred$fit + 1.96 * pred$se.fit
  )

  newdata$variable <- label

  return(newdata)
}


# Predictions

pred_dietinv <- get_predictions(
  model_dietinv,
  data_dietinv,
  "dietinv",
  "Invertebrates"
)

pred_dietvet <- get_predictions(
  model_dietvet,
  data_dietvet,
  "dietvet",
  "Vertebrates"
)

pred_dietplant <- get_predictions(
  model_dietplant,
  data_dietplant,
  "dietplant",
  "Plants"
)


# ============================================================
# 11. OBSERVED DATA
# ============================================================

plot_data <- bind_rows(

  data_dietinv %>%
    mutate(
      diet = dietinv,
      variable = "Invertebrates"
    ),

  data_dietvet %>%
    mutate(
      diet = dietvet,
      variable = "Vertebrates"
    ),

  data_dietplant %>%
    mutate(
      diet = dietplant,
      variable = "Plants"
    )
)


pred_data <- bind_rows(

  pred_dietinv %>%
    mutate(diet = dietinv),

  pred_dietvet %>%
    mutate(diet = dietvet),

  pred_dietplant %>%
    mutate(diet = dietplant)
)


# ============================================================
# 12. THREE-PANEL FIGURE
# ============================================================

figure_diet <- ggplot(
  plot_data,
  aes(x = diet, y = prevalence)
) +

  # 95% CI
  geom_ribbon(
    data = pred_data,
    aes(
      x = diet,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    alpha = 0.20
  ) +

  # Predicted relationship
  geom_line(
    data = pred_data,
    aes(
      x = diet,
      y = fit
    ),
    inherit.aes = FALSE,
    linewidth = 1
  ) +

  # Observed prevalence
  geom_point(
    aes(size = n),
    alpha = 0.70
  ) +

  facet_wrap(
    ~ variable,
    nrow = 1,
    scales = "free_x"
  ) +

  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +

  labs(
    x = "Diet composition (%)",
    y = "Haemoplasma prevalence",
    size = "Number tested"
  ) +

  theme_classic() +

  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )


print(figure_diet)


# ============================================================
# 13. SAVE FIGURE
# ============================================================

ggsave(
  "haemoplasma_prevalence_diet_beta_binomial.png",
  figure_diet,
  width = 11,
  height = 4.5,
  dpi = 300
)

```


# ============================================================
# SPECIES-LEVEL HAEMOPLASMA PREVALENCE
# ~ FORAGING STRATA
# Beta-binomial GLM + LRT + AIC + FDR + predictions
# ============================================================
```
library(dplyr)
library(glmmTMB)
library(ggplot2)


# ============================================================
# 1. PREPARE SPECIES-LEVEL DATA
# ============================================================

species_data_strata <-
  
  data_hemoplasma_stat %>%
  
  group_by(species) %>%
  
  summarise(
    
    n = n(),
    
    n_positive =
      sum(
        hemoplasma == 1,
        na.rm = TRUE
      ),
    
    prevalence =
      n_positive / n,
    
    .groups = "drop"
  ) %>%
  
  left_join(
    
    data_mammal_traits %>%
      
      select(
        species,
        strataAr,
        strataG
      ),
    
    by = "species"
  )


# ------------------------------------------------------------
# Number of species
# ------------------------------------------------------------

cat(
  "\nNumber of species =",
  nrow(species_data_strata),
  "\n"
)


# ------------------------------------------------------------
# Missing values
# ------------------------------------------------------------

cat(
  "\nMissing values:\n"
)

print(
  colSums(
    is.na(
      species_data_strata[
        ,
        c(
          "strataAr",
          "strataG"
        )
      ]
    )
  )
)


# ============================================================
# 2. DATASETS FOR EACH VARIABLE
# ============================================================

data_strataAr <-
  
  species_data_strata %>%
  
  filter(
    !is.na(strataAr)
  )


data_strataG <-
  
  species_data_strata %>%
  
  filter(
    !is.na(strataG)
  )


cat(
  "\nSpecies used for strataAr =",
  nrow(data_strataAr),
  "\n"
)

cat(
  "Species used for strataG  =",
  nrow(data_strataG),
  "\n"
)


# ============================================================
# 3. CHECK FACTOR LEVELS
# ============================================================

cat(
  "\n================ STRATA AR LEVELS ================\n"
)

print(
  levels(
    data_strataAr$strataAr
  )
)

print(
  table(
    data_strataAr$strataAr
  )
)


cat(
  "\n================ STRATA G LEVELS ================\n"
)

print(
  levels(
    data_strataG$strataG
  )
)

print(
  table(
    data_strataG$strataG
  )
)


# Expected:
#
# strataAr:
# 0 = 20 species
# 1 = 20 species
#
# strataG:
# 0 = 14 species
# 1 = 26 species


# ============================================================
# 4. NULL MODELS
# ============================================================

model_strataAr_null <-
  
  glmmTMB(
    
    cbind(
      n_positive,
      n - n_positive
    ) ~ 1,
    
    data = data_strataAr,
    
    family =
      betabinomial(
        link = "logit"
      )
  )


model_strataG_null <-
  
  glmmTMB(
    
    cbind(
      n_positive,
      n - n_positive
    ) ~ 1,
    
    data = data_strataG,
    
    family =
      betabinomial(
        link = "logit"
      )
  )


# ============================================================
# 5. FULL MODELS
# ============================================================

model_strataAr <-
  
  glmmTMB(
    
    cbind(
      n_positive,
      n - n_positive
    ) ~ strataAr,
    
    data = data_strataAr,
    
    family =
      betabinomial(
        link = "logit"
      )
  )


model_strataG <-
  
  glmmTMB(
    
    cbind(
      n_positive,
      n - n_positive
    ) ~ strataG,
    
    data = data_strataG,
    
    family =
      betabinomial(
        link = "logit"
      )
  )


# ============================================================
# 6. LIKELIHOOD-RATIO TESTS
# ============================================================

lrt_strataAr <-
  
  anova(
    model_strataAr_null,
    model_strataAr
  )


lrt_strataG <-
  
  anova(
    model_strataG_null,
    model_strataG
  )


cat(
  "\n================ STRATA AR LRT ================\n"
)

print(
  lrt_strataAr
)


cat(
  "\n================ STRATA G LRT ================\n"
)

print(
  lrt_strataG
)


# ============================================================
# 7. AIC
# ============================================================

cat(
  "\n================ AIC STRATA AR ================\n"
)

print(
  AIC(
    model_strataAr_null,
    model_strataAr
  )
)


cat(
  "\n================ AIC STRATA G ================\n"
)

print(
  AIC(
    model_strataG_null,
    model_strataG
  )
)


# ============================================================
# 8. MODEL SUMMARIES
# ============================================================

cat(
  "\n================ STRATA AR ================\n"
)

print(
  summary(
    model_strataAr
  )
)


cat(
  "\n================ STRATA G ================\n"
)

print(
  summary(
    model_strataG
  )
)


# ============================================================
# 9. EXTRACT RESULTS
# ============================================================

extract_strata_results <-
  
  function(
    model,
    null_model,
    variable,
    lrt
  ) {
    
    
    # --------------------------------------------------------
    # Coefficient table
    # --------------------------------------------------------
    
    coef_table <-
      summary(model)$coefficients$cond
    
    
    # --------------------------------------------------------
    # Find predictor coefficient automatically
    #
    # For factor predictors:
    # strataAr -> strataAr1
    # strataG  -> strataG1
    # --------------------------------------------------------
    
    coefficient_name <-
      
      setdiff(
        rownames(coef_table),
        "(Intercept)"
      )[1]
    
    
    coef_row <-
      
      coef_table[
        coefficient_name,
        ,
        drop = FALSE
      ]
    
    
    # --------------------------------------------------------
    # Coefficient statistics
    # --------------------------------------------------------
    
    estimate <-
      coef_row[
        1,
        "Estimate"
      ]
    
    
    SE <-
      coef_row[
        1,
        "Std. Error"
      ]
    
    
    z <-
      coef_row[
        1,
        "z value"
      ]
    
    
    p <-
      coef_row[
        1,
        "Pr(>|z|)"
      ]
    
    
    # --------------------------------------------------------
    # Likelihood-ratio test
    # --------------------------------------------------------
    
    LRT_chisq <-
      lrt$Chisq[2]
    
    
    LRT_p <-
      lrt$`Pr(>Chisq)`[2]
    
    
    # --------------------------------------------------------
    # AIC
    # --------------------------------------------------------
    
    AIC_null <-
      AIC(
        null_model
      )
    
    
    AIC_model <-
      AIC(
        model
      )
    
    
    # --------------------------------------------------------
    # 95% CI
    # --------------------------------------------------------
    
    CI_low <-
      estimate -
      1.96 * SE
    
    
    CI_high <-
      estimate +
      1.96 * SE
    
    
    # --------------------------------------------------------
    # Odds ratio
    # --------------------------------------------------------
    
    OR <-
      exp(
        estimate
      )
    
    
    OR_low <-
      exp(
        CI_low
      )
    
    
    OR_high <-
      exp(
        CI_high
      )
    
    
    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------
    
    data.frame(
      
      variable =
        variable,
      
      coefficient =
        coefficient_name,
      
      n_species =
        nobs(model),
      
      estimate =
        estimate,
      
      SE =
        SE,
      
      z =
        z,
      
      p_coefficient =
        p,
      
      LRT_chisq =
        LRT_chisq,
      
      LRT_p =
        LRT_p,
      
      AIC_null =
        AIC_null,
      
      AIC_model =
        AIC_model,
      
      delta_AIC =
        AIC_null -
        AIC_model,
      
      CI_low =
        CI_low,
      
      CI_high =
        CI_high,
      
      OR =
        OR,
      
      OR_low =
        OR_low,
      
      OR_high =
        OR_high
    )
  }


# ============================================================
# 10. COMBINE RESULTS
# ============================================================

results_strata <-
  
  bind_rows(
    
    extract_strata_results(
      
      model_strataAr,
      
      model_strataAr_null,
      
      "strataAr",
      
      lrt_strataAr
    ),
    
    
    extract_strata_results(
      
      model_strataG,
      
      model_strataG_null,
      
      "strataG",
      
      lrt_strataG
    )
  )


# ============================================================
# 11. FDR CORRECTION
# ============================================================

results_strata$LRT_p_FDR <-
  
  p.adjust(
    results_strata$LRT_p,
    method = "BH"
  )


results_strata$p_coefficient_FDR <-
  
  p.adjust(
    results_strata$p_coefficient,
    method = "BH"
  )


# ============================================================
# 12. FINAL RESULTS
# ============================================================

cat(
  "\n================ FINAL STRATA RESULTS ================\n"
)

print(
  results_strata,
  row.names = FALSE
)
```
```
Results: 
Neither arboreal foraging nor ground-level foraging was associated with haemoplasma prevalence. Arboreal foraging showed no significant effect (LRT: χ² = 0.73, p = 0.393; OR = 1.58, 95% CI: 0.56–4.45), nor did ground-level foraging (LRT: χ² = 0.81, p = 0.368; OR = 0.59, 95% CI: 0.19–1.81). Neither model improved on the null model (ΔAIC = −1.27 and −1.19, respectively).

Interpretation: 
Foraging stratum was not associated with haemoplasma prevalence among the 40 species with available trait data.
```
```
# ============================================================
# 13. PREDICTION FUNCTION
# ============================================================

get_predictions_binary <-
  
  function(
    model,
    data,
    variable,
    label
  ) {
    
    
    # --------------------------------------------------------
    # Preserve factor levels used by the model
    # --------------------------------------------------------
    
    factor_levels <-
      levels(
        data[[variable]]
      )
    
    
    # --------------------------------------------------------
    # Prediction dataset
    # --------------------------------------------------------
    
    newdata <-
      
      data.frame(
        
        x =
          factor(
            factor_levels,
            levels = factor_levels
          )
      )
    
    
    names(newdata)[1] <-
      variable
    
    
    # --------------------------------------------------------
    # Predictions on response scale
    # --------------------------------------------------------
    
    pred <-
      
      predict(
        
        model,
        
        newdata =
          newdata,
        
        type = "response",
        
        se.fit = TRUE
      )
    
    
    newdata$predicted <-
      as.numeric(
        pred$fit
      )
    
    
    newdata$SE <-
      as.numeric(
        pred$se.fit
      )
    
    
    # --------------------------------------------------------
    # 95% CI on logit scale
    # --------------------------------------------------------
    
    newdata$CI_low <-
      
      plogis(
        
        qlogis(
          newdata$predicted
        ) -
          1.96 *
          newdata$SE
      )
    
    
    newdata$CI_high <-
      
      plogis(
        
        qlogis(
          newdata$predicted
        ) +
          1.96 *
          newdata$SE
      )
    
    
    # --------------------------------------------------------
    # Numeric x for plotting
    # --------------------------------------------------------
    
    newdata$x <-
      
      as.numeric(
        as.character(
          newdata[[variable]]
        )
      )
    
    
    # --------------------------------------------------------
    # Label
    # --------------------------------------------------------
    
    newdata$variable <-
      label
    
    
    newdata
  }


# ============================================================
# 14. PREDICTIONS
# ============================================================

pred_strataAr <-
  
  get_predictions_binary(
    
    model_strataAr,
    
    data_strataAr,
    
    "strataAr",
    
    "Arboreal foraging"
  )


pred_strataG <-
  
  get_predictions_binary(
    
    model_strataG,
    
    data_strataG,
    
    "strataG",
    
    "Ground-level foraging"
  )


# Combine predictions

pred_strata <-
  
  bind_rows(
    
    pred_strataAr,
    
    pred_strataG
  )


cat(
  "\n================ PREDICTIONS ================\n"
)

print(
  pred_strata
)


# ============================================================
# 15. OBSERVED SPECIES-LEVEL PREVALENCE
# ============================================================

# ------------------------------------------------------------
# STRATA AR
# ------------------------------------------------------------

plot_data_strataAr <-
  
  data_strataAr %>%
  
  mutate(
    
    variable =
      "Arboreal foraging",
    
    x =
      as.numeric(
        as.character(
          strataAr
        )
      ),
    
    group =
      factor(
        x,
        levels = c(0, 1),
        labels = c(
          "Non-arboreal",
          "Arboreal"
        )
      )
  ) %>%
  
  select(
    species,
    n,
    n_positive,
    prevalence,
    x,
    group,
    variable
  )


# ------------------------------------------------------------
# STRATA G
# ------------------------------------------------------------

plot_data_strataG <-
  
  data_strataG %>%
  
  mutate(
    
    variable =
      "Ground-level foraging",
    
    x =
      as.numeric(
        as.character(
          strataG
        )
      ),
    
    group =
      factor(
        x,
        levels = c(0, 1),
        labels = c(
          "Non-ground",
          "Ground-level"
        )
      )
  ) %>%
  
  select(
    species,
    n,
    n_positive,
    prevalence,
    x,
    group,
    variable
  )


# ------------------------------------------------------------
# Combine observed data
# ------------------------------------------------------------

plot_data_strata <-
  
  bind_rows(
    
    plot_data_strataAr,
    
    plot_data_strataG
  )


cat(
  "\nNumber of observed species-level data points =",
  nrow(plot_data_strata),
  "\n"
)


# ============================================================
# 16. PREPARE PREDICTION LABELS
# ============================================================

pred_strata <-
  
  pred_strata %>%
  
  mutate(
    
    group =
      
      case_when(
        
        variable ==
          "Arboreal foraging" &
          x == 0 ~
          "Non-arboreal",
        
        variable ==
          "Arboreal foraging" &
          x == 1 ~
          "Arboreal",
        
        variable ==
          "Ground-level foraging" &
          x == 0 ~
          "Non-ground",
        
        variable ==
          "Ground-level foraging" &
          x == 1 ~
          "Ground-level"
      )
  )


# ============================================================
# 17. GRAPH
# ============================================================

p_strata <-
  
  ggplot() +
  
  
  # --------------------------------------------------------
  # Observed species-level prevalence
  # --------------------------------------------------------
  
  geom_jitter(
    
    data =
      plot_data_strata,
    
    aes(
      x = group,
      y = prevalence
    ),
    
    width = 0.08,
    
    height = 0,
    
    alpha = 0.5,
    
    size = 2
  ) +
  
  
  # --------------------------------------------------------
  # Model-predicted prevalence
  # --------------------------------------------------------
  
  geom_point(
    
    data =
      pred_strata,
    
    aes(
      x = group,
      y = predicted
    ),
    
    size = 3
  ) +
  
  
  # --------------------------------------------------------
  # 95% CI
  # --------------------------------------------------------
  
  geom_errorbar(
    
    data =
      pred_strata,
    
    aes(
      x = group,
      ymin = CI_low,
      ymax = CI_high
    ),
    
    width = 0.08
  ) +
  
  
  # --------------------------------------------------------
  # Separate panels
  # --------------------------------------------------------
  
  facet_wrap(
    
    ~ variable,
    
    nrow = 1,
    
    scales = "free_x"
  ) +
  
  
  # --------------------------------------------------------
  # Y axis
  # --------------------------------------------------------
  
  scale_y_continuous(
    
    limits = c(0, 1),
    
    labels =
      scales::percent
  ) +
  
  
  # --------------------------------------------------------
  # Labels
  # --------------------------------------------------------
  
  labs(
    
    x =
      "Foraging stratum",
    
    y =
      "Haemoplasma prevalence"
  ) +
  
  
  # --------------------------------------------------------
  # Theme
  # --------------------------------------------------------
  
  theme_classic() +
  
  theme(
    
    strip.background =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    axis.text.x =
      element_text(
        angle = 0,
        hjust = 0.5
      )
  )


# ============================================================
# 18. DISPLAY FIGURE
# ============================================================

print(
  p_strata
)


# ============================================================
# 19. SAVE FIGURE
# ============================================================

ggsave(
  
  "haemoplasma_prevalence_strata.png",
  
  p_strata,
  
  width = 8,
  
  height = 4,
  
  dpi = 300
)


# ============================================================
# 20. SAVE STATISTICAL RESULTS
# ============================================================

write.csv(
  
  results_strata,
  
  "haemoplasma_prevalence_strata_results.csv",
  
  row.names = FALSE
)
```

#ALTERNATIVE AVEC JUSTE STRATA:
# ============================================================
# SPECIES-LEVEL HAEMOPLASMA PREVALENCE ~ FORAGING STRATUM
# Beta-binomial GLM + LRT + AIC + pairwise comparisons
# + FDR + predictions + figure
#
# strata:
# G  = Ground level, including aquatic foraging
# S  = Scansorial
# Ar = Arboreal
# ============================================================
```
library(dplyr)
library(glmmTMB)
library(ggplot2)
library(emmeans)


# ============================================================
# 1. PREPARE SPECIES-LEVEL DATA
# ============================================================

# G = reference category

data_mammal_traits$strata <- factor(
  data_mammal_traits$strata,
  levels = c("G", "S", "Ar")
)


species_data_strata <- data_hemoplasma_stat %>%
  
  group_by(species) %>%
  
  summarise(
    n = n(),
    n_positive = sum(
      hemoplasma == 1,
      na.rm = TRUE
    ),
    prevalence = n_positive / n,
    .groups = "drop"
  ) %>%
  
  left_join(
    data_mammal_traits %>%
      select(
        species,
        strata
      ),
    by = "species"
  )


# ============================================================
# 2. CHECK DATA
# ============================================================

cat(
  "\nNumber of species =",
  nrow(species_data_strata),
  "\n"
)

cat(
  "\nMissing strata =",
  sum(
    is.na(
      species_data_strata$strata
    )
  ),
  "\n"
)


data_strata <- species_data_strata %>%
  
  filter(
    !is.na(strata)
  )


cat(
  "\nSpecies used in analysis =",
  nrow(data_strata),
  "\n"
)


cat(
  "\nDistribution of foraging strata:\n"
)

print(
  table(
    data_strata$strata
  )
)


# ============================================================
# 3. NULL MODEL
# ============================================================

model_strata_null <- glmmTMB(
  
  cbind(
    n_positive,
    n - n_positive
  ) ~ 1,
  
  data = data_strata,
  
  family = betabinomial(
    link = "logit"
  )
)


# ============================================================
# 4. FULL MODEL
# ============================================================

model_strata <- glmmTMB(
  
  cbind(
    n_positive,
    n - n_positive
  ) ~ strata,
  
  data = data_strata,
  
  family = betabinomial(
    link = "logit"
  )
)


# ============================================================
# 5. LIKELIHOOD-RATIO TEST
# ============================================================

lrt_strata <- anova(
  
  model_strata_null,
  model_strata
  
)


cat(
  "\n================ STRATA LRT ================\n"
)

print(
  lrt_strata
)


# ============================================================
# 6. AIC
# ============================================================

cat(
  "\n================ AIC ================\n"
)

print(
  AIC(
    model_strata_null,
    model_strata
  )
)


# ============================================================
# 7. MODEL SUMMARY
# ============================================================

cat(
  "\n================ STRATA MODEL ================\n"
)

print(
  summary(
    model_strata
  )
)


# ============================================================
# 8. EXTRACT MODEL COEFFICIENTS
# ============================================================

coef_table <-
  summary(
    model_strata
  )$coefficients$cond


# ------------------------------------------------------------
# S vs G
# ------------------------------------------------------------

estimate_S <-
  coef_table[
    "strataS",
    "Estimate"
  ]

SE_S <-
  coef_table[
    "strataS",
    "Std. Error"
  ]

z_S <-
  coef_table[
    "strataS",
    "z value"
  ]

p_S <-
  coef_table[
    "strataS",
    "Pr(>|z|)"
  ]


# ------------------------------------------------------------
# Ar vs G
# ------------------------------------------------------------

estimate_Ar <-
  coef_table[
    "strataAr",
    "Estimate"
  ]

SE_Ar <-
  coef_table[
    "strataAr",
    "Std. Error"
  ]

z_Ar <-
  coef_table[
    "strataAr",
    "z value"
  ]

p_Ar <-
  coef_table[
    "strataAr",
    "Pr(>|z|)"
  ]


# ============================================================
# 9. GLOBAL LRT / AIC
# ============================================================

LRT_chisq <-
  lrt_strata$Chisq[2]

LRT_p <-
  lrt_strata$`Pr(>Chisq)`[2]

AIC_null <-
  AIC(
    model_strata_null
  )

AIC_model <-
  AIC(
    model_strata
  )

delta_AIC <-
  AIC_null - AIC_model


# ============================================================
# 10. ODDS RATIOS + 95% CI
# ============================================================

results_strata <- data.frame(
  
  comparison = c(
    "S vs G",
    "Ar vs G"
  ),
  
  n_species = nobs(
    model_strata
  ),
  
  estimate = c(
    estimate_S,
    estimate_Ar
  ),
  
  SE = c(
    SE_S,
    SE_Ar
  ),
  
  z = c(
    z_S,
    z_Ar
  ),
  
  p_coefficient = c(
    p_S,
    p_Ar
  ),
  
  CI_low = c(
    estimate_S - 1.96 * SE_S,
    estimate_Ar - 1.96 * SE_Ar
  ),
  
  CI_high = c(
    estimate_S + 1.96 * SE_S,
    estimate_Ar + 1.96 * SE_Ar
  )
)


# Odds ratios

results_strata$OR <-
  exp(
    results_strata$estimate
  )

results_strata$OR_low <-
  exp(
    results_strata$CI_low
  )

results_strata$OR_high <-
  exp(
    results_strata$CI_high
  )


# Global LRT and AIC

results_strata$LRT_chisq <-
  LRT_chisq

results_strata$LRT_p <-
  LRT_p

results_strata$AIC_null <-
  AIC_null

results_strata$AIC_model <-
  AIC_model

results_strata$delta_AIC <-
  delta_AIC


# ============================================================
# 11. FDR CORRECTION FOR MODEL COEFFICIENTS
# ============================================================

results_strata$p_coefficient_FDR <-
  p.adjust(
    results_strata$p_coefficient,
    method = "BH"
  )


# ============================================================
# 12. DISPLAY MODEL RESULTS
# ============================================================

cat(
  "\n================ MODEL RESULTS ================\n"
)

print(
  results_strata,
  row.names = FALSE
)


# ============================================================
# 13. ESTIMATED PREVALENCE BY STRATUM
# ============================================================

emm_strata <- emmeans(
  
  model_strata,
  
  ~ strata,
  
  type = "response"
)


cat(
  "\n================ ESTIMATED PREVALENCE ================\n"
)

print(
  emm_strata
)


# ============================================================
# 14. ALL PAIRWISE COMPARISONS
# ============================================================

pairwise_strata <- pairs(
  
  emm_strata,
  
  adjust = "tukey"
)


cat(
  "\n================ PAIRWISE COMPARISONS ================\n"
)

print(
  pairwise_strata
)


# ============================================================
# 15. PAIRWISE RESULTS + 95% CI + FDR
# ============================================================

pairwise_results <- as.data.frame(
  
  summary(
    pairwise_strata,
    infer = c(
      TRUE,
      TRUE
    )
  )
)


# Additional BH correction

pairwise_results$p_FDR <-
  p.adjust(
    pairwise_results$p.value,
    method = "BH"
  )


cat(
  "\n================ PAIRWISE RESULTS ================\n"
)

print(
  pairwise_results,
  row.names = FALSE
)


# ============================================================
# 16. PREDICTED PREVALENCE
# ============================================================

newdata_strata <- data.frame(
  
  strata = factor(
    
    c(
      "G",
      "S",
      "Ar"
    ),
    
    levels = c(
      "G",
      "S",
      "Ar"
    )
  )
)


# Predictions on response scale

newdata_strata$predicted <-
  predict(
    
    model_strata,
    
    newdata = newdata_strata,
    
    type = "response"
  )


# ============================================================
# 17. 95% CI FOR PREDICTIONS
# ============================================================

pred_link <-
  predict(
    
    model_strata,
    
    newdata = newdata_strata,
    
    type = "link",
    
    se.fit = TRUE
  )


newdata_strata$CI_low <-
  plogis(
    
    pred_link$fit -
      1.96 *
      pred_link$se.fit
  )


newdata_strata$CI_high <-
  plogis(
    
    pred_link$fit +
      1.96 *
      pred_link$se.fit
  )


cat(
  "\n================ PREDICTED PREVALENCE ================\n"
)

print(
  newdata_strata
)


# ============================================================
# 18. OBSERVED SPECIES-LEVEL PREVALENCE
# ============================================================

plot_data_strata <- data_strata %>%
  
  mutate(
    
    strata = factor(
      
      strata,
      
      levels = c(
        "G",
        "S",
        "Ar"
      )
    )
  )
```
Results — foraging stratum: 
Haemoplasma prevalence did not differ significantly among foraging strata (beta-binomial GLM: LRT χ²₂ = 0.917, p = 0.632; ΔAIC = +3.08 relative to the null model). Estimated prevalence was 18.0% (95% CI: 9.3–31.9%) for ground-foraging species, 21.9% (7.4–49.7%) for scansorial species, and 28.2% (13.6–49.6%) for arboreal species.
Pairwise comparisons likewise provided no evidence for differences between Ground vs Scansorial (OR = 0.78, Tukey-adjusted p = 0.942), Ground vs Arboreal (OR = 0.56, p = 0.594), or Scansorial vs Arboreal (OR = 0.71, p = 0.903).

Interpretation: 
Overall, foraging stratum was not associated with haemoplasma prevalence. Although prevalence showed a descriptive increase from ground (18%) to scansorial (22%) and arboreal species (28%), the confidence intervals were broad and strongly overlapping, particularly for the smaller scansorial group (n = 6 species). The higher prevalence observed among arboreal species should therefore not be interpreted as evidence of an ecological effect of foraging stratum.
```
# ============================================================
# 19. GRAPH
# ============================================================

p_strata <- ggplot(
  
  plot_data_strata,
  
  aes(
    x = strata,
    y = prevalence
  )
  
) +
  
  # Observed species-level prevalence
  geom_jitter(
    
    width = 0.08,
    
    height = 0,
    
    alpha = 0.5,
    
    size = 2
  ) +
  
  # Predicted prevalence
  geom_point(
    
    data = newdata_strata,
    
    aes(
      x = strata,
      y = predicted
    ),
    
    inherit.aes = FALSE,
    
    size = 3
  ) +
  
  # 95% CI
  geom_errorbar(
    
    data = newdata_strata,
    
    aes(
      x = strata,
      ymin = CI_low,
      ymax = CI_high
    ),
    
    inherit.aes = FALSE,
    
    width = 0.08
  ) +
  
  scale_x_discrete(
    
    labels = c(
      
      "G" = "Ground",
      
      "S" = "Scansorial",
      
      "Ar" = "Arboreal"
    )
  ) +
  
  scale_y_continuous(
    
    limits = c(
      0,
      1
    ),
    
    labels = scales::percent
  ) +
  
  labs(
    
    x = "Foraging stratum",
    
    y = "Haemoplasma prevalence"
  ) +
  
  theme_classic() +
  
  theme(
    
    axis.text.x =
      element_text(
        size = 11
      ),
    
    axis.title =
      element_text(
        size = 12
      )
  )


# ============================================================
# 20. DISPLAY FIGURE
# ============================================================

print(
  p_strata
)


# ============================================================
# 21. SAVE FIGURE
# ============================================================

ggsave(
  
  "haemoplasma_prevalence_strata.png",
  
  p_strata,
  
  width = 6,
  
  height = 5,
  
  dpi = 300
)
```

# ============================================================
# 3. PREPARE DATASETS FOR EACH ACTIVITY VARIABLE
# ============================================================

data_nocturnal <- species_data_activity %>%
  filter(
    !is.na(activitynocturnal)
  )

data_crepuscular <- species_data_activity %>%
  filter(
    !is.na(activitycrepuscular)
  )

data_diurnal <- species_data_activity %>%
  filter(
    !is.na(activitydiurnal)
  )


cat(
  "\n================ SPECIES USED ================\n"
)

cat(
  "Nocturnal =",
  nrow(data_nocturnal),
  "\n"
)

cat(
  "Crepuscular =",
  nrow(data_crepuscular),
  "\n"
)

cat(
  "Diurnal =",
  nrow(data_diurnal),
  "\n"
)


# ============================================================
# 4. NULL MODELS
#
# IMPORTANT:
# Each null model is fitted to exactly the same species
# as its corresponding activity model.
# ============================================================

model_null_nocturnal <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ 1,

  data = data_nocturnal,

  family = betabinomial(
    link = "logit"
  )
)


model_null_crepuscular <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ 1,

  data = data_crepuscular,

  family = betabinomial(
    link = "logit"
  )
)


model_null_diurnal <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ 1,

  data = data_diurnal,

  family = betabinomial(
    link = "logit"
  )
)


# ============================================================
# 5. ACTIVITY MODELS
# ============================================================

model_nocturnal <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ activitynocturnal,

  data = data_nocturnal,

  family = betabinomial(
    link = "logit"
  )
)


model_crepuscular <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ activitycrepuscular,

  data = data_crepuscular,

  family = betabinomial(
    link = "logit"
  )
)


model_diurnal <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ activitydiurnal,

  data = data_diurnal,

  family = betabinomial(
    link = "logit"
  )
)


# ============================================================
# 6. LIKELIHOOD-RATIO TESTS
# ============================================================

lrt_nocturnal <- anova(
  model_null_nocturnal,
  model_nocturnal
)


lrt_crepuscular <- anova(
  model_null_crepuscular,
  model_crepuscular
)


lrt_diurnal <- anova(
  model_null_diurnal,
  model_diurnal
)


cat(
  "\n================ NOCTURNAL LRT ================\n"
)

print(
  lrt_nocturnal
)


cat(
  "\n================ CREPUSCULAR LRT ================\n"
)

print(
  lrt_crepuscular
)


cat(
  "\n================ DIURNAL LRT ================\n"
)

print(
  lrt_diurnal
)


# ============================================================
# 7. AIC
# ============================================================

cat(
  "\n================ AIC ================\n"
)

print(
  AIC(
    model_null_nocturnal,
    model_nocturnal
  )
)

print(
  AIC(
    model_null_crepuscular,
    model_crepuscular
  )
)

print(
  AIC(
    model_null_diurnal,
    model_diurnal
  )
)


# ============================================================
# 8. MODEL SUMMARIES
# ============================================================

cat(
  "\n================ NOCTURNAL MODEL ================\n"
)

print(
  summary(
    model_nocturnal
  )
)


cat(
  "\n================ CREPUSCULAR MODEL ================\n"
)

print(
  summary(
    model_crepuscular
  )
)


cat(
  "\n================ DIURNAL MODEL ================\n"
)

print(
  summary(
    model_diurnal
  )
)


# ============================================================
# 9. FUNCTION TO EXTRACT RESULTS
# ============================================================

extract_activity_results <- function(

  model,
  null_model,
  variable,
  lrt

) {


  coef_table <-
    summary(
      model
    )$coefficients$cond


  # ----------------------------------------------------------
  # Extract binary predictor
  # ----------------------------------------------------------

  predictor_row <-
    setdiff(
      rownames(coef_table),
      "(Intercept)"
    )[1]


  estimate <-
    coef_table[
      predictor_row,
      "Estimate"
    ]


  SE <-
    coef_table[
      predictor_row,
      "Std. Error"
    ]


  z <-
    coef_table[
      predictor_row,
      "z value"
    ]


  p <-
    coef_table[
      predictor_row,
      "Pr(>|z|)"
    ]


  # ----------------------------------------------------------
  # LRT
  # ----------------------------------------------------------

  LRT_chisq <-
    lrt$Chisq[2]


  LRT_p <-
    lrt$`Pr(>Chisq)`[2]


  # ----------------------------------------------------------
  # AIC
  # ----------------------------------------------------------

  AIC_null <-
    AIC(
      null_model
    )


  AIC_model <-
    AIC(
      model
    )


  delta_AIC <-
    AIC_null -
    AIC_model


  # ----------------------------------------------------------
  # 95% CI on log-odds scale
  # ----------------------------------------------------------

  CI_low <-
    estimate -
    1.96 * SE


  CI_high <-
    estimate +
    1.96 * SE


  # ----------------------------------------------------------
  # Odds ratio
  # ----------------------------------------------------------

  OR <-
    exp(
      estimate
    )


  OR_low <-
    exp(
      CI_low
    )


  OR_high <-
    exp(
      CI_high
    )


  # ----------------------------------------------------------
  # Return results
  # ----------------------------------------------------------

  data.frame(

    variable =
      variable,

    coefficient =
      predictor_row,

    n_species =
      nobs(model),

    estimate =
      estimate,

    SE =
      SE,

    z =
      z,

    p_coefficient =
      p,

    LRT_chisq =
      LRT_chisq,

    LRT_p =
      LRT_p,

    AIC_null =
      AIC_null,

    AIC_model =
      AIC_model,

    delta_AIC =
      delta_AIC,

    CI_low =
      CI_low,

    CI_high =
      CI_high,

    OR =
      OR,

    OR_low =
      OR_low,

    OR_high =
      OR_high

  )

}


# ============================================================
# 10. COMBINE RESULTS
# ============================================================

results_activity <-

  bind_rows(

    extract_activity_results(

      model_nocturnal,

      model_null_nocturnal,

      "activitynocturnal",

      lrt_nocturnal

    ),

    extract_activity_results(

      model_crepuscular,

      model_null_crepuscular,

      "activitycrepuscular",

      lrt_crepuscular

    ),

    extract_activity_results(

      model_diurnal,

      model_null_diurnal,

      "activitydiurnal",

      lrt_diurnal

    )

  )


# ============================================================
# 11. FDR CORRECTION
# ============================================================

# Three global LRTs

results_activity$LRT_p_FDR <-

  p.adjust(

    results_activity$LRT_p,

    method = "BH"

  )


# Three coefficient tests

results_activity$p_coefficient_FDR <-

  p.adjust(

    results_activity$p_coefficient,

    method = "BH"

  )


# ============================================================
# 12. FINAL RESULTS
# ============================================================

cat(
  "\n================ FINAL ACTIVITY RESULTS ================\n"
)

print(
  results_activity,
  row.names = FALSE
)

```
Results:
No association was detected between haemoplasma prevalence and foraging activity. Separate beta-binomial models showed no effect of nocturnal (OR = 0.79, 95% CI: 0.09–6.83, LRT p = 0.835), crepuscular (OR = 0.80, 95% CI: 0.23–2.74, LRT p = 0.720), or diurnal activity (OR = 0.80, 95% CI: 0.23–2.74, LRT p = 0.720). All associations remained non-significant after FDR correction (pFDR = 0.835), and activity models had higher AIC than null models.

Interpretation: Temporal foraging activity does not explain interspecific variation in haemoplasma prevalence.
```

# ============================================================
# SPECIES-LEVEL HAEMOPLASMA PREVALENCE
# ~ LIFE-HISTORY TRAITS
#
# Traits:
# bodymass       = Mean adult body mass (g)
# longevity      = Mean longevity (years)
# femalematurity = Mean age at female maturity (days)
# littersize     = Mean litter size (offspring/litter)
#
# Beta-binomial GLM + LRT + AIC + FDR
#
# Each trait is tested separately.
# Continuous predictors are standardized (1 SD increase).
# No (1 | species): one observation per species.
# ============================================================

library(dplyr)
library(glmmTMB)


# ============================================================
# 1. PREPARE SPECIES-LEVEL DATA
# ============================================================

species_data_lifehistory <- data_hemoplasma_stat %>%

  group_by(species) %>%

  summarise(

    n = n(),

    n_positive = sum(
      hemoplasma == 1,
      na.rm = TRUE
    ),

    prevalence = n_positive / n,

    .groups = "drop"

  ) %>%

  left_join(

    data_mammal_traits %>%

      select(
        species,
        bodymass,
        longevity,
        femalematurity,
        littersize
      ),

    by = "species"

  )


# ============================================================
# 2. STANDARDIZE CONTINUOUS TRAITS
# ============================================================

species_data_lifehistory <- species_data_lifehistory %>%

  mutate(

    bodymass_z =
      as.numeric(
        scale(bodymass)
      ),

    longevity_z =
      as.numeric(
        scale(longevity)
      ),

    femalematurity_z =
      as.numeric(
        scale(femalematurity)
      ),

    littersize_z =
      as.numeric(
        scale(littersize)
      )

  )


# ============================================================
# 3. CHECK DATA
# ============================================================

cat(
  "\n================ DATA CHECK ================\n"
)

cat(
  "Number of species =",
  nrow(species_data_lifehistory),
  "\n"
)

print(

  colSums(
    is.na(
      species_data_lifehistory[
        ,
        c(
          "bodymass",
          "longevity",
          "femalematurity",
          "littersize"
        )
      ]
    )
  )

)


# ============================================================
# 4. DATASETS FOR EACH TRAIT
# ============================================================

data_bodymass <- species_data_lifehistory %>%
  filter(!is.na(bodymass_z))

data_longevity <- species_data_lifehistory %>%
  filter(!is.na(longevity_z))

data_femalematurity <- species_data_lifehistory %>%
  filter(!is.na(femalematurity_z))

data_littersize <- species_data_lifehistory %>%
  filter(!is.na(littersize_z))


cat(
  "\n================ SPECIES USED ================\n"
)

cat(
  "Body mass =",
  nrow(data_bodymass),
  "\n"
)

cat(
  "Longevity =",
  nrow(data_longevity),
  "\n"
)

cat(
  "Female maturity =",
  nrow(data_femalematurity),
  "\n"
)

cat(
  "Litter size =",
  nrow(data_littersize),
  "\n"
)


# ============================================================
# 5. NULL MODELS
# ============================================================

model_null_bodymass <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ 1,

  data = data_bodymass,

  family = betabinomial(
    link = "logit"
  )

)


model_null_longevity <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ 1,

  data = data_longevity,

  family = betabinomial(
    link = "logit"
  )

)


model_null_femalematurity <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ 1,

  data = data_femalematurity,

  family = betabinomial(
    link = "logit"
  )

)


model_null_littersize <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ 1,

  data = data_littersize,

  family = betabinomial(
    link = "logit"
  )

)


# ============================================================
# 6. TRAIT MODELS
# ============================================================

model_bodymass <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ bodymass_z,

  data = data_bodymass,

  family = betabinomial(
    link = "logit"
  )

)


model_longevity <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ longevity_z,

  data = data_longevity,

  family = betabinomial(
    link = "logit"
  )

)


model_femalematurity <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ femalematurity_z,

  data = data_femalematurity,

  family = betabinomial(
    link = "logit"
  )

)


model_littersize <- glmmTMB(

  cbind(
    n_positive,
    n - n_positive
  ) ~ littersize_z,

  data = data_littersize,

  family = betabinomial(
    link = "logit"
  )

)


# ============================================================
# 7. LIKELIHOOD-RATIO TESTS
# ============================================================

lrt_bodymass <- anova(
  model_null_bodymass,
  model_bodymass
)

lrt_longevity <- anova(
  model_null_longevity,
  model_longevity
)

lrt_femalematurity <- anova(
  model_null_femalematurity,
  model_femalematurity
)

lrt_littersize <- anova(
  model_null_littersize,
  model_littersize
)


# ============================================================
# 8. FUNCTION TO EXTRACT RESULTS
# ============================================================

extract_lifehistory_results <- function(

  model,
  null_model,
  variable,
  lrt

) {

  coef_table <-
    summary(
      model
    )$coefficients$cond


  predictor_row <-
    setdiff(
      rownames(coef_table),
      "(Intercept)"
    )[1]


  estimate <-
    coef_table[
      predictor_row,
      "Estimate"
    ]


  SE <-
    coef_table[
      predictor_row,
      "Std. Error"
    ]


  z <-
    coef_table[
      predictor_row,
      "z value"
    ]


  p <-
    coef_table[
      predictor_row,
      "Pr(>|z|)"
    ]


  LRT_chisq <-
    lrt$Chisq[2]


  LRT_p <-
    lrt$`Pr(>Chisq)`[2]


  AIC_null <-
    AIC(
      null_model
    )


  AIC_model <-
    AIC(
      model
    )


  delta_AIC <-
    AIC_null -
    AIC_model


  CI_low <-
    estimate -
    1.96 * SE


  CI_high <-
    estimate +
    1.96 * SE


  data.frame(

    variable =
      variable,

    n_species =
      nobs(model),

    estimate =
      estimate,

    SE =
      SE,

    z =
      z,

    p_coefficient =
      p,

    LRT_chisq =
      LRT_chisq,

    LRT_p =
      LRT_p,

    AIC_null =
      AIC_null,

    AIC_model =
      AIC_model,

    delta_AIC =
      delta_AIC,

    CI_low =
      CI_low,

    CI_high =
      CI_high

  )

}


# ============================================================
# 9. COMBINE RESULTS
# ============================================================

results_lifehistory <-

  bind_rows(

    extract_lifehistory_results(
      model_bodymass,
      model_null_bodymass,
      "bodymass",
      lrt_bodymass
    ),

    extract_lifehistory_results(
      model_longevity,
      model_null_longevity,
      "longevity",
      lrt_longevity
    ),

    extract_lifehistory_results(
      model_femalematurity,
      model_null_femalematurity,
      "femalematurity",
      lrt_femalematurity
    ),

    extract_lifehistory_results(
      model_littersize,
      model_null_littersize,
      "littersize",
      lrt_littersize
    )

  )


# ============================================================
# 10. FDR CORRECTION
# ============================================================

results_lifehistory$LRT_p_FDR <-

  p.adjust(
    results_lifehistory$LRT_p,
    method = "BH"
  )


results_lifehistory$p_coefficient_FDR <-

  p.adjust(
    results_lifehistory$p_coefficient,
    method = "BH"
  )


# ============================================================
# 11. DISPLAY RESULTS
# ============================================================

cat(
  "\n================ LIFE-HISTORY RESULTS ================\n"
)

print(
  results_lifehistory,
  row.names = FALSE
)
```
Results : 
No significant association was detected between haemoplasma prevalence and any of the four life-history traits (all FDR-adjusted LRT p > 0.34).
- Body mass: no association (LRT p = 0.754; ΔAIC = −1.90).
- Longevity: positive but non-significant trend (β = 0.41; LRT p = 0.175; ΔAIC = −0.16).
- Female maturity: positive but non-significant trend (β = 0.50; LRT p = 0.159; ΔAIC ≈ 0).
- Litter size: no association (LRT p = 0.547; ΔAIC = −1.64).

Interpretation: 
Haemoplasma prevalence was not significantly related to body mass, longevity, age at female maturity, or litter size.
```































































## Step 4. Variation in hemoplasma infection across mammalian orders (GLMM model 1) 

### Contingency table
```
df_species <- data_hemoplasma_stat %>%
  group_by(species, order) %>%
  summarise(
    n = n(),
    n_infected = sum(as.numeric(as.character(hemoplasma)) == 1, na.rm = TRUE),
    infected = as.integer(n_infected > 0),
    .groups = "drop"
  )
df_order <- df_species %>%
  group_by(order) %>%
  summarise(
    infected_species = sum(infected),
    uninfected_species = n() - sum(infected),
    .groups = "drop"
  )
contingency_table <- as.matrix(df_order[, c("infected_species", "uninfected_species")])
rownames(contingency_table) <- df_order$order
contingency_table
```

Results :
```
Order                infected_species uninfected_species
Carnivora                      3                  3
Cingulata                      1                  1
Didelphimorphia                4                  4
Pilosa                         2                  2
Primates                       2                  3
Rodentia                       8                 11
```

### GLMM preparation
Sampling effort is included as a covariate (log-transformed number of individuals per species) :
```
data_hemoplasma_stat <- data_hemoplasma_stat %>%
  group_by(species) %>%
  mutate(
    n_sampled = n(),
    log_n = log(n_sampled)
  ) %>%
  ungroup() %>%
  mutate(hemoplasma = as.numeric(as.character(hemoplasma)))
```

### GLMM (Model 1)
This model tests whether `hemoplasma` infection probability varies among mammalian orders (`order`) while controlling for differences in sampling effort (`log_n`) and accounting for species-level random effects (`1 | species`).
```
mod1_full <- glmer(
  hemoplasma ~ order * log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)
```

### Model term significance testing
Model terms were evaluated using likelihood ratio tests via single-term deletions (drop1 function with Chi-square tests).
```
res1 <- drop1(mod1_full, test = "Chisq")
res1
```

Results :
```
Single term deletions
Model:
hemoplasma ~ order + log_n + (1 | species)
       npar    AIC     LRT Pr(Chi)  
<none>      414.76                  
order     5 416.94 12.1855 0.03233 *
log_n     1 416.19  3.4387 0.06369 .
```

### Interpretation
Hemoplasma infection probability varied significantly among mammalian orders (χ² test, _p_ = 0.032), indicating a non-random distribution of infection across host taxonomic groups.

A marginal effect of sampling effort (`log_n`) was also detected (_p_ = 0.064), suggesting a weak influence of species sampling intensity on observed prevalence.

### Model comparison with null and univariate models
```
mod1_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

mod1_order <- glmer(
  hemoplasma ~ order + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

mod1_log_n <- glmer(
  hemoplasma ~ log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

anova(mod1_null, mod1_order, test = "Chisq")
anova(mod1_null, mod1_log_n, test = "Chisq")

aics <- AIC(mod1_null, mod1_order, mod1_log_n)
aic_null <- aics["mod1_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
aics[, c("AIC", "delta_AIC_vs_null")]
```

Results : 
```
> anova(mod1_null, mod1_order, test="Chisq")
Data: data_hemoplasma_stat
Models:
mod1_null: hemoplasma ~ 1 + (1 | species)
mod1_order: hemoplasma ~ order + (1 | species)
           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
mod1_null     2 416.38 425.22 -206.19    412.38                       
mod1_order    7 416.19 447.13 -201.10    402.19 10.187  5    0.07011 .

> anova(mod1_null, mod1_log_n, test="Chisq")
Data: data_hemoplasma_stat
Models:
mod1_null: hemoplasma ~ 1 + (1 | species)
mod1_log_n: hemoplasma ~ log_n + (1 | species)
           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod1_null     2 416.38 425.22 -206.19    412.38                     
mod1_log_n    3 416.94 430.20 -205.47    410.94 1.4399  1     0.2302

                AIC delta_AIC_vs_null
mod1_null  416.3804         0.0000000
mod1_order 416.1937        -0.1866964
mod1_log_n 416.9406         0.5601158
```

### Interpretation (model comparison)
The inclusion of `order` slightly improved model fit compared to the null model, although this effect was not statistically significant (_p_ = 0.07), suggesting a weak signal of taxonomic structure in infection probability.

In contrast, `log_n` did not improve model fit (_p_ = 0.23), indicating no detectable effect of sampling effort.

AIC comparisons supported these results, with minimal differences between models (ΔAIC < 1), indicating no strong support for any predictor over the null model.

Overall, results suggest a weak but consistent tendency for variation in hemoplasma infection across mammalian orders.

### Post-hoc analysis of differences between mammalian orders (model-based pairwise comparisons)
We perform post-hoc comparisons to identify which orders differ in hemoplasma infection probability.
```
emm_order <- emmeans(mod1_full, pairwise ~ order, type = "response")
emm_order
```

Results:
```
$emmeans
 order             prob     SE  df asymp.LCL asymp.UCL
 Carnivora       0.7165 0.2960 Inf    0.1267     0.978
 Cingulata       0.2916 0.3570 Inf    0.0138     0.924
 Didelphimorphia 0.1668 0.1250 Inf    0.0334     0.537
 Pilosa          0.1183 0.1300 Inf    0.0115     0.607
 Primates        0.8773 0.1400 Inf    0.3576     0.989
 Rodentia        0.0518 0.0373 Inf    0.0121     0.195
Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 

$contrasts
 contrast                    odds.ratio       SE  df null z.ratio p.value
 Carnivora / Cingulata           6.1399  12.7000 Inf    1   0.878  0.9520
 Carnivora / Didelphimorphia    12.6271  19.6000 Inf    1   1.637  0.5741
 Carnivora / Pilosa             18.8321  38.3000 Inf    1   1.443  0.7005
 Carnivora / Primates            0.3535   0.6030 Inf    1  -0.610  0.9904
 Carnivora / Rodentia           46.2877  65.4000 Inf    1   2.715  0.0722
 Cingulata / Didelphimorphia     2.0566   3.8300 Inf    1   0.387  0.9989
 Cingulata / Pilosa              3.0671   6.6300 Inf    1   0.518  0.9955
 Cingulata / Primates            0.0576   0.1170 Inf    1  -1.402  0.7260
 Cingulata / Rodentia            7.5388  13.4000 Inf    1   1.138  0.8656
 Didelphimorphia / Pilosa        1.4914   2.3300 Inf    1   0.256  0.9999
 Didelphimorphia / Primates      0.0280   0.0410 Inf    1  -2.441  0.1422
 Didelphimorphia / Rodentia      3.6658   3.8600 Inf    1   1.233  0.8206
 Pilosa / Primates               0.0188   0.0346 Inf    1  -2.158  0.2577
 Pilosa / Rodentia               2.4579   3.6800 Inf    1   0.600  0.9911
 Primates / Rodentia           130.9573 177.0000 Inf    1   3.616  0.0041
P value adjustment: tukey method for comparing a family of 6 estimates 
Tests are performed on the log odds ratio scale 
```

### Interpretation
Post-hoc pairwise comparisons showed variation in hemoplasma infection probability among mammalian orders, but most contrasts were not significant after Tukey correction.

A significant difference was found between Primates and Rodentia (_p_ = 0.0041), suggesting higher infection probabilities in Primates compared to Rodents, while all other comparisons were non-significant.

### Create a plot of hemoplasma prevalence by species and mammalian order
```
df_species <- data_hemoplasma_stat %>%
  group_by(species, order) %>%
  summarise(
    n = n(),
    n_infected = sum(as.numeric(as.character(hemoplasma)) == 1, na.rm = TRUE),
    prevalence = n_infected / n,
    .groups = "drop"
  )
mod_glmm <- glmer(
  hemoplasma ~ order + log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 1e5))
)
emm_prob <- emmeans(mod_glmm, ~ order, type = "response")
prob_df <- as.data.frame(emm_prob)
prob_df <- prob_df %>%
  arrange(prob) %>%
  mutate(order = factor(order, levels = order))
df_species <- df_species %>%
  mutate(order = factor(order, levels = levels(prob_df$order)))
order_colors <- c(
  "Primates" = "#0072B2",
  "Pilosa" = "#E69F00",
  "Cingulata" = "#009E73",
  "Rodentia" = "#D55E00",
  "Carnivora" = "#CC79A7",
  "Didelphimorphia" = "#F0E442"
)
set.seed(1)
p <- ggplot() +
  geom_jitter(
    data = df_species,
    aes(x = order, y = prevalence, color = order, size = n),
    width = 0.35,
    height = 0,
    alpha = 0.5
  ) +
  geom_segment(
    data = prob_df,
    aes(
      x = as.numeric(order) - 0.25,
      xend = as.numeric(order) + 0.25,
      y = asymp.LCL,
      yend = asymp.LCL,
      color = order
    ),
    linewidth = 1
  ) +
  geom_segment(
    data = prob_df,
    aes(
      x = as.numeric(order) - 0.25,
      xend = as.numeric(order) + 0.25,
      y = prob,
      yend = prob,
      color = order
    ),
    linewidth = 1.5
  ) +
  geom_segment(
    data = prob_df,
    aes(
      x = as.numeric(order) - 0.25,
      xend = as.numeric(order) + 0.25,
      y = asymp.UCL,
      yend = asymp.UCL,
      color = order
    ),
    linewidth = 1
  ) +
  scale_color_manual(values = order_colors) +
  scale_size_continuous(range = c(2, 10), name = "Sample size") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Mammalian order",
    y = "Hemoplasma prevalence"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )
print(p)
ggsave(
  filename = "Fig_2_Hemoplasma_prevalence_by_order.pdf",
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)
```

## Step 5. Variation in hemoplasma infection across sex within species where infection was detected (GLMM model 2)

### Data preparation
We restricted here the analysis to mammalian species with at least one infected individual, in order to test sex effects within relevant host species.
```
data_sex <- data_hemoplasma_stat[
  complete.cases(data_hemoplasma_stat[, c("hemoplasma", "sex", "species")]),
]

species_infected_sex <- data_sex %>%
  group_by(species) %>%
  summarise(infected = any(hemoplasma == 1, na.rm = TRUE)) %>%
  filter(infected) %>%
  pull(species)

data_inf_sex <- data_sex %>%
  filter(species %in% species_infected_sex) %>%
  mutate(hemoplasma = as.numeric(as.character(hemoplasma))) %>%
  group_by(species) %>%
  mutate(
    n_sampled = n(),
    log_n = log(n_sampled)
  ) %>%
  ungroup()
```

### GLMM (Model 2)
This model tests whether hemoplasma infection probability differs between sexes (`sex`) while controlling for sampling effort (`log_n`) and accounting for species-level random effects (`1 | species`).
```
mod2_full <- glmer(
  hemoplasma ~ sex + log_n + (1 | species),
  family = binomial,
  data = data_inf_sex,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)
```

### Model term significance testing
Model terms were evaluated using likelihood ratio tests via single-term deletions (drop1 function with Chi-square tests).
```
res2 <- drop1(mod2_full, test = "Chisq")
res2
```

Results :
```
Single term deletions
Model:
hemoplasma ~ sex + log_n + (1 | species)
       npar    AIC     LRT Pr(Chi)  
<none>      279.07                
sex       1 279.59 2.52007  0.1124
log_n     1 277.51 0.44085  0.5067
```

### Interpretation
No significant effect of `sex` on `hemoplasma` infection probability was detected (χ² = 2.52, _p_ = 0.11).

Sampling effort (`log_n`) had no detectable effect (p = 0.51).

### Model comparison with null and univariate models
```
mod2_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  family = binomial,
  data = data_inf_sex,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

mod2_sex <- glmer(
  hemoplasma ~ sex + (1 | species),
  family = binomial,
  data = data_inf_sex,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

mod2_log_n <- glmer(
  hemoplasma ~ log_n + (1 | species),
  family = binomial,
  data = data_inf_sex,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

anova(mod2_null, mod2_sex, test = "Chisq")
anova(mod2_null, mod2_log_n, test = "Chisq")

aics <- AIC(mod2_null, mod2_sex, mod2_log_n)
aic_null <- aics["mod2_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
aics[, c("AIC", "delta_AIC_vs_null")]
```

Results :
```
> anova(mod2_null, mod2_sex, test = "Chisq")
Data: data_inf_sex
Models:
mod2_null: hemoplasma ~ 1 + (1 | species)
mod2_sex: hemoplasma ~ sex + (1 | species)
          npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod2_null    2 278.07 285.65 -137.03    274.07                     
mod2_sex     3 277.51 288.89 -135.75    271.51 2.5598  1     0.1096

> anova(mod2_null, mod2_log_n, test = "Chisq")
Data: data_inf_sex
Models:
mod2_null: hemoplasma ~ 1 + (1 | species)
mod2_log_n: hemoplasma ~ log_n + (1 | species)
           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod2_null     2 278.07 285.65 -137.03    274.07                     
mod2_log_n    3 279.59 290.97 -136.79    273.59 0.4806  1     0.4881

                AIC delta_AIC_vs_null
mod2_null  278.0684         0.0000000
mod2_sex   277.5086        -0.5598415
mod2_log_n 279.5878         1.5193843
```

### Interpretation (model comparison)
The inclusion of `sex` or `log_n` did not improve model fit compared to the null model (_p_ = 0.11 and 0.49, respectively).

AIC comparison supported these results, with only a minimal improvement for the model including `sex` (ΔAIC < 1), indicating weak support for this predictor. Overall, results provide no clear evidence for sex-biased hemoplasma infection.

### Post-hoc analysis of differences between sex in infected species (model-based pairwise comparisons)

We perform post-hoc comparisons to identify if sexes differ in hemoplasma infection probability.

```
emm_sex <- emmeans(mod2_full, pairwise ~ sex, type = "response")
emm_sex
```

Results:
```
$emmeans
 sex  prob    SE  df asymp.LCL asymp.UCL
 F   0.371 0.202 Inf    0.0972     0.764
 M   0.506 0.214 Inf    0.1603     0.846
Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 

$contrasts
 contrast odds.ratio    SE  df null z.ratio p.value
 F / M         0.576 0.199 Inf    1  -1.594  0.1109
```

### Interpretation
Predicted infection probabilities were slightly higher in males than in females, but this difference was not statistically significant (_p_ = 0.11).

### Species-level sex bias in hemoplasma infection (Fisher exact tests)
To test whether `hemoplasma` infection differs between males and females within host `species`, we restricted analyses to infected `species` that had at least 20 sexed individuals to ensure minimum statistical power.


### Data preparation
```
data_clean <- data_hemoplasma_stat %>%
  mutate(hemoplasma = as.numeric(as.character(hemoplasma))) %>%
  filter(!is.na(sex))

species_keep <- data_clean %>%
  group_by(species) %>%
  summarise(
    n_sexed = n(),
    any_infected = any(hemoplasma == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_sexed >= 20, any_infected) %>%
  pull(species)

write.csv(species_keep,
          "species_used_fisher_sex_bias.csv",
          row.names = FALSE)

species_sex_summary <- data_clean %>%
  filter(species %in% species_keep) %>%
  group_by(species, sex) %>%
  summarise(
    n_total = n(),
    n_infected = sum(hemoplasma == 1, na.rm = TRUE),
    n_uninfected = sum(hemoplasma == 0, na.rm = TRUE),
    prevalence = n_infected / n_total,
    .groups = "drop"
  ) %>%
  arrange(species, sex)
species_sex_summary %>%
  arrange(species, sex)
```

Results for infected `species` that had at least 20 sexed individuals :
```
# A tibble: 8 × 6
  species               sex   n_total n_infected n_uninfected prevalence
  <fct>                 <fct>   <int>      <int>        <int>      <dbl>
1 Bradypus_tridactylus  F          43          2           41     0.0465
2 Bradypus_tridactylus  M          49          2           47     0.0408
3 Choloepus_didactylus  F          49         39           10     0.796 
4 Choloepus_didactylus  M          34         29            5     0.853 
5 Didelphis_marsupialis F          19          8           10     0.421 
6 Didelphis_marsupialis M          19         13            5     0.684 
7 Saguinus_midas        F          15         15            0     1     
8 Saguinus_midas        M          23         23            0     1     
```

### Fisher exact tests per species
```
fisher_results <- data_fisher %>%
  group_by(species) %>%
  summarise(
    tab = list(table(sex, hemoplasma)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    fisher = list(
      if (all(dim(tab) == c(2, 2))) {
        fisher.test(tab)
      } else {
        NULL
      }
    ),
    p_value = if (!is.null(fisher)) fisher$p.value else NA_real_,
    odds_ratio = if (!is.null(fisher)) as.numeric(fisher$estimate) else NA_real_
  ) %>%
  ungroup() %>%
  select(species, p_value, odds_ratio)
```

### Multiple testing correction (FDR)
```
fisher_results <- fisher_results %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr")
  ) %>%
  arrange(p_adj)

fisher_results
```

Results : 
```
# A tibble: 4 × 4
  species               p_value odds_ratio  p_adj
  <fct>                   <dbl>      <dbl>  <dbl>
1 Didelphis_marsupialis   0.176      3.14   0.527
2 Choloepus_didactylus    0.573      1.48   0.860
3 Bradypus_tridactylus    1          0.874  1    
4 Saguinus_midas         NA         NA     NA    
```

### Interpretation
Fisher exact tests performed separately for each species (restricted to species with ≥20 sexed individuals and at least one infected individual) revealed no significant sex differences in hemoplasma infection after FDR correction. For _Saguinus midas_, the test could not be computed due to insufficient data structure in the contingency table (ie, all individuals are infected).

### Create a plot of hemoplasma prevalence by infected species
```
species_infected <- data_hemoplasma_stat %>%
  mutate(hemoplasma = as.numeric(as.character(hemoplasma))) %>%
  group_by(species) %>%
  summarise(any_infected = any(hemoplasma == 1, na.rm = TRUE)) %>%
  filter(any_infected) %>%
  pull(species)

data_clean <- data_hemoplasma_stat %>%
  filter(
    species %in% species_infected,
    !is.na(sex)
  ) %>%
  mutate(hemoplasma = as.numeric(as.character(hemoplasma)))

species_keep <- data_clean %>%
  group_by(species) %>%
  summarise(n_sexed = n()) %>%
  filter(n_sexed >= 20) %>%
  pull(species)

df_plot <- data_clean %>%
  filter(species %in% species_keep) %>%
  group_by(species, sex) %>%
  summarise(
    n = n(),
    n_infected = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_infected / n,
    .groups = "drop"
  ) %>%
  mutate(
    species_sex = paste(species, sex)
  )

species_order <- unique(df_plot$species)

levels_ordered <- unlist(lapply(species_order, function(sp) {
  c(paste(sp, "M"),
    paste(sp, "F"),
    paste0("gap_", sp))
}))

df_plot$species_sex <- factor(df_plot$species_sex,
                              levels = levels_ordered)

species_cols <- setNames(hue_pal()(length(species_order)), species_order)

fill_cols <- c()
for (sp in species_order) {
  fill_cols[paste(sp, "M")] <- species_cols[sp]
  fill_cols[paste(sp, "F")] <- alpha(species_cols[sp], 0.5)
}

df_plot_plot <- df_plot %>%
  mutate(
    species_sex_plot = ifelse(grepl("gap_", species_sex), NA, species_sex)
  )

p_sex_species <- ggplot(df_plot_plot,
                        aes(x = species_sex_plot,
                            y = prevalence,
                            fill = species_sex)) +
  
  geom_bar(stat = "identity",
           color = "black",
           na.rm = TRUE) +
  
  scale_fill_manual(values = fill_cols,
                    na.translate = FALSE) +
  
  scale_x_discrete(drop = FALSE) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank()
  ) +
  
  labs(
    x = "Species / Sex",
    y = "Hemoplasma infection prevalence",
    fill = "Species & Sex"
  )

p_sex_species

pdf("Fig_1C_hemoplasma_sex_species.pdf", width = 12, height = 5)
print(p_sex_species)
dev.off()
```

## Step 6. Variation in hemoplasma infection across host ecological traits (GLMM Model 3) 

### Data preparation
```
data_hemoplasma_stat <- data_hemoplasma_stat %>%
  group_by(species) %>%
  mutate(
    n_sampled = n(),
    log_n = log(n_sampled)
  ) %>%
  ungroup()
```

### GLMM (Model 3)

This model tests whether hemoplasma infection probability varies according to host ecological traits (`strata`, `activity`, `diet`, `sociality`), while controlling for sampling effort (`log_n`) and accounting for species-level random effects (`1 | species`).
```
mod3_full <- glmer(
  hemoplasma ~ strata + activity + diet + sociality + log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)
```

### Model term significance testing
Model terms were evaluated using likelihood ratio tests via single-term deletions (drop1 function with Chi-square tests).
```
res3 <- drop1(mod3_full, test = "Chisq")
res3
```

Results :
```
Single term deletions
Model:
hemoplasma ~ strata + activity + diet + sociality + log_n + (1 | species)
                 npar    AIC    LRT Pr(Chi)  
<none>                421.76                 
strata    2 418.86 1.0977 0.57760  
activity            1 423.09 3.3279 0.06812 .
diet                3 416.11 0.3554 0.94930  
sociality           1 419.80 0.0445 0.83297
log_n               1 423.37 3.6160 0.05723 .  
```

### Interpretation
No significant effect of host ecological traits was detected on hemoplasma infection probability. Activity and sampling effort showed weak non-significant trends but did not reach statistical significance.


### Model comparison with null and univariate models
```

mod3_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod3_strata <- glmer(
  hemoplasma ~ strata + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod3_activity <- glmer(
  hemoplasma ~ activity + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod3_diet <- glmer(
  hemoplasma ~ diet + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod3_sociality <- glmer(
  hemoplasma ~ sociality + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod3_log_n <- glmer(
  hemoplasma ~ log_n + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

anova(mod3_null, mod3_strata, test = "Chisq")
anova(mod3_null, mod3_activity, test = "Chisq")
anova(mod3_null, mod3_diet, test = "Chisq")
anova(mod3_null, mod3_sociality, test = "Chisq")
anova(mod3_null, mod3_log_n, test = "Chisq")

aics <- AIC(mod3_null,
            mod3_strata,
            mod3_activity,
            mod3_diet,
            mod3_sociality,
            mod3_log_n)

aic_null <- aics["mod3_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
aics[, c("AIC", "delta_AIC_vs_null")]
```

Results : 
```
> anova(mod3_null, mod3_strata, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod3_null: hemoplasma ~ 1 + (1 | species)
mod3_strata: hemoplasma ~ strata + (1 | species)
                      npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod3_null                2 416.38 425.22 -206.19    412.38                     
mod3_strata    4 419.55 437.23 -205.78    411.55 0.8301  2     0.6603

> anova(mod3_null, mod3_activity, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod3_null: hemoplasma ~ 1 + (1 | species)
mod3_activity: hemoplasma ~ activity + (1 | species)
              npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
mod3_null        2 416.38 425.22 -206.19    412.38                       
mod3_activity    3 413.14 426.40 -203.57    407.14 5.2439  1    0.02202 *

> anova(mod3_null, mod3_diet, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod3_null: hemoplasma ~ 1 + (1 | species)
mod3_diet: hemoplasma ~ diet + (1 | species)
          npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod3_null    2 416.38 425.22 -206.19    412.38                     
mod3_diet    5 419.44 441.54 -204.72    409.44 2.9409  3     0.4008

> anova(mod3_null, mod3_sociality, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod3_null: hemoplasma ~ 1 + (1 | species)
mod3_sociality: hemoplasma ~ sociality + (1 | species)
               npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod3_null         2 416.38 425.22 -206.19    412.38                     
mod3_sociality    3 417.58 430.84 -205.79    411.58 0.8027  1     0.3703

> anova(mod3_null, mod3_log_n, test = "Chisq")
Data: data_hemoplasma_stat
Models:
mod3_null: hemoplasma ~ 1 + (1 | species)
mod3_log_n: hemoplasma ~ log_n + (1 | species)
           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod3_null     2 416.38 425.22 -206.19    412.38                     
mod3_log_n    3 416.94 430.20 -205.47    410.94 1.4399  1     0.2302

                           AIC delta_AIC_vs_null
mod3_null             416.3804         0.0000000
mod3_strata 419.5503         3.1698637
mod3_activity         413.1365        -3.2439419
mod3_diet             419.4395         3.0590596
mod3_sociality        417.5777         1.1973107
mod3_log_n            416.9406         0.5601158
```

### Interpretation
Activity showed a marginal effect in the univariate model comparison (_p_ = 0.022), but this pattern was not robust in the full multivariate model (see above).

### Marginal means (estimated infection probabilities)
We estimated marginal means (back-transformed to response scale) for each level of host ecological traits.

```
emm_strata <- emmeans(mod3_full, ~ strata, type = "response")
emm_activity <- emmeans(mod3_full, ~ activity, type = "response")
emm_diet <- emmeans(mod3_full, ~ diet, type = "response")
emm_sociality <- emmeans(mod3_full, ~ sociality, type = "response")

emm_strata
emm_activity
emm_diet
emm_sociality
```

Results :
```
 strata   prob    SE  df asymp.LCL asymp.UCL
 Canopy           0.4112 0.298 Inf   0.05898     0.886
 Ground           0.4723 0.261 Inf   0.10315     0.874
 Mixed            0.0876 0.177 Inf   0.00126     0.879

 activity    prob     SE  df asymp.LCL asymp.UCL
 Diurnal   0.6517 0.3250 Inf   0.10161     0.969
 Nocturnal 0.0757 0.0919 Inf   0.00621     0.518

 diet         prob    SE  df asymp.LCL asymp.UCL
 Carnivore   0.469 0.583 Inf   0.00888     0.989
 Insectivore 0.199 0.317 Inf   0.00498     0.925
 Omnivore    0.212 0.172 Inf   0.03459     0.670
 Phytophage  0.284 0.272 Inf   0.02822     0.845

 sociality  prob    SE  df asymp.LCL asymp.UCL
 Group     0.312 0.307 Inf    0.0268     0.882
 Solitary  0.252 0.226 Inf    0.0313     0.779

Results are averaged over the levels of: strata, activity, diet 
Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 
```

### Post-hoc analysis of differences between levelS of host ecological traits (model-based pairwise comparisons)
```
pairs(emmeans(mod3_full, ~ strata), adjust = "tukey")
pairs(emmeans(mod3_full, ~ activity), adjust = "tukey")
pairs(emmeans(mod3_full, ~ diet), adjust = "tukey")
pairs(emmeans(mod3_full, ~ sociality), adjust = "tukey")
```

Results:
```
> pairs(emmeans(mod3_full, ~ strata), adjust = "tukey")
 contrast        estimate   SE  df z.ratio p.value
 Canopy - Ground   -0.248 1.13 Inf  -0.219  0.9739
 Canopy - Mixed     1.984 2.28 Inf   0.869  0.6596
 Ground - Mixed     2.233 2.20 Inf   1.017  0.5661

> pairs(emmeans(mod3_full, ~ activity), adjust = "tukey")
 contrast            estimate   SE  df z.ratio p.value
 Diurnal - Nocturnal     3.13 1.61 Inf   1.946  0.0517

> pairs(emmeans(mod3_full, ~ diet), adjust = "tukey")
 contrast                 estimate   SE  df z.ratio p.value
 Carnivore - Insectivore    1.2706 2.79 Inf   0.455  0.9687
 Carnivore - Omnivore       1.1872 2.42 Inf   0.490  0.9614
 Carnivore - Phytophage     0.7991 2.63 Inf   0.304  0.9903
 Insectivore - Omnivore    -0.0834 1.76 Inf  -0.047  1.0000
 Insectivore - Phytophage  -0.4716 2.06 Inf  -0.229  0.9958
 Omnivore - Phytophage     -0.3881 1.18 Inf  -0.328  0.9878

> pairs(emmeans(mod3_full, ~ sociality), adjust = "tukey")
 contrast         estimate   SE  df z.ratio p.value
 Group - Solitary    0.298 1.41 Inf   0.211  0.8328
```

### Interpretation
Marginal means showed broadly overlapping confidence intervals across all ecological trait categories. Tukey-adjusted pairwise comparisons revealed no significant differences among levels of vertical stratum, diet, or sociality, and only a marginal diurnal–nocturnal contrast (_p_ = 0.052).

Create a plot of hemoplasma prevalence by species ecological traits :
```
df_species <- data_hemoplasma_stat %>%
  group_by(species, strata, activity, diet, sociality) %>%
  summarise(
    n = n(),
    n_infected = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_infected / n,
    .groups = "drop"
  ) %>%
  mutate(
    prevalence = ifelse(is.nan(prevalence), 0, prevalence)
  )

fix_emm <- function(df, var){

  df <- as.data.frame(df)

  df[[var]] <- as.character(df[[var]])

  df$prob  <- as.numeric(df$prob)
  df$lower <- as.numeric(df$asymp.LCL)
  df$upper <- as.numeric(df$asymp.UCL)

  df
}

prob_stratum  <- fix_emm(prob_stratum, "strata")
prob_activity  <- fix_emm(prob_activity, "activity")
prob_diet      <- fix_emm(prob_diet, "diet")
prob_social    <- fix_emm(prob_social, "sociality")

make_panel <- function(df, pred, var, palette, title){

  df[[var]] <- factor(df[[var]])
  pred[[var]] <- factor(pred[[var]])

  ggplot() +

    geom_jitter(
      data = df,
      aes(x = .data[[var]], y = prevalence,
          color = .data[[var]], size = n),
      width = 0.20,
      alpha = 0.65,
      na.rm = TRUE
    ) +

    geom_segment(
      data = pred,
      aes(
        x = as.numeric(.data[[var]]) - 0.15,
        xend = as.numeric(.data[[var]]) + 0.15,
        y = lower,
        yend = lower,
        color = .data[[var]]
      ),
      linewidth = 0.9,
      na.rm = TRUE
    ) +

    geom_segment(
      data = pred,
      aes(
        x = as.numeric(.data[[var]]) - 0.15,
        xend = as.numeric(.data[[var]]) + 0.15,
        y = prob,
        yend = prob,
        color = .data[[var]]
      ),
      linewidth = 1.3,
      na.rm = TRUE
    ) +

    geom_segment(
      data = pred,
      aes(
        x = as.numeric(.data[[var]]) - 0.15,
        xend = as.numeric(.data[[var]]) + 0.15,
        y = upper,
        yend = upper,
        color = .data[[var]]
      ),
      linewidth = 0.9,
      na.rm = TRUE
    ) +

    scale_color_manual(values = palette, drop = FALSE) +
    scale_size_continuous(range = c(2, 9), name = "Sample size") +

    scale_y_continuous(
      labels = percent_format(accuracy = 1)
    ) +

    labs(
      x = NULL,
      y = "Hemoplasma prevalence",
      title = title
    ) +

    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}
pA <- make_panel(df_species, prob_stratum, "strata", stratum_colors, "A")
pB <- make_panel(df_species, prob_activity, "activity", activity_colors, "B")
pC <- make_panel(df_species, prob_diet, "diet", diet_colors, "C")
pD <- make_panel(df_species, prob_social, "sociality", social_colors, "D")
figure1 <- (pA | pB) / (pC | pD)
figure1
pdf(file = file.path(getwd(), "Figure2_Hemoplasma_ecological_traits.pdf"),
    width = 12, height = 9, useDingbats = FALSE)
print(figure1)
dev.off()
```

## Step 7. Association between hemoplasma infection and co-infecting blood parasites (GLMM Model 4)
### Data preparation
We tested whether `hemoplasma` infection probability is associated with co-infection by other blood parasites, including `anaplasmataceae`, `apicomplexa`, `trypanosoma`, and `filaria`. Species-level sampling effort was accounted for using `log_n`, and species identity was included as a random effect.
```
data_hemoplasma_stat <- data_hemoplasma_stat %>%
  group_by(species) %>%
  mutate(
    n_sampled = n(),
    log_n = log(n_sampled)
  ) %>%
  ungroup()

data_coinf <- data_hemoplasma_stat %>%
  filter(
    !is.na(hemoplasma),
    !is.na(anaplasmataceae),
    !is.na(apicomplexa),
    !is.na(trypanosoma),
    !is.na(filaria)
  )
```

## GLMM (Model 4)
This model tests whether `hemoplasma` infection probability is associated with co-infection by other blood parasites.
```
mod4_full <- glmer(
  hemoplasma ~ anaplasmataceae + apicomplexa + trypanosoma + filaria +
    log_n + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)
```

## Model term significance testing

Model terms were evaluated using likelihood ratio tests via single-term deletions.
```
res4 <- drop1(mod4_full, test = "Chisq")
res4
```

Results :
```
Single term deletions
Model:
hemoplasma ~ anaplasmataceae + apicomplexa + trypanosoma + filaria + 
    log_n + (1 | species)
                npar    AIC    LRT  Pr(Chi)   
<none>               120.13                   
anaplasmataceae    1 122.84 4.7040 0.030093 * 
apicomplexa        1 118.14 0.0109 0.916691   
trypanosoma        1 120.53 2.3976 0.121525   
filaria            1 118.14 0.0042 0.948289   
log_n              1 127.77 9.6327 0.001911 **
```

## Interpretation
`anaplasmataceae` infection was significantly associated with `hemoplasma` infection probability (χ² = 4.70, _p_ = 0.030), with higher infection probability in co-infected hosts. Sampling effort also had a strong effect (χ² = 9.63, _p_ = 0.0019), indicating lower detection probability with increasing sample size. No significant effects were detected for `apicomplexa`, `trypanosoma`, or `filaria` (all _p_ > 0.10).

## Model comparison with null and univariate models
```
mod4_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

mod4_anaplasma <- glmer(
  hemoplasma ~ anaplasmataceae + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

mod4_apicomplexa <- glmer(
  hemoplasma ~ apicomplexa + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

mod4_trypanosoma <- glmer(
  hemoplasma ~ trypanosoma + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

mod4_filaria <- glmer(
  hemoplasma ~ filaria + (1 | species),
  family = binomial,
  data = data_coinf,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

anova(mod4_null, mod4_anaplasma, test = "Chisq")
anova(mod4_null, mod4_apicomplexa, test = "Chisq")
anova(mod4_null, mod4_trypanosoma, test = "Chisq")
anova(mod4_null, mod4_filaria, test = "Chisq")

aics <- AIC(mod4_null,
            mod4_anaplasma,
            mod4_apicomplexa,
            mod4_trypanosoma,
            mod4_filaria)

aic_null <- aics["mod4_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
aics[, c("AIC", "delta_AIC_vs_null")]
```

Results : 
```
> anova(mod4_null, mod4_anaplasma, test = "Chisq")
Data: data_coinf
Models:
mod4_null: hemoplasma ~ 1 + (1 | species)
mod4_anaplasma: hemoplasma ~ anaplasmataceae + (1 | species)
               npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
mod4_null         2 125.55 131.92 -60.777    121.55                       
mod4_anaplasma    3 124.09 133.63 -59.044    118.09 3.4647  1    0.06269 .

> anova(mod4_null, mod4_apicomplexa, test = "Chisq")
Data: data_coinf
Models:
mod4_null: hemoplasma ~ 1 + (1 | species)
mod4_apicomplexa: hemoplasma ~ apicomplexa + (1 | species)
                 npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod4_null           2 125.55 131.92 -60.777    121.55                     
mod4_apicomplexa    3 127.47 137.02 -60.736    121.47 0.0816  1     0.7751

> anova(mod4_null, mod4_trypanosoma, test = "Chisq")
Data: data_coinf
Models:
mod4_null: hemoplasma ~ 1 + (1 | species)
mod4_trypanosoma: hemoplasma ~ trypanosoma + (1 | species)
                 npar    AIC    BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)
mod4_null           2 125.55 131.92 -60.777    121.55                    
mod4_trypanosoma    3 126.22 135.77 -60.112    120.22 1.329  1      0.249

> anova(mod4_null, mod4_filaria, test = "Chisq")
Data: data_coinf
Models:
mod4_null: hemoplasma ~ 1 + (1 | species)
mod4_filaria: hemoplasma ~ filaria + (1 | species)
             npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod4_null       2 125.55 131.92 -60.777    121.55                     
mod4_filaria    3 127.55 137.10 -60.776    121.55 0.0018  1     0.9662

                      AIC delta_AIC_vs_null
mod4_null        125.5536         0.0000000
mod4_anaplasma   124.0889        -1.4646561
mod4_apicomplexa 127.4719         1.9183628
mod4_trypanosoma 126.2246         0.6710417
mod4_filaria     127.5518         1.9981994
```

## Interpretation
Model comparisons confirmed a weak but consistent support for an association between `hemoplasma` infection and `anaplasmataceae` co-infection (ΔAIC ≈ −1.46), whereas no improvement in model fit was observed for `apicomplexa`, `trypanosoma`, or `filaria`.

## Post-hoc analysis of differences between blood parasites (model-based pairwise comparisons)
```
emm_anaplasma <- emmeans(mod4_full, ~ anaplasmataceae, type = "response")
emm_apicomplexa <- emmeans(mod4_full, ~ apicomplexa, type = "response")
emm_trypanosoma <- emmeans(mod4_full, ~ trypanosoma, type = "response")
emm_filaria <- emmeans(mod4_full, ~ filaria, type = "response")

emm_anaplasma
pairs(emmeans(mod4_full, ~ anaplasmataceae), adjust = "tukey")
emm_apicomplexa
pairs(emmeans(mod4_full, ~ apicomplexa), adjust = "tukey")
emm_trypanosoma
pairs(emmeans(mod4_full, ~ trypanosoma), adjust = "tukey")
emm_filaria
pairs(emmeans(mod4_full, ~ filaria), adjust = "tukey")
```

Results : 
```
 anaplasmataceae   prob     SE  df asymp.LCL asymp.UCL
 0               0.0635 0.0688 Inf   0.00695     0.396
 1               0.1878 0.1480 Inf   0.03345     0.607

 contrast                            estimate   SE  df z.ratio p.value
 anaplasmataceae0 - anaplasmataceae1    -1.23 0.61 Inf  -2.013  0.0442

 apicomplexa  prob     SE  df asymp.LCL asymp.UCL
 0           0.117 0.0844 Inf    0.0263     0.396
 1           0.105 0.1370 Inf    0.0068     0.669

 contrast                    estimate   SE  df z.ratio p.value
 apicomplexa0 - apicomplexa1    0.123 1.16 Inf   0.106  0.9158

 trypanosoma   prob     SE  df asymp.LCL asymp.UCL
 0           0.2947 0.1350 Inf   0.10506     0.598
 1           0.0361 0.0583 Inf   0.00141     0.499

 contrast                    estimate  SE  df z.ratio p.value
 trypanosoma0 - trypanosoma1     2.41 1.5 Inf   1.608  0.1079

 filaria  prob    SE  df asymp.LCL asymp.UCL
 0       0.113 0.100 Inf    0.0178     0.475
 1       0.109 0.111 Inf    0.0128     0.537

 contrast            estimate    SE  df z.ratio p.value
 filaria0 - filaria1   0.0426 0.657 Inf   0.065  0.9482

Results are averaged over the levels of: anaplasmataceae, apicomplexa, trypanosoma, log_n 
Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 
```

## Interpretation
`hemoplasma` infection probability was higher in hosts infected by `anaplasmataceae` (0.19 vs 0.06; _p_ = 0.044), while no significant differences were detected for `apicomplexa`, `trypanosoma`, or `filaria` (all _p _> 0.10), with strongly overlapping confidence intervals across groups.

## Create a plot of hemoplasma prevalence by co-infecting blood parasites
A FAIRE !!!



yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy





Class size effects on hemoplasma infection prevalence (GLMM with species random effect, model #3) : 
```
df_body <- data %>%
  group_by(species, body_size) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )
# Full model (body size only)
mod_body <- glmer(
  hemoplasma ~ body_size + (1 | species),
  data = data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

# Null model (only random effect)
mod_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  data = data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

# Summary
summary(mod_body)

# LRT (global test of body_size effect)
anova(mod_null, mod_body, test = "Chisq")

# AIC comparison
AIC_table <- data.frame(
  model = c("body_size + species", "null (species only)"),
  AIC = c(AIC(mod_body), AIC(mod_null))
)
AIC_table$delta_AIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table


# Post-hoc (if needed)
library(emmeans)
emmeans(mod_body, ~ body_size, type = "response")

```

Results are : 
```
Data: data
Models:
mod_null: hemoplasma ~ 1 + (1 | species)
mod_body: hemoplasma ~ body_size + (1 | species)
         npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
mod_null    2 416.38 425.22 -206.19    412.38                     
mod_body    4 417.54 435.22 -204.77    409.54 2.8398  2     0.2417

                model      AIC delta_AIC
1 body_size + species 417.5406  1.160212
2 null (species only) 416.3804  0.000000
```

Predicted infection probabilities  with confidence interval per mammalian class size :
```
emm_body <- emmeans(mod_body, ~ body_size, type = "response")
body_pred <- as.data.frame(emm_body)
body_pred
```

Results are : 
```
body_size      prob        SE  df   asymp.LCL asymp.UCL
 Large     0.5705841 0.5767685 Inf 0.013002934 0.9925935
 Medium    0.0857229 0.0578205 Inf 0.021605423 0.2847425
 Small     0.0258783 0.0243023 Inf 0.003999273 0.1494869
```


XXXXXXXXXXXXXXXXXXXXXXXXXX









## Step 3. Hemoplasma prevalence analysis across host species

Summarize by species:

```
df_species <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),                             
    n_infected = sum(hemoplasma == 1, na.rm = TRUE),  # infected individuals
    prevalence = n_infected / n  
  )
print(df_species, n = Inf)
```

Results are:
```
# A tibble: 44 × 4
   species                       n n_infected prevalence
   <fct>                     <int>      <int>      <dbl>
 1 Alouatta_macconnelli         22         20     0.909 
 2 Bradypus_tridactylus        108          4     0.0370
 3 Cabassous_unicinctus          2          0     0     
 4 Caluromys_philander           5          0     0     
 5 Sapajus_apella                  1          0     0     
 6 Choloepus_didactylus         90         72     0.8   
 7 Coendou_melanurus             1          0     0     
 8 Coendou_prehensilis                    3          1     0.333 
 9 Cyclopes_didactylus           1          0     0     
10 Dasypus_novemcinctus         15          5     0.333 
11 Didelphis_marsupialis        51         22     0.431 
12 Eira_barbara                  4          0     0     
13 Leopardus_wiedii                  1          0     0     
14 Galictis_vittata              4          3     0.75  
15 Holochilus_sciureus           5          1     0.2   
16 Hydrochoerus_hydrochaeris     2          0     0     
17 Hylaeamys_megacephalus       15          0     0     
18 Hylaeamys_yunganus           10          0     0     
19 Lontra_longicaudis            1          1     1     
20 Makalata_didelphoides         8          0     0     
21 Marmosa_lepida                1          1     1     
22 Marmosa_murina               20          2     0.1   
23 Marmosops_parvidens           5          0     0     
24 Mesomys_hispidus             13          0     0     
25 Metachirus_nudicaudatus       5          0     0     
26 Micoureus_demerarae          16          0     0     
27 Mus_musculus                 34          0     0     
28 Neacomys_dubosti              1          0     0     
29 Neacomys_paracou              8          0     0     
30 Nectomys_rattus               4          2     0.5   
31 Oecomys_auyantepui           16          1     0.0625
32 Oecomys_bicolor              16          0     0     
33 Oligoryzomys_fulvescens       7          1     0.143 
34 Philander_opossum            20          8     0.4   
35 Pithecia_pithecia             1          0     0     
36 Potos_flavus                  2          1     0.5   
37 Proechimys_cuvieri           18          1     0.0556
38 Proechimys_guyannensis       20          1     0.05  
39 Puma_yagouaroundi             5          0     0     
40 Rattus_rattus                19          1     0.0526
41 Saguinus_midas               41         41     1     
42 Saimiri_sciureus              1          0     0     
43 Sciurus_aestuans              1          0     0     
44 Tamandua_tetradactyla         3          0     0     
```

Scatter plot
```
ggplot(df_species, aes(x = n, y = prevalence)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(
    x = "Number of sampled individuals per species",
    y = "Hemoplasma prevalence"
  )
```

Binomial GLM with all species
```
glm_model <- glm(
  cbind(n_infected, n - n_infected) ~ n,
  family = binomial,
  data = df_species
)
summary(glm_model)

# Null model
glm_model_null <- glm(
  cbind(n_infected, n - n_infected) ~ 1,
  family = binomial,
  data = df_species
)

# Likelihood ratio test
anova_all <- anova(glm_model_null, glm_model, test = "Chisq")
print(anova_all)
```

Results are:
```
Analysis of Deviance Table
Model 1: cbind(n_infected, n - n_infected) ~ 1
Model 2: cbind(n_infected, n - n_infected) ~ n
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1        43     442.01                          
2        42     428.11  1   13.903 0.0001925 ***
```

Binomial GLM filtering species with n < 5
```
df_species_filtered <- df_species %>%
  filter(n >= 5)
# GLM with filtered data
glm_model_f <- glm(
  cbind(n_infected, n - n_infected) ~ n,
  family = binomial,
  data = df_species_filtered
)
summary(glm_model_f)

# Null model for filtered data
glm_model_null_f <- glm(
  cbind(n_infected, n - n_infected) ~ 1,
  family = binomial,
  data = df_species_filtered
)

# Likelihood ratio test
anova_filtered <- anova(glm_model_null_f, glm_model_f, test = "Chisq")
print(anova_filtered)
```

Results are:
```
Analysis of Deviance Table
Model 1: cbind(n_infected, n - n_infected) ~ 1
Model 2: cbind(n_infected, n - n_infected) ~ n
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1        25     419.11                          
2        24     405.19  1   13.916 0.0001912 ***
```

## Step 5. Hemoplasma prevalence by species

Prepare species-level data and keep only species with at least 1 infected individual
```
df_species <- data_hemoplasma_stat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_infected = sum(as.numeric(as.character(hemoplasma)) == 1, na.rm = TRUE),
    prevalence = n_infected / n,
    .groups = "drop"
  )
df_infected <- df_species %>%
  filter(n_infected > 0)
```

Binomial GLM (for all infected species)
```
glm_all <- glm(
  cbind(n_infected, n - n_infected) ~ species,
  family = binomial,
  data = df_infected
)
glm_null_all <- glm(
  cbind(n_infected, n - n_infected) ~ 1,
  family = binomial,
  data = df_infected
)
anova_all <- anova(glm_null_all, glm_all, test = "Chisq")
print("=== All infected species ===")
print(summary(glm_all))
print(anova_all)
```

Results are:
```
Call:
glm(formula = cbind(n_infected, n - n_infected) ~ species, family = binomial, 
    data = df_infected)
Coefficients:
                                 Estimate Std. Error z value Pr(>|z|)    
(Intercept)                        2.3026     0.7416   3.105  0.00190 ** 
speciesBradypus_tridactylus       -5.5607     0.8998  -6.180 6.41e-10 ***
speciesCholoepus_didactylus       -0.9163     0.7870  -1.164  0.24434    
speciesCoendou_prehensilis                 -2.9957     1.4318  -2.092  0.03641 *  
speciesDasypus_novemcinctus       -2.9957     0.9220  -3.249  0.00116 ** 
speciesDidelphis_marsupialis      -2.5788     0.7937  -3.249  0.00116 ** 
speciesGalictis_vittata           -1.2040     1.3723  -0.877  0.38032    
speciesHolochilus_sciureus        -3.6889     1.3416  -2.750  0.00597 ** 
speciesLontra_longicaudis         21.2635 79462.0195   0.000  0.99979    
speciesMarmosa_lepida             21.2635 79461.9966   0.000  0.99979    
speciesMarmosa_murina             -4.4998     1.0515  -4.280 1.87e-05 ***
speciesNectomys_rattus            -2.3026     1.2450  -1.849  0.06439 .  
speciesOecomys_auyantepui         -5.0106     1.2715  -3.941 8.12e-05 ***
speciesOligoryzomys_fulvescens    -4.0943     1.3102  -3.125  0.00178 ** 
speciesPhilander_opossum          -2.7081     0.8708  -3.110  0.00187 ** 
speciesPotos_flavus               -2.3026     1.5969  -1.442  0.14932    
speciesProechimys_cuvieri         -5.1358     1.2684  -4.049 5.14e-05 ***
speciesProechimys_guyannensis     -5.2470     1.2660  -4.145 3.40e-05 ***
speciesRattus_rattus              -5.1930     1.2671  -4.098 4.16e-05 ***
speciesSaguinus_midas             24.1352 52162.4892   0.000  0.99963    

(Dispersion parameter for binomial family taken to be 1)
    Null deviance: 3.0552e+02  on 19  degrees of freedom
Residual deviance: 5.0351e-10  on  0  degrees of freedom
AIC: 81.765
Number of Fisher Scoring iterations: 22

Analysis of Deviance Table
Model 1: cbind(n_infected, n - n_infected) ~ 1
Model 2: cbind(n_infected, n - n_infected) ~ species
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1        19     305.52                          
2         0       0.00 19   305.52 < 2.2e-16 ***
```

Binomial GLM (for infected species with n >= 5)
```
df_infected_f <- df_infected %>% filter(n >= 5)
glm_filtered <- glm(
  cbind(n_infected, n - n_infected) ~ species,
  family = binomial,
  data = df_infected_f
)
glm_null_filtered <- glm(
  cbind(n_infected, n - n_infected) ~ 1,
  family = binomial,
  data = df_infected_f
)
anova_filtered <- anova(glm_null_filtered, glm_filtered, test = "Chisq")
print("=== Species with n >= 5 ===")
print(summary(glm_filtered))
print(anova_filtered)
```
Results are:
```
Call: glm(formula = cbind(n_infected, n - n_infected) ~ species, family = binomial, 
    data = df_infected_f)
Coefficients:
                                 Estimate Std. Error z value Pr(>|z|)    
(Intercept)                        2.3026     0.7416   3.105  0.00190 ** 
speciesBradypus_tridactylus       -5.5607     0.8998  -6.180 6.41e-10 ***
speciesCholoepus_didactylus       -0.9163     0.7870  -1.164  0.24434    
speciesDasypus_novemcinctus       -2.9957     0.9220  -3.249  0.00116 ** 
speciesDidelphis_marsupialis      -2.5788     0.7937  -3.249  0.00116 ** 
speciesHolochilus_sciureus        -3.6889     1.3416  -2.750  0.00597 ** 
speciesMarmosa_murina             -4.4998     1.0515  -4.280 1.87e-05 ***
speciesOecomys_auyantepui         -5.0106     1.2715  -3.941 8.12e-05 ***
speciesOligoryzomys_fulvescens    -4.0943     1.3102  -3.125  0.00178 ** 
speciesPhilander_opossum          -2.7081     0.8708  -3.110  0.00187 ** 
speciesProechimys_cuvieri         -5.1358     1.2684  -4.049 5.14e-05 ***
speciesProechimys_guyannensis     -5.2470     1.2660  -4.145 3.40e-05 ***
speciesRattus_rattus              -5.1930     1.2671  -4.098 4.16e-05 ***
speciesSaguinus_midas             24.1352 52162.4892   0.000  0.99963    

(Dispersion parameter for binomial family taken to be 1)
    Null deviance: 2.9957e+02  on 13  degrees of freedom
Residual deviance: 2.7045e-10  on  0  degrees of freedom
AIC: 63.069
Number of Fisher Scoring iterations: 22

Analysis of Deviance Table
Model 1: cbind(n_infected, n - n_infected) ~ 1
Model 2: cbind(n_infected, n - n_infected) ~ species
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1        13     299.57                          
2         0       0.00 13   299.57 < 2.2e-16 ***
```

Scatter plot for infected species
```
plot_file <- "hemoplasma_prevalence_species.pdf"
pdf(plot_file, width = 8, height = 5)  # PDF output
ggplot(df_infected, aes(x = species, y = prevalence, size = n)) +
  geom_point(alpha = 0.5, color = "#69b3a2") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Species",
    y = "Hemoplasma prevalence",
    size = "Sample size"
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_size_continuous(range = c(2, 8))
dev.off()
cat("PDF saved to:", plot_file, "\n")
```

## Step 6. Hemoplasma prevalence by species | order

GLMM GLOBAL
```
glmm_order <- glmer(
  hemoplasma ~ order + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(optimizer = "bobyqa")
)
glmm_null <- glmer(
  hemoplasma ~ 1 + (1 | species),
  family = binomial,
  data = data_hemoplasma_stat,
  control = glmerControl(optimizer = "bobyqa")
)
anova(glmm_null, glmm_order, test = "Chisq")
summary(glmm_order)

Results are:
```
Data: data_hemoplasma_stat
Models:
glmm_null: hemoplasma ~ 1 + (1 | species)
glmm_order: hemoplasma ~ order + (1 | species)
           npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)  
glmm_null     2 416.38 425.22 -206.19    412.38                       
glmm_order    7 416.19 447.13 -201.10    402.19 10.187  5    0.07011 .

Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: hemoplasma ~ order + (1 | species)
   Data: data_hemoplasma_stat
Control: glmerControl(optimizer = "bobyqa")

      AIC       BIC    logLik -2*log(L)  df.resid 
    416.2     447.1    -201.1     402.2       607 

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-2.9140 -0.2196 -0.1671  0.1416  4.8727 

Random effects:
 Groups  Name        Variance Std.Dev.
 species (Intercept) 3.633    1.906   
Number of obs: 614, groups:  species, 44

Fixed effects:
                     Estimate Std. Error z value Pr(>|z|)  
(Intercept)            -0.748      1.154  -0.648   0.5167  
orderCingulata         -1.034      2.046  -0.505   0.6132  
orderDidelphimorphia   -1.449      1.439  -1.007   0.3141  
orderPilosa            -1.094      1.661  -0.659   0.5102  
orderPrimates           1.732      1.668   1.039   0.2990  
orderRodentia          -2.950      1.327  -2.222   0.0263 *

Correlation of Fixed Effects:
            (Intr) ordrCn ordrDd ordrPl ordrPr
orderCinglt -0.552                            
ordrDdlphmr -0.788  0.459                     
orderPilosa -0.680  0.402  0.566              
orderPrimts -0.675  0.403  0.567  0.498       
orderRodent -0.853  0.500  0.706  0.617  0.619
```


FIGURE TRAITS
```
library(dplyr)
library(ggplot2)
library(tidyr)

# =========================
# 1. AGGREGATION PAR ESPECE
# =========================
df_species <- data_hemoplasma_stat %>%
  group_by(species, body_size, strata, locomotion,
           activity, diet, sociality) %>%
  summarise(
    n = n(),
    n_infected = sum(hemoplasma == 1, na.rm = TRUE),
    prevalence = n_infected / n,
    .groups = "drop"
  )

# =========================
# 2. FORMAT LONG
# =========================
df_long <- df_species %>%
  pivot_longer(
    cols = c(body_size, strata, locomotion,
             activity, diet, sociality),
    names_to = "trait",
    values_to = "category"
  )

df_long$trait <- factor(df_long$trait,
                        levels = c("body_size",
                                   "strata",
                                   "locomotion",
                                   "activity",
                                   "diet",
                                   "sociality"))

# =========================
# 3. COULEURS PAR TRAIT (SEULEMENT 6)
# =========================
trait_colors <- c(
  body_size = "#A6CEE3",
  strata = "#B2DF8A",
  locomotion = "#FDBF6F",
  activity = "#CAB2D6",
  diet = "#FFFF99",
  sociality = "#FB9A99"
)

# =========================
# 4. PLOT
# =========================
p <- ggplot(df_long, aes(x = category, y = prevalence)) +
  
  geom_point(aes(size = n, fill = trait),
             shape = 21, color = "black", alpha = 0.85) +
  
  scale_fill_manual(values = trait_colors) +
  
  scale_size(range = c(2, 8)) +
  
  facet_wrap(~ trait, scales = "free_x", ncol = 3) +
  
  theme_classic(base_size = 14) +
  
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  
  labs(
    x = NULL,
    y = "Hemoplasma prevalence (per species)"
  )

p
ggsave(
  filename = "hemoplasma_traits_prevalence.pdf",
  plot = p,
  width = 12,
  height = 7
)
```

