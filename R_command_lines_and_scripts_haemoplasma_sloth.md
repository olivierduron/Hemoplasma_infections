# **Hemoplasma veterinary evaluations : R command lines and script**

We analyzed data from 175 wild sloths captured between 1994 and 1995 during the flooding of the Petit Saut Dam (5°03′43″ N, 53°03′00″ O) on the Sinnamary River (French Guiana, South America). The clinical data include the following variables for each examined sloth: 
- `species` : Sloth species (Bt: *Bradypus tridactylus*; Cd: *Choloepus didactylus*)
- `sex` : Sex of the sloth (F: female; M: male)
- `age_class` : Age category (A: adult; J: juvenile)
- `season` : Season of capture (W: wet; D: dry)
- `weight` : Body weight (quantitative variable, in kg)
- `total_length` : Total body length (quantitative variable, in cm)
- `wither_height` : Height at the withers (quantitative variable, in cm)
- `neck_size` : Neck circumference (quantitative variable, in cm)
- `temperature` : Body temperature (quantitative variable, in °C)
- `hematocrit` : Hematocrit level (quantitative variable, in %)
- `health_condition` : Overall health status (G: good; D: deteriorated)
- `female_reproductive_status` : Reproductive state of female individuals, indicating whether they were reproductively active or inactive at the time of sampling (Female non pregnant non lactating / Pregnant female / Female lactating with a young)
- `hemoplasma` : Infection status with hemotropic mycoplasmas (0 = uninfected; 1 = infected)
- `anaplasmataceae` : Infection status with bacteria of the Anaplasmataceae family (here, *Anaplasma amazonensis*) (0 = uninfected; 1 = infected)
- `apicomplexa` : Infection status with piroplasmids (here, *Babesia* sp.) (0 = uninfected; 1 = infected)
  
Details about all the experimental methods and measures are available in the related manuscript.


## Table of contents A CORRIGER!
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

All veterinary clinical data for the two sloth species are available [here](https://github.com/olivierduron/Hemoplasma_infections/blob/main/data_hemoplasma_sloth.csv).

This database will be referred to as `data_hemoplasma` throughout the R command lines and scripts provided below. It corresponds to the dataset provided in Table S1 of the related manuscript.

Load the dataset directly from the GitHub repository to R
```
data_hemoplasma <- read.csv2(
  "https://raw.githubusercontent.com/olivierduron/Hemoplasma_infections/main/data_hemoplasma_sloth.csv",
  na.strings = c("NA", "")
)
data_hemoplasma
```


## Step 2. Prepare the data for analysis

Convert categorical variables into factors
```
data_hemoplasma$species        <- as.factor(data_hemoplasma$species)
data_hemoplasma$season         <- as.factor(data_hemoplasma$season)
data_hemoplasma$sex            <- as.factor(data_hemoplasma$sex)
data_hemoplasma$age            <- as.factor(data_hemoplasma$age)
data_hemoplasma$hemoplasma      <- as.factor(data_hemoplasma$hemoplasma)
data_hemoplasma$anaplasmataceae      <- as.factor(data_hemoplasma$anaplasmataceae)
data_hemoplasma$apicomplexa       <- as.factor(data_hemoplasma$apicomplexa)
```

Load libraries for analysis
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

## Step 3. Calculate `hemoplasma` infection prevalence
Calculate `hemoplasma` infection prevalence and 95% confidence interval for _Bradypus tridactylus_ (Bt) and _Choloepus didactylus_ (Cd)

```
prevalence_results <- data_hemoplasma %>% group_by(species) %>% summarise(n = n(), positives = sum(hemoplasma == 1), prevalence = positives / n, conf_low = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$lower, conf_high = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$upper)
print(prevalence_results)
```

Results :
| species | n | positives | prevalence | conf_low | conf_high |
|:--------|--:|----------:|-----------:|---------:|----------:|
| Bt      | 92 | 4  | 0.0435 | 0.0120 | 0.108 |
| Cd      | 83 | 68 | 0.8190 | 0.7200 | 0.895 |

Test if `hemoplasma` is influenced by sloth `species`:
```
chisq.test(table(data_hemoplasma$hemoplasma, data_hemoplasma$species))
```

-> Results : `hemoplasma` prevalence differed strongly between the two sloth species, from 4.3% (4/92; 95% CI: 1.2–10.8%) in *Bradypus tridactylus* to 81.9% (68/83; 95% CI: 72.0–89.5%) in *Choloepus didactylus* (χ²₁ = 105.27, p < 2.2 × 10⁻¹⁶).

-> Interpretation : This strong interspecific difference suggests that host species and associated ecological or evolutionary traits may constrain `hemoplasma` infection.

## Step 4. Create the `pathogens` variable and sloth `species` subset

Create the `pathogens` variable by merging `anaplasmataceae` and `apicomplexa` (0 = uninfected ; 1 = infected by `anaplasmataceae` and/or `apicomplexa`)
```
data_hemoplasma <- data_hemoplasma %>%
  mutate(
    pathogens = ifelse(
      anaplasmataceae == 1 | apicomplexa == 1,
      1, 0
    ),
    species = factor(species)
  )
data_hemoplasma$pathogens <- as.factor(data_hemoplasma$pathogens)
```

Convert variables:
```
data_hemoplasma$weight <- as.numeric(data_hemoplasma$weight)
data_hemoplasma$total_length <- as.numeric(data_hemoplasma$total_length)
data_hemoplasma$wither_height <- as.numeric(data_hemoplasma$wither_height)
data_hemoplasma$neck_size <- as.numeric(data_hemoplasma$neck_size)
data_hemoplasma$temperature <- as.numeric(data_hemoplasma$temperature)
data_hemoplasma$hematocrit <- as.numeric(data_hemoplasma$hematocrit)
```

Create a subset `data_Bt` containing only records for _Bradypus tridactylus_ (Bt)
```
data_Bt <- subset(data_hemoplasma, species == "Bt")
table(data_Bt$hemoplasma)
```

Create a subset `data_Cd` containing only records for _Choloepus didactylus_ (Cd):
```
data_Cd <- subset(data_hemoplasma, species == "Cd")
table(data_Cd$hemoplasma)
```

## Step 5. Impact of `hemoplasma` infections on Scale Mass Index (SMI) in adult Bt
The Scaled Mass Index (SMI) was used as a body condition indicator that standardizes individual `weight` to `body_length`, using an allometric scaling relationship. SMI was calculated following Peig & Green (2009) (https://doi.org/10.1111/j.1600-0706.2009.17643.x).

Function to calculate SMI for adult Bt
```
data_adult_Bt <- subset(data_Bt, age == "A")
sma_model_Bt <- sma(log(weight) ~ log(total_length), data = data_adult_Bt)
b <- coef(sma_model_Bt)[2]
L0 <- mean(data_adult_Bt$total_length, na.rm = TRUE)
data_adult_Bt$SMI <- data_adult_Bt$weight * (L0 / data_adult_Bt$total_length)^b
```

Fit a GLM to test whether SMI is influenced by interactions among `hemoplasma`, `pathogens`, `sex` and `season` in Bt
```
model_SMIBt_full <- glm(
  SMI ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

summary(model_SMIBt_full)

model_SMIBt_3way <- glm(
  SMI ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

anova(model_SMIBt_3way, model_SMIBt_full, test = "Chisq")

AIC(model_SMIBt_full, model_SMIBt_3way)

model_SMIBt_2way <- glm(
  SMI ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

anova(model_SMIBt_2way, model_SMIBt_3way, test = "Chisq")

AIC(model_SMIBt_3way, model_SMIBt_2way)

model_SMIBt_add <- glm(
  SMI ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

anova(model_SMIBt_add, model_SMIBt_2way, test = "Chisq")

AIC(model_SMIBt_2way, model_SMIBt_add)

model_SMIBt_null <- glm(
  SMI ~ 1,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

anova(model_SMIBt_null, model_SMIBt_add, test = "Chisq")

AIC(model_SMIBt_null, model_SMIBt_add)

AIC_models_SMIBt <- AIC(
  model_SMIBt_full,
  model_SMIBt_3way,
  model_SMIBt_2way,
  model_SMIBt_add,
  model_SMIBt_null
)

AIC_models_SMIBt$delta_AIC <- 
  AIC_models_SMIBt$AIC - min(AIC_models_SMIBt$AIC)

AIC_models_SMIBt
```

-> Results and interpretation : Several interaction terms were not estimable because some combinations of predictors were absent from the dataset. The model was therefore simplified hierarchically by removing higher-order interaction terms. Removing three-way interactions did not significantly reduce model fit (LRT: χ²₁ = 0.31, *p* = 0.305), and removing all two-way interactions likewise resulted in no significant loss of fit (LRT: χ²₄ = 1.12, *p* = 0.437). The additive model had the lowest AIC (141.18) and was strongly supported over the null model (LRT: χ²₄ = 7.80, *p* < 0.001; ΔAIC = 16.19). We therefore retained the additive model for further analyses.

Fit a GLM to test whether SMI is influenced by additive effets of `hemoplasma`, `pathogens`, `sex` and `season` in Bt
```
model_SMIBt_add <- glm(SMI ~ hemoplasma + pathogens + season + sex, data = data_adult_Bt, family = gaussian(link = "identity"))
model_no_season <- update(model_SMIBt_add, . ~ . - season)
anova(model_no_season, model_SMIBt_add, test = "Chisq")
AIC(model_SMIBt_add, model_no_season)
model_SMIBt_reduced1 <- glm(SMI ~ hemoplasma + pathogens + sex, data = data_adult_Bt, family = gaussian(link = "identity"))
model_no_pathogens <- update(model_SMIBt_reduced1, . ~ . - pathogens)
anova(model_no_pathogens, model_SMIBt_reduced1, test = "Chisq")
AIC(model_SMIBt_reduced1, model_no_pathogens)
model_SMIBt_reduced2 <- glm(SMI ~ hemoplasma + sex, data = data_adult_Bt, family = gaussian(link = "identity"))
model_no_hemoplasma <- update(model_SMIBt_reduced2, . ~ . - hemoplasma)
anova(model_no_hemoplasma, model_SMIBt_reduced2, test = "Chisq")
AIC(model_SMIBt_reduced2, model_no_hemoplasma)
model_no_sex <- update(model_SMIBt_reduced2, . ~ . - sex)
anova(model_no_sex, model_SMIBt_reduced2, test = "Chisq")
AIC(model_SMIBt_reduced2, model_no_sex)
model_SMIBt_final <- model_SMIBt_reduced2
summary(model_SMIBt_final)
AIC(model_SMIBt_add, model_SMIBt_reduced1, model_SMIBt_reduced2, model_no_hemoplasma, model_no_sex)
```
-> Results : `season` and other blood-borne `pathogens` were excluded from the final model because their removal did not affect model fit (LRT: χ²₁ = 0.02, *p* = 0.804 and χ²₁ = 0.60, *p* = 0.153, respectively). The final model retained `hemoplasma` and `sex` (AIC = 139.36). SMI was lower in `hemoplasma`-positive individuals (β = −0.92 ± 0.28 SE, *p* = 0.0014) and higher in males (β = 0.44 ± 0.12 SE, *p* < 0.001).

-> Interpretation : `hemoplasma` infection was negatively associated with body condition in adult *B. tridactylus*, independently of `sex`, whereas `season` and other blood-borne `pathogens` showed no detectable effect.

Model diagnostics (diagnostic plots, Shapiro–Wilk test, and Breusch–Pagan test)
```
par(mfrow = c(2, 2))
plot(model_SMIBt_final)
par(mfrow = c(1, 1))
shapiro.test(residuals(model_SMIBt_final))
library(lmtest)
bptest(model_SMIBt_final)
```

-> Results and interpretation :Residuals showed no evidence of departure from normality (Shapiro–Wilk: W = 0.988, p = 0.661) or heteroscedasticity (Breusch–Pagan: BP = 4.10, p = 0.129), supporting the use of the Gaussian model.


Calculation of mean and standard error of SMI by `hemoplasma` infection status and `sex` for Bt
```
data_adult_Bt %>% 
  group_by(sex, hemoplasma) %>% 
  summarise(
    n = sum(!is.na(SMI)),
    mean = mean(SMI, na.rm = TRUE),
    se = sd(SMI, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  ) %>% 
  mutate(SMI = sprintf("%.2f ± %.2f", mean, se))
```

Results :
| sex | hemoplasma | n | mean | se | SMI |
|:----|:-----------|--:|-----:|----:|:----|
| F | 0 | 39 | 4.43 | 0.0749 | 4.43 ± 0.07 |
| F | 1 | 2 | 3.56 | 0.0284 | 3.56 ± 0.03 |
| M | 0 | 40 | 4.88 | 0.0985 | 4.88 ± 0.10 |
| M | 1 | 2 | 3.90 | 0.298 | 3.90 ± 0.30 |


Generate SMI chart for Bt
```
clean_data <- data_adult_Bt %>%
  filter(
    !is.na(weight), !is.na(total_length), !is.na(SMI),
    is.finite(weight), is.finite(total_length), is.finite(SMI)
  ) %>%
  mutate(
    sex_infect = case_when(
      sex == "M" & hemoplasma == 0 ~ "Male, uninfected",
      sex == "M" & hemoplasma == 1 ~ "Male, infected",
      sex == "F" & hemoplasma == 0 ~ "Female, uninfected",
      sex == "F" & hemoplasma == 1 ~ "Female, infected",
      TRUE ~ NA_character_
    )
  )

levels_order <- c(
  "Male, uninfected",
  "Male, infected",
  "Female, uninfected",
  "Female, infected"
)

clean_data <- clean_data %>%
  mutate(
    sex_infect = factor(sex_infect, levels = levels_order),
    point_size = case_when(
      sex_infect %in% c("Male, uninfected", "Male, infected") ~ 3.25,
      TRUE ~ 4
    )
  )

interp_data <- with(clean_data, akima::interp(
  x = weight,
  y = total_length,
  z = SMI,
  duplicate = "mean",
  extrap = FALSE
))

interp_df <- expand.grid(
  x = interp_data$x,
  y = interp_data$y
)

interp_df$z <- as.vector(interp_data$z)

legend_point_sizes <- c(3.25, 3.25, 4, 4) / 2

p_SMI_Bt <- ggplot() +
  geom_contour_filled(
    data = interp_df,
    aes(x = x, y = y, z = z)
  ) +
  geom_point(
    data = clean_data,
    aes(
      x = weight,
      y = total_length,
      shape = sex_infect,
      size = point_size
    ),
    color = "black",
    stroke = 1
  ) +
  scale_fill_brewer(
    palette = "Oranges",
    name = "SMI level"
  ) +
  scale_shape_manual(
    name = expression(paste("Hemoplasma", " infection status")),
    values = c(
      "Male, uninfected" = 0,
      "Male, infected" = 12,
      "Female, uninfected" = 1,
      "Female, infected" = 10
    )
  ) +
  scale_size_identity(guide = "none") +
  guides(
    shape = guide_legend(
      override.aes = list(size = legend_point_sizes)
    )
  ) +
  labs(
    x = "Body mass (kg)",
    y = "Total Length (cm)",
    title = expression(
      paste(
        "Scale Mass Index (SMI) of ",
        italic("Bradypus tridactylus")
      )
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    panel.background = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

p_SMI_Bt

ggsave(
  filename = "SMI_Bradypus_tridactylus.png",
  plot = p_SMI_Bt,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300
)
```

## Step 6. Impact of `hemoplasma` infections on Scale Mass Index (SMI) in adult Cd
Function to calculate SMI for adult Cd
```
data_adult_Cd <- subset(data_Cd, age == "A")
sma_model_Cd <- sma(log(weight) ~ log(total_length), data = data_adult_Cd)
b <- coef(sma_model_Cd)[2]
L0 <- mean(data_adult_Cd$total_length, na.rm = TRUE)
data_adult_Cd$SMI <- data_adult_Cd$weight * (L0 / data_adult_Cd$total_length)^b
```

Fit a GLM to test whether SMI is influenced by interactions among `hemoplasma`, `pathogens`, `sex`, and `season` in Cd
```
model_SMICd_full <- glm(
  SMI ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

summary(model_SMICd_full)

model_SMICd_3way <- glm(
  SMI ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

anova(model_SMICd_3way, model_SMICd_full, test = "Chisq")

AIC(model_SMICd_full, model_SMICd_3way)

model_SMICd_2way <- glm(
  SMI ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

anova(model_SMICd_2way, model_SMICd_3way, test = "Chisq")

AIC(model_SMICd_3way, model_SMICd_2way)

model_SMICd_add <- glm(
  SMI ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

anova(model_SMICd_add, model_SMICd_2way, test = "Chisq")

AIC(model_SMICd_2way, model_SMICd_add)

model_SMICd_null <- glm(
  SMI ~ 1,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

anova(model_SMICd_null, model_SMICd_add, test = "Chisq")

AIC(model_SMICd_null, model_SMICd_add)

AIC_models_SMICd <- AIC(
  model_SMICd_full,
  model_SMICd_3way,
  model_SMICd_2way,
  model_SMICd_add,
  model_SMICd_null
)

AIC_models_SMICd$delta_AIC <- 
  AIC_models_SMICd$AIC - min(AIC_models_SMICd$AIC)

AIC_models_SMICd
```

-> Results: For adult *Choloepus didactylus*, neither the interaction models nor the additive model improved model fit. Higher-order interactions were not supported (LRT: χ²₂ = 2.81, *p* = 0.087), and removing all two-way interactions did not reduce model fit (LRT: χ²₆ = 1.37, *p* = 0.896). The additive model also did not improve on the null model (LRT: χ²₄ = 1.02, *p* = 0.774) and had a higher AIC (136.32 vs. 130.25). The null model was therefore retained.

-> Interpretation : SMI in adult *C. didactylus* showed no detectable association with `hemoplasma` infection, `sex`, `season`, or other blood-borne `pathogens`.

Model diagnostics (Shapiro–Wilk test)
```
shapiro.test(residuals(model_SMICd_null))
```
-> Results and interpretation : Residuals of the null model showed no evidence of departure from normality (Shapiro–Wilk: W = 0.979, p = 0.428).

Generate SMI chart for Cd
```
clean_data <- data_adult_Cd %>%
  filter(
    !is.na(weight), !is.na(total_length), !is.na(SMI),
    is.finite(weight), is.finite(total_length), is.finite(SMI)
  ) %>%
  mutate(
    sex_infect = case_when(
      sex == "M" & hemoplasma == 0 ~ "Male, uninfected",
      sex == "M" & hemoplasma == 1 ~ "Male, infected",
      sex == "F" & hemoplasma == 0 ~ "Female, uninfected",
      sex == "F" & hemoplasma == 1 ~ "Female, infected",
      TRUE ~ NA_character_
    )
  )

levels_order <- c("Male, uninfected", "Male, infected", "Female, uninfected", "Female, infected")

clean_data <- clean_data %>%
  mutate(
    sex_infect = factor(sex_infect, levels = levels_order),
    point_size = case_when(
      sex_infect %in% c("Male, uninfected", "Male, infected") ~ 3.25,
      TRUE ~ 4
    )
  )

interp_data <- with(clean_data, akima::interp(
  x = weight,
  y = total_length,
  z = SMI,
  duplicate = "mean",
  extrap = FALSE
))

interp_df <- expand.grid(
  x = interp_data$x,
  y = interp_data$y
)

interp_df$z <- as.vector(interp_data$z)

legend_point_sizes <- c(3.25, 3.25, 4, 4) / 2

p_SMI_Cd <- ggplot() +
  geom_contour_filled(data = interp_df, aes(x = x, y = y, z = z)) +
  geom_point(
    data = clean_data,
    aes(
      x = weight,
      y = total_length,
      shape = sex_infect,
      size = point_size
    ),
    color = "black",
    stroke = 1
  ) +
  scale_fill_brewer(palette = "YlOrBr", name = "SMI level") +
  scale_shape_manual(
    name = expression(paste("Hemoplasma", " infection status")),
    values = c(
      "Male, uninfected" = 0,
      "Male, infected" = 12,
      "Female, uninfected" = 1,
      "Female, infected" = 10
    )
  ) +
  scale_size_identity(guide = "none") +
  guides(
    shape = guide_legend(override.aes = list(size = legend_point_sizes))
  ) +
  labs(
    x = "Body mass (kg)",
    y = "Total Length (cm)",
    title = expression(paste("Scale Mass Index (SMI) of ", italic("Choloepus didactylus")))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

p_SMI_Cd

ggsave(
  filename = "SMI_Choloepus_didactylus.png",
  plot = p_SMI_Cd,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300
)
```

## Step 7. Impact of `hemoplasma` infections on neck circumference

Fit a GLM to test whether neck circumference is influenced by interactions among `hemoplasma`, `pathogens`, `sex`, and `season` in Bt
```
model_5 <- glm(
  log(neck_size) ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

model_5_3way <- glm(
  log(neck_size) ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

model_5_2way <- glm(
  log(neck_size) ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

model_5_add <- glm(
  log(neck_size) ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

model_5_null <- glm(
  log(neck_size) ~ 1,
  data = data_adult_Bt,
  family = gaussian(link = "identity")
)

anova(model_5_3way, model_5, test = "Chisq")
anova(model_5_2way, model_5_3way, test = "Chisq")
anova(model_5_add, model_5_2way, test = "Chisq")
anova(model_5_null, model_5_add, test = "Chisq")

AIC(model_5, model_5_3way, model_5_2way, model_5_add, model_5_null)
```

-> Results: For adult *Bradypus tridactylus*, none of the interaction models improved model fit (three-way vs. two-way interactions: LRT, χ²₁ = 0.003, *p* = 0.554; two-way vs. additive model: χ²₃ = 0.010, *p* = 0.765). The additive model also did not improve on the null model (χ²₄ = 0.042, *p* = 0.298), and the null model had the lowest AIC (−112.64 vs. −109.75 for the additive model). The null model was therefore retained.

-> Interpretation : Neck size showed no detectable association with `hemoplasma` infection, other blood-borne `pathogens`, `season`, or `sex` in adult *Bradypus tridactylus*.

Model diagnostics (Shapiro–Wilk test)
```
shapiro.test(residuals(model_5_null))
```
-> Results and interpretation : Residuals showed no significant departure from normality (Shapiro–Wilk: W = 0.967, p = 0.093); this supports the Gaussian assumption.






Fit a GLM to test whether neck circumference is influenced by interactions among `hemoplasma`, `pathogens`, `sex`, and `season` in Cd:
```
model_6 <- glm(log(neck_size) ~ hemoplasma * bloodparasite * season * sex, data = data_adult_Cd, family = gaussian(link = "identity"))
```

Fit a GLM to test whether SMI is influenced by additive effects of `hemoplasma`, `bloodparasite`, `sex`, and `season` in Cd:
```
model_6a <- glm(log(neck_size) ~ hemoplasma + bloodparasite + season + sex, data = data_adult_Cd, family = gaussian(link = "identity"))
```

Compare the additive model (model_6a) to the interaction model (model_6) using a likelihood ratio test:
```
anova(model_6a, model_6, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: log(neck_size) ~ hemoplasma + bloodparasite + season + sex
Model 2: log(neck_size) ~ hemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        43    0.55930                     
2        36    0.42106  7  0.13824   0.1067
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_6, model_6a)
```

Results are:
```
         df       AIC
model_6  13 -65.11804
model_6a  6 -65.49074
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_6a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
log(neck_size) ~ hemoplasma + bloodparasite + season + sex
              Df Deviance     AIC scaled dev. Pr(>Chi)
<none>            0.55930 -65.491                     
hemoplasma    1  0.55935 -67.486     0.00426   0.9480
bloodparasite  1  0.56108 -67.338     0.15253   0.6961
season         1  0.57932 -65.802     1.68826   0.1938
sex            1  0.57644 -66.042     1.44873   0.2287
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_6a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                  AIC delta_AIC
<none>        -65.491   0.00000
hemoplasma   -67.486  -1.99574
bloodparasite -67.338  -1.84747
season        -65.802  -0.31174
sex           -66.042  -0.55127
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model6_null <- glm(log(neck_size) ~ 1, data = data_adult_Cd, family = gaussian(link = "identity"))
model6_hemoplasma <- glm(log(neck_size) ~ hemoplasma, data = data_adult_Cd, family = gaussian(link = "identity"))
model6_bloodparasite <- glm(log(neck_size) ~ bloodparasite, data = data_adult_Cd, family = gaussian(link = "identity"))
model6_season <- glm(log(neck_size) ~ season, data = data_adult_Cd, family = gaussian(link = "identity"))
model6_sex <- glm(log(neck_size) ~ sex, data = data_adult_Cd, family = gaussian(link = "identity"))
anova(model6_null, model6_hemoplasma, test="Chisq")
anova(model6_null, model6_bloodparasite, test="Chisq")
anova(model6_null, model6_season, test="Chisq")
anova(model6_null, model6_sex, test="Chisq")
aics <- AIC(model6_null, model6_hemoplasma, model6_bloodparasite, model6_season, model6_sex)
aic_null <- aics["model6_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ hemoplasma
  Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
1        47    0.60752                      
2        46    0.60512  1 0.0024025   0.6691
---
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ bloodparasite
  Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
1        47    0.60752                      
2        46    0.60334  1 0.0041827   0.5723
---
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ season
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        47    0.60752                     
2        46    0.58020  1 0.027321   0.1411
---
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        47    0.60752                     
2        46    0.58267  1 0.024851   0.1613
---
                           AIC delta_AIC_vs_null
model6_null          -69.52144       0.000000000
model6_hemoplasma   -67.71164       1.809799276
model6_bloodparasite -67.85305       1.668385770
model6_season        -69.73015      -0.208715492
model6_sex           -69.52621      -0.004777113
```

Fit a linear model to test the null hypothesis (`neck_size` ~ 1) in adult Cd, assessing model fit and checking residual normality:
```
model_6b <- glm(log(neck_size) ~ 1, data = data_adult_Cd, family = gaussian(link = "identity"))
anova(model_6b, model_6, test = "Chisq")
AIC(model_6b, model_6)
shapiro.test(model_6b$residuals)
```

Results are:
```
Analysis of Deviance Table
Model 1: log(neck_size) ~ 1
Model 2: log(neck_size) ~ hemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        47    0.60752                     
2        36    0.42106 11  0.18645   0.1433
---
         df       AIC
model_6b  2 -69.52144
model_6  13 -65.11804
---
Shapiro-Wilk normality test
data:  model_6b$residuals
W = 0.97242, p-value = 0.3137
```

## Step 9. Impact of _Anaplasma_ infections on hematocrit levels (GLM models 7, 8 and 9)
Fit a GLM to test whether `hematocrit` is influenced by interactions among `hemoplasma`, `bloodparasite`, `sex`, and `season` in Bt:
```
model_7 <- glm(hematocrit ~ hemoplasma * bloodparasite * season * sex, data = data_adult_Bt, family = Gamma(link = "log"))
```

Fit a GLM to test whether `hematocrit` is influenced by additive effects of hemoplasma`, `bloodparasite`, `sex`, and `season` in Bt:
```
model_7a <- glm(hematocrit ~ hemoplasma + bloodparasite + season + sex, data = data_adult_Bt, family = Gamma(link = "log"))
```

Compare the additive model (model_7a) to the interaction model (model_7) using a likelihood ratio test:
```
anova(model_7a, model_7, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: hematocrit ~ hemoplasma + bloodparasite + season + sex
Model 2: hematocrit ~ hemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        79     1.3166                     
2        74     1.1985  5  0.11806   0.2065
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_7, model_7a)
```

Results are:
```
         df      AIC
model_7  11 517.8352
model_7a  6 515.7467
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_7a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
hematocrit ~ hemoplasma + bloodparasite + season + sex
              Df Deviance    AIC scaled dev. Pr(>Chi)  
<none>             1.3166 515.75                       
hemoplasma    1   1.3209 514.00      0.2550  0.61358  
bloodparasite  1   1.3258 514.29      0.5474  0.45937  
season         1   1.3238 514.18      0.4285  0.51273  
sex            1   1.3825 517.67      3.9191  0.04774 *
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_7a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        515.75    0.0000
hemoplasma   514.00   -1.7450
bloodparasite 514.29   -1.4526
season        514.18   -1.5715
sex           517.67    1.9191
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model7_null <- glm(hematocrit ~ 1, data = data_adult_Bt, family = Gamma(link = "log"))
model7_hemoplasma <- glm(hematocrit ~ hemoplasma, data = data_adult_Bt, family = Gamma(link = "log"))
model7_bloodparasite <- glm(hematocrit ~ bloodparasite, data = data_adult_Bt, family = Gamma(link = "log"))
model7_season <- glm(hematocrit ~ season, data = data_adult_Bt, family = Gamma(link = "log"))
model7_sex <- glm(hematocrit ~ sex, data = data_adult_Bt, family = Gamma(link = "log"))
anova(model7_null, model7_hemoplasma, test="Chisq")
anova(model7_null, model7_bloodparasite, test="Chisq")
anova(model7_null, model7_season, test="Chisq")
anova(model7_null, model7_sex, test="Chisq")
aics <- AIC(model7_null, model7_hemoplasma, model7_bloodparasite, model7_season, model7_sex)
aic_null <- aics["model7_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
Analysis of Deviance Table
Model 1: hematocrit ~ 1
Model 2: hematocrit ~ hemoplasma
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        83     1.4015                     
2        82     1.3936  1 0.007882   0.4994
---
Analysis of Deviance Table
Model 1: hematocrit ~ 1
Model 2: hematocrit ~ bloodparasite
  Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
1        83     1.4015                      
2        82     1.3950  1 0.0064542   0.5418
---
Analysis of Deviance Table
Model 1: hematocrit ~ 1
Model 2: hematocrit ~ season
  Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
1        83     1.4015                      
2        82     1.3926  1 0.0089277   0.4709
---
Analysis of Deviance Table
Model 1: hematocrit ~ 1
Model 2: hematocrit ~ sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
1        83     1.4015                       
2        82     1.3423  1 0.059141  0.05883 .
---
                          AIC delta_AIC_vs_null
model7_null          513.0105          0.000000
model7_hemoplasma   514.5354          1.524936
model7_bloodparasite 514.6217          1.611190
model7_season        514.4722          1.461711
model7_sex           511.3790         -1.631481
```

Fit a linear model to test the null hypothesis (`hematocrit` ~ 1) in adult Bt, assessing model fit:
```
model_7b <- glm(hematocrit ~ 1, data = data_adult_Bt, family = Gamma(link = "log"))
anova(model_7b, model_7, test = "Chisq")
AIC(model_7b, model_7)
```

Results are:
```
Analysis of Deviance Table
Model 1: hematocrit ~ 1
Model 2: hematocrit ~ hemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        83     1.4015                     
2        74     1.1985  9  0.20297   0.1932
---
         df      AIC
model_7b  2 513.0105
model_7  11 517.8352
```

Generate QQ plot for model_7b to assess normality of residuals:
```
res_dev <- residuals(model_7b, type = "deviance")
qqnorm(res_dev, main = "QQ Plot (Deviance residuals)")
qqline(res_dev, col = "red")
```

Calculation of mean and standard error of `hematocrit` by `hemoplasma` for Bt:
```
data_adult_Bt %>%
  group_by(hemoplasma) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results are:
```
A tibble: 2 × 3
  hemoplasma mean_hematocrit se_hematocrit
  <fct>                 <dbl>         <dbl>
1 0                      39.0         0.569
2 1                      40.8         3.04 
```

Fit a GLM to test whether `hematocrit` is influenced by interactions among `hemoplasma`, `bloodparasite`, `sex`, and `season` in Cd:
```
model_8 <- glm(hematocrit ~ hemoplasma * bloodparasite * season * sex, data = data_adult_Cd, family = Gamma(link = "log"))
```

Fit a GLM to test whether `hematocrit` is influenced by additive effects of `hemoplasma`, `bloodparasite`, `sex`, and `season` in Cd:
```
model_8a <- glm(hematocrit ~ hemoplasma + bloodparasite + season + sex, data = data_adult_Cd, family = Gamma(link = "log"))
```

Compare the additive model (model_8a) to the interaction model (model_8) using a likelihood ratio test:
```
anova(model_8a, model_8, test = "Chisq")
```

Results are:
```
Analysis of Deviance Table
Model 1: hematocrit ~ hemoplasma + bloodparasite + season + sex
Model 2: hematocrit ~ hemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        55     1.1832                     
2        49     1.1262  6 0.057072   0.8368
```

Compute AIC for both models to evaluate model fit:
```
AIC(model_8, model_8a)
```

Results are:
```
         df      AIC
model_8  12 393.2948
model_8a  6 384.2705
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_8a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
hematocrit ~ hemoplasma + bloodparasite + season + sex
              Df Deviance    AIC scaled dev. Pr(>Chi)   
<none>             1.1832 384.27                        
hemoplasma    1   1.1883 382.54      0.2645 0.607020   
bloodparasite  1   1.2285 384.64      2.3705 0.123646   
season         1   1.3174 389.30      7.0307 0.008013 **
sex            1   1.2199 384.19      1.9221 0.165622   
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_8a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        384.27    0.0000
hemoplasma   382.54   -1.7355
bloodparasite 384.64    0.3705
season        389.30    5.0307
sex           384.19   -0.0779
```

Fit a linear model to test the model_8b (`hematocrit` ~ `season`) in adult Cd, assessing model fit:
```
model_8b <- glm(hematocrit ~ season, data = data_adult_Cd, family = Gamma(link = "log"))
anova(model_8b, model_8, test = "Chisq")
AIC(model_8b, model_8)
```

Results are:
```
Analysis of Deviance Table
Model 1: hematocrit ~ season
Model 2: hematocrit ~ hemoplasma * bloodparasite * season * sex
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1        58     1.2872                     
2        49     1.1262  9  0.16103    0.552
---
         df      AIC
model_8b  3 383.3404
model_8  12 393.2948
```

Generate QQ plot for model_8b to assess normality of residuals:
```
res_dev <- residuals(model_8b, type = "deviance")
qqnorm(res_dev, main = "QQ Plot (Deviance residuals)")
qqline(res_dev, col = "red")
```

Calculation of mean and standard error of `hematocrit` by `hemoplasma` for Cd:
```
data_adult_Cd %>%
  group_by(hemoplasma) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results are:
```
A tibble: 2 × 3
  hemoplasma mean_hematocrit se_hematocrit
  <fct>                 <dbl>         <dbl>
1 0                      39.1         1.48 
2 1                      38.6         0.798
```

Calculation of mean and standard error of `hematocrit` by `season` for Cd:
```
data_adult_Cd %>%
  group_by(season) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results are:
```
A tibble: 2 × 3
  season mean_hematocrit se_hematocrit
  <fct>            <dbl>         <dbl>
1 D                 37.6         0.892
2 W                 41.1         0.944
```

Create violin plots for `hematocrit`
```
label_style <- element_text(size = 28, face = "bold") 

pA <- ggplot(data_adult_Bt, aes(x = factor(hemoplasma, levels = c(0, 1),
                                labels = c("Uninfected", "Infected")),
                                y = hematocrit)) +
  geom_violin(fill = "darkolivegreen3", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = NULL,  
       y = "Hematocrit (%)",   
       title = "A") +
  scale_y_continuous(limits = c(20, 60)) +
  scale_x_discrete(labels = NULL) + 
  theme_minimal() +
  theme(plot.title = label_style)

pB <- ggplot(data_adult_Bt, aes(x = factor(season, levels = c("D", "W"),
                                labels = c("Dry", "Wet")),
                                y = hematocrit)) +
  geom_violin(fill = "darkolivegreen3", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = NULL,
       y = NULL,
       title = "B") +
  scale_y_continuous(limits = c(20, 60)) +
  scale_x_discrete(labels = NULL) +  
  theme_minimal() +
  theme(plot.title = label_style)

pC <- ggplot(data_adult_Bt, aes(x = factor(sex, levels = c("M", "F"),
                                labels = c("Male", "Female")),
                                y = hematocrit)) +
  geom_violin(fill = "darkolivegreen3", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = NULL,
       y = NULL,
       title = "C") +
  scale_y_continuous(limits = c(20, 60)) +
  scale_x_discrete(labels = NULL) +  
  theme_minimal() +
  theme(plot.title = label_style)

pD <- ggplot(data_adult_Cd, aes(x = factor(hemoplasma, levels = c(0, 1),
                                labels = c("Uninfected", "Infected")),
                                y = hematocrit)) +
  geom_violin(fill = "goldenrod1", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = expression(paste("Hemoplasma", " infection status")),
       y = "Hematocrit (%)",
       title = "D") +
  scale_y_continuous(limits = c(10, 60)) +
  theme_minimal() +
  theme(plot.title = label_style)

pE <- ggplot(data_adult_Cd, aes(x = factor(season, levels = c("D", "W"),
                                labels = c("Dry", "Wet")),
                                y = hematocrit)) +
  geom_violin(fill = "goldenrod1", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = "Season",
       y = NULL,
       title = "E") +
  scale_y_continuous(limits = c(10, 60)) +
  theme_minimal() +
  theme(plot.title = label_style)

pF <- ggplot(data_adult_Cd, aes(x = factor(sex, levels = c("M", "F"),
                                labels = c("Male", "Female")),
                                y = hematocrit)) +
  geom_violin(fill = "goldenrod1", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = "Sex",
       y = NULL,
       title = "F") +
  scale_y_continuous(limits = c(10, 60)) +
  theme_minimal() +
  theme(plot.title = label_style)
final_plot <- (pA | pB | pC) / (pD | pE | pF)
print(final_plot)
```

## Step 10. Impact of Hemoplasma infections on body temperature (CLRM models 10 and 11)
Convert `temperature` to numeric, handle left-censored values (<32°C) for analysis in Bt:
```
data_adult_Bt <- data_adult_Bt %>%
  mutate(
    temperature_numeric = as.numeric(ifelse(temperature == "< 32.00", 32, temperature)),
    censored = ifelse(temperature == "< 32.00", TRUE, FALSE)
  )
```

Create a left-censored Surv object (temp) for `temperature` in Bt:
```
temp <- Surv(data_adult_Bt$temperature_numeric,
                  event = !data_adult_Bt$censored,
                  type = "left")
```

Fit Gaussian survival regression models to test the effects of `hemoplasma`, `bloodparasite`, `season`, `sex` on `temperature`, with (model_10) and without (model_10b) interactions in Bt:
```
model_10 <- survreg(temp ~ hemoplasma * bloodparasite * season * sex, data = data_adult_Bt, dist = "gaussian")
model_10a <- survreg(temp ~ hemoplasma + bloodparasite + season + sex, data = data_adult_Bt, dist = "gaussian")

```

Compare models using ANOVA and AIC to evaluate the contribution of interaction terms:
```
anova(model_10a, model_10, test = "Chisq")
AIC(model_10a, model_10)
```

Results are:
```
                                       Terms Resid. Df    -2*LL Test Df  Deviance  Pr(>Chi)
1 hemoplasma + bloodparasite + season + sex        27 79.74336      NA        NA        NA
2 hemoplasma * bloodparasite * season * sex        16 79.05558    = 11 0.6877884 0.9999927
---
          df       AIC
model_10a  6  91.74336
model_10  17 113.05558
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_10a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
temp ~ hemoplasma + bloodparasite + season + sex
              Df    AIC    LRT Pr(>Chi)  
<none>           91.743                  
hemoplasma    1 91.269 1.5261   0.2167  
bloodparasite  1 89.952 0.2084   0.6480  
season         1 90.275 0.5320   0.4658  
sex            1 94.666 4.9228   0.0265 *
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_10a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        91.743   0.00000
hemoplasma   91.269  -0.47387
bloodparasite 89.952  -1.79159
season        90.275  -1.46797
sex           94.666   2.92284
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model10_null <- survreg(temp ~ 1, data = data_adult_Bt, dist = "gaussian")
model10_hemoplasma <- survreg(temp ~ hemoplasma, data = data_adult_Bt, dist = "gaussian")
model10_bloodparasite <- survreg(temp ~ bloodparasite, data = data_adult_Bt, dist = "gaussian")
model10_season <- survreg(temp ~ season, data = data_adult_Bt, dist = "gaussian")
model10_sex <- survreg(temp ~ sex, data = data_adult_Bt, dist = "gaussian")
anova(model10_null, model10_hemoplasma, test="Chisq")
anova(model10_null, model10_bloodparasite, test="Chisq")
anova(model10_null, model10_season, test="Chisq")
anova(model10_null, model10_sex, test="Chisq")
aics <- AIC(model10_null, model10_hemoplasma, model10_bloodparasite, model10_season, model10_sex)
aic_null <- aics["model10_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
        Terms Resid. Df    -2*LL Test Df  Deviance  Pr(>Chi)
1           1        31 85.36015      NA        NA        NA
2 hemoplasma        30 84.81430    =  1 0.5458493 0.4600186
---
          Terms Resid. Df    -2*LL Test Df   Deviance  Pr(>Chi)
1             1        31 85.36015      NA         NA        NA
2 bloodparasite        30 85.27469    =  1 0.08546037 0.7700297
---
   Terms Resid. Df    -2*LL Test Df    Deviance  Pr(>Chi)
1      1        31 85.36015      NA          NA        NA
2 season        30 85.35333    =  1 0.006827965 0.9341446
---
  Terms Resid. Df    -2*LL Test Df Deviance   Pr(>Chi)
1     1        31 85.36015      NA       NA         NA
2   sex        30 81.66861    =  1 3.691544 0.05468897
---
                           AIC delta_AIC_vs_null
model10_null          89.36015          0.000000
model10_hemoplasma   90.81430          1.454151
model10_bloodparasite 91.27469          1.914540
model10_season        91.35333          1.993172
model10_sex           87.66861         -1.691544
```

Generate a QQ-plot of deviance residuals from model11_null to visually assess normality:
```
resid_temp <- residuals(model10_null, type = "deviance")
qqnorm(resid_temp)
qqline(resid_temp, col = "red", lwd = 1)
```

Convert `temperature` to numeric, handle left-censored values (<32°C) for sanalysis in Cd:
```
data_adult_Cd <- data_adult_Cd %>%
  mutate(
    temperature_numeric = as.numeric(ifelse(temperature == "< 32.00", 32, temperature)),
    censored = ifelse(temperature == "< 32.00", TRUE, FALSE)
  )
```

Create a left-censored Surv object (temp) for `temperature` in Cd:
```
temp <- Surv(data_adult_Cd$temperature_numeric,
                  event = !data_adult_Cd$censored,
                  type = "left")
```

Fit Gaussian survival regression models to test the effects of `hemoplasma`, `bloodparasite`, `season`, `sex` on `temperature`, with (model_11) and without (model_11b) interactions in Cd:
```
model_11 <- survreg(temp ~ hemoplasma * bloodparasite * season * sex, data = data_adult_Cd, dist = "gaussian")
model_11a <- survreg(temp ~ hemoplasma + bloodparasite + season + sex, data = data_adult_Cd, dist = "gaussian")
```

Compare models using ANOVA and AIC to evaluate the contribution of interaction terms:
```
anova(model_11a, model_11, test = "Chisq")
AIC(model_11a, model_11)
```

Results are:
```
                                       Terms Resid. Df    -2*LL Test Df Deviance  Pr(>Chi)
1 hemoplasma + bloodparasite + season + sex        13 58.22666      NA       NA        NA
2 hemoplasma * bloodparasite * season * sex         2 54.04089    = 11  4.18577 0.9641642
---
          df      AIC
model_11a  6 70.22666
model_11  17 88.04089
```

Perform drop-one-term analysis on the additive model:
```
res <- drop1(model_11a, test = "Chisq")
```

Results are:
```
Single term deletions
Model:
temp ~ hemoplasma + bloodparasite + season + sex
              Df    AIC      LRT Pr(>Chi)
<none>           70.227                  
hemoplasma    1 68.390 0.163699   0.6858
bloodparasite  1 68.235 0.008139   0.9281
season         1 68.269 0.042462   0.8367
sex            1 68.405 0.178242   0.6729
```

Calculate delta AIC for each term to assess its contribution to model fit:
```
aic_full <- AIC(model_11a)
res$delta_AIC <- res$AIC - aic_full
print(res[, c("AIC", "delta_AIC")])
```

Results are:
```
                 AIC delta_AIC
<none>        70.227    0.0000
hemoplasma   68.390   -1.8363
bloodparasite 68.235   -1.9919
season        68.269   -1.9575
sex           68.405   -1.8218
```

Compare the null model (model_null) to univariate models using likelihood ratio tests and AIC:
```
model11_null <- survreg(temp ~ 1, data = data_adult_Cd, dist = "gaussian")
model11_hemoplasma <- survreg(temp ~ hemoplasma, data = data_adult_Cd, dist = "gaussian")
model11_bloodparasite <- survreg(temp ~ bloodparasite, data = data_adult_Cd, dist = "gaussian")
model11_season <- survreg(temp ~ season, data = data_adult_Cd, dist = "gaussian")
model11_sex <- survreg(temp ~ sex, data = data_adult_Cd, dist = "gaussian")
anova(model11_null, model11_hemoplasma, test="Chisq")
anova(model11_null, model11_bloodparasite, test="Chisq")
anova(model11_null, model11_season, test="Chisq")
anova(model11_null, model11_sex, test="Chisq")
aics <- AIC(model11_null, model11_hemoplasma, model11_bloodparasite, model11_season, model11_sex)
aic_null <- aics["model11_null", "AIC"]
aics$delta_AIC_vs_null <- aics$AIC - aic_null
print(aics[, c("AIC", "delta_AIC_vs_null")])
```

Results are:
```
        Terms Resid. Df    -2*LL Test Df  Deviance  Pr(>Chi)
1           1        17 58.73649      NA        NA        NA
2 hemoplasma        16 58.42303    =  1 0.3134589 0.5755654
---
          Terms Resid. Df    -2*LL Test Df   Deviance  Pr(>Chi)
1             1        17 58.73649      NA         NA        NA
2 bloodparasite        16 58.67821    =  1 0.05828631 0.8092253
---
   Terms Resid. Df    -2*LL Test Df   Deviance  Pr(>Chi)
1      1        17 58.73649      NA         NA        NA
2 season        16 58.64685    =  1 0.08964222 0.7646325
---
  Terms Resid. Df    -2*LL Test Df  Deviance  Pr(>Chi)
1     1        17 58.73649      NA        NA        NA
2   sex        16 58.60173    =  1 0.1347589 0.7135479
---
                           AIC delta_AIC_vs_null
model11_null          62.73649          0.000000
model11_hemoplasma   64.42303          1.686541
model11_bloodparasite 64.67821          1.941714
model11_season        64.64685          1.910358
model11_sex           64.60173          1.865241
```

Generate a QQ-plot of deviance residuals from model11_null to visually assess normality:
```
resid_temp <- residuals(model11_null, type = "deviance")
qqnorm(resid_temp)
qqline(resid_temp, col = "red", lwd = 1)
```

## Step 11. Impact of Hemoplasma infections on general health condition 
Test the association between `hemoplasma` and `health_condition` in Bt:
```
table_health_condition_hemoplasma_Bt <- table(data_Bt$hemoplasma, data_Bt$health_condition)
table_health_condition_hemoplasma_Bt
fisher.test(table_health_condition_hemoplasma_Bt)
```

Results are:
```
     D  G
  0 10 78
  1  0  4
---
Fisher's Exact Test for Count Data
data:  table_health_condition_hemoplasma_Bt
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.0750502       Inf
sample estimates:
odds ratio 
       Inf 
```

Test the association between `hemoplasma` and `health_condition` in Cd:
```
table_health_condition_hemoplasma_Cd <- table(data_Cd$hemoplasma, data_Cd$health_condition)
table_health_condition_hemoplasma_Cd
fisher.test(table_health_condition_hemoplasma_Cd)
```

Results are:
```
     D  G
  0  1 14
  1  4 64
---
Fisher's Exact Test for Count Data
data:  table_health_condition_hemoplasma_Cd
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  0.02167699 12.75597732
sample estimates:
odds ratio 
  1.140982 
```

## Step 12. Impact of Hemoplasma_ infections on female reproductive status 
Test the association between `hemoplasma` and `female_reproductive_status` in Bt:
```
table_hemoplasma_infection_female_Bt <- table(data_Bt$hemoplasma, data_Bt$female_reproductive_status)
table_hemoplasma_infection_female_Bt
fisher.test(table_hemoplasma_infection_female_Bt)
```

Results are:
```
    Female lactating with a young Female non pregnant non lactating Pregnant female
  0                             8                                27               6
  1                             0                                 2               0
---
Fisher's Exact Test for Count Data
data:  table_hemoplasma_infection_female_Bt
p-value = 1
alternative hypothesis: two.sided
```

Test the association between `hemoplasma` and `female_reproductive_status` in Cd:
```
table_hemoplasma_infection_female_Cd <- table(data_Cd$hemoplasma, data_Cd$female_reproductive_status)
table_hemoplasma_infection_female_Cd
fisher.test(table_hemoplasma_infection_female_Cd)
```

Results are:
```
    Female lactating with a young Female non pregnant non lactating Pregnant female
  0                             2                                 7               1
  1                             5                                32               2
---
Fisher's Exact Test for Count Data
data:  table_hemoplasma_infection_female_Cd
p-value = 0.5051
alternative hypothesis: two.sided
```
