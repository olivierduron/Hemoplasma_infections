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


## Table of contents
- [Step 1. Retrieving the data](#step-1-retrieving-the-data)
- [Step 2. Prepare the data for analysis](#step-2-prepare-the-data-for-analysis)
- [Step 3. Calculate hemoplasma infection prevalence](#step-3-calculate-hemoplasma-infection-prevalence)
- [Step 4. Create the pathogens variable and sloth species subset](#step-4-create-the-pathogens-variable-and-sloth-species-subset)
- [Step 5. Impact of hemoplasma infections on Scale Mass Index (SMI) in adult Bt](#step-5-impact-of-hemoplasma-infections-on-scale-mass-index-smi-in-adult-bt)
- [Step 6. Impact of hemoplasma infections on Scale Mass Index (SMI) in adult Cd](#step-6-impact-of-hemoplasma-infections-on-scale-mass-index-smi-in-adult-cd)
- [Step 7. Impact of hemoplasma infections on neck circumference](#step-7-impact-of-hemoplasma-infections-on-neck-circumference)
- [Step 8. Impact of hemoplasma infections on hematocrit levels](#step-8-impact-of-hemoplasma-infections-on-hematocrit-levels)
- [Step 9. Impact of hemoplasma infections on body temperature](#step-9-impact-of-hemoplasma-infections-on-body-temperature)
- [Step 10. Impact of hemoplasma infections on general health_condition](#step-10-impact-of-hemoplasma-infections-on-general-health_condition)
- [Step 11. Impact of hemoplasma infections on female_reproductive_status](#step-11-impact-of-hemoplasma-infections-on-female_reproductive_status)

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

### Convert categorical variables into factors
```
data_hemoplasma$species        <- as.factor(data_hemoplasma$species)
data_hemoplasma$season         <- as.factor(data_hemoplasma$season)
data_hemoplasma$sex            <- as.factor(data_hemoplasma$sex)
data_hemoplasma$age            <- as.factor(data_hemoplasma$age)
data_hemoplasma$hemoplasma      <- as.factor(data_hemoplasma$hemoplasma)
data_hemoplasma$anaplasmataceae      <- as.factor(data_hemoplasma$anaplasmataceae)
data_hemoplasma$apicomplexa       <- as.factor(data_hemoplasma$apicomplexa)
```

### Load libraries for analysis
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
### Calculate `hemoplasma` infection prevalence and 95% confidence interval for _Bradypus tridactylus_ (Bt) and _Choloepus didactylus_ (Cd)

```
prevalence_results <- data_hemoplasma %>% group_by(species) %>% summarise(n = n(), positives = sum(hemoplasma == 1), prevalence = positives / n, conf_low = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$lower, conf_high = binom.confint(positives, n, conf.level = 0.95, methods = "exact")$upper)
print(prevalence_results)
```

Results :
| species | n | positives | prevalence | conf_low | conf_high |
|:--------|--:|----------:|-----------:|---------:|----------:|
| Bt      | 92 | 4  | 0.0435 | 0.0120 | 0.108 |
| Cd      | 83 | 68 | 0.8190 | 0.7200 | 0.895 |

### Test if `hemoplasma` is influenced by sloth `species`:
```
chisq.test(table(data_hemoplasma$hemoplasma, data_hemoplasma$species))
```

-> Results : `hemoplasma` prevalence differed strongly between the two sloth species, from 4.3% (4/92; 95% CI: 1.2–10.8%) in *Bradypus tridactylus* to 81.9% (68/83; 95% CI: 72.0–89.5%) in *Choloepus didactylus* (χ²₁ = 105.27, p < 2.2 × 10⁻¹⁶).

-> Interpretation : This strong interspecific difference suggests that host species and associated ecological or evolutionary traits may constrain `hemoplasma` infection.

## Step 4. Create the `pathogens` variable and sloth `species` subset

### Create the `pathogens` variable by merging `anaplasmataceae` and `apicomplexa` (0 = uninfected ; 1 = infected by `anaplasmataceae` and/or `apicomplexa`)
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

### Convert variables:
```
data_hemoplasma$weight <- as.numeric(data_hemoplasma$weight)
data_hemoplasma$total_length <- as.numeric(data_hemoplasma$total_length)
data_hemoplasma$wither_height <- as.numeric(data_hemoplasma$wither_height)
data_hemoplasma$neck_size <- as.numeric(data_hemoplasma$neck_size)
data_hemoplasma$temperature <- as.numeric(data_hemoplasma$temperature)
data_hemoplasma$hematocrit <- as.numeric(data_hemoplasma$hematocrit)
```

### Create a subset `data_Bt` containing only records for _Bradypus tridactylus_ (Bt)
```
data_Bt <- subset(data_hemoplasma, species == "Bt")
table(data_Bt$hemoplasma)
```

### Create a subset `data_Cd` containing only records for _Choloepus didactylus_ (Cd):
```
data_Cd <- subset(data_hemoplasma, species == "Cd")
table(data_Cd$hemoplasma)
```

## Step 5. Impact of `hemoplasma` infections on Scale Mass Index (SMI) in adult Bt
The Scaled Mass Index (SMI) was used as a body condition indicator that standardizes individual `weight` to `body_length`, using an allometric scaling relationship. SMI was calculated following Peig & Green (2009) (https://doi.org/10.1111/j.1600-0706.2009.17643.x).

### Function to calculate SMI for adult Bt
```
data_adult_Bt <- subset(data_Bt, age == "A")
sma_model_Bt <- sma(log(weight) ~ log(total_length), data = data_adult_Bt)
b <- coef(sma_model_Bt)[2]
L0 <- mean(data_adult_Bt$total_length, na.rm = TRUE)
data_adult_Bt$SMI <- data_adult_Bt$weight * (L0 / data_adult_Bt$total_length)^b
```

### Fit a GLM to test whether SMI is influenced by interactions among `hemoplasma`, `pathogens`, `sex` and `season` in Bt
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

### Fit a GLM to test whether SMI is influenced by additive effets of `hemoplasma`, `pathogens`, `sex` and `season` in Bt
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

### Model diagnostics (diagnostic plots, Shapiro–Wilk test, and Breusch–Pagan test)
```
par(mfrow = c(2, 2))
plot(model_SMIBt_final)
par(mfrow = c(1, 1))
shapiro.test(residuals(model_SMIBt_final))
library(lmtest)
bptest(model_SMIBt_final)
```

-> Results and interpretation :Residuals showed no evidence of departure from normality (Shapiro–Wilk: W = 0.988, p = 0.661) or heteroscedasticity (Breusch–Pagan: BP = 4.10, p = 0.129), supporting the use of the Gaussian model.


### Calculation of mean and standard error of SMI by `hemoplasma` infection status and `sex` for Bt
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

### Generate SMI chart for Bt
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
### Function to calculate SMI for adult Cd
```
data_adult_Cd <- subset(data_Cd, age == "A")
sma_model_Cd <- sma(log(weight) ~ log(total_length), data = data_adult_Cd)
b <- coef(sma_model_Cd)[2]
L0 <- mean(data_adult_Cd$total_length, na.rm = TRUE)
data_adult_Cd$SMI <- data_adult_Cd$weight * (L0 / data_adult_Cd$total_length)^b
```

### Fit a GLM to test whether SMI is influenced by interactions among `hemoplasma`, `pathogens`, `sex`, and `season` in Cd
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

### Model diagnostics (Shapiro–Wilk test)
```
shapiro.test(residuals(model_SMICd_null))
```
-> Results and interpretation : Residuals of the null model showed no evidence of departure from normality (Shapiro–Wilk: W = 0.979, p = 0.428).

### Generate SMI chart for Cd
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

### Fit a GLM to test whether neck circumference is influenced by interactions among `hemoplasma`, `pathogens`, `sex`, and `season` in Bt
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

### Model diagnostics (Shapiro–Wilk test)
```
shapiro.test(residuals(model_5_null))
```
-> Results and interpretation : Residuals showed no significant departure from normality (Shapiro–Wilk: W = 0.967, p = 0.093); this supports the Gaussian assumption.

### Fit a GLM to test whether neck circumference is influenced by interactions among `hemoplasma`, `pathogens`, `sex`, and `season` in Cd
```
model_6 <- glm(
  log(neck_size) ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

model_6_3way <- glm(
  log(neck_size) ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

model_6_2way <- glm(
  log(neck_size) ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

model_6_add <- glm(
  log(neck_size) ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

model_6_null <- glm(
  log(neck_size) ~ 1,
  data = data_adult_Cd,
  family = gaussian(link = "identity")
)

anova(model_6_3way, model_6, test = "Chisq")
anova(model_6_2way, model_6_3way, test = "Chisq")
anova(model_6_add, model_6_2way, test = "Chisq")
anova(model_6_null, model_6_add, test = "Chisq")

AIC(
  model_6,
  model_6_3way,
  model_6_2way,
  model_6_add,
  model_6_null
)
```

-> Results : For adult *Choloepus didactylus*, the four-way interaction model and the three-way interaction model were identical in fit (Δdeviance = 0, df = 0). Adding three-way interactions to the two-way model resulted in a small improvement in fit (LRT: χ²₂ = 0.093, *p* = 0.026), whereas adding two-way interactions to the additive model did not improve fit (χ²₆ = 0.012, *p* = 0.992). The additive model also did not improve on the null model (χ²₄ = 0.055, *p* = 0.370), and the null model had the lowest AIC (−69.52 vs. −66.07 for the additive model). The null model was therefore retained.

-> Interpretation : Neck size showed no overall detectable association with `hemoplasma` infection, other blood-borne `pathogens`, `season`, or `sex` in adult *Choloepus didactylus*, despite a nominal improvement in fit associated with three-way interactions.

### Model diagnostics (Shapiro–Wilk test)
```
shapiro.test(residuals(model_6_null))
```
-> Results and interpretation : Residuals of the null model did not significantly deviate from normality (Shapiro–Wilk test, W = 0.972, p = 0.314), supporting the use of a Gaussian error distribution.

## Step 8. Impact of `hemoplasma` infections on hematocrit levels
### Fit a GLM to test whether `hematocrit` is influenced by interactions among `hemoplasma`, `pathogens`, `sex`, and `season` in Bt
```
model_7 <- glm(
  hematocrit ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

model_7_3way <- glm(
  hematocrit ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

model_7_2way <- glm(
  hematocrit ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

model_7_add <- glm(
  hematocrit ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)

model_7_null <- glm(
  hematocrit ~ 1,
  data = data_adult_Bt,
  family = Gamma(link = "log")
)
anova(model_7_3way, model_7, test = "Chisq")
anova(model_7_2way, model_7_3way, test = "Chisq")
anova(model_7_add, model_7_2way, test = "Chisq")
anova(model_7_null, model_7_add, test = "Chisq")
AIC(model_7, model_7_3way, model_7_2way, model_7_add, model_7_null)
```
-> Results : For adult *Bradypus tridactylus*, adding three-way interactions did not improve model fit over the two-way interaction model (LRT: χ²₁ = 0.0001, *p* = 0.928), and the two-way interaction model did not improve fit over the additive model (χ²₄ = 0.099, *p* = 0.193). The additive model also did not improve on the null model (χ²₄ = 0.079, *p* = 0.321), and the null model had the lowest AIC (513.01 vs. 516.14 for the additive model). The null model was therefore retained.

-> Interpretation : `hematocrit` showed no detectable association with `hemoplasma` infection, other blood-borne `pathogens`, `season`, or `sex` in adult *Bradypus tridactylus*.

### Calculation of mean and standard error of `hematocrit` by `hemoplasma` for Bt
```
data_adult_Bt %>%
  group_by(hemoplasma) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results :
| Hemoplasma infection | Mean hematocrit (%) | SE |
|:---------------------|--------------------:|---:|
| Negative             | 39.0                | 0.57 |
| Positive             | 40.8                | 3.04 |


### Fit a GLM to test whether `hematocrit` is influenced by interactions among `hemoplasma`, `pathogens`, `sex`, and `season` in Cd
```
model_8 <- glm(
  hematocrit ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

model_8_3way <- glm(
  hematocrit ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

model_8_2way <- glm(
  hematocrit ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

model_8_add <- glm(
  hematocrit ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)

model_8_null <- glm(
  hematocrit ~ 1,
  data = data_adult_Cd,
  family = Gamma(link = "log")
)
anova(model_8_3way, model_8, test = "Chisq")
anova(model_8_2way, model_8_3way, test = "Chisq")
anova(model_8_add, model_8_2way, test = "Chisq")
anova(model_8_null, model_8_add, test = "Chisq")
AIC(model_8, model_8_3way, model_8_2way, model_8_add, model_8_null)
```

-> Results : For adult *Choloepus didactylus*, adding three-way interactions did not improve model fit over the two-way interaction model (LRT: χ²₁ = 0.003, *p* = 0.703), and the two-way interaction model did not improve fit over the additive model (χ²₆ = 0.079, *p* = 0.692). In contrast, the additive model significantly improved fit relative to the null model (χ²₄ = 0.209, *p* = 0.028) and had a lower AIC (383.99 vs. 385.81). The additive model was therefore retained.


### Fit a GLM to test whether `hematocrit` is influenced by `hemoplasma`, `pathogens`, `sex`, and `season` in Cd
```
model_8_final <- model_8_add
summary(model_8_final)
AIC(model_8_final)
```
-> Results : For adult *Choloepus didactylus*, the additive model was retained (AIC = 383.99; ΔAIC = 1.82 relative to the null model). `hematocrit` was higher during the wet `season` (β = 0.116 ± 0.041 SE, *p* = 0.006), but was not significantly associated with `hemoplasma` infection, other blood-borne `pathogens`, or `sex` (*p* > 0.11).

-> Interpretation : `hematocrit` showed a seasonal pattern in adult *Choloepus didactylus*, with higher values during the wet `season`, but no detectable association with `hemoplasma` infection.

### Calculation of mean and standard error of `hematocrit` by `hemoplasma` for Cd
```
data_adult_Cd %>%
  group_by(hemoplasma) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results :
| Hemoplasma infection | Mean hematocrit (%) | SE |
|:---------------------|--------------------:|---:|
| Negative             | 39.1                | 1.48 |
| Positive             | 38.6                | 0.80 |

### Calculation of mean and standard error of `hematocrit` by `season` for Cd
```
data_adult_Cd %>%
  group_by(season) %>%
  summarise(
    mean_hematocrit = mean(hematocrit, na.rm = TRUE),
    se_hematocrit = sd(hematocrit, na.rm = TRUE) / sqrt(sum(!is.na(hematocrit)))
  )
```

Results :
| Season | Mean hematocrit (%) | SE |
|:-------|--------------------:|---:|
| Dry    | 37.6                | 0.89 |
| Wet    | 41.1                | 0.94 |


Create violin plots for `hematocrit`
```
label_style <- element_text(size = 28, face = "bold")

pA <- ggplot(data_adult_Bt, aes(x = factor(hemoplasma, levels = c(0, 1),
                                labels = c("Uninfected", "Infected")),
                                y = hematocrit)) +
  geom_violin(fill = "darkorange2", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = NULL, y = "Hematocrit (%)", title = "A") +
  scale_y_continuous(limits = c(20, 60)) +
  scale_x_discrete(labels = NULL) +
  theme_minimal() +
  theme(plot.title = label_style)

pB <- ggplot(data_adult_Bt, aes(x = factor(season, levels = c("D", "W"),
                                labels = c("Dry", "Wet")),
                                y = hematocrit)) +
  geom_violin(fill = "darkorange2", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = NULL, y = NULL, title = "B") +
  scale_y_continuous(limits = c(20, 60)) +
  scale_x_discrete(labels = NULL) +
  theme_minimal() +
  theme(plot.title = label_style)

pC <- ggplot(data_adult_Bt, aes(x = factor(sex, levels = c("M", "F"),
                                labels = c("Male", "Female")),
                                y = hematocrit)) +
  geom_violin(fill = "darkorange2", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = NULL, y = NULL, title = "C") +
  scale_y_continuous(limits = c(20, 60)) +
  scale_x_discrete(labels = NULL) +
  theme_minimal() +
  theme(plot.title = label_style)

pD <- ggplot(data_adult_Cd, aes(x = factor(hemoplasma, levels = c(0, 1),
                                labels = c("Uninfected", "Infected")),
                                y = hematocrit)) +
  geom_violin(fill = "darkorange2", color = "black", alpha = 0.7, trim = FALSE) +
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
  geom_violin(fill = "darkorange2", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = "Season", y = NULL, title = "E") +
  scale_y_continuous(limits = c(10, 60)) +
  theme_minimal() +
  theme(plot.title = label_style)

pF <- ggplot(data_adult_Cd, aes(x = factor(sex, levels = c("M", "F"),
                                labels = c("Male", "Female")),
                                y = hematocrit)) +
  geom_violin(fill = "darkorange2", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(x = "Sex", y = NULL, title = "F") +
  scale_y_continuous(limits = c(10, 60)) +
  theme_minimal() +
  theme(plot.title = label_style)

final_plot <- (pA | pB | pC) / (pD | pE | pF)

print(final_plot)

ggsave(
  "hematocrit_by_species_and_predictors.png",
  plot = final_plot,
  width = 14,
  height = 12,
  units = "in",
  dpi = 300
)
```

## Step 9. Impact of `hemoplasma` infections on body `temperature`
### Convert `temperature` to numeric, handle left-censored values (<32°C) for analysis in Bt
```
data_adult_Bt <- data_adult_Bt %>%
  mutate(
    temperature_numeric = as.numeric(ifelse(temperature == "< 32.00", 32, temperature)),
    censored = ifelse(temperature == "< 32.00", TRUE, FALSE)
  )
```

### Create a left-censored Surv object (temp) for `temperature` in Bt
```
temp <- Surv(data_adult_Bt$temperature_numeric,
                  event = !data_adult_Bt$censored,
                  type = "left")
```

### Convert `temperature` to numeric, handle left-censored values (<32°C) for analysis in Cd
```
data_adult_Cd <- data_adult_Cd %>%
  mutate(
    temperature_numeric = as.numeric(ifelse(temperature == "< 32.00", 32, temperature)),
    censored = ifelse(temperature == "< 32.00", TRUE, FALSE)
  )
```

### Create a left-censored Surv object (temp) for `temperature` in Cd
```
temp <- Surv(data_adult_Cd$temperature_numeric,
                  event = !data_adult_Cd$censored,
                  type = "left")
```

### Fit Gaussian survival regression models to test the interaction effects among `hemoplasma`, `pathogens`, `season`, `sex` on `temperature` in Bt
```
model_10 <- survreg(
  temp ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Bt,
  dist = "gaussian"
)

model_10_3way <- survreg(
  temp ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Bt,
  dist = "gaussian"
)

model_10_2way <- survreg(
  temp ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Bt,
  dist = "gaussian"
)

model_10_add <- survreg(
  temp ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Bt,
  dist = "gaussian"
)

model_10_null <- survreg(
  temp ~ 1,
  data = data_adult_Bt,
  dist = "gaussian"
)

anova(model_10_3way, model_10)
anova(model_10_2way, model_10_3way)
anova(model_10_add, model_10_2way)
anova(model_10_null, model_10_add)

AIC(
  model_10,
  model_10_3way,
  model_10_2way,
  model_10_add,
  model_10_null
)
```

-> Results : For adult *Bradypus tridactylus*, the four-way and three-way interaction models did not differ in fit (LRT: χ²₁ = 0, *p* = 1.00). Adding three-way interactions to the two-way model did not improve fit (χ²₄ = 0.20, *p* = 0.995), nor did adding two-way interactions to the additive model (χ²₆ = 11.14, *p* = 0.084). However, the additive model improved fit relative to the null model (χ²₄ = 11.61, *p* = 0.020) and had the lowest AIC (43.40). The additive model was therefore retained.

### Fit Gaussian survival regression models to test the effects of `hemoplasma`, `pathogens`, `season`, `sex` on `temperature` in Bt
```
model_10_final <- model_10_add
summary(model_10_final)
```

-> Results : For adult *Bradypus tridactylus*, the additive model was retained (AIC = 43.40) and was significantly better than the null model (LRT: χ²₄ = 11.61, *p* = 0.020). Body `temperature` was higher during the wet `season` (β = 0.97 ± 0.40 SE, *p* = 0.016) and lower in `males` (β = −0.87 ± 0.39 SE, *p* = 0.025). `hemoplasma` infection was not associated with body `temperature` (β = −0.40 ± 0.78 SE, *p* = 0.612), nor were other blood-borne `pathogens` (β = 0.59 ± 0.39 SE, *p* = 0.128).

-> Interpretation : Body `temperature` showed significant seasonal and `sex`-related variation in adult *Bradypus tridactylus*, but no detectable association with `hemoplasma` infection.

### Fit Gaussian survival regression models to test the interaction effects among `hemoplasma`, `pathogens`, `season`, `sex` on `temperature` in Cd
```
model_11 <- survreg(
  temp ~ hemoplasma * pathogens * season * sex,
  data = data_adult_Cd,
  dist = "gaussian"
)

model_11_3way <- survreg(
  temp ~ (hemoplasma + pathogens + season + sex)^3,
  data = data_adult_Cd,
  dist = "gaussian"
)

model_11_2way <- survreg(
  temp ~ (hemoplasma + pathogens + season + sex)^2,
  data = data_adult_Cd,
  dist = "gaussian"
)

model_11_add <- survreg(
  temp ~ hemoplasma + pathogens + season + sex,
  data = data_adult_Cd,
  dist = "gaussian"
)

model_11_null <- survreg(
  temp ~ 1,
  data = data_adult_Cd,
  dist = "gaussian"
)

anova(model_11_3way, model_11)
anova(model_11_2way, model_11_3way)
anova(model_11_add, model_11_2way)
anova(model_11_null, model_11_add)

AIC(
  model_11,
  model_11_3way,
  model_11_2way,
  model_11_add,
  model_11_null
)
```

-> Results : For adult *Choloepus didactylus*, none of the interaction models improved model fit (three-way vs. two-way: LRT, χ²₄ = 0, *p* = 1.00; two-way vs. additive: χ²₆ = 3.55, *p* = 0.737). The additive model also did not improve on the null model (χ²₄ = 0.56, *p* = 0.967), and the null model had the lowest AIC (62.74 vs. 70.17 for the additive model). The null model was therefore retained.

-> Interpretation : Body `temperature` showed no detectable association with `hemoplasma` infection, other blood-borne `pathogens`, `season`, or `sex` in adult *Choloepus didactylus*.

## Step 10. Impact of `hemoplasma` infections on general `health_condition` 
### Test the association between `hemoplasma` and `health_condition` in Bt
```
table_health_condition_hemoplasma_Bt <- table(data_Bt$hemoplasma, data_Bt$health_condition)
table_health_condition_hemoplasma_Bt
fisher.test(table_health_condition_hemoplasma_Bt)
```
-> Results : In adult *Bradypus tridactylus*, `hemoplasma` infection was not associated with `health_condition` (Fisher’s exact test, p = 1.00). All four `hemoplasma`-positive individuals were classified as having a good `health_condition`.

-> Interpretation : There was no detectable association between `hemoplasma` infection and `health_condition` in adult *Bradypus tridactylus*. 

### Test the association between `hemoplasma` and `health_condition` in Cd
```
table_health_condition_hemoplasma_Cd <- table(data_Cd$hemoplasma, data_Cd$health_condition)
table_health_condition_hemoplasma_Cd
fisher.test(table_health_condition_hemoplasma_Cd)
```
-> Results :In adult *Choloepus didactylus*, `hemoplasma` infection was not associated with `health_condition` (Fisher’s exact test, *p* = 1.00). Among `hemoplasma`-positive individuals, 64/68 (94.1%) were classified as having a good `health_condition`, compared with 14/15 (93.3%) among `hemoplasma`-negative individuals.

-> Interpretation : There was no detectable association between `hemoplasma` infection and `health_condition` in adult *Choloepus didactylus*.


## Step 11. Impact of `hemoplasma` infections on `female_reproductive_status`
### Test the association between `hemoplasma` and `female_reproductive_status` in Bt
```
table_hemoplasma_infection_female_Bt <- table(data_Bt$hemoplasma, data_Bt$female_reproductive_status)
table_hemoplasma_infection_female_Bt
fisher.test(table_hemoplasma_infection_female_Bt)
```
-> Results : In adult *Bradypus tridactylus* females, `hemoplasma` infection was not associated with `female_reproductive_status` (Fisher’s exact test, *p* = 1.00). The two `hemoplasma`-positive females were both classified as non-pregnant and non-lactating, while no infections were detected among lactating or pregnant females.

-> Interpretation: There was no detectable association between `hemoplasma` infection and `female_reproductive_status` in adult *Bradypus tridactylus*. However, the very small number of `hemoplasma`-positive females (n = 2) limits the power of this comparison.

### Test the association between `hemoplasma` and `female_reproductive_status` in Cd
```
table_hemoplasma_infection_female_Cd <- table(data_Cd$hemoplasma, data_Cd$female_reproductive_status)
table_hemoplasma_infection_female_Cd
fisher.test(table_hemoplasma_infection_female_Cd)
```
-> Results : In adult *Choloepus didactylus* females, `hemoplasma` infection was not associated with `female_reproductive_status` (Fisher’s exact test, *p* = 0.505). Among `hemoplasma`-positive females, 5/39 (12.8%) were lactating with a young, 32/39 (82.1%) were non-pregnant and non-lactating, and 2/39 (5.1%) were pregnant. The corresponding proportions among `hemoplasma`-negative females were 2/10 (20.0%), 7/10 (70.0%), and 1/10 (10.0%), respectively.

-> Interpretation : There was no detectable association between `hemoplasma` infection and `female_reproductive_status` in adult *Choloepus didactylus*.
