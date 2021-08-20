setwd("~/Downloads/") #Set working directory
rm(list = ls()) #Clear working environment
graphics.off() #Clear graphics


###########Overall review of projects###########
###Load packages
library(dplyr) #Load dplyr package for data wrangling 
library(plotrix) #Load plotrix package for calculating standard error
library(ggplot2) #Load ggplot2 package for plotting
library(ggpubr) #Load ggpubr package for arranging multiple ggplots into one figure
library(forcats) #Load package for ggplot x-axis ordering
library(sf) #Load sf package for converting objects to sf object
library(rnaturalearth) #Load rnaturalearth package for world map plotting
library(rnaturalearthdata) #Load rnaturalearthdata package for world map plotting


###Summarising data structure
data <- read.csv("project_study.csv") #Load in dataset
data %>% group_by(continent) %>% summarise(n()) #Number of projects in each continent
data %>% group_by(study_design) %>% summarise(n()) #Number of projects per study design
data %>% filter(publication!="") %>% group_by(publication) %>% summarise(n()) #Number of grey and primary literature studies
data %>% group_by(species) %>% summarise(n()) #Target species of restoration project
data %>% group_by(restoration) %>% summarise(n()) #Number of projects of each restoration type
mean(data$monitoring, na.rm = TRUE) #Mean post-monitoring years
std.error(data$monitoring, na.rm = TRUE) #Standard error of mean post-monitoring years
range(data$monitoring, na.rm = TRUE) # Range of post-monitoring years
sum(data$monitoring<5, na.rm = TRUE) #Number of projects monitoring for less than 5 years
sum(data$monitoring>=10, na.rm = TRUE) #Number of projects monitoring for 10 or more years

###Plotting Map of Sites
coord <- read.csv("study_coordinates.csv") #Read in dataset with latitude and longitude of each project
world <- ne_countries(scale = "medium", returnclass = "sf") #Load world map
coord_sf <- st_as_sf(coord, coords = c("Longitude", "Latitude"),
                     crs = 4236, agr = "constant") #Convert dataframe to sf object using WGS84 projection
map_plot <- ggplot(data = world) +
  geom_sf() +
  geom_sf(data = coord_sf, colour = "black", size = 1.5, pch = 21, fill = "red") +
  theme_classic() #Plot world map with sites
map_plot

###Plotting restoration studies
levels(data$restoration) <- gsub(" ", "\n", levels(data$restoration)) # Add line breaks to x-axis labels in ggplot
restor_plot <- ggplot(data, aes(x = fct_infreq(restoration), fill = species)) + 
  geom_bar(stat = "count", width = 0.8) + 
  theme_classic() + 
  labs(x = "Restoration Action", y = "Number of Projects") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13, vjust = -.75),
        legend.title = element_blank()) + 
  scale_y_continuous(limit = c(0,25),
    expand = c(0,0)) #Barplot of number of restoration projects, coloured by species
restor_plot
map_restor <- ggarrange(map_plot, restor_plot, labels = c("A", "B"),
                        ncol = 1, nrow = 2) #Combine map plot and restoration plot
map_restor
ggexport(map_restor, filename = "figure1.pdf") #Export map_restor as pdf

###########Meta-Analysis###########
rm(list = ls()) #Clear working environment

###Load packages
library(metapower) #Load metapower package for statistical power analysis
library(metafor) #Load metafor package for meta-analysis
library(dmetar) #Load dmetar package for statistical power analysis
library(assertive.base) #Load package for parentheses function
library(psych) #Load package for co-variate correlation analysis

###Load dataset
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis effect size dataset
MA_data <- subset(MA_data, after_sample_size!="1") #Remove projects with one sample size, cannot compute standardised mean difference effect sizes with sample size of one
MA_data$publication_year <- as.numeric(as.character(MA_data$publication_year)) #Set publication yeaar variable to numeric

###Power analysis
n1 <- mean(MA_data$before_sample_size) #Mean sample size of control group
n2 <- mean(MA_data$after_sample_size) #Mean sample size of treatment group
power.analysis(d = 0.85, 
               k = 94, 
               n1 = n1, 
               n2 = n2, 
               p = 0.05,
               heterogeneity = "high") #Power analysis of random-effects model, heterogeneity is assumed to be high in ecological meta-analyses

meta_subgroup_power <- subgroup_power(n_groups = X, effect_sizes = c(X, X), study_size = 6.75, 
                                      i2 = 0.75, k = 73, es_type = "d") #statistical power of subgroup differences
summary(meta_subgroup_power) #summary of power analysis of subgroup differences

###Calculate effect sizes
effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       data = MA_data) #Calculate standardised mean difference effect sizes 

###Meta-analysis of all data
random_effect_model_SMD <- rma.mv(yi, 
                                  vi, 
                                  random = ~ 1| study_ID, #Study ID as random effect
                                  data = effect_sizes) #Multi-level random-effects model with REML estimator
summary.rma(random_effect_model_SMD) #Print model results
effect_sizes %>% filter(yi>0) %>% summarise(n()) #Number of positive effect sizes
effect_sizes %>% filter(yi<0) %>% summarise(n()) #Number of negative effect sizes
confid_interval <- summary.escalc(effect_sizes) #Extract confidence intervals of each effect size
confid_interval %>% filter(zi>1.96) %>% summarise(n()) #Number of significantly positive effect sizes
confid_interval %>% filter(ci.lb<0) %>% filter(ci.ub<0) %>% summarise(n()) #Number of significantly negative effect sizes
res <- rma.uni(yi, vi, data = effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res) #Print model results
qqnorm.rma.uni(res, envelope = FALSE) #Model validation - normality of residuals

###Sensitivity analysis - all effect sizes
#Standardized residuals
rs_overall <- rstudent.rma.mv(random_effect_model_SMD) #Standardized residuals
rs_overall <- as.data.frame(rs_overall) #Set as dataframe 
limit_rs <- 3  #Set limit of standardised deleted residuals to 3
rs_overall$ES_ID <- as.numeric(1:nrow(rs_overall)) #Add a new column with the effect size ID number

#Hat values
hat_overall <- hatvalues.rma.mv(random_effect_model_SMD) #Hat values
hat_overall <- as.data.frame(hat_overall) #Set as dataframe
mean_overall <- mean(hat_overall$hat) #Average hat values
limit_hat_overall <- 2 * mean_overall #Set limit of hat values to twice the average hat values
hat_overall$ES_ID <- as.numeric(1:nrow(hat_overall)) #Add a new column with the effect size ID number

effect_sizes$ES_ID <- as.numeric(1:nrow(effect_sizes)) #Add a new column with the effect size ID number
sensitivity_res <- left_join(hat_overall,rs_overall, by="ES_ID")%>%
  left_join(effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs) %>%
  filter(hat_overall < limit_hat_overall) #Remove outliers

sensitivity_overall <- rma.mv(yi, #Effect size from each row in database
                              vi, #Measure of variance from each row in database
                              random = ~ 1| study_ID, #Study ID as random effect
                              data = sensitivity_res) #Multi-level random-effects model with REML estimator excluding outliers
summary.rma(sensitivity_overall) #Print model results without outliers
summary(random_effect_model_SMD) #Print model with outliers


###Assessment of publication bias
primary <- rma.mv(yi, 
                  vi, 
                  subset = (publication_type=="Journal"),
                  random = ~ 1| study_ID, #Study ID as random effect
                  data = effect_sizes) #Multi-level random-effects model with REML estimator of primary literature
summary.rma(primary) #Print model results
grey <- rma.mv(yi, 
                  vi, 
                  subset = (publication_type=="Grey"),
                  random = ~ 1| study_ID, #Study ID as random effect
                  data = effect_sizes) #Multi-level random-effects model with REML estimator of grey literature
summary.rma(grey) #Print model results
pub_type <- rma.mv(yi, 
                   vi, 
                   mods = ~ publication_type, 
                   random = ~ 1| study_ID,
                   data = effect_sizes) #Multi-level random-effects model with REML with publication type as a moderator variable to assess for bias
summary.rma(pub_type) #Print model results, no difference between effect sizes of grey and primary literature

funnel(random_effect_model_SMD,
       level=c(95, 99), 
       shade=c("white", "gray55", "gray85"),
       refline = 0,
       legend=TRUE, 
       back = "grey85") #Contour-enhanced funnel plot to assess for publication bias, with 95% and 99% significance levels

effect_sizes$se <- confid_interval$sei #Add standard error to effect sizes dataframe
effect_sizes$seinvers <- 1/effect_sizes$se #Inverse of standard error
egger_test_se <- rma.mv(yi, 
                     vi, 
                     mods = ~ se, 
                     random = ~ 1| study_ID,
                     data = effect_sizes) #Adjusted Egger's regression test to test for funnel asymmetry with standard error as moderator
summary.rma(egger_test_se) #Print model results, #Intercept significantly different from zero
#Slope of B1 is significant, so use vi as moderator 
egger_test_vi <- rma.mv(yi,
                     vi, 
                     mods = ~ vi, 
                     random = ~ 1| study_ID,
                     data = effect_sizes) #Adjusted Egger's regression test to test for funnel asymmetry with sampling variance as moderator
summary.rma(egger_test_vi) #Print model results, #Intercept significantly different from zero
#Slope of B1 is significant, so use sample size as moderator 
effect_sizes$sample_size <- (4*effect_sizes$before_sample_size*effect_sizes$after_sample_size) / (effect_sizes$before_sample_size +effect_sizes$after_sample_size) #Calculate effective sample size as per methods of Nakagawa et al. (2021)
effect_sizes$sample_size <- effect_sizes$sample_size / 4 #Calculate n1
effect_sizes$sample_size_invers <- sqrt((1/effect_sizes$sample_size)) #Calculate sqrt of inverse  sample size
egger_test_sample_size <- rma.mv(yi, 
                     vi, 
                     mods = ~ sample_size_invers, 
                     random = ~ 1| study_ID,
                     data = effect_sizes) #Adjusted Egger's regression test to test for funnel asymmetry with inverse sample size as moderator
summary.rma(egger_test_sample_size) #Print model results, #Intercept does not statistically differ from zero
funnel(egger_test_sample_size,
       level=c(95, 99), 
       shade=c("white", "gray55", "gray85"),
       legend=TRUE, 
       refline = 0,
       back = "grey85") #Funnel plot from adjusted Egger's regression test using inverse sample size

time_lag_bias <- rma.mv(yi, 
                        vi, 
                        mods = ~ publication_year, 
                        random = ~ 1| study_ID,
                        data = effect_sizes) #Multi-level random-effects model with REML with publication year as a moderator variable to assess for time-lag bias
summary.rma(time_lag_bias) #Print model results, #No indication of time-lag bias
cor.test(effect_sizes$post_monitoring_years, 
         effect_sizes$publication_year, method = "kendall") #Test for correlation between project age (post monitoring years) and publication year

random_effect_model_species <- rma.mv(yi,
                                      vi, 
                                      mods = ~ species, 
                                      random = ~ 1| study_ID,
                                      data = effect_sizes) #Multi-level random-effects model with REML with species as a moderator variable to test for differences between species
summary.rma(random_effect_model_species) #Print model results, no difference between species

###Sensitivity analysis of species
#Standardized residuals
rs_species <- rstudent.rma.mv(random_effect_model_species)
rs_species <- as.data.frame(rs_species)
limit_rs <- 3
rs_species$ES_ID <- as.numeric(1:nrow(rs_species)) #add a new column with the effect size ID number

#Hat values
hat_species <- hatvalues.rma.mv(random_effect_model_species)
hat_species <- as.data.frame(hat_species)
mean_species <- mean(hat_species$hat)
limit_hat_species <- 2 * mean_species
hat_species$ES_ID <- as.numeric(1:nrow(hat_species)) #add a new column with the effect size ID number

effect_sizes$ES_ID <- as.numeric(1:nrow(effect_sizes))
sensitivity_species <- left_join(hat_species,rs_species, by="ES_ID")%>%
  left_join(effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs) %>%
  filter(hat_species < limit_hat_species)

sensitivity_species_res <- rma.mv(yi, 
                              vi, 
                              mods = ~ factor(species),
                              random = ~ 1| study_ID, #Study ID as random effect
                              data = sensitivity_species) #Multivariate meta-analysis random-effects model with REML
summary.rma(sensitivity_species_res) #Print model results without outliers
summary(random_effect_model_species) #Print model with outliers


###Salmon and trout separate meta-analysis
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis data 
MA_data <- MA_data %>% filter(before_sample_size!=1) #Remove projects with sample size of 1
salmon <- subset(MA_data, species=="Salmon") #Subset salmon data
salmon_effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       data = salmon) #Calculate standardised mean difference effect sizes 

random_effect_model_SMD_Salmon <- rma.mv(yi, 
                                         vi, 
                                         random = ~ 1 | study_ID,
                                         method = "REML",
                                         data = salmon_effect_sizes) #Multi-level random-effects model with REML estimator 
summary.rma(random_effect_model_SMD_Salmon) #Print model results
salmon_effect_sizes %>% filter(yi>0) %>% summarise(n()) #Number of positive effect sizes
salmon_effect_sizes %>% filter(yi<0) %>% summarise(n()) #Number of negative effect sizes
con_interval_salmon <- summary.escalc(salmon_effect_sizes) #Extract confidence intervals of each effect size
con_interval_salmon %>% filter(ci.lb>0) %>% filter(ci.ub>0) %>% summarise(n()) #Number of significantly positive effect sizes
con_interval_salmon %>% filter(ci.lb<0) %>% filter(ci.ub<0) %>% summarise(n()) #Number of significantly negative effect sizes
res_salmon <- rma.uni(yi, vi, data = salmon_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_salmon) #Print model results
qqnorm.rma.uni(res_salmon, envelope = FALSE) #Model validation - normality of residuals

forest(random_effect_model_SMD_Salmon, 
       header = c("Author(s) and Year"),
       slab = salmon_effect_sizes$author, #Label individual studies
       mlab = "RE Model",
       annotate = TRUE, #Display individual effect sizes and 95% confidence intervals
       xlab = "Standardised Mean Difference", # Label for x-axis
       cex = 1.25, # Text side for study labels
       pch = 16, # Shape of effect size in forest plot
       psize = 0.75, #Size of individual effect size bars
       cex.lab = 1.5, # Size of x-axis label
) #Forest plot with overall  effect size

trout <- subset(MA_data, species=="Trout") #Subset trout data
trout_effect_sizes <- escalc("SMD", 
                                 m2i = before_mean,  
                                 sd2i = before_SD,  
                                 n2i = before_sample_size, 
                                 m1i = after_mean, 
                                 sd1i = after_SD, 
                                 n1i = after_sample_size, 
                                 data = trout) #Calculate standardised mean difference effect sizes 
random_effect_model_SMD_Trout <- rma.mv(yi, 
                                        vi,
                                        random = ~ 1 | study_ID,
                                        method = "REML",
                                        data = trout_effect_sizes) #Multi-level random-effects model with REML estimator 
summary.rma(random_effect_model_SMD_Trout) #Print model results
trout_effect_sizes %>% filter(yi>0) %>% summarise(n()) #Number of positive effect sizes
trout_effect_sizes %>% filter(yi<0) %>% summarise(n()) #Number of negative effect sizes
con_interval_trout <- summary.escalc(trout_effect_sizes) #Extract confidence intervals of each effect size
con_interval_trout %>% filter(ci.lb>0) %>% filter(ci.ub>0) %>% summarise(n()) #Number of significantly positive effect sizes
con_interval_trout %>% filter(ci.lb<0) %>% filter(ci.ub<0) %>% summarise(n()) #Number of significantly negative effect sizes
res_trout <- rma.uni(yi, vi, data = trout_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_trout) #Print model results
qqnorm.rma.uni(res_trout, envelope = FALSE) #Model validation - normality of residuals

forest(random_effect_model_SMD_Trout, 
       header = c("Author(s) and Year"),
       top = 0.75,
       slab = trout_effect_sizes$author, #Label individual studies
       mlab = "RE Model",
       annotate = TRUE, #Display individual effect sizes and 95% confidence intervals
       xlab = "Standardised Mean Difference", # Label for x-axis
       cex = 1, # Text side for study labels
       pch = 16, # Shape of effect size in forest plot
       psize = 0.75,
       cex.lab = 1.5, # Size of x-axis label
) #Forest plot with overall grand effect size

species_overall_ES <- data.frame(species = c("Salmon", "Trout"),
                                 k = c(random_effect_model_SMD_Salmon[["k"]], 
                                       random_effect_model_SMD_Trout[["k"]]),
                                 yi = c(random_effect_model_SMD_Salmon[["b"]], 
                                        random_effect_model_SMD_Trout[["b"]]),
                                 ci.lb = c(random_effect_model_SMD_Salmon[["ci.lb"]], 
                                           random_effect_model_SMD_Trout[["ci.lb"]]),
                                 ci.ub = c(random_effect_model_SMD_Salmon[["ci.ub"]], 
                                           random_effect_model_SMD_Trout[["ci.ub"]]),
                                 p = c(random_effect_model_SMD_Salmon[["pval"]], 
                                       random_effect_model_SMD_Trout[["pval"]])) #Salmon and trout number of studies, overall effect size, confidence intervals, p-value
ggplot(species_overall_ES, 
       aes(x = species, y = yi)) +  
  geom_bar(stat = "identity", width = 0.5,
           fill = "grey55") +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.1) +
  theme_classic() +
  xlab("Species") +
  ylab("Effect Size (95% CI)") +
  theme(axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0,0),
                     limit = c(0,1.5)) #Plot salmon and trout overall effect sizes

###Species-specifc responses to restoration sub-groups
salmon %>% group_by(restoration) %>% summarise(n()) #Number of projects per restoration category in salmon dataset
#Minimum sample size of 10 projects for restoration sub-group analysis, analyse salmon response to in-stream structures and liming

Instream_Salmon <- rma.mv(yi,
                          vi,
                          method = "REML", 
                          subset = (restoration=="In-stream_structure"),
                          random = ~ 1| study_ID,
                          data = salmon_effect_sizes) #Multi-level random-effects model with REML for in-stream structures
summary.rma(Instream_Salmon) #Print model results
Liming_Salmon <- rma.mv(yi, 
                        vi, 
                        method = "REML",
                        subset = (restoration=="Liming"),
                        random = ~ 1| study_ID,
                        data = salmon_effect_sizes) #Multi-level random-effects model with REML for liming
summary.rma(Liming_Salmon) #Print model results
Salmon_restoration <- data.frame(species = "Salmon",
                                 restoration = c("Barrier Removal",  
                                                 "In-stream Structures",
                                                 "Liming",  
                                                 "Spawning Gravel"),
                                 k = c(4,
                                       Instream_Salmon[["k"]],
                                       Liming_Salmon[["k"]],
                                       2),
                                 yi = c("NA",
                                        Instream_Salmon[["b"]],
                                        Liming_Salmon[["b"]],
                                        "NA"),
                                 ci.lb = c("NA",
                                           Instream_Salmon[["ci.lb"]],
                                           Liming_Salmon[["ci.lb"]],
                                           "NA"),
                                 ci.ub = c("NA",
                                           Instream_Salmon[["ci.ub"]],
                                           Liming_Salmon[["ci.ub"]],
                                           "NA"),
                                 p = c("NA",
                                       Instream_Salmon[["pval"]],
                                       Liming_Salmon[["pval"]],
                                       "NA")) #Dataframe of salmon responses to restoration type

trout %>% group_by(restoration) %>% summarise(n()) #Number of projects per restoration category in trout dataset
#Minimum sample size of 10 projects for restoration sub-group analysis, analyse salmon response to barrier removal, in-stream structures and spawning gravel

Barrier_Removal_Trout <- rma.mv(yi, 
                                vi, 
                                method = "REML", 
                                subset = (restoration=="Barrier_removal"),
                                random = ~ 1| study_ID,
                                data = trout_effect_sizes) #Multi-level random-effects model with REML for barrier removal
summary.rma(Barrier_Removal_Trout) #Print model results
Instream_Trout <- rma.mv(yi, # effect size from each row in database
                         vi, # measure of variance from each row in database
                         method = "REML", 
                         subset = (restoration=="In-stream_structure"),
                         random = ~ 1| study_ID,
                         data = trout_effect_sizes) #Multi-level random-effects model with REML for in-stream structures
summary.rma(Instream_Trout) #Print model results
SpawningGravel_Trout <- rma.mv(yi, # effect size from each row in database
                               vi, # measure of variance from each row in database
                               method = "REML", 
                               subset = (restoration=="Spawning_gravel"),
                               random = ~ 1| study_ID,
                               data = trout_effect_sizes) #Multi-level random-effects model with REML for spawning gravel
summary.rma(SpawningGravel_Trout) #Print model results
Trout_restoration <- data.frame(species = "Trout",
                                restoration = c("Barrier Removal",  
                                                "In-stream Structures",
                                                "Liming",  
                                                "Spawning Gravel"),
                                k = c(Barrier_Removal_Trout[["k"]],
                                      Instream_Trout[["k"]],
                                      2,
                                      SpawningGravel_Trout[["k"]]),
                                yi = c(Barrier_Removal_Trout[["b"]],
                                       Instream_Trout[["b"]],
                                       "NA",
                                       SpawningGravel_Trout[["b"]]),
                                ci.lb = c(Barrier_Removal_Trout[["ci.lb"]],
                                          Instream_Trout[["ci.lb"]],
                                          "NA",
                                          SpawningGravel_Trout[["ci.lb"]]),
                                ci.ub = c(Barrier_Removal_Trout[["ci.ub"]],
                                          Instream_Trout[["ci.ub"]],
                                          "NA",
                                          SpawningGravel_Trout[["ci.ub"]]),
                                p = c(Barrier_Removal_Trout[["pval"]],
                                      Instream_Trout[["pval"]],
                                      "NA",
                                      SpawningGravel_Trout[["pval"]])) #Dataframe of trout responses to restoration type
Salmon_Trout_Effect_Size_Restoration <- rbind(Salmon_restoration, Trout_restoration) #Combine salmon and trout dataframes

###Plot salmon and trout responses to restoration type 
Salmon_Trout_Effect_Size_Restoration$yi <- as.numeric(as.character(Salmon_Trout_Effect_Size_Restoration$yi)) #Set as numeric
Salmon_Trout_Effect_Size_Restoration$ci.lb <- as.numeric(as.character(Salmon_Trout_Effect_Size_Restoration$ci.lb)) #Set as numeric
Salmon_Trout_Effect_Size_Restoration$ci.ub <- as.numeric(as.character(Salmon_Trout_Effect_Size_Restoration$ci.ub)) #Set as numeric
levels(Salmon_Trout_Effect_Size_Restoration$restoration) <- 
  gsub(" ", "\n", levels(Salmon_Trout_Effect_Size_Restoration$restoration)) # Add line breaks
restor_species_plot <- 
  ggplot(Salmon_Trout_Effect_Size_Restoration, 
         aes(x = restoration, y = yi, ymin = ci.lb, ymax = ci.ub)) +  
  geom_pointrange(aes(color = species, group = species), 
                  position = position_dodge(width = 0.5)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  scale_y_continuous(breaks = c(-1,0,1,2,3), 
                     limit = c(-1,3), expand = c(0,0)) +
  xlab("Restoration Type") +
  ylab("Effect Size (95% CI)") +
  theme(axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 14, vjust = -0.75),
        axis.title.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()) #Point range plot of species-specific responses to restoration type
restor_species_plot



###Restoration sub-group meta-analysis (for both salmonid species)
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis effect size dataset
MA_data <- subset(MA_data, after_sample_size!="1") #Remove projects with one sample size, cannot compute standardised mean difference effect sizes with sample size of one
MA_data$publication_year <- as.numeric(as.character(MA_data$publication_year)) #Set publication yeaar variable to numeric
MA_data %>% group_by(restoration) %>% summarise(n()) #Number of projects per restoration category
#Minimum sample size of 10 projects for restoration sub-group analysis
#Analyse barrier removal, in-stream structures, liming and spawning gravel

effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       data = MA_data) #Calculate standardised mean difference effect sizes 
Barrier_Removal_random_effect_model <- rma.mv(yi, 
                                              vi, 
                                              method = "REML",
                                              subset = (restoration=="Barrier_removal"),
                                              random = ~ 1| study_ID,
                                              data = effect_sizes) #Multi-level random-effects model with REML for barrier removal
summary.rma(Barrier_Removal_random_effect_model) #Print model results
Channel_random_effect_model <- rma.mv(yi, 
                                      vi, 
                                      method = "REML",
                                              subset = (restoration=="Channel"),
                                              random = ~ 1| study_ID,
                                              data = effect_sizes) #Multi-level random-effects model with REML for channel restoration
summary.rma(Channel_random_effect_model) #Print model results
Instream_random_effect_model <- rma.mv(yi, 
                                       vi, 
                                       method = "REML",
                                       subset = (restoration=="In-stream_structure"),
                                       random = ~ 1| study_ID,
                                       data = effect_sizes) #Multi-level random-effects model with REML for in-stream structures
summary.rma(Instream_random_effect_model) #Print model results
Liming_random_effect_model <- rma.mv(yi, 
                                     vi, 
                                     method = "REML",
                                     subset = (restoration=="Liming"),
                                     random = ~ 1| study_ID,
                                     data = effect_sizes) #Multi-level random-effects model with REML for liming
summary.rma(Liming_random_effect_model) #Print model results
Nutrient_random_effect_model <- rma.mv(yi, 
                                       vi, 
                                       method = "REML",
                                       subset = (restoration=="Nutrient"),
                                       random = ~ 1| study_ID,
                                       data = effect_sizes) #Multi-level random-effects model with REML for nutrient addition
summary.rma(Nutrient_random_effect_model) #Print model results
Riparian_random_effect_model <- rma.mv(yi, 
                                     vi, 
                                     method = "REML",
                                     subset = (restoration=="Riparian"),
                                     random = ~ 1| study_ID,
                                     data = effect_sizes) #Multi-level random-effects model with REML for riparian restoration
summary.rma(Riparian_random_effect_model) #Print model results
SpawningGravel_random_effect_model <- rma.mv(yi, 
                                             vi, 
                                             method = "REML",
                                             subset = (restoration=="Spawning_gravel"),
                                             random = ~ 1| study_ID,
                                             data = effect_sizes) #Multi-level random-effects model with REML for spawning gravel
summary.rma(SpawningGravel_random_effect_model) #Print model results

MA_restoration <- data.frame(restoration = c("Barrier Removal", "Channel Restoration", 
                                             "In-stream Structures",
                                             "Liming", "Nutrient Restoration", 
                                             "Spawning Gravel"),
                             k = c(Barrier_Removal_random_effect_model[["k"]], 
                                   Channel_random_effect_model[["k"]],
                                   Instream_random_effect_model[["k"]],
                                   Liming_random_effect_model[["k"]],
                                   Nutrient_random_effect_model[["k"]],
                                   SpawningGravel_random_effect_model[["k"]]),
                             yi = c(Barrier_Removal_random_effect_model[["b"]], 
                                    Channel_random_effect_model[["b"]],
                                    Instream_random_effect_model[["b"]],
                                    Liming_random_effect_model[["b"]],
                                    Nutrient_random_effect_model[["b"]],
                                    SpawningGravel_random_effect_model[["b"]]),
                             ci.lb = c(Barrier_Removal_random_effect_model[["ci.lb"]], 
                                       Channel_random_effect_model[["ci.lb"]],
                                       Instream_random_effect_model[["ci.lb"]],
                                       Liming_random_effect_model[["ci.lb"]],
                                       Nutrient_random_effect_model[["ci.lb"]],
                                       SpawningGravel_random_effect_model[["ci.lb"]]),
                             ci.ub = c(Barrier_Removal_random_effect_model[["ci.ub"]], 
                                       Channel_random_effect_model[["ci.ub"]],
                                       Instream_random_effect_model[["ci.ub"]],
                                       Liming_random_effect_model[["ci.ub"]],
                                       Nutrient_random_effect_model[["ci.ub"]],
                                       SpawningGravel_random_effect_model[["ci.ub"]]),
                             p = c(Barrier_Removal_random_effect_model[["pval"]], 
                                   Channel_random_effect_model[["pval"]],
                                   Instream_random_effect_model[["pval"]],
                                   Liming_random_effect_model[["pval"]],
                                   Nutrient_random_effect_model[["pval"]],
                                   SpawningGravel_random_effect_model[["pval"]])) #Dataframe of salmonid responses to restoration type (overall effect size, number of effect sizes, confidence intervals, p-value)

MA_restoration <- data.frame(restoration = c("Barrier Removal",  
                                             "In-stream Structures",
                                             "Liming",  
                                             "Spawning Gravel"),
                             k = c(Barrier_Removal_random_effect_model[["k"]],
                                   Instream_random_effect_model[["k"]],
                                   Liming_random_effect_model[["k"]],
                                   SpawningGravel_random_effect_model[["k"]]),
                             yi = c(Barrier_Removal_random_effect_model[["b"]],
                                    Instream_random_effect_model[["b"]],
                                    Liming_random_effect_model[["b"]],
                                    SpawningGravel_random_effect_model[["b"]]),
                             ci.lb = c(Barrier_Removal_random_effect_model[["ci.lb"]],
                                       Instream_random_effect_model[["ci.lb"]],
                                       Liming_random_effect_model[["ci.lb"]],
                                       SpawningGravel_random_effect_model[["ci.lb"]]),
                             ci.ub = c(Barrier_Removal_random_effect_model[["ci.ub"]],
                                       Instream_random_effect_model[["ci.ub"]],
                                       Liming_random_effect_model[["ci.ub"]],
                                       SpawningGravel_random_effect_model[["ci.ub"]]),
                             p = c(Barrier_Removal_random_effect_model[["pval"]],
                                   Instream_random_effect_model[["pval"]],
                                   Liming_random_effect_model[["pval"]],
                                   SpawningGravel_random_effect_model[["pval"]])) #Dataframe of salmonid responses to restoration type (overall effect size, number of effect sizes, confidence intervals, p-value)
levels(MA_restoration$restoration) <- 
  gsub(" ", "\n", levels(MA_restoration$restoration)) # Add line breaks 
MA_restoration$k <- as.character(MA_restoration$k) #Set "K" as a charcter
MA_restoration$restoration <- paste(MA_restoration$restoration, parenthesize(MA_restoration$k)) #Add "K" in parantheses to restoration description
MA_restoration$restoration <- as.character(MA_restoration$restoration)
MA_restoration$restoration <- factor(MA_restoration$restoration, 
                                     levels=
                                       c("Spawning\nGravel (17)", 
                                         "In-stream\nStructures (38)", 
                                         "Channel\nRestoration (6)",
                                         "Barrier\nRemoval (17)",
                                         "Nutrient\nRestoration (2)",
                                         "Liming (13)")) #Re-order dataframe for plotting

restoration_sub_group <- 
  ggplot(MA_restoration, aes(x = restoration, y = yi, ymin = ci.lb, ymax = ci.ub)) +  
  geom_pointrange() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  scale_y_continuous(breaks = c(-5,0,5,10,15), 
                     limit = c(-5,15), expand = c(0,0)) +
  xlab("Restoration Type") +
  ylab("Effect Size (95% CI)") +
  theme(axis.text.x = element_text(size = 13, color = "black"), 
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 15, vjust = -0.75),
        axis.title.y = element_text(size = 15),
        axis.ticks.x = element_blank()) #Point range plot of salmonid responses to restoration type
restoration_sub_group


###Moderator analyses (factors influencing restoration)
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis effect size dataset
MA_data <- subset(MA_data, after_sample_size!="1") #Remove projects with one sample size, cannot compute standardised mean difference effect sizes with sample size of one
MA_data$publication_year <- as.numeric(as.character(MA_data$publication_year)) #Set publication yeaar variable to numeric
pairs.panels(MA_data[,c(16:18,28)]) #Test for collinearity between moderators 
effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       data = MA_data) #Calculate standardised mean difference effect sizes 
random_effect_model_latitude <- rma.mv(yi, 
                                       vi, 
                                       mods = ~ Latitude, 
                                       random = ~ 1| study_ID,
                                       data = effect_sizes) #Multi-level random-effects model with REML with latitude as a moderator variable
summary.rma(random_effect_model_latitude) #Print model results
res <- rma.uni(yi, 
               vi, 
               mods = ~ Latitude,
               data = effect_sizes) #Multi-level random-effects model for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality 
random_effect_model_riverdepth <- rma.mv(yi, 
                                         vi, 
                                         mods = ~ river_depth, 
                                         random = ~ 1| study_ID,
                                         data = effect_sizes) #Multi-level random-effects model with REML with river depth as a moderator variable
summary.rma(random_effect_model_riverdepth) #Print model results
res <- rma.uni(yi, 
               vi, 
               mods = ~ river_depth,
               data = effect_sizes) #Multi-level random-effects model for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality
random_effect_model_riverwidth <- rma.mv(yi, 
                                         vi, 
                                         mods = ~ river_width, 
                                         random = ~ 1| study_ID,
                                         data = effect_sizes) #Multi-level random-effects model with REML with river width as a moderator variable
summary.rma(random_effect_model_riverwidth) #print model results
res <- rma.uni(yi, 
               vi,
               mods = ~ river_width,
               data = effect_sizes) #Multi-level random-effects model for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality
random_effect_model_landuse <- rma.mv(yi, 
                                           vi, 
                                           mods = ~ factor(land_use), 
                                           random = ~ 1| study_ID,
                                           data = effect_sizes) #Multi-level random-effects model with REML with reach land-use as a moderator variable
summary.rma(random_effect_model_landuse) #Print model results
res <- rma.uni(yi, 
               vi,
               mods = ~ land_use,
               data = effect_sizes) #Multi-level random-effects model for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality
random_effect_model_stocking <- rma.mv(yi, # effect size from each row in database
                                       vi, # measure of variance from each row in database
                                       mods = ~ factor(stocking), 
                                       random = ~ 1| study_ID,
                                       data = effect_sizes) #Multi-level meta-analysis random-effects model with REML with presence of stocking as a moderator variable
summary.rma(random_effect_model_stocking) #Print model results
res <- rma.uni(yi, # effect size from each row in database
               vi, # measure of variance from each row in database
               mods = ~ stocking,
               data = effect_sizes) #Multi-level random-effects model for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality
random_effect_model_studydesign <- rma.mv(yi, 
                                          vi, 
                                          mods = ~ factor(study_design), 
                                          random = ~ 1 | study_ID, #study ID as random effect
                                          data = effect_sizes) #Multi-level meta-analysis random-effects model with REML with study design as a moderator variable
summary.rma(random_effect_model_studydesign) #print model results
res <- rma.uni(yi, 
               vi, 
               mods = ~ study_design,
               data = effect_sizes) #Multi-level random-effects model for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality
random_effect_model_projectage <- rma.mv(yi, 
                                         vi, 
                                         mods = ~ post_monitoring_years, 
                                         random = ~ 1| study_ID,
                                         data = effect_sizes) #Multi-level random-effects model with REML with project age as a moderator variable
summary.rma(random_effect_model_projectage) #print model results
res <- rma.uni(yi, 
               vi, 
               mods = ~ post_monitoring_years,
               data = effect_sizes) #Multi-level random-effects model for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality
project_age_plot <- 
  ggplot(effect_sizes, aes(x = post_monitoring_years, y = yi)) + 
  geom_point() +
  geom_smooth(method = lm) + 
  theme_classic() +
  xlab("Post-restoration monitoring length (years)") +
  ylab("Effect Size (Standardised Mean Difference)") +
  theme(axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) #Scatterplot with effect size as a function of project age
project_age_plot

random_effect_model_restoration <- rma.mv(yi, 
                                          vi, 
                                          mods = ~ factor(restoration),
                                          random = ~ 1 | study_ID, #study ID as random effect
                                          data = effect_sizes) #Random-effects model with REML with project age as a moderator variable
summary.rma(random_effect_model_restoration) #Print model results
res <- rma.uni(yi, 
               vi, 
               mods = ~ factor(restoration),
               data = effect_sizes) #Multi-level random-effects model for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality

restoration_mods <- subset(effect_sizes, restoration!="Channel") #Exclude restoration types with fewer than 10 effect sizes
restoration_mods <- subset(restoration_mods, restoration!="Nutrient") #Exclude restoration types with fewer than 10 effect sizes
restoration_mods <- subset(restoration_mods, restoration!="Riparian") #Exclude restoration types with fewer than 10 effect sizes
random_effect_model_restoration_sub <- rma.mv(yi, 
                                          vi, # 
                                          mods = ~ factor(restoration), 
                                          random = ~ 1 | study_ID, #study ID as random effect
                                          data = restoration_mods) #Multi-level random-effects model with REML with restoration as a moderator variable
summary.rma(random_effect_model_restoration_sub) #Print model
summary(glht(random_effect_model_restoration, 
             linfct=cbind(contrMat(rep(1,7), type="Tukey"))), test=adjusted("none")) #Post-hoc Tukey test for pairwise comparisons
random_effect_model_restoration <- rma.mv(yi, 
                                          vi, 
                                          mods = ~ factor(restoration)*species, 
                                          random = ~ 1 | study_ID, #study ID as random effect
                                          data = restoration_mods) #Multi-level random-effects model with REML with restoration as a moderator variable
summary.rma(random_effect_model_restoration_sub) #Print model
anova(random_effect_model_restoration) 
summary(glht(random_effect_model_restoration, 
             linfct=cbind(contrMat(rep(1,4), type="Tukey"))), test=adjusted("none")) #Post-hoc Tukey test for pairwise comparisons
res <- rma.uni(yi, # effect size from each row in database
               vi, # measure of variance from each row in database
               mods = ~ factor(restoration),
               data = restoration_mods) #Multi-level random-effects model for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality



#########Moderator sensitivity analyses
####Latitude
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis effect size dataset
MA_data <- subset(MA_data, after_sample_size!="1") #Remove projects with one sample size, cannot compute standardised mean difference effect sizes with sample size of one
MA_data$publication_year <- as.numeric(as.character(MA_data$publication_year)) #Set publication yeaar variable to numeric
effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       data = MA_data) #Calculate standardised mean difference effect sizes 
effect_sizes <- subset(effect_sizes, Latitude!="NA") #Subset latitude data
random_effect_model_latitude <- rma.mv(yi, 
                                       vi, 
                                       mods = ~ Latitude, 
                                       random = ~ 1| study_ID,
                                       data = effect_sizes) #Multi-level random-effects model with REML with latitude as a moderator variable
#Standardized residuals
rs <- rstudent.rma.mv(random_effect_model_latitude)
rs <- as.data.frame(rs) #Residuals in dataframe
limit_rs <- 3 #Set limit of standardised deleted residuals to 3
rs$ES_ID <- as.numeric(1:nrow(rs)) #Add a new column with the effect size ID number

#Hat values
hat <- hatvalues.rma.mv(random_effect_model_latitude)
hat <- as.data.frame(hat)
mean <- mean(hat$hat) #Average hat value
limit_hat <- 2 * mean #Set limit of hat values to twice the average hat value
hat$ES_ID <- as.numeric(1:nrow(hat)) #Add a new column with the effect size ID number

effect_sizes$ES_ID <- as.numeric(1:nrow(effect_sizes)) #Add a new column with the effect size ID number
sensitivity_res_latitude <- left_join(hat,rs, by="ES_ID")%>%
  left_join(effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat<limit_hat)

sensitivity_latitude <- rma.mv(yi, 
                               vi,
                               mods = ~ Latitude,
                               random = ~ 1| study_ID, #Study ID as random effect
                               data = sensitivity_res_latitude) #Multi-level random-effects model with REML with latitude as a moderator variable
summary.rma(sensitivity_latitude) #Print model results without outliers
summary(random_effect_model_latitude) #Print model with outliers


####River depth
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis effect size dataset
MA_data <- subset(MA_data, after_sample_size!="1") #Remove projects with one sample size, cannot compute standardised mean difference effect sizes with sample size of one
MA_data$publication_year <- as.numeric(as.character(MA_data$publication_year)) #Set publication yeaar variable to numeric
effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       data = MA_data) #Calculate standardised mean difference effect sizes 
effect_sizes <- subset(effect_sizes, river_depth!="NA")
random_effect_model_riverdepth <- rma.mv(yi, # effect size from each row in database
                                         vi, # measure of variance from each row in database
                                         mods = ~ river_depth, 
                                         random = ~ 1| study_ID,
                                         data = effect_sizes) #Multivariate meta-analysis random-effects model with REML with river depth as a moderator variable
#Standardized residuals
rs <- rstudent.rma.mv(random_effect_model_riverdepth)
rs <- as.data.frame(rs) #Residuals in dataframe
limit_rs <- 3 #Set limit of standardised deleted residuals to 3
rs$ES_ID <- as.numeric(1:nrow(rs)) #Add a new column with the effect size ID number

#Hat values
hat <- hatvalues.rma.mv(random_effect_model_riverdepth)
hat <- as.data.frame(hat)
mean <- mean(hat$hat)
limit_hat <- 2 * mean
hat$ES_ID <- as.numeric(1:nrow(hat)) #add a new column with the effect size ID number

effect_sizes$ES_ID <- as.numeric(1:nrow(effect_sizes))
sensitivity_res_depth <- left_join(hat,rs, by="ES_ID")%>%
  left_join(effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat<limit_hat)

sensitivity_depth <- rma.mv(yi, #Effect size from each row in database
                            vi, #Measure of variance from each row in database
                            mods = ~ river_depth,
                            random = ~ 1| study_ID, #Study ID as random effect
                            data = sensitivity_res_depth) #Multivariate meta-analysis random-effects model with REML
summary.rma(sensitivity_depth) #Print model results without outliers
summary(random_effect_model_riverdepth) #Print model with outliers


####River width
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis effect size dataset
MA_data <- subset(MA_data, after_sample_size!="1") #Remove projects with one sample size, cannot compute standardised mean difference effect sizes with sample size of one
MA_data$publication_year <- as.numeric(as.character(MA_data$publication_year)) #Set publication yeaar variable to numeric
effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       data = MA_data) #Calculate standardised mean difference effect sizes 
effect_sizes <- subset(effect_sizes, river_width!="NA")
random_effect_model_riverwidth <- rma.mv(yi, # effect size from each row in database
                                         vi, # measure of variance from each row in database
                                         mods = ~ river_width, 
                                         random = ~ 1| study_ID,
                                         data = effect_sizes) #Multivariate meta-analysis random-effects model with REML with river width as a moderator variable
#Standardized residuals
rs <- rstudent.rma.mv(random_effect_model_riverwidth)
rs <- as.data.frame(rs) #Residuals in dataframe
limit_rs <- 3 #Set limit of standardised deleted residuals to 3
rs$ES_ID <- as.numeric(1:nrow(rs)) #Add a new column with the effect size ID number

#Hat values
hat <- hatvalues.rma.mv(random_effect_model_riverwidth)
hat <- as.data.frame(hat)
mean <- mean(hat$hat)
limit_hat <- 2 * mean
hat$ES_ID <- as.numeric(1:nrow(hat)) #add a new column with the effect size ID number

effect_sizes$ES_ID <- as.numeric(1:nrow(effect_sizes))
sensitivity_res_width <- left_join(hat,rs, by="ES_ID")%>%
  left_join(effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat<limit_hat)

sensitivity_width <- rma.mv(yi, #Effect size from each row in database
                            vi, #Measure of variance from each row in database
                            mods = ~ river_width,
                            random = ~ 1| study_ID, #Study ID as random effect
                            data = sensitivity_res_width) #Multivariate meta-analysis random-effects model with REML
summary.rma(sensitivity_width) #Print model results without outliers
summary(random_effect_model_riverwidth) #Print model with outliers


####Land use
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis effect size dataset
MA_data <- subset(MA_data, after_sample_size!="1") #Remove projects with one sample size, cannot compute standardised mean difference effect sizes with sample size of one
MA_data$publication_year <- as.numeric(as.character(MA_data$publication_year)) #Set publication yeaar variable to numeric
effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       data = MA_data) #Calculate standardised mean difference effect sizes
effect_sizes <- subset(effect_sizes, land_use!="NA")
random_effect_model_landuse <- rma.mv(yi, # effect size from each row in database
                                      vi, # measure of variance from each row in database
                                      mods = ~ factor(land_use), 
                                      random = ~ 1| study_ID,
                                      data = effect_sizes) #Multivariate meta-analysis random-effects model with REML with reach land-use as a moderator variable

#Standardized residuals
rs <- rstudent.rma.mv(random_effect_model_landuse)
rs <- as.data.frame(rs) #Residuals in dataframe
limit_rs <- 3 #Set limit of standardised deleted residuals to 3
rs$ES_ID <- as.numeric(1:nrow(rs)) #Add a new column with the effect size ID number

#Hat values
hat <- hatvalues.rma.mv(random_effect_model_landuse)
hat <- as.data.frame(hat)
mean <- mean(hat$hat)
limit_hat <- 2 * mean
hat$ES_ID <- as.numeric(1:nrow(hat)) #add a new column with the effect size ID number

effect_sizes$ES_ID <- as.numeric(1:nrow(effect_sizes))
sensitivity_res_landuse <- left_join(hat,rs, by="ES_ID")%>%
  left_join(effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat<limit_hat)

sensitivity_landuse <- rma.mv(yi, #Effect size from each row in database
                              vi, #Measure of variance from each row in database
                              mods = ~ factor(land_use),
                              random = ~ 1| study_ID, #Study ID as random effect
                              data = sensitivity_res_landuse) #Multivariate meta-analysis random-effects model with REML
summary.rma(sensitivity_landuse) #Print model results without outliers
summary(random_effect_model_landuse) #Print model with outliers



####Stocking
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis effect size dataset
MA_data <- subset(MA_data, after_sample_size!="1") #Remove projects with one sample size, cannot compute standardised mean difference effect sizes with sample size of one
MA_data$publication_year <- as.numeric(as.character(MA_data$publication_year)) #Set publication yeaar variable to numeric
effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       data = MA_data) #Calculate standardised mean difference effect sizes 
effect_sizes <- subset(effect_sizes, stocking!="U") #Compare stocked and unstocked rivers, not unknown rivers
random_effect_model_stocking <- rma.mv(yi, # effect size from each row in database
                                       vi, # measure of variance from each row in database
                                       mods = ~ factor(stocking), 
                                       random = ~ 1| study_ID,
                                       data = effect_sizes) #Multivariate meta-analysis random-effects model with REML with presence of stocking as a moderator variable
#Standardized residuals
rs <- rstudent.rma.mv(random_effect_model_stocking)
rs <- as.data.frame(rs) #Residuals in dataframe
limit_rs <- 3 #Set limit of standardised deleted residuals to 3
rs$ES_ID <- as.numeric(1:nrow(rs)) #Add a new column with the effect size ID number

#Hat values
hat <- hatvalues.rma.mv(random_effect_model_stocking)
hat <- as.data.frame(hat)
mean <- mean(hat$hat)
limit_hat <- 2 * mean
hat$ES_ID <- as.numeric(1:nrow(hat)) #add a new column with the effect size ID number

effect_sizes$ES_ID <- as.numeric(1:nrow(effect_sizes))
sensitivity_res_stocking <- left_join(hat,rs, by="ES_ID")%>%
  left_join(effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat<limit_hat)

sensitivity_stocking <- rma.mv(yi, #Effect size from each row in database
                               vi, #Measure of variance from each row in database
                               mods = ~ factor(stocking),
                               random = ~ 1| study_ID, #Study ID as random effect
                               data = sensitivity_res_stocking) #Multivariate meta-analysis random-effects model with REML
summary.rma(sensitivity_stocking) #Print model results without outliers
summary(random_effect_model_stocking) #Print model with outliers


####Project age
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis effect size dataset
MA_data <- subset(MA_data, after_sample_size!="1") #Remove projects with one sample size, cannot compute standardised mean difference effect sizes with sample size of one
MA_data$publication_year <- as.numeric(as.character(MA_data$publication_year)) #Set publication yeaar variable to numeric
effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       data = MA_data) #Calculate standardised mean difference effect sizes 
effect_sizes <- subset(effect_sizes, post_monitoring_years!="NA") 
random_effect_model_projectage <- rma.mv(yi, # effect size from each row in database
                                         vi, # measure of variance from each row in database
                                         mods = ~ post_monitoring_years, 
                                         random = ~ 1| study_ID,
                                         data = effect_sizes) #Multivariate meta-analysis random-effects model with REML with project age as a moderator variable
#Standardized residuals
rs <- rstudent.rma.mv(random_effect_model_projectage)
rs <- as.data.frame(rs) #Residuals in dataframe
limit_rs <- 3 #Set limit of standardised deleted residuals to 3
rs$ES_ID <- as.numeric(1:nrow(rs)) #Add a new column with the effect size ID number

#Hat values
hat <- hatvalues.rma.mv(random_effect_model_projectage)
hat <- as.data.frame(hat)
mean <- mean(hat$hat)
limit_hat <- 2 * mean
hat$ES_ID <- as.numeric(1:nrow(hat)) #add a new column with the effect size ID number

effect_sizes$ES_ID <- as.numeric(1:nrow(effect_sizes))
sensitivity_res_projectage <- left_join(hat,rs, by="ES_ID")%>%
  left_join(effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat<limit_hat)

sensitivity_projectage <- rma.mv(yi, #Effect size from each row in database
                                 vi, #Measure of variance from each row in database
                                 mods = ~ post_monitoring_years,
                                 random = ~ 1| study_ID, #Study ID as random effect
                                 data = sensitivity_res_projectage) #Multivariate meta-analysis random-effects model with REML
summary.rma(sensitivity_projectage) #Print model results without outliers
summary(random_effect_model_projectage) #Print model with outliers

####Restoration
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis effect size dataset
MA_data <- subset(MA_data, after_sample_size!="1") #Remove projects with one sample size, cannot compute standardised mean difference effect sizes with sample size of one
MA_data$publication_year <- as.numeric(as.character(MA_data$publication_year)) #Set publication yeaar variable to numeric
effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       data = MA_data) #Calculate standardised mean difference effect sizes
#Standardized residuals
rs <- rstudent.rma.mv(random_effect_model_restoration)
rs <- as.data.frame(rs) #Residuals in dataframe
limit_rs <- 3 #Set limit of standardised deleted residuals to 3
rs$ES_ID <- as.numeric(1:nrow(rs)) #Add a new column with the effect size ID number

#Hat values
hat <- hatvalues.rma.mv(random_effect_model_restoration)
hat<- as.data.frame(hat)
mean <- mean(hat$hat)
limit_hat <- 2 * mean
hist(hat$hat)
hat$ES_ID <- as.numeric(1:nrow(hat)) #Add a new column with the effect size ID number

#Plot hat values agains residual values
plot(x=hat$hat, y= rs$resid,
     ylab= "Standardized residuals", xlab= "Hat values")
text(x=hat$hat, y= rs$resid, labels = effect_sizes$ES_ID, cex= 1, pos = 2)
abline(v = limit_hat)
abline(h = -limit_rs)
abline(h = limit_rs)

restoration_mods$ES_ID <- as.numeric(1:nrow(restoration_mods))
sensitivity_res_restoration <- left_join(hat,rs, by="ES_ID")%>%
  left_join(restoration_mods, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat<limit_hat)

sensitivity_restoration <- rma.mv(yi, #Effect size from each row in database
                                  vi, #Measure of variance from each row in database
                                  mods = ~ factor(restoration),
                                  random = ~ 1| study_ID, #Study ID as random effect
                                  data = sensitivity_res_restoration) #Multivariate meta-analysis random-effects model with REML
summary.rma(sensitivity_restoration) #Print model results without outliers
summary.rma(random_effect_model_restoration) #Print model results with outliers
summary(glht(sensitivity_restoration, 
             linfct=cbind(contrMat(rep(1,4), type="Tukey"))), test=adjusted("bonferroni")) #Post hoc Tukey test for pairwise comparisons




###Stocking random-effects models
N_stock <- subset(MA_data, MA_data$stocking=="N") #Subset without stocking data
N_stock_effect_sizes <- escalc("SMD", 
                               m2i = before_mean,  
                               sd2i = before_SD,  
                               n2i = before_sample_size, 
                               m1i = after_mean, 
                               sd1i = after_SD, 
                               n1i = after_sample_size, 
                               data = N_stock) #Calculate standardised mean difference effect sizes 
N_stock_random_effect_model <- rma.mv(yi, 
                                    vi, 
                                    method = "REML", 
                                    random = ~ 1| study_ID,
                                    data = N_stock_effect_sizes) #Multi-level random-effects model with REML estimator for rivers without stocking
summary.rma(N_stock_random_effect_model) #Print model results

Y_stock <- subset(MA_data, MA_data$stocking=="Y") #Subset with stocking data
Y_stock_effect_sizes <- escalc("SMD", 
                               m2i = before_mean,  
                               sd2i = before_SD,  
                               n2i = before_sample_size, 
                               m1i = after_mean, 
                               sd1i = after_SD, 
                               n1i = after_sample_size, 
                               data = Y_stock) #Calculate standardised mean difference effect sizes
Y_stock_random_effect_model <- rma.mv(yi, 
                                      vi, 
                                      method = "REML", 
                                      random = ~ 1| study_ID,
                                      data = Y_stock_effect_sizes) #Multi-level random-effects model with REML estimator for rivers with stocking
summary.rma(Y_stock_random_effect_model) #Print model results

###Land-use subset analysis
Agri <- subset(MA_data, MA_data$land_use=="Agriculture") #Subset rivers in agricultural areas
Agri_effect_sizes <- escalc("SMD", 
                               m2i = before_mean,  
                               sd2i = before_SD,  
                               n2i = before_sample_size, 
                               m1i = after_mean, 
                               sd1i = after_SD, 
                               n1i = after_sample_size, 
                               data = Agri) #Calculate standardised mean difference effect sizes
Agri_random_effect_model <- rma.mv(yi, 
                                      vi,
                                      method = "REML", 
                                      random = ~ 1| study_ID,
                                      data = Agri_effect_sizes) #Multi-level random-effects model with REML estimator for rivers in agricultural areas
summary.rma(Agri_random_effect_model) #Print model results

Natural <- subset(MA_data, MA_data$land_use=="Natural") #Subset rivers in natural areas (i.e., forest, meadow, mountains)
Natural_effect_sizes <- escalc("SMD", 
                            m2i = before_mean,  
                            sd2i = before_SD,  
                            n2i = before_sample_size, 
                            m1i = after_mean, 
                            sd1i = after_SD, 
                            n1i = after_sample_size, 
                            data = Natural) #Calculate standardised mean difference effect sizes 
Natural_random_effect_model <- rma.mv(yi, 
                                   vi, 
                                   method = "REML", 
                                   random = ~ 1| study_ID,
                                   data = Natural_effect_sizes) #Multi-level random-effects model with REML estimator for rivers in natural areas
summary.rma(Natural_random_effect_model) #Print model results
