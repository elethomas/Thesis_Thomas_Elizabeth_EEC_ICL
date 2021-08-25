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
sum(data$monitoring>=5, na.rm = TRUE) #Number of projects monitoring for less than 5 years
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

###Calculate effect sizes
effect_sizes <- escalc("SMD", 
                       m2i = before_mean,  
                       sd2i = before_SD,  
                       n2i = before_sample_size, 
                       m1i = after_mean, 
                       sd1i = after_SD, 
                       n1i = after_sample_size, 
                       vtype = "AV",
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
res <- rma.uni(yi, 
               vi, 
               data = sensitivity_res) #Heterogeneity of sensitivity model
res


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

###Sensitivity analysis of species as a moderator
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
       xlab = "Effect size (Hedge's g)", # Label for x-axis
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
       slab = trout_effect_sizes$author, #Label individual studies
       mlab = "RE Model",
       annotate = TRUE, #Display individual effect sizes and 95% confidence intervals
       xlab = "Effect size (Hedge's g)", # Label for x-axis
       cex = 1.25, # Text side for study labels
       pch = 16, # Shape of effect size in forest plot
       psize = 0.75,
       cex.lab = 1.75, # Size of x-axis label
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
           fill = "grey65") +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.08) +
  theme_classic() +
  xlab("Species") +
  ylab("Effect Size (95% CI)") +
  theme(axis.text.x = element_text(size = 16, color = "black"), 
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0,0),
                     limit = c(0,1.5)) #Plot salmon and trout overall effect sizes

####Sensitvity analysis trout and salmon
##Salmon
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
                              data = salmon) 
random_effect_model_SMD_Salmon <- rma.mv(yi, # effect size from each row in database
                                         vi, # measure of variance from each row in database
                                         random = ~ 1 | study_ID,
                                         method = "REML",
                                         data = salmon_effect_sizes) #Random-effects model with REML estimator 
summary.rma(random_effect_model_SMD_Salmon) #Print model results
#Standardized residuals
rs_salmon <- rstudent.rma.mv(random_effect_model_SMD_Salmon)
rs_salmon <- as.data.frame(rs_salmon)
limit_rs <- 3
rs_salmon$ES_ID <- as.numeric(1:nrow(rs_salmon)) #add a new column with the effect size ID number

#Hat values
hat_salmon <- hatvalues.rma.mv(random_effect_model_SMD_Salmon)
hat_salmon <- as.data.frame(hat_salmon)
mean_salmon <- mean(hat_salmon$hat_salmon)
limit_hat_salmon <- 2 * mean_salmon
hat_salmon$ES_ID <- as.numeric(1:nrow(hat_salmon)) #add a new column with the effect size ID number

salmon_effect_sizes$ES_ID <- as.numeric(1:nrow(salmon_effect_sizes))
sensitivity_res_salmon <- left_join(hat_salmon,rs_salmon, by="ES_ID")%>%
  left_join(salmon_effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat_salmon<limit_hat_salmon)
sensitivity_salmon <- rma.mv(yi, #Effect size from each row in database
                             vi, #Measure of variance from each row in database
                             random = ~ 1| study_ID, #Study ID as random effect
                             data = sensitivity_res_salmon) #Multivariate meta-analysis random-effects model with REML
summary.rma(sensitivity_salmon) #Print model results without outliers
summary(random_effect_model_SMD_Salmon) #Print model with outliers
res_salmon <- rma.uni(yi, 
                      vi, 
                      data = sensitivity_res_salmon) #Heterogeneity of sensitivity model
res_salmon


####Trout
trout <- subset(MA_data, species=="Trout") #Subset trout data
trout_effect_sizes <- escalc("SMD", 
                             m2i = before_mean,  
                             sd2i = before_SD,  
                             n2i = before_sample_size, 
                             m1i = after_mean, 
                             sd1i = after_SD, 
                             n1i = after_sample_size, 
                             data = trout) #Calculate standardised mean difference effect sizes for each row in dataset
random_effect_model_SMD_Trout <- rma.mv(yi, # effect size from each row in database
                                        vi, # measure of variance from each row in database
                                        random = ~ 1 | study_ID,
                                        method = "REML",
                                        data = trout_effect_sizes) #Random-effects model with REML estimator 
#Standardized residuals
rs_trout <- rstudent.rma.mv(random_effect_model_SMD_Trout)
rs_trout <- as.data.frame(rs_trout)
limit_rs <- 3
rs_trout$ES_ID <- as.numeric(1:nrow(rs_trout)) #add a new column with the effect size ID number

#Hat values
hat_trout <- hatvalues.rma.mv(random_effect_model_SMD_Trout)
hat_trout<- as.data.frame(hat_trout)
mean_trout <- mean(hat_trout$hat_trout)
limit_hat_trout <- 2 * mean_trout
hat_trout$ES_ID <- as.numeric(1:nrow(hat_trout)) #add a new column with the effect size ID number

trout_effect_sizes$ES_ID <- as.numeric(1:nrow(trout_effect_sizes))
sensitivity_res_trout <- left_join(hat_trout,rs_trout, by="ES_ID")%>%
  left_join(trout_effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat_trout<limit_hat_trout)
sensitivity_trout <- rma.mv(yi, #Effect size from each row in database
                            vi, #Measure of variance from each row in database
                            random = ~ 1| study_ID, #Study ID as random effect
                            data = sensitivity_res_trout) #Multivariate meta-analysis random-effects model with REML
summary.rma(sensitivity_trout) #Print model results without outliers
summary(random_effect_model_SMD_Trout) #Print model with outliers
res_trout <- rma.uni(yi, 
                     vi, 
                     data = sensitivity_res_trout) #Heterogeneity of sensitivity model
res_trout

###Salmon and trout specific responses to restoration type
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis data 
MA_data <- MA_data %>% filter(before_sample_size!=1) #Remove projects with sample size of 1
salmon <- subset(MA_data, species=="Salmon") #Subset salmon data
Salmon_SMD_effect_sizes <- escalc("SMD", 
                                  m2i = before_mean, 
                                  sd2i = before_SD,  
                                  n2i = before_sample_size, 
                                  m1i = after_mean, 
                                  sd1i = after_SD, 
                                  n1i = after_sample_size, 
                                  data = salmon) #Calculate standardised mean difference effect size 
salmon %>% group_by(restoration) %>% summarise(n()) #Number of projects per restoration category in salmon dataset
#Minimum sample size of 10 projects for restoration sub-group analysis

Barrier_Removal_Salmon <- rma.mv(yi, 
                                 vi, 
                                 method = "REML", 
                                 subset = (restoration=="Barrier_removal"),
                                 random = ~ 1| study_ID,
                                 slab = paste(author, publication_year, sep = " "),
                                 data = Salmon_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(Barrier_Removal_Salmon) #Print model results
res_barrier_salmon <- rma.uni(yi, vi, subset = (restoration=="Barrier_removal"), 
                              data = Salmon_SMD_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_barrier_salmon) #Print model results
qqnorm.rma.uni(res_barrier_salmon, envelope = FALSE) #Model validation - normality of residuals

Channel_Salmon <- rma.mv(yi, # effect size from each row in database
                         vi, # measure of variance from each row in database
                         method = "REML", # specifies REML estimator for random effects model
                         subset = (restoration=="Channel"),
                         random = ~ 1| study_ID,
                         slab = paste(author, publication_year, sep = " "),
                         data = Salmon_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(Channel_Salmon) #Print model results
res_channel_salmon <- rma.uni(yi, vi, subset = (restoration=="Channel"), 
                              data = Salmon_SMD_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_channel_salmon) #Print model results
qqnorm.rma.uni(res_channel_salmon, envelope = FALSE) #Model validation - normality of residuals

Instream_Salmon <- rma.mv(yi, # effect size from each row in database
                          vi, # measure of variance from each row in database
                          method = "REML", # specifies REML estimator for random effects model
                          subset = (restoration=="In-stream_structure"),
                          random = ~ 1| study_ID,
                          slab = paste(author, publication_year, sep = " "),
                          data = Salmon_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(Instream_Salmon) #Print model results
res_instream_salmon <- rma.uni(yi, vi, subset = (restoration=="In-stream_structure"), 
                               data = Salmon_SMD_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_instream_salmon) #Print model results
qqnorm.rma.uni(res_instream_salmon, envelope = FALSE) #Model validation - normality of residuals

Liming_Salmon <- rma.mv(yi, # effect size from each row in database
                        vi, # measure of variance from each row in database
                        method = "REML", # specifies REML estimator for random effects model
                        subset = (restoration=="Liming"),
                        random = ~ 1| study_ID,
                        slab = paste(author, publication_year, sep = " "),
                        data = Salmon_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(Liming_Salmon) #Print model results
res_liming_salmon <- rma.uni(yi, vi, subset = (restoration=="Liming"), 
                             data = Salmon_SMD_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_liming_salmon) #Print model results
qqnorm.rma.uni(res_liming_salmon, envelope = FALSE) #Model validation - normality of residuals

Nutrient_Salmon <- rma.mv(yi, # effect size from each row in database
                          vi, # measure of variance from each row in database
                          method = "REML", # specifies REML estimator for random effects model
                          subset = (restoration=="Nutrient"),
                          random = ~ 1| study_ID,
                          slab = paste(author, publication_year, sep = " "),
                          data = Salmon_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(Nutrient_Salmon) #Print model results
res_nutrient_salmon <- rma.uni(yi, vi, subset = (restoration=="Nutrient"), 
                               data = Salmon_SMD_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_nutrient_salmon) #Print model results
qqnorm.rma.uni(res_nutrient_salmon, envelope = FALSE) #Model validation - normality of residuals

Riparian_Salmon <- rma.mv(yi, # effect size from each row in database
                          vi, # measure of variance from each row in database
                          method = "REML", # specifies REML estimator for random effects model
                          subset = (restoration=="Riparian"),
                          random = ~ 1| study_ID,
                          slab = paste(author, publication_year, sep = " "),
                          data = Salmon_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(Riparian_Salmon) #Print model results, too few studies

SpawningGravel_Salmon <- rma.mv(yi, # effect size from each row in database
                                vi, # measure of variance from each row in database
                                method = "REML", # specifies REML estimator for random effects model
                                subset = (restoration=="Spawning_gravel"),
                                random = ~ 1| study_ID,
                                slab = paste(author, publication_year, sep = " "),
                                data = Salmon_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(SpawningGravel_Salmon) #Print model results
res_spawning_salmon <- rma.uni(yi, vi, subset = (restoration=="Spawning_gravel"), 
                               data = Salmon_SMD_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_spawning_salmon) #Print model results
qqnorm.rma.uni(res_spawning_salmon, envelope = FALSE) #Model validation - normality of residuals

Salmon_restoration <- data.frame(species = "Salmon",
                                 restoration = c("Barrier removal", "Channel restoration", 
                                                 "In-stream structures",
                                                 "Liming", "Nutrient restoration", 
                                                 "Riparian restoration", "Spawning gravel"),
                                 k = c(Barrier_Removal_Salmon[["k"]], 
                                       Channel_Salmon[["k"]],
                                       Instream_Salmon[["k"]],
                                       Liming_Salmon[["k"]],
                                       Nutrient_Salmon[["k"]],
                                       0,
                                       SpawningGravel_Salmon[["k"]]),
                                 yi = c(Barrier_Removal_Salmon[["b"]], 
                                        Channel_Salmon[["b"]],
                                        Instream_Salmon[["b"]],
                                        Liming_Salmon[["b"]],
                                        Nutrient_Salmon[["b"]],
                                        "NA",
                                        SpawningGravel_Salmon[["b"]]),
                                 ci.lb = c(Barrier_Removal_Salmon[["ci.lb"]], 
                                           Channel_Salmon[["ci.lb"]],
                                           Instream_Salmon[["ci.lb"]],
                                           Liming_Salmon[["ci.lb"]],
                                           Nutrient_Salmon[["ci.lb"]],
                                           "NA",
                                           SpawningGravel_Salmon[["ci.lb"]]),
                                 ci.ub = c(Barrier_Removal_Salmon[["ci.ub"]], 
                                           Channel_Salmon[["ci.ub"]],
                                           Instream_Salmon[["ci.ub"]],
                                           Liming_Salmon[["ci.ub"]],
                                           Nutrient_Salmon[["ci.ub"]],
                                           "NA",
                                           SpawningGravel_Salmon[["ci.ub"]]),
                                 p = c(Barrier_Removal_Salmon[["pval"]], 
                                       Channel_Salmon[["pval"]],
                                       Instream_Salmon[["pval"]],
                                       Liming_Salmon[["pval"]],
                                       Nutrient_Salmon[["pval"]],
                                       "NA",
                                       SpawningGravel_Salmon[["pval"]])) #Dataframe of salmon responses to restoration type

trout <- subset(MA_data, species=="Trout") #Subset trout data
Trout_SMD_effect_sizes <- escalc("SMD", 
                                 m2i = before_mean, 
                                 sd2i = before_SD,  
                                 n2i = before_sample_size, 
                                 m1i = after_mean, 
                                 sd1i = after_SD, 
                                 n1i = after_sample_size, 
                                 data = trout) #Calculate standardised mean difference effect size 

trout %>% group_by(restoration) %>% summarise(n()) #Number of projects per restoration category in trout dataset
#Minimum sample size of 10 projects for restoration sub-group analysis

Barrier_Removal_Trout <- rma.mv(yi, # effect size from each row in database
                                vi, # measure of variance from each row in database
                                method = "REML", # specifies REML estimator for random effects model
                                subset = (restoration=="Barrier_removal"),
                                random = ~ 1| study_ID,
                                slab = paste(author, publication_year, sep = " "),
                                data = Trout_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(Barrier_Removal_Trout) #Print model results
res_barrier_trout <- rma.uni(yi, vi, subset = (restoration=="Barrier_removal"), 
                             data = Trout_SMD_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_barrier_trout) #Print model results
qqnorm.rma.uni(res_barrier_trout, envelope = FALSE) #Model validation - normality of residuals

channel <- subset(Trout_SMD_effect_sizes, yi < 5)
Channel_Trout <- rma.mv(yi, # effect size from each row in database
                        vi, # measure of variance from each row in database
                        method = "REML", # specifies REML estimator for random effects model
                        subset = (restoration=="Channel"),
                        random = ~ 1| study_ID,
                        slab = paste(author, publication_year, sep = " "),
                        data = channel) #Multi-level random-effects model with REML
summary.rma(Channel_Trout) #Print model results
res_channel_trout <- rma.uni(yi, vi,
                             data = channel, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_channel_trout) #Print model results
qqnorm.rma.uni(res_channel_trout, envelope = FALSE) #Model validation - normality of residuals

Instream_Trout <- rma.mv(yi, # effect size from each row in database
                         vi, # measure of variance from each row in database
                         method = "REML", # specifies REML estimator for random effects model
                         subset = (restoration=="In-stream_structure"),
                         random = ~ 1| study_ID,
                         slab = paste(author, publication_year, sep = " "),
                         data = Trout_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(Instream_Trout) #Print model results
res_instream_trout <- rma.uni(yi, vi, subset = (restoration=="In-stream_structure"), 
                              data = Trout_SMD_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_instream_trout) #Print model results
qqnorm.rma.uni(res_instream_trout, envelope = FALSE) #Model validation - normality of residuals

Liming_Trout <- rma.mv(yi, # effect size from each row in database
                       vi, # measure of variance from each row in database
                       method = "REML", # specifies REML estimator for random effects model
                       subset = (restoration=="Liming"),
                       random = ~ 1| study_ID,
                       slab = paste(author, publication_year, sep = " "),
                       data = Trout_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(Liming_Trout) #Print model results
res_liming_trout <- rma.uni(yi, vi, subset = (restoration=="Liming"), 
                            data = Trout_SMD_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_liming_trout) #Print model results
qqnorm.rma.uni(res_liming_trout, envelope = FALSE) #Model validation - normality of residuals

Nutrient_Trout <- rma.mv(yi, # effect size from each row in database
                         vi, # measure of variance from each row in database
                         method = "REML", # specifies REML estimator for random effects model
                         subset = (restoration=="Nutrient"),
                         random = ~ 1| study_ID,
                         slab = paste(author, publication_year, sep = " "),
                         data = Trout_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(Nutrient_Trout) #Print model results, too few studies

Riparian_Trout <- rma.mv(yi, # effect size from each row in database
                         vi, # measure of variance from each row in database
                         method = "REML", # specifies REML estimator for random effects model
                         subset = (restoration=="Riparian"),
                         random = ~ 1| study_ID,
                         slab = paste(author, publication_year, sep = " "),
                         data = Trout_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(Riparian_Trout) #Print model results, too few studies

SpawningGravel_Trout <- rma.mv(yi, # effect size from each row in database
                               vi, # measure of variance from each row in database
                               method = "REML", # specifies REML estimator for random effects model
                               subset = (restoration=="Spawning_gravel"),
                               random = ~ 1| study_ID,
                               slab = paste(author, publication_year, sep = " "),
                               data = Trout_SMD_effect_sizes) #Multi-level random-effects model with REML
summary.rma(SpawningGravel_Trout) #Print model results
res_spawning_trout <- rma.uni(yi, vi, subset = (restoration=="Spawning_gravel"), 
                              data = Trout_SMD_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_spawning_trout) #Print model results
qqnorm.rma.uni(res_spawning_trout, envelope = FALSE) #Model validation - normality of residuals

Trout_restoration <- data.frame(species = "Trout",
                                restoration = c("Barrier removal", "Channel restoration", 
                                                "In-stream structures",
                                                "Liming", "Nutrient restoration",
                                                "Riparian restoration", "Spawning gravel"),
                                k = c(Barrier_Removal_Trout[["k"]], 
                                      Channel_Trout[["k"]],
                                      Instream_Trout[["k"]],
                                      Liming_Trout[["k"]],
                                      0,
                                      1,
                                      SpawningGravel_Trout[["k"]]),
                                yi = c(Barrier_Removal_Trout[["b"]], 
                                       Channel_Trout[["b"]],
                                       Instream_Trout[["b"]],
                                       Liming_Trout[["b"]],
                                       "NA",
                                       "NA",
                                       SpawningGravel_Trout[["b"]]),
                                ci.lb = c(Barrier_Removal_Trout[["ci.lb"]], 
                                          Channel_Trout[["ci.lb"]],
                                          Instream_Trout[["ci.lb"]],
                                          Liming_Trout[["ci.lb"]],
                                          "NA",
                                          "NA",
                                          SpawningGravel_Trout[["ci.lb"]]),
                                ci.ub = c(Barrier_Removal_Trout[["ci.ub"]], 
                                          Channel_Trout[["ci.ub"]],
                                          Instream_Trout[["ci.ub"]],
                                          Liming_Trout[["ci.ub"]],
                                          "NA",
                                          "NA",
                                          SpawningGravel_Trout[["ci.ub"]]),
                                p = c(Barrier_Removal_Trout[["pval"]], 
                                      Channel_Trout[["pval"]],
                                      Instream_Trout[["pval"]],
                                      Liming_Trout[["pval"]],
                                      "NA",
                                      "NA",
                                      SpawningGravel_Trout[["pval"]])) #Dataframe of trout responses to restoration type

Salmon_Trout_Effect_Size_Restoration <- rbind(Salmon_restoration, Trout_restoration) #Combine restoration effect sizes for both salmon and trout
Salmon_Trout_Effect_Size_Restoration$yi <- as.numeric(as.character(Salmon_Trout_Effect_Size_Restoration$yi)) #Set as numeric
Salmon_Trout_Effect_Size_Restoration$ci.lb <- as.numeric(as.character(Salmon_Trout_Effect_Size_Restoration$ci.lb)) #Set as numeric
Salmon_Trout_Effect_Size_Restoration$ci.ub <- as.numeric(as.character(Salmon_Trout_Effect_Size_Restoration$ci.ub)) #Set as numeric
Salmon_Trout_Effect_Size_Restoration$order <- c(9,7,3,13,11,5,1,10,8,4,14,12,6,2) #Order effect sizes for plotting purposes
Salmon_Trout_Effect_Size_Restoration <- 
  Salmon_Trout_Effect_Size_Restoration %>% arrange(-desc(order)) %>%
  filter(restoration!="Riparian restoration") #Order effect sizes for plotting purposes and remove riparian restoration (neither effect size for salmon or trout)

restor_species_plot <- Salmon_Trout_Effect_Size_Restoration %>%
  mutate(restoration = fct_reorder(restoration, order)) %>%
  ggplot(aes(x = restoration, y = yi, ymin = ci.lb, ymax = ci.ub)) +  
  geom_pointrange(aes(color = species, group = species), 
                  position = position_dodge(width = 0.5),
                  size = 1) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limit = c(-3,9),
                     breaks = c(-3,0,3,6,9)) +
  xlab("Restoration Type") +
  ylab("Effect Size (95% CI)") + 
  scale_fill_discrete(labels = c("Atlantic salmon", "Brown trout")) +
  theme(axis.text.x = element_text(size = 20, color = "black"), 
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 25, vjust = 0),
        axis.title.y = element_text(size = 25),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, "cm")) +
  scale_x_discrete(labels=c("Barrier removal" = "Barrier\nremoval\n(4,13)", 
                            "Channel restoration" = "Channel\nrestoration\n(2,3)",
                            "In-stream structures" = "In-stream\nstructures\n(11,27)",
                            "Liming" = "Liming\n(11,2)",
                            "Nutrient restoration" = "Nutrient\nrestoration\n(2,0)",
                            "Spawning gravel" = "Spawning\ngravel\n(2,15)")) #Plot species-specific responses to restoration type 
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
res_barrier <- rma.uni(yi, vi, subset = (restoration=="Barrier_removal"), 
                       data = effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_barrier) #Print model results
qqnorm.rma.uni(res_barrier, envelope = FALSE) #Model validation - normality of residuals

Channel_random_effect_model <- rma.mv(yi, 
                                      vi, 
                                      method = "REML",
                                      subset = (restoration=="Channel"),
                                      random = ~ 1| study_ID,
                                      data = effect_sizes) #Multi-level random-effects model with REML for channel restoration
summary.rma(Channel_random_effect_model) #Print model results
res_channel <- rma.uni(yi, vi, subset = (restoration=="Channel"), 
                       data = effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_channel) #Print model results
qqnorm.rma.uni(res_channel, envelope = FALSE) #Model validation - normality of residuals

Instream_random_effect_model <- rma.mv(yi, 
                                       vi, 
                                       method = "REML",
                                       subset = (restoration=="In-stream_structure"),
                                       random = ~ 1| study_ID,
                                       data = effect_sizes) #Random-effects model with REML estimator for in-stream structures
summary.rma(Instream_random_effect_model) #Print model results
res_instream <- rma.uni(yi, vi, subset = (restoration=="In-stream_structure"), 
                        data = effect_sizes, method = "REML") #Random-effects model with REML estimator to give I^2 statistic and test for normality
summary.rma(res_instream) #Print model results
qqnorm.rma.uni(res_instream, envelope = FALSE) #Model validation - normality of residuals

Liming_random_effect_model <- rma.mv(yi, 
                                     vi, 
                                     method = "REML",
                                     subset = (restoration=="Liming"),
                                     random = ~ 1| study_ID,
                                     data = effect_sizes) #Random-effects model with REML estimator for liming
summary.rma(Liming_random_effect_model) #Print model results
res_liming <- rma.uni(yi, vi, subset = (restoration=="Liming"), 
                      data = effect_sizes, method = "REML") #Random-effects model with REML estimator to give I^2 statistic and test for normality
summary.rma(res_liming) #Print model results
qqnorm.rma.uni(res_liming, envelope = FALSE) #Model validation - normality of residuals

Nutrient_random_effect_model <- rma.mv(yi, 
                                       vi, 
                                       method = "REML",
                                       subset = (restoration=="Nutrient"),
                                       random = ~ 1| study_ID,
                                       data = effect_sizes) #Random-effects model with REML estimator for nutrient addition
summary.rma(Nutrient_random_effect_model) #Print model results
res_nutrient <- rma.uni(yi, vi, subset = (restoration=="Nutrient"), 
                        data = effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_nutrient) #Print model results
qqnorm.rma.uni(res_nutrient, envelope = FALSE) #Model validation - normality of residuals

Riparian_random_effect_model <- rma.mv(yi, 
                                       vi, 
                                       method = "REML",
                                       subset = (restoration=="Riparian"),
                                       random = ~ 1| study_ID,
                                       data = effect_sizes) #Random-effects model with REML estimator for riparian restoration
summary.rma(Riparian_random_effect_model) #Print model results

SpawningGravel_random_effect_model <- rma.mv(yi, 
                                             vi, 
                                             method = "REML",
                                             subset = (restoration=="Spawning_gravel"),
                                             random = ~ 1| study_ID,
                                             data = effect_sizes) #Random-effects model with REML estimator for spawning gravel
summary.rma(SpawningGravel_random_effect_model) #Print model results
res_spawning <- rma.uni(yi, vi, subset = (restoration=="Spawning_gravel"), 
                        data = effect_sizes, method = "REML") #Random-effects model with REML estimator to give I^2 statistic and test for normality
summary.rma(res_spawning) #Print model results
qqnorm.rma.uni(res_spawning, envelope = FALSE) #Model validation - normality of residuals

####Sensitivity analysis restoration type
####Barrier
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
barrier <- subset(effect_sizes,restoration=="Barrier_removal")
barrier_effect_sizes <- escalc("SMD", 
                               m2i = before_mean,  
                               sd2i = before_SD,  
                               n2i = before_sample_size, 
                               m1i = after_mean, 
                               sd1i = after_SD, 
                               n1i = after_sample_size, 
                               data = barrier)
Barrier_Removal_random_effect_model <- rma.mv(yi, 
                                              vi, 
                                              method = "REML",
                                              random = ~ 1| study_ID,
                                              data = barrier_effect_sizes) #Random-effects model with REML estimator
summary.rma(Barrier_Removal_random_effect_model) #Print model results
#Standardized residuals
rs_barrier <- rstudent.rma.mv(Barrier_Removal_random_effect_model)
rs_barrier <- as.data.frame(rs_barrier)
limit_rs <- 3 
rs_barrier$ES_ID <- as.numeric(1:nrow(rs_barrier)) #Add a new column with the effect size ID number

#Hat values
hat_barrier <- hatvalues.rma.mv(Barrier_Removal_random_effect_model)
hat_barrier<- as.data.frame(hat_barrier)
mean_barrier <- mean(hat_barrier$hat_barrier)
limit_hat_barrier <- 2 * mean_barrier
hat_barrier$ES_ID <- as.numeric(1:nrow(hat_barrier)) #Add a new column with the effect size ID number

barrier_effect_sizes$ES_ID <- as.numeric(1:nrow(barrier_effect_sizes))
sensitivity_res_barrier <- left_join(hat_barrier,rs_barrier, by="ES_ID")%>%
  left_join(barrier_effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat_barrier<limit_hat_barrier)
sensitivity_barrier <- rma.mv(yi, #Effect size from each row in database
                              vi, #Measure of variance from each row in database
                              random = ~ 1| study_ID, #Study ID as random effect
                              data = sensitivity_res_barrier) #Random-effects model with REML estimator
summary.rma(sensitivity_barrier) #Print model results without outliers
summary(Barrier_Removal_random_effect_model) #Print model with outliers
res_barrier <- rma.uni(yi, 
                       vi, 
                       data = sensitivity_res_barrier) #Heterogeneity of sensitivity model
res_barrier

channel <- subset(effect_sizes,restoration=="Channel")
channel_effect_sizes <- escalc("SMD", 
                               m2i = before_mean,  
                               sd2i = before_SD,  
                               n2i = before_sample_size, 
                               m1i = after_mean, 
                               sd1i = after_SD, 
                               n1i = after_sample_size, 
                               data = channel)
Channel_random_effect_model <- rma.mv(yi, # effect size from each row in database
                                      vi, # measure of variance from each row in database
                                      method = "REML",
                                      random = ~ 1| study_ID,
                                      data = channel_effect_sizes) #Random-effects model with REML estimator
summary.rma(Channel_random_effect_model) #Print model results
####Channel
#Standardized residuals
rs_channel <- rstudent.rma.mv(Channel_random_effect_model)
rs_channel <- as.data.frame(rs_channel)
limit_rs <- 3
rs_channel$ES_ID <- as.numeric(1:nrow(rs_channel)) #add a new column with the effect size ID number

#Hat values
hat_channel <- hatvalues.rma.mv(Channel_random_effect_model)
hat_channel<- as.data.frame(hat_channel)
mean_channel <- mean(hat_channel$hat_channel)
limit_hat_channel <- 2 * mean_channel
hat_channel$ES_ID <- as.numeric(1:nrow(hat_channel)) #add a new column with the effect size ID number

channel_effect_sizes$ES_ID <- as.numeric(1:nrow(channel_effect_sizes))
sensitivity_res_channel <- left_join(hat_channel,rs_channel, by="ES_ID")%>%
  left_join(channel_effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat_channel<limit_hat_channel)
sensitivity_channel <- rma.mv(yi, #Effect size from each row in database
                              vi, #Measure of variance from each row in database
                              random = ~ 1| study_ID, #Study ID as random effect
                              data = sensitivity_res_channel) #Random-effects model with REML estimator
summary.rma(sensitivity_channel) #Print model results without outliers
summary(Channel_random_effect_model) #Print model with outliers
res_channel <- rma.uni(yi, 
                       vi, 
                       data = sensitivity_res_channel) #Heterogeneity of sensitivity model
res_channel


instream <- subset(effect_sizes,restoration=="In-stream_structure")
instream_effect_sizes <- escalc("SMD", 
                                m2i = before_mean,  
                                sd2i = before_SD,  
                                n2i = before_sample_size, 
                                m1i = after_mean, 
                                sd1i = after_SD, 
                                n1i = after_sample_size, 
                                data = instream)
Instream_random_effect_model <- rma.mv(yi, # effect size from each row in database
                                       vi, # measure of variance from each row in database
                                       method = "REML",
                                       random = ~ 1| study_ID,
                                       data = instream_effect_sizes) #Random-effects model with REML estimator
summary.rma(Instream_random_effect_model) #Print model results
####Instream structures
#Standardized residuals
rs_instream <- rstudent.rma.mv(Instream_random_effect_model)
rs_instream <- as.data.frame(rs_instream)
limit_rs <- 3
rs_instream$ES_ID <- as.numeric(1:nrow(rs_instream)) #add a new column with the effect size ID number

#Hat values
hat_instream <- hatvalues.rma.mv(Instream_random_effect_model)
hat_instream<- as.data.frame(hat_instream)
mean_instream <- mean(hat_instream$hat_instream)
limit_hat_instream <- 2 * mean_instream
hat_instream$ES_ID <- as.numeric(1:nrow(hat_instream)) #add a new column with the effect size ID number

instream_effect_sizes$ES_ID <- as.numeric(1:nrow(instream_effect_sizes))
sensitivity_res_instream <- left_join(hat_instream,rs_instream, by="ES_ID")%>%
  left_join(instream_effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat_instream<limit_hat_instream)
sensitivity_instream <- rma.mv(yi, #Effect size from each row in database
                               vi, #Measure of variance from each row in database
                               random = ~ 1| study_ID, #Study ID as random effect
                               data = sensitivity_res_instream) #Random-effects model with REML estimator
summary.rma(sensitivity_instream) #Print model results without outliers
summary(Instream_random_effect_model) #Print model with outliers
res_instream <- rma.uni(yi, 
                        vi, 
                        data = sensitivity_res_instream) #Heterogeneity of sensitivity model
res_instream


liming <- subset(effect_sizes,restoration=="Liming")
liming_effect_sizes <- escalc("SMD", 
                              m2i = before_mean,  
                              sd2i = before_SD,  
                              n2i = before_sample_size, 
                              m1i = after_mean, 
                              sd1i = after_SD, 
                              n1i = after_sample_size, 
                              data = liming)
Liming_random_effect_model <- rma.mv(yi, # effect size from each row in database
                                     vi, # measure of variance from each row in database
                                     method = "REML",
                                     random = ~ 1| study_ID,
                                     data = liming_effect_sizes) #Random-effects model with REML estimator
summary.rma(Liming_random_effect_model) #Print model results
####Liming
#Standardized residuals
rs_liming <- rstudent.rma.mv(Liming_random_effect_model)
rs_liming <- as.data.frame(rs_liming)
limit_rs <- 3
rs_liming$ES_ID <- as.numeric(1:nrow(rs_liming)) #add a new column with the effect size ID number

#Hat values
hat_liming <- hatvalues.rma.mv(Liming_random_effect_model)
hat_liming<- as.data.frame(hat_liming)
mean_liming <- mean(hat_liming$hat_liming)
limit_hat_liming <- 2 * mean_liming
hat_liming$ES_ID <- as.numeric(1:nrow(hat_liming)) #add a new column with the effect size ID number

liming_effect_sizes$ES_ID <- as.numeric(1:nrow(liming_effect_sizes))
sensitivity_res_liming <- left_join(hat_liming,rs_liming, by="ES_ID")%>%
  left_join(liming_effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat_liming<limit_hat_liming)
sensitivity_liming <- rma.mv(yi, #Effect size from each row in database
                             vi, #Measure of variance from each row in database
                             random = ~ 1| study_ID, #Study ID as random effect
                             data = sensitivity_res_liming) #Random-effects model with REML estimator
summary.rma(sensitivity_liming) #Print model results without outliers
summary(Liming_random_effect_model) #Print model with outliers
res_liming <- rma.uni(yi, 
                      vi, 
                      data = sensitivity_res_liming) #Heterogeneity of sensitivity model
res_liming


nutrient <- subset(effect_sizes,restoration=="Nutrient")
nutrient_effect_sizes <- escalc("SMD", 
                                m2i = before_mean,  
                                sd2i = before_SD,  
                                n2i = before_sample_size, 
                                m1i = after_mean, 
                                sd1i = after_SD, 
                                n1i = after_sample_size, 
                                data = nutrient)
Nutrient_random_effect_model <- rma.mv(yi, # effect size from each row in database
                                       vi, # measure of variance from each row in database
                                       method = "REML",
                                       random = ~ 1| study_ID,
                                       data = nutrient_effect_sizes) #Random-effects model with REML estimator
summary.rma(Nutrient_random_effect_model) #Print model results
####Nutrient
#Standardized residuals
rs_nutrient <- rstudent.rma.mv(Nutrient_random_effect_model)
rs_nutrient <- as.data.frame(rs_nutrient)
limit_rs <- 3
rs_nutrient$ES_ID <- as.numeric(1:nrow(rs_nutrient)) #add a new column with the effect size ID number

#Hat values
hat_nutrient <- hatvalues.rma.mv(Nutrient_random_effect_model)
hat_nutrient<- as.data.frame(hat_nutrient)
mean_nutrient <- mean(hat_nutrient$hat_nutrient)
limit_hat_nutrient <- 2 * mean_nutrient
hat_nutrient$ES_ID <- as.numeric(1:nrow(hat_nutrient)) #add a new column with the effect size ID number

nutrient_effect_sizes$ES_ID <- as.numeric(1:nrow(nutrient_effect_sizes))
sensitivity_res_nutrient <- left_join(hat_nutrient,rs_nutrient, by="ES_ID")%>%
  left_join(nutrient_effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat_nutrient<limit_hat_nutrient)
sensitivity_nutrient <- rma.mv(yi, #Effect size from each row in database
                               vi, #Measure of variance from each row in database
                               random = ~ 1| study_ID, #Study ID as random effect
                               data = sensitivity_res_nutrient) #Random-effects model with REML estimator
summary.rma(sensitivity_nutrient) #Print model results without outliers
summary(Nutrient_random_effect_model) #Print model with outliers
res_nutrient <- rma.uni(yi, 
                        vi, 
                        data = sensitivity_res_nutrient) #Heterogeneity of sensitivity model
res_nutrient



spawning <- subset(effect_sizes,restoration=="Spawning_gravel")
spawning_effect_sizes <- escalc("SMD", 
                                m2i = before_mean,  
                                sd2i = before_SD,  
                                n2i = before_sample_size, 
                                m1i = after_mean, 
                                sd1i = after_SD, 
                                n1i = after_sample_size, 
                                data = spawning)
SpawningGravel_random_effect_model <- rma.mv(yi, # effect size from each row in database
                                             vi, # measure of variance from each row in database
                                             method = "REML",
                                             random = ~ 1| study_ID,
                                             data = spawning_effect_sizes) #Random-effects model with REML estimator
summary.rma(SpawningGravel_random_effect_model) #Print model results
####Nutrient
#Standardized residuals
rs_spawning <- rstudent.rma.mv(SpawningGravel_random_effect_model)
rs_spawning <- as.data.frame(rs_spawning)
limit_rs <- 3
rs_spawning$ES_ID <- as.numeric(1:nrow(rs_spawning)) #add a new column with the effect size ID number

#Hat values
hat_spawning <- hatvalues.rma.mv(SpawningGravel_random_effect_model)
hat_spawning<- as.data.frame(hat_spawning)
mean_spawning <- mean(hat_spawning$hat_spawning)
limit_hat_spawning <- 2 * mean_spawning
hat_spawning$ES_ID <- as.numeric(1:nrow(hat_spawning)) #add a new column with the effect size ID number

spawning_effect_sizes$ES_ID <- as.numeric(1:nrow(spawning_effect_sizes))
sensitivity_res_spawning <- left_join(hat_spawning,rs_spawning, by="ES_ID")%>%
  left_join(spawning_effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat_spawning<limit_hat_spawning)
sensitivity_spawning <- rma.mv(yi, #Effect size from each row in database
                               vi, #Measure of variance from each row in database
                               random = ~ 1| study_ID, #Study ID as random effect
                               data = sensitivity_res_spawning) #Random-effects model with REML estimator
summary.rma(sensitivity_spawning) #Print model results without outliers
summary(SpawningGravel_random_effect_model) #Print model with outliers
res_spawning <- rma.uni(yi, 
                        vi, 
                        data = sensitivity_res_spawning) #Heterogeneity of sensitivity model
res_spawning



MA_restoration <- data.frame(restoration = c("Barrier\nremoval\n(17)", 
                                             "Channel\nrestoration\n(5)", 
                                             "In-stream\nstructures\n(38)",
                                             "Liming\n(13)", 
                                             "Nutrient\nrestoration\n(2)", 
                                             "Spawning\ngravel\n(17)"),
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
MA_restoration$order <- c(4,3,2,6,5,1) #Order restoration type for plotting
MA_restoration <- 
  MA_restoration %>% arrange(-desc(order)) #Order restoration type for plotting

restoration_sub_group <- MA_restoration %>%
  mutate(restoration = fct_reorder(restoration, order)) %>%
  ggplot(aes(x = restoration, y = yi, ymin = ci.lb, ymax = ci.ub)) +  
  geom_pointrange() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  scale_y_continuous(breaks = c(-2,-1,0,1,2,3), 
                     limit = c(-2,3), expand = c(0,0)) +
  xlab("Restoration Type") +
  ylab("Effect Size (95% CI)") +
  theme(axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 16, vjust = -0.5),
        axis.title.y = element_text(size = 16),
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
                                       data = effect_sizes) #Random-effects model with REML estimator with latitude as a moderator variable
summary.rma(random_effect_model_latitude) #Print model results
res <- rma.uni(yi, 
               vi, 
               mods = ~ Latitude,
               data = effect_sizes) #Random-effects model with REML estimator for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality 

random_effect_model_riverdepth <- rma.mv(yi, 
                                         vi, 
                                         mods = ~ river_depth, 
                                         random = ~ 1| study_ID,
                                         data = effect_sizes) #Random-effects model with REML estimator with river depth as a moderator variable
summary.rma(random_effect_model_riverdepth) #Print model results
res <- rma.uni(yi, 
               vi, 
               mods = ~ river_depth,
               data = effect_sizes) #Random-effects model with REML estimator for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality

random_effect_model_riverwidth <- rma.mv(yi, 
                                         vi, 
                                         mods = ~ river_width, 
                                         random = ~ 1| study_ID,
                                         data = effect_sizes) #Random-effects model with REML estimator with river width as a moderator variable
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
                                      data = effect_sizes) #Random-effects model with REML estimator with reach land-use as a moderator variable
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
                                       data = effect_sizes) #Random-effects model with REML estimator with presence of stocking as a moderator variable
summary.rma(random_effect_model_stocking) #Print model results
res <- rma.uni(yi, # effect size from each row in database
               vi, # measure of variance from each row in database
               mods = ~ stocking,
               data = effect_sizes) #Random-effects model with REML estimator for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality

random_effect_model_studydesign <- rma.mv(yi, 
                                          vi, 
                                          mods = ~ factor(study_design), 
                                          random = ~ 1 | study_ID, #study ID as random effect
                                          data = effect_sizes) #Random-effects model with REML estimator with study design as a moderator variable
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
                                         data = effect_sizes) #Random-effects model with REML estimator with project age as a moderator variable
summary.rma(random_effect_model_projectage) #print model results
res <- rma.uni(yi, 
               vi, 
               mods = ~ post_monitoring_years,
               data = effect_sizes) #Multi-level random-effects model for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality
project_age_ <- effect_sizes %>% filter(yi < 5)
random_effect_model_projectage <- rma.mv(yi, 
                                         vi, 
                                         mods = ~ post_monitoring_years, 
                                         random = ~ 1| study_ID,
                                         data = project_age_) #Random-effects model with REML estimator with project age as a moderator variable
summary.rma(random_effect_model_projectage) #print model results
res <- rma.uni(yi, 
               vi, 
               mods = ~ post_monitoring_years,
               data = project_age_) #Multi-level random-effects model for I^2 and normality testing
res #Print model
qqnorm.rma.uni(res, envelope = FALSE) #Model validation of normality
project_age_plot <- 
  ggplot(project_age_, aes(x = post_monitoring_years, y = yi)) + 
  geom_point() +
  geom_smooth(method = lm) + 
  theme_classic() +
  xlab("Project age (years)") +
  ylab("Effect Size") +
  theme(axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) #Scatterplot with effect size as a function of project age
project_age_plot

restoration_mods <- subset(effect_sizes, restoration!="Riparian") #Exclude riparian (only one study)
random_effect_model_restoration <- rma.mv(yi, 
                                          vi, # 
                                          mods = ~ factor(restoration), 
                                          random = ~ 1 | study_ID, #study ID as random effect
                                          data = restoration_mods) #Random-effects model with REML estimator with restoration as a moderator variable
summary.rma(random_effect_model_restoration) #Print model
summary(glht(random_effect_model_restoration, 
             linfct=cbind(contrMat(rep(1,6), type="Tukey"))), test=adjusted("none")) #Post-hoc Tukey test for pairwise comparisons
res <- rma.uni(yi, # effect size from each row in database
               vi, # measure of variance from each row in database
               mods = ~ factor(restoration),
               data = restoration_mods) #Random-effects model with REML estimator for I^2 and normality testing
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
                                       data = effect_sizes) #Random-effects model with REML estimator with latitude as a moderator variable
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
                               data = sensitivity_res_latitude) #Random-effects model with REML estimator with latitude as a moderator variable
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
                                         data = effect_sizes) #Random-effects model with REML estimator with river depth as a moderator variable
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
                            data = sensitivity_res_depth) #Random-effects model with REML estimator
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
                                         data = effect_sizes) #Random-effects model with REML estimator with river width as a moderator variable
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
                            data = sensitivity_res_width) #Random-effects model with REML estimator
summary.rma(sensitivity_width) #Print model results without outliers
summary(random_effect_model_riverwidth) #Print model with outliers


####Land use
rm(list = ls()) #Clear working environment
MA_data <- read.csv("meta_analysis_salmo_data.csv") #Load in meta-analysis effect size dataset
MA_data <- subset(MA_data, after_sample_size!="1") #Remove projects with one sample size, cannot compute standardised mean difference effect sizes with sample size of one
MA_data$publication_year <- as.numeric(as.character(MA_data$publication_year)) #Set publication year variable to numeric
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
                                      data = effect_sizes) #Random-effects model with REML estimator with reach land-use as a moderator variable

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
                              data = sensitivity_res_landuse) #Random-effects model with REML estimator
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
                                       data = effect_sizes) #Random-effects model with REML estimator with presence of stocking as a moderator variable
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
                               data = sensitivity_res_stocking) #Random-effects model with REML estimator
summary.rma(sensitivity_stocking) #Print model results without outliers
summary(random_effect_model_stocking) #Print model with outliers

####Study design
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
random_effect_model_studydesign <- rma.mv(yi, # effect size from each row in database
                                          vi, # measure of variance from each row in database
                                          mods = ~ factor(study_design), 
                                          random = ~ 1| study_ID,
                                          data = effect_sizes) #Random-effects model with REML estimator with presence of stocking as a moderator variable
#Standardized residuals
rs <- rstudent.rma.mv(random_effect_model_studydesign)
rs <- as.data.frame(rs) #Residuals in dataframe
limit_rs <- 3 #Set limit of standardised deleted residuals to 3
rs$ES_ID <- as.numeric(1:nrow(rs)) #Add a new column with the effect size ID number

#Hat values
hat <- hatvalues.rma.mv(random_effect_model_studydesign)
hat <- as.data.frame(hat)
mean <- mean(hat$hat)
limit_hat <- 2 * mean
hat$ES_ID <- as.numeric(1:nrow(hat)) #add a new column with the effect size ID number

effect_sizes$ES_ID <- as.numeric(1:nrow(effect_sizes))
sensitivity_res_study <- left_join(hat,rs, by="ES_ID")%>%
  left_join(effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat<limit_hat)

sensitivity_study <- rma.mv(yi, #Effect size from each row in database
                            vi, #Measure of variance from each row in database
                            mods = ~ factor(study_design),
                            random = ~ 1| study_ID, #Study ID as random effect
                            data = sensitivity_res_study) #Random-effects model with REML estimator
summary.rma(sensitivity_study) #Print model results without outliers
summary(random_effect_model_studydesign) #Print model with outliers


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
effect_sizes <- effect_sizes %>% filter(yi < 5) #Filter out effect sizes that are larger than 5
effect_sizes <- subset(effect_sizes, post_monitoring_years!="NA") #Filter out NAs
random_effect_model_projectage <- rma.mv(yi, # effect size from each row in database
                                         vi, # measure of variance from each row in database
                                         mods = ~ post_monitoring_years, 
                                         random = ~ 1| study_ID,
                                         data = effect_sizes) #Random-effects model with REML estimator with project age as a moderator variable
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
effect_sizes <- subset(effect_sizes, restoration!="Riparian") 
random_effect_model_restoration <- rma.mv(yi, # effect size from each row in database
                                          vi, # measure of variance from each row in database
                                          mods = ~ factor(restoration), 
                                          random = ~ 1| study_ID,
                                          data = effect_sizes) #Random-effects model with REML estimator with restoration as a moderator variable
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
hat$ES_ID <- as.numeric(1:nrow(hat)) #Add a new column with the effect size ID number

effect_sizes$ES_ID <- as.numeric(1:nrow(effect_sizes))
sensitivity_res_restoration <- left_join(hat,rs, by="ES_ID")%>%
  left_join(effect_sizes, by ="ES_ID") %>%
  filter(resid < limit_rs)%>%
  filter(hat<limit_hat)

sensitivity_restoration <- rma.mv(yi, #Effect size from each row in database
                                  vi, #Measure of variance from each row in database
                                  mods = ~ factor(restoration),
                                  random = ~ 1| study_ID, #Study ID as random effect
                                  data = sensitivity_res_restoration) #Random-effects model with REML estimator
summary.rma(sensitivity_restoration) #Print model results without outliers
summary.rma(random_effect_model_restoration) #Print model results with outliers
summary(glht(sensitivity_restoration, 
             linfct=cbind(contrMat(rep(1,4), type="Tukey"))), test=adjusted("bonferroni")) #Post hoc Tukey test for pairwise comparisons




###Stocking random-effects models
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
                                      data = N_stock_effect_sizes) #Random-effects model with REML estimator for rivers without stocking
summary.rma(N_stock_random_effect_model) #Print model results
res_N_stock <- rma.uni(yi, vi, 
                       data = N_stock_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_N_stock) #Print model results
qqnorm.rma.uni(res_N_stock, envelope = FALSE) #Model validation - normality of residuals

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
                                      data = Y_stock_effect_sizes) #Random-effects model with REML estimator for rivers with stocking
summary.rma(Y_stock_random_effect_model) #Print model results
res_Y_stock <- rma.uni(yi, vi, 
                       data = Y_stock_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(res_Y_stock) #Print model results
qqnorm.rma.uni(res_Y_stock, envelope = FALSE) #Model validation - normality of residuals

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
                                   data = Agri_effect_sizes) #Random-effects model with REML estimator for rivers in agricultural areas
summary.rma(Agri_random_effect_model) #Print model results
Agri_res <- rma.uni(yi, vi, 
                       data = Agri_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(Agri_res) #Print model results
qqnorm.rma.uni(Agri_res, envelope = FALSE) #Model validation - normality of residuals

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
                                      data = Natural_effect_sizes) #Random-effects model with REML estimator for rivers in natural areas
summary.rma(Natural_random_effect_model) #Print model results
Natural_res <- rma.uni(yi, vi, 
                    data = Natural_effect_sizes, method = "REML") #Random-effects model to give I^2 statistic and test for normality
summary.rma(Natural_res) #Print model results
qqnorm.rma.uni(Natural_res, envelope = FALSE) #Model validation - normality of residuals
