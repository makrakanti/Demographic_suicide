data_age<-read.csv(file="data_age.csv", h=T)

library(lme4)
library(DHARMa)
library(ggplot2)

bm <- glm(Eggs ~ BM, family = poisson, data = data_age)
scl <- glm(Eggs ~ SCL,   family = poisson, data = data_age)
Age <- glm(Eggs ~ Age,   family = poisson, data = data_age)
Age_locality <- glm(Eggs ~ Age + Locality, family = poisson, data = data_age)
Age_Xlocality <- glm(Eggs ~ Age * Locality, family = poisson, data = data_age)
bm_Xlocality <- glm(Eggs ~ BM * Locality,   family = poisson, data = data_age)
bm_locality <- glm(Eggs ~ BM + Locality,  family = poisson, data = data_age)
scl_locality <- glm(Eggs ~ SCL + Locality,      family = poisson, data = data_age)
scl_Xlocality <- glm(Eggs ~ SCL * Locality, family = poisson, data = data_age)

# perform model selection
library(MuMIn)
model.set <- model.sel(bm, scl, bm_locality, scl_locality, scl_Xlocality, bm_Xlocality, Age, Age_locality, Age_Xlocality)
write.csv(as.data.frame(model.set), "model_selection_table.csv", row.names = TRUE)

# check diagnostics of best model
simulationOutput <- simulateResiduals(scl)
plot(simulationOutput)

#change family and include locality as random factor to address underdispersion
library(glmmTMB)

com_poisson_bm <- glmmTMB(Eggs ~ BM+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_scl <- glmmTMB(Eggs ~ SCL+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_Age <- glmmTMB(Eggs ~ Age+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_Age_Locality <- glmmTMB(Eggs ~ Age + Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_Age_XLocality <- glmmTMB(Eggs ~ Age * Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_bm_XLocality <- glmmTMB(Eggs ~ BM * Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_bm_Locality <- glmmTMB(Eggs ~ BM + Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_scl_Locality <- glmmTMB(Eggs ~ SCL + Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_scl_XLocality <- glmmTMB(Eggs ~ SCL * Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)

#scale predictors due to numeric troubles
data_age$SCL_scaled <- scale(data_age$SCL)
data_age$BM_scaled <- scale(data_age$BM)

com_poisson_bm <- glmmTMB(Eggs ~ BM_scaled+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_scl <- glmmTMB(Eggs ~ SCL_scaled+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_Age <- glmmTMB(Eggs ~ Age+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_Age_Locality <- glmmTMB(Eggs ~ Age + Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_Age_XLocality <- glmmTMB(Eggs ~ Age * Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_bm_XLocality <- glmmTMB(Eggs ~ BM_scaled * Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_bm_Locality <- glmmTMB(Eggs ~ scale(BM) + Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_scl_Locality <- glmmTMB(Eggs ~ SCL_scaled + Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)
com_poisson_scl_XLocality <- glmmTMB(Eggs ~ SCL_scaled * Locality+ (1 | Locality), family = compois(link = "log"), data = data_age)

#new model selection
model.set <- model.sel(com_poisson_bm, com_poisson_scl, com_poisson_bm_Locality, com_poisson_scl_Locality, com_poisson_scl_XLocality, com_poisson_bm_XLocality, com_poisson_Age, com_poisson_Age_Locality, com_poisson_Age_XLocality)
write.csv(as.data.frame(model.set), "com_poisson_model_selection_table.csv", row.names = TRUE)

#model diagnostics of best-fit model
com_ps_simulationOutput <- simulateResiduals(com_poisson_bm_Locality)
plot(com_ps_simulationOutput)

#to predict from best model
# Create a data frame with every combination of body mass and locality
bm_range <- seq(from = 500, to = 3000, by = 50)  
all_localities <- unique(data_age$Locality)
newdata <- expand.grid(
  BM = bm_range,
  Locality = all_localities
)
#predict with SE
predictions_se <- predict(com_poisson_bm_Locality,
                          newdata = newdata,
                          type = "response",
                          re.form = NA,
                          se.fit = TRUE)

# Add the predictions and standard errors to the newdata data frame and calculate CIs
newdata$predicted_clutch <- predictions_se$fit
newdata$se_fit <- predictions_se$se.fit
newdata$lower_ci <- newdata$predicted_clutch - 1.96 * newdata$se_fit
newdata$upper_ci <- newdata$predicted_clutch + 1.96 * newdata$se_fit

#plot predictions and extract
plot<-ggplot(newdata, aes(x = BM, y = predicted_clutch, color = Locality, fill = Locality)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("#e69f00", "#56b4e9", "#cc79a7")) +  # apply custom line colors
  scale_fill_manual(values = c("#e69f00", "#56b4e9", "#cc79a7")) + 
  coord_cartesian(ylim = c(0, 10)) +
  theme_minimal() +
  labs(x = "Body Mass (g)",
       y = "Predicted egg numbers")

plot

library(Cairo)
ggsave("BM-eggs.png", plot = plot, width = 6, height = 4, dpi = 300, type = "cairo")
ggsave("BM-eggs.pdf", plot = plot, width = 6, height = 4)
