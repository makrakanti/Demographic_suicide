library(emmeans)

datasr<-read.csv("M_F_loc.csv", h=T)

model <- glm(cbind(num_males, num_females) ~ locality * year, 
             data = datasr, 
             family = binomial)

# Get the p-values for each term
anova(model, test = "Chisq") 

emmeans_locality <- emmeans(your_glm_model, ~ locality, type = "response")

# Calculate the estimated mean proportion of males for each year
emmeans_year <- emmeans(your_glm_model, ~ year, type = "response")

# Print the results to see the estimated proportions and their confidence intervals
print(emmeans_locality)
print(emmeans_year)