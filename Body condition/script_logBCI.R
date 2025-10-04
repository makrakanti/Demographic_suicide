library(nlme)
library(emmeans)
library(MuMIn)
library(multcomp)

datam<-read.csv(file='datam.csv', h=T)
datam$BCI_log <- log(datam$BCI)
global_model <- lme(log_BCI ~ loc * sex + loc * season + as.factor(year),
                    random = ~1 | ind,
                    data = datam,
                    weights = varIdent(form = ~1 | loc),
                    method = "ML")

# --- Perform automated model selection with dredge ---

options(na.action = "na.fail") # This is required for dredge()
model_selection <- dredge(global_model, extra = "R^2", fixed = "~1 | ind")
print(model_selection)
write.csv(as.data.frame(model_selection), "model_selection_table.csv", row.names = FALSE)

# --- Get the best-fitting model ---
best_model <- get.models(model_selection, subset = 1)[[1]]
print(best_model)

# --- Fit the chosen best-fitting model with REML ---
best_model <- lme(log_BCI ~ loc * sex  + as.factor(year) + season,
                  random = ~1 | ind,
                  data = datam,
                  na.action = na.omit,
                  weights = varIdent(form = ~1 | loc),
                  method = "REML")

# --- Diagnostics for Model Assumptions ---
# Plot 1: Normality of Residuals (Q-Q Plot)
qqnorm(best_model$residuals)
qqline(best_model$residuals)

# Plot 2: Homoscedasticity (Residuals vs. Fitted Plot)
plot(best_model)


# --- Extract model results.
summary(best_model)
model_summary_df <- as.data.frame(summary(best_model)$tTable)
model_summary_df$term <- rownames(model_summary_df)
write.csv(model_summary_df, "best_model_summary.csv")

# After extracting the summary fit the best-fitting model with a new cohort variable that now represents all of the 'loc:sex' interactions
# in order to perform post-hoc tests on a single factor as 'glht' cannot directly handle an interaction term.
datam$cohort <- interaction(datam$loc, datam$sex, sep = " ", lex.ord = TRUE)
best_model_cohort <- lme(log_BCI ~ cohort + season + as.factor(year),
                  random = ~1 | ind,
                  data = datam,
                  na.action = na.omit,
                  weights = varIdent(form = ~1 | loc),
                  method = "REML")
posthoc_comparisons <- glht(best_model_cohort, mcp(cohort = "Tukey"))
summary(posthoc_comparisons)
confint(posthoc_comparisons)

em_model <- emmeans(best_model_cohort, "cohort")
posthoc_table <- pairs(em_model, adjust = "tukey")
posthoc_df <- as.data.frame(posthoc_table)
write.csv(posthoc_df, "cohort_comparisons.csv", row.names = FALSE)
