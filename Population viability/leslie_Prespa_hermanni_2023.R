library(popbio)
library(ggplot2)
library(dplyr)

# Define populations to be simulated
populations <- c("Plateau", "Beach", "Konjsko")

# Helper function to sample within CIs
sample_from_CI <- function(CI) {
  runif(1, min = CI[1], max = CI[2])
}

# --- Historical Population Data (2009-2022) ---
historical_years <- 2009:2022

# Population Plateau
historical_Plateau_males <- c(915, 904, 905, 897, 890, 889, 866, 871, 757, 818, 812, 789, 755, 734)
historical_Plateau_females <- c(53, 35, 33, 41, 43, 36, 52, 39, 38, 34, 35, 35, 32, 28) 

# Population Beach
historical_Beach_males <- c(55, 56, 55, 54, 69, 58, 86, 71, 46, 59, 55, 53, 48, 44)
historical_Beach_females <- c(23, 11, 15, 12, 17, 22, 20, 22, 14, 23, 24, 25, 21, 15)

# Population Konjsko
historical_Konjsko_males <- c(228, 228, 228, 273, 259, 268, 234, 226, 249, 190, 211, 412, 283, 156)
historical_Konjsko_females <- c(245, 245, 245, 291, 667, 215, 270, 394, 398, 422, 345, 187, 373, 500)


# --- Population Plateau Parameters ---
# Male survival CIs
s0mPlateau <- c(0.19, 0.49)
s1mPlateau <- c(0.19, 0.49)
s2mPlateau <- c(0.41, 0.68)
s3mPlateau <- c(0.64, 0.85)
s4mPlateau <- c(0.89, 0.99)
s5mPlateau <- c(0.93, 0.98)
s6mPlateau <- c(0.95, 0.99)
s7mPlateau <- c(0.95, 0.99)
s8mPlateau <- c(0.95, 0.99)
s9mPlateau <- c(0.92, 0.99)
s10mPlateau <- c(0.92, 0.99)

# Female survival CIs
s0fPlateau <- c(0.19, 0.49)
s1fPlateau <- c(0.19, 0.49)
s2fPlateau <- c(0.41, 0.68)
s3fPlateau <- c(0.64, 0.85)
s4fPlateau <- c(0.89, 0.99)
s5fPlateau <- c(0.79, 0.91)
s6fPlateau <- c(0.79, 0.91)
s7fPlateau <- c(0.79, 0.91)
s8fPlateau <- c(0.79, 0.91)
s9fPlateau <- c(0.86, 0.92)
s10fPlateau <- c(0.86, 0.92)

# Fecundity and sex ratio
fecundityPlateau <- c(0.00, 0.89)
srPlateau <- 0.5

# Initial population vector (as provided)
N_Plateau <- c(rep(0, 9), 836, rep(0, 9), 38)

# --- Population Beach Parameters ---
# Male survival CIs
s0mBeach <- c(0.19, 0.49)
s1mBeach <- c(0.19, 0.49)
s2mBeach <- c(0.41, 0.68)
s3mBeach <- c(0.64, 0.85)
s4mBeach <- c(0.89, 0.99)
s5mBeach <- c(0.93, 0.98)
s6mBeach <- c(0.95, 0.99)
s7mBeach <- c(0.95, 0.99)
s8mBeach <- c(0.95, 0.99)
s9mBeach <- c(0.90, 0.94)
s10mBeach <- c(0.90, 0.94)

# Female survival CIs
s0fBeach <- c(0.19, 0.49)
s1fBeach <- c(0.19, 0.49)
s2fBeach <- c(0.41, 0.68)
s3fBeach <- c(0.64, 0.85)
s4fBeach <- c(0.89, 0.99)
s5fBeach <- c(0.93, 0.98)
s6fBeach <- c(0.95, 0.99)
s7fBeach <- c(0.95, 0.99)
s8fBeach <- c(0.95, 0.99)
s9fBeach <- c(0.83, 0.90)
s10fBeach <- c(0.83, 0.90)

# Fecundity and sex ratio
fecundityBeach <- c(1.27, 3.02)
srBeach <- 0.5

# Initial population vector (as provided)
N_Beach <- c(rep(0, 9), 47, rep(0, 9), 18)


# --- Population Konjsko Parameters ---
# Male survival CIs
s0mKonjsko <- c(0.19, 0.49)
s1mKonjsko <- c(0.19, 0.49)
s2mKonjsko <- c(0.41, 0.68)
s3mKonjsko <- c(0.64, 0.85)
s4mKonjsko <- c(0.89, 0.99)
s5mKonjsko <- c(0.93, 0.98)
s6mKonjsko <- c(0.95, 0.99)
s7mKonjsko <- c(0.95, 0.99)
s8mKonjsko <- c(0.95, 0.99)
s9mKonjsko <- c(0.89, 0.93)
s10mKonjsko <- c(0.89, 0.93)

# Female survival CIs
s0fKonjsko <- c(0.19, 0.49)
s1fKonjsko <- c(0.19, 0.49)
s2fKonjsko <- c(0.41, 0.68)
s3fKonjsko <- c(0.64, 0.85)
s4fKonjsko <- c(0.89, 0.99)
s5fKonjsko <- c(0.93, 0.98)
s6fKonjsko <- c(0.95, 0.99)
s7fKonjsko <- c(0.95, 0.99)
s8fKonjsko <- c(0.95, 0.99)
s9fKonjsko <- c(0.92, 0.95)
s10fKonjsko <- c(0.92, 0.95)

# Fecundity and sex ratio
fecundityKonjsko <- c(4.43, 6.70)
srKonjsko <- 0.5

# Initial population vector (as provided)
N_Konjsko <- c(rep(0, 9), 260, rep(0, 9), 349)


# Simulation parameters
iterations <- 100
num_sims <- 1000  # Running 1000 simulations for robust result

# Lists to store the median population sizes for later combined plotting
median_male_list <- list()
median_female_list <- list()
years <- 2023:(2023 + iterations - 1)

# List to store final summary results
results_summary <- list()

# Main loop to run simulation for each population
for (pop_name in populations) {
  # Initialize matrices to store results for the current population
  male_pop_matrix <- matrix(NA, nrow = iterations, ncol = num_sims)
  female_pop_matrix <- matrix(NA, nrow = iterations, ncol = num_sims)
  
  # Set parameters based on the current population
  if (pop_name == "Plateau") {
    N <- N_Plateau
    s_m_list <- list(s0mPlateau, s1mPlateau, s2mPlateau, s3mPlateau, s4mPlateau, s5mPlateau, s6mPlateau, s7mPlateau, s8mPlateau, s9mPlateau, s10mPlateau)
    s_f_list <- list(s0fPlateau, s1fPlateau, s2fPlateau, s3fPlateau, s4fPlateau, s5fPlateau, s6fPlateau, s7fPlateau, s8fPlateau, s9fPlateau, s10fPlateau)
    fecundity_CI <- fecundityPlateau
    sr <- srPlateau
  } else if (pop_name == "Beach") {
    N <- N_Beach
    s_m_list <- list(s0mBeach, s1mBeach, s2mBeach, s3mBeach, s4mBeach, s5mBeach, s6mBeach, s7mBeach, s8mBeach, s9mBeach, s10mBeach)
    s_f_list <- list(s0fBeach, s1fBeach, s2fBeach, s3fBeach, s4fBeach, s5fBeach, s6fBeach, s7fBeach, s8fBeach, s9fBeach, s10fBeach)
    fecundity_CI <- fecundityBeach
    sr <- srBeach
  } else if (pop_name == "Konjsko") {
    N <- N_Konjsko
    s_m_list <- list(s0mKonjsko, s1mKonjsko, s2mKonjsko, s3mKonjsko, s4mKonjsko, s5mKonjsko, s6mKonjsko, s7mKonjsko, s8mKonjsko, s9mKonjsko, s10mKonjsko)
    s_f_list <- list(s0fKonjsko, s1fKonjsko, s2fKonjsko, s3fKonjsko, s4fKonjsko, s5fKonjsko, s6fKonjsko, s7fKonjsko, s8fKonjsko, s9fKonjsko, s10fKonjsko)
    fecundity_CI <- fecundityKonjsko
    sr <- srKonjsko
  }
  
  for (sim_run in 1:num_sims) {
    pop_sizes <- matrix(NA, nrow = iterations, ncol = length(N))
    pop_sizes[1, ] <- N
    A <- matrix(0, nrow = 20, ncol = 20)
    
    for (t in 2:iterations) {
      # Introduce environmental stochasticity (5%)
      environmental_effect <- 1 + rnorm(1, mean = 0, sd = 0.05)
      
      # Update Leslie matrix with survival values sampled from CIs and adjust with environmental effect
      for (i in 1:9) {
        A[i + 1, i] <- sample_from_CI(s_m_list[[i+1]]) * environmental_effect
        A[i + 11, i + 10] <- sample_from_CI(s_f_list[[i+1]]) * environmental_effect
      }
      A[10, 10] <- sample_from_CI(s_m_list[[11]]) * environmental_effect
      A[20, 20] <- sample_from_CI(s_f_list[[11]]) * environmental_effect
      
      # Sample fecundity within its CI and adjust with environmental effect
      fec <- sample_from_CI(fecundity_CI) * environmental_effect
      
      # Update fecundity values in Leslie matrix
      for (j in 18:20) {
        A[1, j] <- sr * fec * sample_from_CI(s_m_list[[1]]) * environmental_effect
        A[11, j] <- (1 - sr) * fec * sample_from_CI(s_f_list[[1]]) * environmental_effect
      }
      
      # Project population to the next time step
      pop_sizes[t, ] <- A %*% pop_sizes[t - 1, ]
    }
    
    # Extract total male and female populations from the full matrix for this simulation run
    total_male_pop <- rowSums(pop_sizes[, 1:10])
    total_female_pop <- rowSums(pop_sizes[, 11:20])
    
    # Save the results of this simulation run
    male_pop_matrix[, sim_run] <- total_male_pop
    female_pop_matrix[, sim_run] <- total_female_pop
  }
  
  # Calculate the median for each sex and store it in the lists
  median_male_list[[pop_name]] <- apply(male_pop_matrix, 1, median)
  median_female_list[[pop_name]] <- apply(female_pop_matrix, 1, median)
  
  # --- New Calculations for Reporting ---
  
  # Calculate Stochastic Lambda (lambda_s)
  total_pop_matrix <- male_pop_matrix + female_pop_matrix
  growth_rates <- total_pop_matrix[2:iterations, ] / total_pop_matrix[1:(iterations-1), ]
  stochastic_lambda <- exp(mean(log(growth_rates)))
  
  # Calculate Final Median Population Size
  final_median_male <- median(male_pop_matrix[iterations, ])
  final_median_female <- median(female_pop_matrix[iterations, ])
  
  # Calculate Extinction Probability
  extinction_count <- 0
  for (sim_run in 1:num_sims) {
    if (any(female_pop_matrix[, sim_run] <= 0.4)) {
      extinction_count <- extinction_count + 1
    }
  }
  extinction_prob <- extinction_count / num_sims
  
  # Store all new results in a summary list
  results_summary[[pop_name]] <- list(
    stochastic_lambda = stochastic_lambda,
    final_median_male = final_median_male,
    final_median_female = final_median_female,
    extinction_prob = extinction_prob
  )
}


# --- Data Wrangling for Plotting ---
# Combine all median projected data into one data frame
projected_data_df <- bind_rows(
  data.frame(
    Year = years,
    Population_Size = median_male_list[["Plateau"]],
    Population = "Plateau",
    Sex = "Male"
  ),
  data.frame(
    Year = years,
    Population_Size = median_female_list[["Plateau"]],
    Population = "Plateau",
    Sex = "Female"
  ),
  data.frame(
    Year = years,
    Population_Size = median_male_list[["Beach"]],
    Population = "Beach",
    Sex = "Male"
  ),
  data.frame(
    Year = years,
    Population_Size = median_female_list[["Beach"]],
    Population = "Beach",
    Sex = "Female"
  ),
  data.frame(
    Year = years,
    Population_Size = median_male_list[["Konjsko"]],
    Population = "Konjsko",
    Sex = "Male"
  ),
  data.frame(
    Year = years,
    Population_Size = median_female_list[["Konjsko"]],
    Population = "Konjsko",
    Sex = "Female"
  )
)

# Combine all historical data into a single data frame
historical_data_df <- bind_rows(
  data.frame(
    Year = historical_years,
    Population_Size = historical_Plateau_males,
    Population = "Plateau",
    Sex = "Male"
  ),
  data.frame(
    Year = historical_years,
    Population_Size = historical_Plateau_females,
    Population = "Plateau",
    Sex = "Female"
  ),
  data.frame(
    Year = historical_years,
    Population_Size = historical_Beach_males,
    Population = "Beach",
    Sex = "Male"
  ),
  data.frame(
    Year = historical_years,
    Population_Size = historical_Beach_females,
    Population = "Beach",
    Sex = "Female"
  ),
  data.frame(
    Year = historical_years,
    Population_Size = historical_Konjsko_males,
    Population = "Konjsko",
    Sex = "Male"
  ),
  data.frame(
    Year = historical_years,
    Population_Size = historical_Konjsko_females,
    Population = "Konjsko",
    Sex = "Female"
  )
)


# --- Plotting medians with ggplot2 ---
combined_plot <- ggplot() +
  # Plot projected median lines
  geom_line(data = projected_data_df, aes(x = Year, y = Population_Size, color = Population, linetype = Sex), linewidth = 1.0) +
  # Plot historical data points
  geom_line(data = historical_data_df, aes(x = Year, y = Population_Size, color = Population, linetype = Sex), linewidth = 1.0) +
  scale_x_continuous(breaks = c(2009, 2023, 2030, 2040, 2050, 2060, 2070, 2083, 2090, 2100, 2110, 2123)) +
  
  # Customization
  labs(
    x = "Time (year)",
    y = "Population Size",
    color = "Population",
    linetype = "Sex",
  ) +
  scale_color_manual(values = c("Plateau" = "#cc79a7", "Beach" = "#ecb94c", "Konjsko" = "#a3daf9")) +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed")) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
  geom_vline(xintercept = 2083, color = "black", alpha = 0.5, linewidth=6) +
  theme_minimal() +
  theme(text = element_text(size = 10), axis.title = element_text(size=12), plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = "none"
  )

print(combined_plot)


# --- set a max y to plot a shorter clearer graph in ggplot2 ----
max_y_limit <- bind_rows(projected_data_df, historical_data_df) %>%
  filter(Year <= 2090) %>%
  pull(Population_Size) %>%
  max() * 1.05  # Add a 5% buffer

max_y_limit_2 <- bind_rows(projected_data_df, historical_data_df) %>%
  filter(Year <= 2050) %>%
  pull(Population_Size) %>%
  max() * 1.05  # Add a 5% buffer

# --- plot graph with shorter x and y ---
combined_plot_2090 <- ggplot() +
  # Plot projected median lines
  geom_line(data = projected_data_df, aes(x = Year, y = Population_Size, color = Population, linetype = Sex), linewidth = 1.0) +
  # Plot historical data points
  geom_line(data = historical_data_df, aes(x = Year, y = Population_Size, color = Population, linetype = Sex), linewidth = 1.0) +
  
  # Customization
  labs(
    x = "Time (year)",
    y = "Population Size",
    color = "Population",
    linetype = "Sex",
  ) +
  scale_color_manual(values = c("Plateau" = "#cc79a7", "Beach" = "#ecb94c", "Konjsko" = "#a3daf9")) +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed")) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
  geom_vline(xintercept = 2083, color = "black", alpha = 0.5, linewidth=6) +
  theme_minimal() +
  theme(text = element_text(size = 10), axis.title = element_text(size=12), plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = "none"
  ) +
  scale_x_continuous(limits = c(min(historical_years), 2090), breaks = c(2009, 2023, 2030, 2040, 2050, 2060, 2070, 2083, 2090)) +
  ylim(0, max_y_limit)

print(combined_plot_2090)

# --- plot graph with shorter x and y ---
combined_plot_2090_2 <- ggplot() +
  # Plot projected median lines
  geom_line(data = projected_data_df, aes(x = Year, y = Population_Size, color = Population, linetype = Sex), linewidth = 1.0) +
  # Plot historical data points
  geom_line(data = historical_data_df, aes(x = Year, y = Population_Size, color = Population, linetype = Sex), linewidth = 1.0) +
  
  # Customization
  labs(
    x = "Time (year)",
    y = "Population Size",
    color = "Population",
    linetype = "Sex",
  ) +
  scale_color_manual(values = c("Plateau" = "#cc79a7", "Beach" = "#ecb94c", "Konjsko" = "#a3daf9")) +
  scale_linetype_manual(values = c("Male" = "solid", "Female" = "dashed")) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
  geom_vline(xintercept = 2083, color = "black", alpha = 0.5, linewidth=6) +
  theme_minimal() +
  theme(text = element_text(size = 10), axis.title = element_text(size=12), plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = "none"
  ) +
  scale_x_continuous(limits = c(min(historical_years), 2090), breaks = c(2009, 2023, 2030, 2040, 2050, 2060, 2070, 2083, 2090)) +
  ylim(0, max_y_limit_2)

print(combined_plot_2090_2)

# --- print plots ---
library(Cairo)
ggsave("plot_projections.png", plot = combined_plot, width = 6, height = 4, dpi = 300, type = "cairo")
ggsave("plot_projections.pdf", plot = combined_plot, width = 6, height = 4)
ggsave("plot_projections_2090.png", plot = combined_plot_2090, width = 6, height = 4, dpi = 300, type = "cairo")
ggsave("plot_projections_2090.pdf", plot = combined_plot_2090, width = 6, height = 4)
ggsave("plot_projections_2090_2.png", plot = combined_plot_2090_2, width = 6, height = 4, dpi = 300, type = "cairo")
ggsave("plot_projections_2090_2.pdf", plot = combined_plot_2090_2, width = 6, height = 4)


# --- Print Final Summary of Results ---
cat("\n--- Population Model Summary (from 1,000 simulations) ---\n")

for (pop_name in populations) {
  cat("\n--", pop_name, "--\n")
  cat("Stochastic Growth Rate (Lambda_s):", round(results_summary[[pop_name]]$stochastic_lambda, 4), "\n")
  cat("Final Median Male Population Size (2123):", round(results_summary[[pop_name]]$final_median_male, 2), "\n")
  cat("Final Median Female Population Size (2123):", round(results_summary[[pop_name]]$final_median_female, 2), "\n")
  cat("Probability of Extinction (females <= 0.4):", round(results_summary[[pop_name]]$extinction_prob * 100, 2), "%\n")
}

