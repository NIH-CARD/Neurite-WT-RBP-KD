#Load necessary libraries
install.packages("readxl")
library(readxl)
install.packages("dplyr")
library(dplyr)
install.packages("ggplot2")
library(ggplot2)
install.packages("emmeans")
library(emmeans)

#Load data
data <- read_excel("RBPKDFUNCAT3.xlsx")

# Step 1: Summarize by image, well, and group for small dots (intensity per image)
df_small <- data %>%
  dplyr::select(Replicates, Stain, sgRNA, Well, MeanMeanObjIntensityNeuriteSG) %>%
  mutate(Category = interaction(Stain, sgRNA, sep = " - ", lex.order = TRUE))

# Step 2: Summarize by well for medium dots (average intensity per well)
df_medium <- data %>%
  dplyr::group_by(Replicates, Stain, sgRNA, Well) %>%
  summarise(AvgIntensity = mean(MeanMeanObjIntensityNeuriteSG)) %>%
  mutate(Category = interaction(Stain, sgRNA, sep = " - ", lex.order = TRUE))

# Step 3: Summarize by Group and Treatment for large dots (median intensity per replicate, treatment, and group)
df_large <- df_medium %>%
  group_by(Replicates, Stain, sgRNA) %>%
  summarise(AvgIntensity = median(AvgIntensity)) %>%
  mutate(Category = interaction(Stain, sgRNA, sep = " - ", lex.order = TRUE))

# Step 4: Create a factor for Group and Treatment to control ordering
df_small <- df_small %>%
  mutate(Stain = factor(Stain, levels = c("Met", "CHX", "Aniso", "AHA")))
df_small <- df_small %>%
  mutate(Category = factor(paste(sgRNA, Stain, sep = " - "), 
                           levels = as.vector(outer(unique(df_small$sgRNA), 
                                                    levels(df_small$Stain),                                                     paste, sep = " - "))))

df_medium <- df_medium %>%
  mutate(Stain = factor(Stain, levels = c("Met", "CHX", "Aniso", "AHA")))
df_medium <- df_medium %>%
  mutate(Category = factor(paste(sgRNA, Stain, sep = " - "), 
                           levels = as.vector(outer(unique(df_medium$sgRNA), 
                                                    levels(df_medium$Stain), 
                                                    paste, sep = " - "))))
df_large <- df_large %>%
  mutate(Stain = factor(Stain, levels = c("Met", "CHX", "Aniso", "AHA")))
df_large <- df_large %>%
  mutate(Category = factor(paste(sgRNA, Stain, sep = " - "), 
                           levels = as.vector(outer(unique(df_large$sgRNA), 
                                                    levels(df_large$Stain), 
                                                    paste, sep = " - "))))

df_large_overall <- df_medium %>%
  group_by(Stain, sgRNA) %>%
  summarise(AvgIntensity = mean(AvgIntensity, na.rm = TRUE)) %>%
  mutate(Category = interaction(Stain, sgRNA, sep = " - ", lex.order = TRUE))
df_large_overall <- df_large_overall %>%
  mutate(Category = factor(paste(sgRNA, Stain, sep = " - "), 
                           levels = as.vector(outer(unique(df_large_overall$sgRNA), 
                                                    levels(df_large_overall$Stain), 
                                                    paste, sep = " - "))))



# Create a new column for color mapping
df_medium <- df_medium %>%
  mutate(ColorGroup = case_when(
    Stain == "AHA" ~ paste(Stain, sgRNA, sep = "_"),  
    TRUE ~ Stain  
  ))

df_large <- df_large %>%
  mutate(ColorGroup = case_when(
    Stain == "AHA" ~ paste(Stain, sgRNA, sep = "_"),  
    TRUE ~ Stain  
  ))

# Define colors for each Stain and sgRNA combination
color_mapping <- c(
  "Met" = "#005492",        
  "CHX" = "#932191",  
  "Aniso" = "forestgreen", 
  "AHA_NT" = "black",     
  "AHA_TDP-43 KD" = "#D682FF",    
  "AHA_hnRNPA1 KD" = "#0095FF"      
)

# Extract the y-value for the "AHA NT" line
aha_nt_y <- df_large_overall %>%
  filter(Stain == "AHA" & sgRNA == "NT") %>%
  pull(AvgIntensity)

# Calculate the mean intensity for all Met sgRNAs
mean_met_y <- df_medium %>%
  filter(Stain == "Met") %>%
  summarise(mean_y = mean(AvgIntensity, na.rm = TRUE)) %>%
  pull(mean_y)
mean2 <- mean(mean_met_y)

# ANOVA model for comparing AHA NT, TDP-43 KD, and hnRNPA1 KD
anova_model <- aov(AvgIntensity ~ sgRNA, 
                   data = df_medium %>%
                     filter(Stain == "AHA" & sgRNA %in% c("NT", "TDP-43 KD", "hnRNPA1 KD")))

# View the ANOVA table
summary(anova_model)

# Create pairwise comparisons for NT vs TDP-43 KD, NT vs hnRNPA1 KD
emmeans_model <- emmeans(anova_model, pairwise ~ sgRNA)

# View all pairwise comparisons
emmeans_model$contrasts

# Extract specific p-values for NT vs TDP-43 KD, and NT vs hnRNPA1 KD
contrast_results <- contrast(emmeans_model, method = "pairwise", adjust = "tukey")
contrast_results


# Plot the data
ggplot() +
  geom_hline(yintercept = aha_nt_y, linetype = "dotted", color = "black") +
  geom_hline(yintercept = mean2, linetype = "dotted", color = "black") +
  geom_point(data = df_small, aes(x = Category, y = MeanMeanObjIntensityNeuriteSG), 
             color = "darkgray", size = 0.15, position = position_jitter(width = 0.2), alpha = 0.5) +
  geom_point(data = df_medium, aes(x = Category, y = AvgIntensity, color = ColorGroup, alpha = 0), 
             size = 0.75, position = position_jitter(width = 0.2)) +
  geom_point(data = df_large, aes(x = Category, y = AvgIntensity, fill = ColorGroup), 
            color = "black", size = 2, shape = 24) + 
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  geom_errorbarh(data = df_large_overall, 
                 aes(xmin = Category, xmax = Category, y = AvgIntensity), 
                 height = 2, color = "black", size = 10, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        panel.background = element_rect(fill = "white"),   # Background of the plot area
        plot.background = element_rect(fill = "white"),    # Background of the entire plot
        panel.grid.major = element_blank(),                # Remove major grid lines
        panel.grid.minor = element_blank(), 
        plot.margin = margin(10, 10, 10, 10)) +  
  labs(x = "", y = "Intensity (A.U.)") +
  scale_x_discrete(labels = rep(unique(df_medium$sgRNA), length(unique(df_medium$Stain)))) +
  ylim(0, 1000) +
  guides(color = "none") +
  guides(fill = "none") +
  guides(alpha = "none") +
  annotate("segment", x = which(levels(df_large$Category) == "NT - AHA"), 
           xend = which(levels(df_large$Category) == "TDP-43 KD - AHA"), 
           y = max(df_large$AvgIntensity) + 250, yend = max(df_large$AvgIntensity) + 250, 
           color = "black", linetype = "solid") +
  annotate("segment", x = which(levels(df_large$Category) == "NT - AHA"), 
           xend = which(levels(df_large$Category) == "NT - AHA"), 
           y = max(df_large$AvgIntensity), yend = max(df_large$AvgIntensity) + 250, 
           color = "black", linetype = "solid") +
  annotate("segment", x = which(levels(df_large$Category) == "TDP-43 KD - AHA"), 
           xend = which(levels(df_large$Category) == "TDP-43 KD - AHA"), 
           y = max(df_large$AvgIntensity) + 150, yend = max(df_large$AvgIntensity) + 250, 
           color = "black", linetype = "solid") +
  annotate("text", x = (which(levels(df_large$Category) == "NT - AHA") + 
                          which(levels(df_large$Category) == "TDP-43 KD - AHA")) / 2, 
           y = max(df_large$AvgIntensity) + 280, label = paste("0.0064"), 
           size = 3, color = "black") +
  annotate("segment", x = which(levels(df_large$Category) == "NT - AHA"), 
           xend = which(levels(df_large$Category) == "hnRNPA1 KD - AHA"), 
           y = max(df_large$AvgIntensity) + 320, yend = max(df_large$AvgIntensity) + 320, 
           color = "black", linetype = "solid") +
  annotate("segment", x = which(levels(df_large$Category) == "NT - AHA"), 
           xend = which(levels(df_large$Category) == "NT - AHA"), 
           y = max(df_large$AvgIntensity) + 300, yend = max(df_large$AvgIntensity) + 320, 
           color = "black", linetype = "solid") +
  annotate("segment", x = which(levels(df_large$Category) == "hnRNPA1 KD - AHA"), 
           xend = which(levels(df_large$Category) == "hnRNPA1 KD - AHA"), 
           y = max(df_large$AvgIntensity + 30), yend = max(df_large$AvgIntensity) + 320, 
           color = "black", linetype = "solid") +
  annotate("text", x = (which(levels(df_large$Category) == "NT - AHA") + 
                          which(levels(df_large$Category) == "hnRNPA1 KD - AHA")) / 2, 
           y = max(df_large$AvgIntensity) + 350, label = paste("0.6127"), 
           size = 3, color = "black")
