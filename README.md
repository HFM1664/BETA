# BETA
BETA DIVERSITY

########Code attempt 1

########## Beta diversity calculation

##  load packages we need for the session
require(vegan)

## read in the data with read.csv
community <- read.csv("DIS/Data/PEAR1.csv", header = TRUE)
meta <- read.csv("DIS/Data/meta_data.csv", header = T)

head(community)

setdiff(meta$Sample.ID, community$Sample.ID)


community <- na.omit(community)
colnames(community)[1] <- "Sample.ID"
meta_names <- meta %>% select (Sample.ID)
community <- left_join(meta_names, community, by = 'Sample.ID')



com_test_site <- adonis2(community[,-1] ~ site, data = meta, method="bray")
print(com_test_site)
## PERMANOVA TEST ##
## run the permanova analysis - here we are running separate tests for each of our factors
## try running a test with multiple factors too if you want to see what that looks like
## site

com_test_xdepth <- adonis2(community[,-1] ~ site$Depth, data = meta, method="bray")

com_test_xdepth1 <- adonis2(community[,-1] ~ Depth, data = meta, method = "bray")
com_test_multifactor <- adonis2(community[,-1] ~ site + Depth, data = meta, method = "bray")


nmds <- metaMDS(community[,-1], distance = "bray")
plot(nmds, type = "t")
ordiellipse(nmds, groups = meta$site, col = 1:3)

summary(meta)  # Provides min, max, mean, and NA counts for all variables
str(meta)      # Displays data types and structure



 



########## Non-Metric Multidimensional Scaling (NMDS)
ord <- metaMDS(community[,-1], distance = "bray")
plot(ord, type = "t")



# Run pairwise comparisons
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
require(cluster)
library(pairwiseAdonis)

# Check metadata structure
head(meta)

# Ensure 'site' is a factor
meta$site <- as.factor(meta$site)
# Exclude Sample.ID and metadata columns to focus on community data
community_matrix <- community[ , !(colnames(community) %in% c("Sample.ID", "site"))]

bray_dist <- vegdist(community_matrix, method = "bray")
pairwise_results <- pairwise.adonis(bray_dist, community$site, sim.method = "bray")

# View results
print(pairwise_results)





pairwise_results <- pairwise.adonis(bray_dist, community$site, sim.method = "bray")




######## Excperimenting with site dependent plot
# Step 1: PERMANOVA for beta diversity by site
com_test_site <- adonis2(community[,-1] ~ site, data = meta, method = "bray")
site_permanova_results <- data.frame(
  Factor = "Site",
  R2 = com_test_site$R2[1],
  p_value = com_test_site$`Pr(>F)`[1]
)

# Optional: Add depth as a second factor if needed
com_test_xdepth <- adonis2(community[,-1] ~ depth, data = meta, method = "bray")
depth_permanova_results <- data.frame(
  Factor = "Depth",
  R2 = com_test_xdepth$R2[1],
  p_value = com_test_xdepth$`Pr(>F)`[1]
)

# Step 2: Pairwise beta diversity tests
bray_dist <- vegdist(community_matrix, method = "bray")
pairwise_results <- pairwise.adonis(bray_dist, meta$site, sim.method = "bray")

# Clean pairwise results for a summary table
pairwise_results_table <- data.frame(
  Comparison = pairwise_results$pairs,
  R2 = pairwise_results$R2,
  p_value = pairwise_results$p.value
)

# Step 3: Combine PERMANOVA and pairwise results into a single table
beta_diversity_table <- list(
  PERMANOVA = rbind(site_permanova_results, depth_permanova_results),
  Pairwise_Comparisons = pairwise_results_table
)

# Step 4: Export tables for review
write.csv(beta_diversity_table$PERMANOVA, "~/DIS/Data/beta_diversity_permanova.csv", row.names = FALSE)
write.csv(beta_diversity_table$Pairwise_Comparisons, "~/DIS/Data/beta_diversity_pairwise.csv", row.names = FALSE)

########Code attempt 2

#### Top taxa coefficient plot
### produced using code from microbiome 

# Load necessary libraries
library(vegan)
library(ggplot2)
library(dplyr)

community <- read.csv("DIS/Data/PEAR1.csv", header = TRUE)
meta <- read.csv("DIS/Data/meta_data.csv", header = T)

# Ensure data is clean and prepared
community <- na.omit(community)
colnames(community)[1] <- "Sample.ID"

# Join the community data with metadata to ensure alignment
community <- left_join(meta %>% select(Sample.ID, site), community, by = "Sample.ID")

# Compute Bray-Curtis dissimilarity
bray_dist <- vegdist(community[,-c(1, 2)], method = "bray")

# Perform constrained ordination (dbRDA)
dbRDA_result <- capscale(bray_dist ~ site, data = community)

coef <- dbRDA_result$CCA$biplot[, 1] # Get the biplot scores for the first axis

top_coef <- sort(head(coef[rev(order(abs(coef)))], 20))

top_coef_df <- data.frame(
  Taxa = factor(names(top_coef), levels = names(top_coef)), # Factor for ordered plotting
  Coefficient = top_coef
)

top_taxa_plot <- ggplot(top_coef_df, aes(x = Coefficient, y = Taxa)) +
  geom_bar(stat = "identity") +
  labs(x = "Coefficient", y = "Taxa", title = "Top 20 Taxa by Coefficients (dbRDA)") +
  theme_minimal()

print(top_taxa_plot)

########Code attempt 3

# Load necessary libraries
library(vegan)
library(ggplot2)
library(dplyr)

community <- read.csv("DIS/Data/PEAR1.csv", header = TRUE)
meta <- read.csv("DIS/Data/meta_data.csv", header = T)


# Clean and prepare data
community <- na.omit(community)
colnames(community)[1] <- "Sample.ID"

# Join community and metadata
community <- left_join(meta %>% select(Sample.ID, site), community, by = "Sample.ID")

# Function to analyze a single site
analyze_site <- function(site_name, community, num_top_taxa = 20) {
  # Subset data for the specific site
  site_data <- community %>% filter(site == site_name)
  
  # Ensure the taxa columns are numeric
  taxa_data <- site_data[,-c(1, 2)]
  taxa_data <- as.data.frame(sapply(taxa_data, as.numeric))
  
  # Check if there is enough variation
  if (nrow(taxa_data) < 2 || all(rowSums(taxa_data) == 0)) {
    warning(paste("Insufficient data for site:", site_name, "- Skipping analysis."))
    return(NULL)
  }
  
  # Compute Bray-Curtis dissimilarity
  bray_dist <- vegdist(taxa_data, method = "bray")
  
  # Perform constrained ordination (dbRDA)
  dbRDA_result <- capscale(bray_dist ~ 1) # No predictors
  
  # Check if biplot coefficients exist
  if (is.null(dbRDA_result$CCA$biplot)) {
    warning(paste("No biplot coefficients found for site:", site_name, "- Skipping."))
    return(NULL)
  }
  
  # Extract coefficients for taxa
  coef <- dbRDA_result$CCA$biplot[, 1] # First axis coefficients
  
  # Ensure coefficients are numeric
  coef <- as.numeric(coef)
  
  # Identify top taxa with the largest absolute coefficients
  top_coef <- sort(head(coef[rev(order(abs(coef)))], num_top_taxa))
  
  # Create a data frame for plotting
  data.frame(
    Taxa = factor(names(top_coef), levels = names(top_coef)),
    Coefficient = top_coef,
    Site = site_name
  )
}

# List of unique sites
sites <- unique(community$site)

# Analyze each site and combine results
results_list <- lapply(sites, function(site) analyze_site(site, community))
results <- do.call(rbind, results_list[!sapply(results_list, is.null)]) # Remove NULL results

# Plot the results if any data exists
if (!is.null(results) && nrow(results) > 0) {
  top_taxa_plot <- ggplot(results, aes(x = Coefficient, y = Taxa, fill = Site)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~Site, scales = "free_y") +
    labs(x = "Coefficient", y = "Taxa", title = "Top Taxa by Site") +
    theme_minimal()
  
  print(top_taxa_plot)
} else {
  print("No valid results to plot.")
}

