# Load libraries
library(tidyverse)

# Define directories
dir <- "/Users/alzorrcarri/projects_genyo/ml_farmacogen/data/snps_filtered/"  # Carpeta con los archivos filtrados (.tsv)
output_file <- paste0(dir, "joined.tsv")
# Get list of .tsv files in folder
files <- list.files(dir, pattern = "*.tsv", full.names = TRUE)
files

# Load dataframes
file_path <- files[1]
file_name <- basename(file_path)
df1 <- read_tsv(file_path)
df1 <- as.data.frame(df1)


file_path <- files[2]
file_name <- basename(file_path)
df2 <- read_tsv(file_path)
df2 <- as.data.frame(df2)


file_path <- files[3]
file_name <- basename(file_path)
df3 <- read_tsv(file_path)
df3 <- as.data.frame(df3)

df <- cbind(df1, df2, df3)


rownames(df) <- df$ID
df <- df[, -which(names(df) == 'ID')]

# Transpose matrix
df <- as.data.frame(t(df))

# Add class type information
df$Class <- c(rep("Good", times = dim(df1)[2]-1), 
              rep("Control", times = dim(df2)[2]-1), 
              rep("Bad", times = dim(df3)[2]-1))

# Error check and fix in data
df <- df %>%
  mutate(across(where(is.character), ~ gsub("\\|", "/", .)))

df <- df %>%
  mutate(across(where(is.character), ~ na_if(., "./.")))

# Change from character to factor type
sum(apply(df, 2, is.character)) == dim(df)[2]

df <- df %>%
  mutate(across(where(is.character), as.factor))

summary(df$rs3795261)

# Add sample information
df$Samples.ID <- rownames(df)

# Save final table
write_tsv(df, output_file)
