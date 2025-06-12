# Libraries
library(tidyverse)
library(ggplot2)
library(snpStats)
library(readr)

data <- read.csv("/Users/alzorrcarri/projects_genyo/ml_farmacogen/data/joined.tsv",
                 sep = "\t", header = TRUE)

# Separate information from the dataframe
class <- data$Class
sample_id <- data$Samples.ID
snps_data <- data[-c(which(colnames(data) == "Class"), which(colnames(data) == "Samples.ID"))]

rownames(snps_data) <- sample_id

# Detect missing values
sum(is.na(snps_data)) # 861 NAs 
sum(rowSums(is.na(snps_data)) > 0) # 271 rows (samples) with NAs
sum(colSums(is.na(snps_data)) > 0) # 65 columns (SNPs) with NAs


# Study NAs by SNP and delete those with a high percentage of missing values
na_by_snps <- colSums(is.na(snps_data))
na_by_snps <- data.frame(Column = names(na_by_snps), NA_Count = na_by_snps)

ggplot(na_by_snps, aes(x = Column, y = NA_Count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Valores NA por columna en snps_data",
       x = "Columnas",
       y = "Cantidad de NA") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Looking the graph, delete columns with mote than 15 NAs
sum(colSums(is.na(snps_data)) >= 15)

snps_data <- snps_data[, colSums(is.na(snps_data)) < 15]

# Write deleted columns in a file
# write(colnames(snps_data)[colSums(is.na(snps_data)) < 15], file = "filteredMutationsByNAs.txt", sep = '\n')

# Idem for samples (rows)
na_by_register <- rowSums(is.na(snps_data))
na_by_register <- data.frame(Row = sample_id, NA_Count = na_by_register)

ggplot(na_by_register, aes(x = Row, y = NA_Count)) +
  geom_bar(stat = "identity", fill = "red", color = "black") +
  labs(title = "Valores NA por fila en snps_data",
       x = "Filas",
       y = "Cantidad de NA") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Delete files with mote than 2 NAs
sum(rowSums(is.na(snps_data)) >= 2)

snps_data <- snps_data[rowSums(is.na(snps_data)) < 2, ]

# Impute the rest of NAs with snpStats
snpsMatrix <- snps_data
snpsMatrix[snpsMatrix == "0/0"] <- 0
snpsMatrix[snpsMatrix == "0/1"] <- 1
snpsMatrix[snpsMatrix == "1/1"] <- 2
snpsMatrix[snpsMatrix == "./."] <- NA

snpsMatrix <- apply(snpsMatrix, 2, function(x) as.numeric(x))

sum(is.na(snpsMatrix)) # 47 cells to impute

snpsMatrix <- as(as.matrix(snpsMatrix), "SnpMatrix")
rownames(snpsMatrix@.Data) <- rownames(snps_data)
dim(snpsMatrix@.Data)

snpsMatrix@.Data[which(is.na(snps_data))] # NAs coded as "00"

select <- colSums(is.na(snps_data)) > 0
training <- rowSums(is.na(snps_data)) == 0
missing <- as(snpsMatrix@.Data[training, select], "SnpMatrix")
present <- as(snpsMatrix@.Data[training,!select], "SnpMatrix")

target <- as(snpsMatrix@.Data[!training, !select], "SnpMatrix")

# Information about the selected SNPs
snps_selected <- colnames(snps_data)[select]
snps_metadata <- read_tsv("/Users/alzorrcarri/projects_genyo/ml_farmacogen/data/rs_info.tsv",
                          col_names = TRUE)
colnames(snps_metadata) <- c("Chr", "Pos", "ID", "Ref", "Alt", "Gene")
snps_selected_metadata <- snps_metadata %>% filter(ID %in% snps_selected)
pos.missing <- snps_selected_metadata$Pos

snps_selected <- colnames(snps_data)[!select]
snps_selected_metadata <- snps_metadata %>% filter(ID %in% snps_selected)
pos.present <- snps_selected_metadata$Pos

# Imputation rules
rules <- snp.imputation(present, missing, pos.present, pos.missing)
rules
rules[colnames(missing)[c(1, 2)]]
summary(rules)
plot(rules)
write(print(rules), file="reglas.txt")

# Missing values imputation
imputed <- impute.snps(rules, target, as.numeric=TRUE)
imputed_rounded <- round(imputed)
imputed_rounded <- as.data.frame(apply(imputed_rounded, 2, as.integer))
rownames(imputed_rounded) <- rownames(target)

# Correctly formatting the dataframe
imputed_rounded[imputed_rounded == 0] <- "0/0"
imputed_rounded[imputed_rounded == 1] <- "0/1"
imputed_rounded[imputed_rounded == 2] <- "1/1"

# Substitute the imputed value from imputed_rounded df into the NAs in snps_data df
snps_data_imputed <- snps_data
na_positions <- which(is.na(snps_data_imputed), arr.ind = TRUE)
for (i in 1:nrow(na_positions)) {
  row <- na_positions[i, 1]
  col <- na_positions[i, 2]
  row_name <- rownames(snps_data_imputed)[row]
  col_name <- colnames(snps_data_imputed)[col]  
  snps_data_imputed[row, col] <- imputed_rounded[row_name, col_name]
}

# Save the imputed matrix in a .tsv file
data_aux <- data[, colSums(is.na(data)) < 15]
data_aux <- data_aux[rowSums(is.na(data_aux)) < 2, ]

new_data <- snps_data_imputed
new_data$class <- data_aux$Class
new_data$sample_id <- data_aux$Samples.ID

write_tsv(new_data, file = "nas_imputed_joined.tsv")


