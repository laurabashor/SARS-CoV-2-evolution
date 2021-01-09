#load libraries
library(tidyverse)
library(data.table)
library(stringr)
library(openxlsx)
library(readxl)

#raw output of pipeline is variant_table.txt
#save as .csv 
#import variant table.csv
df <- read.csv("variant_table.csv")

#clean up column names
colnames(df) <-sapply(strsplit(names(df), ".variant_alleles.txt"), `[[`, 1)
as.character(names(df))

#add one variant column, bring it to the front, remove additional columns
df$variant <- paste0(df$ref_aa, df$codon_number)
df$variant <- paste0(df$variant, df$alt_aa)

df <- df %>%
  select(variant, everything())

df <- df[, !names(df) %in% c("ref_aa","codon_number","alt_aa")]

#make sure data is ordered by position in genome
df <- arrange(df, position)

##export this as an excel sheet, variant table containing all replicates and all v2 samples

wb <- createWorkbook("variant_table_cleaned_all_replicates.xlsx")
addWorksheet(wb, "variant_table")
writeData(wb, "variant_table", df,borders="all")
saveWorkbook(wb, "variant_table_cleaned_all_replicates.xlsx", overwrite = TRUE)

####combining replicates -- average them, unless one is NA, then don't keep it

#subset df to leave v2 samples we don't need out
#just keeping C70 and C58
df2 <- select(df, !c("C43_v2", "C56_v2", "C57_v2", "C69_v2", "SARS2_P3_v2"))

#delete _v2 string so these samples can be processed with the rest
names(df2) <- gsub("_v2", "", names(df2))

#make data long and skinny
df2 <- melt(setDT(df2), id.vars=c("variant", "position","cds", "N_or_S"), variable.name="dataset_ID")

#make a new column labeling replicates (anything with _R is replicate 2, without it gets a 1)
df2$replicate = str_extract(df2$dataset_ID, "_R")
df2$replicate[is.na(df2$replicate)] <- 1
df2$replicate[df2$replicate == "_R"] <- 2

#now delete _R from sample names because that data is encoded in replicate column
as.character(df2$dataset_ID)
df2$dataset_ID <- gsub("_R", "", df2$dataset_ID)

#calculate mean of variants, if there is an NA, mean function will always output NA (this is what we want! if we didn't detect the variant in one replicate, we don't want to keep the value for the other replicate)
df3 <- df2 %>%
 group_by(variant,position, cds, N_or_S, dataset_ID) %>%
  summarise(mean = mean(value))

#recreate variant table by pivoting wide based on dataset_ID
df_wide <- df3 %>% pivot_wider(id_cols = c(variant, position,cds,N_or_S), 
                              names_from=dataset_ID, 
                              values_from=c(mean),
                              names_sort = T)

df_wide <- df_wide %>% arrange(position)

#rename df_wide to df.all so you can sort all variants by species
df.all <- df_wide

#subset into data frames for each animal
df.cats <- df.all %>% 
  select(variant, position, cds, N_or_S, starts_with("C"))

df.dogs <- df.all %>% 
  select(variant, position, cds, N_or_S, starts_with("D"))

df.hamsters <- df.all %>% 
  select(variant, position, cds, N_or_S, starts_with("H"))

df.ferrets <- df.all %>% 
  select(variant, position, cds, N_or_S, starts_with("F"))

#use a function to delete any rows in which there are NAs for all cats
delete.na <- function(df.cats, n=0) {
  df.cats[rowSums(is.na(df.cats)) <= n,]
}
df.cats <- delete.na(df.cats,5) #5 is max NAs allowed for 6 cats

#do it for dogs
delete.na.dogs <- function(df.dogs, n=0) {
  df.dogs[rowSums(is.na(df.dogs)) <= n,]
}
df.dogs <- delete.na.dogs(df.dogs, 2) #2 is max NAs allowed for 3 dogs

#do it for hamsters
delete.na.ham <- function(df.hamsters, n=0) {
  df.hamsters[rowSums(is.na(df.hamsters)) <= n,]
}
df.hamsters <- delete.na.ham(df.hamsters, 2) #2 is max NAs allowed for 3 hamsters

#ferrets
delete.na.fer <- function(df.ferrets, n=0) {
  df.ferrets[rowSums(is.na(df.ferrets)) <= n,]
}
df.ferrets <- delete.na.fer(df.ferrets) #no NAs wanted at all

#save variant table as excel sheet, with sheets for each species
wb2 <- createWorkbook("variant_table_cleaned_combined.xlsx")
addWorksheet(wb2, "variant_table_all_means")
addWorksheet(wb2, "cats")
addWorksheet(wb2, "dogs")
addWorksheet(wb2, "hamsters")
addWorksheet(wb2, "ferrets")
writeData(wb2, "variant_table_all_means", df_wide,borders="all")
writeData(wb2, "cats", df.cats,borders="all")
writeData(wb2, "dogs", df.dogs,borders="all")
writeData(wb2, "hamsters", df.hamsters,borders="all")
writeData(wb2, "ferrets", df.ferrets,borders="all")
saveWorkbook(wb2, "variant_table_cleaned_combined.xlsx", overwrite = TRUE)
