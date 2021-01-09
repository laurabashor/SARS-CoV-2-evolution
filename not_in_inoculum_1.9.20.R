library(tidyverse)
library(openxlsx)
library(readxl)

#read in data
df <- as.data.frame(variant_table_cleaned_combined)

#filter out all data where SARS2 columns have NAs, and name this df.all
df.filtered <-
  df %>%
  filter_at(vars(starts_with("SARS2")), all_vars(is.na(.)))

df.all <- df.filtered

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
df.hamsters <- delete.na.ham(df.hamsters, 2)

#ferrets
delete.na.fer <- function(df.ferrets, n=0) {
  df.ferrets[rowSums(is.na(df.ferrets)) <= n,]
}

df.ferrets <- delete.na.fer(df.ferrets) #for ferrets we want no NAs at all

wb <- createWorkbook("SNVs_not_in_inoculum.xlsx")

addWorksheet(wb, "not_in_inoculum_all")
addWorksheet(wb, "cats")
addWorksheet(wb, "dogs")
addWorksheet(wb, "hamsters")
addWorksheet(wb, "ferrets")

writeData(wb, "not_in_inoculum_all", df.all,borders="all")
writeData(wb, "cats", df.cats,borders="all")
writeData(wb, "dogs", df.dogs,borders="all")
writeData(wb, "hamsters", df.hamsters,borders="all")
writeData(wb, "ferrets", df.ferrets,borders="all")
saveWorkbook(wb, "SNVs_not_in_inoculum.xlsx", overwrite = TRUE)
