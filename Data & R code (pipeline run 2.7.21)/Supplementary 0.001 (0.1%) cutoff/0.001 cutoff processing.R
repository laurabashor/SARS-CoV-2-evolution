#load libraries
library(tidyverse)
library(openxlsx)
library(readxl)

#read in pipeline output variant summary table
#pipeline was run with a 0.1% frequency cutoff instead of 3% cutoff for calling variants
df <- read.xlsx("variant_summary_0.001_cutoff.xlsx")

####combining replicates -- average them, unless one is NA, then don't keep it

#make data long and skinny
#names(df)
df2 <- df %>%
  pivot_longer(!c(reference_sequence, position, gene, codon,indel,variant, reference_base, variant_base, effect,featureid), 
               names_to= "dataset_ID", 
               values_to ="frequency")

#make a new column labeling replicates (anything with _R is replicate 2, without it gets a 1)
df2$replicate = str_extract(df2$dataset_ID, "_R")
df2$replicate[is.na(df2$replicate)] <- 1
df2$replicate[df2$replicate == "_R"] <- 2

#important note: Cat 5 (Cat 70 originally amplified with ARTIC v2 primers) is Cat_5 and there is no Cat_5_R
#also, Cat 4 (Cat 58) has one replicate sequenced with ARTIC v2 primers (Cat_4) and one with ARTIC v3 primers (Cat_4_R)
#all other samples, both replicates were sequenced with ARTIC v3 primers

#now delete _R from sample names because that data is encoded in replicate column
df2$dataset_ID <- gsub("_R", "", df2$dataset_ID)

#right now, there is a 0 if we didn't detect the variant at >3% frequency, and an NA if we didn't have enough read coverage to know if we actually didn't detect it, or it didn't show up due to low coverage 
#for the purpose of combining replicates, we will treat the 0s and NAs the same, as NA
df2[df2 == "0"] <- NA

#what we want is to throw out any data where we didn't detect the variant in both replicates, whatever the reason for this may be
#this calculates mean of variants, if there is an NA in either value, the mean function will always output NA 
df3 <- df2 %>%
  group_by(reference_sequence, position, gene,codon, indel,variant, reference_base, variant_base, effect,featureid,dataset_ID) %>%
  summarise(mean = mean(frequency))

#now we can recreate the variant table by pivoting wider based on dataset_ID
df_wide <- df3 %>% pivot_wider(id_cols = c(reference_sequence, position, gene, indel,variant, reference_base, variant_base, effect), 
                               names_from=dataset_ID, 
                               values_from=c(mean),
                               names_sort = T)

#rename df_wide to df.all before sorting variants by species
df.all <- df_wide

#for our final excel spreadsheet, it will be nice to be able to look at species separately
#subset into data frames for each animal
df.cats <- df.all %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, starts_with("C"))

df.dogs <- df.all %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, starts_with("D"))

df.hamsters <- df.all %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, starts_with("H"))

df.ferrets <- df.all %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect,starts_with("F"))

df.inoc <- df.all %>%
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect,starts_with("P"))

#use a function to delete any rows in which there are NAs for all cats
delete.na <- function(df.cats, n=0) {
  df.cats[rowSums(is.na(df.cats)) <= n,]
}
df.cats <- delete.na(df.cats,5) #5 is max NAs allowed for 6 cats

# #do it for dogs
delete.na.dogs <- function(df.dogs, n=0) {
  df.dogs[rowSums(is.na(df.dogs)) <= n,]
}
df.dogs <- delete.na.dogs(df.dogs, 2) #2 is max NAs allowed for 3 dogs

#do it for hamsters
delete.na.ham <- function(df.hamsters, n=0) {
  df.hamsters[rowSums(is.na(df.hamsters)) <= n,]
}
df.hamsters <- delete.na.ham(df.hamsters, 2) #2 is max NAs allowed for 3 hamsters

#and do it for inoculums
delete.na.inoc <- function(df.inoc, n=0) {
  df.inoc[rowSums(is.na(df.inoc)) <= n,]
}
df.inoc <- delete.na.inoc(df.inoc, 2) #2 is max NAs allowed for 3 passages

#no NAs wanted at all for ferrets
df.ferrets <- na.omit(df.ferrets)

#filter  out just cats 5 & 6 for a table

df.contact <- df.all %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, Cat_5, Cat_6)

delete.na.contact <- function(df.contact, n=0) {
  df.contact[rowSums(is.na(df.contact)) <= n,]
}
df.contact <- delete.na.contact(df.contact, 1)


##now filter for variants that were not detected in the inoculum viral stock across the three passages
#filter to keep all data where Passage 1, 2 and 3 viral stock columns have NAs, and name this df.all
df.filtered <-
  df.all %>%
  filter_at(vars(starts_with("Passage")), all_vars(is.na(.)))

#save variant table as excel sheet, with sheets for each species
wb <- createWorkbook("variant_summary_0.001_processed.xlsx")
addWorksheet(wb, "variant_table_all_means")
addWorksheet(wb, "cats")
addWorksheet(wb, "dogs")
addWorksheet(wb, "hamsters")
addWorksheet(wb, "ferrets")
addWorksheet(wb, "viral_stock_inoculum")
addWorksheet(wb, "all_variants_not_in_inoculum")
addWorksheet(wb, "contact_cats")
writeData(wb, "variant_table_all_means", df_wide,borders="all")
writeData(wb, "cats", df.cats,borders="all")
writeData(wb, "dogs", df.dogs,borders="all")
writeData(wb, "hamsters", df.hamsters,borders="all")
writeData(wb, "ferrets", df.ferrets,borders="all")
writeData(wb, "viral_stock_inoculum", df.inoc, borders="all")
writeData(wb, "all_variants_not_in_inoculum", df.filtered,borders="all")
writeData(wb, "contact_cats",df.contact, borders="all")
saveWorkbook(wb, "variant_summary_0.001_processed.xlsx", overwrite = TRUE)


###further analysis to compare results of 0.1% cutoff pipeline run to 3% cutoff
#variant table was imported as df
#make zeros NAs
df[df == "0"] <- NA

#make a data frame with just inoculum columns
df.inoc <- df %>%
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect,starts_with("P"))

#use a function to delete any rows in which there are NAs for all passage replicates

delete.na.inoc <- function(df.inoc, n=0) {
  df.inoc[rowSums(is.na(df.inoc)) <= n,]
}
df.inoc <- delete.na.inoc(df.inoc, 5) #5 is max NAs allowed for 3 passages x 2 replicates

#now we want to compare the variants in this dataframe to variants we are including in our analysis that were present at >3% but not in the inoculum
#this has been saved as a worksheet in the variant_summary_processed.xlsx file output after initial processing
#save that sheet as a csv and load it in here 

df.not <- read.csv("not_in_inoculum.csv")

#for some reason R doesn't like the indel column
df.not <- df.not %>%
  select(-indel)

##compare df.not to df.inoc
anti <- anti_join(df.not, df.inoc, by="variant") 
#64 variants are not found in the inoculum above 0.1%
#including notable variants like N501T

#now let's filter down to remake Table 2 for variants not detected above 0.1%

#get a df with all variants greater than 50%
anti_long <- anti %>%
  pivot_longer(!c(reference_sequence, position, gene,variant, reference_base, variant_base, effect), 
                names_to= "dataset_ID", 
                values_to ="frequency")

anti0.5 <- anti_long %>%
  filter(frequency >= 0.5)
  
table2 <- anti0.5 %>%
  pivot_wider(id_cols = c(reference_sequence, position, gene, variant, reference_base, variant_base, effect), 
              names_from=dataset_ID, 
              values_from=c(frequency),
              names_sort = T)

#now get dfs with only variants present in every single individual of a species
df.cats <- na.omit(df.cats)
df.dogs <- na.omit(df.dogs)
df.hamsters <- na.omit(df.hamsters)
#none of these variants qualify! meaning only variants in anti0.5 are variants not found in the inoculum above 0.1%

#that means out of the variants in Table 2, only E195G, N501T, D614G and S43Y were not detected in the inoculum above 0.1%
#de novo?? maybe, maybe not

  
  

