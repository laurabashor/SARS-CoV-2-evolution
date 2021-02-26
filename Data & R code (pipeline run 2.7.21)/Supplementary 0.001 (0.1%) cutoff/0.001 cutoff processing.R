#load libraries
library(tidyverse)
library(openxlsx)
library(readxl)

#import variant table.csv
df <- read.csv("variant_summary_0.001_cutoff.csv")

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

#now let's filter down to remake Table 2

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

#save as an excel workbook
wb <- createWorkbook("not in inoculum above 0.1%.xlsx")
addWorksheet(wb, "not_in_inoculum_0.001")
addWorksheet(wb, "not_in_inoculum_0.001_table2")
writeData(wb, "not_in_inoculum_0.001", anti, borders="all")
writeData(wb, "not_in_inoculum_0.001_table2", table2, borders="all")
saveWorkbook(wb, "not in inoculum above 0.1%.xlsx", overwrite=TRUE)  

  
  

