##visualizing TC variant emergence and reversion

#load libraries
library(tidyverse)
library(data.table)

#load in cleaned data with replicates combined
df <- as.data.frame(variant_table_cleaned_combined)

#NAs to 0s for plotting
df[is.na(df)] <- 0

#first graph just the inoculums
df1 <- df[, c("variant", "cds", "SARS2_P1","SARS2_P2","SARS2_P3")]

#rename passages as 1, 2, 3
names(df1) <- gsub("SARS2_P", "", names(df1))

#pull out just TC mutations that went to fixation
df1 <- df1 %>%
  filter(variant %in% c("D135E", "D215H", "R685H", "T7I","S194T"))

#make data long and skinny for plotting
df1.long <- df1 %>%
  pivot_longer(!c(variant, cds), names_to= "passage", values_to ="frequency")

#plot
pdf(file="TC_variants.pdf")
ggplot(df1.long, aes(x=passage, y=frequency, group=variant, color= variant)) +
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=FALSE) +
  facet_wrap(vars(variant))
dev.off()

##want to show the reversion of these variants in vivo 
#working with just inoculums and dogs for now
df2 <- df[, c("variant", "cds", "D46", "D47", "D48", "SARS2_P1","SARS2_P2","SARS2_P3")]

#rename passages as 1, 2, 3
names(df2) <- gsub("SARS2_P", "", names(df2))

#pull out just TC mutations that went to fixation
df2 <- df2 %>%
  filter(variant %in% c("D135E", "D215H", "R685H", "T7I","S194T"))

#make all dogs into one passage
names(df2) <- gsub("D46|D47|D48","Dogs",names(df2))

#make data long and skinny for plotting, rename value as frequency
#can do this with data.table melt or tidyr pivot_longer

df2.long <- df2 %>%
  pivot_longer(!c(variant, cds), names_to= "passage", values_to ="frequency")

df2.long <- melt(setDT(df2), id.vars=c("variant", "cds"), variable.name="passage", value.name="frequency")

#get the order right for passages
df2.long <- df2.long %>%mutate(passage = fct_relevel(passage, "1","2","3", "Dogs"))

#plot
pdf(file="TC_variants_and_dogs.pdf")
ggplot(df2.long, aes(x=passage, y=frequency, group=variant)) +
  geom_point(show.legend=FALSE) +
  geom_point(data=df_highlight, aes(x=passage, y=frequency,color=passage),show.legend=FALSE) +
  facet_wrap(vars(variant))
dev.off()

#pull out dogs in case you want to highligh as a different color
#df_highlight <- df2.long %>%
#filter(passage=="Dogs")

#pull out inoculums to connect with a line/color
df_SARS2 <- df2.long %>%
  filter(passage!="Dogs")

#plot
pdf(file="TC_variants_and_dogs.pdf")
ggplot(df2.long, aes(x=passage, y=frequency, group=variant)) +
  geom_point(show.legend=FALSE) +
  geom_point(data=df_highlight, aes(x=passage,y=frequency), 
             color='red')+
  geom_line(data=df_SARS2)+
  facet_wrap(vars(variant))
dev.off()