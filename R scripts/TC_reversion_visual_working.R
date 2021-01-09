##visualizing TC variant emergence and reversion

#load libraries
library(tidyverse)
library(data.table)

#load in cleaned data with replicates combined
df <- as.data.frame(variant_table_cleaned_combined)

#pull out just inoculums and dogs for now
df2 <- df[, c("variant", "cds", "D46", "D47", "D48", "SARS2_P1","SARS2_P2","SARS2_P3")]
#just inoculums df2 <- df[, c("variant", "cds", "SARS2_P1","SARS2_P2","SARS2_P3")]

#rename passages as 1, 2, 3
names(df2) <- gsub("SARS2_P", "", names(df2))

#pull out just TC mutations that went to fixation
df2 <- df2 %>%
  filter(variant %in% c("D135E", "D215H", "R685H", "T7I","S194T"))

#make all dogs into one passage
names(df2) <- gsub("D46|D47|D48","Dogs",names(df2))

#make data long and skinny for plotting, rename value as frequency
df2.long <- melt(setDT(df2), id.vars=c("variant", "cds"), variable.name="passage")
names(df2.long)[names(df2.long) == 'value'] <- 'frequency'

#NAs to 0s for plotting
df2.long[is.na(df2.long)] <- 0

#get the order right for passages
df2.long <- df2.long %>%mutate(passage = fct_relevel(passage, "1","2","3", "Dogs"))

#plot
ggplot(df2.long, aes(x=passage, y=frequency, group=variant, color= variant)) +
  geom_line(show.legend=FALSE) +
  geom_point(show.legend=FALSE) +
  facet_wrap(vars(variant))
