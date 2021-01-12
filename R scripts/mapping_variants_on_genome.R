###mapping variants along genome in ggplot

#load packages
library(tidyverse)
library(ggplot2)

#start with data frame of cleaned data, combined replicates
df <- as.data.frame(variant_table_cleaned_combined)

#make data long and skinny 
df.long <- df %>%
  pivot_longer(!c(variant,position,cds,N_or_S), names_to= "dataset", values_to ="frequency")

#make datasets a factor, then reorder
factor(df.long$dataset, levels=c("SARS2_P1","SARS2_P2","SARS2_P3","C43","C56","C57","C58","C69","C70", "D46", "D47","D48", "F2678","H63","H65","H68"))

#plot first try
pdf(file="mapping_variants.pdf")
ggplot(df.long %>% filter(!is.na(frequency))) +
  geom_point(aes(y=dataset, x=position, fill=frequency), shape=21, size=2.5, alpha=0.8, color="black", stroke=0.2)  +
  scale_fill_gradient(low="#FFDFDF", high="#FF0000", breaks=c(0,1)) +
  scale_x_continuous(limits=c(0,30000)) +
  theme_classic() +
  theme(text = element_text(size=10, family="sans")) + 
  xlab ("Position in genome (nt)") + 
  ylab ("Dataset") 
dev.off()

#plot as line segments instead, doesn't work yet
# ggplot(df.long %>% filter(!is.na(frequency))) +
#   geom_segment(aes(x=position, xend=position, y=dataset, yend=dataset+1.5, color=frequency), size=1.5) +
#   scale_color_gradient(low="#FFDFDF", high="#FF0000", breaks=c(0,1)) +
#   scale_x_continuous(limits=c(0,30000)) +
#   theme_classic() +
#   theme(text = element_text(size=10, family="sans")) + 
#   xlab ("Position in genome (nt)") + 
#   ylab ("Dataset") 
