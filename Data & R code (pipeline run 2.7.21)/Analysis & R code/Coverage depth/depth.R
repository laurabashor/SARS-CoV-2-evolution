##testing depth script edits

#load libraries
library(tidyverse)
library(openxlsx)
library(pdftools)

#load data, rename columns
depth_df <- read.delim("all.depth", sep="\t", header=FALSE)
colnames(depth_df) <- c("dataset", "virus", "position", "depth")

# rename datasets 
#source(paste0(r_bindir, "/process_dataset_names.R"))
depth_df <- process_dataset_names(depth_df)

# higlight coverage below a certain limit
# TODO: Parameterize this
min_depth_highlight <- 100
depth_df <- depth_df %>% mutate(above_highlight = if_else(depth > min_depth_highlight, TRUE, FALSE))

# max position for drawing a red rectangle showing low coverage
max_position = max(depth_df$position)

######### calculate average depth in windows
# %/% is the integer division operator
window_size = 10
depth_df <- depth_df %>% mutate (window = position %/% window_size)

# calculate average coverage depth in each window
df_windowed <- depth_df %>% 
  group_by(dataset, virus, window)  %>% 
  summarize(depth = mean(depth), .groups = "drop") %>% 
  mutate(position = (window*window_size) + 1) %>% 
  ungroup()

##now plot coverage data on multiple pdf pages

datasets <- depth_df %>% group_by(dataset) %>% summarize(.groups="drop") %>% pull()

plot_some_datasets <- function(datasets){

  datasets_per_page <- 12

  page_number <- 1

  # iterate through the viruses, doing so many per page
  for (i in seq(1, length(datasets), datasets_per_page)) {

    # plots_per_page at a time
    subset_datasets <- datasets[i:(i+(datasets_per_page-1))]

    # output to console which ones we're doing
    print(paste0(subset_datasets))

    pdf_name <- paste0("coverage_plot_page_", page_number, ".pdf")

    # generate & print plot
    plot_datasets(subset_datasets, pdf_name)

    page_number = page_number + 1
    
  }

}

plot_datasets <- function(dataset_names, pdf_name){
  
  # subset the main dataframes to get the data just for these viruses
  subset_df <- df_windowed %>% filter(dataset %in% dataset_names)
  
  # output the virus names to the console
  print(paste0(dataset_names))
  
  # p returned here is a ggplot object
  p <- ggplot(subset_df) + 
    geom_line(aes(x=position, y=depth), size=0.5) +
    geom_area(aes(x=position, y=depth), fill="lightgrey") +
    annotate("rect", xmin=0, xmax=max_position, ymin=1, ymax=min_depth_highlight, alpha=0.05, fill="red") +
    scale_color_manual(values = c("red", "black")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_y_log10() +
    xlab("genome position (nt)") +
    ylab (paste0("coverage depth\n", "coverage below ", min_depth_highlight, " in red"))+
    facet_grid(dataset~virus, scales="free", space="free_x")
  
  ggsave(pdf_name, p, height=10.5, width=7.5, units="in")

    # print p will make the plots appear on viewer
    # print(p)

  }

plot_some_datasets(datasets)

pdf_combine(c(list.files(pattern="pdf$")), output="coverage_plot.pdf")

# calculated median depth of (total) coverage for each virus in each dataset and store it in a new df
median_depths <- depth_df %>% 
  group_by(dataset, virus) %>% 
  summarize(median_depth = median(depth))

median_depth_all_together <- depth_df %>%
  summarize(median_depth = median(depth))

mean_depths <- depth_df %>% 
  group_by(dataset, virus) %>% 
  summarize(mean_depth = mean(depth))

mean_depth_all_together <- depth_df %>%
  summarize(mean_depth = mean(depth))

wb <- createWorkbook("Median_depths.xlsx")
addWorksheet(wb, "median_depth")
addWorksheet(wb, "median_depth_all_together")
addWorksheet(wb, "mean_depth")
addWorksheet(wb, "mean_depth_all_together")
writeData(wb, "median_depth", median_depths)
writeData(wb, "median_depth_all_together", median_depth_all_together)
writeData(wb, "mean_depth", mean_depths)
writeData(wb, "mean_depth_all_together", mean_depth_all_together)
saveWorkbook(wb, "Median_depths.xlsx", overwrite = TRUE)
