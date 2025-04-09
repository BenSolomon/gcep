require(readxl)
require(dplyr)
require(tidyr)
require(readr)
require(here)

df_iuis <- readxl::read_xlsx(here("data/raw_data/IUIS-IEI-list-for-web-site-July-2024V2.xlsx")) %>% 
  select(gene = `Genetic defect`, group = `Major category`, subgroup = Subcategory) %>% 
  separate(group, into = c(NA,"table","table_description"), extra = "merge", sep = " ") %>% 
  separate(subgroup, into = c(NA,"subtable","subtable_description"), extra = "merge", sep = " ") %>% 
  mutate(subtable = ifelse(is.na(subtable), 0, subtable)) %>% 
  mutate(subtable_description = ifelse(is.na(subtable_description), "None", subtable_description)) %>% 
  mutate(gene = gsub(" .*","",gene))

write_csv(df_iuis, here("data/iuis_table.csv"))