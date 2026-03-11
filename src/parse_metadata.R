library(tidyverse)

xlsx <- "info_from_therese/Master_Health_Survey_Nov 26th_2025_TM.xlsx"

cols <- c(
  "sample_id",
  "family_id",
  "subject_id",
  "group",
  "sex",
  "age",
  "antibiotics"
)

dat <- readxl::read_excel(xlsx, skip = 2) %>%
  select(1:3, 5:6, 9, 33) %>%
  setNames(cols) %>%
  arrange(age)

dat$age <- dat$age %>%
  lapply(function(x) if (grepl("Months", x)) {
    as.numeric(gsub(" .*", "", x)) / 12
  } else {
    as.numeric(x)
  }) %>%
  unlist()

d230 <- readxl::read_excel("info_from_therese/Master_FFQ_Nov 26th_2025_TM.xlsx", skip = 2) %>%
  select(1:3, 5) %>%
  setNames(cols[1:4]) %>%
  filter(!sample_id %in% dat$sample_id) %>%
  mutate(sex = NA, age = NA, antibiotics = NA)

dat <- rbind(dat, d230) %>%
  mutate(sample_id = sample_id %>%
           str_replace("D", "") %>%
           str_pad(3, "left", "0") %>%
           str_c("D", .))

write_tsv(dat, "data/metadata.tsv")
