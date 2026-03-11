import_data <- function() {
  
  metadata <- read_delim("data/metadata.tsv") %>%
    drop_na() %>%
    arrange(sample_id)
  
  ## two samples failed on the first sequencing run
  drop <- c("Dsbiome_D61_output", "Dsbiome_138_output")
  
  mpa <- read_delim("data/merged_metaphlan_tables.txt", skip = 1) %>%
    rename(species = 1) %>%
    filter(grepl("s__", species) & !grepl("t__", species)) %>%
    mutate(species = gsub(".*\\|", "", species)) %>%
    select(!contains("Neg"), -Dsbiome_D61_output, -Dsbiome_138_output) %>%
    pivot_longer(!species, names_to = "sample_id", values_to = "relab") %>%
    pivot_wider(names_from = "species", values_from = "relab")
  
  mpa$sample_id <- mpa$sample_id %>%
    str_replace("Dsbiome_", "") %>%
    str_replace("DS_BIOME_", "") %>%
    str_replace("_output", "") %>%
    str_replace("D", "") %>%
    str_pad(3, "left", "0") %>%
    str_c("D", .)
  
  mpa <- arrange(mpa, sample_id) %>%
    filter(sample_id %in% metadata$sample_id) %>%
    arrange(metadata)
  
  metadata <- arrange(metadata, sample_id)
  
  list("mpa" = mpa, "metadata" = metadata)
  
}
