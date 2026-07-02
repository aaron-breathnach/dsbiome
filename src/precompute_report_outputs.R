#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(lmerTest)
  library(glmmTMB)
  library(broom.mixed)
  library(compositions)
  library(ape)
})

source("src/import_data.R")

if (!dir.exists("tables/report")) dir.create("tables/report", recursive = TRUE)
if (!dir.exists("plots/report")) dir.create("plots/report", recursive = TRUE)
if (!dir.exists("tables")) dir.create("tables", recursive = TRUE)

pal <- c("Control" = "#55bce4", "Down Syndrome" = "#661f7b")

# -----------------------------------------------------------------------------
# Data import
# -----------------------------------------------------------------------------

data <- import_data()

metadata <- data$metadata %>%
  filter(group %in% c("Control", "Down Syndrome")) %>%
  mutate(
    group = factor(group, levels = c("Control", "Down Syndrome")),
    family_id = factor(family_id),
    sex = factor(sex),
    antibiotics = factor(antibiotics)
  )

mpa <- data$mpa %>%
  filter(sample_id %in% metadata$sample_id) %>%
  arrange(match(sample_id, metadata$sample_id))

humann <- data$humann %>%
  filter(sample_id %in% metadata$sample_id) %>%
  arrange(match(sample_id, metadata$sample_id))

# -----------------------------------------------------------------------------
# Study design tables
# -----------------------------------------------------------------------------

sample_overview <- metadata %>%
  count(group, name = "n_samples")

household_overview <- tibble(
  n_samples = nrow(metadata),
  n_families = n_distinct(metadata$family_id),
  n_species_features = ncol(mpa) - 1,
  n_humann_pathways = ncol(humann) - 1
)

family_size <- metadata %>%
  count(family_id, name = "n_family_members") %>%
  count(n_family_members, name = "n_families")

age_summary <- metadata %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median_age = median(age, na.rm = TRUE),
    iqr_age = IQR(age, na.rm = TRUE),
    min_age = min(age, na.rm = TRUE),
    max_age = max(age, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(sample_overview, "tables/report/sample_overview.tsv")
readr::write_tsv(household_overview, "tables/report/household_overview.tsv")
readr::write_tsv(family_size, "tables/report/family_size.tsv")
readr::write_tsv(age_summary, "tables/report/age_summary.tsv")

# -----------------------------------------------------------------------------
# Alpha diversity
# -----------------------------------------------------------------------------

alpha <- mpa %>%
  column_to_rownames("sample_id") %>%
  vegan::diversity(index = "shannon") %>%
  as.data.frame() %>%
  rename(shannon_index = 1) %>%
  rownames_to_column("sample_id") %>%
  inner_join(metadata, by = "sample_id")

alpha_mod <- lmerTest::lmer(shannon_index ~ group * age + (1 | family_id), data = alpha)

alpha_stats <- anova(alpha_mod) %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  rename(p.value = `Pr(>F)`) %>%
  select(term, NumDF, DenDF, `F value`, p.value)

readr::write_tsv(alpha_stats, "tables/report/alpha_stats.tsv")

p_alpha <- ggplot(alpha, aes(x = group, y = shannon_index)) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_jitter(aes(colour = group, fill = group), pch = 21, width = 0.25, alpha = 0.6, show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  theme_classic(base_size = 12.5) +
  theme(axis.title = element_text(face = "bold")) +
  labs(x = "Group", y = "Shannon diversity")

ggsave("plots/report/alpha_diversity.png", p_alpha, width = 9, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Beta diversity and household distances
# -----------------------------------------------------------------------------

mpa_mat <- mpa %>%
  column_to_rownames("sample_id")

bc <- vegan::vegdist(mpa_mat, method = "bray")

meta_for_vegan <- metadata %>%
  filter(sample_id %in% labels(bc)) %>%
  column_to_rownames("sample_id")

permanova <- vegan::adonis2(
  bc ~ group + age + sex + antibiotics,
  data = meta_for_vegan,
  strata = meta_for_vegan$family_id,
  by = "margin"
) %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  select(term, Df, SumOfSqs, R2, F, `Pr(>F)`)

readr::write_tsv(permanova, "tables/report/permanova.tsv")

pcoa <- ape::pcoa(bc)
eig_val <- round(100 * pcoa$values[1:2, 2], 2)

pcoa_points <- pcoa$vectors[, 1:2] %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  inner_join(metadata, by = "sample_id")

p_pcoa <- ggplot(pcoa_points, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(fill = group), pch = 21, size = 3, alpha = 0.85) +
  scale_fill_manual(values = pal) +
  theme_classic(base_size = 12.5) +
  theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold")) +
  labs(
    x = paste0("PCoA1 [", eig_val[1], "%]"),
    y = paste0("PCoA2 [", eig_val[2], "%]"),
    fill = "Group"
  )

ggsave("plots/report/pcoa_bray_curtis.png", p_pcoa, width = 9, height = 6, dpi = 300)

families <- metadata %>%
  select(sample_id, family_id)

pairwise_dists <- as.matrix(bc) %>%
  as.data.frame() %>%
  rownames_to_column("sample1") %>%
  pivot_longer(-sample1, names_to = "sample2", values_to = "distance") %>%
  filter(sample1 < sample2) %>%
  inner_join(families, by = c("sample1" = "sample_id")) %>%
  inner_join(families, by = c("sample2" = "sample_id"), suffix = c("1", "2")) %>%
  mutate(comparison = if_else(family_id1 == family_id2, "Within family", "Between families"))

household_distance_summary <- pairwise_dists %>%
  group_by(comparison) %>%
  summarise(
    n_pairs = n(),
    median_distance = median(distance),
    mean_distance = mean(distance),
    .groups = "drop"
  )

readr::write_tsv(household_distance_summary, "tables/report/household_distance_summary.tsv")

p_household <- ggplot(pairwise_dists, aes(x = distance, colour = comparison)) +
  geom_density(linewidth = 1) +
  theme_classic(base_size = 12.5) +
  theme(axis.title = element_text(face = "bold"), legend.title = element_blank()) +
  labs(x = "Bray-Curtis distance", y = "Density")

ggsave("plots/report/household_distance_density.png", p_household, width = 9, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Species-level mixed models
# -----------------------------------------------------------------------------

get_top_spp <- function(mpa_in, abundance = 0.01, prevalence = 0.1) {
  mpa_in %>%
    pivot_longer(!sample_id, names_to = "species", values_to = "relab") %>%
    filter(relab > abundance) %>%
    group_by(species) %>%
    tally() %>%
    filter(n > prevalence * nrow(mpa_in)) %>%
    pull(species)
}

get_species_dat <- function(sp, mat, metadata_in, value_name = "relab") {
  mat %>%
    select(sample_id, all_of(sp)) %>%
    rename(!!value_name := 2) %>%
    inner_join(metadata_in, by = "sample_id")
}

top_spp <- get_top_spp(mpa)

clr <- mpa %>%
  column_to_rownames("sample_id") %>%
  compositions::clr() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

run_lmer_species <- function(sp) {
  dat <- get_species_dat(sp, clr, metadata, "relab")
  lmerTest::lmer(relab ~ group * age + (1 | family_id), data = dat) %>%
    anova() %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    rename(p.value = `Pr(>F)`) %>%
    mutate(feature = sp, method = "CLR abundance") %>%
    select(feature, method, term, p.value)
}

run_glmm_species <- function(sp) {
  dat <- get_species_dat(sp, mpa, metadata, "relab") %>%
    mutate(pres_abs = if_else(relab > 0.01, 1, 0))

  glmmTMB::glmmTMB(
    pres_abs ~ group * age + (1 | family_id),
    data = dat,
    family = "binomial"
  ) %>%
    broom.mixed::tidy() %>%
    filter(effect == "fixed", term != "(Intercept)") %>%
    transmute(
      feature = sp,
      method = "Prevalence",
      term = str_replace_all(term, c("groupDown Syndrome" = "group", "groupDown Syndrome:age" = "group:age")),
      p.value = p.value
    )
}

species_stats <- bind_rows(
  purrr::map_dfr(top_spp, purrr::possibly(run_lmer_species, otherwise = tibble())),
  purrr::map_dfr(top_spp, purrr::possibly(run_glmm_species, otherwise = tibble()))
) %>%
  group_by(method, term) %>%
  mutate(fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(fdr)

species_summary <- species_stats %>%
  group_by(method, term) %>%
  summarise(
    n_features = n(),
    min_fdr = min(fdr, na.rm = TRUE),
    n_fdr_0_05 = sum(fdr < 0.05, na.rm = TRUE),
    n_fdr_0_10 = sum(fdr < 0.10, na.rm = TRUE),
    n_fdr_0_25 = sum(fdr < 0.25, na.rm = TRUE),
    .groups = "drop"
  )

top_species_group <- species_stats %>%
  filter(term == "group") %>%
  arrange(fdr) %>%
  group_by(method) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  select(method, feature, p.value, fdr)

readr::write_tsv(species_stats, "tables/report/species_stats.tsv")
readr::write_tsv(species_summary, "tables/report/species_summary.tsv")
readr::write_tsv(top_species_group, "tables/report/top_species_group.tsv")

# -----------------------------------------------------------------------------
# HUMAnN mixed models
# -----------------------------------------------------------------------------

run_lmer_humann <- function(pwy, humann_in, metadata_in) {
  dat <- humann_in %>%
    select(sample_id, all_of(pwy)) %>%
    rename(cpm = 2) %>%
    mutate(log_cpm = log1p(cpm)) %>%
    inner_join(metadata_in, by = "sample_id")

  lmerTest::lmer(log_cpm ~ group * age + (1 | family_id), data = dat) %>%
    anova() %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    rename(p.value = `Pr(>F)`) %>%
    mutate(feature = pwy) %>%
    select(feature, term, p.value)
}

humann_stats <- purrr::map_dfr(colnames(humann)[-1], function(pwy) {
  run_lmer_humann(pwy, humann, metadata)
}) %>%
  group_by(term) %>%
  mutate(fdr = p.adjust(p.value, "BH")) %>%
  ungroup() %>%
  arrange(p.value)

humann_summary <- humann_stats %>%
  group_by(term) %>%
  summarise(
    n_pathways = n(),
    min_p = min(p.value, na.rm = TRUE),
    min_fdr = min(fdr, na.rm = TRUE),
    n_fdr_0_05 = sum(fdr < 0.05, na.rm = TRUE),
    n_fdr_0_10 = sum(fdr < 0.10, na.rm = TRUE),
    n_fdr_0_25 = sum(fdr < 0.25, na.rm = TRUE),
    .groups = "drop"
  )

humann_top <- humann_stats %>%
  arrange(fdr) %>%
  group_by(term) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  select(term, feature, p.value, fdr)

readr::write_tsv(humann_stats, "tables/stats_humann.tsv")
readr::write_tsv(humann_summary, "tables/report/humann_summary.tsv")
readr::write_tsv(humann_top, "tables/report/humann_top.tsv")

# -----------------------------------------------------------------------------
# Strain-sharing summary tables for report rendering
# -----------------------------------------------------------------------------

strain_summary <- readr::read_tsv("figures/strain_sharing/strain_sharing_summary.tsv", show_col_types = FALSE)
shared_by_species <- readr::read_tsv("figures/strain_sharing/shared_strains_by_species.tsv", show_col_types = FALSE)

strain_summary_report <- strain_summary %>%
  mutate(proportion_shared = scales::percent(proportion_shared, accuracy = 0.1))

readr::write_tsv(strain_summary_report, "tables/report/strain_sharing_summary.tsv")
readr::write_tsv(shared_by_species %>% arrange(desc(n_shared_strains)), "tables/report/strain_sharing_by_species.tsv")

message("Report outputs precomputed successfully.")
