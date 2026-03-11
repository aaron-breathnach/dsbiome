get_top_spp <- function(mpa, abundance = 0.01, prevalence = 0.1) {
  mpa %>%
    pivot_longer(!sample_id, names_to = "species", values_to = "relab") %>%
    filter(relab > abundance) %>%
    group_by(species) %>%
    tally() %>%
    filter(n > prevalence * nrow(mpa)) %>%
    pull(species)
}

get_dat <- function(sp, mpa, metadata) {
  
  mpa %>%
    select(sample_id, all_of(sp)) %>%
    rename(relab = 2) %>%
    inner_join(metadata, by = "sample_id") %>%
    filter(group %in% c("Control", "Down Syndrome") & !is.na(family_id))
  
}

run_lmer <- function(sp, clr, metadata) {
  
  dat <- get_dat(sp, clr, metadata)
  
  lmerTest::lmer(relab ~ group * age + (1 | family_id), data = dat) %>%
    anova() %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    rename(p.value = 7) %>%
    mutate(species = sp) %>%
    select(species, term, p.value)
  
}

run_glmm <- function(sp, mpa, metadata) {
  
  dat <- get_dat(sp, mpa, metadata) %>%
    mutate(pres_abs = ifelse(relab > 0.01, 1, 0))
  
  pval <- glmmTMB::glmmTMB(
    pres_abs ~ group * age + (1 | family_id),
    data = dat,
    family = "binomial"
  ) %>%
    broom.mixed::tidy() %>%
    select(4, 8) %>%
    filter(!grepl("\\(", term)) %>%
    mutate(term = gsub("Down Syndrome", "", term)) %>%
    mutate(species = sp) %>%
    select(species, term, p.value)
  
}

## visualisation

plot_res <- function(sp, mpa, metadata, res, pal) {
  
  print(sp)
  
  tmp <- res %>%
    filter(species == sp)
  
  lmm  <- tmp[[1, 2]]
  glmm <- tmp[[1, 3]]
  
  subtitle1 <- sprintf(
    "LMM *q*-value=%s",
    ifelse(lmm < 0.001, scales::scientific(lmm), round(lmm, 3))
  )
  
  subtitle2 <- sprintf(
    "GLMM *q*-value=%s",
    ifelse(glmm < 0.001, scales::scientific(glmm), round(glmm, 3))
  )
  
  df1 <- get_dat(sp, mpa, metadata)
  
  p1 <- ggplot(df1, aes(x = group, y = relab)) +
    geom_jitter(
      aes(colour = group, fill = group),
      pch = 21,
      alpha = 0.5,
      width = 0.25,
      show.legend = FALSE
    ) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    scale_y_continuous(trans = "log1p") +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    ggtitle(sp) +
    labs(x = "Group", y = "log(relative abundance [%])", subtitle = subtitle1) +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      plot.subtitle = ggtext::element_markdown()
    )
  
  df2 <- df1 %>%
    mutate(pres_abs = ifelse(relab > 0, 1, 0)) %>%
    group_by(group) %>%
    summarise(perc = 100 * sum(pres_abs) / n())
  
  p2 <- ggplot(df2, aes(x = group, y = perc)) +
    geom_bar(
      stat = "identity",
      aes(fill = group), colour = "black",
      width = 0.5,
      show.legend = FALSE
    ) +
    geom_text(
      aes(label = paste0(round(perc, 2), "%")),
      vjust = -0.5
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = get_pal("group")) +
    labs(x = "Group", y = "Prevalence [%]", subtitle = subtitle2) +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      plot.subtitle = ggtext::element_markdown()
    )
  
  patchwork::wrap_plots(p1, p2)
  
}

plot_lmm <- function(sp, mpa, metadata, res, pal) {
  
  print(sp)
  
  tmp <- res %>%
    filter(species == sp)
  
  lmm  <- tmp[[1, 2]]
  
  subtitle1 <- sprintf(
    "LMM *q*-value=%s",
    ifelse(lmm < 0.001, scales::scientific(lmm), round(lmm, 3))
  )
  
  df1 <- get_dat(sp, mpa, metadata)
  
  ggplot(df1, aes(x = group, y = relab)) +
    geom_jitter(
      aes(colour = group, fill = group),
      pch = 21,
      alpha = 0.5,
      width = 0.25,
      show.legend = FALSE
    ) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    scale_y_continuous(trans = "log1p") +
    scale_colour_manual(values = get_pal("group")) +
    scale_fill_manual(values = get_pal("group")) +
    ggtitle(sp) +
    labs(x = "Group", y = "log(relative abundance [%])", subtitle = subtitle1) +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      plot.subtitle = ggtext::element_markdown()
    )
  
}

plot_glmm <- function(sp, mpa, metadata, res, pal) {
  
  dat <- get_dat(sp, mpa, metadata) %>%
    mutate(pres_abs = ifelse(relab > 0, 1, 0))
  
  mod <- glmmTMB::glmmTMB(
    pres_abs ~ group + (1|family_id),
    data = dat,
    family = "binomial"
  )
  
  pred <- ggeffects::ggpredict(
    mod,
    terms = "group",
    bias_correction = TRUE
  ) %>%
    as_tibble() %>%
    select(1:5) %>%
    rename(group = 1)
  
  p_value <- res %>%
    filter(species == sp) %>%
    pull(GLMM)
  
  subtitle <- sprintf(
    "*p*=%s",
    ifelse(p_value < 0.001, scales::scientific(p_value), round(p_value, 3))
  )
  
  ggplot(pred, aes(x = group, y = predicted)) +
    geom_segment(
      aes(
        x = group,
        xend = group,
        y = conf.low,
        yend = conf.high,
        colour = group
      ),
      position = position_nudge(x = c(-0.1, 0.1), y = 0),
      size = 1,
      show.legend = FALSE
    ) +
    geom_point(
      aes(colour = group),
      position = position_nudge(x = c(-0.1, 0.1), y = 0),
      size = 5,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = get_pal("group")) +
    labs(
      title = sp,
      subtitle = subtitle,
      x = "Group",
      y = "Prevalence [%]"
    ) +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.subtitle = ggtext::element_markdown()
    )
  
}

######################
## run the analysis ##
######################

data <- import_data()
mpa <- data$mpa
metadata <- data$metadata

top_spp <- get_top_spp(mpa)

clr <- mpa %>%
  column_to_rownames("sample_id") %>%
  compositions::clr() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

res_lmm <- bind_rows(lapply(top_spp, function(x) run_lmer(x, clr, metadata))) %>%
  mutate(method = "LMM") %>%
  select(method, species, term, p.value)

res_glmm <- bind_rows(lapply(top_spp, function(x) run_glmm(x, mpa, metadata))) %>%
  mutate(method = "GLMM") %>%
  select(method, species, term, p.value)

res <- bind_rows(res_lmm, res_glmm) %>%
  group_by(method, term) %>%
  mutate(fdr = p.adjust(p, "BH")) %>%
  arrange(species) %>%
  select(species, method, fdr) %>%
  pivot_wider(names_from = "method", values_from = "fdr")

sig_spp <- res %>%
  filter(LMM <= 0.25 | GLMM <= 0.25) %>%
  pull(species)

spp_lmm <- res %>%
  filter(LMM <= 0.25) %>%
  pull(species)

p_lmm <- lapply(spp_lmm, function(x) plot_lmm(x, mpa, metadata, res, pal)) %>%
  patchwork::wrap_plots()

spp_glmm <- res %>%
  filter(GLMM <= 0.25) %>%
  pull(species)

p_glmm <- lapply(spp_glmm, function(x) plot_glmm(x, mpa, metadata, res, pal)) %>%
  patchwork::wrap_plots()

filenames <- sprintf(
  "plots/differentially_%s_species.png",
  c("abundant", "prevalent")
)

ggsave(filenames[1], p_lmm, width = 10, height = 5)
ggsave(filenames[2], p_glmm, width = 7.5, height = 7.5)
