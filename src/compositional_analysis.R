library(tidyverse)

get_alpha_div <- function(mat, metadata) {
  
  mpa %>%
    column_to_rownames("sample_id") %>%
    vegan::diversity() %>%
    as.data.frame() %>%
    rename(shannon_index = 1) %>%
    rownames_to_column("sample_id") %>%
    inner_join(metadata, by = "sample_id") %>%
    filter(group %in% c("Control", "Down Syndrome")) %>%
    mutate(group = str_to_sentence(group)) %>%
    as_tibble()
  
}

make_alpha_div_plot <- function(mpa, metadata) {
  
  alpha <- get_alpha_div(mpa, metadata)
  
  title <- lmerTest::lmer(shannon_index ~ group * age + (1 | family_id), data = alpha) %>%
    anova() %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term == "group") %>%
    pull(7) %>%
    round(2) %>%
    sprintf("*p*=%s", .)
  
  ggplot(alpha, aes(x = group, y = shannon_index)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
      aes(colour = group, fill = group),
      show.legend = FALSE,
      pch = 21,
      width = 0.25,
      alpha = 0.5
    ) +
    scale_colour_manual(values = get_pal("group")) +
    scale_fill_manual(values = get_pal("group")) +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      plot.title = ggtext::element_markdown()
    ) +
    ggtitle(title) +
    labs(x = "Group", y = "Shannon index")
  
}

viz_alpha_div_corr <- function(mpa, metadata) {
  
  dat <- get_alpha_div(mpa, metadata) %>%
    select(family_id, group, shannon_index) %>%
    pivot_wider(names_from = "group", values_from = "shannon_index") %>%
    rename(down_syndrome = 2, control = 3) %>%
    drop_na()
  
  cor <- rstatix::cor_test(dat[-1], method = "spearman")
  rval <- cor$cor
  pval <- cor$p
  title <- sprintf("\U03C1=%s; *p*=%s", rval, pval)
  
  ggplot(dat, aes(x = down_syndrome, y = control)) +
    geom_smooth(method = "lm", colour = "#BD3027") +
    geom_point(pch = 21, fill = "#90D4CC") +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      plot.title = ggtext::element_markdown()
    ) +
    ggtitle(title) +
    labs(
      x = "Shannon index of Down syndrome family member",
      y = "Shannon index of control family member"
    )
  
}

make_pcoa_plot <- function(veg_dis, metadata) {
  
  dat <- metadata %>%
    filter(sample_id %in% labels(veg_dis)) %>%
    drop_na() %>%
    column_to_rownames("sample_id")
  
  permanova <- vegan::adonis2(
    veg_dis ~ group + age + sex + antibiotics,
    dat,
    strata = dat$family_id,
    by = "margin"
  ) %>%
    as.data.frame() %>%
    rownames_to_column("covariate") %>%
    select(covariate, 4, 6) %>%
    drop_na() %>%
    arrange(R2) %>%
    mutate(covariate = factor(covariate, levels = .$covariate))

  title <- sprintf(
    "PERMANOVA: R<sup>2</sup>=%s; *p*=%s",
    round(filter(permanova, covariate == "group")[[1, 2]], 2),
    round(filter(permanova, covariate == "group")[[1, 3]], 2)
  )
  
  pcoa <- ape::pcoa(veg_dis)
  
  points <- pcoa$vectors[,1:2] %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    inner_join(metadata, by = "sample_id") %>%
    as_tibble() %>%
    filter(group %in% c("Control", "Down Syndrome")) %>%
    mutate(group = str_to_sentence(group))
  
  eig_val <- round(100 * pcoa$values[1:2, 2], 2)
  xlab <- paste0("PCoA1 [", eig_val[1], "%]")
  ylab <- paste0("PCoA2 [", eig_val[2], "%]")
  
  ggplot(points, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = group), pch = 21, size = 3) +
    scale_fill_manual(values = get_pal("group")) +
    ggtitle(title) +
    labs(fill = "Group", x = xlab, y = ylab) +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.title = ggtext::element_markdown()
    )
  
}

viz_bc_dists <- function(veg_dis, metadata) {
  
  families <- metadata %>%
    select(sample_id, family_id)
  
  dists <- as.matrix(veg_dis) %>%
    as.data.frame() %>%
    rownames_to_column("sample1") %>%
    pivot_longer(-sample1, names_to = "sample2", values_to = "dist") %>%
    filter(sample1 < sample2) %>%
    inner_join(families, by = c("sample1" = "sample_id")) %>%
    inner_join(families, by = c("sample2" = "sample_id")) %>%
    rename(family1 = 4, family2 = 5)
  
  intra_dists <- dists %>%
    filter(sample1 != sample2 & family1 == family2) %>%
    mutate(intra_inter = "Intra-family")
  
  inter_dists <- dists %>%
    filter(sample1 != sample2 & family1 != family2) %>%
    mutate(intra_inter = "Inter-family")
  
  df <- rbind(intra_dists, inter_dists)
  
  medians <- df %>%
    group_by(intra_inter) %>%
    summarise(median = median(dist))
  
  ggplot(df, aes(x = dist)) +
    geom_vline(
      xintercept = medians$median,
      colour = wesanderson::wes_palettes$Cavalcanti1[1:2],
      linetype = "dashed") +
    geom_density(aes(colour = intra_inter), linewidth = 1) +
    labs(x = "Bray-Curtis distance", y = "Density", colour = "") +
    theme_classic(base_size = 12.5) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.position = c(0.2, 0.9)
    ) +
    scale_colour_manual(values = wesanderson::wes_palettes$Cavalcanti1)
  
}

run_diversity_analysis <- function() {
  
  data <- import_data()
  mpa <- data$mpa
  metadata <- data$metadata
  
  veg_dis <- mpa %>%
    column_to_rownames("sample_id") %>%
    vegan::vegdist()
  
  p1 <- make_alpha_div_plot(mpa, metadata)
  p2 <- viz_alpha_div_corr(mpa, metadata)
  p3 <- make_pcoa_plot(veg_dis, metadata)
  p4 <- viz_bc_dists(veg_dis, metadata)
  
  plot_list <- list(p2, p1, p3, p4)
  
  p <- patchwork::wrap_plots(plot_list)
  
  ggsave("plots/diversity.png", p, height = 10, width = 12.5)
  
}
