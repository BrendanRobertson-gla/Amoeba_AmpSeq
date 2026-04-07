suppressPackageStartupMessages({
  library(dplyr) 
  library(tidyr) 
  library(readr)
  library(readxl)
  library(ggplot2) 
  library(stringr)
})

# ---------- paths ----------
base_dir  <- "/Users/brendan/Documents/Documents - Brendan’s MacBook Pro/University of Glasgow/PhD/Geographical Study/Write up/VCFs/VCFs_forbrendan/Amoeba_GLST"
meta_xlsx <- file.path(base_dir, "MetaData.xlsx")
counts_f  <- file.path(base_dir, "target_counts.tsv")
plot_dir  <- file.path(base_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- options ----------
hit_threshold <- 50

# ---------- helpers ----------
norm_id <- function(x) str_trim(as.character(x))
base_from_sample <- function(x) sub("_.*$", "", norm_id(x))

# ---------- load ----------
stopifnot(file.exists(meta_xlsx), file.exists(counts_f))

metadata <- read_excel(meta_xlsx) %>%
  rename_with(~ str_replace_all(., "\\s+", "_")) %>%
  mutate(across(everything(), as.character)) %>%
  filter(!is.na(Sequencing_ID) & str_trim(Sequencing_ID) != "") %>%
  filter(is.na(Country) | Country != "Ireland")   # remove Ireland

counts <- read_tsv(counts_f, show_col_types = FALSE, col_types = "ccc") %>%
  mutate(
    Sample     = norm_id(Sample),
    SampleBase = base_from_sample(Sample),
    Target     = norm_id(Target),
    Reads      = as.numeric(Reads)
  )

# ---------- metadata join ----------
metadata_keys <- metadata %>%
  transmute(
    Sequencing_ID = norm_id(Sequencing_ID),
    Country       = if ("Country" %in% names(.)) Country else NA_character_,
    Site          = if ("Site" %in% names(.)) Site else NA_character_
  ) %>%
  mutate(
    SID_underscore = str_replace_all(Sequencing_ID, "-", "_"),
    SID_hyphen     = str_replace_all(Sequencing_ID, "_", "-")
  ) %>%
  pivot_longer(c(Sequencing_ID, SID_underscore, SID_hyphen),
               values_to = "KeyID") %>%
  distinct(KeyID, .keep_all = TRUE)

counts <- counts %>%
  left_join(metadata_keys, by = c("SampleBase" = "KeyID"))

# ---------- clean country + site ----------
counts <- counts %>%
  mutate(
    Country = ifelse(is.na(Country) | Country == "", "Unknown", Country),
    Site    = ifelse(is.na(Site) | Site == "", SampleBase, Site)
  ) %>%
  filter(Country != "Ireland")

# ---------- totals ----------
sample_totals <- counts %>%
  group_by(SampleBase, Country, Site) %>%
  summarise(TotalReads = sum(Reads, na.rm = TRUE), .groups = "drop")

site_labels <- sample_totals %>%
  group_by(Site) %>%
  mutate(SiteLabel = ifelse(n() == 1, Site, paste0(Site, " (", row_number(), ")"))) %>%
  ungroup() %>%
  mutate(
    SiteLabel = make.unique(SiteLabel),
    Total_log10 = log10(TotalReads + 1)
  )

site_labels$Country <- factor(site_labels$Country,
                              levels = c("Scotland","Norway","Unknown"))

sample_order <- site_labels %>%
  arrange(Country, desc(TotalReads)) %>%
  distinct(SiteLabel) %>%
  pull(SiteLabel)

# ---------- target ordering ----------
hit_stats <- counts %>%
  group_by(Target) %>%
  summarise(
    HitRate = mean(Reads > hit_threshold, na.rm = TRUE),
    TargetTotal = sum(Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(HitRate), desc(TargetTotal))

target_order <- unique(hit_stats$Target)

# ---------- heatmap data ----------
tmp <- counts %>%
  left_join(
    site_labels %>%
      select(SampleBase, SiteLabel, Country_site = Country, TotalReads, Total_log10),
    by = "SampleBase"
  ) %>%
  mutate(Country = coalesce(Country_site, Country)) %>%
  transmute(SiteLabel, Target, Reads, TotalReads, Total_log10, Country)

country_lookup <- tmp %>%
  distinct(SiteLabel, Country)

heat_dat <- tmp %>%
  select(-Country) %>%
  complete(SiteLabel = sample_order, Target = target_order,
           fill = list(Reads = 0, TotalReads = 0, Total_log10 = 0)) %>%
  left_join(country_lookup, by = "SiteLabel") %>%
  mutate(
    SiteLabel = factor(SiteLabel, levels = sample_order),
    Target    = factor(Target, levels = target_order),
    logReads  = log10(Reads + 1),
    Hit       = Reads > hit_threshold
  )

# ---------- HEATMAP ONLY ----------
p_heat <- ggplot(heat_dat, aes(x = Target, y = SiteLabel, fill = logReads)) +
  geom_tile() +
  geom_tile(data = subset(heat_dat, Hit),
            color = "white", linewidth = 0.15, fill = NA) +
  scale_fill_viridis_c(name = "log10(reads+1)", option = "C") +
  labs(
    x = "Target",
    y = "Site",
    title = "Targets hit & reads per target"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
    panel.grid  = element_blank(),
    plot.margin = margin(6, 6, 6, 6)
  )

print(p_heat)

#ggsave(file.path(plot_dir, "heatmap_targets_reads.png"), p_heat, width = 14, height = 12, dpi = 300)
#ggsave(file.path(plot_dir, "heatmap_targets_reads.pdf"), p_heat, width = 14, height = 12)




# =========================================================
# Custering
# =========================================================

suppressPackageStartupMessages({
  library(tibble)
})

# ---------- parameters ----------
min_prevalence <- 0.05   # keep targets present in ≥5% of sites
agg_fun <- "sum"         # "sum" or "mean"

# ---------- collapse to site level ----------
site_level <- heat_dat %>%
  group_by(SiteLabel, Target, Country) %>%
  summarise(
    Reads = if (agg_fun == "sum") sum(Reads, na.rm = TRUE) else mean(Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(logReads = log10(Reads + 1))

# ---------- filter low-prevalence targets ----------
target_prev <- site_level %>%
  group_by(Target) %>%
  summarise(prev = mean(Reads > 0), .groups = "drop")

keep_targets <- target_prev %>%
  filter(prev >= min_prevalence) %>%
  pull(Target)

site_level <- site_level %>%
  filter(Target %in% keep_targets)

# ---------- create matrix for clustering ----------
mat <- site_level %>%
  select(SiteLabel, Target, logReads) %>%
  pivot_wider(names_from = Target, values_from = logReads, values_fill = 0) %>%
  column_to_rownames("SiteLabel") %>%
  as.matrix()

# ---------- hierarchical clustering ----------
row_order <- rownames(mat)[hclust(dist(mat))$order]
col_order <- colnames(mat)[hclust(dist(t(mat)))$order]

# ---------- apply ordering ----------
site_level <- site_level %>%
  mutate(
    SiteLabel = factor(SiteLabel, levels = row_order),
    Target    = factor(Target, levels = col_order)
  )

# ---------- clustered heatmap ----------
p_clustered <- ggplot(site_level, aes(x = Target, y = SiteLabel, fill = logReads)) +
  geom_tile() +
  scale_fill_viridis_c(name = "log10(reads+1)", option = "C") +
  labs(
    x = "Target",
    y = "Site",
    title = "Clustered heatmap (site-level, filtered targets)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
    panel.grid  = element_blank(),
    plot.margin = margin(6, 6, 6, 6)
  )

# ---------- facet by country ----------
p_clustered_facet <- p_clustered +
  facet_grid(Country ~ ., scales = "free_y", space = "free_y")

print(p_clustered_facet)

# ---------- save outputs ----------
#ggsave(file.path(plot_dir, "heatmap_clustered.png"),
 #      p_clustered, width = 14, height = 10, dpi = 300)

#ggsave(file.path(plot_dir, "heatmap_clustered.pdf"),
  #     p_clustered, width = 14, height = 10)

#ggsave(file.path(plot_dir, "heatmap_clustered_by_country.png"),
   #    p_clustered_facet, width = 14, height = 12, dpi = 300)

#ggsave(file.path(plot_dir, "heatmap_clustered_by_country.pdf"),
 #      p_clustered_facet, width = 14, height = 12)

print(p_clustered)
print(p_clustered_facet)


