#R Code for North Atlantic Map
# Sample counts in labels, no size scaling, corrected sample types

# Load required packages
library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(readxl)
library(ggrepel)

# Set working directory and file paths
base_dir  <- "/Users/brendan/Documents/Documents - Brendan's MacBook Pro/University of Glasgow/PhD/Geographical Study/Write up/VCFs/VCFs_forbrendan"


# Read the metadata file
metadata <- read_excel("MetaData.xlsx")

# Clean and prepare the data with corrected sample types
map_data <- metadata %>%
  rename(Latitude = Lattitude) %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  # Group by site to get sample counts and types
  group_by(Site, Country, Latitude, Longitude) %>%
  summarise(
    Sample_Types = paste(unique(Sample_Type), collapse = " & "),
    Sample_Count = n(),
    .groups = "drop"
  ) %>%
  # Create corrected sample type categories
  mutate(
    Sample_Type_Category = case_when(
      Site %in% c("Fosså", "Salvågvika") ~ "Gill Tissue",
      TRUE ~ "Gill Swab"  # All other sites are swab only
    ),
    # Create site labels with sample counts
    Site_Label = paste0(Site, " (n = ", Sample_Count, ")")
  )

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define map boundaries (back to wider view for context)
lon_limits <- c(-12, 15)  
lat_limits <- c(52, 66)

# Create the final map
p <- ggplot() +
  # Add country polygons
  geom_sf(data = world, 
          fill = "lightgray", 
          color = "white", 
          size = 0.3) +
  
  # Add white outline around each point for contrast
  geom_point(data = map_data,
             aes(x = Longitude, y = Latitude, 
                 shape = Sample_Type_Category),
             color = "white",
             size = 5,            # Much smaller size
             alpha = 1,
             stroke = 2) +
  
  # Add the colored points on top
  geom_point(data = map_data,
             aes(x = Longitude, y = Latitude, 
                 color = Country, 
                 shape = Sample_Type_Category),
             size = 4,            # Much smaller size
             alpha = 0.95,
             stroke = 1.5) +
  
  # Add site labels with sample counts - improved positioning
  geom_text_repel(data = map_data,
                  aes(x = Longitude, y = Latitude, 
                      label = Site_Label),
                  size = 4,
                  fontface = "bold",
                  color = "black",
                  bg.color = "white",
                  bg.r = 0.25,
                  box.padding = 2,      # Much more padding
                  point.padding = 1,    # More space from points
                  segment.color = "black",
                  segment.size = 0.8,
                  segment.alpha = 0.9,
                  arrow = arrow(length = unit(0.03, "npc"), 
                                type = "closed"),
                  max.overlaps = Inf,
                  force = 0.75,            # Stronger force to separate labels
                  force_pull = 0.5,     # Less pull back to points
                  min.segment.length = 0.1,
                  xlim = c(-12, 15),    # Allow labels to extend beyond map
                  ylim = c(52, 66),     # Allow labels to extend beyond map
                  direction = "both") + # Allow movement in both directions
  
  # Set coordinate limits
  coord_sf(xlim = lon_limits, ylim = lat_limits, expand = FALSE) +
  
  # Color scheme: Red for Norway, Blue for Scotland
  scale_color_manual(values = c("Norway" = "#CC0000",     # Darker red
                                "Scotland" = "#0066CC"),   # Darker blue
                     name = "Country") +
  
  # Shape scheme for sample types (only two types now)
  scale_shape_manual(values = c("Gill Tissue" = 16,    # solid circle
                                "Gill Swab" = 17),     # triangle
                     name = "Sample Type") +
  
  # Labels and title
  labs(title = "",
       x = "Longitude",
       y = "Latitude") +
  
  # Theme
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 20),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 15, r = 15, b = 50, l = 15)
  ) +
  
  # Simplified guides (no size legend)
  guides(
    color = guide_legend(override.aes = list(size = 3.5),
                         title.position = "top",
                         title.hjust = 0.5,
                         order = 1),
    shape = guide_legend(override.aes = list(size = 3.5),
                         title.position = "top", 
                         title.hjust = 0.5,
                         order = 2)
  )

# Print the map
print(p)

# Print summary to verify the sample types are correct
cat("\n=== SAMPLING SITE SUMMARY ===\n")
cat("Total sites:", nrow(map_data), "\n")
print(map_data %>% 
        select(Site, Country, Sample_Type_Category, Sample_Count, Site_Label) %>%
        arrange(Country, Site))

cat("\nSample types by site:\n")
gill_tissue_sites <- map_data %>% filter(Sample_Type_Category == "Gill Tissue")
cat("Gill Tissue sites:", paste(gill_tissue_sites$Site, collapse = ", "), "\n")
cat("All other sites: Gill Swab\n")

cat("\nFinal map saved as:\n")
cat("- north_atlantic_sampling_sites_FINAL.png\n")
cat("- north_atlantic_sampling_sites_FINAL.pdf\n")
