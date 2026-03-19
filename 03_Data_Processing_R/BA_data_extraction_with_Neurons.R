suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# ------------------------------
# 1. LOAD FILES
# ------------------------------

base_path <- "/Users/izgivarol/Documents/Heid Uni/BA/Results/ROI Analysis/Cell Profiler/Analysis"

syp      <- read_csv(file.path(base_path,"Overlap_BASYP_puncta.csv"), show_col_types = FALSE)
strict   <- read_csv(file.path(base_path,"Overlap_BASYP_on_BDA_strict.csv"), show_col_types = FALSE)
loose    <- read_csv(file.path(base_path,"Overlap_BASYP_on_BDA_loose.csv"), show_col_types = FALSE)
neurons  <- read_csv(file.path(base_path,"Overlap_BANeurons.csv"), show_col_types = FALSE)
image_df <- read_csv(file.path(base_path,"Overlap_BAImage.csv"), show_col_types = FALSE)

# ------------------------------
# 2. CLEAN METADATA
# ------------------------------

clean_meta <- function(df){
  
  cohort_col <- grep("Metadata_Cohort", names(df), value=TRUE)[1]
  animal_col <- grep("Metadata_Animal", names(df), value=TRUE)[1]
  sample_col <- grep("Metadata_Sample", names(df), value=TRUE)[1]
  
  df %>%
    mutate(
      Cohort = as.character(.data[[cohort_col]]),
      Animal = as.character(.data[[animal_col]]),
      Sample = as.character(.data[[sample_col]])
    )
  
}

syp     <- clean_meta(syp)
strict  <- clean_meta(strict)
loose   <- clean_meta(loose)
neurons <- clean_meta(neurons)

# ------------------------------
# 3. COUNT OBJECTS PER IMAGE
# ------------------------------

total_syp <- syp %>%
  group_by(ImageNumber,Cohort,Animal,Sample) %>%
  summarise(total_syp=n(), .groups="drop")

strict_syp <- strict %>%
  group_by(ImageNumber,Cohort,Animal,Sample) %>%
  summarise(strict_syp_on_bda=n(), .groups="drop")

loose_syp <- loose %>%
  group_by(ImageNumber,Cohort,Animal,Sample) %>%
  summarise(loose_syp_on_bda=n(), .groups="drop")

neun_summary <- neurons %>%
  group_by(ImageNumber,Cohort,Animal,Sample) %>%
  summarise(neun_cells=n(), .groups="drop")

# ------------------------------
# 4. BUILD IMAGE SUMMARY TABLE
# ------------------------------

image_summary <- total_syp %>%
  full_join(strict_syp, by=c("ImageNumber","Cohort","Animal","Sample")) %>%
  full_join(loose_syp,  by=c("ImageNumber","Cohort","Animal","Sample")) %>%
  left_join(neun_summary, by=c("ImageNumber","Cohort","Animal","Sample")) %>%
  mutate(
    
    strict_syp_on_bda = replace_na(strict_syp_on_bda,0),
    loose_syp_on_bda  = replace_na(loose_syp_on_bda,0),
    neun_cells        = replace_na(neun_cells,0),
    
    strict_fraction = ifelse(total_syp>0, strict_syp_on_bda/total_syp, NA),
    loose_fraction  = ifelse(total_syp>0, loose_syp_on_bda/total_syp, NA),
    
    strict_synapses_per_neuron = ifelse(neun_cells>0, strict_syp_on_bda/neun_cells, NA),
    loose_synapses_per_neuron  = ifelse(neun_cells>0, loose_syp_on_bda/neun_cells, NA)
    
  ) %>%
  arrange(Cohort,as.numeric(Animal),Sample,ImageNumber)

image_summary <- image_summary %>%
  distinct(Cohort, Animal, Sample, .keep_all = TRUE)

# ------------------------------
# 6. CHECK RESULTS
# ------------------------------

View(image_summary)

image_summary %>%
  count(Cohort,Animal,Sample)

# ------------------------------
# 7. SAVE FINAL TABLES
# ------------------------------

write_csv(image_summary, file.path(base_path,"FINAL_image_summary.csv"))
