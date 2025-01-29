
#To determine the missing sites of each clocks:
#First download the clocks cpg sites, available in Biolearn: https://bio-learn.github.io/clocks.html

# Load required packages
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)  # Official hg38 annotation

#anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)


# Modified probe retrieval function with EPICv2 normalization
get_array_probes <- function(array_type) {
  switch(array_type,
         "450k" = {
           anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
           unique(c(
             anno$Name,
             anno@listData$Methyl27_Loci
           ))
         },
         "EPICv1" = {
           anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
           unique(anno$Name)
           unique(c(
             anno$Name,
             anno@listData$Methyl450_Loci,
             anno@listData$Methyl27_Loci
           ))
         },
         "EPICv2" = {
           # Get EPICv2 probes and normalize names
           anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
           unique(c(
             gsub("_.*$", "", anno$Name),
             anno@listData$EPICv1_Loci,
             anno@listData$Methyl450_Loci,
             anno@listData$Methyl27_Loci
           ))
         })
}

#anno@listData$EPICv1_Loci



array_probes <- list(
  "450k" = get_array_probes("450k"),
  "EPICv1" = get_array_probes("EPICv1"),
  "EPICv2" = get_array_probes("EPICv2")
)

# Function to check CpG presence
check_cpg_presence <- function(cpg_sites, array_probes) {
  # Remove non-CpG coefficient in Biolearn
  non_cpg_terms <- c("(Intercept)", "Age", "Female", "Intercept", "m_age", 
                     "m_cox", "sd_age", "sd_cox" ,"DNAmADM", "DNAmB2M","DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPACKYRS","DNAmPAI1", "DNAmTIMP1", "DNAmlogA1C", "DNAmlogCRP"
  )
  
  cpg_sites <- setdiff(cpg_sites, non_cpg_terms)
  total <- length(cpg_sites)
  missing <- setdiff(cpg_sites, array_probes)
  
  list(
    total = total,
    present = total - length(missing),
    missing_count = length(missing),
    missing_percent = round(length(missing)/total * 100, 2),
    missing_sites = missing
  )
}


# List of clocks (update with your actual data)
clocks <- list(
  "Horvath1" = Horvath1_coef$CpGmarker,
  "Hannum" = Hannum_coef$CpGmarker,
  "PhenoAge" = PhenoAge_coef$CpGmarker,
  "GrimAge" = GrimAgeV1_coef$var,
  "Horvath2" = Horvath2_coef$CpGmarker,
  "DNAmTL" = DNAmTL_coef$CpGmarker,
  "DNAmFitAge" = DNAmFitAge_coef,
  "GrimAge2" = GrimAgeV2_coef$var,
  "DunedinPoAm" = DunedinPoAm38_coef$CpGmarker,
  "DunedinPACE" = DunedinPACE_coef$CpGmarker,
  "YingCausAge" = YingCaus_coef$CpGmarker,
  "YingAdaptAge" = YingAdapt_coef$CpGmarker,
  "YingDamAge" = YingDam_coef$CpGmarker
)

# Generate missing CpG report
generate_missing_report <- function(clocks, array_probes) {
  report <- list()
  
  for (clock_name in names(clocks)) {
    clock_cpgs <- clocks[[clock_name]]
    clock_report <- list()
    
    for (array_name in names(array_probes)) {
      result <- check_cpg_presence(clock_cpgs, array_probes[[array_name]])
      
      clock_report[[array_name]] <- list(
        missing_count = result$missing_count,
        missing_percent = result$missing_percent,
        missing_sites = result$missing_sites
      )
    }
    
    report[[clock_name]] <- clock_report
  }
  
  return(report)
}

# Generate full report
missing_report <- generate_missing_report(clocks, array_probes)

# Print summary table
cat("\nMissing CpGs Summary:\n")
cat("-----------------------------------------\n")
cat("Clock Name            | 450k | EPICv1 | EPICv2\n")
cat("-----------------------------------------\n")
for (clock in names(missing_report)) {
  counts <- sapply(missing_report[[clock]], function(x) paste0(x$missing_count, " (", x$missing_percent, "%)"))
  cat(sprintf("%-20s | %-8s | %-8s | %-8s\n", 
              clock, 
              counts["450k"],
              counts["EPICv1"],
              counts["EPICv2"]))
}
cat("-----------------------------------------\n")






