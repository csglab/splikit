source("/home/arsham79/Splicing-Pipeline/R/gedi_input_matrices.R")
suppressPackageStartupMessages({
  
  library(dplyr)
  library(Matrix)
  library(data.table)
  
})


sample_metadata <- readRDS("/home/arsham79/scratch/Sandrine-s-project/data/meta.rds") |> as.data.table()

sample_ids <- c("GCRC2076_2_1of2",
                "GCRC2076_2_2of2",
                "GCRC2076_Rb1_1of2",
                "GCRC2076_Rb1_2of2",
                "GCRC2076_RD_1of2", 
                "GCRC2076_RD_2of2",
                "GCRC2076_RES_1of2",
                "GCRC2076_RES_2of2")

barcode <- list()
for(sample in sample_ids){
  temp_bar <- sample_metadata[donor_sample == sample, Barcode]
  barcode[[sample]] <- sub("(^.{16})-1", "\\1", temp_bar)
  
}

junction_ab_object <- multigedi_make_junction_ab(STARsolo_SJ_dirs = paste0(list.files("/home/arsham79/scratch/Sandrine-s-project/results/star-output-new-whitelist/", full.names = TRUE), 
                                                                           "/Solo.out/SJ/") |> as.list(), 
                                                 
                                                 sample_ids = sample_ids |> as.list(), 
                                                 white_barcode_lists =  barcode, 
                                                 use_internal_whitelist = FALSE)


m1 <- multigedi_make_m1(junction_ab_object = junction_ab_object)

m1_old <- readRDS("/home/arsham79/scratch/m1_merged.rds")
m1_old <- m1_old[rownames(m1$m1_inclusion_matrix), colnames(m1$m1_inclusion_matrix)]

identical(m1_old, m1$m1_inclusion_matrix)


multigedi_countsplit(m1_inclusion_matrix = m1$m1_inclusion_matrix)
m2 <- multigedi_make_m2(m1_inclusion_matrix = m1$m1_inclusion_matrix, eventdata = m1$event_data)



m2_old <- readRDS("/home/arsham79/scratch/m2_merged.rds")
m2_old <- m2_old[rownames(m2), colnames(m2)]
identical(m2, m2_old)

