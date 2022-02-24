### Prepare data for NCBI submission 

files <- list.files(paste("//NRONP6AwvFSP001/EnvChem/project_data_Venier/AshNet_biodiversity/sequence.data/raw.sequence.data"), recursive = T, full.names=T) %>%
        str_subset("\\.zip$")

files <- sapply(files, function(x){unzip(x, list=TRUE)})

files <- unlist(files) %>% str_subset("q.gz$")

NCBI <- data.frame(file=files, set = ifelse(grepl("16S", files), "16S V4-V5", "18S, ITS, BR5, F230"))

d16S_sampID <-  read.csv(paste("//NRONP6AwvFSP001/EnvChem/project_data_Venier/AshNet_biodiversity/sequence.data", "AshNet.16S.sample.ID.csv", sep="/"))

NCBI$Label <- gsub("Caroline Emilson 16S V4-V5\\/|Pinkn*ey-Ashnet-", "", NCBI$file) %>% 
        str_extract(".+_S[[:digit:]]+") %>% 
        str_remove("_S[[:digit:]]+$")

NCBI$Label <- ifelse(NCBI$Label %in% d16S_sampID$Sample.no, d16S_sampID$Sample.name[as.numeric(NCBI$Label)], NCBI$Label) 

load("data/AshNet_Clean_meta_func.RData")

Sam_data4join <- Sample_metadata %>% left_join(metadata %>% select(Site_name, Province), by="Site_name")

Sam_data4join$Label4join <- str_replace(Sam_data4join$Label, "[[:alpha:]]$", "") 

NCBI$Label4join <- str_replace(NCBI$Label, "[[:alpha:]]$", "")

NCBI$sample_name <- NCBI$Label4join

Sam_data4join <- unique(Sam_data4join[, c("Label4join", "Site", "ash_amt", "Replicate",  "ash_type", "Province", "Soil Sample Type", "Forest.floor.cm", "X_Longitude", "Y_Latitiude", "Z_Elevation")])

NCBI <- left_join(NCBI, Sam_data4join, by= "Label4join")
        
NCBIbiosample <- NCBI %>%
        mutate(lat_long = paste(paste(round(abs(Y_Latitiude), 4), ifelse(Y_Latitiude < 0, "S", "N")), paste(round(abs(X_Longitude), 4), ifelse(X_Longitude < 0, "W", "E"))), 
               depth = ifelse(`Soil Sample Type` == "MIN", "0-10 cm", paste(str_replace(Forest.floor.cm, "cm| ", ""), "cm"))) %>%
        rename(elev = Z_Elevation) %>%
        mutate(organism = "soil metagenome", 
               collection_date = 2017, 
               env_broad_scale = "soil[ENVO:00001998]", 
               env_local_scale = "forest soil[ENVO:00002261]",
               env_medium = ifelse(`Soil Sample Type` == "MIN", "mineral horizon[ENVO:03600011]",
                                   ifelse(grepl("FH", `Soil Sample Type`),
                                          "organic horizon[ENVO:03600018]", 
                                          "litter layer[ENVO:01000338]")),
               pool_dna_extracts = 3, 
               geo_loc_name = paste("Canada", Province, sep=": " 
                                    ),  
               bioproject_accession = "",
               agrochem_addition = ifelse(!is.na(ash_type), paste(ash_amt, "t/ha", ash_type, "wood ash"), ""), 
               sample_title = "",
               Replicate = paste("replicate ", Replicate)) %>%
        select(sample_name, sample_title, bioproject_accession, organism, collection_date, depth, elev, env_broad_scale, env_local_scale, env_medium, geo_loc_name, lat_long, agrochem_addition, Replicate) %>%
        unique()

write.csv(NCBIbiosample, "data/NCBIbs_data.csv")

NCBISRA <- NCBI %>% mutate(library_id = gsub("Caroline Emilson 16S V4-V5\\/|_L001_R._001.fastq.gz", "", file)) %>%
        select(sample_name, library_id) %>% unique()
        
write.csv(NCBISRA, "data/NCBIsra_data.csv")



files <- list.files(paste("//NRONP6AwvFSP001/EnvChem/project_data_Venier/AshNet_biodiversity/sequence.data/raw.sequence.data"), recursive = T, full.names=T) %>%
        str_subset("\\.zip$")
dir.create("temp")
for(i in files){
        unzip(i, exdir="temp")
}
