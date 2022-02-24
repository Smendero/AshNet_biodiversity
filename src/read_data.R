## Script to read in datasets from the ashnet dataset from the EnvChem Drive##

library(readxl)
library(tidyverse)
source("src/helper_functions.R")

#### Compile site characteristics data ####
path_to_files <- "E://project_data_Venier/AshNet_biodiversity/site.metadata"
metadata <- read_excel(paste(path_to_files, "AshNet.metadata.xlsx", sep="/"), skip=8)
colnames(metadata)[1] <- "Site_name"
## cut down to columns of interest
metadata <- metadata[, !colnames(metadata) %in% c("Lead Contributor",
                                                  "Associated contacts", 
                                                  "Other.treatment", 
                                                  "Treatment", 
                                                  "Ash.type", 
                                                  "Ash.Ca", 
                                                  "Ca.kg/ha*application rate", 
                                                  "Ash.application.rate", 
                                                  "...39", 
                                                  "Coordinates", 
                                                  "Elevation.m")]

metadata$Year.ash.applied <- as.character(ifelse(grepl("[[:alpha:]]", 
                                          metadata$Year.ash.applied), 
                                    metadata$Year.ash.applied, 
                                    as.character(as.Date(as.numeric(metadata$Year.ash.applied), 
                                                         origin="1899-12-30"))))
metadata$Year.ash.applied <- str_extract(metadata$Year.ash.applied, "[[:digit:]]{4}")
metadata$Years_ash <- 2017 - as.numeric(metadata$Year.ash.applied)

metadata$Stand <- ifelse(grepl("spruce", tolower(metadata$`Stand type`)), 
                               "Spruce", 
                               ifelse(grepl("pine", tolower(metadata$`Stand type`)),
                                            "Pine", 
                                            ifelse(grepl("poplar", tolower(metadata$`Stand type`)), 
                                                  "Poplar", 
                                                  "Mixed deciduous")))

metadata$Oldest_Stand_age <- c(2017-1997, 2017-1991, 2017-1995, 2017-2015, 2017-2012, 2017-2013, 80, 80, 2017-2012)

#### Pull Climate Data ####
climate_data_2016 <- read_excel(paste(path_to_files, "AshNet - climate data - Sept 2016 - from Pia.xls", sep ="/"), skip=2)
colnames(climate_data_2016) <- read_excel(paste(path_to_files, "AshNet - climate data - Sept 2016 - from Pia.xls", sep ="/"))[1,]

climate_data_2016 <- climate_data_2016[c(1, 8:14),]

climate_data_2016$Site_name <- c("Eastern Township Sugar Maple", "Haliburton", "Island Lake", "25th sideroad (Lakehead)", "Pineland", "Mistik (Burness)", "Aleza Lake S", "Aleza Lake N")

climate_data_2016 <- climate_data_2016[, !colnames(climate_data_2016) == "ID"]

## manually Went down "leaves" of grouped variables 
## until pulled single variable that is highly correlated with related variables (.9 or higher).
## use spearmans correlation because of the small sample #'s and skewed distribution of variables.

finalgroups <- Hclust_dimension_reduc(climate_data_2016[, !colnames(climate_data_2016) %in% c("Site_name", "ID", "X_Longitude", "Y_Latitiude", "Z_Elevation")], 
                                      0.9, 
                                      "intermediate/Climate Representative variables for analyses.csv", method="spearman")

## Subset climate data to new groups
climate_data_2016 <- climate_data_2016[, colnames(climate_data_2016) %in% c("Site_name", "ID", "X_Longitude", "Y_Latitiude", "Z_Elevation", unique(finalgroups$rep_var))]

rm(finalgroups, group_num, i, g, group_cor_chk, groups, climate_groups, climate_cor, climate_clust, captured, climate_dist, temporary_groups, threshold)


### Grab Sample info ####
path_to_files <- "E://project_data_Venier/AshNet_biodiversity/project.planning.files"

Samples <- read_excel(paste(path_to_files, "Master List for AshNet Metabarcoding Study.xlsx", sep="/"), skip =4)

Samples <- Samples[, c(1:5)]

#### Populate full names ####
Samples$Site_name <- sapply(as.character(Samples$Site), function(x){
        if(x == "ALN"){return("Aleza Lake N")}
        if(x == "ALS"){return("Aleza Lake S")}
        if(x == "ETM"){return("Eastern Township Sugar Maple")}
        if(x == "SRD"){return("25th sideroad (Lakehead)")}
        if(x == "ILK"){return("Island Lake")}
        if(x == "HLB"){return("Haliburton")}
        if(x == "PLD"){return("Pineland")}
        if(x == "MSK"){return("Mistik (Burness)")}
})
## populate soil types
Samples$Soil_type <- sapply(as.character(Samples$`Soil Sample Type`), function(x){
        if(x == "LM"){return("surface litter &/or moss")}
        if(x == "FH"){return("FH-layer forest floor")}
        if(x == "LFH"){return("surface litter &/or moss with FH layer")}
        if(x == "MIN"){return("0-10 cm mineral soil")}
})
## Populate ash amounts
Samples$ash_amt <- sapply(1:nrow(Samples), function(x){
        if(Samples$Treatment[x] == "C"){return(0)}
        if(grepl("AL", Samples$Site[x]) & !Samples$Treatment[x] == "C"){
                return(5)
        }
        if(Samples$Site[x] == "PLD" & Samples$Treatment[x] == "A"){
                return(1.5)
        }
        if(Samples$Site[x] == "HLB" & Samples$Treatment[x] %in% c("A1", "A4")){
                return(1)
        }
        if(Samples$Site[x] == "HLB" & Samples$Treatment[x] %in% c("A2", "A5")){
                return(4)
        }
        if(Samples$Site[x] == "HLB" & Samples$Treatment[x] %in% c("A3", "A6")){
                return(8)
        }
        if(Samples$Site[x] == "ILK" & Samples$Treatment[x] %in% c("A1")){
                return(0.7)
        }
        if(Samples$Site[x] == "ILK" & Samples$Treatment[x] %in% c("A2")){
                return(1.4)
        }
        if(Samples$Site[x] == "ILK" & Samples$Treatment[x] %in% c("A3")){
                return(2.8)
        }
        if(Samples$Site[x] == "ILK" & Samples$Treatment[x] %in% c("A4")){
                return(5.6)
        }
        if(Samples$Site[x] == "SRD" & Samples$Treatment[x] %in% c("A1")){
                return(1)
        }
        if(Samples$Site[x] == "SRD" & Samples$Treatment[x] %in% c("A2")){
                return(10)
        }
        if(Samples$Site[x] == "ETM" & Samples$Treatment[x] %in% c("A")){
                return(20)
        }
        if(Samples$Site[x] == "MSK" & Samples$Treatment[x] %in% c("A1")){
                return(1)
        }
        if(Samples$Site[x] == "MSK" & Samples$Treatment[x] %in% c("A2")){
                return(5)
        }
})
## Populate Calcium Application rate (kg/ha)
Samples$Ca_application_amt <- sapply(1:nrow(Samples), function(x){
        if(Samples$Treatment[x] == "C"){return(0)}
        if(grepl("AL", Samples$Site[x]) & Samples$Treatment[x] == "A1"){
                return(335)
        }
        if(grepl("AL", Samples$Site[x]) & Samples$Treatment[x] == "A2"){
                return(970)
        }
        if(Samples$Site[x] == "PLD" & Samples$Treatment[x] == "A"){
                return(273.495)
        }
        if(Samples$Site[x] == "HLB" & Samples$Treatment[x] %in% c("A1")){
                return(43.6)
        }
        if(Samples$Site[x] == "HLB" & Samples$Treatment[x] %in% c("A4")){
                return(101.1)
        }
        if(Samples$Site[x] == "HLB" & Samples$Treatment[x] %in% c("A2")){
                return(174.4)
        }
        if(Samples$Site[x] == "HLB" & Samples$Treatment[x] %in% c("A5")){
                return(404.4)
        }
        if(Samples$Site[x] == "HLB" & Samples$Treatment[x] %in% c("A3")){
                return(348.8)
        }
        if(Samples$Site[x] == "HLB" & Samples$Treatment[x] %in% c("A6")){
                return(808.8)
        }
        if(Samples$Site[x] == "ILK" & Samples$Treatment[x] %in% c("A1")){
                return(50)
        }
        if(Samples$Site[x] == "ILK" & Samples$Treatment[x] %in% c("A2")){
                return(100)
        }
        if(Samples$Site[x] == "ILK" & Samples$Treatment[x] %in% c("A3")){
                return(200)
        }
        if(Samples$Site[x] == "ILK" & Samples$Treatment[x] %in% c("A4")){
                return(400)
        }
        if(Samples$Site[x] == "SRD" & Samples$Treatment[x] %in% c("A1")){
                return(170.28)
        }
        if(Samples$Site[x] == "SRD" & Samples$Treatment[x] %in% c("A2")){
                return(1702.8)
        }
        if(Samples$Site[x] == "MSK" & Samples$Treatment[x] %in% c("A1")){
                return(16.84)
        }
        if(Samples$Site[x] == "MSK" & Samples$Treatment[x] %in% c("A2")){
                return(84.2)
        }
        if(grepl("ETM", Samples$Site[x]) & Samples$Treatment[x] == "A"){
                return(3516)
        }
})
## Populate ash Type information
Samples$ash_type <- sapply(1:nrow(Samples), function(x){
        if(Samples$Treatment[x] == "C"){return(NA)}
        if(grepl("AL", Samples$Site[x]) & Samples$Treatment[x] == "A1"){
                return("UNBC Bottom")
        }
        if(grepl("AL", Samples$Site[x]) & Samples$Treatment[x] == "A2"){
                return("CPLP Bottom")
        }
        if(Samples$Site[x] == "PLD" & Samples$Treatment[x] %in% c("A")){
                return("Mixed")
        }
        if(Samples$Site[x] == "HLB" & Samples$Treatment[x] %in% c("A1", "A2", "A3")){
                return("Bottom")
        }
        if(Samples$Site[x] == "HLB" & Samples$Treatment[x] %in% c("A4", "A5", "A6")){
                return("Fly")
        }
        if(Samples$Site[x] == "MSK" & Samples$Treatment[x] %in% c("A1", "A2")){
                return("Bottom")
        }
        if(Samples$Site[x] == "SRD" & Samples$Treatment[x] %in% c("A1", "A2")){
                return("Fly") ### Assuming just fly ash
        }
        if(Samples$Site[x] == "ETM" & !Samples$Treatment[x] %in% c("C")){
                return("Bottom") 
        }
        if(Samples$Site[x] == "ILK" & !Samples$Treatment[x] %in% c("C")){
                return("Bottom") 
        }
})

Samples$`Soil Sample Type` <- ifelse(Samples$`Soil Sample Type` == "LM", "L", Samples$`Soil Sample Type`)
#### Pull in Ash information ####
path_to_files <- "E://project_data_Venier/AshNet_biodiversity/wood.ash.chemistry.data"

Ashchem <- read_excel(paste(path_to_files, "AshNet.Ash.chemistries.xlsx", sep="/"))[1:11,]

Ashchem$Ash.type[1:2] <- Ashchem$Ash.type[3:4]

colnames(Ashchem)[1:2] <- c("Site_name", "ash_type")

Ashchem$Site_name[grep("25th Sideroad", Ashchem$Site_name)] <- "25th sideroad (Lakehead)" 
Ashchem$Site_name[grep("Mistik", Ashchem$Site_name)] <- "Mistik (Burness)" 

Ashchem[, !colnames(Ashchem) %in% c("Source", "Site_name", "ash_type", "Ash.Mg.g.kg")] <- apply(Ashchem[, !colnames(Ashchem) %in% c("Source", "Site_name", "ash_type", "Ash.Mg.g.kg")],
                                                                                                2,                                                                                       function(x) as.numeric(trimws(gsub("\\(|\\)|-|[[:alpha:]]", "", x))))
Ashchem <- Ashchem[, !colnames(Ashchem) %in% c("Ash.Al.g.kg", "Ash.Cs.g.kg", "Source", "Ash.Mg.g.kg")]

## recalculate as g/Mg and replace
Ashchem[, grep("Ash", colnames(Ashchem))] <- Ashchem[, grep("Ash", colnames(Ashchem))] * 1000

colnames(Ashchem)[grep("Ash", colnames(Ashchem))] <- gsub("g\\.kg", "g.Mg", colnames(Ashchem[, grep("Ash", colnames(Ashchem))]))

### Merge datasets

Sample_metadata <- left_join(Samples, metadata[, c("Site_name", "Years_ash", "Soil.type", "Forest.floor.cm", "Oldest_Stand_age", "Stand")], by = "Site_name")
Sample_metadata <- left_join(Sample_metadata, climate_data_2016, by="Site_name")
Sample_metadata <- left_join(Sample_metadata, Ashchem, by=c("Site_name", "ash_type"))

Sample_metadata$Years_ash <- ifelse(Sample_metadata$Treatment=='C', 0, Sample_metadata$Years_ash)
## recalculate as kg element added and replace
Sample_metadata[, grepl("Ash", colnames(Sample_metadata)) & !colnames(Sample_metadata) == "Ash.feedstock"] <- Sample_metadata[, grepl("Ash", colnames(Sample_metadata)) & !colnames(Sample_metadata) == "Ash.feedstock"] /1000 * Sample_metadata$ash_amt

colnames(Sample_metadata)[grep("Ash", colnames(Sample_metadata))] <- gsub("g\\.Mg", "kg.applied", colnames(Sample_metadata[, grep("Ash", colnames(Sample_metadata))]))

rm(Samples, climate_data_2016, path_to_files, Hclust_dimension_reduc)


#### Read in Chemistry Data #####
path_to_files <- "E://project_data_Venier/AshNet_biodiversity/soil.chemistry.data"
TotalChemistry <- read_excel(paste(path_to_files, "Ruth's data compiled.xlsx", sep="/"), sheet='total carbon and nitrogen')
TotalChemistry <- TotalChemistry[, c("Label", "C%", "N%", "C:N")]

MBCChemistry <- read_excel(paste(path_to_files, "Ruth's data compiled.xlsx", sep="/"), sheet='MBC', col_types = c("numeric", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
MBCChemistry <- MBCChemistry[, c("Label", "MBC (mg C/g soil)")]

Moisture <- read_excel(paste(path_to_files, "Ruth's data compiled.xlsx", sep="/"), sheet='moisture')
Moisture <- Moisture[, c("Label", "Moisture content (% DW)")]

TotalChemistry <- left_join(TotalChemistry, MBCChemistry, by="Label")
TotalChemistry <- left_join(TotalChemistry, Moisture, by="Label")

rm(Moisture, MBCChemistry)

## Combine with Sample meta-data, applying to replicates in the sample metadata dataset
chemcomb_term <- "[[:alpha:]]+-[[:alnum:]]+-[[:alpha:]]+-[[:digit:]]"

Sample_metadata$Chem_combine <- str_extract(Sample_metadata$Label, chemcomb_term)

TotalChemistry <- TotalChemistry %>%
        mutate(Chem_combine = str_extract(Label, chemcomb_term) )%>%
        select(-Label) %>%
        data.frame()

Sample_metadata <- left_join(Sample_metadata, TotalChemistry, by = "Chem_combine") %>%
        select(-Chem_combine) %>%
        mutate(pH_combine = ifelse(grepl("ILK", Label), Label, str_extract(Label, chemcomb_term)))

pHData <- read_excel(paste(path_to_files, "AshNet_pHdata.xlsx", sep="/"), skip=3)
colnames(pHData) <- c("Sample_Name", "Sample_number", "Lab_ID", "Weight", "Volume", "pH", "Notes")
pHData <- pHData[, c("Sample_Name", "pH")]

pHData <- pHData %>%
        mutate(pH_combine = ifelse(grepl("ILK", Sample_Name), Sample_Name, str_extract(Sample_Name, chemcomb_term)))%>%
        select(-Sample_Name) %>%
        data.frame()

Sample_metadata <- left_join(Sample_metadata, pHData, by = "pH_combine") %>%
        select(-pH_combine)

#### Read in Enzyme Data ####
path_to_files <- "E://project_data_Venier/AshNet_biodiversity/enzyme.data"
Hydrolases <- read.csv(paste(path_to_files, "AshNet.hydrolases.csv", sep="/"))
Hydrolases <- Hydrolases[!is.na(Hydrolases$sample.name), ]
Hydrolases$Site <- gsub("\\*", "", str_extract(Hydrolases$sample.name, "^.[[:alpha:]]+"))
Hydrolases$Treatment <- gsub(".[[:alpha:]]+-", "", str_extract(Hydrolases$sample.name, "^.[[:alpha:]]+-[[:alnum:]]+"))

ggplot(Hydrolases, aes(Treatment, activity.nmol.h.g..wet.weight., color=Treatment))+
        geom_boxplot()+
        geom_rug()+
        facet_grid(enzyme~Site, scales = 'free')

Lignases <- read.csv(paste(path_to_files, "AshNet.lignases.csv", sep="/"))
Lignases$Site <- gsub("\\*", "", str_extract(Lignases$sample.name, "^.[[:alpha:]]+"))
Lignases$Treatment <- gsub(".[[:alpha:]]+-", "", str_extract(Lignases$sample.name, "^.[[:alpha:]]+-[[:alnum:]]+"))

ggplot(Lignases, aes(Treatment, activity.umol.h.g..wet.weight., color=Treatment))+
        geom_boxplot()+
        geom_rug()+
        facet_grid(enzyme~Site, scales = 'free')

## High amount of results negative, not including at this time.

hyd_lign_soil_moist <- read_excel(paste(path_to_files, "AshNet_2017-18_SoilMoisture.xlsx", sep="/"), skip = 1)
hyd_lign_soil_moist <- hyd_lign_soil_moist %>%
        mutate(Label = gsub("\\*",  "", `Sample Name`))%>%
        select(Label, `Moisture content (% DW)`, `% Soil Moisture (wet mass - dry mass) / wet mass *100`) %>%
        data.frame()

moist_compare <- left_join(Sample_metadata[, c("Label", "Moisture.content....DW.")], hyd_lign_soil_moist, by = "Label")

ggplot(moist_compare, aes(Moisture.content....DW..x, Moisture.content....DW..y))+
        geom_point()+
        stat_smooth(method="lm")
## Results highly related to results from higher sample

## Correct for soil moisture - get activities as dry weight
colnames(hyd_lign_soil_moist) <- c("Label", "Moisture_DW", "Soil_Moisture")

Hydrolases <- Hydrolases %>% 
        filter(!enzyme == "BG") %>% ## High amount of results negative, not including at this time.
        mutate(Label = gsub("\\*", "", sample.name)) %>%
        select(-sample.name, - Site, -Treatment) %>%
        left_join(hyd_lign_soil_moist, by="Label") %>%
        mutate(activity_nmol.h.g_g.dry.wt = activity.nmol.h.g..wet.weight. * (1-(Soil_Moisture/100))) %>%
        data.frame()

Hydrolases_comb <- Hydrolases %>% 
        filter(!enzyme == "BG") %>% ## High amount of results negative, not including at this time.
        pivot_wider(id_cols = Label, names_from = enzyme, values_from = activity_nmol.h.g_g.dry.wt, values_fn = mean) %>%
        data.frame()

Sample_metadata <- left_join(Sample_metadata, Hydrolases_comb, by="Label") ## need to recalculate hydrolases based on dry weight

rm(Lignases, Hydrolases_comb, moist_compare, hyd_lign_soil_moist, pHData, TotalChemistry, chemcomb_term, path_to_files)


#### Read in ITS Data ####
path_to_files <- "E://project_data_Venier/AshNet_biodiversity/sequence.data"

ITS <- read.csv(paste(path_to_files, "AshNet.biodiversity.ITS.ESVtaxr.csv", sep="/")) ## only using rarefied table with taxonomies
ITS <- ITS %>%
        rename(zOTU = X)
ITSraw <- read.csv(paste(path_to_files, "AshNet.biodiversity.ITS.results.csv", sep="/")) ## only using rarefied table with taxonomies

#### Read in 18S Data ####
d18S <-  read.csv(paste(path_to_files, "AshNet.biodiversity.18S.ESVtaxr.csv", sep="/")) ## only using rarefied table with taxonomies
d18S <- d18S %>%
        rename(zOTU = X)
d18Sraw <- read.csv(paste(path_to_files, "AshNet.biodiversity.18S.results.csv", sep="/")) ## only using rarefied table with taxonomies

#### Read in F230 Data ####
F230 <-  read.csv(paste(path_to_files, "AshNet.biodiversity.F230.ESVtaxr.csv", sep="/")) ## only using rarefied table with taxonomies
F230 <- F230 %>%
        rename(zOTU = X)
F230raw <- read.csv(paste(path_to_files, "AshNet.biodiversity.F230.results.csv", sep="/")) ## only using rarefied table with taxonomies

#### Read in 16S Data ####
d16S <-  read.csv(paste(path_to_files, "AshNet.biodiversity.16S.ESVtaxr.csv", sep="/")) ## only using rarefied table with taxonomies
d16S_sampID <-  read.csv(paste(path_to_files, "AshNet.16S.sample.ID.csv", sep="/")) ## only using rarefied table with taxonomies
d16S <- d16S %>%
        rename(zOTU = X)
colnames(d16S)[2:172] <- d16S_sampID$Sample.name
d16Sraw <- read.csv(paste(path_to_files, "AshNet.biodiversity.16S.results.csv", sep="/")) ## only using rarefied table with taxonomies
d16Sraw$SampleName <- d16S_sampID$Sample.name[d16Sraw$SampleName]

#### Save Datasets to RData object

save(Hydrolases, metadata, Sample_metadata, file="data/AshNet_Clean_meta_func.RData")
save(ITS, ITSraw, file="data/AshNet_ITS.RData")
save(d18S, d18Sraw, F230, F230raw, file="data/AshNet_Arthro_Seq.RData")
save(d16S, d16Sraw, file="data/AshNet_16S.RData")
