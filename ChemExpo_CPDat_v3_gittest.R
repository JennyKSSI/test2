# AUTHOR: Jenny Kay, with support from Anaga Dinesh
# PURPOSE: classify sources of exposure to breast cancer-relevant chemicals 
#          identified in Kay et al EHP 2024 using EPA's Chemicals and Products Database (CPDat)
#          Bulk data downloaded from ChemExpo (beta v0.1): https://comptox.epa.gov/chemexpo/get_data/
# STARTED: 2024-01-08
# written in version: R version 4.2.2 (2022-10-31 ucrt)

#UPDATES Feb 2024, per discussion with ToxTeam: 
#       make construction/bldg materials its own category, not consumer products or industrial
#       pesticides include sanitizing products & bug spray
#       move non-FDA regulated OTC products like supplements, Airborne, & lozenges to pharma
#       sunscreen is both pharma and consumer

library(tidyverse)
library(readxl)
library(stringr)


# Assign the folder where the R script lives to working directory
workingdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workingdir)

options(stringsAsFactors = FALSE)
options(scipen=999)


## Combined MC+E/P+ER list of BC-relevant chemicals created in 3_BCrelevant_chem_effects.R
BCrelList <- read.csv("outputs/BCRelList.csv") 

BCrelchems <- BCrelList %>% 
  select(CASRN:preferred_name)


## glossary matching CASRNs to DTXSIDs and chem names
# downloaded version 2021r1 on Jan 19, 2022 from CompTox Dashboard 
# https://epa.figshare.com/articles/dataset/Chemistry_Dashboard_Data_DSSTox_Identifiers_Mapped_to_CAS_Numbers_and_Names/5588566
chemids <- read.csv("./inputs/DSSTox_Identifiers_and_CASRN_2021r1.csv") %>%
  rename(CASRN = casrn, DTXSID = dtxsid, preferred_name = preferredName) %>% 
  select(CASRN, DTXSID, preferred_name)%>% 
  # some characters show up weird for 1,2-benzenediol
  mutate(preferred_name = ifelse(CASRN == "120-80-9", "1,2-Benzenediol", preferred_name))


# need to make chemids glossary all lowercase for matching to FDA drug database
lowercase_chemids <- chemids %>% 
  mutate(preferred_name = tolower(preferred_name))



#### ChemExpo List Presence #### 
# Downloaded List Presence Data from https://comptox.epa.gov/chemexpo/get_data/ Sep 28, 2023
listpres_orig <- read.csv("inputs/ChemExpo_bulk_list_presence_chemicals_20230928.csv")


#make a dataframe using the helpful bits from listpres_problems
listpres <- listpres_orig %>% 
  mutate(preferred_name = ifelse(Curated.Chemical.Name == "", Raw.Chemical.Name, Curated.Chemical.Name),
         #some chems are missing CASRNs/DTXSIDs, but their names are all over the place with capitalization
         #    convert to lowercase so we can fill in blanks w/ lowercase_chemids 
         preferred_name = tolower(preferred_name)) %>% 
  left_join(lowercase_chemids, by = "preferred_name") %>% 
  mutate(CASRN = case_when(is.na(CASRN) & Curated.CAS == "" ~ Raw.CAS,
                           is.na(CASRN) ~ Curated.CAS,
                           TRUE ~ CASRN)) %>% 
  left_join(chemids, by = "CASRN") %>% 
  #clean things up. Note that if something has a CASRN but nothing for DTXSID or preferred_name.y, 
  #    that means EPA doesn't have identifiers for it
  #Note also that there are thousands of things with names but not CASRNs, and those also 
  #    can't be matched to EPA identifiers (all names are problematic for different reasons)
  select(CASRN, DTXSID, preferred_name.y, Reported.Function, Harmonized.Function, 
         Keyword.Set, Organization, Data.Document.Title, Data.Document.Subtitle) %>% 
  filter(CASRN != "" & !is.na(preferred_name.y)) %>% 
  # Get rid of nondetects and banned chems - Kristin Isaacs confirmed that 
  #     "nondetect" keyword means the chemical was not present in the media listed
  # also delete illegal, withdrawn, and veterinary drugs
  filter(!str_detect(Keyword.Set, "nondetect|banned|prohibited") & 
           !str_detect(Data.Document.Title, "Withdrawn|Illicit|EPARs") &
           # exclude things measured in UK, Belgian, & Canadian drinking water. 
           #    All other international stuff could be exported
           !str_detect(Organization, "Inspectorate|Belgaqua|Health Canada"))   


#narrow down to just BCRCs in listpres to make categorization easier/more focused
listpres_bcrcsonly <- inner_join(BCrelchems, listpres, by = c("CASRN", "DTXSID"))


#I did a bunch of iteration to develop categorization schema and checked validity w/ Kristin Isaacs from EPA
bcrc_listpres_cats <- listpres_bcrcsonly %>% 
  mutate(consumer = ifelse(
    str_detect(Keyword.Set, "consumer|Personal|[Tt]oy|children|care|Arts") |
      str_detect(Keyword.Set, "appliance|Furniture|SCIL|[Hh]ome|smoking|Batteries") |
      str_detect(Keyword.Set, "Packaging|[Vv]ehicle|Cleaning") |
      
      str_detect(Data.Document.Title, "[Hh]air [Dd]ye|[Tt]extile|Cigarette|Hobby|toy"), 
    "list_consumer", NA), 
    
    diet = ifelse(
      str_detect(Keyword.Set, "animal_products|baby_food|general_foods|Food contact") |
        str_detect(Keyword.Set, "additives food|food_additive|EAFUS|dairy|drinking") |
        str_detect(Keyword.Set, "grain|fruits|legume|nuts|CEDI") |
        
        str_detect(Harmonized.Function, "[Ff]lavo"), 
      "list_diet", NA),
    
    #adding construction & building material category, removed "building" keyword from consumer
    construction = ifelse(
      str_detect(Keyword.Set, "Construction"), 
      "list_construction", NA), 
    
    industrial = ifelse(
      str_detect(Keyword.Set, "fracking|Raw|plastic|fuel") |
        str_detect(Reported.Function, "Industrial|plastici[sz]er|Intermediate"), 
      "list_industrial", NA), 
    
    pesticide = ifelse(str_detect(Keyword.Set, "[Pp]esticide") |
                         str_detect(Reported.Function, "icide"), "list_pest", NA),
    
    pharma = ifelse(
      (str_detect(Keyword.Set, "[Pp]harma") &
         #exclude drugs used from diff countries, since this paper is about the US and 
         #    any international drugs would have to go through FDA approval
         !str_detect(Keyword.Set, "Canada|Europe")) |
        str_detect(Harmonized.Function, "[Pp]harma"), 
      "list_pharma", NA), 
    
    water = case_when(str_detect(Keyword.Set, "drinking") ~ "list_drinking_water", 
                      str_detect(Keyword.Set, "surface_water") ~ "list_surface_water", 
                      str_detect(Keyword.Set, "ground_water") ~ "list_groundwater", 
                      str_detect(Keyword.Set, "wastewater") ~ "list_wastewater",
                      TRUE ~ NA)) %>% 
  distinct() 

#save this df so we don't have to re-create it every time we want to look at specific listings
write_csv(bcrc_listpres_cats, "outputs/bcrc_listpres.csv")

bcrc_listpres_cats_clean <- bcrc_listpres_cats %>% 
  select(CASRN:preferred_name, consumer:water) %>% 
  pivot_longer(cols = c(consumer:water), names_to = "expsource", values_drop_na = TRUE) %>% 
  select(-expsource) %>% 
  distinct() %>% 
  pivot_wider(names_from = value, values_from = value) 




#### ChemExpo Functional Use data ####
# Downloaded Function Data from https://comptox.epa.gov/chemexpo/get_data/ Oct 12, 2023
functional_orig <- read.csv("inputs/ChemExpo_bulk_functional_uses_20231012.csv") 

functional <- functional_orig %>% 
  mutate(preferred_name = ifelse(Curated.Chemical.Name == "", Raw.Chemical.Name, Curated.Chemical.Name),
         #some chems are missing CASRNs/DTXSIDs, but their names are all over the place with capitalization
         #    convert to lowercase so we can fill in some blanks w/ lowercase_chemids 
         preferred_name = tolower(preferred_name)) %>% 
  left_join(lowercase_chemids, by = "preferred_name") %>% 
  mutate(CASRN = case_when(is.na(CASRN) & Curated.CAS == "" ~ Raw.CAS,
                           is.na(CASRN) ~ Curated.CAS,
                           TRUE ~ CASRN)) %>% 
  left_join(chemids, by = "CASRN") %>% 
  #filter out things w/o EPA-recognized identifiers or that are prohibited/discontinued/not used
  filter(CASRN != "" & !is.na(preferred_name.y) & 
           !str_detect(preferred_name.y, "PROHIB") & 
           !str_detect(Data.Document.Title, "Prohib|do not use|[Dd]iscont|DISCONT")) %>%
  select(CASRN, DTXSID, preferred_name.y, Reported.Functional.Use:Component, 
         Data.Source, Data.Document.Title, Data.Type) 


functional_bcrcsonly <- inner_join(BCrelchems, functional, by = c("CASRN", "DTXSID")) %>% 
  select(-preferred_name.y)


#set up categorization scheme from cats and orgs notes
bcrc_functional_cats <- functional_bcrcsonly %>% 
  mutate(consumer = ifelse(
    #these are companies that sell products that can be purchased by consumers in stores or online
    #   or databases of chems used in products
    str_detect(Data.Source, "Cleaning|Church|Clorox|CosIng|Danish|Environ|IFRA") |
      str_detect(Data.Source, "Johnson|Method|Palmolive|Proctor|Published|RB|Safer")|
      str_detect(Data.Source, "Seventh|Simple|Sun|Turtle|Unilever|Vi-Jon") 
      
      #determined next two lines not necessary
      # | str_detect(Data.Document.Title, "[Cc]hild|Cigarette|[Cc]osmetic|[Tt]oy") |
      # str_detect(Reported.Functional.Use, "[Cc]osmetic")
    , 
    "func_consumer", NA)) %>% 
  
  mutate(diet = ifelse(
    str_detect(Data.Document.Title, "Fl@vis|Food|[Dd]rink|DRINK|Community") |
      str_detect(Harmonized.Functional.Use, "Food") |
      #next line for pesticides used on food
      (str_detect(Data.Source, "USDA") & !str_detect(Data.Document.Title, "Ground|GROUND|Untreated")), 
    "func_diet", NA)) %>% 
  
  mutate(construction = ifelse(
    #Health = Health Product Declaration Collaborative and Healthy Stuff (only HS BCRC entries are PVC carpet backing)
    str_detect(Data.Source, "Health|Tarkett"),
    "func_construction", NA)) %>% 
  
  mutate(industrial = ifelse(
    str_detect(Data.Source, "Chemours|FracFocus|Graham|Kimball|Mayzo|SpecialChem|Sprayon")|
      str_detect(Reported.Functional.Use, "Industrial"), 
    "func_industrial", NA)) %>% 
  
  mutate(pesticide = ifelse(
    str_detect(Harmonized.Functional.Use, "Biocide") |
      str_detect(Reported.Functional.Use, "icide|insect") |
      str_detect(Data.Document.Title, "icide|Raid|OFF!|DISINFECT|[Dd]isinfect")|
      str_detect(Data.Document.Title, "SANITIZ|[Ss]anitiz"),
    "func_pest", NA)) %>% 
  
  mutate(pharma = ifelse(
    str_detect(Harmonized.Functional.Use, "Pharma") |
      str_detect(Reported.Functional.Use, "Dietary Supplement") |
      #these OTC supplements are not FDA-approved
      str_detect(Data.Document.Title, "AIRBORNE|CEPACOL|DIGESTIVE|SCHIFF") |
      #these OTC drugs ARE FDA-regulated
      str_detect(Data.Document.Title, "MUCINEX|DELSYM|STREPSILS|LANACANE"),
    "func_pharma", NA)) %>% 
  
  mutate(water = case_when(
    str_detect(Data.Document.Title, "[Ss]urface [Ww]ater|[Rr]iver|[Ss]tream|Lake") ~ "func_surface_water",
    str_detect(Data.Document.Title, "[Gg]round|GROUND|Aquifer") ~ "func_groundwater", 
    str_detect(Data.Document.Title, "[Dd]rink|DRINK|Community") ~ "func_drinking_water", 
    TRUE ~ NA)) %>% 
  distinct() 

#save this df so we don't have to re-create it every time we want to look at specific listings
write_csv(bcrc_functional_cats, "outputs/bcrc_functional.csv")


bcrc_functional_cats_clean <- bcrc_functional_cats %>% 
  select(CASRN:preferred_name, consumer:water) %>%
  pivot_longer(cols = c(consumer:water), names_to = "expsource", values_drop_na = TRUE) %>% 
  select(-expsource) %>% 
  distinct() %>% 
  pivot_wider(names_from = value, values_from = value)
  


#### ChemExpo Composition data ####
# Downloaded Composition Data from https://comptox.epa.gov/chemexpo/get_data/ Oct 12, 2023
composition_orig <- read.csv("inputs/ChemExpo_bulk_composition_chemicals_20231012.csv") 

composition <- composition_orig %>% 
  mutate(preferred_name = ifelse(Curated.Chemical.Name == "", Raw.Chemical.Name, Curated.Chemical.Name),
         #some chems are missing CASRNs/DTXSIDs, but their names are all over the place with capitalization
         #    convert to lowercase so we can fill in blanks w/ lowercase_chemids 
         preferred_name = tolower(preferred_name)) %>% 
  left_join(lowercase_chemids, by = "preferred_name") %>% 
  mutate(CASRN = case_when(is.na(CASRN) & Curated.CAS == "" ~ Raw.CAS,
                           is.na(CASRN) ~ Curated.CAS,
                           TRUE ~ CASRN)) %>% 
  left_join(chemids, by = "CASRN") %>%
  select(CASRN, DTXSID, preferred_name.y, Data.Source, Data.Document.Title, Data.Document.Subtitle, 
         Product.Name:PUC.Classification.Method, Component) 


composition_bcrcsonly <- inner_join(BCrelchems, composition, by = c("CASRN", "DTXSID")) %>% 
  #filter out things that are discontinued/not used
  filter(!str_detect(Data.Document.Title, "[Dd]iscont|DISCONT") & 
           #filter out things with "automatic" PUC classification, these were not curated and contain 
           #a bunch of errors - like cosmetics with the word "camera" classified as camera batteries
           #TO DO/CHECK NOTE: Manual Batch has lots of issues too... but automatic seems worse
           !str_detect(PUC.Classification.Method, "Automatic") &
           
           # Laboratory supplies PUC makes lots of categorizations difficult... remove?
           !str_detect(PUC.General.Category, "Laboratory")) %>% 
  select(-c(DTXSID, preferred_name.y)) %>% 
  unique()

#this commented block of code just helped to sort through things
# #start narrowing down with PUCs
# composition_bcrcsonly_PUCs <- composition_bcrcsonly %>% 
#   select(PUC.Kind:PUC.Product.Type) %>% 
#   unique()
# #finding that lots of PUCs are miscoded - like calling "asphalt black" personal care products construction materials 
# #more accurate to base categories off of Data.Source
# compositiondatasources <- composition_bcrcsonly %>% 
#   select(Data.Source) %>% unique()
# #figure out what chems I'm missing without pucs
# comp_nopucs <- composition_bcrcsonly %>% 
#   filter(PUC.General.Category == "") %>% 
#   select(preferred_name) %>% unique()

bcrc_comp_cats <- composition_bcrcsonly %>% 
  
    mutate(consumer = ifelse(
    #these companies all make consumer products. 
    str_detect(Data.Source, "Church|Clorox|Environ|Johnson|Method|Palmolive|Proctor") |
      str_detect(Data.Source, "RB|Seventh|Simple|Sun|Turtle|Unilever|Vi-Jon") |
      
      #set conditions for these companies to only get consumer-grade products (others are industrial-grade)
      (str_detect(Data.Source, "Gojo") & PUC.Kind == "Formulation") |
      (str_detect(Data.Source, "Zep") & str_detect(PUC.General.Category, "Vehicle|Home|house")) |
      
      #ID PUCs specific to consumer products
      str_detect(PUC.General.Category, "Arts|house|Home|Furniture|Personal|Sports|yard") |
      #Pet PUC only works if DataSource is not SIRI, those are animal agriculture & euthanasia
      (str_detect(PUC.General.Category, "Pet") & Data.Source != "SIRI") |
      
      #car in product family for boat & auto care, carpets, and body care 
      str_detect(PUC.Product.Family, "car|computer|ink") |
      
      #these products are home-use pesticides or bug sprays
      str_detect(Product.Name, "scotts|[Ff]lea") |
      #roach killers
      (str_detect(Product.Name, "[Rr]oach") & !str_detect(Product.Name, "approach")),
    "comp_consumer", NA)) %>% 
  
  mutate(construction = ifelse(
    #Some personal care products from Environ. Working Group miscoded as PUC.Gen.Cat Construction, otherwise ok
    (str_detect(PUC.General.Category, "Construction") & !str_detect(Data.Source, "Environmental")) |
      
      #many "Home maintenance" things are construction materials that can be bought at hardware stores
      #but many personal care products are miscoded due to product name keywords, so exclude EWG entries
      (str_detect(PUC.General.Category, "Home maintenance") & 
         !str_detect(Data.Source, "Environmental") & 
         str_detect(Product.Family, "paint|caulk|adhesive|insulation|finish|concrete")) |
      
      str_detect(Data.Source, "Health Product|Tarkett") |
      #exclude Declare/Living Future furniture entries, all remaining are construction 
      (str_detect(Data.Source, "Declare") & !str_detect(PUC.General.Category, "Furniture")) , 
    "comp_construction", NA)) %>% 
  
  mutate(industrial = ifelse(
    #these Data sources are chemical or industrial manufacturers or industrial databases
    #   some SIRI entries may seem odd (eg. drugs); SIRI encompasses any type of manufacturing
    str_detect(Data.Source, "3M|British|Chemours|CRC|Kimball|Mayzo|McKinley") |
      str_detect(Data.Source, "Safety|Schaeffer|SIRI|Sprayon|Stepan|Williams") |
      
      #these PUCs capture anything not included by sources above
      str_detect(PUC.General.Category, "[Ii]ndustrial|Raw") | #raw = raw materials
      str_detect(PUC.Product.Family, "engine|insulation")|
      str_detect(PUC.Product.Type, "[Ii]ndustrial") |
      #filter out industrial-sized personal care product entries from EWG
      (str_detect(Product.Name, "[Ii]ndustrial|INDUSTRIAL") & !str_detect(Data.Source, "Environmental")),
    "comp_industrial", NA)) %>%
  
  mutate(pesticide = ifelse(
    str_detect(Data.Source, "ECRI") | #all ECRI entries are pesticides
      str_detect(PUC.General.Category, "Pesticides") |
      str_detect(PUC.Product.Family, "icide") |
      str_detect(PUC.Product.Type, "icide|disinfectant|sanitizer|antiseptic") |
      #next terms are for biocide products not coded as such in PUCs, excluding lab standards 
      (str_detect(Product.Name, "icid") & !str_detect(Product.Name, "standard")) |
      str_detect(Product.Name, "[Ff]ung|[Dd]isinfect|DISINFECT|asept"),
    "comp_pest", NA)) %>% 
  
  mutate(pharma = ifelse(
    # Note, Cannot categorize by PUC.General.Category = Medical/dental. Mostly things used in dentistry
    #     (e.g. orthodontia, making tooth molds) and lots of errors (like gasoline and body oil) 
    (str_detect(PUC.General.Category, "drug") & 
       #Exclude things miscoded with pharmaceutical PUC.Product.Family 
       #1,2-dichloroethane, BPA, BPAF, trifluralin, benfluralin, fenarimol, 1-methoxy-2-propanol, 
       #    N-isopropyl-N'-phenyl-p-phenylenediamine, oryzalin, and styrene are 
       #    NOT used in pharmaceuticals
       !str_detect(preferred_name, "1,2-Dichloroethane|Bisphenol|fluralin|Fenarimol|1-Methoxy|phenylene|Oryzalin|Styrene")), 
    "comp_pharma", NA)) %>% 
  distinct()

#save this df so we don't have to re-create it every time we want to look at specific listings
write_csv(bcrc_comp_cats, "outputs/bcrc_composition.csv")


bcrc_comp_cats_clean <- bcrc_comp_cats %>% 
  select(CASRN:preferred_name, consumer:pharma) %>%
  pivot_longer(cols = c(consumer:pharma), names_to = "expsource", values_drop_na = TRUE) %>% 
  select(-expsource) %>% 
  distinct() %>% 
  pivot_wider(names_from = value, values_from = value)


#### combine final clean versions of list presence, functional use, and composition ####
# into a single file with chem identifiers and the 6 categories we binned things into 

#make a list of the three data frames: list presence, functional use, and composition of BCRCs
combo_list <- list(bcrc_listpres_cats_clean, bcrc_functional_cats_clean, bcrc_comp_cats_clean)
#merge together via full_join and reduce
combo <- combo_list %>% reduce(full_join, by = "CASRN") %>% 
  mutate(DTXSID = coalesce(DTXSID.x, DTXSID.y), 
         preferred_name = coalesce(preferred_name.x, preferred_name.y, preferred_name))%>% 

  #create single column for each type of exposure source 
  # mutate(consumer = ifelse((!is.na(list_consumer) | !is.na(func_consumer) | !is.na(comp_consumer)),
  #                          "CPDat", NA),
  #        diet = ifelse((!is.na(list_diet) | !is.na(func_diet) | !is.na(comp_diet)),
  #                      "CPDat", NA),
  #        industrial = ifelse((!is.na(list_industrial) | !is.na(func_industrial) | !is.na(comp_industrial)),
  #                            "CPDat", NA),
  #        pesticide = ifelse((!is.na(list_pest) | !is.na(func_pest) | !is.na(comp_pest)),
  #                           "CPDat", NA), 
  #        pharma = ifelse((!is.na(list_pharma) | !is.na(func_pharma) | !is.na(comp_pharma)),
  #                        "CPDat", NA)) %>% 
  unite(consumer, list_consumer, func_consumer, comp_consumer, sep = ", ", na.rm = TRUE) %>% 
  unite(diet, list_diet, func_diet, sep = ", ", na.rm = TRUE) %>% 
  unite(industrial, list_industrial, func_industrial, comp_industrial, sep = ", ", na.rm = TRUE) %>% 
  unite(pest, list_pest, func_pest, comp_pest, sep = ", ", na.rm = TRUE) %>% 
  unite(pharma, list_pharma, func_pharma, comp_pharma, sep = ", ", na.rm = TRUE) %>% 
  unite(construction, list_construction, func_construction, comp_construction, sep = ", ", na.rm = TRUE) %>% 
  
  mutate(drinking_water = ifelse(!is.na(list_drinking_water) | !is.na(func_drinking_water), "drinking_water", NA),
         surface_water = ifelse(!is.na(list_surface_water) | !is.na(func_surface_water), "surface_water", NA),
         groundwater = ifelse(!is.na(list_groundwater) | !is.na(func_groundwater), "groundwater", NA),
         wastewater = ifelse(!is.na(list_wastewater), "wastewater", NA)) %>% 
  unite(water, drinking_water:wastewater, sep = ", ", na.rm = TRUE) %>% 
  mutate(water = ifelse(water == "", NA, water)) %>% 
  
  select(CASRN, DTXSID, preferred_name, consumer, diet, construction, pest, industrial, pharma, water) %>% 
  distinct()


write_csv(combo, "outputs/CPDatcombined_list_func_comp.csv")

