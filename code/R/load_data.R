
#linelist for descriptive analysis
start_date <- as.Date("2022-05-01")
end_date <- as.Date("2023-04-30") 


#FDOH dataset
fdoh_linelist_raw <- xlsx::read.xlsx("data/final_set/cases for paper.xlsx",
                                                 password="Forrest line list",
                                                 sheetIndex = 1
)

fdoh_linelist <- fdoh_linelist_raw %>%
                    select(
                      case_id=Case..,
                      Onset.Date,
                      County,
                      Imported.Status,
                      Origin,
                      original_serotype=Serotype,
                      age=Age,
                      sex=Gender,
                      race=Race,
                      ethnicity=Ethnicity
                      )%>%
  mutate(travel.status=ifelse(Imported.Status=="Imported",
                              "Travel associated",
                              "Locally acquired"
                              )) %>%
  mutate(case_week=floor_date(Onset.Date,"weeks"))%>%
  filter(Onset.Date>=start_date) %>%
  filter(Onset.Date<=end_date) %>%
  mutate(denvserotype=ifelse(original_serotype %in% c("DENV-1",
                                                      "DENV-2",
                                                      "DENV-3",
                                                      "DENV-4"
                                                      ),
                             original_serotype,"Unknown")) %>%
  mutate(countryoforigin=case_when(
    Origin %in% c("Caribbean",
                  "Central America",
                  "Cuba/Central America",
                  "Unknown"
                  ) ~ "Unknown",
    Origin %in% c("Broward","Collier",
                  "Miami-Dade","Volusia") ~ NA,
    TRUE ~ Origin
  ))


denv3_only <-  fdoh_linelist %>% 
  filter(denvserotype=="DENV-3") 


#arbonet data
arboNet_linelist_raw<- read_excel("data/FL_2023-04-27_analyses.arbonet.data.2010.2023.xlsx") 
arboNet_linelist <- arboNet_linelist_raw%>%
  mutate(onsetdate =as.Date(onsetdate)) %>%
  mutate(case_week=floor_date(onsetdate,"weeks"))%>%
  mutate(age=as.numeric(age)) %>%
  mutate(travel.status2=case_when(
      travel.status == "Travel associated" & countryoforigin=="Cuba"~ "Traveler - Cuba",
      travel.status == "Travel associated"  & countryoforigin!="Cuba" ~ "Traveler - Not Cuba",
      travel.status == "Locally acquired" ~ "Locally acquired"
  ))%>%
  filter(case_week>=start_date) %>%
  filter(case_week<end_date)

arboNet_denv3_only <- arboNet_linelist %>% 
  filter(denvserotype=="DENV-3") 


# WGS data
genome_size <- 10170

#list of sequences
sequences <- read_excel("data/Dengue 3 positives for phylogeographic analysis.xlsx") %>%
  mutate(`Onset Date`=as.Date(`Onset Date`))%>%
  mutate(case_week=floor_date(`Onset Date`,"weeks"))%>%
  mutate(travel.status=case_when(
    `Import Status (US)`== "Acquired in Florida" ~ "Locally acquired",
    `Import Status (US)`== "Acquired outside the US" ~ "Travel associated",
    TRUE~as.character(`Import Status (US)`))) %>%
    mutate(seq_id=`GenBank acc#`)



# genetic distance matrix
pairwise <- read_excel("data/paiwise_distance_FLzoomTree_v3.xlsx",
                       col_names = FALSE
                       ) %>%
            bind_cols(NA)



pairwise_ids_original <- pairwise[,1] %>% unlist() %>% unlist()
pairwise_ids <- str_extract(pairwise_ids_original,"O[A-Z][0-9]{6}") %>% unlist()
pairwise_ids[is.na(pairwise_ids)] <-pairwise_ids_original[is.na(pairwise_ids)]

pairwise_matrix <- pairwise[1:nrow(pairwise),2:(nrow(pairwise)+1)] %>% as.matrix()
colnames(pairwise_matrix) <-NULL 
pairwise_matrix <- pairwise_matrix *genome_size
pairwise_matrix <- Matrix::forceSymmetric(pairwise_matrix,uplo="L") %>% as.matrix() 
diag(pairwise_matrix) <- 0

#data for vimes analysis
vimes_analysis <- sequences %>% filter(seq_id %in% pairwise_ids) %>%
                      arrange(`Onset Date`) %>% filter(travel.status=="Locally acquired") 

# different=setdiff(sequences$seq_id, pairwise_ids)
# miss=sequences[which(sequences$seq_id %in% different),]
# write.csv(miss, "output/missing_cases.csv", row.names=F)
# pairwisediff=setdiff(pairwise_ids, sequences$seq_id)
# write.csv(pairwisediff, "output/pairwise_extra.csv", row.names=F)

D_dates <- dist(vimes_analysis$`Onset Date`)
D_dna <-   pairwise_matrix[match(vimes_analysis$seq_id,pairwise_ids),
                           match(vimes_analysis$seq_id,pairwise_ids)] %>%
              as.dist()

D_all <- vimes::vimes_data(dates = D_dates, dna = D_dna)


# data for transmisson tree_analysis
tt_analysis <- sequences %>% filter(seq_id %in% pairwise_ids) %>%
  arrange(`Onset Date`) 
gendist <- pairwise_matrix[match(tt_analysis$seq_id,pairwise_ids),
                match(tt_analysis$seq_id,pairwise_ids)]
write.csv(tt_analysis, "data/cases_sorted.csv", row.names = F)
write.csv(gendist, "data/gendist_mat.csv", row.names = F)



# census_api_key("33b6456dae2fa1564107e17b476071389575effa",
#                install = TRUE)
# #download census data
# pop20 <- get_decennial(
#       geography = "county",
#       state="Florida",
#       variables = "P1_001N",
#       year = 2020,
#       geometry = TRUE
# )

# write_rds(pop20,"data/census/pop20.rds")

pop20 <- read_rds("data/census/pop20.rds") %>%
  mutate(County =str_remove(NAME,", Florida"))


# zipcodes <- tigris::zctas(cb=TRUE, year=2020)
# 
# miami_zip_shape <-   st_intersection(zipcodes,
#                        filter(pop20,County=="Miami-Dade County"))
# 
# write_rds(miami_zip_shape,"data/miami_zip_shape.rds")

miami_zip_shape <- read_rds("data/miami_zip_shape.rds")


# 
# 
# miami_zip %>% ggplot()+
#     geom_sf()+
#     geom_sf(data=filter(pop20,County=="Miami-Dade County"),
#             alpha=0.5,fill="red"
#             )
# 
# 
# data.frame(area=miami_zip %>% st_area) %>% mutate(side=sqrt(area)) %>%
#       mutate(radius=sqrt(area/pi)) %>%
#       mutate(area=area/(1000)^2) %>%
#       summary()




