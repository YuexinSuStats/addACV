library(TCGAbiolinks)
library(dplyr)

###prepare data---------------
clinical<-TCGAbiolinks::GDCquery_clinic(project= "TCGA-BRCA", type = "clinical")

query_biotab <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

GDCdownload(query_biotab)
biotab_data <- GDCprepare(query_biotab)
patient_data <- biotab_data[["clinical_patient_brca"]]

# er_status_by_ihc (ER)
# pr_status_by_ihc (PR)
# her2_status_by_ihc (HER2)

er_data <- patient_data %>%
  dplyr::select(bcr_patient_barcode, er_status_by_ihc) %>%
  dplyr::mutate(
    ER_Status = dplyr::case_when(
      er_status_by_ihc == "Positive" ~ 1,
      er_status_by_ihc == "Negative" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::filter(!is.na(ER_Status))

clinical <- merge(clinical, er_data, by = "bcr_patient_barcode", all.x = TRUE)


# - a preparation
sdata.pre <- data.frame(
  yobs   = ifelse(!is.na(clinical[,'days_to_death']),
                  clinical[,'days_to_death'],
                  clinical[,'days_to_last_follow_up']) / 365,
  delta  = ifelse(!is.na(clinical[,'days_to_death']), 1, 0),
  
  Age    = clinical[,'age_at_diagnosis'] / 365,
  
  Race   = ifelse(clinical[,'race'] == 'black or african american', 'black',
                  ifelse(clinical[,'race'] == 'white', 'white', 'other')),
  
  Gender = ifelse(clinical[,'gender'] == 'female', 'Female',
                  ifelse(clinical[,'gender'] == 'male', 'Male', NA)),
  Laterality = ifelse(clinical[,'laterality'] == 'Left','Left','Right'),
  Stage  = sapply(
    clinical[,'ajcc_pathologic_stage'],
    function(x, pattern = 'Stage X|Stage IV|Stage [I]*') {
      ifelse(grepl(pattern, x),
             regmatches(x, regexpr(pattern, x)),
             NA)
    },
    USE.NAMES = FALSE
  ),
  
  ERstatus = clinical[,'ER_Status'],
  Day = clinical[,'days_to_last_follow_up']
  
)


# - extract covariates and remove undesiable subjects and NA
sdata.TCGA <- na.omit(
  sdata.pre[
    sdata.pre[,'yobs'] > 0
    & sdata.pre[,'Age'] <= 75
    & sdata.pre[,'Gender'] == "Female"
    & sdata.pre[,'Race'] %in% c('white')
    & sdata.pre[,'Laterality'] %in% c('Right')
    & sdata.pre[,'Stage'] %in% c('Stage I','Stage II','Stage III')
    & sdata.pre[,'Day'] > 0 ,
    c('yobs','delta','Age','ERstatus'),
  ]
)

rownames(sdata.TCGA) <- NULL
save(
  sdata.TCGA,
  file = "sdata.TCGA.RData",
  compress = "xz"
)