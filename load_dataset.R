rm(list=ls())
library(readxl)
library(lubridate)
library(tidyr)
library(dplyr)

## Define the reference date for week number 1
reference_date <- as.Date("2002-01-01", format = "%Y-%m-%d") #as.Date("2002-07-20", format = "%Y-%m-%d")

vec_weeks <- seq.Date(from = reference_date, to = as.Date("2022-10-26"), by = "weeks")
vec_years <- format(vec_weeks - 19*7, "%Y")
head(data.frame(wks = vec_weeks, yrs = vec_years), 80)
names(vec_years) <- 1:length(vec_years)

# Convert weeks to years
week_to_year <- function(week) {
  if (week > 0 && week <= length(vec_years)) {
    return(vec_years[as.character(week)])
  }
  return(NA)
}

# Identify the weeks corresponding to January 1st of each year
january_first_dates <- seq.Date(from = as.Date("2003-05-06"), to = as.Date("2022-05-06"), by = "year")
weeks_january_first <- sapply(january_first_dates, function(date){
  which(vec_weeks <= date & (vec_weeks + 6) >= date)
})

# Check the identified weeks to ensure correctness
print(vec_weeks[weeks_january_first])

#######################
### DATA HORSES     ###
#######################

## Load horse data
tab.cv=as.data.frame(read_excel("./BDD_cv_propre_2.xlsx"))
source("./fix_data_horses.R")
tab.cv2=tab.cv
tab.cv2$`DATE PRELEVEMENT`=as.Date(tab.cv2$`DATE PRELEVEMENT`)

tab.cv2$wk_y= floor(as.numeric(difftime(tab.cv2$`DATE PRELEVEMENT`, reference_date), units="weeks")+1 )

tab.cv2=dplyr::filter(tab.cv2,`RESULTAT  CONF`%in% c("Détecté","Non détecté","Douteux"))
tab.cv2$`RESULTAT  CONF`=ifelse(tab.cv2$`RESULTAT  CONF`=="Détecté",1,0)
tab.coord.cv <- unique(tab.cv2[,c("CLUB","Lon","Lat")])

# Calcul des intervalles en années et segmentation des groupes
tab_inter_an_chev <- tab.cv2 %>%
  # arrange(NOM, `DATE PRELEVEMENT`) %>%
  group_by(NOM) %>%
  mutate(
    Year = year(`DATE PRELEVEMENT`),
    Year_Diff = Year - lag(Year, default=Year[1]),  # Différence en années entre les prélèvements
    Group_ID = cumsum(Year_Diff > 5 | is.na(Year_Diff))  # Nouveau groupe si l'intervalle est > 1 an
  ) %>%
  ungroup()

filtered_data_chev <- tab_inter_an_chev %>%
  group_by(NOM, Group_ID) %>%
  mutate(
    Count = n(),  # Total de prélèvements par groupe
    Echantillon_Num = row_number()  # Numéro séquentiel de chaque prélèvement dans le groupe
  ) %>%
  filter(Count > 1) %>%
  ungroup()

# Pivoter les données pour une visualisation large des résultats sérologiques par échantillon 
wide_data_filtrer_chev <- filtered_data_chev %>%
  unite("ID", NOM, Group_ID, sep = "_") %>%
  pivot_wider(
    names_from = Echantillon_Num,
    values_from = `RESULTAT  CONF`,
    id_cols = ID
  )

# Pivoter les données pour une visualisation large des semaines d'échantillonnage 
wide_data_filtrer_chev_week <- filtered_data_chev %>%
  unite("ID", NOM, Group_ID, sep = "_") %>%
  pivot_wider(
    names_from = Echantillon_Num,
    values_from = wk_y,
    id_cols = ID
  )

####remplacer les NA des 2 matrices par des -1
wide_data_filtrer_chev[is.na(wide_data_filtrer_chev)] <- - 1
wide_data_filtrer_chev_week[is.na(wide_data_filtrer_chev_week)] <- - 1

wide_data_filtrer_chev1=wide_data_filtrer_chev
wide_data_filtrer_chev_week1=wide_data_filtrer_chev_week

# Conversion du tibble en data frame
wide_data_filtrer_chev1 <- as.data.frame(wide_data_filtrer_chev1)
wide_data_filtrer_chev_week1 <- as.data.frame(wide_data_filtrer_chev_week1)

# Assignation des noms de lignes
rownames(wide_data_filtrer_chev1) <- wide_data_filtrer_chev1$ID
wide_data_filtrer_chev1 = dplyr::select(wide_data_filtrer_chev1, -ID)
rownames(wide_data_filtrer_chev_week1) <- wide_data_filtrer_chev_week1$ID
wide_data_filtrer_chev_week1 = dplyr::select(wide_data_filtrer_chev_week1, -ID)

## Création du vecteur avec le nombre d'échantillon par cheval
last_ech_chev=rowSums(wide_data_filtrer_chev1 != -1)

## Création de la liste de vecteur avec en entré les années et pour chaque année le nombre de chevaux 

# Creating a list to store the horse names per year
liste_annees_chevaux2 <- vector("list", length(unique(vec_years)))
names(liste_annees_chevaux2) <- unique(vec_years)

# Loop through each row of the weekly sampling matrix
for (i in 1:nrow(wide_data_filtrer_chev_week)) {
  # Extract the horse's name
  cheval <- wide_data_filtrer_chev_week[i, 1]
  
  # Loop through the weeks (columns) for the current horse
  for (j in 2:ncol(wide_data_filtrer_chev_week)) {
    # Extract the week number
    semaine <- wide_data_filtrer_chev_week[i, j]
    
    # If the week number is valid (not -1)
    if (semaine != -1) {
      # Find the corresponding year
      annee <- vec_years[as.numeric(semaine)]
      
      # If the list for this year doesn't exist, create it
      if (is.null(liste_annees_chevaux2[[as.character(annee)]])) {
        liste_annees_chevaux2[[as.character(annee)]] <- c()
      }
      
      # Add the horse's name to the list for the corresponding year
      liste_annees_chevaux2[[as.character(annee)]] <- c(liste_annees_chevaux2[[as.character(annee)]], cheval)
    }
  }
}

## Create a data frame to store the pairs of samples
paires_chv_df <- data.frame(pair_number = NA,
                        semaine_1 = NA,
                        semaine_2 = NA)

paires_etat_serologique_chv <- data.frame(semaine_1 = integer(),
                                      semaine_2 = integer(),
                                      etat_serologique_1 = character(),
                                      etat_serologique_2 = character())

## Create table where 1 row = 1 pair. First column is individual number, second column is number of the pair in this individual
tab_indiv_pair_chev <- data.frame(indiv = integer(),
                                  pair_in_indiv = integer())

# Initialize pair number counter
pair_counter_chv <- 0

# Loop over each row in the dataframe
for (i in 1:nrow(wide_data_filtrer_chev_week1)) {
  # Get the sample weeks (excluding the name column)
  semaine_i <- as.numeric(wide_data_filtrer_chev_week1[i, ])
  etat_serologique_i <- wide_data_filtrer_chev1[i, ]
  
  # Filter out invalid weeks
  semaine_i <- semaine_i[semaine_i != -1]
  etat_serologique_i <- etat_serologique_i[etat_serologique_i != -1]
  
  # Iterate over pairs of consecutive weeks
  for (j in 1:(length(semaine_i) - 1)) {
    
    semaine1 <- semaine_i[j]
    semaine2 <- semaine_i[j + 1]
    
    # Ajoutez les informations de la paire de semaines et d'états sérologiques au data frame
    paires_etat_serologique_chv <- rbind(paires_etat_serologique_chv,
                                         data.frame(semaine_1=semaine1, semaine_2=semaine2, etat_serologique_1=etat_serologique_i[j], etat_serologique_2=etat_serologique_i[j+1]))
    
    tab_indiv_pair_chev <- rbind(tab_indiv_pair_chev,
                                 data.frame(indiv = i, pair_in_indiv = j))
    
    annee1 <- week_to_year(semaine1)
    annee2 <- week_to_year(semaine2)
    
    # Create a vector for this pair
    paire_j <- c(pair_number=pair_counter_chv, semaine_1=semaine1, semaine_2=semaine2)
    
    # Identify the weeks of year change
    semaines_changement_annee <- c()
    if (annee1 != annee2) {
      # Check if any of the weeks are January 1st weeks
      semaines_changement_annee <- weeks_january_first[weeks_january_first >= semaine1-1 & weeks_january_first < semaine2]
      
      # Add columns for each year change week
      for (k in 1:length(semaines_changement_annee)) {
        paire_j[paste0("changement_annee_", k)] <- semaines_changement_annee[k]
      }
    }
    
    # Ensure all rows have the same columns by adding NA to missing columns
    missing_cols <- setdiff(names(paires_chv_df), names(paire_j))
    for (col in missing_cols) {
      paire_j[[col]] <- NA
    }
    
    extra_cols <- setdiff(names(paire_j), names(paires_chv_df))
    for (col in extra_cols) {
      paires_chv_df[[col]] <- NA
    }
    
    # Bind the new row to the dataframe
    paires_chv_df <- rbind(paires_chv_df, paire_j)
    
    # Increment pair counter
    pair_counter_chv <- pair_counter_chv + 1
  }
}
paires_chv_df <- paires_chv_df[-1,]

###remplacer les NA par la dernière valeur de chaque colonne###
paires_chv_df1 = paires_chv_df
for(row_i in 1:nrow(paires_chv_df1)){
  max_no_NA <- max(which(!is.na(paires_chv_df1[row_i,])))
  paires_chv_df1[row_i, is.na(paires_chv_df1[row_i,])] <- paires_chv_df1[row_i,max_no_NA]
}

# Order columns
rownames(paires_chv_df1) <- paires_chv_df1$pair_number
paires_chv_df1 <- paires_chv_df1[, c("semaine_1", sort(setdiff(names(paires_chv_df), c("pair_number","semaine_1","semaine_2"))), "semaine_2")]

# Calculer le nombre total de paires
nb_total_paires <- sum(apply(wide_data_filtrer_chev_week1, 1, function(x) sum(x != -1) - 1))

paires_etat_serologique_chv1=paires_etat_serologique_chv

#######################
### DATA CHICKENS   ###
#######################

tab.pl=as.data.frame(read_excel("./BDD_pl_propre_2.xlsx"))
source("./fix_data_chickens.R")
tab.pl$`DATE PRELEVEMENT`=as.Date(tab.pl$`DATE PRELEVEMENT`)

tab.pl2=tab.pl
tab.pl2$anne<-format(as.Date(tab.pl2$`DATE PRELEVEMENT`,
                             format="%d/%m/%Y"),"%Y")

# Convert DATE_PRELEVEMENT to POSIXct format
tab.pl2$`DATE PRELEVEMENT` <- as.Date(tab.pl2$`DATE PRELEVEMENT`, format = "%d/%m/%Y")

# Calculate the week number relative to the reference date
tab.pl2$wk_y <- floor(as.numeric(difftime(tab.pl2$`DATE PRELEVEMENT`, reference_date), units = "weeks") + 1)

# table(tab.pl2$`RESULTAT  CONF`)
tab.pl2 = dplyr::filter(tab.pl2,`Type d'Enquete`%in% c("Sentinelle"),
                        ! anne %in% c("2019","2020","2021","2022"),
                        `RESULTAT  CONF`%in% c("Détecté","Non détecté","DOUTEUX"))
tab.pl2$`RESULTAT  CONF`=ifelse(tab.pl2$`RESULTAT  CONF`=="Détecté",1,0)
tab.coord.pl <- unique(tab.pl2[,c("NOM Eleveur","Lon","Lat")])

###comme y a des poulets qui sont testé 2 fois de suite la même semaine donc on prend un seul résultat par semaine
#pour éviter qu'il soit en double dans les tables wide_data_poul
tab.pl2 = tab.pl2 %>%
  group_by(`DATE PRELEVEMENT`,DESCRIPTION,`NOM Eleveur`,anne,wk_y) %>%
  summarise(resultat=max(`RESULTAT  CONF`)) %>%
  mutate(des_elveur=paste(DESCRIPTION, `NOM Eleveur`, sep="_"))

# Enleve poulet anormal (trop longue longevite)
tab.pl2 <- filter(tab.pl2, des_elveur != "B150_Mr Lada")

# Calcul des intervalles en années et segmentation des groupes
tab_inter_an_poul <- tab.pl2 %>%
  arrange(des_elveur, `DATE PRELEVEMENT`) %>%
  group_by(des_elveur) %>%
  mutate(
    Year = year(`DATE PRELEVEMENT`),
    Year_Diff = Year - lag(Year, default=Year[1]),  # Différence en années entre les prélèvements
    Group_ID = cumsum(Year_Diff > 5 | is.na(Year_Diff))  # Nouveau groupe si l'intervalle est > 1 an
  ) %>%
  ungroup()


filtered_data_poul <- tab_inter_an_poul %>%
  group_by(des_elveur, Group_ID) %>%
  mutate(
    Count = n(),  # Total de prélèvements par groupe
    Echantillon_Num = row_number()  # Numéro séquentiel de chaque prélèvement dans le groupe
  ) %>%
  filter(Count > 1) %>%
  ungroup()



# Pivoter les données pour une visualisation large des résultats sérologiques par échantillon 
wide_data_filtrer_poul <- filtered_data_poul %>%
  unite("ID", des_elveur, Group_ID, sep = "_") %>%
  pivot_wider(
    names_from = Echantillon_Num,
    values_from = resultat,
    id_cols = ID
  )

# Pivoter les données pour une visualisation large des semaines d'échantillonnage 
wide_data_filtrer_poul_week <- filtered_data_poul %>%
  unite("ID", des_elveur, Group_ID, sep = "_") %>%
  pivot_wider(
    names_from = Echantillon_Num,
    values_from = wk_y,
    id_cols = ID
  )


####remplacer les NA des 2 matrices par des -1
wide_data_filtrer_poul[is.na(wide_data_filtrer_poul)] <- - 1
wide_data_filtrer_poul_week[is.na(wide_data_filtrer_poul_week)] <- - 1


wide_data_filtrer_poul1=wide_data_filtrer_poul
wide_data_filtrer_poul_week1=wide_data_filtrer_poul_week

# Conversion du tibble en data frame
wide_data_filtrer_poul1 <- as.data.frame(wide_data_filtrer_poul1)
wide_data_filtrer_poul_week1 <- as.data.frame(wide_data_filtrer_poul_week1)

# Assignation des noms de lignes
rownames(wide_data_filtrer_poul1) <- wide_data_filtrer_poul1$ID
wide_data_filtrer_poul1 = dplyr::select(wide_data_filtrer_poul1, -ID)
rownames(wide_data_filtrer_poul_week1) <- wide_data_filtrer_poul_week1$ID
wide_data_filtrer_poul_week1 = dplyr::select(wide_data_filtrer_poul_week1, -ID)

## Création du vecteur avec le nombre d'échantillon par poulet
last_ech_poul = rowSums(wide_data_filtrer_poul1 != -1)

# Creating a list to store the chicken IDs per year
liste_annees_poulet2 <- vector("list", length(unique(vec_years)))
names(liste_annees_poulet2) <- unique(vec_years)

# Loop through each row of the weekly sampling matrix
for (i in 1:nrow(wide_data_filtrer_poul_week)) {
  # Extract the horse's name
  poulet <- wide_data_filtrer_poul_week[i, 1]
  
  # Loop through the weeks (columns) for the current horse
  for (j in 2:ncol(wide_data_filtrer_poul_week)) {
    # Extract the week number
    semaine <- wide_data_filtrer_poul_week[i, j]
    
    # If the week number is valid (not -1)
    if (semaine != -1) {
      # Find the corresponding year
      annee <- vec_years[as.numeric(semaine)]
      
      # If the list for this year doesn't exist, create it
      if (is.null(liste_annees_poulet2[[as.character(annee)]])) {
        liste_annees_poulet2[[as.character(annee)]] <- c()
      }
      
      # Add the horse's name to the list for the corresponding year
      liste_annees_poulet2[[as.character(annee)]] <- c(liste_annees_poulet2[[as.character(annee)]], poulet)
    }
  }
}

# Calculate the number of unique chickens per year
unique_poulet_per_year <- sapply(liste_annees_poulet2, function(x) length(unique(x)))
print(unique_poulet_per_year)

# Identify years with less than 3 unique chickens
years_with_few_poulet <- names(unique_poulet_per_year[unique_poulet_per_year < 3])
print(years_with_few_poulet)

## Create a data frame to store the pairs of samples
paires_pl_df <- data.frame(pair_number = NA,
                            semaine_1 = NA,
                            semaine_2 = NA)

paires_etat_serologique_pl <- data.frame(semaine_1 = integer(),
                                          semaine_2 = integer(),
                                          etat_serologique_1 = character(),
                                          etat_serologique_2 = character())

## Create table where 1 row = 1 pair. First column is individual number, second column is number of the pair in this individual
tab_indiv_pair_poul <- data.frame(indiv = integer(),
                                  pair_in_indiv = integer())

# Initialize pair number counter
pair_counter_pl <- 0

# Loop over each row in the dataframe
for (i in 1:nrow(wide_data_filtrer_poul_week1)) {
  # Get the sample weeks (excluding the name column)
  semaine_i <- as.numeric(wide_data_filtrer_poul_week1[i, ])
  etat_serologique_i <- wide_data_filtrer_poul1[i, ]
  
  # Filter out invalid weeks
  semaine_i <- semaine_i[semaine_i != -1]
  etat_serologique_i <- etat_serologique_i[etat_serologique_i != -1]
  
  # Iterate over pairs of consecutive weeks
  for (j in 1:(length(semaine_i) - 1)) {
    
    semaine1 <- semaine_i[j]
    semaine2 <- semaine_i[j + 1]
    
    # Ajoutez les informations de la paire de semaines et d'états sérologiques au data frame
    paires_etat_serologique_pl <- rbind(paires_etat_serologique_pl,
                                        data.frame(semaine_1=semaine1, semaine_2=semaine2, etat_serologique_1=etat_serologique_i[j], etat_serologique_2=etat_serologique_i[j+1]))
    
    tab_indiv_pair_poul <- rbind(tab_indiv_pair_poul,
                                 data.frame(indiv = i, pair_in_indiv = j))
    
    annee1 <- week_to_year(semaine1)
    annee2 <- week_to_year(semaine2)
    
    # Create a vector for this pair
    paire_j <- c(pair_number=pair_counter_pl, semaine_1=semaine1, semaine_2=semaine2)
    
    # Identify the weeks of year change
    semaines_changement_annee <- c()
    if (annee1 != annee2) {
      # Check if any of the weeks are January 1st weeks
      semaines_changement_annee <- weeks_january_first[weeks_january_first >= semaine1-1 & weeks_january_first < semaine2]
      
      # Add columns for each year change week
      for (k in 1:length(semaines_changement_annee)) {
        paire_j[paste0("changement_annee_", k)] <- semaines_changement_annee[k]
      }
    }
    
    # Ensure all rows have the same columns by adding NA to missing columns
    missing_cols <- setdiff(names(paires_pl_df), names(paire_j))
    for (col in missing_cols) {
      paire_j[[col]] <- NA
    }
    
    extra_cols <- setdiff(names(paire_j), names(paires_pl_df))
    for (col in extra_cols) {
      paires_pl_df[[col]] <- NA
    }
    
    # Bind the new row to the dataframe
    paires_pl_df <- rbind(paires_pl_df, paire_j)
    
    # Increment pair counter
    pair_counter_pl <- pair_counter_pl + 1
  }
}
paires_pl_df <- paires_pl_df[-1,]

###remplacer les NA par la dernière valeur de chaque colonne###
paires_pl_df1 = paires_pl_df
for(row_i in 1:nrow(paires_pl_df1)){
  max_no_NA <- max(which(!is.na(paires_pl_df1[row_i,])))
  paires_pl_df1[row_i, is.na(paires_pl_df1[row_i,])] <- paires_pl_df1[row_i,max_no_NA]
}

# Order columns
rownames(paires_pl_df1) <- paires_pl_df1$pair_number
paires_pl_df1 <- paires_pl_df1[, c("semaine_1", sort(setdiff(names(paires_pl_df), c("pair_number","semaine_1","semaine_2"))), "semaine_2")]

# Calculer le nombre total de paires
nb_total_paires <- sum(apply(wide_data_filtrer_poul_week1, 1, function(x) sum(x != -1) - 1))

paires_etat_serologique_pl1=paires_etat_serologique_pl





###### préparer le vecteur du numéro d'année#####

year_num <- case_when(vec_years == "2001" ~ 0,
                      vec_years == "2002" ~ 1,
                      vec_years == "2003" ~ 2,
                      vec_years == "2004" ~ 3,
                      vec_years == "2005" ~ 4,
                      vec_years == "2006" ~ 5,
                      vec_years == "2007" ~ 6,
                      vec_years == "2008" ~ 7,
                      vec_years == "2009" ~ 8,
                      vec_years == "2010" ~ 9,
                      vec_years == "2011" ~ 10,
                      vec_years == "2012" ~ 11,
                      vec_years == "2013" ~ 12,
                      vec_years == "2014" ~ 13,
                      vec_years == "2015" ~ 14,
                      vec_years == "2016" ~ 15,
                      vec_years == "2017" ~ 16,
                      vec_years == "2018" ~ 17,
                      vec_years == "2019" ~ 18,
                      vec_years == "2020" ~ 19,
                      vec_years == "2021" ~ 20,
                      vec_years == "2022" ~ 21,
                      TRUE ~ NA_integer_)

##### Nombre de paire d'échantillon chez les poulets et chevaux
N_pairs_chev <- nrow(paires_etat_serologique_chv1)
N_pairs_poul <- nrow(paires_etat_serologique_pl1)

