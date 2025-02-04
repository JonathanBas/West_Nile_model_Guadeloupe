library(ggplot2)
library(dplyr)

### Plot des différents modèles:

v_week <- 1:500
week_peak <- 0
prop_min <- 0.2
foi_year_1 <- rep(0.4, 11)
foi_year_2 <- c(0.3, 0.1, 0.2, 0.5, 0.3, 0.3, 0.01, 0.2, 0.5, 0.3, 0.01)
year_num <- 1 + (v_week - week_peak + 26) %/% 52

plot_mod <- rbind(data.frame(mod = "FlatStable",
                             sem = v_week,
                             foi = (foi_year_1[year_num[v_week]]/2)),
                  data.frame(mod = "FlatVary",
                             sem = v_week,
                             foi = (foi_year_2[year_num[v_week]]/2)),
                  data.frame(mod = "SeasoStable",
                             sem = v_week,
                             foi = (foi_year_1[year_num[v_week]]/2) * (1-prop_min) * (cos((v_week - week_peak)*3.142/26) + 1) + prop_min * foi_year_1[year_num[v_week]]),
                  data.frame(mod = "SeasoVary",
                             sem = v_week,
                             foi = (foi_year_2[year_num[v_week]]/2) * (1-prop_min) * (cos((v_week - week_peak)*3.142/26) + 1) + prop_min * foi_year_2[year_num[v_week]]))

pexmod = ggplot(data=plot_mod, aes(x=sem, y=foi)) +
  geom_line(linewidth=1) +
  facet_wrap(~mod, ncol=2, scales='free_x') +
  xlab("Time") + ylab(~paste("Force of infection ", lambda["i"], "(t)")) +
  scale_x_continuous(breaks=c(seq(0, 500, 52)), labels=c(paste0("201",(1:9)),"2020")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5))

grDevices::cairo_pdf("./fig_ex_mod.pdf", width=6, height=4.5)
pexmod
dev.off()

### Figure sampling

tab_samp <- rbind(data.frame(Species = "Horses", date = filtered_data_chev$`DATE PRELEVEMENT`, res = filtered_data_chev$`RESULTAT  CONF`),
                  data.frame(Species = "Chickens", date = filtered_data_poul$`DATE PRELEVEMENT`, res = filtered_data_poul$resultat))

psamp = ggplot(tab_samp, aes(x = date, fill = as.factor(res))) +
  geom_histogram(col="black", bins=75) +
  scale_x_date(name="Sampling dates", date_labels="%m-%Y", date_breaks="year") +
  scale_fill_manual(name="Serological result:", labels=c("0"="Negative", "1"="Positive"), values=c("0"="#1b9e77", "1"="#d95f02")) +
  ylab("Number of samples") +
  facet_wrap(~Species, ncol=1, scales="free_y") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        legend.position=c(0.3, 0.8),
        legend.box.background=element_rect(colour="black", linewidth=1))

pdf(file="./fig_samp_dates.pdf", width=5, height=5)
psamp
dev.off()

### Statistiques descriptives donnees sero

length(unique(filtered_data_chev$NOM))
nrow(filtered_data_chev)
unique(filtered_data_chev$CLUB)
summary(filtered_data_chev$`DATE PRELEVEMENT`)
delais_ech_chv <- as.data.frame(filtered_data_chev %>%
                                  group_by(NOM) %>%
                                  mutate(dates_diff = `DATE PRELEVEMENT` - lag(`DATE PRELEVEMENT`)))
summary(as.numeric(delais_ech_chv$dates_diff))

length(unique(filtered_data_poul$des_elveur))
nrow(filtered_data_poul)
unique(filtered_data_poul$`NOM Eleveur`)
summary(filtered_data_poul$`DATE PRELEVEMENT`)
delais_ech_pl <- as.data.frame(filtered_data_poul %>%
                                  group_by(des_elveur) %>%
                                  mutate(dates_diff = `DATE PRELEVEMENT` - lag(`DATE PRELEVEMENT`)))
summary(as.numeric(delais_ech_pl$dates_diff))

### Tableau états sérologiques consécutifs

table(paires_etat_serologique_chv1[,c("etat_serologique_2", "etat_serologique_1")])
nrow(paires_etat_serologique_chv1)

table(paires_etat_serologique_pl1[,c("etat_serologique_2", "etat_serologique_1")])
nrow(paires_etat_serologique_pl1)

