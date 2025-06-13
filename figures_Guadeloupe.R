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
  scale_x_continuous(breaks=c(seq(0, 500, 52)), labels=c(paste0("Year ",(1:9)),"Year 10")) +
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
all_trim <- levels(interaction(sep="-", substr(min(tab_samp$date),1,4):substr(max(bdd_mous_aggr$sem_pvmt),1,4), c("Q1","Q2","Q3","Q4")))
all_trim <- expand.grid(Species = unique(tab_samp$Species), year_trim = all_trim)

tab_samp_2 <- as.data.frame(tab_samp %>%
                              mutate(year_trim = paste(substr(date,1,4), quarters(date), sep="-")) %>%
                              group_by(Species, year_trim) %>%
                              summarise(Positive = sum(res),
                                        Negative = sum(res == 0)) %>%
                              merge(y=all_trim, by=c("Species","year_trim"), all=T) %>%
                              pivot_longer(-c(Species,year_trim), values_to="nump_samp", names_to="resul") %>%
                              replace_na(list(nump_samp=0)) %>%
                              mutate(resul = factor(resul, levels=c("Negative","Positive")),
                                     year_trim = factor(year_trim, levels=unique(year_trim))))

psamp = ggplot(tab_samp_2, aes(x=year_trim, y=nump_samp, fill=resul)) +
  geom_bar(col="black", stat="identity", position="stack") +
  # scale_x_date(name="Sampling dates", date_labels="%m-%Y", date_breaks="year") +
  scale_fill_manual(name="Serological result:", values=c("Negative"="#1b9e77", "Positive"="#d95f02")) +
  ylab("Number of samples") + xlab("Sampling dates") +
  facet_wrap(~Species, ncol=1, scales="free_y") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, size=c(8,0), colour=c("gray30","white")),
        legend.position=c(0.3, 0.8),
        legend.box.background=element_rect(colour="black", linewidth=1))

pdf(file="./fig_samp_dates.pdf", width=7, height=5)
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
nb_samp_per_chev <- as.data.frame(filtered_data_chev %>%
                                    group_by(NOM) %>%
                                    summarise(nb_samp = n()))
summary(as.numeric(nb_samp_per_chev$nb_samp))

length(unique(filtered_data_poul$des_elveur))
nrow(filtered_data_poul)
unique(filtered_data_poul$`NOM Eleveur`)
summary(filtered_data_poul$`DATE PRELEVEMENT`)
delais_ech_pl <- as.data.frame(filtered_data_poul %>%
                                  group_by(des_elveur) %>%
                                  mutate(dates_diff = `DATE PRELEVEMENT` - lag(`DATE PRELEVEMENT`)))
summary(as.numeric(delais_ech_pl$dates_diff))
nb_samp_per_poul <- as.data.frame(filtered_data_poul %>%
                                    group_by(des_elveur) %>%
                                    summarise(nb_samp = n()))
summary(as.numeric(nb_samp_per_poul$nb_samp))

pdistrnb1 <- ggplot() +
  geom_bar(aes(x = nb_samp_per_chev$nb_samp), stat="count", col="black", fill="#1b9e77") +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks=seq(0,100,2), limits=c(0,NA)) +
  theme_bw()
pdistrnb2 <- ggplot() +
  geom_bar(aes(x = nb_samp_per_poul$nb_samp), stat="count", col="black", fill="#1b9e77") +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks=seq(0,100,5), limits=c(0,NA)) +
  theme_bw()
pdistrnb <- ggarrange(pdistrnb1, pdistrnb2, align="hv", ncol=2, nrow=1, labels=c("Horses", "Chickens"), label.x=0.4, label.y=0.97)
pdistrnb <- annotate_figure(pdistrnb, left = text_grob("Count", size=11, rot=90))
pdistrnb <- annotate_figure(pdistrnb, bottom = text_grob("Number of samples per individual", size=11))

png(file="./fig_distr_nb_samp_per_indiv.png", width=16, height=8, units="cm", res=300)
pdistrnb
dev.off()

### Tableau états sérologiques consécutifs

table(paires_etat_serologique_chv1[,c("etat_serologique_2", "etat_serologique_1")])
nrow(paires_etat_serologique_chv1)

table(paires_etat_serologique_pl1[,c("etat_serologique_2", "etat_serologique_1")])
nrow(paires_etat_serologique_pl1)

### Cartes

library(sf)
library(terra)
library(ggspatial)

guad = read_sf("./gadm41_GLP_shp/gadm41_GLP_2.shp")

comm_cv_agr <- read.table("./comm_cv_agr.csv", sep=";", header=T)
comm_pl_agr <- read.table("./comm_pl_agr.csv", sep=";", header=T)
comm_mosq_agr <- read.table("./comm_mosq_agr.csv", sep=";", header=T)

guad_plot <- guad %>%
  full_join(comm_cv_agr, join_by(NAME_2==Commune)) %>%
  full_join(comm_pl_agr, join_by(NAME_2==Commune)) %>%
  full_join(comm_mosq_agr, join_by(NAME_2==Commune)) %>%
  replace_na(list(n_samp_cv=0, n_samp_pl=0, n_sit_cv=0, n_sit_pl=0, n_sit_mosq=0))

mapguad <- function(var_plot, plot_tit=NULL, leg_tit=NULL, max_scale=NA, y_lab=T, plot_leg=T){
  p <- ggplot() +
    ggtitle(plot_tit) +
    geom_sf(data=guad_plot, aes(fill = .data[[var_plot]])) +
    scale_fill_gradient(name=leg_tit, low="#fff7bc", high="#d95f02", limits=c(0,max_scale)) +
    coord_sf(crs = crs(guad_plot)) +
    annotation_scale(location="bl", width_hint=0.5)+
    theme_minimal()
  
  if(y_lab){
    p <- p + labs(x="Longitude", y="Latitude")
  }else{
    p <- p + labs(x="Longitude", y=NULL)
  }
  
  if(plot_leg | is.na(max_scale)){
    p <- p + theme(legend.position=c(0.83, 0.3), legend.background=element_rect(colour="black"))
  }else{
    p <- p + theme(legend.position = "none")
  }
  
  p
}

pmap1 <- mapguad("n_sit_cv", plot_tit="Horses", max_scale=3, y_lab=T, plot_leg=F)
pmap2 <- mapguad("n_sit_pl", plot_tit="Chickens", max_scale=3, y_lab=F, plot_leg=F)
pmap3 <- mapguad("n_sit_mosq", plot_tit="Mosquitoes", leg_tit="Number of\nsites sampled:", max_scale=3, y_lab=F, plot_leg=T)
# mapguad("n_samp_cv", "Number of\nsamples\nin horses:")
# mapguad("n_samp_pl", "Number of\nsamples\nin chickens:")

pmap <- ggarrange(pmap1, pmap2, pmap3, nrow=1, align="v", labels=c("A","C","D"), label.x=0.03, label.y=0.07)

all_coun = read_sf("./geoBoundariesCGAZ_ADM0/geoBoundariesCGAZ_ADM0.shp")
pcar <- ggplot() +
  geom_sf(data=all_coun, fill="grey80") +
  geom_segment(aes(x=-67, y=13, xend=-62.5, yend=15.75), col="#d95f02", linewidth=0.8, arrow=arrow(length=unit(0.15,"cm"))) +
  # geom_rect(aes(xmin=-62.75, xmax=-60.25, ymin=15.5, ymax=17), col="#d95f02", alpha=0) +
  geom_text(aes(label="B", x=-58.5, y=28.5), size=5, fontface='bold') +
  coord_sf(xlim = c(-115,-55), ylim = c(1,31)) +
  # geom_text(aes(label="B", x=-62, y=29), size=5, fontface='bold') +
  # coord_sf(xlim = c(-84,-58), ylim = c(6,32)) +
  theme_void() +
  theme(panel.background = element_rect(fill="white"),
        panel.border = element_rect(colour="black", fill=NA, linewidth=1))

pdf(file="./fig_map.pdf", width=15, height=5)
print(pmap)
print(pcar, vp=viewport(width = 0.22, height = 0.22, x = 0.27, y = 0.82))
dev.off()

pmap4 <- mapguad("prop_pos_cv", plot_tit="Horses", leg_tit="Proportion of\npositive samples:", y_lab=T, plot_leg=T)
pmap5 <- mapguad("prop_pos_pl", plot_tit="Chickens", leg_tit="Proportion of\npositive samples:", y_lab=F, plot_leg=T)
pmapbis <- ggarrange(pmap4, pmap5, nrow=1, align="v", labels=c("A","B"), label.x=0.06, label.y=0.08)

png(file="./fig_map_res.png", width=26, height=13, units="cm", res=300)
pmapbis
dev.off()

### Comparison of parameter values depending on model

to_mcmc.list <- function(samples1, to_burn, which_par = "all") {
  nchains = length(as.mcmc.list(samples1[[1]]))
  
  if (which_par == "all") {
    var_mod = names(samples1)[which(!(names(samples1) %in% c("pD")))]
  } else {
    var_mod = names(samples1)[which(!grepl("prob_sample|serop_true_sample|pD", names(samples1)))]
  }
  
  Mch = as.list(rep(NA, nchains))
  
  for (chain_i in 1:nchains) {
    for (na in var_mod) {
      post_var_chain_i = as.mcmc.list(samples1[[na]])[[chain_i]]
      Mch[[chain_i]] = cbind(Mch[[chain_i]], post_var_chain_i)
    }
    
    included_in_chain = (1:nrow(Mch[[chain_i]]))[! (1:nrow(Mch[[chain_i]])) %in% to_burn]
    Mch[[chain_i]] = mcmc(Mch[[chain_i]][included_in_chain, -1])
  }
  
  output_Mch = as.mcmc.list(Mch)
  output_combined_Mch = do.call("rbind", Mch)
  
  list(output_var_mod = var_mod, output_Mch = output_Mch, output_combined_Mch = output_combined_Mch)
}

summ_par <- function(param_i, Mark_Ch, comb_Mark_Ch, modlab){
  return(c(name_var = rownames(summary(Mark_Ch)[[2]])[param_i],
           median_post = summary(Mark_Ch)[[2]][param_i,"50%"],
           inf_post = HPDinterval(as.mcmc(comb_Mark_Ch))[param_i,1],
           sup_post = HPDinterval(as.mcmc(comb_Mark_Ch))[param_i,2],
           mod_lab = modlab))
}

load(file="./sample_model_SeasoVary.RDATA")
transformed_chain_main <- to_mcmc.list(samples1, to_burn = 1:10, which_par = "not_all")
combined_Mch_main = transformed_chain_main[["output_combined_Mch"]]
Mch_main = transformed_chain_main[["output_Mch"]]
tab_param_main <- as.data.frame(t(sapply(X=1:nrow(summary(Mch_main)[[2]]), FUN=summ_par, modlab="SeasoVary", Mark_Ch=Mch_main, comb_Mark_Ch=combined_Mch_main)))

load(file="./sample_model_FlatVary.RDATA")
transformed_chain_flat <- to_mcmc.list(samples1, to_burn = 1:10, which_par = "not_all")
combined_Mch_flat = transformed_chain_flat[["output_combined_Mch"]]
Mch_flat = transformed_chain_flat[["output_Mch"]]
tab_param_flat <- as.data.frame(t(sapply(X=1:nrow(summary(Mch_flat)[[2]]), FUN=summ_par, modlab="FlatVary", Mark_Ch=Mch_flat, comb_Mark_Ch=combined_Mch_flat)))

labvar = c("eta" = bquote(eta), "spec" = bquote(psi), "nu" = bquote(nu), "fact_pl" = bquote(beta), "NPV_1st_samp" = bquote(NPV[1]), "PPV_1st_samp" = bquote(PPV[1]),
           "min_prop" = bquote(epsilon), "week_peak" = bquote(delta), "param_P" =  bquote(P[1]), "param_alpha_1"=bquote(alpha[1]), "param_alpha_2"=bquote(alpha[2]))

tab_comb_par <- rbind(filter(tab_param_main, name_var %in% c("NPV_1st_samp","PPV_1st_samp","eta","spec","fact_pl","nu","param_P","param_alpha_1","param_alpha_2")), #,"min_prop","week_peak",paste0("foi_max[",1:16,"]"))),
                      filter(tab_param_flat, name_var %in% c("NPV_1st_samp","PPV_1st_samp","eta","spec","fact_pl","nu","param_P","param_alpha_1","param_alpha_2"))) #,"min_prop","week_peak",paste0("foi_max[",1:16,"]"))))
tab_comb_par$median_post <- as.numeric(tab_comb_par$median_post)
tab_comb_par$inf_post <- as.numeric(tab_comb_par$inf_post)
tab_comb_par$sup_post <- as.numeric(tab_comb_par$sup_post)
tab_comb_par$mod_lab <- factor(tab_comb_par$mod_lab, levels = unique(tab_comb_par$mod_lab))
tab_comb_par$lab_var <- sapply(X = tab_comb_par$name_var, FUN = function(x){labvar[x]})
tab_comb_par$name_var <- factor(tab_comb_par$name_var, levels=unique(tab_comb_par$name_var), labels=unique(tab_comb_par$lab_var))

p_compar_mod <- ggplot(data=tab_comb_par, aes(x=mod_lab, col=mod_lab, shape=mod_lab, y=median_post, ymin=inf_post, ymax=sup_post)) +
  geom_pointrange() +
  scale_color_brewer(name="Model:", type="qual", palette=2) +
  scale_shape_manual(name="Model:", values=c(15,16,17)) +
  facet_wrap(~name_var, scales="free_y", labeller = label_parsed) +
  xlab(NULL) + ylab("Posterior distribution") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "right", #c(0.88,0.2),
        legend.direction = "vertical")

png(file="./fig_compar_mod.png", width=16, height=10, units="cm", res=300)
p_compar_mod
dev.off()





