source("./load_dataset.R")

library(rjags)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pROC)
library(gridGraphics)

nb_sem_tab_chev <- ncol(paires_chv_df1)
nb_sem_tab_poul <- ncol(paires_pl_df1)

# Fonction de conversion
to_mcmc.list <- function(samples1, to_burn, which_par = "all") {
  nchains = length(as.mcmc.list(samples1[[1]]))
  
  if (which_par == "all") {
    var_mod = names(samples1)[which(!(names(samples1) %in% c("pD")))]
  } else {
    var_mod = names(samples1)[which(!grepl("prob_sample2_chev_egal_1|serop_true_sample|pD", names(samples1)))]
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

## Prepare donnees Chevaux pour plot

tab_chev <- data.frame(pair_id = as.character(NULL),
                       res_ech_1 = as.numeric(NULL),
                       res_ech_2 = as.numeric(NULL),
                       sem = as.numeric(NULL),
                       col = as.numeric(NULL))

for(pair_i in 1:nrow(paires_etat_serologique_chv1)){
  cat(paste0("\r", pair_i, "/", nrow(paires_etat_serologique_chv1)))
  
  res_ech_1 <- paires_etat_serologique_chv1[pair_i,"etat_serologique_1"]
  res_ech_2 <- paires_etat_serologique_chv1[pair_i,"etat_serologique_2"]
  nb_sem <- paires_etat_serologique_chv1[pair_i,"semaine_2"] - paires_etat_serologique_chv1[pair_i,"semaine_1"] + 1
  for(sem_i in 0:(nb_sem)){
    tab_chev <- rbind(tab_chev, data.frame(pair_id = pair_i,
                                           res_ech_1 = res_ech_1,
                                           res_ech_2 = res_ech_2,
                                           sem = paires_etat_serologique_chv1[pair_i,"semaine_1"] + sem_i,
                                           col = res_ech_1 + (res_ech_2-res_ech_1)*(sem_i/(nb_sem-1))))
  }
}

ordre_paire_id <- order(paires_etat_serologique_chv1$etat_serologique_1, paires_etat_serologique_chv1$etat_serologique_2, paires_etat_serologique_chv1$semaine_2, paires_etat_serologique_chv1$semaine_1)
tab_chev$order <- sapply(X = tab_chev$pair_id, FUN = function(x){which(ordre_paire_id == x)})

## Prepare donnees Poulets pour plot

tab_pl <- data.frame(pair_id = as.character(NULL),
                     res_ech_1 = as.numeric(NULL),
                     res_ech_2 = as.numeric(NULL),
                     sem = as.numeric(NULL),
                     col = as.numeric(NULL))

for(pair_i in 1:nrow(paires_etat_serologique_pl1)){
  cat(paste0("\r", pair_i, "/", nrow(paires_etat_serologique_pl1)))
  
  res_ech_1 <- paires_etat_serologique_pl1[pair_i,"etat_serologique_1"]
  res_ech_2 <- paires_etat_serologique_pl1[pair_i,"etat_serologique_2"]
  nb_sem <- paires_etat_serologique_pl1[pair_i,"semaine_2"] - paires_etat_serologique_pl1[pair_i,"semaine_1"] + 1
  for(sem_i in seq(0, nb_sem-0.5, 0.5)){
    tab_pl <- rbind(tab_pl, data.frame(pair_id = pair_i,
                                       res_ech_1 = res_ech_1,
                                       res_ech_2 = res_ech_2,
                                       sem = paires_etat_serologique_pl1[pair_i,"semaine_1"] + sem_i,
                                       col = res_ech_1 + (res_ech_2-res_ech_1)*(sem_i/(nb_sem-1))))
  }
}

ordre_paire_id <- order(paires_etat_serologique_pl1$etat_serologique_1, paires_etat_serologique_pl1$etat_serologique_2, paires_etat_serologique_pl1$semaine_2, paires_etat_serologique_pl1$semaine_1)
tab_pl$order <- sapply(X = tab_pl$pair_id, FUN = function(x){which(ordre_paire_id == x)})

## Simule modele pour plot

modplot <- function(mod, n_simu=1000){
  load(file = paste0("./sample_model_", mod, ".RDATA"))
  
  transformed_chain <- to_mcmc.list(samples1, to_burn = 1:10, which_par = "not_all")
  
  var_mod = transformed_chain[["output_var_mod"]]
  combined_Mch = transformed_chain[["output_combined_Mch"]]
  
  simu_foi <- data.frame(simu = as.numeric(NULL),
                         sem = as.numeric(NULL),
                         valfoi = as.numeric(NULL))
  vsem <- min(min(tab_chev$sem), min(tab_pl$sem)):(which(vec_years == 2019)[1])
  roc_tab_chev = roc_tab_poul = data.frame(simu = as.numeric(NULL),
                                           sensi = as.numeric(NULL),
                                           speci = as.numeric(NULL))
  AUC_moy_chev = AUC_moy_poul = 0
  cumul_incid_rate <- data.frame()
  
  for(simu_i in 1:n_simu){
    cat(paste0("\r", simu_i, "/", n_simu))
    
    rand_i <- sample(1:nrow(combined_Mch), 1)
    valnu <- combined_Mch[rand_i, "nu"]
    valfacpl <- combined_Mch[rand_i, "fact_pl"]
    
    if(grepl("FlatStable", mod)){
      vfoi <- combined_Mch[rand_i, "foi"]
      
      simu_foi <- rbind(simu_foi, data.frame(simu = simu_i,
                                             sem = vsem,
                                             valfoi = vfoi))

    }else if(grepl("FlatVary", mod)){
      vfoi <- combined_Mch[rand_i, paste0("foi[",year_num[vsem],"]")]

      simu_foi <- rbind(simu_foi, data.frame(simu = simu_i,
                                             sem = vsem,
                                             valfoi = vfoi))
      
      foi_cum_chev <- rep(0, nrow(paires_chv_df1))
      for (sem_j in 1:(nb_sem_tab_chev - 1)) {
        foi_cum_chev <- foi_cum_chev + combined_Mch[rand_i, paste0("foi[",year_num[paires_chv_df1[,sem_j]],"]")] * (paires_chv_df1[,sem_j + 1] - paires_chv_df1[,sem_j])
      }

      foi_cum_poul <- rep(0, nrow(paires_pl_df1))
      for (sem_j in 1:(nb_sem_tab_poul - 1)) {
        foi_cum_poul <- foi_cum_poul + valfacpl * combined_Mch[rand_i, paste0("foi[",year_num[paires_pl_df1[,sem_j]],"]")] * (paires_pl_df1[,sem_j + 1] - paires_pl_df1[,sem_j])
      }

    }else if(grepl("SeasoStable", mod)){
      vfoimax <- combined_Mch[rand_i, "foi_max"]
      valminprop <- combined_Mch[rand_i, "min_prop"]
      valweekpeak <- combined_Mch[rand_i, "week_peak"]

      simu_foi <- rbind(simu_foi, data.frame(simu = simu_i,
                                             sem = vsem,
                                             valfoi = (vfoimax/2) * (1-valminprop) * (cos((vsem - valweekpeak)*3.142/26) + 1) + valminprop*vfoimax))
      
      foi_cum_chev <- ((1-valminprop)*vfoimax/2) * ((26 / 3.14) * (sin((paires_etat_serologique_chv1$semaine_2 - valweekpeak) * 3.14 / 26) - sin((paires_etat_serologique_chv1$semaine_1 - valweekpeak) * 3.14 / 26)) + (paires_etat_serologique_chv1$semaine_2 - paires_etat_serologique_chv1$semaine_1)) + valminprop*vfoimax*(paires_etat_serologique_chv1$semaine_2 - paires_etat_serologique_chv1$semaine_1)

      foi_cum_poul <- (valfacpl*(1-valminprop)*vfoimax/2) * ((26 / 3.14) * (sin((paires_etat_serologique_pl1$semaine_2 - valweekpeak) * 3.14 / 26) - sin((paires_etat_serologique_pl1$semaine_1 - valweekpeak) * 3.14 / 26)) + (paires_etat_serologique_pl1$semaine_2 - paires_etat_serologique_pl1$semaine_1)) + valfacpl*valminprop*vfoimax*(paires_etat_serologique_pl1$semaine_2 - paires_etat_serologique_pl1$semaine_1)

    }else if(grepl("SeasoVary", mod)){
      vfoimax <- combined_Mch[rand_i, paste0("foi_max[",year_num[vsem],"]")]
      valminprop <- combined_Mch[rand_i, "min_prop"]
      valweekpeak <- combined_Mch[rand_i, "week_peak"]

      simu_foi <- rbind(simu_foi, data.frame(simu = simu_i,
                                             sem = vsem,
                                             valfoi = (vfoimax/2) * (1-valminprop) * (cos((vsem - valweekpeak)*3.142/26) + 1) + valminprop*vfoimax))
      
      cumul_incid_rate <- rbind(cumul_incid_rate,
                                data.frame(cum_inc_2002 = ((1-valminprop)*combined_Mch[rand_i,"foi_max[1]"]/2) * ((26 / 3.14) * (sin((72 - valweekpeak) * 3.14 / 26) - sin((20 - valweekpeak) * 3.14 / 26)) + (72 - 20)) + valminprop*combined_Mch[rand_i,"foi_max[1]"]*(72 - 20),
                                           cum_inc_2007 = ((1-valminprop)*combined_Mch[rand_i,"foi_max[6]"]/2) * ((26 / 3.14) * (sin((332 - valweekpeak) * 3.14 / 26) - sin((281 - valweekpeak) * 3.14 / 26)) + (332 - 281)) + valminprop*combined_Mch[rand_i,"foi_max[6]"]*(332 - 281),
                                           cum_inc_2012 = ((1-valminprop)*combined_Mch[rand_i,"foi_max[11]"]/2) * ((26 / 3.14) * (sin((593 - valweekpeak) * 3.14 / 26) - sin((542 - valweekpeak) * 3.14 / 26)) + (593 - 542)) + valminprop*combined_Mch[rand_i,"foi_max[11]"]*(593 - 542)))
      
      foi_cum_chev <- rep(0, nrow(paires_chv_df1))
      for (sem_j in 1:(nb_sem_tab_chev - 1)) {
        v_foimax_j <- combined_Mch[rand_i, paste0("foi_max[",year_num[paires_chv_df1[,sem_j]],"]")]
        foi_cum_chev <- foi_cum_chev +  ((1-valminprop)*v_foimax_j/2) * ((26 / 3.14) * (sin((paires_chv_df1[,sem_j + 1] - valweekpeak) * 3.14 / 26) - sin((paires_chv_df1[,sem_j] - valweekpeak) * 3.14 / 26)) + (paires_chv_df1[,sem_j + 1] - paires_chv_df1[,sem_j])) + valminprop*v_foimax_j*(paires_chv_df1[,sem_j + 1] - paires_chv_df1[,sem_j])
      }

      foi_cum_poul <- rep(0, nrow(paires_pl_df1))
      for (sem_j in 1:(nb_sem_tab_poul - 1)) {
        v_foimax_j <- combined_Mch[rand_i, paste0("foi_max[",year_num[paires_pl_df1[,sem_j]],"]")]
        foi_cum_poul <- foi_cum_poul +  (valfacpl*(1-valminprop)*v_foimax_j/2) * ((26 / 3.14) * (sin((paires_pl_df1[,sem_j + 1] - valweekpeak) * 3.14 / 26) - sin((paires_pl_df1[,sem_j] - valweekpeak) * 3.14 / 26)) + (paires_pl_df1[,sem_j + 1] - paires_pl_df1[,sem_j])) + valfacpl*valminprop*v_foimax_j*(paires_pl_df1[,sem_j + 1] - paires_pl_df1[,sem_j])
      }
    }
  }

  simu_foi_aggr <- as.data.frame(simu_foi %>%
                                   # filter(sem <= max(max(tab_chev$sem), max(tab_pl$sem))) %>%
                                   filter(sem < which(vec_years == 2019)[1], sem >= min(min(tab_chev$sem), min(tab_pl$sem))) %>%
                                   group_by(sem) %>%
                                   summarise(pred_val = median(valfoi),
                                             pred_025 = quantile(valfoi, probs=0.025),
                                             pred_10 = quantile(valfoi, probs=0.10),
                                             pred_25 = quantile(valfoi, probs=0.25),
                                             pred_75 = quantile(valfoi, probs=0.75),
                                             pred_90 = quantile(valfoi, probs=0.90),
                                             pred_975 = quantile(valfoi, probs=0.975)))
  
  ### Cree plot
  
  p1 = ggplot(data=simu_foi_aggr, aes(x=sem)) +
    geom_ribbon(aes(ymin=pred_10, ymax=pmin(pred_90, 0.075)), fill="gray65") +
    # geom_ribbon(aes(ymin=pred_25, ymax=pmin(pred_75, 0.075)), fill="gray80") +
    # geom_ribbon(aes(ymin=pred_025, ymax=pmin(pred_975, 0.075)), fill="gray50") +
    geom_line(linewidth=1, aes(y=pred_val)) +
    ylab("Predicted force of infection") + xlab(NULL) +
    ylim(0, NA) + #min(0.055, max(simu_foi_aggr$pred_sup))) +
    scale_x_continuous(limits = c(min(simu_foi_aggr$sem), max(simu_foi_aggr$sem)),
                       breaks = seq(min(simu_foi_aggr$sem)+1, max(simu_foi_aggr$sem), 52),
                       labels = format(vec_weeks[seq(min(simu_foi_aggr$sem)+1, max(simu_foi_aggr$sem), 52)], "%Y-%m")) +
    theme_bw() +
    # theme(axis.text.x = element_text(angle=90, vjust=0.4))
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 = ggplot() +
    # geom_line(data=filter(tab_chev, res_ech_1==0, res_ech_2==0), aes(x=sem, y=-0.0001*order, group=pair_id, col=col)) +
    geom_line(data=filter(tab_chev, res_ech_1==0), aes(x=sem, y=-0.0001*order, group=pair_id, col=col), linewidth=0.25) +
    geom_ribbon(data=data.frame(fill="Negative"), aes(x=0, ymin=0, ymax=0, fill=fill)) +
    geom_ribbon(data=data.frame(fill="Positive"), aes(x=0, ymin=0, ymax=0, fill=fill)) +
    scale_color_gradient(low="#1b9e77", high="#d95f02", guide="none", limits=c(0,1)) +
    scale_fill_manual(name="Serology:", labels=c("Negative","Positive"), values=c("#1b9e77","#d95f02")) +
    ylab("Horses\nSerological data") + xlab(NULL) +
    scale_x_continuous(limits = c(min(simu_foi_aggr$sem), max(simu_foi_aggr$sem)),
                       breaks = seq(min(simu_foi_aggr$sem)+1, max(simu_foi_aggr$sem), 52),
                       labels = format(vec_weeks[seq(min(simu_foi_aggr$sem)+1, max(simu_foi_aggr$sem), 52)], "%Y-%m")) +
    scale_y_continuous(breaks = -0.0001 * max(filter(tab_chev, res_ech_1==0)$order) /2,
                       labels = "Pairs of samples") +
    theme_bw() +
    theme(axis.text.y = element_text(angle=90, hjust=0.5),
          axis.ticks.y = element_blank(),
          legend.position = c(0.85, 0.7),
          # axis.text.x = element_text(angle=90, vjust=0.4))
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p3 = ggplot() +
    # geom_line(data=filter(tab_pl, res_ech_1==0, res_ech_2==1), aes(x=sem, y=-0.0001*order, group=pair_id, col=col)) +
    geom_line(data=filter(tab_pl, res_ech_1==0), aes(x=sem, y=-0.0001*order, group=pair_id, col=col), linewidth=0.25) +
    scale_color_gradient(low="#1b9e77", high="#d95f02", guide="none", limits=c(0,1)) +
    ylab("Chickens\nSerological data") + xlab("Time") +
    scale_x_continuous(limits = c(min(simu_foi_aggr$sem), max(simu_foi_aggr$sem)),
                       breaks = seq(min(simu_foi_aggr$sem)+1, max(simu_foi_aggr$sem), 52),
                       labels = format(vec_weeks[seq(min(simu_foi_aggr$sem)+1, max(simu_foi_aggr$sem), 52)], "%Y-%m")) +
    scale_y_continuous(breaks = -0.0001 * max(filter(tab_pl, res_ech_1==0)$order) /2,
                       labels = "Pairs of samples") +
    
    geom_rect(aes(xmin=600, xmax=670, ymin=-0.0001*3200, ymax=-0.0001*3415), col="black", alpha=0, linewidth=0.5) +
    geom_segment(aes(x=600, y=-0.0001*3415, xend=398, yend=-0.0001*2050), col="black", linewidth=0.5, linetype="dashed") +
    geom_segment(aes(x=670, y=-0.0001*3200, xend=527, yend=-0.0001*240), col="black", linewidth=0.5, linetype="dashed") +
  
    theme_bw() +
    theme(axis.text.y = element_text(angle=90, hjust=0.5),
          axis.ticks.y = element_blank(),
          legend.position = c(0.85, 0.8),
          axis.text.x = element_text(angle=90, vjust=0.4))
  
  tab_seroconv_pl <- filter(tab_pl, res_ech_1==0, res_ech_2==1)
  p3bis = ggplot() +
    # geom_line(data=filter(tab_pl, res_ech_1==0, res_ech_2==1), aes(x=sem, y=-0.0001*order, group=pair_id, col=col)) +
    geom_line(data=tab_seroconv_pl, aes(x=sem, y=-0.0001*order, group=pair_id, col=col), linewidth=3) +
    ylim(min(-0.0001*tab_seroconv_pl$order)-0.0001, max(-0.0001*tab_seroconv_pl$order)+0.0001) +
    scale_color_gradient(low="#1b9e77", high="#d95f02", guide="none", limits=c(0,1)) +
    ylab(NULL) + xlab(NULL) +
    scale_x_continuous(limits = c(621, 648),
                       breaks = seq(621+1, 648, 5),
                       labels = vec_weeks[seq(621+1, 648, 5)]) +
    theme_void() + #theme_bw() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.background = element_rect(fill="white"),
          panel.border = element_rect(colour="black", fill=NA, size=0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
          # axis.text.x = element_text(angle=90, vjust=0.4), text = element_text(size = 8))
  
  pfit = ggarrange(p1, p2, p3, ncol=1, align="v", labels=c("A","B","C"))
  
  pdf(file=paste0("./fig_fit_", mod, ".pdf"), width=6, height=7)
  print(pfit)
  print(p3bis, vp=viewport(width = 0.12, height = 0.1, x = 0.55, y = 0.25))
  dev.off()
  
  cumul_incid_rate <- 1-exp(-cumul_incid_rate)
  return(cumul_incid_rate)
}

modplot("FlatStable")
modplot("FlatVary")
modplot("SeasoStable")
cumul_inc_rate <- modplot("SeasoVary", n_simu=5000)

summary(cumul_inc_rate)

