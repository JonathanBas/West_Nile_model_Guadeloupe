
####################################
### Load data
####################################

source("./load_dataset.R")

library(rjags)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

###définir nb_sem_tab pour le modèle JAGS 
nb_sem_tab_chev <- ncol(paires_chv_df1)
nb_sem_tab_poul <- ncol(paires_pl_df1)

####################################
### Run MCMC
####################################

load.module("dic")

### MODELE FlatStable

jags1 <- jags.model("./FlatStable.BUGS",
                    data = list(
                      serop_pair_chev = paires_etat_serologique_chv1,
                      serop_pair_poul = paires_etat_serologique_pl1,
                      N_pairs_chev = N_pairs_chev,
                      N_pairs_poul = N_pairs_poul,
                      list_sem_chev = paires_chv_df1,
                      list_sem_poul = paires_pl_df1,
                      nb_sem_tab_chev = nb_sem_tab_chev,
                      nb_sem_tab_poul = nb_sem_tab_poul,
                      tab_indiv_pair_chev = tab_indiv_pair_chev,
                      tab_indiv_pair_poul = tab_indiv_pair_poul),
                    n.chains = 3,
                    n.adapt = 0)

samples1 <- jags.samples(jags1,
                         c("foi","nu","deviance","fact_pl","eta","param_alpha_1","param_alpha_2","param_P","NPV_1st_samp","pD"),
                         # c("foi","nu","deviance","fact_pl","eta","param_alpha_1","param_alpha_2","param_P","NPV_1st_samp","spec","PPV_1st_samp","pD"),
                         40000,
                         thin = 200)
save(samples1, file="./sample_model_FlatStable.RDATA")

### MODELE FlatVary

jags1 <- jags.model("./FlatVary.BUGS",
                    data = list(
                      serop_pair_chev = paires_etat_serologique_chv1,
                      serop_pair_poul = paires_etat_serologique_pl1,
                      N_pairs_chev = N_pairs_chev,
                      N_pairs_poul = N_pairs_poul,
                      list_sem_chev = paires_chv_df1,
                      list_sem_poul = paires_pl_df1,
                      year_num = year_num,
                      year_max = max(year_num),
                      nb_sem_tab_chev = nb_sem_tab_chev,
                      nb_sem_tab_poul = nb_sem_tab_poul,
                      tab_indiv_pair_chev = tab_indiv_pair_chev,
                      tab_indiv_pair_poul = tab_indiv_pair_poul),
                    n.chains = 3,
                    n.adapt = 0)

samples1 <- jags.samples(jags1,
                         c("foi","nu","deviance","fact_pl","eta","param_alpha_1","param_alpha_2","param_P","NPV_1st_samp","pD"),
                         40000,
                         thin = 200)
save(samples1, file="./sample_model_FlatVary.RDATA")

### MODELE SeasoStable

jags1 <- jags.model("./SeasoStable.BUGS",
                    data = list(serop_pair_chev = paires_etat_serologique_chv1,
                                serop_pair_poul = paires_etat_serologique_pl1,
                                N_pairs_chev = N_pairs_chev,
                                N_pairs_poul = N_pairs_poul,
                                list_sem_chev = paires_chv_df1,
                                list_sem_poul = paires_pl_df1,
                                nb_sem_tab_chev = nb_sem_tab_chev,
                                nb_sem_tab_poul = nb_sem_tab_poul,
                                tab_indiv_pair_chev = tab_indiv_pair_chev,
                                tab_indiv_pair_poul = tab_indiv_pair_poul),
                    n.chains = 3,
                    n.adapt = 0)

samples1 <- jags.samples(jags1,
                         c("foi_max","nu","deviance","fact_pl","week_peak","min_prop","eta","param_alpha_1","param_alpha_2","param_P","NPV_1st_samp","pD"),
                         40000,
                         thin = 200)
save(samples1, file="./sample_model_SeasoStable.RDATA")

### MODELE SeasoVary

jags1 <- jags.model("./SeasoVary.BUGS",
                    data = list(
                      serop_pair_chev = paires_etat_serologique_chv1,
                      serop_pair_poul = paires_etat_serologique_pl1,
                      N_pairs_chev = N_pairs_chev,
                      N_pairs_poul = N_pairs_poul,
                      list_sem_chev = paires_chv_df1,
                      list_sem_poul = paires_pl_df1,
                      year_num = year_num,
                      year_max = max(year_num),
                      nb_sem_tab_chev = nb_sem_tab_chev,
                      nb_sem_tab_poul = nb_sem_tab_poul,
                      tab_indiv_pair_chev = tab_indiv_pair_chev,
                      tab_indiv_pair_poul = tab_indiv_pair_poul),
                    n.chains = 3,
                    n.adapt = 0)

samples1 <- jags.samples(jags1,
                         c("foi_max","nu","deviance","fact_pl","week_peak","min_prop","eta","param_alpha_1","param_alpha_2","param_P","NPV_1st_samp","pD"), #,"prob_sample2_chev_egal_1","prob_sample2_poul_egal_1"),
                         # c("foi_max","nu","deviance","fact_pl","week_peak","min_prop","eta","param_alpha_1","param_alpha_2","param_P","NPV_1st_samp","spec","PPV_1st_samp","pD"),
                          40000,
                         thin = 200)
save(samples1, file="./sample_model_SeasoVary.RDATA")

####################################
### Analyze MCMC outputs
####################################

load(file="./sample_model_FlatStable.RDATA") ; mod = "FlatStable"
load(file="./sample_model_FlatVary.RDATA") ; mod = "FlatVary"
load(file="./sample_model_SeasoStable.RDATA") ; mod = "SeasoStable"
load(file="./sample_model_SeasoVary.RDATA") ; mod = "SeasoVary"

####CALCUL DU DIC : 
summary(samples1$deviance, mean)$stat + summary(samples1$pD, mean)$stat
DIC_all_mod = function(which_mod, path = "./c.hamouche/"){
  DIC_mod = rep(NA, length(which_mod))
  names(DIC_mod) = as.character(which_mod)
  
  for (n_mod in which_mod){
    load(paste0(path, "sample_model_", n_mod, ".RDATA"))
    DIC_mod[as.character(n_mod)] = summary(samples1$deviance, mean)$stat + summary(samples1$pD, mean)$stat
  }
  DIC_mod
}
DIC_all_mod(path = "./", which_mod = c("FlatStable", "FlatVary", "SeasoStable", "SeasoVary"))

# Fonction de conversion
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

# Convertir les échantillons JAGS en mcmc.list en spécifiant les itérations à brûler
transformed_chain <- to_mcmc.list(samples1, to_burn = 1:10, which_par = "not_all")
var_mod = transformed_chain[["output_var_mod"]]
combined_Mch = transformed_chain[["output_combined_Mch"]]
Mch = transformed_chain[["output_Mch"]]
plot(Mch)

tab_param <- t(sapply(X=1:nrow(summary(Mch)[[2]]), FUN=function(param_i, rou=3){
  return(c(rownames(summary(Mch)[[2]])[param_i], paste0(round(summary(Mch)[[2]][param_i,"50%"], rou), " [", round(HPDinterval(as.mcmc(combined_Mch))[param_i,1], rou), " ; ", round(HPDinterval(as.mcmc(combined_Mch))[param_i,2], rou), "]")))
}))
tab_param

par(mfrow = c(3, 3)) ; for(i in var_mod){hist(combined_Mch[,which(grepl(i, colnames(combined_Mch)))], nclass=20, main = i, col = "skyblue")} ; par(mfrow = c(1, 1))
par(mfrow = c(1, 1))
par(mfrow = c(3, 3)) ; for(i in var_mod){acf(combined_Mch[,which(grepl(i, colnames(combined_Mch)))], lag.max = 30, ylim = c(-1,1), main = i)} ; par(mfrow = c(1, 1))
print(cor(combined_Mch))
ggplot(data = melt(cor(combined_Mch)), aes(x = Var1, y = Var2, fill = value, col = (abs(value) >0.4))) + geom_tile(linewidth = 1.5) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + scale_colour_manual(values = c("TRUE" = "forestgreen", "FALSE" = "white"))
print(effectiveSize(Mch))
print(gelman.diag(Mch))
gelman.plot(Mch) # should be < 1.05 (or 1.1 in other sources)

fit_poster_prior <- function(param, name_param, prior_dis, plot_leg=F, fit_prior=T){
  poster_dis <- combined_Mch[,param]
  
  distr = data.frame(Distribution = c(rep("Prior",length(prior_dis)), rep("Posterior",length(poster_dis))),
                     val = c(prior_dis, poster_dis))
  distr$Distribution = factor(distr$Distribution, levels=c("Prior", "Posterior"))
  
  if(! fit_prior){distr <- filter(distr, Distribution=="Posterior")}
  
  distrib_plot <- ggplot(data = distr, aes(x=val, fill=Distribution, alpha=Distribution)) +
    geom_density(color="black") +
    xlab(name_param) + ylab(NULL) +
    scale_fill_manual(name="Distribution:", values=c("Prior"="#d95f02", "Posterior"="#1b9e77")) +
    scale_alpha_manual(name="Distribution:", values=c("Prior"=0.3, "Posterior"=0.8)) +
    theme_bw() +
    theme(legend.position="none", axis.text.x=element_text(angle=90, vjust=0.4))
  if(plot_leg){distrib_plot = distrib_plot + theme(legend.position = c(0.35, 0.65))}
  
  return(distrib_plot)
}

p_param_fit <- ggarrange(fit_poster_prior("eta", expression(eta), runif(100000, 0, 1), plot_leg=T, fit_prior = T),
                         fit_poster_prior("min_prop", expression(epsilon), rbeta(100000, 3.80, 22.43), fit_prior = T),
                         fit_poster_prior("week_peak", expression(delta), rnorm(100000, 45.26, 0.988), fit_prior = T),
                         fit_poster_prior("nu", expression(mu), runif(100000, 0, 0.5), fit_prior = F),
                         fit_poster_prior("fact_pl", expression(beta), runif(100000, 0, 10), fit_prior = T),
                         fit_poster_prior("NPV_1st_samp", bquote(NPV[1]), runif(100000, 0, 1), fit_prior = F),
                         fit_poster_prior("foi_max[1]", expression(Lambda["2002"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[2]", expression(Lambda["2003"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[3]", expression(Lambda["2004"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[4]", expression(Lambda["2005"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[5]", expression(Lambda["2006"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[6]", expression(Lambda["2007"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[7]", expression(Lambda["2008"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[8]", expression(Lambda["2009"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[9]", expression(Lambda["2010"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[10]", expression(Lambda["2011"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[11]", expression(Lambda["2012"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[12]", expression(Lambda["2013"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[13]", expression(Lambda["2014"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[14]", expression(Lambda["2015"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[15]", expression(Lambda["2016"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         fit_poster_prior("foi_max[16]", expression(Lambda["2017"]), exp(runif(100000, -20, 3)), fit_prior = F),
                         align="hv", ncol=4, nrow=6)
p_param_fit = annotate_figure(p_param_fit, left = text_grob("Density", size=11, rot=90))

png(filename="./fig_estim_param.png",pointsize=6,res=300,width = 22.5, height = 33, units = "cm")
p_param_fit
dev.off()

