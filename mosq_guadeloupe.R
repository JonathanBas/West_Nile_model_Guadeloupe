library(readxl)
library(rjags)
library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(ggspatial)

to_mcmc.list = function(mcarray_obj, to_burn, which_par = "all"){
  nchains = length(as.mcmc.list(mcarray_obj[[1]]))
  
  if(which_par == "all"){
    var_mod = names(mcarray_obj)[which(!(names(mcarray_obj) %in% c("pD")))]
  }else{
    var_mod = names(mcarray_obj)[which(!grepl("pD|site_eff[1]", names(mcarray_obj)))]
  }
  
  Mch = as.list(rep(NA, nchains))
  for (chain_i in 1:nchains){
    for (na in var_mod){
      post_var_chain_i = as.mcmc.list(mcarray_obj[[na]])[[chain_i]]
      Mch[[chain_i]] = cbind(Mch[[chain_i]], post_var_chain_i)
    }
    included_in_chain = (1:nrow(Mch[[chain_i]]))[! (1:nrow(Mch[[chain_i]])) %in% to_burn]
    Mch[[chain_i]] = mcmc(Mch[[chain_i]][included_in_chain,-1])
  }
  output_Mch = as.mcmc.list(Mch)
  output_combined_Mch = do.call("rbind", Mch)
  
  list(output_var_mod = var_mod, output_Mch=output_Mch, output_combined_Mch=output_combined_Mch)
}

### Importe, prepare et plot les donnees Moustiques

bdd_mous_1 <- as.data.frame(read_excel("./Field_Diagnose_v5_Roxane_Moustiques 2015-2021.xlsx",
                                     sheet="Field_Samples", skip=1) %>%
                            filter(Habitat == "Farm") %>%
                              rename(lati = X...11,
                                     longi = Y...12) %>%
                            mutate(Date = as.Date(Date, format = "%Y%m%d")) %>%
                            mutate(sem_pvmt = format(Date, format = "%Y-%U")))

ttes_semaines = unique(format(seq.Date(from=min(bdd_mous_1$Date), to=max(bdd_mous_1$Date), by=1), "%Y-%U"))
any(ttes_semaines != sort(ttes_semaines)) # must be FALSE

bdd_mous_2 <- as.data.frame(read_excel("./Field_Diagnose_v5_Roxane_Moustiques 2015-2021.xlsx",
                                       sheet="Field_Samples_Diagnose") %>%
                              filter(grepl("Cx.|Culex", Species)) %>%
                              dplyr::select(ID_FS, Individuals) %>%
                              mutate(Individuals = as.numeric(Individuals)))

ttes_sem = expand.grid(sem_pvmt = ttes_semaines, Site = unique(bdd_mous_1$Site))
bdd_mous_aggr <- as.data.frame(bdd_mous_1 %>%
                                 merge(bdd_mous_2, by.x="ID_F", by.y="ID_FS", all.x=T) %>%
                                 replace_na(list(Individuals = 0)) %>%
                                 group_by(sem_pvmt, Site) %>%
                                 summarise(abun_mosq = mean(Individuals)) %>%
                                 merge(ttes_sem, by=c("sem_pvmt", "Site"), all.y=T))

### Ajuste modele aux donnees Moustiques
## AVEC effet fixe "Site"

bdd_mous_aggr_2 <- as.data.frame(bdd_mous_aggr %>% pivot_wider(id_cols=sem_pvmt, values_from=abun_mosq, names_from=Site))
jags <- jags.model("./mod_longitu_mosquito_Guad_2.BUGS",
                   data = list(foi_week = bdd_mous_aggr_2[,(2:5)],
                               last_week = nrow(bdd_mous_aggr_2)),
                   n.chains = 4,
                   n.adapt = 0)
samples <- jags.samples(jags,
                        c("foi_max", "week_peak","stand_dev","min_prop","site_eff"),
                        50000,
                        thin = 100)

save(samples, file="./samples_fit_mosqui_data_Guadeloupe.RDATA")

### Charge les estimations du modele

load(file="./samples_fit_mosqui_data_Guadeloupe.RDATA")

transformed_chain = to_mcmc.list(mcarray_obj = samples, to_burn = 1:200, which_par = "not_all")
var_mod = transformed_chain[["output_var_mod"]]
Mch = transformed_chain[["output_Mch"]]
combined_Mch = transformed_chain[["output_combined_Mch"]]
plot(Mch) ; par(mfrow = c(1, 1))
summary(Mch)[[2]][,"50%"]
HPDinterval(as.mcmc(combined_Mch))

print(cor(combined_Mch))
print(effectiveSize(Mch))
print(gelman.diag(Mch))
gelman.plot(Mch) # should be < 1.05 (or 1.1 in other sources)
par(mfrow = c(3, 3)) ; for(i in colnames(combined_Mch)){acf(combined_Mch[,i], lag.max = 30, ylim = c(-1,1), main = i)} ; par(mfrow = c(1, 1))

### Simule et plot l'abondance de moustiques observee et predite par le modele

foi_wk = data.frame(Site=as.character(NULL), simu=as.numeric(NULL), week=as.numeric(NULL), week_name=as.character(NULL), foi=as.numeric(NULL))
n_wks = nrow(bdd_mous_aggr_2)
n_simu = 500
for(simu_i in 1:n_simu){
  cat(paste0(simu_i, "/", n_simu, "\r"))
  part_i = sample(x=1:nrow(combined_Mch), size=1)
  param_i = combined_Mch[part_i,]
  
  for(site_i in 1:4){
    foi_week_esp <- param_i[paste0("site_eff[", site_i, "]")] * ((param_i["foi_max"]/2) * (1-param_i["min_prop"]) * (cos(((1:n_wks) - param_i["week_peak"])*3.142/26) + 1) + param_i["min_prop"]*param_i["foi_max"])
    
    foi_wk = rbind(foi_wk,
                   data.frame(Site = rep(colnames(bdd_mous_aggr_2)[2:5][site_i], n_wks),
                              simu = rep(simu_i, n_wks),
                              week = 1:n_wks,
                              week_name = bdd_mous_aggr_2$sem_pvmt,
                              foi = foi_week_esp))
    
  }
}

foi_wk_aggr <- as.data.frame(foi_wk %>%
                               group_by(Site, week, week_name) %>%
                               summarise(inf = quantile(foi, 0.025),
                                         val = median(foi),
                                         sup = quantile(foi, 0.975)) %>%
                               mutate(week = as.numeric(week)) %>%
                               merge(y=bdd_mous_aggr, by.x=c("Site","week_name"), by.y=c("Site","sem_pvmt")))

p_fit_abun = ggplot() +
  geom_bar(data=foi_wk_aggr, aes(x=week_name, y=abun_mosq), stat="identity", width=1) +
  geom_ribbon(data=foi_wk_aggr, aes(x=week_name, ymin=inf, ymax=sup, group=1), linetype="dashed", alpha=0, col=1) +
  geom_line(data=foi_wk_aggr, aes(x=week_name, y=val, group=1)) +
  xlab("Weeks") + ylab("Mosquito abundance") +
  facet_wrap(~Site, ncol=2, labeller = as_labeller(c("GP001"="Site 1", "GP002"="Site 2", "GP009"="Site 3", "GP010"="Site 4"))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size=c(8,rep(0, 20)), colour=c("black",rep("white", 20))),
        axis.ticks.x = element_blank())
p_fit_abun

pdf(file="./fig_fit_abun_mosqui.pdf", width=6, height=4.5)
p_fit_abun
dev.off()

### Bout de code pour ajuster une distribution Beta a la posterior du parametre "min_prop"

library(bbmle)

min_prop_val <- combined_Mch[,"min_prop"]

ll_beta <- function(alpha, beta){
  return(-sum(dbeta(x=min_prop_val, shape1=alpha, shape2=beta, log=T), na.rm=T))
}

fit_mle2 = mle2(minuslogl = ll_beta, start = list(alpha=0.5, beta=0.5), lower=c(alpha=0, beta=0))
fit_coef = coef(fit_mle2)
print(fit_coef)

### Plot les distributions de deux parametres apres ajustement aux donnees Moustiques

week_peak_val <- combined_Mch[,"week_peak"]
p_week_peak = ggplot() +
  geom_density(aes(x=week_peak_val), fill="#69b3a2", color="#69b3a2", alpha=0.8) +
  geom_density(aes(x = rnorm(100000,
                             45.26,
                             0.988)), linewidth=1) +
  xlab(expression(delta)) + ylab(NULL) +
  theme_bw()
p_week_peak

min_prop_val <- combined_Mch[,"min_prop"]
p_min_prop = ggplot() +
  geom_density(aes(x=min_prop_val), fill="#69b3a2", color="#69b3a2", alpha=0.8) +
  geom_density(aes(x = rbeta(100000,
                              shape1 = 3.80,
                              shape2 =  22.43)), linewidth=1) +
  xlab(expression(epsilon)) + ylab(NULL) +
  theme_bw()
p_min_prop

library(ggpubr)
pdf(file="./fig_prior_two_param.pdf", width=4, height=2)
ggarrange(p_week_peak, p_min_prop, ncol=2, labels=c("A", "B"), label.x=0.8, label.y=0.95)
dev.off()


