library(DAMOCLES)

gammas = c(0.01, 0.05, 0.5, 1) 
mus = c(0, 0.05, 0.5, 1)

##### Estimate parameters for colonization and local extinction #####
gammas <- seq(from = 0.1, to = 1.0, by = 0.15)
mus <- seq(from = 0.1, to = 1.0, by = 0.15)

## Simulate a yule phylogeny
set.seed(12345)
trsYule <- geiger::sim.bdtree(b = 0.2, d = 0, stop = 'taxa', n = 100)

## Simulate local assemblages under know rates of colonization and local extinction

sim0.1 <- list()
sim0.25 <- list()
sim0.4 <- list()
sim0.55 <- list()
sim0.7 <- list()
sim0.85 <- list()
sim1.0 <- list()
simMu <- list()
for(i in 1:length(gammas)){
  svMisc::progress(i, max.value = length(gammas))
  sim0.1[[i]] <- DAMOCLES_sim(trsYule, gamma_0 = gammas[i], gamma_td = 0, mu = mus[1], sigma = 0, psiBranch = 0, 
                            psiTrait = 0, z = 10, phi = 0, traitOpt = 1, br0 = 0.1, br_td = -0.1,
                            nTdim = 2, root.state = 1, root.trait.state = 0, plotit = FALSE,
                            keepExtinct = FALSE)[[1]]
  sim0.25[[i]] <- DAMOCLES_sim(trsYule, gamma_0 = gammas[i], gamma_td = 0, mu = mus[2], sigma = 0, psiBranch = 0, 
                                psiTrait = 0, z = 10, phi = 0, traitOpt = 1, br0 = 0.1, br_td = -0.1,
                                nTdim = 2, root.state = 1, root.trait.state = 0, plotit = FALSE,
                                keepExtinct = FALSE)[[1]]
  sim0.4[[i]] <- DAMOCLES_sim(trsYule, gamma_0 = gammas[i], gamma_td = 0, mu = mus[3], sigma = 0, psiBranch = 0, 
                                psiTrait = 0, z = 10, phi = 0, traitOpt = 1, br0 = 0.1, br_td = -0.1,
                                nTdim = 2, root.state = 1, root.trait.state = 0, plotit = FALSE,
                                keepExtinct = FALSE)[[1]]
  sim0.55[[i]] <- DAMOCLES_sim(trsYule, gamma_0 = gammas[i], gamma_td = 0, mu = mus[4], sigma = 0, psiBranch = 0, 
                                psiTrait = 0, z = 10, phi = 0, traitOpt = 1, br0 = 0.1, br_td = -0.1,
                                nTdim = 2, root.state = 1, root.trait.state = 0, plotit = FALSE,
                                keepExtinct = FALSE)[[1]]
  sim0.7[[i]] <- DAMOCLES_sim(trsYule, gamma_0 = gammas[i], gamma_td = 0, mu = mus[5], sigma = 0, psiBranch = 0, 
                                psiTrait = 0, z = 10, phi = 0, traitOpt = 1, br0 = 0.1, br_td = -0.1,
                                nTdim = 2, root.state = 1, root.trait.state = 0, plotit = FALSE,
                                keepExtinct = FALSE)[[1]]
  sim0.85[[i]] <- DAMOCLES_sim(trsYule, gamma_0 = gammas[i], gamma_td = 0, mu = mus[6], sigma = 0, psiBranch = 0, 
                                psiTrait = 0, z = 10, phi = 0, traitOpt = 1, br0 = 0.1, br_td = -0.1,
                                nTdim = 2, root.state = 1, root.trait.state = 0, plotit = FALSE,
                                keepExtinct = FALSE)[[1]]
  sim1.0[[i]] <- DAMOCLES_sim(trsYule, gamma_0 = gammas[i], gamma_td = 0, mu = mus[7], sigma = 0, psiBranch = 0, 
                                psiTrait = 0, z = 10, phi = 0, traitOpt = 1, br0 = 0.1, br_td = -0.1,
                                nTdim = 2, root.state = 1, root.trait.state = 0, plotit = FALSE,
                                keepExtinct = FALSE)[[1]]
}


pa1 = matrix(c(phy$tip.label, patable$state), nrow = length(phy$tip.label), ncol = 2)
## Presence absence matrices under know parameters
PAs0.1 = list()
PAs0.25 = list()
PAs0.4 = list()
PAs0.55 = list()
PAs0.7 = list()
PAs0.85 = list()
PAs1.0 = list()

for(i in 1:length(sim0.1)){
  PAs0.1[[i]] <- matrix(c(trsYule$tip.label, sim0.1[[i]]$state), nrow = length(trsYule$tip.label), ncol = 2)
  PAs0.25[[i]] <- matrix(c(trsYule$tip.label, sim0.25[[i]]$state), nrow = length(trsYule$tip.label), ncol = 2)
  PAs0.4[[i]] <- matrix(c(trsYule$tip.label, sim0.4[[i]]$state), nrow = length(trsYule$tip.label), ncol = 2)
  PAs0.55[[i]] <- matrix(c(trsYule$tip.label, sim0.55[[i]]$state), nrow = length(trsYule$tip.label), ncol = 2)
  PAs0.7[[i]] <- matrix(c(trsYule$tip.label, sim0.7[[i]]$state), nrow = length(trsYule$tip.label), ncol = 2)
  PAs0.85[[i]] <- matrix(c(trsYule$tip.label, sim0.85[[i]]$state), nrow = length(trsYule$tip.label), ncol = 2)
  PAs1.0[[i]] <- matrix(c(trsYule$tip.label, sim1.0[[i]]$state), nrow = length(trsYule$tip.label), ncol = 2)
}

## Estimate parameters of colonization and local extinction based on know paraterms
out0.1 <- list()
out0.25 <- list()
out0.4 <- list()
out0.55 <- list()
out0.7 <- list()
out0.85 <- list()
out1.0 <- list()

for(i in 1:length(PAs0.1)){
  svMisc::progress(i, max.value = length(PAs0.1))
  out0.1[[i]] = DAMOCLES_bootstrap(phy = trsYule, pa = PAs0.1[[i]], initparsopt = c(0.01, 1.8), 
                           idparsopt = c(1, 2), pars2 = c(1E-3, 1E-4, 1E-5, 10000), 
                           pchoice = 1, runs = 100, estimate_pars = TRUE,
                           conf.int = 0.95)[[1]]
  out0.25[[i]] = DAMOCLES_bootstrap(phy = trsYule, pa = PAs0.25[[i]], initparsopt = c(0.01, 1.8), 
                                  idparsopt = c(1, 2), pars2 = c(1E-3, 1E-4, 1E-5, 10000), 
                                  pchoice = 1, runs = 100, estimate_pars = TRUE,
                                  conf.int = 0.95)[[1]]
  out0.4[[i]] = DAMOCLES_bootstrap(phy = trsYule, pa = PAs0.4[[i]], initparsopt = c(0.01, 1.8), 
                                  idparsopt = c(1, 2), pars2 = c(1E-3, 1E-4, 1E-5, 10000), 
                                  pchoice = 1, runs = 100, estimate_pars = TRUE,
                                  conf.int = 0.95)[[1]]
  out0.55[[i]] = DAMOCLES_bootstrap(phy = trsYule, pa = PAs0.55[[i]], initparsopt = c(0.01, 1.8), 
                                  idparsopt = c(1, 2), pars2 = c(1E-3, 1E-4, 1E-5, 10000), 
                                  pchoice = 1, runs = 100, estimate_pars = TRUE,
                                  conf.int = 0.95)[[1]]
  out0.7[[i]] = DAMOCLES_bootstrap(phy = trsYule, pa = PAs0.7[[i]], initparsopt = c(0.01, 1.8), 
                                  idparsopt = c(1, 2), pars2 = c(1E-3, 1E-4, 1E-5, 10000), 
                                  pchoice = 1, runs = 100, estimate_pars = TRUE,
                                  conf.int = 0.95)[[1]]
  out0.85[[i]] = DAMOCLES_bootstrap(phy = trsYule, pa = PAs0.85[[i]], initparsopt = c(0.01, 1.8), 
                                  idparsopt = c(1, 2), pars2 = c(1E-3, 1E-4, 1E-5, 10000), 
                                  pchoice = 1, runs = 100, estimate_pars = TRUE,
                                  conf.int = 0.95)[[1]]
  out1.0[[i]] = DAMOCLES_bootstrap(phy = trsYule, pa = PAs1.0[[i]], initparsopt = c(0.01, 1.8), 
                                  idparsopt = c(1, 2), pars2 = c(1E-3, 1E-4, 1E-5, 10000), 
                                  pchoice = 1, runs = 100, estimate_pars = TRUE,
                                  conf.int = 0.95)[[1]]
}

## Save results
save.image("DamoclesSimulations")

## Extract results
load("DamoclesSimulations")

out0.1par <- list()
out0.25par <- list()
out0.4par <- list()
out0.55par <- list()
out0.7par <- list()
out0.85par <- list()
out1.0par <- list()

for(i in 1:length(out0.1)){
  out0.1par <- out0.1[[i]][2:3, ]
  out0.25par <- out0.25[[i]][2:3, ]
  out0.4par <- out0.4[[i]][2:3, ]
  out0.55par <- out0.55[[i]][2:3, ]
  out0.7par <- out0.7[[i]][2:3, ]
  out0.85par <- out0.85[[i]][2:3, ]
  out1.0par <- out1.0[[i]][2:3, ] 
}
## Make final tables
library(tidyr)

tab0.1 <- list()
tab0.25 <- list()
tab0.4 <- list()
tab0.55 <- list()
tab0.7 <- list()
tab0.85 <- list()
tab1.0 <- list()

for(i in 1:length(out0.1)){
  tab0.1 <- separate(data = out0.1par[[i]], col = value, into = c("OBS", "inter75"), sep = "\\(")
  tab0.25 <- separate(data = out0.25par[[i]], col = value, into = c("OBS", "inter75"), sep = "\\(")
  tab0.4 <- separate(data = out0.4par[[i]], col = value, into = c("OBS", "inter75"), sep = "\\(")
  tab0.55 <- separate(data = out0.55par[[i]], col = value, into = c("OBS", "inter75"), sep = "\\(")
  tab0.7 <- separate(data = out0.7par[[i]], col = value, into = c("OBS", "inter75"), sep = "\\(")
  tab0.85 <- separate(data = out0.85par[[i]], col = value, into = c("OBS", "inter75"), sep = "\\(")
  tab1.0 <- separate(data = out1.0par[[i]], col = value, into = c("OBS", "inter75"), sep = "\\(")
}

tab0.1 <- do.call(rbind, tab0.1)
tab0.25 <- do.call(rbind, tab0.25)
tab0.4 <- do.call(rbind, tab0.4)
tab0.55 <- do.call(rbind, tab0.55)
tab0.7 <- do.call(rbind, tab0.7)
tab0.85 <- do.call(rbind, tab0.85)
tab1.0 <- do.call(rbind, tab1.0)

## Merge all tables
tablesParameters <- rbind(tab0.1, tab0.25, tab0.4, tab0.55, tab0.7, tab0.85, tab1.0)
tablesParameters$Mu <- rep(mus, each = 7)
tablesParameters$Gamma <- rep(gammas, each = 1)

head(tablesParameters)

write.csv(tablesParameters, "tableDAMOCLESsimulation.csv")

