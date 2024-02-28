# Load library
library(Matrix) # tidyverse not loading without it
library(tidyverse) 
library(cowplot)
library(colorspace)
library(brms)
library(rstan)
library(ggeffects)
library(ape)
library(phytools)
library(diversitree)
library(corHMM)

# Functions
mytheme <- function() {
  theme_bw() + 
    theme(panel.border          = element_rect(fill = NA, colour = "black"), # set border around plot.
          panel.grid.major      = element_blank(), # remove major grid lines
          panel.grid.minor      = element_blank(), # remove minor grid lines
          axis.line             = element_blank(), # remove axis lines
          axis.ticks            = element_line(colour = "black"),
          axis.text             = element_text(size = 10, colour = "black"), # axis text size
          axis.title            = element_text(size = 10), # axis title size
          axis.title.y          = element_text(vjust = 3), # increase distance from the y-axis
          axis.title.x          = element_text(vjust = -1), # increase distance from the x-axis
          panel.background      = element_rect(fill = NA),
          plot.background       = element_rect(fill = NA, color = NA), # remove background colour
          plot.margin           = unit(c(0.2, 0.2, 0.2, 0.2), units = , "cm"), 
          legend.background     = element_rect(fill = NA, color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = NA, color = NA), # get rid of legend panel bg
          strip.text.x          = element_text(size = 10, color = "black", face = "bold"), # for facet plots
          strip.background      = element_rect(fill = NA, color = NA)
    )
} # set up plot theme

# Set directory
setwd('/Users/nicholaswu/Dropbox/Eye malar study')
setwd('C:/Users/nicho/Dropbox (Personal)/Eye malar study')

clean_dat <- read.csv("malar_snake_phylogeny.csv") %>%
  dplyr::mutate(active_mode = factor(active_mode), 
                malar       = ifelse(malar == "Yes",1,0),
                active_bin  = case_when(
                  active_mode %in% "nocturnal" ~ 0,
                  active_mode %in% "underground" ~ 0,
                  active_mode %in% "both" ~ 0.5,
                  active_mode %in% "diurnal" ~ 1),
                state       = case_when(
                  malar %in% 0 & active_mode %in% "nocturnal" ~ 1, # state 1 = no malar and nocturnal
                  malar %in% 0 & active_mode %in% "underground" ~ 1, # state 1 = no malar and nocturnal
                  malar %in% 1 & active_mode %in% "nocturnal" ~ 2, # state 2 = malar and nocturnal
                  malar %in% 0 & active_mode %in% "diurnal" ~ 3, # state 3 = no malar and diurnal
                  malar %in% 1 & active_mode %in% "diurnal" ~ 4)) %>% # state 4 = malar and diurnal
  dplyr::filter(!is.na(active_mode)) %>%
  drop_na(active_mode)

## DATA SUMMARY ##--------------------------------------------------------------------------------------------------------
clean_dat %>%
  dplyr::group_by(malar) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n) * 100)

clean_phylo_dat %>%
  dplyr::group_by(active_mode, malar) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n) * 100) %>%
  dplyr::ungroup() %>%
  tidyr::complete(active_mode, malar, fill = list(n = 0, freq = 0)) %>% # Turns implicit missing values into explicit missing values
  data.frame() 

## SNAKE PHYLOGENY ##-----------------------------------------------------------------------------------------------------
phylo_tree     <- ape::read.tree("squam_shl_new.tre") # Tree with molecular supermatrix 
tree_tip_label <- phylo_tree$tip.label # extract tree tip names
species_list   <- clean_dat$species_tree # extract species name from diet analysis

pruned_tree     <- ape::drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, species_list)) # prune phylo_tree to keep species from the diet analysis
clean_phylo_dat <- subset(clean_dat, species_tree %in% pruned_tree$tip.label) # remove species from dataset not available in the pruned_tree
clean_phylo_dat$species_tree <- factor(clean_phylo_dat$species_tree, levels = c(tree_tip_label))

# check if lengths match for both data and tree
length(unique(pruned_tree$tip.label)) # 1535
length(unique(clean_phylo_dat$species_tree)) # 1535

# reorder factor levels and add species names to rowname
clean_phylo_dat$species_tree <- factor(clean_phylo_dat$species_tree, levels = pruned_tree$tip.label)
rownames(clean_phylo_dat)    <- NULL
trait_dat                    <- tibble::column_to_rownames(clean_phylo_dat, var = 'species_tree')

#plot.phylo(pruned_tree, show.tip.label = FALSE)
ape::is.binary(pruned_tree) # check if polytomies exist

# Time-calibrated tree via Penalised Likelihood (http://phylobotanist.blogspot.com/2018/04/time-calibrated-or-at-least-ultrametric.html)
tree_calibration <- ape::makeChronosCalib(pruned_tree, node = "root", age.min = 81.1, age.max = 137.0) # set 95% HPDI estimated root age (myr) from Tonini et al 2016

# Run all 3 clock models
timetree_corr  <- ape::chronos(pruned_tree, lambda = 1, model = "correlated", calibration = tree_calibration, control = chronos.control()) # adjacent parts of the phylogeny are not allowed to evolve at rates that are very different
timetree_disc  <- ape::chronos(pruned_tree, lambda = 1, model = "discrete", calibration = tree_calibration, control = chronos.control(nb.rate.cat = 1)) # models different parts of the tree as evolving at different rates
timetree_relax <- ape::chronos(pruned_tree, lambda = 1, model = "relaxed", calibration = tree_calibration, control = chronos.control()) # rates to vary freely

# log-likelihood higher is better, PHIIC lower is better
# timetree_corr = log-Lik = -358.4109, PHIIC = 9920.83
# timetree_disc = log-Lik = -363.7136, PHIIC = 3797.43 
# timetree_relax = log-Lik = -363.676, PHIIC = 9954.26 

pruned_tree_cal <- timetree_corr # correlated model best

#saveRDS(pruned_tree_cal, "pruned_tree_cal.rds")
pruned_tree_cal  <- readRDS("pruned_tree_cal.rds")

plot(pruned_tree_cal, cex = 0.3)
axisPhylo()

# check tree is ultrametric
is.ultrametric(pruned_tree_cal) # TRUE

# Create correlation matrix for analysis
phylo_cor <- vcv(pruned_tree_cal, cor = T)
#saveRDS(phylo_cor, "phylo_cor.rds")
phylo_cor <- readRDS("phylo_cor.rds")

# Pearson correlation of subsetted data
family_n_subset <- clean_phylo_dat %>% dplyr::group_by(family) %>% summarise(n_subset = n())
family_n_full   <- clean_dat %>% dplyr::group_by(family) %>% summarise(n_full = n())
family_n        <- merge(family_n_full, family_n_subset, by = "family", all = TRUE)
family_cor      <- cor.test(family_n$n_full, family_n$n_subset, method = "pearson", use = "complete.obs")
family_cor$estimate
family_cor$p.value

as.data.frame(family_n %>% group_by(family) %>% summarise(n_subset / n_full * 100)) # proportion of species represented

ggplot(family_n, aes(x = log(n_full + 1), y = log(n_subset + 1))) + geom_point()+ mytheme() + xlab("Species per family (ln(n); full)") + ylab("Species per family (ln(n); strict)")

## PHYLOGENETIC SIGNAL AND ANCESTRAL STATE RECONSTRUCTION ## ---------------------------------------------------------------------------------------
## ASR of eye malar ##
malar_matrix  <- setNames(trait_dat$malar, rownames(trait_dat))
active_matrix <- setNames(trait_dat$active_bin, rownames(trait_dat))

phytools::phylosig(pruned_tree, malar_matrix, method = "lambda", test = TRUE, nsim = 1000)
phytools::phylosig(pruned_tree, active_matrix, method = "lambda", test = TRUE, nsim = 1000)

# Plot ASR
active_obj <- phytools::contMap(pruned_tree_cal, active_matrix, plot = FALSE)
active_obj <- phytools::setMap(active_obj,  c("black","grey", "#879FDB"))

# Fig 1a - internal (original pruned_tree - no ultrametric)
plot(active_obj, type = "fan", legend = 0.7 * max(nodeHeights(pruned_tree_cal)), lwd = 0.8, outline = FALSE, ftype = "off")


col_range <- colorRampPalette(c("white", "#002F70"))
malar_col <- col_range(length(unique(trait_dat$malar)))
names(malar_col) <- levels(factor(trait_dat$malar))

# Fig 3a
#active_col <- viridis::viridis(4) # set 4 discrete colours
active_col <- c("#5F1415","grey", "#002F70")
names(active_col)  <- levels(factor(trait_dat$active_bin)) # set 4 discrete colours

diversitree::trait.plot(pruned_tree_cal, trait_dat,
                        cols = list(malar = malar_col, active_bin = active_col),
                        cex.lab = 0.00001, w = 0.06)
plot(active_obj$tree, 
     colors  = active_obj$cols, 
     type    = "fan",
     add     = TRUE, 
     ftype   = "off", 
     lwd     = 0.5,
     outline = FALSE,
     xlim    = get("last_plot.phylo", envir =.PlotPhyloEnv)$x.lim,
     ylim    = get("last_plot.phylo", envir =.PlotPhyloEnv)$y.lim)

## MALAR ANALYSIS ## ---------------------------------------------------------------------------------------------------------------------------
# Set options in Rstan
rstan::rstan_options(auto_write = TRUE) # translate to STAN platform for running Bayesian model
options(mc.cores = parallel::detectCores()) # detects how many cores available to use
set.seed(10)

malar_phylo_model <- brms::brm(malar ~ -1 + active_mode + (1 | gr(species_tree, cov = A)), 
                               data    = clean_phylo_dat, 
                               family  = 'bernoulli',
                               prior   = prior(normal(0, 3), "b"),
                               data2   = list(A = phylo_cor),
                               iter    = 5e3, warmup = 2.5e3, chains = 4, cores = 4,
                               control = list(adapt_delta = 0.99, max_treedepth = 15))

summary(malar_phylo_model)
brms::pp_check(malar_phylo_model)

#saveRDS(malar_phylo_model, "malar_phylo_model.rds")
malar_phylo_model <- readRDS("malar_phylo_model.rds")

# Phylogenetic signal via hypothesis method
#parnames(diet_phylo_model)

(hyp <- hypothesis(malar_phylo_model, "sd_species_tree__Intercept^2 / (sd_species_tree__Intercept^2 + (3.141593^2/3)) = 0", class = NULL))

# PLOT FIG 3b
# Extract estimates and 95% CI
malar_eff <- data.frame(ggeffects::ggpredict(malar_phylo_model, "active_mode")) %>%
  dplyr::rename(active_mode = x) %>%
  dplyr::mutate(active_mode = forcats::fct_reorder(active_mode, predicted),
                predicted = predicted * 100,
                conf.low  = conf.low * 100,
                conf.high = conf.high * 100)

malar_eff %>%
  ggplot(aes(x = reorder(active_mode, predicted), y = predicted, colour = active_mode)) +
  geom_point(size = 2, show.legend = FALSE) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), size = 1, width = 0.1, show.legend = FALSE) + 
  scale_colour_manual(values = c("black", "black", "grey", "#879FDB")) +
  xlab(NULL) + ylab("Probability of malar stripe (%)") +
  mytheme() +
  coord_flip()

## ANCESTRAL STATE INFERENCE FOR EYE MALAR VIA corHMM2.1 - ALL ## --------------------------------------------------------------------------
malar_dat <- clean_phylo_dat %>% dplyr::select(species_tree, malar)
geiger::name.check(pruned_tree$tree.label, malar_dat$species_tree, data.names=NULL) #using geiger

## ER: Equal rates model ##
malar_ER <- corHMM::corHMM(pruned_tree_cal, malar_dat, model = c("ER"), 
                           node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100)
#saveRDS(malar_ER, "malar_ER.rds")
malar_ER <- readRDS("malar_ER.rds")

# Visualise model
corHMM::plotMKmodel(malar_ER) 

## ARD: All rates different model ##
malar_ARD <- corHMM::corHMM(pruned_tree_cal, malar_dat, model = c("ARD"),
                            node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100) 
#saveRDS(malar_ARD, "malar_ARD.rds")
malar_ARD <- readRDS("malar_ARD.rds")

# Plot tree
malar_ARD_model       <- malar_ARD$solution
malar_ARD_model[is.na(malar_ARD_model)] <- 0
diag(malar_ARD_model) <- -rowSums(malar_ARD_model)
malar_ARD_simmap      <- corHMM::makeSimmap(tree = pruned_tree_cal, data = malar_dat, model = malar_ARD_model, rate.cat = 1, nSim = 1, nCores = 1)
malar_ARD_cols        <- setNames(c("light grey", "black"), c("1","2"))
phytools::plotSimmap(malar_ARD_simmap[[1]], fsize = 0.001, lwd = 0.6, type = "fan", colors = malar_ARD_cols)

## Model ARD/ARD: Hidden Markov model using ARD model in both matrices ##
malar_ARD_ARD <- corHMM::corHMM(pruned_tree_cal, malar_dat, model = c("ARD"), 
                                node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 2, get.tip.states = TRUE, nstarts = 100)
#saveRDS(malar_ARD_ARD, "malar_ARD_ARD.rds")
malar_ARD_ARD <- readRDS("malar_ARD_ARD.rds")

## Model ER/ARD: Hidden Markov model using ARD model and ER model as the two rate matrices ##
# Construct two 'within' rate category rate.mat objects (R1 and R) - state-dependent processes
# Used to index the unique parameters to be estimated
malar_R1_ER <- corHMM::getStateMat4Dat(malar_dat, model = "ER")$rate.mat # R1 assume a drift-like hypothesis where all transition rates are equal
malar_R2_ER <- corHMM::getStateMat4Dat(malar_dat, model = "ER")$rate.mat # R2
malar_R2_ARD <- corHMM::getStateMat4Dat(malar_dat, model = "ARD")$rate.mat # R2 assume differences in transition rate

# getRateCatMat function to generate matrix for transitions 'among' the different rate classes - parameter process
RateClassMat  <- corHMM::getRateCatMat(2)

# Create list where first element corresponds to R1 and the second to R2 to create full model using getFullMat
malar_FullMat <- corHMM::getFullMat(list(malar_R1_ER, malar_R2_ARD), RateClassMat)

# Visualize model setup
#corHMM::plotMKmodel(ophio_FullMat, rate.cat = 2, text.scale = 0.7) # check that the model is correct

malar_ER_ARD <- corHMM::corHMM(pruned_tree_cal, malar_dat, rate.cat = 2, rate.mat = malar_FullMat, 
                               node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), get.tip.states = TRUE, nstarts = 100) 
#saveRDS(malar_ER_ARD, "malar_ER_ARD.rds")
malar_ER_ARD <- readRDS("malar_ER_ARD.rds")

## Model ER/ER: Hidden Markov model using two ER models as the two rate matrices ##
malar_FullMat2 <- corHMM::getFullMat(list(malar_R1_ER, malar_R2_ER), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal

# Visualize model setup
#corHMM::plotMKmodel(ophio_FullMat2, rate.cat = 2, text.scale = 0.7)

malar_ER_ER <- corHMM::corHMM(pruned_tree_cal, malar_dat, rate.cat = 2, rate.mat = malar_FullMat2, 
                             node.states  = "marginal", lewis.asc.bias = FALSE,  root.p = c(1,0), get.tip.states = TRUE, nstarts = 100) 
#saveRDS(malar_ER_ER, "malar_ER_ER.rds")
malar_ER_ER <- readRDS("malar_ER_ER.rds")

## ANCESTRAL STATE INFERENCE FOR ACTIVE MODE VIA corHMM2.1 - ALL ## --------------------------------------------------------------------------
active_dat <- clean_phylo_dat %>% dplyr::select(species_tree, active_bin)

## ER: Equal rates model ##
active_ER <- corHMM::corHMM(pruned_tree_cal, active_dat, model = c("ER"), 
                           node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100)
saveRDS(active_ER, "active_ER.rds")
active_ER <- readRDS("active_ER.rds")

## ARD: All rates different model ##
active_ARD <- corHMM::corHMM(pruned_tree_cal, active_dat, model = c("ARD"),
                            node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100) 
#saveRDS(active_ARD, "active_ARD.rds")
active_ARD <- readRDS("active_ARD.rds")

## Model ARD/ARD: Hidden Markov model using ARD model in both matrices ##
active_ARD_ARD <- corHMM::corHMM(pruned_tree_cal, active_dat, model = c("ARD"), 
                                node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 2, get.tip.states = TRUE, nstarts = 100)
#saveRDS(active_ARD_ARD, "active_ARD_ARD.rds")
active_ARD_ARD <- readRDS("active_ARD_ARD.rds")

## Model ER/ARD: Hidden Markov model using ARD model and ER model as the two rate matrices ##
# Construct two 'within' rate category rate.mat objects (R1 and R) - state-dependent processes
# Used to index the unique parameters to be estimated
active_R1_ER  <- corHMM::getStateMat4Dat(active_dat, model = "ER")$rate.mat  # R1
active_R2_ER  <- corHMM::getStateMat4Dat(active_dat, model = "ER")$rate.mat  # R2
active_R2_ARD <- corHMM::getStateMat4Dat(active_dat, model = "ARD")$rate.mat # R2 

# Create list where first element corresponds to R1 and the second to R2 to create full model using getFullMat
active_FullMat <- corHMM::getFullMat(list(active_R1_ER, active_R2_ARD), RateClassMat)
active_ER_ARD <- corHMM::corHMM(pruned_tree_cal, active_dat, rate.cat = 2, rate.mat = active_FullMat, 
                               node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), get.tip.states = TRUE, nstarts = 100) 
#saveRDS(active_ER_ARD, "active_ER_ARD.rds")
active_ER_ARD <- readRDS("active_ER_ARD.rds")

## Model ER/ER: Hidden Markov model using two ER models as the two rate matrices ##
active_FullMat2 <- corHMM::getFullMat(list(active_R1_ER, active_R2_ER), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal
active_ER_ER <- corHMM::corHMM(pruned_tree_cal, active_dat, rate.cat = 2, rate.mat = active_FullMat2, 
                             node.states  = "marginal", lewis.asc.bias = FALSE,  root.p = c(1,0), get.tip.states = TRUE, nstarts = 100) 
saveRDS(active_ER_ER, "active_ER_ER.rds")
active_ER_ER <- readRDS("active_ER_ER.rds")

##Plotting the tree and model
par(mfrow = c(1,2))
# Malar pattern
corHMM::plotRECON(pruned_tree_cal, malar_ER_ER$states, 
                  pie.cex = 0.3, 
                  piecolors = (c("black","orange", "grey","#DF536B")), 
                  show.tip.label = FALSE,
                  cex = 0.25, 
                  adj = 0.25) 
add.scale.bar(pruned_tree_cal)
cols <- setNames(c("grey", "#DF536B"), levels(malar_dat$malar))
malar_dat2 <- malar_dat %>% dplyr::arrange(factor(species_tree, levels = levels(malar_dat$species_tree))) # reorder whole dataset to match phylogeny
ape::tiplabels(pie = to.matrix(malar_dat2$malar, sort(unique(malar_dat2$malar))), piecol = cols, cex = 0.3)

# Active mode
corHMM::plotRECON(pruned_tree_cal, active_ARD_ARD$states, 
                  pie.cex = 0.3, 
                  piecolors = (c("black","grey", "#53d1df", "dark grey", "light grey","#5384df")),
                  show.tip.label = FALSE,
                  cex = 0.25, 
                  adj = 0.25) 
add.scale.bar(pruned_tree_cal)
cols <- setNames(c("dark grey", "light grey", "#5384df"), levels(active_dat$active_bin))
active_dat2 <- active_dat %>% dplyr::arrange(factor(species_tree, levels = levels(active_dat$species_tree))) # reorder whole dataset to match phylogeny
ape::tiplabels(pie = to.matrix(active_dat2$active_bin, sort(unique(active_dat2$active_bin))), piecol = cols, cex = 0.3)
par(mfrow = c(1,1))

## ANCESTRAL STATE INFERENCE FOR EYE MALAR AND ACTIVE MODE corHMM2.1 - ALL ## --------------------------------------------------------------------------
malar_active_dat <- clean_phylo_dat %>% dplyr::select(species_tree, malar, active_bin)

## ER: Equal rates model ##
malar_active_ER <- corHMM::corHMM(pruned_tree_cal, malar_active_dat, model = c("ER"), 
                          node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0,0,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100)
#saveRDS(malar_active_ER, "malar_active_ER.rds")
malar_active_ER <- readRDS("malar_active_ER.rds")

#corHMM::plotMKmodel(band_ophio_ER) 

## SYM: All rates different model ##
malar_active_SYM <- corHMM::corHMM(pruned_tree_cal, malar_active_dat, model = c("SYM"),
                                 node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0,0,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100) 
#saveRDS(malar_active_SYM, "malar_active_SYM.rds")
malar_active_SYM <- readRDS("malar_active_SYM.rds")

## ARD: All rates different model ##
malar_active_ARD <- corHMM::corHMM(pruned_tree_cal, malar_active_dat, model = c("ARD"),
                           node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0,0,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100) 
#saveRDS(malar_active_ARD, "malar_active_ARD.rds")
malar_active_ARD <- readRDS("malar_active_ARD.rds")

# Construct two 'within' rate category rate.mat objects (R1 and R) - state-dependent processes
# Used to index the unique parameters to be estimated
malar_active_R1_ER  <- corHMM::getStateMat4Dat(malar_active_dat, model= "ER")$rate.mat # R1 assume a drift-like hypothesis where all transition rates are equal
malar_active_R2_ER  <- corHMM::getStateMat4Dat(malar_active_dat, model= "ER")$rate.mat # R2
malar_active_R1_SYM <- corHMM::getStateMat4Dat(malar_active_dat, model = "SYM")$rate.mat # R1 
malar_active_R2_SYM <- corHMM::getStateMat4Dat(malar_active_dat, model = "SYM")$rate.mat # R2 
malar_active_R1_ARD <- corHMM::getStateMat4Dat(malar_active_dat, model = "ARD")$rate.mat # R1 
malar_active_R2_ARD <- corHMM::getStateMat4Dat(malar_active_dat, model = "ARD")$rate.mat # R2 

## Model ER/ARD: Hidden Markov model using ARD model and ER model as the two rate matrices ##
# Create list where first element corresponds to R1 and the second to R2 to create full model using getFullMat
RateClassMat             <- corHMM::getRateCatMat(2)
RateClassMat1            <- corHMM::equateStateMatPars(RateClassMat, c(1, 2))
malar_active_FullMat_mix <- corHMM::getFullMat(list(malar_active_R1_ER, malar_active_R2_ARD), RateClassMat1)

# Visualize model setup
#corHMM::plotMKmodel(band_ophio_FullMat, rate.cat = 2, text.scale = 0.7) # check that the model is correct

malar_active_ER_ARD <- corHMM::corHMM(pruned_tree_cal, malar_active_dat, rate.cat = 2, rate.mat = malar_active_FullMat_mix, 
                              node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0,0,0), get.tip.states = TRUE, nstarts = 100) 
#saveRDS(malar_active_ER_ARD, "malar_active_ER_ARD.rds")
malar_active_ER_ARD <- readRDS("malar_active_ER_ARD.rds")

## Model ER/ER: Hidden Markov model using two ER models as the two rate matrices ##
malar_active_FullMat_ER <- corHMM::getFullMat(list(malar_active_R1_ER, malar_active_R2_ER), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal
malar_active_ER_ER <- corHMM::corHMM(pruned_tree_cal, malar_active_dat, rate.cat = 2, rate.mat = malar_active_FullMat_ER, 
                             node.states  = "marginal", lewis.asc.bias = FALSE,  root.p = c(1,0,0,0), get.tip.states = TRUE, nstarts = 100) 
#saveRDS(malar_active_ER_ER, "malar_active_ER_ER.rds")
malar_active_ER_ER <- readRDS("malar_active_ER_ER.rds")

## Model SYM/SYM: Hidden Markov model using two ARD models as the two rate matrices ##
malar_active_FullMat_SYM <- corHMM::getFullMat(list(malar_active_R1_SYM, malar_active_R2_SYM), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal
malar_active_SYM_SYM <- corHMM::corHMM(pruned_tree_cal, malar_active_dat, rate.cat = 2, rate.mat = malar_active_FullMat_SYM, 
                                     node.states  = "marginal", lewis.asc.bias = FALSE,  root.p = c(1,0,0,0), get.tip.states = TRUE, nstarts = 100) 
saveRDS(malar_active_SYM_SYM, "malar_active_SYM_SYM.rds")
malar_active_SYM_SYM <- readRDS("malar_active_SYM_SYM.rds")

## Model ARD/ARD: Hidden Markov model using two ARD models as the two rate matrices ##
malar_active_FullMat_ARD <- corHMM::getFullMat(list(malar_active_R1_ARD, malar_active_R2_ARD), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal
malar_active_ARD_ARD <- corHMM::corHMM(pruned_tree_cal, malar_active_dat, rate.cat = 2, rate.mat = malar_active_FullMat_ARD, 
                             node.states  = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0,0,0), get.tip.states = TRUE, nstarts = 100) 
saveRDS(malar_active_ARD_ARD, "malar_active_ARD_ARD.rds")
malar_active_ARD_ARD <- readRDS("malar_active_ARD_ARD.rds")

malar_active_R1_ARD_dual      <- corHMM::getStateMat4Dat(malar_active_dat, model = "ARD", dual = TRUE)$rate.mat #R1
malar_active_R2_ARD_dual      <- corHMM::getStateMat4Dat(malar_active_dat, model = "ARD", dual = TRUE)$rate.mat #R2
malar_active_FullMat_ARD_dual <- corHMM::getFullMat(list(malar_active_R1_ARD_dual, malar_active_R2_ARD_dual), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal
malar_active_ARD_ARD_dual     <- corHMM::corHMM(pruned_tree_cal, malar_active_dat, rate.cat = 2, rate.mat = malar_active_FullMat_ARD_dual, 
                                     node.states  = "marginal", lewis.asc.bias = FALSE,  root.p = c(1,0,0,0), get.tip.states = TRUE, nstarts = 100) 
saveRDS(malar_active_ARD_ARD_dual, "malar_active_ARD_ARD_dual.rds")
malar_active_ARD_ARD_dual <- readRDS("malar_active_ARD_ARD_dual.rds")

corHMM::plotMKmodel(malar_active_ARD_ARD) 


## ALTERNATIVE HYPOTHESIS ## ------------------------------------------------------------------------------------ 
# Frequency of banded stripes with ophiophagy, reptiles diet, mammal diet
as.data.frame(diet_dat %>%
                dplyr::group_by(pattern_type) %>%
                dplyr::summarise(ophio = sum(ophiophagy_bin == 1),
                                 reptile = sum(rept_special == 1),
                                 mammal = sum(mammal == 1)) %>%
                dplyr::mutate(ophio_freq = ophio / sum(ophio) * 100,
                              reptile_freq = reptile / sum(reptile) * 100,
                              mammal_freq = mammal / sum(mammal) * 100) %>%
                dplyr::ungroup() %>%
                tidyr::complete(pattern_type, fill = list(n = 0, freq = 0))) # Turns implicit missing values into explicit missing values) 


              
