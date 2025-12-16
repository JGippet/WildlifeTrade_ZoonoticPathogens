## R code for Gippet et al. (2026): 
## Wildlife trade drives animal-to-human pathogen transmission over 40 years



#############################################################################
###########      Part 0 - Data preparation & Core functions     ############# 
#############################################################################

# Load packages
library(ComplexHeatmap)
library(DHARMa)
library(dplyr)
library(effectsize)
library(eulerr)
library(ggeffects)
library(ggmosaic)
library(ggnewscale)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(ggtreeExtra)
library(glmmTMB)
library(grid)
library(gridExtra)
library(multcomp)
library(performance)
library(phytools)
library(piecewiseSEM)
library(purrr)
library(PVR)
library(rlang)
library(sjPlot)
library(tibble)
library(tidyr)
library(treeio)


####>>  Functions

### Cross-referencing to  MDD
crossRefTaxa <- function(df_toCR= NULL, df_base=NULL){
  print(cat(paste0(dim(df_toCR)[1], " species in the dataframe to cross-reference \n")))
  df_toCR$MDD_id <- NA
  
  df_toCR$sciName = stringr::str_replace_all(tolower(df_toCR$sciName), " ", "_")
  df_base$sciName = stringr::str_replace_all(tolower(df_base$sciName), " ", "_")

  df_base <- df_base %>% 
    mutate(
      PossibleNames_MDD = purrr::map(PossibleNames_MDD, ~ {
        x <- .x
        
        # if it's already a length>1 character vector, keep it
        if (length(x) > 1) {
          real_names <- x
          
          # if it's a single string that *looks* like "c(...)" then parse it
        } else if (str_detect(x, "^c\\(")) {
          # remove leading "c(" and trailing ")"
          inner <- str_remove_all(x, "^c\\(|\\)$")
          # split on comma, then strip off any quotes and whitespace
          real_names <- str_split(inner, ",\\s*")[[1]] %>% 
            str_remove_all('^"|"$')
          
          # else it's a plain single name
        } else {
          real_names <- x
        }
        
        # now you have a true chr vector; lower‐case + replace spaces
        real_names %>% 
          tolower() %>% 
          str_replace_all(" ", "_")
      })
    )
  
  df_base <- df_base %>% 
    mutate(
      all_possibleNames = purrr::map(all_possibleNames, ~ {
        x <- .x
        
        # if it's already a length>1 character vector, keep it
        if (length(x) > 1) {
          real_names <- x
          
          # if it's a single string that *looks* like "c(...)" then parse it
        } else if (str_detect(x, "^c\\(")) {
          # remove leading "c(" and trailing ")"
          inner <- str_remove_all(x, "^c\\(|\\)$")
          # split on comma, then strip off any quotes and whitespace
          real_names <- str_split(inner, ",\\s*")[[1]] %>% 
            str_remove_all('^"|"$')
          
          # else it's a plain single name
        } else {
          real_names <- x
        }
        
        # now you have a true chr vector; lower‐case + replace spaces
        real_names %>% 
          tolower() %>% 
          str_replace_all(" ", "_")
      })
    )
  
  
  # first round
  for (i in 1:dim(df_toCR)[1]){
    if(df_toCR$sciName[i] %in% df_base$sciName){
      df_toCR$MDD_id[i] <- df_base$MDD_id[which(df_base$sciName==df_toCR$sciName[i])]
    }
  }
  
  print(cat(paste0(sum(is.na(df_toCR$MDD_id)), " species not matched after first round \n"))) 
  
  # second round
  for (i in 1:dim(df_toCR)[1]){
    if(is.na(df_toCR$MDD_id[i])){
      listIDs <- NULL
      for (p in 1:dim(df_base)[1]){
        if(df_toCR$sciName[i] %in% df_base$PossibleNames_MDD[p][[1]]){
          listIDs <- c(listIDs, df_base$MDD_id[p])
        }
      }
      df_toCR$MDD_id[i] <- list(listIDs)
    }
    svMisc::progress(i, dim(df_toCR)[1])
  }
  
  w = 0
  for (i in 1:dim(df_toCR)[1]){
    if(is.null(df_toCR$MDD_id[i][[1]])){
      w=w+1
    }
  }
  
  print(cat(paste0(w, " species not matched after second round \n"))) 
  
  
  # third round
  for (i in 1:dim(df_toCR)[1]){
    if(is.null(df_toCR$MDD_id[i][[1]])){
      listIDs <- NULL
      for (p in 1:dim(df_base)[1]){
        if(df_toCR$sciName[i] %in% df_base$all_possibleNames[p][[1]]){
          listIDs <- c(listIDs, df_base$MDD_id[p])
        }
      }
      df_toCR$MDD_id[i] <- list(listIDs)
    }
    svMisc::progress(i, dim(df_toCR)[1])
  }
  
  w = 0
  for (i in 1:dim(df_toCR)[1]){
    if(is.null(df_toCR$MDD_id[i][[1]])){
      w=w+1
    }
  }
  
  print(cat(paste0(w, " species not matched after third round \n"))) 
  
  df_toCR$sciName = stringr::str_replace_all(tolower(df_toCR$sciName), "_", " ")
  return(df_toCR)
}


### Dropping tree tips for multiPhylo object
# taken here: http://blog.phytools.org/2020/06/pruning-tips-from-multiphylo-object.html
drop.tip.multiPhylo_jg <- function(phy, tip, ...){
  if(!inherits(phy,"multiPhylo"))
    stop("phy is not an object of class \"multiPhylo\".")
  else {
    trees<-lapply(phy, ape::drop.tip,tip=tip,...)
    class(trees)<-"multiPhylo"
  }
  trees
}


### Pruning phylogenetic trees from invalid and duplicated taxa and changing tips labels to MDD_id
ProcessTree <- function(tree= NULL, DB_ref= dataURall_generalTrends, tips_to_rm= NULL, checked_names= NULL){
  if(is.null(tips_to_rm)){
    # finding invalid taxa
    phylo_species <- tree[1][[1]]$tip.label %>% as.data.frame
    colnames(phylo_species) <- "species"
    phylo_species$sciName <- tolower(stringr::str_replace_all(phylo_species$species, "_", " "))
    # Cross-referencing to MDD
    phylo_species_CrossRef <- crossRefTaxa(df_toCR = phylo_species, df_base = DB_ref)
    # checking correspondences
    phylo_species_CrossRef$MDD_id_check1 <- NA
    for (i in 1:dim(phylo_species_CrossRef)[1]){
      phylo_species_CrossRef$MDD_id_check1[i] <- length(phylo_species_CrossRef$MDD_id[i][[1]])
    }
    print(phylo_species_CrossRef$MDD_id_check1 %>% table(useNA="always"))
    
    # finding duplicated taxa
    print(phylo_species_CrossRef %>% filter(MDD_id_check1==1) %>% dplyr::select("MDD_id") %>% distinct() %>% dim)
    phylo_species_CrossRef_v2 <- phylo_species_CrossRef %>% filter(MDD_id_check1==1) %>% filter(!duplicated(MDD_id))
    phylo_species_CrossRef_v2$MDD_id <- unlist(phylo_species_CrossRef_v2$MDD_id)
    
    # removing invalid tips and duplicates from input tree
    tips_to_rm <- phylo_species_CrossRef$species[which(phylo_species_CrossRef$species %notin% phylo_species_CrossRef_v2$species)]
    tree_prun <- drop.tip.multiPhylo_jg(phy=tree, tip=tips_to_rm)
    
    # assigning MDD_id as tip label
    for(t in 1:length(tree_prun)){
      print(t)
      for (i in 1:length(tree_prun[t][[1]]$tip.label)){
        tree_prun[t][[1]]$tip.label[i] <- phylo_species_CrossRef_v2$MDD_id[which(phylo_species_CrossRef_v2$species==tree_prun[t][[1]]$tip.label[i])]
      }
    }
  }
  return(tree_prun)
}


### Helper that builds one tree for a given column name
make_tree <- function(var, tree = treeData_all,
                      palette_discrete = c("no" = "grey80", "yes" = "black"),
                      add_new_fill = FALSE) {
  # Build once to expose fortified data columns
  p0 <- ggtree(tree,
               layout = "circular",
               branch.length = "branch.length",
               size = 0.2)
  
  # Validate the column against fortified data
  if (!var %in% names(p0$data)) {
    stop(sprintf("Column '%s' not found in ggtree(tree)$data. Available: %s",
                 var, paste(names(p0$data), collapse = ", ")))
  }
  
  vec <- p0$data[[var]]
  is_cont <- is.numeric(vec)
  
  # Map color by var on branches
  p <- p0 + geom_tree(aes(color = !!sym(var)), size = 0.2)
  
  if (is_cont) {
    p <- p + scale_color_viridis_c(direction = -1, option = "D", end = 0.95, begin = 0)
  } else {
    lv  <- if (is.factor(vec)) levels(vec) else sort(unique(as.character(vec)))
    pal <- if (length(lv) <= length(palette_discrete)) {
      x <- palette_discrete[lv]; names(x) <- lv; x
    } else {
      setNames(rep_len(gray.colors(max(3, length(lv))), length(lv)), lv)
    }
    p <- p + scale_color_manual(values = pal, na.translate = TRUE)
  }
  
  p <- p +
    ggtitle(var) +
    theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          legend.title = element_blank())
  
  if (add_new_fill) p <- p + ggnewscale::new_scale_fill()
  
  p
}




####>> Data preparation

custom_colors <- c( 
  "timeTraded_all_scale"= "#8e0062ff",
  "Synanthropy_level_simpleyes"="#2cacdfff", 
  "WildMeatyes" = "#e3d585ff",  
  "Traded_allyes" = "black",
  "Traded_legalyes" = "#adc178ff", 
  "Traded_illegalyes"= "#38a3a5ff", 
  "pubs_X_scale"= "grey", 
  "Traded_all_productyes" = "#7e4618ff", 
  "Traded_all_liveyes"= "#ffb057ff")


Dataset_complete_main <- readRDS("Dataset_complete.RDS") %>% as.data.frame

Dataset_complete_main %>% dim # 6596  90
Dataset_complete_main %>% colnames

Dataset_complete_main <- Dataset_complete_main %>% 
  mutate(order = as.factor(order),
         family = as.factor(family),
         genus = as.factor(genus),
         BioRealm = as.factor(BioRealm),
         CITES_appendix = as.factor(CITES_appendix),
         pubs_X_scale = as.numeric(scale(pubs_X)),
         PhyloDist_Homosapiens_mean_scale= as.numeric(scale(PhyloDist_Homosapiens_mean)),
         Traded_all = as.factor(Traded_all),
         Traded_legal = as.factor(Traded_legal),
         Traded_all_CITES = as.factor(Traded_all_CITES),
         Traded_all_LEMIS2025 = as.factor(Traded_all_LEMIS2025),
         Traded_illegal = as.factor(Traded_illegal),
         Traded_all_live = as.factor(Traded_all_live),
         Traded_all_product = as.factor(Traded_all_product),
         Traded_live_CITES = as.factor(Traded_live_CITES),
         Traded_product_CITES = as.factor(Traded_product_CITES),
         Traded_live_LEMIS2025 = as.factor(Traded_live_LEMIS2025),
         Traded_product_LEMIS2025 = as.factor(Traded_product_LEMIS2025),
         Synanthropy_level_simple = as.factor(Synanthropy_level_simple),
         WildMeat = as.factor(WildMeat),
         Zoonose_vector_1 = as.factor(Zoonose_vector_1),
         Zoonose_virus_vector_1 = as.factor(Zoonose_virus_vector_1),
         noEID2_Zoonose_vector_1 = as.factor(noEID2_Zoonose_vector_1),
         VIRION_Zoonose_vector_1 = as.factor(VIRION_Zoonose_vector_1),
         noPREDICT_VIRION_Zoonose_vector_1 = as.factor(noPREDICT_VIRION_Zoonose_vector_1)
  )

Dataset_complete_main %>% dim # 6596  92



####>> loading phylogenetic tree
# Phylogenetic trees from Vertlife (https://vertlife.org/data/mammals/)
# Download the tree from VertLife and upload it here
treesMammals <- ape::read.nexus("output.nex") # 100 trees
treesMammals

#  Pruning phylogenetic trees from invalid and duplicated taxa and changing tips labels to MDD_id
treesMammals_cleaned <- ProcessTree(tree= treesMammals, DB_ref= Dataset_complete_main) # takes a few minutes

#  the tree has already been cleaned and pruned for the all dataset
treesALLMammals_cleaned_pruned <- treesMammals_cleaned[1:20] # to limit computing time and RAM overload

#  consensus tree and rooting
system.time(treesALLMammals_cleaned_pruned_c <- phytools::consensus.edges(treesALLMammals_cleaned_pruned, method="least.squares")) # takes a ~10 minutes
treesALLMammals_cleaned_pruned_c2 <- phytools::midpoint.root(treesALLMammals_cleaned_pruned_c) # 5658 tips (species)

# Make the consensus tree ultrametric
treesALLMammals_cleaned_pruned_c2_ultra <- phytools::force.ultrametric(treesALLMammals_cleaned_pruned_c2)



####>> Grafting species not in the tree
genera_of_missing_species <- Dataset_complete_main %>% filter(MDD_id %notin% treesALLMammals_cleaned_pruned_c2$tip.label) %>% pull(genus) %>% unique()
genera_of_species_inTree <- Dataset_complete_main %>% filter(MDD_id %in% treesALLMammals_cleaned_pruned_c2$tip.label) %>% pull(genus) %>% unique()

(genera_of_missing_species %in% genera_of_species_inTree) %>% sum
# 318 out of 331 genera of species missing from the tree are in the tree, so easy to graft

dat  <- Dataset_complete_main
treeINIT <- treesALLMammals_cleaned_pruned_c2_ultra

# Maps (character because tip labels are character)
id2bin  <- dat %>% distinct(MDD_id, sciName) %>% 
  transmute(MDD_id = as.character(MDD_id), sciName) %>% tibble::deframe()
bin2id  <- dat %>% distinct(MDD_id, sciName) %>% 
  transmute(sciName, MDD_id = as.character(MDD_id)) %>% tibble::deframe()

# Relabel tips: MDD_id -> sciName
tree_tmp <- treeINIT
tree_tmp$tip.label <- unname(id2bin[treeINIT$tip.label])

# Species to add = in data but not in tree, and whose genus is in the tree
all_binom   <- unique(dat$sciName)
missing_bin <- setdiff(all_binom, tree_tmp$tip.label)
genera_in_tree <- unique(sub("_.*", "", tree_tmp$tip.label))
to_add <- missing_bin[ sub("_.*", "", missing_bin) %in% genera_in_tree ]

# Random within-genus graft
set.seed(42)
tree_aug <- tree_tmp
for (sp in to_add) tree_aug <- add.species.to.genus(tree_aug, sp)

# Relabel back: sciName -> MDD_id
tree_aug$tip.label <- unname(bin2id[tree_aug$tip.label]) 
tree_aug # 6582 tips



####>> Computing eigenvectors to account for phylogenetic non-independence between mammal species

Dataset_complete_main %>% dim # 6596 

# --- 1) Make tree match your dataset (drop species not in Dataset_complete_main) ---
Dataset_tipset <- unique(Dataset_complete_main$MDD_id)

Dataset_complete_main_TREE <- drop.tip(
  tree_aug,
  setdiff(tree_aug$tip.label, Dataset_tipset)
)

# --- 2) Compute PVR eigenvectors from the pruned tree ---
pvr_01 <- PVRdecomp(Dataset_complete_main_TREE)
EV_01  <- pvr_01@Eigen$vectors
# Annotate rows with species names (same order as the EV rows)
rownames(EV_01) <- pvr_01@phylo$tip.label

eigvals <- pvr_01@Eigen$values
cumsum(eigvals) / sum(eigvals) 
plot(eigvals, type="b", main="Phylogenetic eigenvalues", ylab="Eigenvalue", xlim=c(0,100))
# We retained the first five phylogenetic eigenvectors (PEV1–PEV5), which together explained ~40 % of total phylogenetic variation. This number was chosen to capture the major axes of phylogenetic structure while limiting the number of additional parameters and avoiding overfitting from numerous low-variance eigenvectors

k <- 10  # number of phylogenetic eigenvectors to keep >> 54% of all variance
# --- 3) Prepare eigenvectors as a joinable data.frame (character key) ---
EV_df <- as.data.frame(EV_01[, 1:k, drop = FALSE])
EV_df$MDD_id <- rownames(EV_df)
names(EV_df)[1:k] <- paste0("PEV", 1:k)

EV_df_keyed <- EV_df %>%
  mutate(MDD_id_chr = trimws(as.character(MDD_id))) %>%
  dplyr::select(MDD_id_chr, starts_with("PEV")) %>%
  distinct(MDD_id_chr, .keep_all = TRUE) # ensure unique keys on EV side

# --- 4) Prepare your dataset with a matching (character) key and filter to tips present in EV ---
Dataset_complete_main_keyed <- Dataset_complete_main %>%
  mutate(MDD_id_chr = trimws(as.character(MDD_id))) %>%
  filter(MDD_id_chr %in% rownames(EV_01))

# --- 5) Join eigenvectors to the filtered dataset ---
Dataset_complete_main_2 <- Dataset_complete_main_keyed %>%
  left_join(EV_df_keyed, by = "MDD_id_chr")

# --- 6) Scale the PEV columns (optional but common) ---
pev_cols <- paste0("PEV", 1:k)
Dataset_complete_main_2[pev_cols] <- scale(Dataset_complete_main_2[pev_cols])


Dataset_complete_main_2 %>% dim
Dataset_complete_main_2 %>% colnames


####>> Filter domesticated and extinct mammals (and humans
Dataset_01 <- Dataset_complete_main_2 %>% filter(DOMESTICATED_large=="no" & extinct==0)
Dataset_01 %>% dim # 6446  103
Dataset_01 %>% colnames

#################################################################
###########   Part 1 - Analyses linked to Figure 1    ########### 
#################################################################

###>> Figure 1a

# removing tips of domesticated, extinct and humans and all other species not in Dataset_01
tips_to_remove_1 <- c(Dataset_complete_main_TREE$tip.label[which(Dataset_complete_main_TREE$tip.label %notin% as.character(Dataset_01$MDD_id))] )
treesALLMammals_cleaned_pruned_c3 <- ape::drop.tip(phy=Dataset_complete_main_TREE, tip=tips_to_remove_1) # 6446 tips remaining

# rename nodes labels
treesALLMammals_cleaned_pruned_c3$node.label <- paste0("Node_", seq(1:treesALLMammals_cleaned_pruned_c3$Nnode))

# setting a new column to map only main orders
mainOrders_all <- c(Dataset_01$order %>% table %>% sort(decreasing = T) %>% names)[1:10]
Dataset_01 <- Dataset_01 %>% mutate(order2= ifelse(order %in% mainOrders_all, order, NA))

Dataset_01$nodeType <- "tip"

Dataset_01 %>% dplyr::select(c("MDD_id", 
                               "Traded_all", 
                               "pubs_pubmed_all_log",
                               "pubs_pubmed_all_log_rescaled", 
                               "pubs_ALLdisease_PROP_logit",
                               "pubs_ALLdisease_PROP_logit_rescaled",
                               "Synanthropy_level_simple",
                               "PEV1",
                               "PEV2",
                               "PEV3",
                               "PEV4",
                               "PEV5",
                               "PEV6",
                               "PEV7",
                               "PEV8",
                               "PEV9",
                               "PEV10",
                               "nodeType"))

# add nodes to dataframe
Dataset_01_withNodes <- rbind(Dataset_01[,c("MDD_id", 
                                            "Traded_all", 
                                            "pubs_pubmed_all_log",
                                            "pubs_pubmed_all_log_rescaled", 
                                            "pubs_ALLdisease_PROP_logit",
                                            "pubs_ALLdisease_PROP_logit_rescaled",
                                            "Synanthropy_level_simple",
                                            "PEV1",
                                            "PEV2",
                                            "PEV3",
                                            "PEV4",
                                            "PEV5",
                                            "PEV6",
                                            "PEV7",
                                            "PEV8",
                                            "PEV9",
                                            "PEV10",
                                            "nodeType")], 
                              data.frame(MDD_id= paste0("Node_", 1:treesALLMammals_cleaned_pruned_c3$Nnode), 
                                         Traded_all=NA, 
                                         pubs_pubmed_all_log=NA,
                                         pubs_pubmed_all_log_rescaled=NA, 
                                         pubs_ALLdisease_PROP_logit=NA,
                                         pubs_ALLdisease_PROP_logit_rescaled=NA,
                                         Synanthropy_level_simple=NA,
                                         PEV1=NA,
                                         PEV2=NA,
                                         PEV3=NA,
                                         PEV4=NA,
                                         PEV5=NA,
                                         PEV6=NA,
                                         PEV7=NA,
                                         PEV8=NA,
                                         PEV9=NA,
                                         PEV10=NA,
                                         nodeType="node"))

# separate the IDs of tips and nodes
index_firstNode = TreeTools::NTip(treesALLMammals_cleaned_pruned_c3) + 1
index_lastNode = TreeTools::NTip(treesALLMammals_cleaned_pruned_c3) + treesALLMammals_cleaned_pruned_c3$Nnode

# add tips id to the dataframe
Dataset_01_withNodes$tree_id <- NA
for (j in 1:(index_firstNode-1)){
  Dataset_01_withNodes$tree_id[which(as.character(Dataset_01_withNodes$MDD_id) == geiger::tips(treesALLMammals_cleaned_pruned_c3, j))] <- j
}

# add nodes id to the dataframe and mean values of trade
Dataset_01_withNodes$tree_id[index_firstNode:index_lastNode] <- index_firstNode:index_lastNode

# all traded
for (n in index_firstNode:index_lastNode){
  Dataset_01_withNodes$pubs_pubmed_all_log[n] <- median(Dataset_01_withNodes$pubs_pubmed_all_log[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$pubs_ALLdisease_PROP_logit[n] <- median(Dataset_01_withNodes$pubs_ALLdisease_PROP_logit[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$pubs_pubmed_all_log_rescaled[n] <- median(Dataset_01_withNodes$pubs_pubmed_all_log_rescaled[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$pubs_ALLdisease_PROP_logit_rescaled[n] <- median(Dataset_01_withNodes$pubs_ALLdisease_PROP_logit_rescaled[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  
  Dataset_01_withNodes$PEV1[n] <- median(Dataset_01_withNodes$PEV1[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$PEV2[n] <- median(Dataset_01_withNodes$PEV2[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$PEV3[n] <- median(Dataset_01_withNodes$PEV3[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$PEV4[n] <- median(Dataset_01_withNodes$PEV4[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$PEV5[n] <- median(Dataset_01_withNodes$PEV5[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$PEV6[n] <- median(Dataset_01_withNodes$PEV6[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$PEV7[n] <- median(Dataset_01_withNodes$PEV7[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$PEV8[n] <- median(Dataset_01_withNodes$PEV8[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$PEV9[n] <- median(Dataset_01_withNodes$PEV9[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
  Dataset_01_withNodes$PEV10[n] <- median(Dataset_01_withNodes$PEV10[which(Dataset_01_withNodes$MDD_id %in% geiger::tips(treesALLMammals_cleaned_pruned_c3, n) & Dataset_01_withNodes$nodeType == "tip")] )
}

# attribute tree ids to rownames
row.names(Dataset_01_withNodes) <- Dataset_01_withNodes$tree_id
Dataset_01_withNodes$node <- Dataset_01_withNodes$tree_id

# Join the data to the tree!
# from: https://mran.microsoft.com/snapshot/2018-07-11/web/packages/tidytree/vignettes/tidytree.html
treesALLMammals_cleaned_pruned_c3_DF <- tibble::as_data_frame(treesALLMammals_cleaned_pruned_c3)
treeData_all <- tidytree::as.treedata(full_join(treesALLMammals_cleaned_pruned_c3_DF, Dataset_01_withNodes, by = 'node'))

treeData_all <- treeData_all %>%
  mutate(pubs_all = scales::rescale(pubs_pubmed_all_log),
         pubs_propDisease = scales::rescale(pubs_ALLdisease_PROP_logit)) %>%
  mutate(pubs_X = pubs_pubmed_all_log_rescaled*pubs_ALLdisease_PROP_logit_rescaled) %>% 
  mutate(pubs_X_logit = scales::rescale(car::logit(pubs_X), to=c(0,1))
  ) 


# set color scheme for tree branches
pal_RE <- paletteer::paletteer_c("palr::bathy_deep_pal", 100, direction = -1)
pal_trunc_RE <- pal_RE[10:90]  

## plotting legend inset for research effort index
fig1a_inset <-
  ggplot(treeData_all, aes(pubs_pubmed_all_log, pubs_ALLdisease_PROP_logit)) +
  geom_point(aes(color = pubs_X_logit), 
             size=2.5, alpha=1, 
             position=position_jitter(width=0.025, height=0.025)
  ) +
  scale_x_continuous(breaks=log(c(0,10,100,1000,10000)+1),
                     labels=c(0,10,100,1000,10000)) +
  scale_y_continuous(breaks=car::logit(c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)),
                     labels=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
                     limits = car::logit(c(0.024, 0.925))
  ) +
  scale_color_gradientn(
    colours = pal_trunc_RE,
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1)
  ) +
  xlab("Total number of publications") +
  ylab("Proportion of publications on pathogens") +
  theme_bw() +
  theme(legend.position= "none")
fig1a_inset


system.time(ggsave("./outputs/fig1a_inset.pdf", fig1a_inset, width = 12, height = 12, units = "cm", limitsize = FALSE))



### Fig1a - main

## --- 0) Safety: kill any stray 'mapping' object that can confuse ggplot
if (exists("mapping", inherits = FALSE)) rm(mapping)

## --- 1) Base tree (continuous branch colour)
a1 <- ggtree(treeData_all, layout = "circular",
             branch.length = "branch.length", size = 0.2) +
  geom_tree(aes(color = pubs_X_logit), size = 0.2) +
  scale_color_gradientn(colours = pal_trunc_RE,
                        breaks = c(0, 0.5, 1),
                        labels = c(0, 0.5, 1))
a1$mapping <- ggplot2::aes()   # belt & suspenders

## --- 2) Prepare annotations (rows = tips, columns = 1; plain matrix)
Dataset_01 <- as.data.frame(Dataset_01)

# use character IDs to match tree tip labels
Dataset_01$label <- as.character(Dataset_01$MDD_id_chr)
rownames(Dataset_01) <- Dataset_01$label

# tidy up values so factors don’t become NA unexpectedly
norm_yn <- function(x) {
  v <- tolower(trimws(as.character(x)))
  v[v %in% c("1","true","y","yes")] <- "yes"
  v[v %in% c("0","false","n","no")]  <- "no"
  v[v == "."] <- NA_character_
  v
}
Dataset_01$Traded_all               <- norm_yn(Dataset_01$Traded_all)
Dataset_01$Synanthropy_level_simple <- norm_yn(Dataset_01$Synanthropy_level_simple)
Dataset_01$WildMeat                 <- norm_yn(Dataset_01$WildMeat)
# keep original labels for zoonosis column
Dataset_01$Zoonose_vector_1 <- dplyr::case_when(
  is.na(Dataset_01$Zoonose_vector_1) ~ NA_character_,
  TRUE ~ as.character(Dataset_01$Zoonose_vector_1)
)

tips <- as.character(a1$data$label)

mk_ann <- function(df, colname, levels = NULL) {
  ann <- df[, colname, drop = FALSE]
  # align rows to tip order; introduces NA where metadata is missing
  ann <- ann[match(tips, rownames(ann)), , drop = FALSE]
  rownames(ann) <- tips
  if (!is.null(levels)) ann[[1]] <- factor(ann[[1]], levels = levels)
  as.matrix(ann)
}

ann_traded <- mk_ann(Dataset_01, "Traded_all",               levels = c("no","yes"))
ann_syn    <- mk_ann(Dataset_01, "Synanthropy_level_simple", levels = c("no","yes"))
ann_wm     <- mk_ann(Dataset_01, "WildMeat",                 levels = c("no","yes"))
ann_zoon   <- mk_ann(Dataset_01, "Zoonose_vector_1",         levels = c("NotZoonotic","Zoonotic"))

## --- 3) Add heatmaps cumulatively + reset fill between them
p <- a1

p <- ggtree::gheatmap(p, ann_traded, offset = 0.5,  width = 0.05,
                      colnames = FALSE, color = NA) +
  ggplot2::scale_fill_manual(values = c(no="white", yes="black"),
                             na.value = "white", name = "Traded")
p$mapping <- ggplot2::aes()
p <- p + ggnewscale::new_scale_fill()

p <- ggtree::gheatmap(p, ann_syn, offset = 12.5, width = 0.05,
                      colnames = FALSE, color = NA) +
  ggplot2::scale_fill_manual(values = c(no="white", yes="#5d8891ff"),
                             na.value = "white", name = "Synanthropic")
p$mapping <- ggplot2::aes()
p <- p + ggnewscale::new_scale_fill()

p <- ggtree::gheatmap(p, ann_wm, offset = 24.5, width = 0.05,
                      colnames = FALSE, color = NA) +
  ggplot2::scale_fill_manual(values = c(no="white", yes="#e3d585"),
                             na.value = "white", name = "Wild meat")
p$mapping <- ggplot2::aes()
p <- p + ggnewscale::new_scale_fill()

p <- ggtree::gheatmap(p, ann_zoon, offset = 36.5, width = 0.05,
                      colnames = FALSE, color = NA) +
  ggplot2::scale_fill_manual(values = c(NotZoonotic="white", Zoonotic="#780000"),
                             na.value = "white", name = "Zoonotic host")
p$mapping <- ggplot2::aes()

## --- 4) Clade labels (use character tip IDs + constant colour)
o <- 0
color_orderrings <- c("grey65","grey85","grey65","grey65","grey65",
                      "grey85","grey85","grey85","grey85","grey65")
for (ORDER in mainOrders_all) {
  tips1 <- subset(Dataset_01, order2 == ORDER)[["label"]]   # character IDs
  if (!length(tips1)) next
  node1 <- castor::get_mrca_of_set(treesALLMammals_cleaned_pruned_c3, tips1)
  if (is.null(node1) || !length(node1) || is.na(node1)) next
  o <- o + 1
  p <- p + ggtree::geom_cladelabel(
    node = node1, label = ORDER,
    offset = 54, offset.text = 1, align = TRUE,
    barsize = 1.5, fontsize = 1, color = color_orderrings[o]
  )
}

## --- 5) Final tweaks
fig1a <- p %>%
  ggtree::open_tree(20) %>%
  ggtree::rotate_tree(-80) +
  ggplot2::xlim(-20, NA) +
  ggtree::geom_treescale(x = 0, y = 100, width = 20, offset = -45,
                         linesize = 1, fontsize = 3, color = "black") +
  ggplot2::theme(legend.position = c(0.8, 0.2),
                 legend.key.size = grid::unit(5, "mm"))

fig1a

ggsave("./outputs/fig1a.pdf", fig1a, width = 18, height = 18, units = "cm", limitsize = FALSE)
########################


###>> Figure S1

# CITES LEMIS DSW
euler_Trade_S1a <- list(
  CITES = Dataset_01 %>% filter(Traded_all_CITES =="yes") %>% pull(MDD_id),
  LEMIS = Dataset_01 %>% filter(Traded_all_LEMIS2025 =="yes") %>% pull(MDD_id),
  DSW = Dataset_01 %>% filter(Traded_illegal =="yes") %>% pull(MDD_id),
  All = Dataset_01  %>% pull(MDD_id)
)

figS1a <- plot(euler(euler_Trade_S1a, shape = "ellipse"),
               fill = alpha(c(viridis::viridis(3, option = "D"), "white"), alpha= 0.75), # , "grey99", "grey99"
               lwd=1.5,
               edges = "black", #, "black", "grey99"
               labels = list(cex = 0.5), 
               quantities = list(cex = 0.75))


# live-product
euler_Trade_S1b <- list(
  Live = Dataset_01 %>% filter(Traded_all_live =="yes") %>% pull(MDD_id),
  Product = Dataset_01 %>% filter(Traded_all_product =="yes") %>% pull(MDD_id),
  All = Dataset_01 %>% pull(MDD_id)
)

figS1b <- plot(euler(euler_Trade_S1b, shape = "ellipse"),
               fill = c("#ffb057ff", "#7e4618ff", "white"),
               lwd=1.5,
               edges = c('black'),
               labels = list(cex = 0.5), 
               quantities = list(cex = 0.75))

# legal-illegal
euler_Trade_S1c <- list(
  Legal = Dataset_01 %>% filter(Traded_legal =="yes") %>% pull(MDD_id),
  Illegal = Dataset_01 %>% filter(Traded_illegal =="yes") %>% pull(MDD_id),
  All = Dataset_01 %>% pull(MDD_id)
)

figS1c <- plot(euler(euler_Trade_S1c, shape = "ellipse"),
               fill = c("#adc178ff", "#38a3a5ff", "white"),
               lwd=1.5,
               edges = c('black'),
               labels = list(cex = 0.5), 
               quantities = list(cex = 0.75))

figS1 <- ggpubr::ggarrange(figS1a, figS1b, figS1c, ncol=3)
figS1

system.time(ggsave("./outputs/figS1.png", figS1, width = 12, height = 12, units = "cm", limitsize = FALSE))
########################



###>> Figure S2

# pick the variables you want to plot
pev_vars <- paste0("PEV", 1:10)

# build all plots
tree_plots <- purrr::map(setNames(pev_vars, pev_vars),
                         ~ make_tree(.x, tree = treeData_all))

tree_plots_2 <- purrr::map(tree_plots, ~ .x %>% 
                             open_tree(20) %>% 
                             rotate_tree(-80) +
                             xlim(-20, NA))

# arrange in a grid (2 rows x 5 cols) with one shared legend on the right
fig_S2 <- ggarrange(
  plotlist = tree_plots_2,
  nrow = 2, ncol = 5,
  common.legend = TRUE, legend = "right",
  align = "hv"
)

fig_S2
ggsave("./outputs/fig_S2.png", fig_S2, width = 16, height = 7, dpi = 300)
########################



###>> Figure 1b 

# All mammals
euler_Trade_0 <- list(
  Zoonotic = Dataset_01 %>% filter(Zoonose_vector_1=="Zoonotic") %>% pull(MDD_id),
  Synanthropic = Dataset_01 %>% filter(Synanthropy_level_simple =="yes") %>% pull(MDD_id),
  WildMeat = Dataset_01 %>% filter(WildMeat =="yes") %>% pull(MDD_id),
  All = Dataset_01  %>% pull(MDD_id)
)

# Traded mammals
euler_Trade_1 <- list(
  Zoonotic = Dataset_01 %>% filter(Traded_all =="yes" & Zoonose_vector_1=="Zoonotic") %>% pull(MDD_id),
  Synanthropic = Dataset_01 %>% filter(Traded_all =="yes" & Synanthropy_level_simple =="yes") %>% pull(MDD_id),
  WildMeat = Dataset_01 %>% filter(Traded_all =="yes" & WildMeat =="yes") %>% pull(MDD_id),
  All_traded = Dataset_01 %>% filter(Traded_all =="yes")  %>% pull(MDD_id),
  All = Dataset_01  %>% pull(MDD_id)
)

# Not traded mammals
euler_Trade_2 <- list(
  Zoonotic = Dataset_01 %>% filter(Traded_all =="no" & Zoonose_vector_1=="Zoonotic") %>% pull(MDD_id),
  Synanthropic = Dataset_01 %>% filter(Traded_all =="no" & Synanthropy_level_simple =="yes") %>% pull(MDD_id),
  WildMeat = Dataset_01 %>% filter(Traded_all =="no" & WildMeat =="yes") %>% pull(MDD_id),
  All_nottraded = Dataset_01 %>% filter(Traded_all =="no")  %>% pull(MDD_id),
  All = Dataset_01  %>% pull(MDD_id)
)


fig1b <- ggpubr::ggarrange(
  plot(euler(euler_Trade_0, shape = "ellipse"),
       fill = alpha(c('#780000', '#2cacdf9a', "#e3d585ff"), alpha= 0.75), # , "grey99", "grey99"
       lwd=2.5,
       edges = c('grey25', 'grey25', 'grey25'), #, "black", "grey99"
       labels = list(cex = 0.5), 
       quantities = list(cex = 0.75)),
  plot(euler(euler_Trade_1, shape = "ellipse"),
       fill = alpha(c('#780000', '#2cacdf9a', "#e3d585ff"), alpha= 0.75), # , "grey99", "grey99"
       lwd=2.5,
       edges = c('grey25', 'grey25', 'grey25'), #, "black", "grey99"
       labels = list(cex = 0.5), 
       quantities = list(cex = 0.75)),
  plot(euler(euler_Trade_2, shape = "ellipse"),
       fill = alpha(c('#780000', '#2cacdf9a', "#e3d585ff"), alpha= 0.75), # , "grey99", "grey99"
       lwd=2.5,
       edges = c('grey25', 'grey25', 'grey25'), #, "black", "grey99"
       labels = list(cex = 0.5), 
       quantities = list(cex = 0.75)),
  ncol=3
)

fig1b

system.time(ggsave("./outputs/fig1b.pdf", fig1b, width = 12, height = 6, units = "cm", limitsize = FALSE))
########################


###>> Main model

Dataset_01_m3 <- Dataset_01 
Dataset_01_m3 %>% head

M1_01 <- glmmTMB(Zoonose_vector_1  ~ 1 
                 + Traded_all
                 
                 + Synanthropy_level_simple
                 + WildMeat
                 
                 + pubs_X_scale
                 
                 + PEV1
                 + PEV2
                 + PEV3
                 + PEV4
                 + PEV5
                 + PEV6
                 + PEV7
                 + PEV8
                 + PEV9
                 + PEV10
                 
                 + (1|BioRealm)
                 ,
                 
                 family=binomial(link="logit"),
                 
                 data= Dataset_01_m3,
                 control = glmmTMBControl(parallel = 8)
)


summary(M1_01)
plot(simulateResiduals(M1_01))
performance::r2(M1_01)

performance::check_singularity(M1_01)
performance::check_collinearity(M1_01)

# effect sizes
plot_model(M1_01, type="std2") + theme_bw()
sdteff.M1_01 <- get_model_data(M1_01, type="std2")
sdteff.M1_01



### Compute models effects
## Trade
M1_01_effect_Traded_all <- ggaverage(M1_01, 
                                     c("Traded_all"),
                                     type="response")

M1_01_effect_Traded_all <- M1_01_effect_Traded_all %>% as.data.frame
M1_01_effect_Traded_all$x <- factor(M1_01_effect_Traded_all$x, levels=c("no", "yes"))
colnames(M1_01_effect_Traded_all)[1] <- "Traded"

# marginal risk ratio
M1_01_effect_Traded_all$predicted[2]/M1_01_effect_Traded_all$predicted[1] # 1.49

## Synanthropy
M1_01_effect_Synanthropy_level_simple <- ggaverage(M1_01, 
                                                   c("Synanthropy_level_simple"),
                                                   type="response")

M1_01_effect_Synanthropy_level_simple <- M1_01_effect_Synanthropy_level_simple %>% as.data.frame
M1_01_effect_Synanthropy_level_simple$x <- factor(M1_01_effect_Synanthropy_level_simple$x, levels=c("no", "yes"))
colnames(M1_01_effect_Synanthropy_level_simple)[1] <- "Traded"

# marginal risk ratio
M1_01_effect_Synanthropy_level_simple$predicted[2]/M1_01_effect_Synanthropy_level_simple$predicted[1] # 1.2

## Wildmeat
M1_01_effect_WildMeat <- ggaverage(M1_01, 
                                   c("WildMeat"),
                                   type="response")

M1_01_effect_WildMeat <- M1_01_effect_WildMeat %>% as.data.frame
M1_01_effect_WildMeat$x <- factor(M1_01_effect_WildMeat$x, levels=c("no", "yes"))
colnames(M1_01_effect_WildMeat)[1] <- "Traded"

# marginal risk ratio
M1_01_effect_WildMeat$predicted[2]/M1_01_effect_WildMeat$predicted[1] # 1.1



## Figure 1 panels d, e and f
fig1def <- ggpubr::ggarrange(
  ggplot(data=M1_01_effect_Traded_all,
         aes(x=Traded, y=predicted, 
             fill= Traded,
             col= Traded,
             ymin=conf.low, ymax=conf.high, 
             group=group)) +
    geom_linerange(position = position_dodge(.25), 
                   size=2.5,
                   col="black"
    ) +
    scale_fill_manual(values=c("no" = "grey", 
                               "yes" = "black"),
                      na.value = "grey",
                      name="Traded") +
    scale_colour_manual(values=c("no" = "grey", 
                                 "yes" = "black"),
                        na.value = "grey",
                        name="Traded") +
    geom_point(position=position_dodge(width=0.25),
               shape=22, 
               size=7,
               stroke=2,
               col="black"
    ) +
    scale_y_continuous(limits=c(0.1, 0.25)) +
    xlab("Traded") +
    ylab('Probability of sharing at least one pathogen with humans') +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank()),
  
  ggplot(data=M1_01_effect_Synanthropy_level_simple,
         aes(x=Traded, y=predicted, 
             fill= Traded,
             col= Traded,
             ymin=conf.low, ymax=conf.high, 
             group=group)) +
    geom_linerange(position = position_dodge(.25), 
                   size=2.5,
                   col="black"
    ) +
    geom_point(position=position_dodge(width=0.25),
               shape=22, 
               size=7,
               stroke=2,
               col="black"
    ) +
    scale_fill_manual(values=c("no" = "grey", 
                               "yes" = "#2cacdfff"),
                      na.value = "grey",
                      name="Synanthropic") +
    scale_colour_manual(values=c("no" = "grey", 
                                 "yes" = "#2cacdfff"),
                        na.value = "grey",
                        name="Synanthropic") +
    scale_y_continuous(limits=c(0.1, 0.25)) +
    xlab("Synanthropic") +
    ylab('Probability of sharing at least one pathogen with humans') +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank()),
  
  ggplot(data=M1_01_effect_WildMeat,
         aes(x=Traded, y=predicted, 
             fill= Traded,
             col= Traded,
             ymin=conf.low, ymax=conf.high, 
             group=group)) +
    geom_linerange(position = position_dodge(.25), 
                   size=2.5,
                   col="black"
    ) +
    geom_point(position=position_dodge(width=0.25),
               shape=22, 
               size=7,
               stroke=2,
               col="black"
    ) +
    scale_fill_manual(values=c("no" = "grey", 
                               "yes" = "#e3d585"),
                      na.value = "grey",
                      name="Wild meat") +
    scale_colour_manual(values=c("no" = "grey", 
                                 "yes" = "#e3d585"),
                        na.value = "grey",
                        name="Wild meat") +
    scale_y_continuous(limits=c(0.1, 0.25)) +
    xlab("Wild meat") +
    ylab('Probability of sharing at least one pathogen with humans') +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank()),
  ncol=3
)

fig1def
system.time(ggsave("./outputs/fig1def.pdf", fig1def, width = 12, height = 8, units = "cm", limitsize = FALSE))
########################



### Testing for phylogenetic signal in residuals
sim_res_M1_01 <- simulateResiduals(fittedModel = M1_01)
residuals_M1_01 <- residuals(sim_res_M1_01)
names(residuals_M1_01) <- Dataset_01_m3$MDD_id

# Match residuals to tree
residuals_M1_01.2 <- residuals_M1_01[names(residuals_M1_01) %in% tree_aug$tip.label]
consTree_M01_01 <- drop.tip(tree_aug, setdiff(tree_aug$tip.label, names(residuals_M1_01.2))) 

# Test for phylogenetic signal
phySignal_M1_01.2 <- phylosig(consTree_M01_01, residuals_M1_01.2, method = "K", test=TRUE) # can take a several minutes
phySignal_M1_01.2
########################


####>> Check consistence with phylogenetic model
tips_to_remove_2 <- c(tree_aug$tip.label[which(tree_aug$tip.label %notin% as.character(Dataset_01_m3$MDD_id))] )
treesALLMammals_cleaned_pruned_c4 <- ape::drop.tip(phy=tree_aug, tip=tips_to_remove_2) 

Dataset_01_m3.phylolm <- Dataset_01_m3 %>% 
  filter(MDD_id %in% treesALLMammals_cleaned_pruned_c4$tip.label) %>% 
  mutate(Zoonose_vector_1_bin = ifelse(Zoonose_vector_1=="Zoonotic", 1, 0)) %>% 
  as.data.frame()
rownames(Dataset_01_m3.phylolm) <- Dataset_01_m3.phylolm$MDD_id

# Phylogenetic model
M1_01_phylolm = phylolm::phyloglm(Zoonose_vector_1_bin  ~ 1  
                                  + Traded_all
                                  
                                  + Synanthropy_level_simple
                                  + WildMeat
                                  
                                  + pubs_X_scale
                                  
                                  ,
                                  method="logistic_MPLE",
                                  phy=treesALLMammals_cleaned_pruned_c4,
                                  data=Dataset_01_m3.phylolm,
                                  boot=100,
                                  btol= 50)

summary(M1_01_phylolm)
tab_model(M1_01_phylolm)

# Getting model's R2
mod_full <- M1_01_phylolm

# Null model (intercept only)
mod_null <- phylolm::phyloglm(Zoonose_vector_1_bin ~ 1,
                              method = "logistic_MPLE",
                              phy = treesALLMammals_cleaned_pruned_c4,
                              data = Dataset_01_m3.phylolm,
                              boot = 100,
                              btol = 50)

rr2::R2_lik(mod_full, mod_null)
########################




###>> Structural Equation Model (fig. 1d)

# adapt variables to ensure good functioning of SEM
Dataset_01_m3_SEM <- Dataset_01_m3 %>% mutate(var_T = ifelse(Traded_all=="yes", 1, 0),
                                                  var_S = ifelse(Synanthropy_level_simple=="yes", 1, 0),
                                                  var_Z = ifelse(Zoonose_vector_1=="Zoonotic", 1, 0),
                                                  var_W = ifelse(WildMeat=="yes", 1, 0))

# model to predict research effort
M1_01_RE_mod <-
  glmmTMB(pubs_X_scale  ~ 1
          + var_S
          + var_T
          + var_W
          
          + (1|BioRealm)
          
          + PEV1
          + PEV2
          + PEV3
          + PEV4
          + PEV5
          + PEV6
          + PEV7
          + PEV8
          + PEV9
          + PEV10
          ,
          family=gaussian(link="identity"),
          data= Dataset_01_m3_SEM,
          control = glmmTMBControl(parallel = 8))
summary(M1_01_RE_mod)

# model to predict the probability of trade
M1_01_T_mod <- 
  glmmTMB(var_T  ~ 1
          + var_S
          + var_W
          
          + (1|BioRealm)
          
          + PEV1
          + PEV2
          + PEV3
          + PEV4
          + PEV5
          + PEV6
          + PEV7
          + PEV8
          + PEV9
          + PEV10
          ,
          family=binomial(link="logit"),
          data= Dataset_01_m3_SEM,
          control = glmmTMBControl(parallel = 8))       
summary(M1_01_T_mod)

# model to predict the probability of sharing zoonotic pathogens with humans
M1_01_Z_mod <- 
  glmmTMB(var_Z  ~ 1
          + var_T
          + var_S
          + var_W
          
          + pubs_X_scale
          
          + (1|BioRealm)
          
          + PEV1
          + PEV2
          + PEV3
          + PEV4
          + PEV5
          + PEV6
          + PEV7
          + PEV8
          + PEV9
          + PEV10
          ,
          family=binomial(link="logit"),
          data= Dataset_01_m3_SEM,
          control = glmmTMBControl(parallel = 8))
summary(M1_01_Z_mod)

M1_01_model_full <- psem(M1_01_RE_mod,
                         M1_01_T_mod,
                         M1_01_Z_mod)
M1_01_model_full
plot(M1_01_model_full)

M1_01_sumSEM <- summary(M1_01_model_full, conserve="TRUE")
M1_01_sumSEM


## Getting direct and indirect effect strenghts
# Get standardized coefficients for all paths
std_coefs <- piecewiseSEM::coefs(M1_01_model_full, standardize = "scale") %>%
  as.data.frame()
colnames(std_coefs) <- c("Response", "Predictor", "Estimate", "Std.Error",
                         "DF", "Crit.Value", "P.Value", "Std.Estimate", "stars")

# Keep only key paths
std_coefs %>%
  dplyr::select(Response, Predictor, Std.Estimate, P.Value) %>%
  filter(Response %in% c("var_T", "pubs_X_scale", "var_Z"),
         Predictor %in% c("var_S", "var_W", "var_T", "pubs_X_scale"))


# Extract relevant standardized effects
b_S_T <- std_coefs %>% filter(Response=="var_T", Predictor=="var_S") %>% pull(Std.Estimate)
b_W_T <- std_coefs %>% filter(Response=="var_T", Predictor=="var_W") %>% pull(Std.Estimate)

b_T_R <- std_coefs %>% filter(Response=="pubs_X_scale", Predictor=="var_T") %>% pull(Std.Estimate)
b_S_R <- std_coefs %>% filter(Response=="pubs_X_scale", Predictor=="var_S") %>% pull(Std.Estimate)
b_W_R <- std_coefs %>% filter(Response=="pubs_X_scale", Predictor=="var_W") %>% pull(Std.Estimate)

b_R_Z <- std_coefs %>% filter(Response=="var_Z", Predictor=="pubs_X_scale") %>% pull(Std.Estimate)
b_T_Z <- std_coefs %>% filter(Response=="var_Z", Predictor=="var_T") %>% pull(Std.Estimate)
b_S_Z <- std_coefs %>% filter(Response=="var_Z", Predictor=="var_S") %>% pull(Std.Estimate)
b_W_Z <- std_coefs %>% filter(Response=="var_Z", Predictor=="var_W") %>% pull(Std.Estimate)

# compute indirect and total effects 
effects <- tibble(
  path = c("Trade → Zoonotic", 
           "Synanthropy → Zoonotic", 
           "Wild meat → Zoonotic"),
  
  # direct effects (from the Z model)
  direct = c(b_T_Z, b_S_Z, b_W_Z),
  
  # indirect effects (sum of all possible mediations)
  indirect = c(
    b_T_R * b_R_Z,                                      # Trade → Research → Zoonotic
    b_S_T*b_T_Z + b_S_T*b_T_R*b_R_Z + b_S_R*b_R_Z,      # Synanthropy → Trade(+Research) → Zoonotic
    b_W_T*b_T_Z + b_W_T*b_T_R*b_R_Z + b_W_R*b_R_Z       # WildMeat → Trade(+Research) → Zoonotic
  )
) %>%
  mutate(total = direct + indirect) %>%
  arrange(desc(total))

effects
#####################



#################################################################
########                SENSITIVITY ANALYSES              #######

### Sensitivity on wildlife trade dataset
# legal trade (CITES)
M1_01_CITES <- glmmTMB(Zoonose_vector_1  ~ 1 
                       + Traded_all_CITES
                       + Synanthropy_level_simple
                       + WildMeat
                       
                       + pubs_X_scale
                       
                       + PEV1
                       + PEV2
                       + PEV3
                       + PEV4
                       + PEV5
                       + PEV6
                       + PEV7
                       + PEV8
                       + PEV9
                       + PEV10
                       
                       + (1|BioRealm)
                       
                       ,
                       
                       family=binomial(link="logit"),
                       
                       data= Dataset_01_m3,
                       control = glmmTMBControl(parallel = 8)
)

summary(M1_01_CITES)
tab_model(M1_01_CITES)


##>> legal trade (LEMIS)
M1_01_LEMIS <- glmmTMB(Zoonose_vector_1 ~ 1
                       + Traded_all_LEMIS2025
                       + Synanthropy_level_simple
                       + WildMeat
                       
                       + pubs_X_scale
                       
                       + PEV1
                       + PEV2
                       + PEV3
                       + PEV4
                       + PEV5
                       + PEV6
                       + PEV7
                       + PEV8
                       + PEV9
                       + PEV10
                       
                       + (1|BioRealm)
                       , 
                       
                       family=binomial(link="logit"),
                       
                       data= Dataset_01_m3,
                       control = glmmTMBControl(parallel = 8)
)

summary(M1_01_LEMIS)
tab_model(M1_01_LEMIS)


## Effect sizes
# Full model
effectsize_M1_01_BOTH_raw_0.95 <- get_model_data(M1_01, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_01_BOTH_raw_0.99 <- get_model_data(M1_01, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_01_BOTH_raw <- effectsize_M1_01_BOTH_raw_0.95 %>% 
  left_join(effectsize_M1_01_BOTH_raw_0.99, by="term") %>% 
  mutate(dataset= "BOTH")

# CITES only model
effectsize_M1_01_CITES_raw_0.95 <- get_model_data(M1_01_CITES, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_01_CITES_raw_0.99 <- get_model_data(M1_01_CITES, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_01_CITES_raw <- effectsize_M1_01_CITES_raw_0.95 %>% 
  left_join(effectsize_M1_01_CITES_raw_0.99, by="term") %>% 
  mutate(dataset= "CITES")

effectsize_M1_01_CITES_raw$term = as.character(effectsize_M1_01_CITES_raw$term)
effectsize_M1_01_CITES_raw$term[1] = "Traded_allyes"
effectsize_M1_01_CITES_raw$term = as.factor(effectsize_M1_01_CITES_raw$term)


# LEMIS only
effectsize_M1_01_LEMIS_raw_0.95 <- get_model_data(M1_01_LEMIS, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_01_LEMIS_raw_0.99 <- get_model_data(M1_01_LEMIS, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_01_LEMIS_raw <- effectsize_M1_01_LEMIS_raw_0.95 %>% 
  left_join(effectsize_M1_01_LEMIS_raw_0.99, by="term") %>% 
  mutate(dataset= "LEMIS")

effectsize_M1_01_LEMIS_raw$term = as.character(effectsize_M1_01_LEMIS_raw$term)
effectsize_M1_01_LEMIS_raw$term[1] = "Traded_allyes"
effectsize_M1_01_LEMIS_raw$term = as.factor(effectsize_M1_01_LEMIS_raw$term)

# combine
effectsize_M1_01_wildlifetradeDataset_raw <- rbind(effectsize_M1_01_BOTH_raw, 
                                                   effectsize_M1_01_CITES_raw,
                                                   effectsize_M1_01_LEMIS_raw)

effectsize_M1_01_wildlifetradeDataset_cond <- effectsize_M1_01_wildlifetradeDataset_raw %>% mutate(direction = ifelse(estimate>=1, "pos", "neg"))

# no manipulation for the conditional model
effectsize_M1_01_wildlifetradeDataset_cond <- effectsize_M1_01_wildlifetradeDataset_cond %>% 
  mutate(estimate_rescaled = estimate,
         conf.low.95_rescaled = conf.low.95,
         conf.high.95_rescaled = conf.high.95,
         conf.low.99_rescaled = conf.low.99,
         conf.high.99_rescaled = conf.high.99)

# then creating plotting values for visualisation
effectsize_M1_01_wildlifetradeDataset_cond <- effectsize_M1_01_wildlifetradeDataset_cond %>% 
  mutate(estimate_rescaled_2 = ifelse(direction=="pos", estimate_rescaled, (-1/estimate_rescaled)+2),
         conf.low.95_rescaled_2 = ifelse(direction=="pos", conf.low.95_rescaled, (-1/conf.low.95_rescaled)+2),
         conf.high.95_rescaled_2 = ifelse(direction=="pos", conf.high.95_rescaled, (-1/conf.high.95_rescaled)+2),
         conf.low.99_rescaled_2 = ifelse(direction=="pos", conf.low.99_rescaled, (-1/conf.low.99_rescaled)+2),
         conf.high.99_rescaled_2 = ifelse(direction=="pos", conf.high.99_rescaled, (-1/conf.high.99_rescaled)+2))


effectsize_M1_01_wildlifetradeDataset <- effectsize_M1_01_wildlifetradeDataset_cond
effectsize_M1_01_wildlifetradeDataset %>% colnames

effectsize_M1_01_wildlifetradeDataset %>% filter(term=="Traded_allyes") %>% dplyr::select(dataset, estimate, p.value, p.stars) 
effectsize_M1_01_wildlifetradeDataset$term %>% table

# Desired terms
term_levels_all <- c(
  "pubs_X_scale", 
  "WildMeatyes", 
  "Synanthropy_level_simpleyes",
  "Traded_allyes")
dataset_levels_all_01 <- c("LEMIS", "CITES", "BOTH")

# Create full grid
ghost_grid_all<- tidyr::expand_grid(term = term_levels_all,
                                    dataset = dataset_levels_all_01)

# Merge with actual data
effectsize_complete_all <- ghost_grid_all %>%
  left_join(effectsize_M1_01_wildlifetradeDataset, by = c("term", "dataset")) %>%
  mutate(term = factor(term, levels = term_levels_all),
         dataset= factor(dataset, levels = dataset_levels_all_01),
         signif = ifelse(p.value<0.05, "yes", "no"),
         alpha_value = ifelse(p.value > 0.05 | is.na(p.value), 0.1, 
                              ifelse(p.value > 0.01, 0.6, 1)))  # ensure correct order

ci_min_all <- -12
dodge_width_all <- 0.6

# Conditional
figS3a <- ggplot(data = effectsize_complete_all,
                 aes(x = estimate, y = term, fill = term, shape= dataset),
                 alpha=alpha_value) +
  geom_vline(xintercept = 1, linetype = 1, color = "black", size = 1) +
  geom_errorbarh(aes(xmin = conf.low.95_rescaled_2, xmax = conf.high.95_rescaled_2, color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width_all),
                 size = 1.25, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low.99_rescaled_2, xmax = conf.high.99_rescaled_2, color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width_all),
                 size = 0.5, na.rm = TRUE) +
  geom_point(aes(alpha = alpha_value),
             stroke=1,
             #shape = 21, 
             size = 2.5, 
             position = position_dodge(width = dodge_width_all),
             na.rm = TRUE) +
  scale_x_continuous(limits = c(0.5, 7.25),
                     breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7),
                     labels = c(3, 2, 1, 2, 3, 4, 5, 6, 7)) +
  scale_shape_manual(values=c(21, 22, 23),
                     breaks=c("BOTH", "CITES", "LEMIS"),
                     labels= c("CITES + LEMIS", "CITES only", "LEMIS only"),
                     name= "Wildlife trade \ndataset used") +
  # annotate("text", x=-Inf, y=Inf, label= "d", size=8, vjust = 1.5, hjust = -1) +
  theme_bw() +
  scale_y_discrete(labels = c("Traded_allyes" = "Traded",
                              "Synanthropy_level_simpleyes" = "Synanthropic",
                              "WildMeatyes" = "Used as\nwild meat",
                              "pubs_X_scale" = "Research\neffort")) +
  scale_colour_manual(values = custom_colors,
                      labels = c("Traded_allyes" = "Traded",
                                 "Synanthropy_level_simpleyes" = "Synanthropic",
                                 "WildMeatyes" = "Used as\nwild meat",
                                 "pubs_X_scale" = "Research\neffort")) +
  scale_fill_manual(values = custom_colors,
                    labels = c("Traded_allyes" = "Traded",
                               "Synanthropy_level_simpleyes" = "Synanthropic",
                               "WildMeatyes" = "Used as\nwild meat",
                               "pubs_X_scale" = "Research\neffort")) +
  scale_alpha_continuous(breaks= c(0.1, 0.5, 1),
                         labels= c("0.1" = "p > 0.05",
                                   "0.5" = "0.05 > p > 0.01",
                                   "1" = "p < 0.01"),
                         name= "Significance level") +
  labs(#title= "Response: Probability of sharing\nat least one pathogen with humans",
    x= "Std effect size (odds ratio)") +
  guides(colour = "none",
         fill = "none")  +
  theme(legend.position = c(0.8, 0.75),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(size=10), 
        title = element_text(size=10), 
        # axis.text.y = element_blank(), 
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
  )

figS3a


### Sensitivity on host-pathogen interaction dataset 

# Model M01 using CLOVER virus only
M1_01_clovervirus <- glmmTMB(Zoonose_virus_vector_1  ~ 1 
                             + Traded_all
                             
                             + Synanthropy_level_simple
                             + WildMeat
                             
                             + pubs_X_scale
                             
                             + PEV1
                             + PEV2
                             + PEV3
                             + PEV4
                             + PEV5
                             + PEV6
                             + PEV7
                             + PEV8
                             + PEV9
                             + PEV10
                             
                             + (1|BioRealm)
                             
                             ,
                             
                             family=binomial(link="logit"),
                             
                             data= Dataset_01_m3,
                             control = glmmTMBControl(parallel = 8)
)


summary(M1_01_clovervirus)
plot(simulateResiduals(fittedModel = M1_01_clovervirus))
tab_model(M1_01_clovervirus)


# Model M01 using CLOVER no EID2
Dataset_01_m3$noEID2_ZOONOTIC_NB_all %>% table
Dataset_01_m3 <- Dataset_01_m3 %>% 
  mutate(Zoonose_vector_1_noEID2= as.factor(ifelse(noEID2_ZOONOTIC_NB_all>0, "Zoonotic", "NotZoonotic")))
Dataset_01_m3$Zoonose_vector_1_noEID2 %>% table

M1_01_clovernoEID2 <- glmmTMB(Zoonose_vector_1_noEID2  ~ 1 
                              + Traded_all
                              
                              + Synanthropy_level_simple
                              + WildMeat
                              
                              + pubs_X_scale
                              
                              + PEV1
                              + PEV2
                              + PEV3
                              + PEV4
                              + PEV5
                              + PEV6
                              + PEV7
                              + PEV8
                              + PEV9
                              + PEV10
                              
                              + (1|BioRealm)
                              
                              ,
                              
                              family=binomial(link="logit"),
                              
                              data= Dataset_01_m3,
                              control = glmmTMBControl(parallel = 8)
)


summary(M1_01_clovernoEID2)
plot(simulateResiduals(fittedModel = M1_01_clovernoEID2))
tab_model(M1_01_clovernoEID2)



# Model M01 using VIRION
M1_01_virion <- glmmTMB(VIRION_Zoonose_vector_1  ~ 1 
                        + Traded_all
                        
                        + Synanthropy_level_simple
                        + WildMeat
                        
                        + pubs_X_scale
                        
                        + PEV1
                        + PEV2
                        + PEV3
                        + PEV4
                        + PEV5
                        + PEV6
                        + PEV7
                        + PEV8
                        + PEV9
                        + PEV10
                        
                        + (1|BioRealm)
                        
                        ,
                        
                        family=binomial(link="logit"),
                        
                        data= Dataset_01_m3,
                        control = glmmTMBControl(parallel = 8)
)

summary(M1_01_virion)
plot(simulateResiduals(fittedModel = M1_01_virion))
tab_model(M1_01_virion)


# Model M01 using VIRION no PREDICT
Dataset_01_m3$noPREDICT_VIRION_ZOONOTIC_NB_viruses %>% table
Dataset_01_m3 <- Dataset_01_m3 %>% 
  mutate(VIRION_Zoonose_vector_1_noPREDICT= as.factor(ifelse(noPREDICT_VIRION_ZOONOTIC_NB_viruses>0, "Zoonotic", "NotZoonotic")))
Dataset_01_m3$VIRION_Zoonose_vector_1_noPREDICT %>% table

M1_01_virionnoPREDICT <- glmmTMB(VIRION_Zoonose_vector_1_noPREDICT  ~ 1 
                                 + Traded_all
                                 
                                 + Synanthropy_level_simple
                                 + WildMeat
                                 
                                 + pubs_X_scale
                                 
                                 + PEV1
                                 + PEV2
                                 + PEV3
                                 + PEV4
                                 + PEV5
                                 + PEV6
                                 + PEV7
                                 + PEV8
                                 + PEV9
                                 + PEV10
                                 
                                 + (1|BioRealm)
                                 
                                 ,
                                 
                                 family=binomial(link="logit"),
                                 
                                 data= Dataset_01_m3,
                                 control = glmmTMBControl(parallel = 8)
)

summary(M1_01_virionnoPREDICT)
plot(simulateResiduals(fittedModel = M1_01_virionnoPREDICT))
tab_model(M1_01_virionnoPREDICT)


# Basic model, full CLOVER dataset
# Full model
effectsize_M1_01_CLOVER_raw_0.95 <- get_model_data(M1_01, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_01_CLOVER_raw_0.99 <- get_model_data(M1_01, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_01_CLOVER_raw <- effectsize_M1_01_CLOVER_raw_0.95 %>% 
  left_join(effectsize_M1_01_CLOVER_raw_0.99, by="term") %>% 
  mutate(pathogen_dataset= "CLOVER") 

# CLOVER viruses model
effectsize_M1_01_CLOVERvirus_raw_0.95 <- get_model_data(M1_01_clovervirus, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_01_CLOVERvirus_raw_0.99 <- get_model_data(M1_01_clovervirus, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_01_CLOVERvirus_raw <- effectsize_M1_01_CLOVERvirus_raw_0.95 %>% 
  left_join(effectsize_M1_01_CLOVERvirus_raw_0.99, by="term") %>% 
  mutate(pathogen_dataset= "CLOVER virus")

effectsize_M1_01_CLOVERvirus_raw$term = as.character(effectsize_M1_01_CLOVERvirus_raw$term)
effectsize_M1_01_CLOVERvirus_raw$term[1] = "Traded_allyes"
effectsize_M1_01_CLOVERvirus_raw$term = as.factor(effectsize_M1_01_CLOVERvirus_raw$term)

# CLOVER no EID2 model
effectsize_M1_01_CLOVERnoEID2_raw_0.95 <- get_model_data(M1_01_clovernoEID2, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_01_CLOVERnoEID2_raw_0.99 <- get_model_data(M1_01_clovernoEID2, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_01_CLOVERnoEID2_raw <- effectsize_M1_01_CLOVERnoEID2_raw_0.95 %>% 
  left_join(effectsize_M1_01_CLOVERnoEID2_raw_0.99, by="term") %>% 
  mutate(pathogen_dataset= "CLOVER noEID2")

effectsize_M1_01_CLOVERnoEID2_raw$term = as.character(effectsize_M1_01_CLOVERnoEID2_raw$term)
effectsize_M1_01_CLOVERnoEID2_raw$term[1] = "Traded_allyes"
effectsize_M1_01_CLOVERnoEID2_raw$term = as.factor(effectsize_M1_01_CLOVERnoEID2_raw$term)


# VIRION
effectsize_M1_01_VIRION_raw_0.95 <- get_model_data(M1_01_virion, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_01_VIRION_raw_0.99 <- get_model_data(M1_01_virion, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_01_VIRION_raw <- effectsize_M1_01_VIRION_raw_0.95 %>% 
  left_join(effectsize_M1_01_VIRION_raw_0.99, by="term") %>% 
  mutate(pathogen_dataset= "VIRION")

effectsize_M1_01_VIRION_raw$term = as.character(effectsize_M1_01_VIRION_raw$term)
effectsize_M1_01_VIRION_raw$term[1] = "Traded_allyes"
effectsize_M1_01_VIRION_raw$term = as.factor(effectsize_M1_01_VIRION_raw$term)


# VIRION no PREDICT
effectsize_M1_01_VIRIONnoPREDICT_raw_0.95 <- get_model_data(M1_01_virionnoPREDICT, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_01_VIRIONnoPREDICT_raw_0.99 <- get_model_data(M1_01_virionnoPREDICT, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_01_VIRIONnoPREDICT_raw <- effectsize_M1_01_VIRIONnoPREDICT_raw_0.95 %>% 
  left_join(effectsize_M1_01_VIRIONnoPREDICT_raw_0.99, by="term") %>% 
  mutate(pathogen_dataset= "VIRION noPREDICT")

effectsize_M1_01_VIRIONnoPREDICT_raw$term = as.character(effectsize_M1_01_VIRIONnoPREDICT_raw$term)
effectsize_M1_01_VIRIONnoPREDICT_raw$term[1] = "Traded_allyes"
effectsize_M1_01_VIRIONnoPREDICT_raw$term = as.factor(effectsize_M1_01_VIRIONnoPREDICT_raw$term)

# combine
effectsize_M1_01_pathogenDatasets_raw <- rbind(effectsize_M1_01_CLOVER_raw, 
                                               effectsize_M1_01_CLOVERvirus_raw,
                                               effectsize_M1_01_CLOVERnoEID2_raw,
                                               effectsize_M1_01_VIRION_raw,
                                               effectsize_M1_01_VIRIONnoPREDICT_raw)

effectsize_M1_01_pathogenDatasets_cond <- effectsize_M1_01_pathogenDatasets_raw %>% mutate(direction = ifelse(estimate>=1, "pos", "neg"))

# no manipulation for the conditional model
effectsize_M1_01_pathogenDatasets_cond <- effectsize_M1_01_pathogenDatasets_cond %>% 
  mutate(estimate_rescaled = estimate,
         conf.low.95_rescaled = conf.low.95,
         conf.high.95_rescaled = conf.high.95,
         conf.low.99_rescaled = conf.low.99,
         conf.high.99_rescaled = conf.high.99)

# then creating plotting values for visualisation
effectsize_M1_01_pathogenDatasets_cond <- effectsize_M1_01_pathogenDatasets_cond %>% 
  mutate(estimate_rescaled_2 = ifelse(direction=="pos", estimate_rescaled, (-1/estimate_rescaled)+2),
         conf.low.95_rescaled_2 = ifelse(direction=="pos", conf.low.95_rescaled, (-1/conf.low.95_rescaled)+2),
         conf.high.95_rescaled_2 = ifelse(direction=="pos", conf.high.95_rescaled, (-1/conf.high.95_rescaled)+2),
         conf.low.99_rescaled_2 = ifelse(direction=="pos", conf.low.99_rescaled, (-1/conf.low.99_rescaled)+2),
         conf.high.99_rescaled_2 = ifelse(direction=="pos", conf.high.99_rescaled, (-1/conf.high.99_rescaled)+2))


effectsize_M1_01_pathogenDatasets <- effectsize_M1_01_pathogenDatasets_cond
effectsize_M1_01_pathogenDatasets %>% colnames

effectsize_M1_01_pathogenDatasets %>% filter(term=="Traded_allyes") %>% dplyr::select(pathogen_dataset, estimate, p.value, p.stars) 
effectsize_M1_01_pathogenDatasets$term %>% table


pathogen_dataset_levels_all_02 <- c("VIRION noPREDICT", "VIRION", "CLOVER virus", "CLOVER noEID2", "CLOVER")

# Create full grid
ghost_grid_pathogenDatasets_all <- tidyr::expand_grid(term = term_levels_all,
                                                      pathogen_dataset = pathogen_dataset_levels_all_02)

# Merge with actual data
effectsize_complete_pathogenDatasets_all <- ghost_grid_pathogenDatasets_all %>%
  left_join(effectsize_M1_01_pathogenDatasets, by = c("term", "pathogen_dataset")) %>%
  mutate(term = factor(term, levels = term_levels_all),
         pathogen_dataset= factor(pathogen_dataset, levels = pathogen_dataset_levels_all_02),
         signif = ifelse(p.value<0.05, "yes", "no"),
         alpha_value = ifelse(p.value > 0.05 | is.na(p.value), 0.1, 
                              ifelse(p.value > 0.01, 0.6, 1)))  # ensure correct order

dodge_width_all <- 0.6

# Conditional
figS3b <- ggplot(data = effectsize_complete_pathogenDatasets_all,
                 aes(x = estimate, y = term, fill = term, shape= pathogen_dataset),
                 alpha=alpha_value) +
  geom_vline(xintercept = 1, linetype = 1, color = "black", size = 1) +
  geom_errorbarh(aes(xmin = conf.low.95_rescaled_2, xmax = conf.high.95_rescaled_2, color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width_all),
                 size = 1.25, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low.99_rescaled_2, xmax = conf.high.99_rescaled_2, color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width_all),
                 size = 0.5, na.rm = TRUE) +
  geom_point(aes(alpha = alpha_value),
             stroke=1,
             #shape = 21, 
             size = 2.5, 
             position = position_dodge(width = dodge_width_all),
             na.rm = TRUE) +
  scale_x_continuous(limits = c(0.5, 7.25),
                     breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7),
                     labels = c(3, 2, 1, 2, 3, 4, 5, 6, 7)) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25),
                     breaks=c("CLOVER", "CLOVER noEID2", "CLOVER virus", "VIRION", "VIRION noPREDICT"),
                     labels= c("CLOVER - Full", "CLOVER - No EID2", "CLOVER - Viruses only", "VIRION - Full", "VIRION - No PREDICT"),
                     name= "Host-pathogen \ndataset used") +
  # annotate("text", x=-Inf, y=Inf, label= "d", size=8, vjust = 1.5, hjust = -1) +
  theme_bw() +
  scale_y_discrete(labels = c("Traded_allyes" = "Traded",
                              "Synanthropy_level_simpleyes" = "Synanthropic",
                              "WildMeatyes" = "Used as\nwild meat",
                              "pubs_X_scale" = "Research\neffort")) +
  scale_colour_manual(values = custom_colors,
                      labels = c("Traded_allyes" = "Traded",
                                 "Synanthropy_level_simpleyes" = "Synanthropic",
                                 "WildMeatyes" = "Used as\nwild meat",
                                 "pubs_X_scale" = "Research\neffort")) +
  scale_fill_manual(values = custom_colors,
                    labels = c("Traded_allyes" = "Traded",
                               "Synanthropy_level_simpleyes" = "Synanthropic",
                               "WildMeatyes" = "Used as\nwild meat",
                               "pubs_X_scale" = "Research\neffort")) +
  scale_alpha_continuous(breaks= c(0.1, 0.5, 1),
                         labels= c("0.1" = "p > 0.05",
                                   "0.5" = "0.05 > p > 0.01",
                                   "1" = "p < 0.01"),
                         name= "Significance level") +
  labs(#title= "Probability of sharing\nat least one pathogen with humans",
    x= "Std effect size (odds ratio)") +
  guides(colour = "none",
         fill = "none")  +
  theme(legend.position = c(0.75, 0.75),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(size=10), 
        title = element_text(size=10), 
        # axis.text.y = element_blank(), 
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
  )

figS3b

figS3 <- ggpubr::ggarrange(figS3a,# + theme(legend.position = "none"), 
                           figS3b, # + theme(legend.position = "none"),
                           ncol=2)
figS3
system.time(ggsave("./outputs/figS3.pdf", figS3, width = 24, height = 18, units = "cm", limitsize = FALSE))




#################################################################
###########   Part 2 - Analyses linked to Figure 2    ########### 
#################################################################

fill_markets <- c(Legal = "#adc178ff",
                  Illegal = "#38a3a5ff",
                  Product = "#7e4618ff", 
                  Live = "#ffb057ff")

####> Figure 2a
euler_Trade_live_illegal <- list(
  Zoonotic = Dataset_01 %>% filter(Zoonose_vector_1 == "Zoonotic") %>% pull(MDD_id),
  Traded = Dataset_01 %>% filter(Traded_all =="yes") %>% pull(MDD_id),
  Live = Dataset_01 %>% filter(Traded_all_live =="yes") %>% pull(MDD_id),
  Illegal = Dataset_01 %>% filter(Traded_illegal =="yes") %>% pull(MDD_id),
  All = Dataset_01 %>% pull(MDD_id)
)

Fig2a <- plot(euler(euler_Trade_live_illegal, shape = "ellipse"),
              fill = alpha(c("#780000cc", "grey50", fill_markets[c(4,2)], "grey90"), alpha= 0.85),
              lwd=2.5,
              edges = c('black'),
              labels = list(cex = 0.5), 
              quantities = list(cex = 0.75))
Fig2a



####>>  Add illegal trade and live market to M1_01
M1_02 <- glmmTMB(Zoonose_vector_1 ~ 1 
                 + Traded_all
                 
                 + Traded_illegal
                 + Traded_all_live
                 
                 + Synanthropy_level_simple
                 + WildMeat
                 
                 + pubs_X_scale
                 
                 + (1|BioRealm)
                 
                 + PEV1
                 + PEV2
                 + PEV3
                 + PEV4
                 + PEV5
                 + PEV6
                 + PEV7
                 + PEV8
                 + PEV9
                 + PEV10
                 ,
                 
                 family=binomial(link="logit"),
                 
                 data= Dataset_01_m3, 
                 control = glmmTMBControl(parallel = 8)
)

summary(M1_02)
plot(simulateResiduals(M1_02))
performance::r2(M1_02)
performance::check_singularity(M1_02)
performance::check_collinearity(M1_02)

tab_model(M1_02)

# compare this expanded model to the previous one
anova(M1_01, M1_02)


# effect sizes
sdteff_M1_02_0.95 <- get_model_data(M1_02, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
sdteff_M1_02_0.99 <- get_model_data(M1_02, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

sdteff_M1_02 <- sdteff_M1_02_0.95 %>% 
  left_join(sdteff_M1_02_0.99, by="term")

stars_M1_02  <- sdteff_M1_02$p.stars[1:2]



###########. Phylogenetic analyses on model M1_02
##>>  testing for phylogenetic signal in residuals
sim_res_M1_02 <- simulateResiduals(fittedModel = M1_02)
residuals_M1_02 <- residuals(sim_res_M1_02)
names(residuals_M1_02) <- Dataset_01_m3$MDD_id 

# Match residuals to tree
residuals_M1_02.2 <- residuals_M1_02[names(residuals_M1_02) %in% tree_aug$tip.label]
consTree_M1_02 <- drop.tip(tree_aug, setdiff(tree_aug$tip.label, names(residuals_M1_02.2)))

# Test for phylogenetic signal
phySignal_M1_02.2 <- phylosig(consTree_M1_02, residuals_M1_02.2, method = "K", test=TRUE)  
phySignal_M1_02.2



### Plotting effect sizes
sdteff_M1_02 <- sdteff_M1_02 %>% 
  slice(1:6) %>% 
  mutate(signif = ifelse(p.value<0.05, "yes", "no"),
         alpha_value = ifelse(p.value > 0.05 | is.na(p.value), 0.1, 
                              ifelse(p.value > 0.01, 0.5, 1))) 
# plotting effect sizes
dodge_width=2

fig2_supp <- ggplot(data= sdteff_M1_02,
                    aes(x = estimate, y = term, alpha = alpha_value, color=term, fill= term))  +
  geom_vline(xintercept = 1, linetype = 1, color = "black", size = 1) +
  geom_errorbarh(aes(xmin = conf.low.95, xmax = conf.high.95),
                 height = 0, 
                 color="black",
                 position = position_dodge(width = dodge_width),
                 size = 1.25, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low.99, xmax = conf.high.99),
                 height = 0, 
                 color="black",
                 position = position_dodge(width = dodge_width),
                 size = 0.5, na.rm = TRUE) +
  geom_point(stroke=1,
             color="black",
             shape = 21, 
             size = 2.5, 
             na.rm = TRUE) +
  scale_x_continuous(limits=c(0.75, max(sdteff_M1_02$conf.high.99)+0.25),
                     breaks=c(-1, 0, 1, 2, 3, 4, 5, 6, 7),
                     labels=c(3, 2, 1, 2, 3, 4, 5, 6, 7)) +    
  scale_y_discrete(labels = c("Traded_allyes" = "Traded",
                              "Traded_illegalyes" = "Occurs in\nillegally trade",
                              "Traded_all_liveyes" = "Occurs in\nlive animal market",
                              "Synanthropy_level_simpleyes" = "Synanthropic",
                              "WildMeatyes" = "Used as\nwild meat",
                              "pubs_X_scale" = "Research\neffort")) +
  scale_fill_manual(values = custom_colors,
                    labels = c("Traded_allyes" = "Traded",
                               "Traded_illegalyes" = "Occurs in\nillegally trade",
                               "Traded_all_liveyes" = "Occurs in\nlive animal market",
                               "Synanthropy_level_simpleyes" = "Synanthropic",
                               "WildMeatyes" = "Used as\nwild meat",
                               "pubs_X_scale" = "Research\neffort")) +
  scale_alpha_continuous(range= c(0.1, 1),
                         breaks= c(0.1, 0.5, 1),
                         labels= c("0.1" = "p > 0.05",
                                   "0.5" = "0.05 > p > 0.01",
                                   "1" = "p < 0.01"),
                         name= "Significance level") +
  theme_bw() +
  theme(legend.position = c(0.75, 0.7),
        legend.key.size = unit(0.01, "cm"),
        legend.text = element_text(size=8),
        axis.text.x = element_text(size=10), 
        title = element_text(size=10), 
        # axis.text.y = element_blank(), 
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
  )
fig2_supp

system.time(ggsave("./outputs/figSupp_optional.pdf", fig2_supp, width = 6, height = 6, units = "cm", limitsize = FALSE))


###### >> Plot average predictions 

### Compute models effects
## Trade
M1_02_effect_trade <- ggaverage(M1_02, 
                                c("Traded_all"), 
                                type="response") %>% 
  as.data.frame
colnames(M1_02_effect_trade) <- c("traded", "predicted",  "std.error",  "conf.low", "conf.high")
M1_02_effect_illegal$traded <- factor(M1_02_effect_trade$traded, levels=c("no", "yes"))

# marginal risk ratio
marginalRR_traded <- M1_02_effect_trade$predicted[2]/M1_02_effect_trade$predicted[1] # 1.36

## Illegal trade
M1_02_effect_illegal <- ggaverage(M1_02, 
                                  c("Traded_illegal"), 
                                  type="response") %>% 
  as.data.frame
colnames(M1_02_effect_illegal) <- c("traded_illegal", "predicted",  "std.error",  "conf.low", "conf.high")
M1_02_effect_illegal$traded_illegal <- factor(M1_02_effect_illegal$traded_illegal, levels=c("no", "yes"))

# marginal risk ratio
marginalRR_illegal <- M1_02_effect_illegal$predicted[2]/M1_02_effect_illegal$predicted[1] # 1.11

## Live market
M1_02_effect_live <- ggaverage(M1_02, 
                               c("Traded_all_live"), 
                               type="response") %>% 
  as.data.frame
colnames(M1_02_effect_live) <- c("traded_all_live", "predicted",  "std.error",  "conf.low", "conf.high")
M1_02_effect_live$traded_live <- factor(M1_02_effect_live$traded_all_live, levels=c("no", "yes"))

# marginal risk ratio
marginalRR_live <- M1_02_effect_live$predicted[2]/M1_02_effect_live$predicted[1] # 1.34


# Plot model's effects
Fig2bc <-
  ggpubr::ggarrange(
    ggplot(data=M1_02_effect_live,
           aes(x=traded_all_live, y=predicted, ymin=conf.low, ymax=conf.high, fill=traded_all_live, col=traded_all_live)) +
      geom_linerange(position = position_dodge(.25), 
                     color="black",
                     size=2) +
      geom_point(position=position_dodge(width=0.25),
                 shape=22, 
                 size=6,
                 stroke=1,
                 col="black") +
      scale_fill_manual(values=c("grey75", as.character(fill_markets[4]))) +
      scale_colour_manual(values=c("grey75", as.character(fill_markets[4]))) +
      labs(fill="Live market", color="Live market") +
      ylab('Probability of sharing at least one pathogen with humans') +
      ylim(c(0.14, 0.28)) +
      geom_segment(aes(x=1, xend=2, y=0.275,yend=0.275), color="black", size=1) +
      annotate(geom="text", x=1.5, y=0.28, label=paste0(round(marginalRR_live, 1), stars_M1_02[1]),
               color="black") +
      theme_bw() +
      theme(legend.position = "none",
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text = element_blank()
      ),
    
    ggplot(data=M1_02_effect_illegal,
           aes(x=traded_illegal, y=predicted, ymin=conf.low, ymax=conf.high, fill=traded_illegal, col=traded_illegal)) +
      geom_linerange(position = position_dodge(.25), 
                     color="black",
                     size=2) +
      geom_point(position=position_dodge(width=0.25),
                 shape=22, 
                 size=6,
                 stroke=1,
                 col="black") +
      scale_fill_manual(values=c("grey75", as.character(fill_markets[2]))) +
      scale_colour_manual(values=c("grey75", as.character(fill_markets[2]))) +
      labs(fill="Illegal trade", color="Illegal trade") +
      ylab('Probability of sharing at least one pathogen with humans') +
      ylim(c(0.14, 0.28)) +
      geom_segment(aes(x=1, xend=2, y=0.275,yend=0.275), color="black", size=1) +
      annotate(geom="text", x=1.5, y=0.28, label=paste0(round(marginalRR_illegal, 1), stars_M1_02[2]),
               color="black") +
      theme_bw() +
      theme(legend.position = "none",
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text = element_blank()
      )
    
    ,
    
    ncol=2)

Fig2bc

blankP <- ggplot() + theme_void()

Fig2bc_2 <- ggpubr::ggarrange(blankP, Fig2bc, blankP, nrow=3, heights = c(0.1, 0.3, 0.1))

Fig2 <- ggpubr::ggarrange(Fig2a, Fig2bc_2, ncol=2)
Fig2

system.time(ggsave("./outputs/Fig2.pdf", Fig2, width = 15, height = 12, units = "cm", limitsize = FALSE))
############################




#################################################################
########             SENSITIVITY ANALYSES                 #######

### Sensitivity to the wildlife trade dataset
# CITES only
M1_02_CITES <- glmmTMB(Zoonose_vector_1 ~ 1 
                       + Traded_all_CITES
                       
                       + Traded_illegal
                       + Traded_live_CITES
                       
                       + Synanthropy_level_simple
                       + WildMeat
                       
                       + pubs_X_scale
                       
                       + (1|BioRealm)
                       
                       + PEV1
                       + PEV2
                       + PEV3
                       + PEV4
                       + PEV5
                       + PEV6
                       + PEV7
                       + PEV8
                       + PEV9
                       + PEV10
                       ,
                       
                       family=binomial(link="logit"),
                       
                       data= Dataset_01_m3,
                       control = glmmTMBControl(parallel = 8)
)

summary(M1_02_CITES)
plot(simulateResiduals(M1_02_CITES))
performance::check_singularity(M1_02_CITES)
performance::check_collinearity(M1_02_CITES)

tab_model(M1_02_CITES)

# LEMIS only
M1_02_LEMIS <- glmmTMB(Zoonose_vector_1 ~ 1 
                       + Traded_all_LEMIS2025
                       
                       + Traded_illegal
                       + Traded_live_LEMIS2025
                       
                       + Synanthropy_level_simple
                       + WildMeat
                       
                       + pubs_X_scale
                       
                       + (1|BioRealm)
                       
                       + PEV1
                       + PEV2
                       + PEV3
                       + PEV4
                       + PEV5
                       + PEV6
                       + PEV7
                       + PEV8
                       + PEV9
                       + PEV10
                       ,
                       
                       family=binomial(link="logit"),
                       
                       data= Dataset_01_m3,
                       control = glmmTMBControl(parallel = 8)
)

summary(M1_02_LEMIS)
plot(simulateResiduals(M1_02_LEMIS))
plot_model(M1_02_LEMIS, type="std2")
performance::check_singularity(M1_02_LEMIS)
performance::check_collinearity(M1_02_LEMIS)

tab_model(M1_02_LEMIS)


### plotting effect sizes for both submodels >>. !!!!! Re-run the models before running the following lines
effectsize_M1_02_BOTH_raw_0.95 <- get_model_data(M1_02, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_02_BOTH_raw_0.99 <- get_model_data(M1_02, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_02_BOTH_raw <- effectsize_M1_02_BOTH_raw_0.95 %>% 
  left_join(effectsize_M1_02_BOTH_raw_0.99, by="term") %>% 
  mutate(dataset= "BOTH")

# CITES only
effectsize_M1_02_CITES_raw_0.95 <- get_model_data(M1_02_CITES, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_02_CITES_raw_0.99 <- get_model_data(M1_02_CITES, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_02_CITES_raw <- effectsize_M1_02_CITES_raw_0.95 %>% 
  left_join(effectsize_M1_02_CITES_raw_0.99, by="term") %>% 
  mutate(dataset= "CITES")
effectsize_M1_02_CITES_raw$term = as.character(effectsize_M1_02_CITES_raw$term)
effectsize_M1_02_CITES_raw$term[1] = "Traded_allyes"
effectsize_M1_02_CITES_raw$term[3] = "Traded_all_liveyes"
effectsize_M1_02_CITES_raw$term = as.factor(effectsize_M1_02_CITES_raw$term)

# LEMIS only
effectsize_M1_02_LEMIS_raw_0.95 <- get_model_data(M1_02_LEMIS, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_02_LEMIS_raw_0.99 <- get_model_data(M1_02_LEMIS, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_02_LEMIS_raw <- effectsize_M1_02_LEMIS_raw_0.95 %>% 
  left_join(effectsize_M1_02_LEMIS_raw_0.99, by="term")%>% 
  mutate(dataset= "LEMIS")
effectsize_M1_02_LEMIS_raw$term = as.character(effectsize_M1_02_LEMIS_raw$term)
effectsize_M1_02_LEMIS_raw$term[1] = "Traded_allyes"
effectsize_M1_02_LEMIS_raw$term[3] = "Traded_all_liveyes"
effectsize_M1_02_LEMIS_raw$term = as.factor(effectsize_M1_02_LEMIS_raw$term)

effectsize_M1_02_raw <- rbind(effectsize_M1_02_BOTH_raw, 
                              effectsize_M1_02_CITES_raw,
                              effectsize_M1_02_LEMIS_raw)

effectsize_M1_02_cond <- effectsize_M1_02_raw %>% mutate(direction = ifelse(estimate>=1, "pos", "neg"))

# no manipulation for the conditional model
effectsize_M1_02_cond <- effectsize_M1_02_cond %>% 
  mutate(estimate_rescaled = estimate,
         conf.low.95_rescaled = conf.low.95,
         conf.high.95_rescaled = conf.high.95,
         conf.low.99_rescaled = conf.low.99,
         conf.high.99_rescaled = conf.high.99)

# then creating plotting values for visualisation
effectsize_M1_02_cond <- effectsize_M1_02_cond %>% 
  mutate(estimate_rescaled_2 = ifelse(direction=="pos", estimate_rescaled, (-1/estimate_rescaled)+2),
         conf.low.95_rescaled_2 = ifelse(direction=="pos", conf.low.95_rescaled, (-1/conf.low.95_rescaled)+2),
         conf.high.95_rescaled_2 = ifelse(direction=="pos", conf.high.95_rescaled, (-1/conf.high.95_rescaled)+2),
         conf.low.99_rescaled_2 = ifelse(direction=="pos", conf.low.99_rescaled, (-1/conf.low.99_rescaled)+2),
         conf.high.99_rescaled_2 = ifelse(direction=="pos", conf.high.99_rescaled, (-1/conf.high.99_rescaled)+2))


effectsize_M1_02 <- effectsize_M1_02_cond

effectsize_M1_02 %>% colnames

effectsize_M1_02 %>% filter(term=="Traded_all_liveyes") %>% dplyr::select(dataset, estimate, p.value, p.stars) 
effectsize_M1_02 %>% filter(term=="Traded_illegalyes") %>% dplyr::select(dataset, estimate, p.value, p.stars) 
effectsize_M1_02$term %>% table


# Desired terms (already in order)
term_levels_M1_02 <- c(
  "pubs_X_scale", 
  "WildMeatyes", 
  "Synanthropy_level_simpleyes",
  "Traded_all_liveyes",
  "Traded_illegalyes",
  "Traded_allyes")
dataset_levels_M1_02 <- c("LEMIS", "CITES", "BOTH")


# Create full grid
ghost_grid_M1_02 <- tidyr::expand_grid(term = term_levels_M1_02,
                                       dataset = dataset_levels_M1_02)

# Merge with actual data
effectsize_complete_M1_02 <- ghost_grid_M1_02 %>%
  left_join(effectsize_M1_02, by = c("term", "dataset")) %>%
  mutate(term = factor(term, levels = term_levels_M1_02),
         dataset= factor(dataset, levels = dataset_levels_M1_02),
         signif = ifelse(p.value<0.05, "yes", "no"),
         alpha_value = ifelse(p.value > 0.05 | is.na(p.value), 0.1, 
                              ifelse(p.value > 0.01, 0.5, 1)))  # ensure correct order

ci_min_M1_02 <- -12
dodge_width_M1_02 <- 0.6

# 
FigS4a <- ggplot(data = effectsize_complete_M1_02,
                 aes(x = estimate_rescaled_2, y = term, fill = term, shape= dataset, alpha= signif)) +
  geom_vline(xintercept = 1, linetype = 1, color = "black", size = 1) +
  geom_errorbarh(aes(xmin = conf.low.95_rescaled_2, xmax = conf.high.95_rescaled_2, color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width_M1_02),
                 size = 1.25, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low.99_rescaled_2, xmax = conf.high.99_rescaled_2, color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width_M1_02),
                 size = 0.5, na.rm = TRUE) +
  geom_point(aes(alpha = alpha_value),
             stroke=1,
             #shape = 21, 
             size = 2.5, 
             position = position_dodge(width = dodge_width_M1_02),
             na.rm = TRUE) +
  scale_x_continuous(limits = c(-0.75, 7),
                     breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7),
                     labels = c(3, 2, 1, 2, 3, 4, 5, 6, 7)) +
  scale_shape_manual(values=c(21, 22, 23),
                     breaks=c("BOTH", "CITES", "LEMIS"),
                     labels= c("CITES + LEMIS", "CITES only", "LEMIS only"),
                     name= "Wildlife trade \ndataset used") +
  # annotate("text", x=-Inf, y=Inf, label= "d", size=8, vjust = 1.5, hjust = -1) +
  theme_bw() +
  scale_y_discrete(labels = c("Traded_allyes" = "Traded",
                              "Traded_illegalyes" = "Occurs in\nillegally trade",
                              "Traded_all_liveyes" = "Occurs in\nlive animal market",
                              "Synanthropy_level_simpleyes" = "Synanthropic",
                              "WildMeatyes" = "Used as\nwild meat",
                              "pubs_X_scale" = "Research\neffort")) +
  scale_colour_manual(values = custom_colors,
                      labels = c("Traded_allyes" = "Traded",
                                 "Traded_illegalyes" = "Occurs in\nillegally trade",
                                 "Traded_all_liveyes" = "Occurs in\nlive animal market",
                                 "Synanthropy_level_simpleyes" = "Synanthropic",
                                 "WildMeatyes" = "Used as\nwild meat",
                                 "pubs_X_scale" = "Research\neffort")) +
  scale_fill_manual(values = custom_colors,
                    labels = c("Traded_allyes" = "Traded",
                               "Traded_illegalyes" = "Occurs in\nillegally trade",
                               "Traded_all_liveyes" = "Occurs in\nlive animal market",
                               "Synanthropy_level_simpleyes" = "Synanthropic",
                               "WildMeatyes" = "Used as\nwild meat",
                               "pubs_X_scale" = "Research\neffort")) +
  scale_alpha_continuous(breaks= c(0.1, 0.5, 1),
                         labels= c("0.1" = "p > 0.05",
                                   "0.5" = "0.05 > p > 0.01",
                                   "1" = "p < 0.01"),
                         name= "Significance level") +
  labs(#title= "Response: Probability of sharing\nat least one pathogen with humans",
    x= "Std effect size (odds ratio)") +
  guides(colour = "none",
         fill = "none")  +
  theme(legend.position = c(0.8, 0.75),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(size=10), 
        title = element_text(size=10), 
        # axis.text.y = element_blank(), 
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
  )

FigS4a



### Sensitivity on host-pathogen interaction dataset 
# Model M1_04 using CLOVER virus only
M1_02_clovervirus <- glmmTMB(Zoonose_virus_vector_1  ~ 1 
                             + Traded_all
                             
                             + Traded_illegal
                             + Traded_all_live
                             
                             + Synanthropy_level_simple
                             + WildMeat
                             
                             + pubs_X_scale
                             
                             + PEV1
                             + PEV2
                             + PEV3
                             + PEV4
                             + PEV5
                             + PEV6
                             + PEV7
                             + PEV8
                             + PEV9
                             + PEV10
                             
                             + (1|BioRealm)
                             
                             ,
                             
                             family=binomial(link="logit"),
                             
                             data= Dataset_01_m3,
                             control = glmmTMBControl(parallel = 8)
)

summary(M1_02_clovervirus)
plot(simulateResiduals(fittedModel = M1_02_clovervirus))
tab_model(M1_02_clovervirus)


# Model M1_04 using CLOVER virus only
M1_02_clovernoEID2 <- glmmTMB(noEID2_Zoonose_vector_1  ~ 1 
                              + Traded_all
                              
                              + Traded_illegal
                              + Traded_all_live
                              
                              + Synanthropy_level_simple
                              + WildMeat
                              
                              + pubs_X_scale
                              
                              + PEV1
                              + PEV2
                              + PEV3
                              + PEV4
                              + PEV5
                              + PEV6
                              + PEV7
                              + PEV8
                              + PEV9
                              + PEV10
                              
                              + (1|BioRealm)
                              
                              ,
                              
                              family=binomial(link="logit"),
                              
                              data= Dataset_01_m3,
                              control = glmmTMBControl(parallel = 8)
)


summary(M1_02_clovernoEID2)
plot(simulateResiduals(fittedModel = M1_02_clovernoEID2))
tab_model(M1_02_clovernoEID2)


# Model M1_03 using VIRION
M1_02_virion <- glmmTMB(VIRION_Zoonose_vector_1 ~ 1 
                        + Traded_all
                        
                        + Traded_illegal
                        + Traded_all_live
                        
                        + Synanthropy_level_simple
                        + WildMeat
                        
                        + pubs_X_scale
                        
                        + PEV1
                        + PEV2
                        + PEV3
                        + PEV4
                        + PEV5
                        + PEV6
                        + PEV7
                        + PEV8
                        + PEV9
                        + PEV10
                        
                        + (1|BioRealm)
                        , 
                        
                        family=binomial(link="logit"),
                        
                        data= Dataset_01_m3,
                        control = glmmTMBControl(parallel = 8)
)



summary(M1_02_virion)
plot(simulateResiduals(fittedModel = M1_02_virion))
tab_model(M1_02_virion)


# Model M1_03 using VIRION
M1_02_virionnoPREDICT <- glmmTMB(noPREDICT_VIRION_Zoonose_vector_1 ~ 1 
                                 + Traded_all
                                 
                                 + Traded_illegal
                                 + Traded_all_live
                                 
                                 + Synanthropy_level_simple
                                 + WildMeat
                                 
                                 + pubs_X_scale
                                 
                                 + PEV1
                                 + PEV2
                                 + PEV3
                                 + PEV4
                                 + PEV5
                                 + PEV6
                                 + PEV7
                                 + PEV8
                                 + PEV9
                                 + PEV10
                                 
                                 + (1|BioRealm)
                                 , 
                                 
                                 family=binomial(link="logit"),
                                 
                                 data= Dataset_01_m3,
                                 control = glmmTMBControl(parallel = 8)
)



summary(M1_02_virionnoPREDICT)
plot(simulateResiduals(fittedModel = M1_02_virionnoPREDICT))
tab_model(M1_02_virionnoPREDICT)


## Effect sizes
# Full model
effectsize_M1_02_CLOVER_raw_0.95 <- get_model_data(M1_02, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_02_CLOVER_raw_0.99 <- get_model_data(M1_02, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_02_CLOVER_raw <- effectsize_M1_02_CLOVER_raw_0.95 %>% 
  left_join(effectsize_M1_02_CLOVER_raw_0.99, by="term") %>% 
  mutate(pathogen_dataset= "CLOVER")

# CLOVER viruses model
effectsize_M1_02_CLOVERvirus_raw_0.95 <- get_model_data(M1_02_clovervirus, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_02_CLOVERvirus_raw_0.99 <- get_model_data(M1_02_clovervirus, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_02_CLOVERvirus_raw <- effectsize_M1_02_CLOVERvirus_raw_0.95 %>% 
  left_join(effectsize_M1_02_CLOVERvirus_raw_0.99, by="term") %>% 
  mutate(pathogen_dataset= "CLOVER virus")

effectsize_M1_02_CLOVERvirus_raw$term = as.character(effectsize_M1_02_CLOVERvirus_raw$term)
effectsize_M1_02_CLOVERvirus_raw$term = as.factor(effectsize_M1_02_CLOVERvirus_raw$term)


# CLOVER no EID2 model
effectsize_M1_02_CLOVERnoEID2_raw_0.95 <- get_model_data(M1_02_clovernoEID2, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_02_CLOVERnoEID2_raw_0.99 <- get_model_data(M1_02_clovernoEID2, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_02_CLOVERnoEID2_raw <- effectsize_M1_02_CLOVERnoEID2_raw_0.95 %>% 
  left_join(effectsize_M1_02_CLOVERnoEID2_raw_0.99, by="term") %>% 
  mutate(pathogen_dataset= "CLOVER noEID2")

effectsize_M1_02_CLOVERnoEID2_raw$term = as.character(effectsize_M1_02_CLOVERnoEID2_raw$term)
effectsize_M1_02_CLOVERnoEID2_raw$term = as.factor(effectsize_M1_02_CLOVERnoEID2_raw$term)


# VIRION
effectsize_M1_02_VIRION_raw_0.95 <- get_model_data(M1_02_virion, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_02_VIRION_raw_0.99 <- get_model_data(M1_02_virion, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_02_VIRION_raw <- effectsize_M1_02_VIRION_raw_0.95 %>% 
  left_join(effectsize_M1_02_VIRION_raw_0.99, by="term") %>% 
  mutate(pathogen_dataset= "VIRION")

effectsize_M1_02_VIRION_raw$term = as.character(effectsize_M1_02_VIRION_raw$term)
effectsize_M1_02_VIRION_raw$term = as.factor(effectsize_M1_02_VIRION_raw$term)


# VIRION no PREDICT
effectsize_M1_02_VIRIONnoPREDICT_raw_0.95 <- get_model_data(M1_02_virionnoPREDICT, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M1_02_VIRIONnoPREDICT_raw_0.99 <- get_model_data(M1_02_virionnoPREDICT, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M1_02_VIRIONnoPREDICT_raw <- effectsize_M1_02_VIRIONnoPREDICT_raw_0.95 %>% 
  left_join(effectsize_M1_02_VIRIONnoPREDICT_raw_0.99, by="term") %>% 
  mutate(pathogen_dataset= "VIRION noPREDICT")

effectsize_M1_02_VIRIONnoPREDICT_raw$term = as.character(effectsize_M1_02_VIRIONnoPREDICT_raw$term)
effectsize_M1_02_VIRIONnoPREDICT_raw$term = as.factor(effectsize_M1_02_VIRIONnoPREDICT_raw$term)


# combine dataframes
effectsize_M1_02_pathogenDatasets_raw <- rbind(effectsize_M1_02_CLOVER_raw, 
                                               effectsize_M1_02_CLOVERvirus_raw,
                                               effectsize_M1_02_CLOVERnoEID2_raw,
                                               effectsize_M1_02_VIRION_raw,
                                               effectsize_M1_02_VIRIONnoPREDICT_raw)

effectsize_M1_02_pathogenDatasets_cond <- effectsize_M1_02_pathogenDatasets_raw %>% mutate(direction = ifelse(estimate>=1, "pos", "neg"))

# no manipulation for the conditional model
effectsize_M1_02_pathogenDatasets_cond <- effectsize_M1_02_pathogenDatasets_cond %>% 
  mutate(estimate_rescaled = estimate,
         conf.low.95_rescaled = conf.low.95,
         conf.high.95_rescaled = conf.high.95,
         conf.low.99_rescaled = conf.low.99,
         conf.high.99_rescaled = conf.high.99)

# then creating plotting values for visualisation
effectsize_M1_02_pathogenDatasets_cond <- effectsize_M1_02_pathogenDatasets_cond %>% 
  mutate(estimate_rescaled_2 = ifelse(direction=="pos", estimate_rescaled, (-1/estimate_rescaled)+2),
         conf.low.95_rescaled_2 = ifelse(direction=="pos", conf.low.95_rescaled, (-1/conf.low.95_rescaled)+2),
         conf.high.95_rescaled_2 = ifelse(direction=="pos", conf.high.95_rescaled, (-1/conf.high.95_rescaled)+2),
         conf.low.99_rescaled_2 = ifelse(direction=="pos", conf.low.99_rescaled, (-1/conf.low.99_rescaled)+2),
         conf.high.99_rescaled_2 = ifelse(direction=="pos", conf.high.99_rescaled, (-1/conf.high.99_rescaled)+2))


effectsize_M1_02_pathogenDatasets <- effectsize_M1_02_pathogenDatasets_cond
effectsize_M1_02_pathogenDatasets %>% colnames

effectsize_M1_02_pathogenDatasets %>% filter(term=="Traded_all_liveyes") %>% dplyr::select(pathogen_dataset, estimate, p.value, p.stars) 
effectsize_M1_02_pathogenDatasets %>% filter(term=="Traded_illegalyes") %>% dplyr::select(pathogen_dataset, estimate, p.value, p.stars) 
effectsize_M1_02_pathogenDatasets$term %>% table


pathogen_dataset_levels_all_02 <- c("VIRION noPREDICT", "VIRION", "CLOVER virus", "CLOVER noEID2", "CLOVER")

# Create full grid
ghost_grid_M1_02_pathogenDatasets_all <- tidyr::expand_grid(term = term_levels_M1_02,
                                                            pathogen_dataset = pathogen_dataset_levels_all_02)

# Merge with actual data
effectsize_complete_M1_02_pathogenDatasets_all <- ghost_grid_M1_02_pathogenDatasets_all %>%
  left_join(effectsize_M1_02_pathogenDatasets, by = c("term", "pathogen_dataset")) %>%
  mutate(term = factor(term, levels = term_levels_M1_02),
         pathogen_dataset= factor(pathogen_dataset, levels = pathogen_dataset_levels_all_02),
         signif = ifelse(p.value<0.05, "yes", "no"),
         alpha_value = ifelse(p.value > 0.05 | is.na(p.value), 0.1, 
                              ifelse(p.value > 0.01, 0.6, 1)))  

dodge_width_M1_02 <- 0.6

# 
FigS4b <- ggplot(data = effectsize_complete_M1_02_pathogenDatasets_all,
                 aes(x = estimate, y = term, fill = term, shape= pathogen_dataset),
                 alpha=alpha_value) +
  geom_vline(xintercept = 1, linetype = 1, color = "black", size = 1) +
  geom_errorbarh(aes(xmin = conf.low.95_rescaled_2, xmax = conf.high.95_rescaled_2, color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width_M1_02),
                 size = 1.25, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low.99_rescaled_2, xmax = conf.high.99_rescaled_2, color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width_M1_02),
                 size = 0.5, na.rm = TRUE) +
  geom_point(aes(alpha = alpha_value),
             stroke=1,
             #shape = 21, 
             size = 2.5, 
             position = position_dodge(width = dodge_width_M1_02),
             na.rm = TRUE) +
  scale_x_continuous(limits = c(0.5, 7.25),
                     breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7),
                     labels = c(3, 2, 1, 2, 3, 4, 5, 6, 7)) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25),
                     breaks=c("CLOVER", "CLOVER noEID2", "CLOVER virus", "VIRION", "VIRION noPREDICT"),
                     labels= c("CLOVER - Full", "CLOVER - No EID2", "CLOVER - Viruses only", "VIRION - Full", "VIRION - No PREDICT"),
                     name= "Host-pathogen \ndataset used") +
  # annotate("text", x=-Inf, y=Inf, label= "d", size=8, vjust = 1.5, hjust = -1) +
  theme_bw() +
  scale_y_discrete(labels = c("Traded_allyes" = "Traded",
                              "Traded_illegalyes" = "Occurs in\nillegally trade",
                              "Traded_all_liveyes" = "Occurs in\nlive animal market",
                              "Synanthropy_level_simpleyes" = "Synanthropic",
                              "WildMeatyes" = "Used as\nwild meat",
                              "pubs_X_scale" = "Research\neffort")) +
  scale_colour_manual(values = custom_colors,
                      labels = c("Traded_allyes" = "Traded",
                                 "Traded_illegalyes" = "Occurs in\nillegally trade",
                                 "Traded_all_liveyes" = "Occurs in\nlive animal market",
                                 "Synanthropy_level_simpleyes" = "Synanthropic",
                                 "WildMeatyes" = "Used as\nwild meat",
                                 "pubs_X_scale" = "Research\neffort")) +
  scale_fill_manual(values = custom_colors,
                    labels = c("Traded_allyes" = "Traded",
                               "Traded_illegalyes" = "Occurs in\nillegally trade",
                               "Traded_all_liveyes" = "Occurs in\nlive animal market",
                               "Synanthropy_level_simpleyes" = "Synanthropic",
                               "WildMeatyes" = "Used as\nwild meat",
                               "pubs_X_scale" = "Research\neffort")) +
  scale_alpha_continuous(breaks= c(0.1, 0.5, 1),
                         labels= c("0.1" = "p > 0.05",
                                   "0.5" = "0.05 > p > 0.01",
                                   "1" = "p < 0.01"),
                         name= "Significance level") +
  labs(#title= "Response: Probability of sharing\nat least one pathogen with humans",
    x= "Std effect size (odds ratio)") +
  guides(colour = "none",
         fill = "none")  +
  theme(legend.position = c(0.75, 0.75),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(size=10), 
        title = element_text(size=10), 
        # axis.text.y = element_blank(), 
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
  )
FigS4b

FigS4 <- ggpubr::ggarrange(FigS4a, FigS4b, ncol=2)

system.time(ggsave("./outputs/FigS4.pdf", FigS4, width = 24, height = 18, units = "cm", limitsize = FALSE))
#################################







#################################################################
###########   Part 3 - Analyses linked to Figure 3    ########### 
#################################################################

Dataset_01 %>% colnames
### Representing which species were traded when based on CITES data, and different metrics for describing temporal trade trends
Dataset_02.1 <- Dataset_01 %>% 
  filter(timeTraded_all_1980_2019>0 & CITES_appendix %in% c("I", "II") & First_time_CITES<1981) %>% 
  mutate(across(starts_with("Traded_all_CITES_"),
                ~ ifelse(tolower(.) == "yes", 1, 0))) %>% 
  
  mutate(across(starts_with("Traded_live_CITES_"),
                ~ ifelse(tolower(.) == "yes", 1, 0))) %>% 
  mutate(across(starts_with("Traded_product_CITES_"),
                ~ ifelse(tolower(.) == "yes", 1, 0))) 
Dataset_02.1 %>% dim # 583 species
Dataset_02.1 %>% colnames


# All trade
Dataset_02.1$timeTraded_all <- rowSums(Dataset_02.1[, 37:76], na.rm = T)
# # timeTraded_all should be equal to timeTraded_all_1980_2019
# ggplot(data= Dataset_02.1) +
#   geom_point(aes(x=timeTraded_all, y=timeTraded_all_1980_2019))


Dataset_02.4_m3 <- Dataset_02.1 %>% mutate(
  ZOONOTIC_NB_all_log = log(ZOONOTIC_NB_all +1),
  timeTraded_all_scale = as.numeric(scale(timeTraded_all)),
  pubs_X_scale = as.numeric(scale(pubs_X)),
  PhyloDist_Homosapiens_mean_scale = as.numeric(scale(PhyloDist_Homosapiens_mean))
)



### Phylogenetic eigenvectors on the rediced dataset
DatasetP3_tipset <- unique(Dataset_02.4_m3$MDD_id)
DatasetP3_complete_TREE <- drop.tip(
  tree_aug,
  setdiff(tree_aug$tip.label, DatasetP3_tipset)
)

#
pvr_02 <- PVRdecomp(DatasetP3_complete_TREE)
EV_02  <- pvr_02@Eigen$vectors
rownames(EV_02) <- pvr_02@phylo$tip.label

eigvals2 <- pvr_02@Eigen$values
cumsum(eigvals2) / sum(eigvals2) 
plot(eigvals2, type="b", main="Phylogenetic eigenvalues", ylab="Eigenvalue", xlim=c(0,25))

k <- 10 
EV_df2 <- as.data.frame(EV_02[, 1:k, drop = FALSE])
EV_df2$MDD_id <- rownames(EV_df2)
names(EV_df2)[1:k] <- paste0("PEVbis", 1:k)

EV_df_keyed2 <- EV_df2 %>%
  mutate(MDD_id_chr = trimws(as.character(MDD_id))) %>%
  dplyr::select(MDD_id_chr, starts_with("PEVbis")) %>%
  distinct(MDD_id_chr, .keep_all = TRUE) 

DatasetP3_complete_keyed <- Dataset_02.4_m3 %>%
  mutate(MDD_id_chr = trimws(as.character(MDD_id))) %>%
  filter(MDD_id_chr %in% rownames(EV_02))

Dataset_02.4_m3 <- DatasetP3_complete_keyed %>%
  left_join(EV_df_keyed2, by = "MDD_id_chr")

pev_cols2 <- paste0("PEVbis", 1:k)
Dataset_02.4_m3[pev_cols2] <- scale(Dataset_02.4_m3[pev_cols2])
###############




####>   Figure 3a  
# Data preparation for plotting
mainOrders_22 <- Dataset_02.4_m3$order %>% table %>% sort(decreasing = T) %>% as.data.frame() %>% filter(Freq>50) %>% pull(".")
Dataset_02.4_m3 <- Dataset_02.4_m3 %>% mutate(order2= ifelse(order %in% mainOrders_22, order, "other"))

Dataset_02.4_m3$Zoonose_vector_1 <- factor(Dataset_02.4_m3$Zoonose_vector_1, levels=rev(levels(Dataset_02.4_m3$Zoonose_vector_1)))

# calculating first year in trade for sorting species in fig 3a
Dataset_02.4_m3$timeFirstTrade <- NA
for (species in Dataset_02.4_m3$MDD_id){
  Dataset_02.4_m3$timeFirstTrade[Dataset_02.4_m3$MDD_id==species] = 1979+min(which(Dataset_02.4_m3[Dataset_02.4_m3$MDD_id==species, 37:76]==1))
}

Dataset_02.4_m3$timeFirstTrade %>% hist

Dataset_02.4_m3_2 <- Dataset_02.4_m3 %>% 
  arrange(
    desc(timeTraded_all),
    timeFirstTrade,
    family,
    genus
  ) %>% 
  mutate(Synanthropy_bin = ifelse(Synanthropy_level_simple=="yes", 1, 0),
         WildMeat_bin = ifelse(WildMeat=="yes", 1, 0),
         Traded_illegal_bin = ifelse(Traded_illegal=="yes", 1, 0),
         Traded_all_live_bin = ifelse(Traded_all_live=="yes", 1, 0))

allTrade_matrix_2 <- as.matrix(Dataset_02.4_m3_2[,37:76], ncol=ncol(Dataset_02.4_m3_2[,37:76]))
colnames(allTrade_matrix_2) <- 1980:2019
rownames(allTrade_matrix_2) <- Dataset_02.4_m3_2$MDD_id

# colors barplot nb years
color_map_TimeinTrade <- paletteer::paletteer_c("grDevices::Rocket", length(rowSums(allTrade_matrix_2)/40), direction = -1)
colBars_TimeinTrade <- color_map_TimeinTrade[cut(rowSums(allTrade_matrix_2)/40, length(rowSums(allTrade_matrix_2)))]

# colors barplot nb zoonoses
color_map_NbZoonoses <- paletteer::paletteer_c("grDevices::OrRd", 9, direction = -1)
colBars_NbZoonoses = circlize::colorRamp2(seq(from=0, to=log(40), length.out=9), color_map_NbZoonoses)


fig3a <- Heatmap(allTrade_matrix_2,
                   col=c("grey99", "black"),
                   column_order= colnames(allTrade_matrix_2),
                   column_names_gp = gpar(fontsize = 4),
                   show_heatmap_legend = FALSE,
                   column_names_rot=70,
                   cluster_rows = F,
                   row_labels = rep('', nrow(allTrade_matrix_2)),
                   cluster_row_slices = FALSE,
                   row_title_gp  = gpar(col = c("#780000", "black"),
                                        fontsize = c(10, 10)),
                   right_annotation = rowAnnotation(
                     NbZoo = anno_simple(Dataset_02.4_m3_2$ZOONOTIC_NB_all_log,
                                         width = unit(0.4, "cm"),
                                         col=colBars_NbZoonoses),
                     Synanthropy = anno_simple(Dataset_02.4_m3_2$Synanthropy_bin,
                                               width = unit(0.2, "cm"),
                                               col=c("1"="#2cacdf9a", "0"="white")),
                     WilMeat = anno_simple(Dataset_02.4_m3_2$WildMeat_bin,
                                           width = unit(0.2, "cm"),
                                           col=c("1"="#e3d585ff", "0"="white")),
                     TradedIllegal = anno_simple(Dataset_02.4_m3_2$Traded_illegal_bin,
                                                 width = unit(0.2, "cm"),
                                                 col=c("1"="#38a3a5ff", "0"="white")),
                     TradedLive = anno_simple(Dataset_02.4_m3_2$Traded_all_live_bin,
                                              width = unit(0.2, "cm"),
                                              col=c("1"="#ffb057", "0"="white")),
                     TimeInTrade = anno_barplot(rowSums(allTrade_matrix_2),
                                                gp = gpar(col = colBars_TimeinTrade),
                                                width = unit(1, "cm"))
                   )
)

anno_legend <- Legend(title = "NbZoo", 
                      col_fun = colBars_NbZoonoses, 
                      at = log(c(0, 1, 2, 5, 10, 20, 50)+1),
                      labels = c(0, 1, 2, 5, 10, 20, 50),
                      labels_gp = gpar(fontsize = 10))

draw(fig3a, annotation_legend_list = list(anno_legend))
#####################




####>> Negative binomial model

# 1) FULL CLOVER
M2_01_clover <- glmmTMB(ZOONOTIC_NB_all ~ 1
                        + timeTraded_all_scale
                        
                        + Synanthropy_level_simple
                        + WildMeat
                        + Traded_all_live
                        + Traded_illegal
                        
                        + pubs_X_scale
                        
                        + PEVbis1
                        + PEVbis2
                        + PEVbis3
                        + PEVbis4
                        + PEVbis5
                        + PEVbis6
                        + PEVbis7
                        + PEVbis8
                        + PEVbis9
                        + PEVbis10
                        
                        # + (1|BioRealm)
                        ,
                        
                        family= nbinom2(link="log"), 
                        
                        data= Dataset_02.4_m3,
                        control = glmmTMBControl(parallel = 8)
)


summary(M2_01_clover)

# MuMIn::r.squaredGLMM(M2_01_clover)

M2_01_clover_res <- simulateResiduals(M2_01_clover)
plot(M2_01_clover_res)

performance::check_collinearity(M2_01_clover)
performance::check_singularity(M2_01_clover)
performance::check_overdispersion(M2_01_clover)
performance::check_zeroinflation(M2_01_clover)

plot_model(M2_01_clover, type="std2")

tab_model(M2_01_clover)


##>>  testing for phylogenetic signal in residuals
residuals_M2_01_clover <- residuals(M2_01_clover_res)
names(residuals_M2_01_clover) <- Dataset_02.4_m3$MDD_id  

# Match residuals to tree
residuals_M2_01_clover <- residuals_M2_01_clover[names(residuals_M2_01_clover) %in% tree_aug$tip.label]
consTree_M2_01_clover <- drop.tip(tree_aug, setdiff(tree_aug$tip.label, names(residuals_M2_01_clover))) # 583 tips

# Test for phylogenetic signal
phySignal_M2_01_clover <- phylosig(consTree_M2_01_clover, residuals_M2_01_clover, method = "K", test=TRUE)  # Blomberg's K
phySignal_M2_01_clover



# 2) no EID2 CLOVER
M2_01_clover_noEID2 <- glmmTMB(noEID2_ZOONOTIC_NB_all ~ 1
                               + timeTraded_all_scale
                               
                               + Synanthropy_level_simple
                               + WildMeat
                               + Traded_all_live
                               + Traded_illegal
                               
                               + pubs_X_scale
                               
                               + PEVbis1
                               + PEVbis2
                               + PEVbis3
                               + PEVbis4
                               + PEVbis5
                               + PEVbis6
                               + PEVbis7
                               + PEVbis8
                               + PEVbis9
                               + PEVbis10
                               
                               # + (1|BioRealm)
                               ,

                               family= nbinom2(link="log"), 
                               
                               data= Dataset_02.4_m3,
                               control = glmmTMBControl(parallel = 8)
)


summary(M2_01_clover_noEID2)

M2_01_clover_noEID2_res <- simulateResiduals(M2_01_clover_noEID2)
plot(M2_01_clover_noEID2_res)

performance::check_collinearity(M2_01_clover_noEID2)
performance::check_singularity(M2_01_clover_noEID2)
performance::check_overdispersion(M2_01_clover_noEID2)
performance::check_zeroinflation(M2_01_clover_noEID2)

plot_model(M2_01_clover_noEID2, type="std2")


# 3) VIRUS-ONLY CLOVER
M2_01_clover_viruses <- glmmTMB(ZOONOTIC_NB_viruses ~ 1
                                + timeTraded_all_scale
                                
                                + Synanthropy_level_simple
                                + WildMeat
                                + Traded_all_live
                                + Traded_illegal
                                
                                + pubs_X_scale
                                
                                + PEVbis1
                                + PEVbis2
                                + PEVbis3
                                + PEVbis4
                                + PEVbis5
                                + PEVbis6
                                + PEVbis7
                                + PEVbis8
                                + PEVbis9
                                + PEVbis10
                                
                                # + (1|BioRealm)
                                ,
                                
                                family= nbinom2(link="log"), 
                                
                                data= Dataset_02.4_m3,
                                control = glmmTMBControl(parallel = 8)
)


summary(M2_01_clover_viruses)

M2_01_clover_viruses_res <- simulateResiduals(M2_01_clover_viruses)
plot(M2_01_clover_viruses_res)

performance::check_collinearity(M2_01_clover_viruses)
performance::check_singularity(M2_01_clover_viruses)
performance::check_overdispersion(M2_01_clover_viruses)
performance::check_zeroinflation(M2_01_clover_viruses)

plot_model(M2_01_clover_viruses, type="std2")


# 4) FULL VIRION
M2_01_virion <- glmmTMB(VIRION_ZOONOTIC_NB_viruses ~ 1
                        + timeTraded_all_scale
                        
                        + Synanthropy_level_simple
                        + WildMeat
                        + Traded_all_live
                        + Traded_illegal
                        
                        + pubs_X_scale
                        
                        + PEVbis1
                        + PEVbis2
                        + PEVbis3
                        + PEVbis4
                        + PEVbis5
                        + PEVbis6
                        + PEVbis7
                        + PEVbis8
                        + PEVbis9
                        + PEVbis10
                        
                        # + (1|BioRealm)
                        ,
                        
                        family= nbinom2(link="log"), 
                        
                        data= Dataset_02.4_m3, 
                        control = glmmTMBControl(parallel = 8)
)

summary(M2_01_virion)

M2_01_virion_res <- simulateResiduals(M2_01_virion)
plot(M2_01_virion_res)

performance::check_collinearity(M2_01_virion)
performance::check_singularity(M2_01_virion)
performance::check_overdispersion(M2_01_virion)
performance::check_zeroinflation(M2_01_virion)

plot_model(M2_01_virion, type="std2")


##>>  testing for phylogenetic signal in residuals
residuals_M2_01_virion <- residuals(M2_01_virion_res)
names(residuals_M2_01_virion) <- Dataset_02.4_m3$MDD_id 

# Match residuals to tree
residuals_M2_01_virion <- residuals_M2_01_virion[names(residuals_M2_01_virion) %in% tree_aug$tip.label]
consTree_M2_01_virion <- drop.tip(tree_aug, setdiff(tree_aug$tip.label, names(residuals_M2_01_virion))) # 583 tips

# Test for phylogenetic signal
phySignal_M2_01_virion <- phylosig(consTree_M2_01_virion, residuals_M2_01_virion, method = "K", test=TRUE)  # Blomberg's K
phySignal_M2_01_virion


# 5) FULL VIRION
M2_01_noPREDICT_virion <- glmmTMB(noPREDICT_VIRION_ZOONOTIC_NB_viruses ~ 1
                                  + timeTraded_all_scale
                                  
                                  + Synanthropy_level_simple
                                  + WildMeat
                                  + Traded_all_live
                                  + Traded_illegal
                                  
                                  + pubs_X_scale
                                  
                                  + PEVbis1
                                  + PEVbis2
                                  + PEVbis3
                                  + PEVbis4
                                  + PEVbis5
                                  + PEVbis6
                                  + PEVbis7
                                  + PEVbis8
                                  + PEVbis9
                                  + PEVbis10
                                  
                                  # + (1|BioRealm)
                                  ,
                                  
                                  family= nbinom2(link="log"), 
                                  
                                  data= Dataset_02.4_m3,
                                  control = glmmTMBControl(parallel = 8)
)


summary(M2_01_noPREDICT_virion)

M2_01_noPREDICT_virion_res <- simulateResiduals(M2_01_noPREDICT_virion)
plot(M2_01_noPREDICT_virion_res)

performance::check_collinearity(M2_01_noPREDICT_virion)
performance::check_singularity(M2_01_noPREDICT_virion)
performance::check_overdispersion(M2_01_noPREDICT_virion)
performance::check_zeroinflation(M2_01_noPREDICT_virion)

plot_model(M2_01_noPREDICT_virion, type="std2")



###   Figure 3, panels b, c, d, e and supp associated

## plotting effect sizes for both submodels >>. !!!!! Re-run the models before running the following lines
# CLOVER full
effectsize_M2_01_clover_raw_0.95 <- get_model_data(M2_01_clover, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M2_01_clover_raw_0.99 <- get_model_data(M2_01_clover, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M2_01_clover_raw <- effectsize_M2_01_clover_raw_0.95 %>% 
  left_join(effectsize_M2_01_clover_raw_0.99, by="term") %>% 
  mutate(dataset= "CLOVER")

# CLOVER - no EID2
effectsize_M2_01_clover_noEID2_0.95 <- get_model_data(M2_01_clover_noEID2, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M2_01_clover_noEID2_0.99 <- get_model_data(M2_01_clover_noEID2, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M2_01_clover_noEID2_raw <- effectsize_M2_01_clover_noEID2_0.95 %>% 
  left_join(effectsize_M2_01_clover_noEID2_0.99, by="term") %>% 
  mutate(dataset= "CLOVER-noEID2")


# CLOVER - viruses only
effectsize_M2_01_clover_viruses_0.95 <- get_model_data(M2_01_clover_viruses, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M2_01_clover_viruses_0.99 <- get_model_data(M2_01_clover_viruses, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M2_01_clover_viruses_raw <- effectsize_M2_01_clover_viruses_0.95 %>% 
  left_join(effectsize_M2_01_clover_viruses_0.99, by="term") %>% 
  mutate(dataset= "CLOVER-viruses")


# VIRION full
effectsize_M2_01_virion_raw_0.95 <- get_model_data(M2_01_virion, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M2_01_virion_raw_0.99 <- get_model_data(M2_01_virion, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M2_01_virion_raw <- effectsize_M2_01_virion_raw_0.95 %>% 
  left_join(effectsize_M2_01_virion_raw_0.99, by="term") %>% 
  mutate(dataset= "VIRION")

# VIRION - no PREDICT
effectsize_M2_01_virion_noPredict_raw_0.95 <- get_model_data(M2_01_noPREDICT_virion, type="std2", ci.lvl = 0.95) %>% 
  dplyr::rename(conf.low.95  = conf.low,
                conf.high.95 = conf.high)
effectsize_M2_01_virion_noPredict_raw_0.99 <- get_model_data(M2_01_noPREDICT_virion, type="std2", ci.lvl = 0.99)  %>% 
  dplyr::rename(conf.low.99  = conf.low,
                conf.high.99 = conf.high) %>% 
  dplyr::select(term, conf.low.99, conf.high.99)

effectsize_M2_01_virion_noPredict_raw <- effectsize_M2_01_virion_noPredict_raw_0.95 %>% 
  left_join(effectsize_M2_01_virion_noPredict_raw_0.99, by="term") %>% 
  mutate(dataset= "VIRION-noPredict")


effectsize_M2_01_raw <- rbind(effectsize_M2_01_clover_raw, 
                              effectsize_M2_01_clover_viruses_raw,
                              effectsize_M2_01_clover_noEID2_raw,
                              effectsize_M2_01_virion_raw,
                              effectsize_M2_01_virion_noPredict_raw)


effectsize_M2_01_raw <- effectsize_M2_01_raw %>% mutate(direction = ifelse(estimate>=1, "pos", "neg"))

# no manipulation for the conditional model
effectsize_M2_01_cond <- effectsize_M2_01_raw %>% 
  mutate(estimate_rescaled = estimate,
         conf.low.95_rescaled = conf.low.95,
         conf.high.95_rescaled = conf.high.95,
         conf.low.99_rescaled = conf.low.99,
         conf.high.99_rescaled = conf.high.99)

# then creating plotting values for visualisation
effectsize_M2_01_cond <- effectsize_M2_01_cond %>% 
  mutate(estimate_rescaled_2 = ifelse(direction=="pos", estimate_rescaled, (-1/estimate_rescaled)+2),
         conf.low.95_rescaled_2 = ifelse(direction=="pos", conf.low.95_rescaled, (-1/conf.low.95_rescaled)+2),
         conf.high.95_rescaled_2 = ifelse(direction=="pos", conf.high.95_rescaled, (-1/conf.high.95_rescaled)+2),
         conf.low.99_rescaled_2 = ifelse(direction=="pos", conf.low.99_rescaled, (-1/conf.low.99_rescaled)+2),
         conf.high.99_rescaled_2 = ifelse(direction=="pos", conf.high.99_rescaled, (-1/conf.high.99_rescaled)+2))

effectsize_M2_01 <- effectsize_M2_01_cond 

effectsize_M2_01 %>% filter(term=="timeTraded_all_scale") %>% dplyr::select(dataset, estimate, p.value, p.stars) 
effectsize_M2_01$term %>% table


# Desired terms (already in order)
term_levels <- c("pubs_X_scale", "WildMeatyes", "Synanthropy_level_simpleyes", 
                 "Traded_all_liveyes", "Traded_illegalyes", "timeTraded_all_scale")
dataset_levels <- c("VIRION-noPredict", "VIRION", "CLOVER-viruses", "CLOVER-noEID2", "CLOVER")

# Datasets and facets
datasets <- rev(dataset_levels)

# Create full grid
ghost_grid <- tidyr::expand_grid(term = term_levels,
                                 dataset = datasets)

# Merge with actual data
effectsize_complete <- ghost_grid %>%
  left_join(effectsize_M2_01, by = c("term", "dataset")) %>%
  mutate(term = factor(term, levels = term_levels),
         dataset= factor(dataset, levels = dataset_levels),
         signif = ifelse(p.value<0.05, "yes", "no"),
         alpha_value = ifelse(p.value > 0.05 | is.na(p.value), 0.1, 
                              ifelse(p.value > 0.01, 0.5, 1)))  

ci_min <- -12  
dodge_width <- 0.6

# Conditional
figS5 <- ggplot(data = effectsize_complete,
                aes(x = estimate_rescaled_2, y = term, fill = term, shape= dataset)) +
  geom_vline(xintercept = 1, linetype = 1, color = "black", size = 1) +
  geom_errorbarh(aes(xmin = conf.low.95_rescaled_2, xmax = conf.high.95_rescaled_2, 
                     color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width),
                 size = 1.25, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low.99_rescaled_2, xmax = conf.high.99_rescaled_2, 
                     color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width),
                 size = 0.5, na.rm = TRUE) +
  geom_point(aes(alpha = alpha_value),
             stroke=1,
             #shape = 21, 
             size = 2.5, 
             position = position_dodge(width = dodge_width),
             na.rm = TRUE) +
  scale_x_continuous(limits = c(-0.1, 3.75),
                     breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
                     labels = c(3.5, 3, 2.5, 2, 1.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25),
                     breaks=c("CLOVER", "CLOVER-noEID2", "CLOVER-viruses", "VIRION", "VIRION-noPredict"),
                     labels= c("CLOVER - Full", "CLOVER - No EID2", "CLOVER - Viruses only", "VIRION - Full", "VIRION - No PREDICT"),
                     name= "Dataset used") +
  theme_bw() +
  scale_y_discrete(labels = c("timeTraded_all_scale" = "Time in \ntrade",
                              "Traded_illegalyes" = "Occurs in\nillegally trade",
                              "Traded_all_liveyes" = "Occurs in\nlive animal market",
                              "Synanthropy_level_simpleyes" = "Synanthropic",
                              "WildMeatyes" = "Used as\nwild meat",
                              "pubs_X_scale" = "Research\neffort")) +
  scale_colour_manual(values = custom_colors,
                      labels = c("timeTraded_all_scale" = "Time in \ntrade",
                                 "Traded_illegalyes" = "Occurs in\nillegally trade",
                                 "Traded_all_liveyes" = "Occurs in\nlive animal market",
                                 "Synanthropy_level_simpleyes" = "Synanthropic",
                                 "WildMeatyes" = "Used as\nwild meat",
                                 "pubs_X_scale" = "Research\neffort")) +
  scale_fill_manual(values = custom_colors,
                    labels = c("timeTraded_all_scale" = "Time in \ntrade",
                               "Traded_illegalyes" = "Occurs in\nillegally trade",
                               "Traded_all_liveyes" = "Occurs in\nlive animal market",
                               "Synanthropy_level_simpleyes" = "Synanthropic",
                               "WildMeatyes" = "Used as\nwild meat",
                               "pubs_X_scale" = "Research\neffort")) +
  scale_alpha_continuous(range= c(0.1, 1),
                         breaks= c(0.1, 0.5, 1),
                         labels= c("0.1" = "p > 0.05",
                                   "0.5" = "0.05 > p > 0.01",
                                   "1" = "p < 0.01"),
                         name= "Significance level") +
  labs(# title= "Response: Number of pathogens\nshared with humans",
    x= "Standardized effect size") +
  guides(colour = "none",
         fill = "none")  +
  theme(legend.position = c(0.75, 0.35),
        legend.key.size = unit(0.01, "cm"),
        legend.text = element_text(size=8),
        axis.text.x = element_text(size=10), 
        title = element_text(size=10), 
        axis.title.y = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
  )

figS5

system.time(ggsave("../outputs/figS5.pdf", figS5, width = 12, height = 15, units = "cm", limitsize = FALSE))


#####################


##>   Figure 3  <<## ## CLOVER model only ##

# Main predictor
avg_effects_clover_timeTraded_all_scale <- ggeffects::ggaverage(M2_01_clover, terms=c('timeTraded_all_scale [all]'), type="response")
avg_effects_clover_timeTraded_all_scale_df <-  avg_effects_clover_timeTraded_all_scale %>% data.frame()
avg_effects_clover_timeTraded_all_scale_df
avg_effects_clover_timeTraded_all_scale_df <- avg_effects_clover_timeTraded_all_scale_df %>% mutate(x_unscaled = scales::rescale(x, to=c(1, 40)))

# Live trade effect
avg_effects_clover_Traded_live <- ggeffects::ggaverage(M2_01_clover, terms=c('Traded_all_live [all]'), type="response")
avg_effects_clover_Traded_live_df <-  avg_effects_clover_Traded_live %>% data.frame()
avg_effects_clover_Traded_live_df
colnames(avg_effects_clover_Traded_live_df)[1] <- "Traded_live"
marginalRR_avg_effects_clover_Traded_live_df <- avg_effects_clover_Traded_live_df$predicted[2]/avg_effects_clover_Traded_live_df$predicted[1] # 1.49

# illegal trade effect
avg_effects_clover_Trade_illegal <- ggeffects::ggaverage(M2_01_clover, terms=c('Traded_illegal [all]'), type="response")
avg_effects_clover_Trade_illegal_df <-  avg_effects_clover_Trade_illegal %>% data.frame()
avg_effects_clover_Trade_illegal_df
colnames(avg_effects_clover_Trade_illegal_df)[1] <- "Trade_illegal"
marginalRR_avg_effects_clover_Trade_illegal_df <- avg_effects_clover_Trade_illegal_df$predicted[2]/avg_effects_clover_Trade_illegal_df$predicted[1] # 1.37



## Plotting the panels
fig3b <- ggplot(data = filter(effectsize_complete %>% filter(dataset=="CLOVER")),
                aes(x = estimate_rescaled_2, y = term, fill = term)) +
  geom_vline(xintercept = 1, linetype = 1, color = "black", size = 1) +
  geom_errorbarh(aes(xmin = conf.low.95_rescaled_2, xmax = conf.high.95_rescaled_2, color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width),
                 size = 1.25, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low.99_rescaled_2, xmax = conf.high.99_rescaled_2, color = term, alpha= alpha_value),
                 height = 0, 
                 position = position_dodge(width = dodge_width),
                 size = 0.5, na.rm = TRUE) +
  geom_point(aes(alpha = alpha_value),
             stroke=1,
             shape = 21, 
             size = 2.5, 
             position = position_dodge(width = dodge_width),
             na.rm = TRUE) +
  scale_x_continuous(limits = c(-0.1, 3.25),
                     breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3),
                     labels = c(3.5, 3, 2.5, 2, 1.5, 1, 1.5, 2, 2.5, 3)) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  annotate("text", x=Inf, y=Inf, label= "b", size=8, vjust = 1.5, hjust = 2) +
  theme_bw() +
  scale_y_discrete(labels = c("timeTraded_all_scale" = "Time in \n trade",
                              "Synanthropy_level_simpleyes" = "Synanthropy",
                              "WildMeatyes" = "Used as \n wild meat",
                              "Traded_illegalyes" = "Illegal \n trade",
                              "Traded_all_liveyes" = "Live \n market",
                              "pubs_X_scale" = "Research \n effort")) +
  scale_alpha_continuous(breaks= c(0.1, 0.5, 1),
                         labels= c("0.1" = "p > 0.05",
                                   "0.5" = "0.05 > p > 0.01",
                                   "1" = "p < 0.01"),
                         name= "Significance level") +
  guides(colour = "none",
         fill = "none")  +
  theme(legend.position = c(0.75, 0.7),
        legend.key.size = unit(0.01, "cm"),
        legend.text = element_text(size=8),
        axis.text.x = element_text(size=10), 
        title = element_text(size=10), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
  )

fig3b

# For visualisation purpose, make the zero values a bit closer to the non-zero values 
# when the axis is log transformed
Dataset_02.4_m4 <- Dataset_02.4_m3 %>% 
  mutate(ZOONOTIC_NB_all_2 = ifelse(ZOONOTIC_NB_all==0, 0.25, ZOONOTIC_NB_all))

alpha_datapoint = 0.7
size_datapoint = 1.75

set.seed(2)
fig3c <- 
  ggplot() +
  geom_jitter(data=Dataset_02.4_m4, 
              aes(x= timeTraded_all, 
                  y= ZOONOTIC_NB_all_2+1,
                  fill= ZOONOTIC_NB_all+1,
                  # stroke= 1/(10+ZOONOTIC_NB_all), # ifelse(Zoonose_vector_1 == "Zoonotic", 0, 0.1)
                  # color= 1/(10+ZOONOTIC_NB_all)
              ),
              shape=21,
              color="black",
              stroke= 0.1,
              size=size_datapoint,
              alpha=alpha_datapoint,
              width=1, 
              height=0.15) +
  
  geom_ribbon(data=avg_effects_clover_timeTraded_all_scale_df, 
              aes(x=x_unscaled, y=predicted+1,
                  ymax=conf.high+1, ymin=conf.low+1),
              fill=custom_colors[1],
              col="black",
              alpha=0.7) + 
  geom_line(data=avg_effects_clover_timeTraded_all_scale_df, 
            aes(x=x_unscaled, y=predicted+1),
            size=1,
            col= "black") +
  
  annotate("text", x=-Inf, y=Inf, label= "c", size=8, vjust = 1.5, hjust = -1) +
  scale_fill_gradientn(colours=color_map_NbZoonoses,
                       trans = "log",
                       breaks = c(1, 2, 5, 10, 20, 50)
  ) +
  scale_colour_gradientn(colours= (c("white", "grey2"))) +
  scale_y_continuous(trans="log",
                     breaks = c(0.25, 1, 2, 4, 8, 16, 32, 64, 128)+1,
                     labels = c(0, 1, 2, 4, 8, 16, 32, 64, 128),
                     name='Number of pathgens\nshared with humans') +
  
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)*40,
                     labels = c(0, 0.25, 0.5, 0.75, 1)*40,
                     limits=c(-0.015, 1.025)*40,
                     name='Time in trade (years)') +
  theme_bw() +
  theme(legend.position= "none",
        axis.text.x = element_text(size=10), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )

# fig3c



fig3d <- ggplot() +
  PupillometryR::geom_flat_violin(data=Dataset_02.4_m4, aes(y= ZOONOTIC_NB_all_2+1, 
                                                            x=Traded_illegal, 
                                                            fill=Traded_illegal,
                                                            col=Traded_illegal),
                                  position = position_nudge(x = 0.3, y = 0),
                                  adjust =1.5,
                                  size=0.7,
                                  col="black",
                                  alpha = 0.75) +
  geom_point(data=Dataset_02.4_m4, aes(y= ZOONOTIC_NB_all_2+1, 
                                       x=Traded_illegal, 
                                       fill=Traded_illegal),
             alpha=alpha_datapoint-0.2,
             col="black",
             stroke=0.1,
             shape = 21,
             size = size_datapoint,
             position = position_jitter(seed = 1, width = 0.15, height = 0.1)) + 
  
  scale_colour_gradientn(colours= paletteer::paletteer_c("ggthemes::Classic Gray", 9, direction = -1)) +
  ggnewscale::new_scale_colour() +
  geom_linerange(data=avg_effects_clover_Trade_illegal_df,
                 aes(x=Trade_illegal, 
                     y=predicted+1, 
                     ymin=conf.low+1, ymax=conf.high+1, 
                     group=group),
                 position = position_dodge(.25), 
                 size=1.5) +
  geom_point(data=avg_effects_clover_Trade_illegal_df,
             aes(x=Trade_illegal, y=predicted+1,
                 fill= Trade_illegal,
                 group=group),
             position=position_dodge(width=0.25),
             shape=22, 
             size=2.5,
             stroke=1) +
  annotate("text", x=-Inf, y=Inf, label= "d", size=8, vjust = 1.5, hjust = -1) +
  #coord_flip() +
  
  scale_fill_manual(values= (c("grey90", custom_colors[["Traded_illegalyes"]]))) +
  scale_colour_manual(values= (c("grey90", custom_colors[["Traded_illegalyes"]]))) +
  
  scale_x_discrete(name='Occurs in the illegal trade') +
  scale_y_continuous(trans="log",
                     breaks = c(0.25, 1, 2, 4, 8, 16, 32, 64, 128)+1,
                     labels = c(0, 1, 2, 4, 8, 16, 32, 64, 128),
                     name='Number of pathogens \nshared with humans') +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()
  )

# fig3d

fig3e <- ggplot() +
  PupillometryR::geom_flat_violin(data= Dataset_02.4_m4, aes(y= ZOONOTIC_NB_all_2+1, 
                                                             x=Traded_all_live, 
                                                             fill=Traded_all_live,
                                                             col=Traded_all_live),
                                  size=0.7,
                                  position = position_nudge(x = 0.3, y = 0),
                                  adjust =1.5,
                                  col="black",
                                  alpha = 0.75) +
  geom_point(data= Dataset_02.4_m4, aes(y= ZOONOTIC_NB_all_2+1, 
                                        x=Traded_all_live, 
                                        fill=Traded_all_live),
             alpha=alpha_datapoint-0.2,
             col="black",
             stroke=0.1,
             shape = 21,
             size = size_datapoint,
             position = position_jitter(seed = 1, width = 0.15, height = 0.1)) + 
  geom_linerange(data=avg_effects_clover_Traded_live_df,
                 aes(x=Traded_live, 
                     y=predicted+1, 
                     ymin=conf.low+1, 
                     ymax=conf.high+1,
                     group=group),
                 position = position_dodge(.25), 
                 size=1.5) +
  geom_point(data=avg_effects_clover_Traded_live_df,
             aes(x=Traded_live, 
                 y=predicted+1,
                 fill= Traded_live,
                 group=group),
             position=position_dodge(width=0.25),
             shape=22, 
             size=2.5,
             stroke=1) +
  scale_fill_manual(values=c("grey90", custom_colors[["Traded_all_liveyes"]])) +
  scale_colour_manual(values=c("grey90", custom_colors[["Traded_all_liveyes"]])) +
  
  annotate("text", x=-Inf, y=Inf, label= "e", size=8, vjust = 1.5, hjust = -1) +
  #coord_flip() +
  
  scale_x_discrete(name='Occurs in the live market') +
  scale_y_continuous(trans="log",
                     breaks = c(0.25, 1, 2, 4, 8, 16, 32, 64, 128)+1,
                     labels = c(0, 1, 2, 4, 8, 16, 32, 64, 128),
                     name='Number of pathogens \nshared with humans') +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()
  )

# fig3e


### Assembling figure 3 - option n°1
blankP <- ggplot() + theme_void()

fig3a_grob <- grid::grid.grabExpr(draw(fig3a))
fig3bc <- ggpubr::ggarrange(blankP, fig3b, blankP, fig3c, ncol=4, widths=c(0.1, 0.3, 0.05, 0.35))
fig3bc_2 <- ggpubr::ggarrange(fig3bc, blankP, nrow=2, heights=c(0.9,0.1))
fig3de <- ggpubr::ggarrange(blankP, fig3d, fig3e, ncol=3, widths=c(0.1, 0.35, 0.35))
fig3de_2 <- ggpubr::ggarrange(fig3de, blankP, nrow=2, heights=c(0.9,0.1))
fig3bcde <- ggpubr::ggarrange(fig3bc_2, fig3de_2, blankP, nrow=3, heights=c(0.45, 0.45, 0.1))

fig3 <- grid.arrange(fig3a_grob, 
                     #nullGrob(),
                     fig3bcde, 
                     ncol = 2, 
                     widths = c(0.5, 0.8)
)


system.time(ggsave("./outputs/fig3_vf88.pdf", fig3, width = 20, height = 20, units = "cm", limitsize = FALSE))



#### Derivative to calculate the numbers of years needed to share one more pathogen with humans
# Extract predictor and effect values
x_clover <- avg_effects_clover_timeTraded_all_scale_df$x_unscaled
y_clover <- avg_effects_clover_timeTraded_all_scale_df$predicted
ymin_clover <- avg_effects_clover_timeTraded_all_scale_df$conf.low
ymax_clover <- avg_effects_clover_timeTraded_all_scale_df$conf.high

# Compute derivative using finite differences
dy_dx_clover <- diff(y_clover) / diff(x_clover)
dymin_dx_clover <- diff(ymin_clover) / diff(x_clover)
dymax_dx_clover <- diff(ymax_clover) / diff(x_clover)

# Align derivative results with predictor values (midpoints of x_unscaled-intervals)
x_mid_clover <- head(x_clover, -1) + diff(x_clover) / 2

# Create a data frame with the results
derivative_data_clover <- data.frame(x_mid = x_mid_clover, 
                                     dy_dx = dy_dx_clover, 
                                     dymin_dx = dymin_dx_clover, 
                                     dymax_dx = dymax_dx_clover)

derivative_data_clover$rate_y= 1/derivative_data_clover$dy_dx
derivative_data_clover$rate_ymin= 1/derivative_data_clover$dymin_dx
derivative_data_clover$rate_ymax= 1/derivative_data_clover$dymax_dx

derivative_data_clover$rate_y %>% summary
derivative_data_clover$rate_ymin %>% summary
derivative_data_clover$rate_ymax %>% summary


# Plot the derivative 
figST1_derivative <-
  ggplot(data=derivative_data_clover, aes(x=x_mid)) +
  geom_hline(yintercept=mean(derivative_data_clover$rate_y), # average of mean = 17.41
             linetype="solid", color = "indianred", size=1.25) +
  geom_hline(yintercept=mean(derivative_data_clover$rate_ymin), # average of upper CI bound = 20.43
             linetype="dashed", color = "indianred", size=0.75) +
  geom_hline(yintercept=mean(derivative_data_clover$rate_ymax), # average of loser CI bound = 15.150
             linetype="dashed", color = "indianred", size=0.75) +
  geom_ribbon(aes(ymax=rate_ymin, ymin=rate_ymax),
              fill="#ac76a47f") +
  geom_line(aes(y=rate_y),
            size=1.5) +
  scale_y_continuous(limits=c(0, 20)) +
  # annotate("text", x=-Inf, y=Inf, label= "c", size=8, vjust = 1.5, hjust = -1) +
  xlab("Time in trade (years)") +
  ylab("Number of years in trade \n to share one more pathogen with humans") +
  ggtitle("") +
  theme_bw()

figST1_derivative


                                  ## THE END ##

