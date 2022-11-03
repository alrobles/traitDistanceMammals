if(!require(tidyverse)) install.packages("tidyverse") 
if(!require(readxl)) install.packages("readxl") 
if(!require(janitor)) install.packages("janitor") 
if(!require(ape)) install.packages("ape") # example imputing
if(!require(Rphylopars)) install.packages("Rphylopars") # example imputing



# We read the raw excel data
traits <- readxl::read_excel("data-raw/MAMMALS_traitsToPublish_MS2b_2022.xlsx", sheet = 1)

# We keep only terrestrial extant mammals
traits <- traits %>% 
  filter(MarineOrNot != 1) %>% 
  filter(`extinct?` == 0)

# We are goint to impute the traits according to mammals phylogeny in Upham 2019

traits_diet <- traits %>% 
  select(tiplabel,
         `Diet-Inv`,
         `Diet-Vend`,
         `Diet-Vect`,
         `Diet-Vfish`,
         `Diet-Vunk`,
         `Diet-Scav`,
         `Diet-Fruit`,
         `Diet-Nect`,
         `Diet-Seed`,
         `Diet-PlantO`,
  ) 

# We clean the variable names to have consistence
traits_diet <- traits_diet %>% janitor::clean_names() 

# We change all diet features as numeric values
traits_diet <- traits_diet %>% 
  mutate_at(2:ncol(traits_diet), as.numeric)

# Because diet values are form 0 to 100 we change rescale from 0 to 1
traits_diet <- traits_diet %>% 
  mutate_at(2:ncol(traits_diet), function(x) x/100) 

# We rename tiplabel according to Rphylopars package
traits_diet <- traits_diet %>% 
  rename(species = tiplabel) 

# We read the consensus tree
mammal_tree <- ape::read.nexus("https://raw.githubusercontent.com/n8upham/MamPhy_v1/master/_DATA/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")

# We prune the tree keeping onlye species in our trait matrix
mammal_tree_prune <- ape::keep.tip(mammal_tree, traits_diet$species)


# We create an iteration variable to map the imputation accross all the features 
iter <- 2:ncol(traits_diet)


# We map an anonymus function that impute the feature value
# according a brownian motion model on the phylogenetic tree
# via purrr package

reconstruction <- purrr::map(iter, function(x){
  p_BM <- phylopars(trait_data = select(traits_diet, c(1, all_of(x))),
                    tree = mammal_tree_prune) 
  
  as_tibble(p_BM$anc_recon, rownames = "species")
})

# We keep only the ancestral reconstruction but is piossible to store 
# the variance and the means.

# for a parallel version of this imputation we use the package furrr

# We collect all the outputs in a data frame
reconstruction_all <- reconstruction %>% reduce(inner_join) 

# We separate the cases with NA to mark as imputation case
traits_diet_NA <- traits_diet %>% 
  filter(is.na(diet_inv)) %>% 
  select(species) %>% 
  mutate(imputation = 1)

# We join with the cases that have information
traits_diet_species <- traits_diet %>%
  select(species) %>% 
  left_join(., traits_diet_NA) %>% 
  replace(is.na(.), 0)

# Finally, we join the species information with the traits reconstruction 
# and round the imputed value

traits_diet_species_imputed <- inner_join(traits_diet_species, reconstruction_all) %>% 
  mutate_at(3:12, function(x) round(x, digits = 1) ) 

# We add species name
traits_diet_species_imputed <- select(traits, sciName, tiplabel) %>% 
  rename(species = tiplabel) %>% 
  inner_join(traits_diet_species_imputed) %>% 
  rename(tiplabel = species)

# We write the final data frame as csv file and as rds file
write_csv(traits_diet_species_imputed, "data-raw/traits_diet_species_BM_imputed.csv")
write_rds(traits_diet_species_imputed, "data/traits_diet_species_BM_imputed.rds")


