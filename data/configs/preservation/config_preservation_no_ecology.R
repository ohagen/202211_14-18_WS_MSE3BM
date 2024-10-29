########################
### General settings ###
########################
# set the random seed for the simulation
random_seed = 0001

# set the starting time step or leave NA to use the earliest/highest timestep
start_time = NA

# set the end time step or leave as NA to use the lates/lowest timestep (0)
end_time = NA

# maximum total number of species in the simulation before it is aborted
max_number_of_species = 10000

# maximum number of species within one cell before the simulation is aborted
max_number_of_coexisting_species = 1000

# a list of traits to include with each species
# a "dispersion" trait is implictly added in any case
#trait_names = c("t_min", "a_min", "competition", "dispersion")
trait_names = c("temp",  "dispersal", "foss", "niche_wd") # "arid",

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# lsited with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
environmental_ranges = list() #"temp" = c(-45, 55), "area"=c(6895.094, 196948.4), "arid"=c(1,0.5))

# a place to inspect the internal state of the simulation and collect additional information if desired
end_of_timestep_observer = function(data, vars, config){
  #plot
  par(mfrow=c(1,2))
  plot_richness(data$all_species, data$landscape)
  plot_species_presence(data$all_species[[1]], data$landscape)
  #save
  save_species()
  save_landscape()
}



######################
### Initialization ###
######################
# the intial abundace of a newly colonized cell, both during setup and later when colonizing a cell during the dispersal
initial_abundance = 1
# place species within rectangle, our case entire globe
create_ancestor_species <- function(landscape, config) {
  range <- c(-180, 180, -90, 90)
  co <- landscape$coordinates
  selection <- co[, "x"] >= range[1] &
    co[, "x"] <= range[2] &
    co[, "y"] >= range[3] &
    co[, "y"] <= range[4]
  initial_cells <- rownames(co)[selection]
  new_species <- create_species(initial_cells, config)
  #set local adaptation to max optimal temp equals local temp
  new_species$traits[ , "temp"] <- landscape$environment[,"temp"]
  #set fossilization to be randomly selected in a uniform prior
  new_species$traits[ , "foss"] <- runif(1,0,1)
  new_species$traits[ , "dispersal"] <- 1
  new_species$traits[, "niche_wd"] <- 5
  return(list(new_species))
}


#################
### Dispersal ###
#################
# returns n dispersal values
get_dispersal_values <- function(n, species, landscape, config) {
  values <- rweibull(n, shape =3, scale =781) ### VARY
  return(values)
}


##################
### Speciation ###
##################
# threshold for genetic distance after which a speciation event takes place
divergence_threshold=5

# factor by which the divergence is increased between geographicaly isolated population
# can also be a matrix between the different population clusters
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  return(1)
}


################
### Mutation ###
################
# mutate the traits of a species and return the new traits matrix
apply_evolution <- function(species, cluster_indices, landscape, config) {
  # trait_evolutionary_power <-0.085 ### VARY
  traits <- species[["traits"]]
  # cells <- rownames(traits)
  # #homogenize trait based on abundance
  # for(cluster_index in unique(cluster_indices)){
  #   # cluster_index <- 1
  #   cells_cluster <- cells[which(cluster_indices == cluster_index)]
  #   # hist(traits[cells_cluster, "temp"], main="before")
  #   mean_abd <- mean(species$abundance[cells_cluster])
  #   weight_abd <- species$abundance[cells_cluster]/mean_abd
  #   traits[cells_cluster, "temp"] <- mean(traits[cells_cluster, "temp"]*weight_abd)
  #   # hist(traits[cells_cluster, "temp"], main="after")
  # }
  # #mutations
  # mutation_deltas <-rnorm(length(traits[, "temp"]), mean=0, sd=trait_evolutionary_power)
  # traits[, "temp"] <- traits[, "temp"] + mutation_deltas
  # mutation_deltas <-rnorm(length(traits[, "foss"]), mean=0, sd=trait_evolutionary_power)
  # traits[, "foss"] <- traits[, "foss"] + mutation_deltas
  # # rang between 0 and 1
  # tf <- traits[,"foss"]
  # tf[tf>1] <- 1
  # tf[tf<0] <- 0
  # traits[,"foss"] <- tf
  return(traits)
}


###############
### Ecology ###
###############
# called for every cell with all occuring species, this functin calculates the who survives in the current cells
# returns a vector of abundances
# set the abundance to 0 for every species supposed to die
apply_ecology <- function(abundance, traits, landscape, config, abundance_scale = 10, abundance_threshold = 1) {
  
  # fg <- function(x,a,b,c){
  #   v <- (a/c)*exp(-((x-b)^2/(2*c^2)))
  #   return(v)
  # }
  # #
  # plot(fg(x=0, a=10, b=seq(-50,50,1), c=5), type='l') # c ranges from 0.001 to 0.3 (very wide niche)
  # # abline(h=1)
  
  # gaussian
  
  
  
  # abundance <- ((abundance_scale/traits[,"niche_wd"])*exp(-((traits[, "temp"] - landscape[, "temp"])**2/(traits[,"niche_wd"]**2))))*(landscape[,"arid"])
  # #abundance thhreashold
  # abundance[abundance<abundance_threshold] <- 0
  
  return(abundance)
}
