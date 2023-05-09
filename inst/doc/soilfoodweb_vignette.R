## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# Load the package
library(soilfoodwebs)

## -----------------------------------------------------------------------------
properties = data.frame(ID = c("Predator", "Prey1", "Prey2", "Microbe", "Detritus"), # Name
                        d = c(1,3,0.5,1.2,0), # Death rate
                    a = c(0.61,0.65,0.45,1,1), # Assimilation efficiency
                    p = c(0.5,0.4,0.3,0.3,1), # Production efficiency for carbon (nitrogen production efficiency is assumed to be 1)
                    B = c(0.1,8,5,9000,1380000), # Biomass
                    CN = c(4.5,4.8,5,8,30), # Carbon to nitrogen ratio
                    DetritusRecycling = c(0,0,0,0,1), # proportion of detritus recycling
                    isDetritus = c(0,0,0,0,1), # Boolean: Is this pool detritus?
                    isPlant = c(0,0,0,0,0), # Boolean: Is this pool a plant?
                    canIMM = c(0,0,0,1,0)) # Boolean: Can the pool immobilize inorganic nitrogen?

## -----------------------------------------------------------------------------
# Create a list of feeding relationships:
feedinglist <- data.frame(Predator = c("Predator", "Predator","Prey1","Prey1", "Prey2", "Microbe"),
           Prey = c("Prey1", "Prey2", "Microbe", "Detritus", "Detritus", "Detritus"),
           Preference = c(1,1,1,1,1,1))

## -----------------------------------------------------------------------------
# Make the community

our_comm <- build_foodweb(feeding = feedinglist,
                          properties = properties)

# Clean up our working files
rm(feedinglist)

## -----------------------------------------------------------------------------
sp_names = c("Predator", "Prey1", "Prey2", "Microbe", "Detritus")
feedingmatrix = matrix(c(0,1,1,0,0,
                         0,0,0,1,1,
                         0,0,0,0,1,
                         0,0,0,0,1,
                         0,0,0,0,0), 
                       dimnames = list(sp_names, sp_names),
                       ncol = 5, nrow = 5, byrow = TRUE)

## -----------------------------------------------------------------------------
# Make the community
our_comm <- list(imat = feedingmatrix,
                 prop = properties)
# Clean up our working files
rm(feedingmatrix, sp_names, properties)

## -----------------------------------------------------------------------------
# Analysis
ana1 <- comana(our_comm)

#Returns the following results:

# 1. Consumption rates for each species or input rates into the system for any species/pool that has no prey
ana1$consumption

# 2. Carbon mineralized by each species
ana1$Cmin

# 3. Nitrogen mineralization by each species
ana1$Nmin

# 3.1. Contribution of each prey species to nitrogen mineralized when feeding on each of their prey
ana1$Nminmat

# 4. Consumption rates for each species on each of their prey items
ana1$fmat

# 5. The community returned back
ana1$usin


## ---- fig.dim = c(6,6)--------------------------------------------------------
# Produce a plot if you want
old.par <- par(no.readonly = TRUE)
par(mfrow = c(2,2), mar = c(4,4,2,2))
ana1 <- comana(our_comm, mkplot = T, BOX.SIZE = 0.15, BOX.CEX = 0.7, edgepos = c(0.15, 0.85))
par(old.par)

## ---- warning=F---------------------------------------------------------------
# Draw parameter uncertainty:

# Build function inputs to indicate distribution, error measure, and error type using an empty matrix.
guidemat <- matrix(NA, nrow = length(our_comm$prop$ID), ncol = 2, dimnames = list(our_comm$prop$ID, c("B","a")))

# Replicate this matrix across the inputs.
distributionin = guidemat
errormeasurein = guidemat
errortypein = guidemat

# Test two most useful distributions (normal often produces negative values).
distributionin[,1] = "gamma"
distributionin[,2] = "uniform"

# Gamma error uses coefficient of variation. Uniform distribution MUST be the minimum (maximum is the parameter value in the input community).
errortypein[,1] = "CV"
errortypein[,2] = "Min"

# Set these values.
errormeasurein[,1] = 0.2
errormeasurein[,2] = 0

# Run the calculation.
oot <- parameter_uncertainty(usin = our_comm, distribution = distributionin, errormeasure = errormeasurein, errortype = errortypein, fcntorun = "whomineralizes", parameters = c("B", "a"), replicates = 10)

# Bind together the result
oot <- do.call("rbind", oot)

# Summarize the direct and indirect effects of Predators across 10 parameter draws:
summary(subset(oot, ID == "Predator")[,-1])
# NOTE: You can subset to either include or exclude any trophic species from the result by filtering on the ID column.

## -----------------------------------------------------------------------------
our_comm2 <- corrstoich(our_comm)
our_comm2

## -----------------------------------------------------------------------------
# Print the new nitrogen mineralization rates.
comana(our_comm2)$Nmin

## -----------------------------------------------------------------------------
# Set diet limits
DL <- our_comm$imat 
# Note: This only works because all feeding is currently set as 1--the same as 100% of the diet can be from that resource if necessary.
DL["Prey1", "Microbe"] = 0.2 # Microbes only 20% of the diet.
corrstoich(our_comm, dietlimits = DL)


## -----------------------------------------------------------------------------
# Combine the most similar trophic species
comtrosp(our_comm)

# Combine the Predator and Prey1
comtrosp(our_comm, selected = c("Predator", "Prey1"))

## -----------------------------------------------------------------------------
# Delete cannibalism and reset feeding preferences
comtrosp(our_comm, selected = c("Predator", "Prey1"), deleteCOMBOcannibal = T, allFEEDING1 = T)

## -----------------------------------------------------------------------------
# Delete trophic species
removenodes(our_comm, "Predator")

## -----------------------------------------------------------------------------
# Add a new trophic species
newnode(our_comm, "NewNode", prey = c(Detritus = 1), predator = c(Predator = 2, Prey1 = 0.1), newprops = c(d = 1, a = 0.1, p = 0.1, B = 10, CN = 10, DetritusRecycling = 0, isDetritus = 0, isPlant = 0, canIMM = 0))

## ---- fig.dim = c(6, 6)-------------------------------------------------------
# Create function to modify fungal production efficiency
ccxf <- function(p, ccx){
  
  ccx$prop$p[4] = p[1]
  ccx$prop$p[5] = p[2]
  
  return(c(Cmin = sum(comana(corrstoich(ccx))$Cmin), Nmin = sum(comana(corrstoich(ccx))$Nmin)))
}

# Return carbon and nitrogen mineralization across these gradients
res1 = expand.grid(1:10/10,1:10/10)
res2 = res1

for(i in 1:100){
  res2[i,] = ccxf(as.numeric(res1[i,]), ccx = intro_comm)
}

res = cbind(res1, res2)
colnames(res) = c("p1", "p2", "Cmin", "Nmin")

# Create color gradient to work with
res$p1col = palette.colors(10, "Polychrome 36")[res$p1*10]
res$p2col = palette.colors(10, "Polychrome 36")[res$p2*10]

# Plot the results
old.par <- par(no.readonly = TRUE)
par(mfrow = c(2,2))
plot(Cmin~p1, data = res, col = p2col, xlab = "Prod. eff. Fungi 1", pch = 19)
abline(h = sum(comana(corrstoich(intro_comm))$Cmin), lty = 2)

plot(Nmin~p1, data = res, col = p2col, xlab = "Prod. eff. Fungi 1", pch = 19)
abline(h = sum(comana(corrstoich(intro_comm))$Nmin), lty = 2)

plot(Nmin~Cmin, data = res, col = p2col, bg = p1col, pch = 21)
abline(h = sum(comana(corrstoich(intro_comm))$Nmin), lty = 2)
abline(v = sum(comana(corrstoich(intro_comm))$Cmin), lty = 2)
plot.new()
legend("topright", legend = seq(0.1, 1, by = 0.1), 
       col = palette.colors(10, "Polychrome 36"),
       pch = 19, cex = 0.5, title = "Production efficiency scale", ncol = 2, xpd = T)
par(old.par)

## ---- fig.dim = c(6, 6)-------------------------------------------------------

temp_moist_microbe <- function(temp, moist, comm){
  
  # Set functional range:
  if(temp < 0 | moist < 0.05) stop("Function not coded for these conditions.")
  
  # Functions and parameters from Xu et al. 2014
  fmm = ifelse(moist < 0.6,log10(0.05/moist)/log10(0.05/moist),1)
  
  fmt = 2.5^((temp - 25)/10)

  # New microbial death rate
  comm$prop[comm$prop$ID == "Microbe", "d"] = 
    comm$prop[comm$prop$ID == "Microbe", "d"]*fmm*fmt
  
  # New microbial assimilation efficiency
  comm$prop[comm$prop$ID == "Microbe", "a"] = 
    (0.43 - (-0.012*(temp-15)))*
    (comm$prop[comm$prop$ID == "Microbe", "CN"]/comm$prop[comm$prop$ID == "Detritus", "CN"])^0.6
  
  # Ignoring changes in maintenance respiration here. See the cited paper for the details.

  return(comm)
}

# Explore the carbon mineralization rate across temperature as an example

res <- data.frame(temp = seq(1,30, length = 10))

res$Cmin = NA

for(i in 1:10){
  res$Cmin[i] = 
    sum( # Sum all C min rates
      comana( # Calculate Cmin
        corrstoich( # Correct stoichiometry
      temp_moist_microbe(res[i, "temp"], 0.7, our_comm)
      )
      )$Cmin
      )
}

plot(Cmin~temp, data = res, xlab = "Temperature", ylab = "C. mineralization", type = "l", lwd = 2)


## -----------------------------------------------------------------------------
# Calculate the mineralization contributions
whomineralizes(our_comm)

## -----------------------------------------------------------------------------
inout <- calculate_inputs(our_comm)
inout

## -----------------------------------------------------------------------------
# Add the new node
our_comm_necro = newnode(our_comm, "Necromass", prey = NA, predator = c(Prey1 = 1, Prey2 = 1, Microbe = 1), newprops = c(d = 0, a = 1, p = 1, B = 100, CN = 5, DetritusRecycling = 1, isDetritus = 1, isPlant = 0, canIMM = 0))
# Modify the DetritusRecycling of the old detritus pool
our_comm_necro$prop$DetritusRecycling = c(0,0,0,0,0,1)

# Check out the community:
our_comm_necro

# New inputs
inout2 = calculate_inputs(our_comm_necro)

## -----------------------------------------------------------------------------
stab1 <- stability(our_comm)

# The Jacobian:
stab1$J

# The eigen decomposition:
stab1$eigen

# The largest eignvalue:
stab1$rmax # The community is stable if this is negative.

# Use the function calc_smin to recover smin:
calc_smin(our_comm)

# Add density-dependence using stability2:
stability2(our_comm, densitydependence = c(1,1,1,0,0))

# stability and stability2 should produce similar results when used in default:
stability(our_comm)$rmax
stability2(our_comm)$rmax

## ---- fig.dim = c(6, 6)-------------------------------------------------------
# Simulate equilibirum over 100 time steps
sim1 <- CNsim(our_comm)
# Modify predator biomass to 80% of equilibrium value
sim2 <- CNsim(our_comm, start_mod = c(0.8,1,1,1,1))

# Modify predator biomass to 80% of equilibrium value with density-dependence.
sim3 <- CNsim(our_comm, start_mod = c(0.8,1,1,1,1), densitydependence = c(1,0,0,0,0))

# Plot the results:
old.par <- par(no.readonly = TRUE)
par(mfrow= c(2,2), mar = c(5,4,2,2))
plot(Predator_Carbon~Day, data = sim2, type = "l", col = "blue")
points(Predator_Carbon~Day, data = sim3, type = "l", col = "orange")
points(Predator_Carbon~Day, data = sim1, type = "l")

plot(Prey1_Carbon~Day, data = sim2, type = "l", col = "blue")
points(Prey1_Carbon~Day, data = sim3, type = "l", col = "orange")
points(Prey1_Carbon~Day, data = sim1, type = "l")

plot(Prey2_Carbon~Day, data = sim2, type = "l", col = "blue")
points(Prey2_Carbon~Day, data = sim3, type = "l", col = "orange")
points(Prey2_Carbon~Day, data = sim1, type = "l")

plot(Microbe_Carbon~Day, data = sim2, type = "l", col = "blue")
points(Microbe_Carbon~Day, data = sim3, type = "l", col = "orange")
points(Microbe_Carbon~Day, data = sim1, type = "l")
legend("bottomright", legend = c("Base", "80% predator", "w/ density-dependence"), col = c("black", "blue", "orange"), lty = 1)
par(old.par)

## ---- fig.dim = c(6, 4)-------------------------------------------------------
# The function decompexpt 
decompres = decompexpt(our_comm, overtime = 10)

# The decomposition constants for each detritus pool identified in the community
decompres$basedecomp

# A table of the direct and indirect effects 
decompres$decompeffects

# Plot the decomposition with and without Microbe
decompdata = decompres$overtime$Detritus
plot(Original~Day, data = decompdata, type = "l", ylim = c(0.5,1))
points(Microbe~Day, data = decompdata, col = "red", type = "l")
legend("bottomleft", legend = c("Original", "No Microbe"), col = c("black", "red"), lty = 1)


## ---- fig.dim = c(6, 4)-------------------------------------------------------
# Simulate detritus experiment at equilibirum over 100 time steps.
# NOTE: This is equivalent to the overtime option presented for decompexpt for the original community although numerical errors in the solver makes it diverge slightly.
sim1 <- CNsim(our_comm, DETEXPT = "Detritus", TIMES = 1:50)

# Modify microbial biomass to 80% of equilibrium value
sim2 <- CNsim(our_comm, start_mod = c(1,1,1,0.5,1), DETEXPT = "Detritus", TIMES = 1:50)

# Modify microbial biomass to 80% of equilibrium value with density-dependence.
sim3 <- CNsim(our_comm, DETEXPT = "Detritus", start_mod = c(1,1,1,0.5,1), densitydependence = c(0,0,0,1,0), TIMES = 1:50)

# Plot the results:
plot(DetExpt~Day, data = sim2, type = "l", col = "blue")
points(DetExpt~Day, data = sim3, type = "l", col = "orange")
points(DetExpt~Day, data = sim1, type = "l")
legend("topright", legend = c("Base", "80% microbial biomass", "w/ density-dependence"), col = c("black", "blue", "orange"), lty = 1)

