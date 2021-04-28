## Eularian Particle Sinking and Remineralization Functions

## Load in dependency functions

source("UVP_2017_library.R")

## Global Settings

lb_vec <- c(0.102, 0.128, 0.161, 0.203, 0.256, 0.323, 0.406, 0.512, 0.645, 
            0.813, 1.02, 1.29, 1.63, 2.05, 2.58, 3.25, 4.1, 5.16, 6.5, 8.19, 
            10.3, 13, 16.4, 20.6, 26)

binsize_vec <- c(0.026, 0.033, 0.042, 0.053, 0.067, 0.083, 0.106, 0.133, 0.168, 
                 0.207, 0.27, 0.34, 0.42, 0.53, 0.67, 0.85, 1.06, 1.34, 1.69, 
                 2.11, 2.7, 3.4, 4.2, 5.4, 6)

little_lb <- lb_vec[1] - (lb_vec[2] - lb_vec[1])/2 # size of the particle that the UVP can't see anymore. Eg, things actually shrink to this size but then they vanish from the UVP's view. Just leting it be like, the difference in size of the smallest two bins smaller than the smallest bin.

# properties of a 1mm particle
m1mm = 3.3 * 10^-6; #%g % Alldgedge 1998 % mass of 1mm particle
w1mm = 2; #% m/day # Alldredge and Gotschalk, methinks % sinking speed of 1mm particle
micron = 1e-6; # unit conversion

# C_f = 10.51 # for now
# alpha = .52
# gamma = .26
# ag = alpha + gamma

## Describe relationship between particle size and sinking
# Flux = C_f * D ^ (alpha + gamma)
# Mass = C_m * D ^ alpha
# Sinking Speed = C_w * D ^ gamma
# Where D is sinking speed

C_f_global <- C_f
alpha_global <- alpha
gamma_global <- gamma
ag_global <- ag

Cm = m1mm
Cw = w1mm

## Vectors of the masses, sinking speed and flux contributions of single particles of different sizes
m_vec =  Cm * lb_vec ^ alpha;
w_vec = Cw * lb_vec ^ gamma;
f_vec = C_f_global * lb_vec ^ alpha;


## Functions

#' Apply remineralization dynamics to a particle spectrum, given a flux attenuation
#'remin_shuffle is the core function that runs particle remineralization. It takes an abundance profile,
#'and a percentage change in flux, and some other parameters,
#'and returns a the changes in the number of particles in each size bin accross the size profile
#'and some additonal outputs
#' @param abun_in An abundance profile. This is a one dimensional vector,
#'  with each value corresponding to particle numbers one smaller.
#' @param DFpct The percentage change in flux. 
#' Value is a fraction of one so a value of 0.9 would correspond to 90% of flux retaioned
#' One gets negative values if this is too large. Recommend setting to
#' < 99% and using remin_smooth_shuffle to handle larger values.
#' @param DeltaZ the depth between this depth and the one before it. 
#' DeltaZ affects `Cr` but cancells itself out in such a way that the main delta_nj values are not affected.
#' @param Cm A mass constant, defaults to 3.3 * 10^-6; #%g % Alldgedge 1998 % mass of 1mm particle
#' @param CW A sinking speed constant defaults to 2 #% m/day # Alldredge and Gotschalk, methinks % sinking speed of 1mm particle
#' @param lbv a vector of lower bounds of particle sizes for each size bin
#' should have same length as abun_in
#' @param mv a vector of the masses of the particles in each bin
#' @param wv a vector of the sinking speeds of the particles in each bin
#' @param alpha the particle mass fractal dimension
#' @param gamma the particle sinking speed fractal dimension
#' @param llb the size that particles become when they are no longer considered by the model
#' @return A list with different parameter outputs
#' \itimize{
#' list(
#' Cr - Effective partical remineralization rate (1/day)
#' \item phi - a placeholder vector, used to calculate changes in particle numbers
#' \item Delta_nj_net - the net change in particle numbers in each size bin (Delta_nj_in - Delta_nj_out)
#' \item Delta_nj_in - the numbers of particles lost from each size bin due to remineralization
#' \item Delta_nj_out - the numbers of particles added to each size bin due to remineralization
#' )
#' }
remin_shuffle <- function(abun_in, DFpct, DeltaZ = 10, Cm = m1mm, Cw = w1mm, lbv = lb_vec, mv = m_vec, wv = w_vec,
                          alpha = alpha_global, gamma = gamma_global, llb = little_lb){
 ail = length(abun_in)
 rn = abun_in * lbv
 ran = abun_in * lbv ^ alpha
 srn = sum(rn)
 sran = sum(ran)
 
 # For C_r only, didn't help
 ran2 = ran[2:ail]
 sran2 = sum(ran2)
 
 omega = lbv[1] ^ (2 * alpha)  * abun_in[1]/(lbv[1] ^ alpha - llb ^ alpha) # flipped the denominater so its positive

 #omega = lbv[1] ^ (2 * alpha)  * abun_in[1]/(llb ^ alpha - lbv[1] ^ alpha) ## Omega should calculate like this but the flux attenuation is closer to expected if I use the equation above
 
 nmw = abun_in * mv * wv
 F1 = sum(nmw)
 DeltaF = (F1 * DFpct) - F1 # should be negative
 
 #Cr = DeltaF/ (Cm*(1+gamma/alpha) * DeltaZ * (sran + omega));
 #Cr = DeltaF/ (Cm*  DeltaZ * ((1+gamma/alpha) * sran + omega)); # Possible correction (latest)
 Cr = DeltaF/ (Cm*  DeltaZ * (1+gamma/alpha) * (sran + omega)); #J burchfield modification
 #Cr = DeltaF/ (Cm*  DeltaZ * ((1+gamma/alpha) * sran2 + omega)); # Avoid double dip (worse fit)
 #Cr = DeltaF/ (Cm*(1+gamma/alpha) * DeltaZ * (sran)); # Test
 CrCw = Cr/Cw
 
 phi = CrCw * (1+gamma/alpha) * ran * 
   DeltaZ/
   (c(llb, lbv)[1:length(lbv)]^alpha - lbv ^ alpha) # added extra parentheses
 
 Delta_nj_out = phi/lbv ^ gamma
 Delta_nj_in = c(phi[2:length(phi)],0)/c(lbv[2:length(phi)],1) ^ gamma
 #Delta_nj_in = -c(phi[2:length(phi)],0)/lbv ^ gamma
 Delta_nj_net = Delta_nj_in - Delta_nj_out # positive because out is negative
   
  return(list(Cr = Cr, phi = phi, dnet = Delta_nj_net, din =  Delta_nj_in, dout = Delta_nj_out))
}

#' runs remin shuffle and uses it to provide an expected number of particles in each bin, rather than
#' changes in particle numbers in each bin and other parameters
#' inherits all input variables from remin_shuffle
remin_shuffle_spec <- function(abun_in, DFpct, lbv = lb_vec, mv = m_vec, wv = w_vec, llb = little_lb, ...){
  core <- remin_shuffle(abun_in, DFpct, lbv = lbv, mv = mv, wv = wv, llb = llb, ...)
  abun_in + core$dnet
}

#' Iterate over remin_shuffle so that changes behave in a compound way,
#' rather than extrapolating changes from the particle size distribution as it was
#' prevents negative abunadance values
#' All values inherit from remin_shuffle except
#' @param Ipct the flux loss from each step of the function
#' @return abun_est vector of particle abundances in each size class
remin_smooth_shuffle <- function(abun_in, DFpct, Ipct = 0.9999, lbv = lb_vec, mv = m_vec, wv = w_vec, llb = little_lb, ...){
  # DFpct: Fractional mass retained between depths
  # Ipct: Fractional mass retained between iterations
  # ...: Passed to remin_shuffle
  IMirror <- 2  - Ipct
  
  abun_est = abun_in
  Fpct = DFpct # gets overwritten if we are iterating
  
  # If we are loossing flux, and we loose loose more flux than Ipct
  # Iterate remin shuffle only keeping ipct each time
  if(DFpct < Ipct){
    iters = floor(log(DFpct)/log(Ipct)) # why is this a ratio of log transformed values? 
    # I need to understand this before I can address the case where DFpct > IMirror
    iterFlux = Ipct^iters
    Fpct = 1-(iterFlux-DFpct)
    
    for (i in 1:iters){
      abun_est = remin_shuffle_spec(abun_in = abun_est, DFpct = Ipct, lbv = lbv, llb = llb, mv = mv, wv = wv, ...)
    }
  }
  
  # If we are gaining flux, and we gain more than IMirror, iterate
  if(DFpct > IMirror){
    iters = floor(log(DFpct)/log(IMirror))
    iterFlux = IMirror^iters
    Fpct = 1 - (iterFlux - DFpct)
    
    for (i in 1:iters){
      abun_est = remin_shuffle_spec(abun_in = abun_est, DFpct = IMirror, lbv = lbv, llb = llb, mv = mv, wv = wv, ...)
    }
  }
  
  # Deal with remainder. In the case where the loss is less than ipct, or greater than 2-ipct (Imirror), just do this part
  abun_est = remin_shuffle_spec(abun_in = abun_est, DFpct = Fpct, lbv = lbv, llb = llb, mv = mv, wv = wv, ...)
  
  abun_est
}


#' Calculate the rmse of the difference between the flux loss from remin_smooth_shuffle and the
#'  intended flux loss
#'  Inherits all parametrs from remin_smooth shuffle plus
#'  @param DFpct_toRemin Percentage of flux retained, sent to remin_smooth shuffle,
#'   varied by optimization functions
#'   @param DFpct_target The amount of flux actually lost
shuffle_tune <- function(DFpct_toRemin, abun_in,  DFpct_target, llb = little_lb, mv = m_vec, wv = w_vec, lbv = lb_vec,...){
  abun_out <- remin_smooth_shuffle(abun_in = abun_in, DFpct = DFpct_toRemin, llb = llb, mv = mv, wv = wv, lbv = lbv, ...)
  flux_in <- sum(abun_in * C_f_global * lbv ^ ag_global)
  flux_out <- sum(abun_out * C_f_global * lbv ^ ag_global)
  DFPct_actually_happened <- flux_out/flux_in
  rmse <- (DFPct_actually_happened - DFpct_target)^2
  rmse * 100
}

#' Determine the DFP value that actually provides the preferred flux loss. 
#' @return The DFP value which if passed to remin_smooth_shuffle, will
#'  actually warrant the DFpct value passed to this function
optFun <- function(abun_in, DFpct, lbv = lb_vec, llb = little_lb, mv = m_vec, wv = w_vec, ...){
  opt <- optimize(shuffle_tune, c(0, 2), abun_in = abun_in, DFpct_target = DFpct, llb = llb, mv = mv, wv = wv, lbv = lbv, ...)
  opt$minimum
}