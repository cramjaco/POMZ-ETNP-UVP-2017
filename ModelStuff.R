## Global Settings

lb_vec <- c(0.102, 0.128, 0.161, 0.203, 0.256, 0.323, 0.406, 0.512, 0.645, 
            0.813, 1.02, 1.29, 1.63, 2.05, 2.58, 3.25, 4.1, 5.16, 6.5, 8.19, 
            10.3, 13, 16.4, 20.6, 26)

binsize_vec <- c(0.026, 0.033, 0.042, 0.053, 0.067, 0.083, 0.106, 0.133, 0.168, 
                 0.207, 0.27, 0.34, 0.42, 0.53, 0.67, 0.85, 1.06, 1.34, 1.69, 
                 2.11, 2.7, 3.4, 4.2, 5.4, 6)

little_lb <- lb_vec[1] - (lb_vec[2] - lb_vec[1])/2 # size of the particle that the UVP can't see anymore. Eg, things actually shrink to this size but then they vanish from the UVP's view. Just leting it be like, the difference in size of the smallest two bins smaller than the smallest bin.

# mass of a 1mm particle
m1mm = 3.3 * 10^-6; #%g % Alldgedge 1998 % mass of 1mm particle
w1mm = 2; #% m/day # Alldredge and Gotschalk, methinks % sinking speed of 1mm particle
micron = 1e-6;

C_f = 10.51 # for now
alpha = .52
gamma = .26
ag = alpha + gamma

C_f_global <- C_f
alpha_global <- alpha
gamma_global <- gamma
ag_global <- ag

Cm = m1mm
Cw = w1mm
m_vec =  Cm * lb_vec ^ alpha;
w_vec = Cw * lb_vec ^ gamma;
f_vec = C_f_global * lb_vec ^ alpha;


## Disagg functions

remin_shuffle <- function(abun_in, DFpct, DeltaZ = 10, Cm = m1mm, Cw = w1mm, lbv = lb_vec, mv = m_vec, wv = w_vec,
                          alpha = 0.52, gamma = 0.26, llb = little_lb){
 rn = abun_in * lbv
 ran = abun_in * lbv ^ alpha
 srn = sum(rn)
 sran = sum(ran)
 
 omega = lbv[1] ^ (2 * alpha)  * abun_in[1]/(lbv[1] ^ alpha - llb ^ alpha)
 #omega = lbv[1] ^ (2 * alpha)  * abun_in[1]/(llb ^ alpha - lbv[1] ^ alpha) # possible correction
 
 nmw = abun_in * mv * wv
 F1 = sum(nmw)
 DeltaF = (F1 * DFpct) - F1 # should be negative
 
 #Cr = DeltaF/ (Cm*(1+gamma/alpha) * DeltaZ * (sran + omega));
 Cr = DeltaF/ (Cm*  DeltaZ * ((1+gamma/alpha) * sran + omega)); # Possible correction
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

remin_shuffle_spec <- function(abun_in, DFpct, llb = lb_vec, mv = m_vec, wv = w_vec, ...){
  core <- remin_shuffle(abun_in, DFpct, llb = llb, mv = mv, wv = wv, ...)
  abun_in + core$dnet
}

remin_smooth_shuffle <- function(abun_in, DFpct, Ipct = 0.9999, llb = lb_vec, mv = m_vec, wv = w_vec, ...){
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
      abun_est = remin_shuffle_spec(abun_in = abun_est, DFpct = Ipct, llb = llb, mv = mv, wv = wv, ...)
    }
  }
  
  # If we are gaining flux, and we gain more than IMirror, iterate
  if(DFpct > IMirror){
    iters = floor(log(DFpct)/log(IMirror))
    iterFlux = IMirror^iters
    Fpct = 1 - (iterFlux - DFpct)
    
    for (i in 1:iters){
      abun_est = remin_shuffle_spec(abun_in = abun_est, DFpct = IMirror, llb = llb, mv = mv, wv = wv, ...)
    }
  }
  
  # Deal with remainder. In the case where the loss is less than ipct, or greater than 2-ipct (Imirror), just do this part
  abun_est = remin_shuffle_spec(abun_in = abun_est, DFpct = Fpct, llb = llb, mv = mv, wv = wv, ...)
  
  abun_est
}

shuffle_tune <- function(DFpct_toRemin, abun_in,  DFpct_target, llb = lb_vec, mv = m_vec, wv = w_vec,...){
  abun_out <- remin_smooth_shuffle(abun_in = abun_in, DFpct = DFpct_toRemin, llb = llb, mv = mv, wv = wv, ...)
  flux_in <- sum(abun_in * C_f_global * llb ^ ag_global)
  flux_out <- sum(abun_out * C_f_global * llb ^ ag_global)
  DFPct_actually_happened <- flux_out/flux_in
  rmse <- (DFPct_actually_happened - DFpct_target)^2
  rmse
}

optFun <- function(abun_in, DFpct, llb = lb_vec, mv = m_vec, wv = w_vec, ...){
  opt <- optimize(shuffle_tune, c(0, 2), abun_in = abun_in, DFpct_target = DFpct, llb = llb, mv = mv, wv = wv, ...)
  opt$minimum
}