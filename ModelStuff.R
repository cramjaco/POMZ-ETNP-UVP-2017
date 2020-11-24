## Global Settings

C_f = 10.51 # for now
alpha = .52
gamma = .26
ag = alpha + gamma

C_f_global <- C_f
alpha_global <- alpha
gamma_global <- gamma
ag_global <- ag


## Disagg functions

remin_shuffle <- function(abun_in, DFpct, DeltaZ = 10, size = lb_vec, Cm = m1mm, Cw = w1mm, lbv = lb_vec, mv = m_vec, wv = w_vec,
                          alpha = 0.52, gamma = 0.26, llb = little_lb){
 rn = abun_in * lb_vec
 ran = abun_in * lb_vec ^ alpha
 srn = sum(rn)
 sran = sum(ran)
 
 omega = lb_vec[1] ^ (2 * alpha)  * abun_in[1]/(lb_vec[1] ^ alpha - llb ^ alpha)
 #omega = lb_vec[1] ^ (2 * alpha)  * abun_in[1]/(llb ^ alpha - lb_vec[1] ^ alpha) # possible correction
 
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

remin_shuffle_spec <- function(abun_in, ...){
  core <- remin_shuffle(abun_in, ...)
  abun_in + core$dnet
}

remin_smooth_shuffle <- function(abun_in, DFpct, Ipct = 0.9999, ...){
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
      abun_est = remin_shuffle_spec(abun_in = abun_est, Ipct, ...)
    }
  }
  
  # If we are gaining flux, and we gain more than IMirror, iterate
  if(DFpct > IMirror){
    iters = floor(log(DFpct)/log(IMirror))
    iterFlux = IMirror^iters
    Fpct = 1 - (iterFlux - DFpct)
    
    for (i in 1:iters){
      abun_est = remin_shuffle_spec(abun_in = abun_est, IMirror, ...)
    }
  }
  
  # Deal with remainder. In the case where the loss is less than ipct, or greater than 2-ipct (Imirror), just do this part
  abun_est = remin_shuffle_spec(abun_in = abun_est, Fpct, ...)
  
  abun_est
}

shuffle_tune <- function(DFpct_toRemin, abun_in,  DFpct_target,...){
  abun_out <- remin_smooth_shuffle(abun_in, DFpct_toRemin, ...)
  flux_in <- sum(C_f_global * abun_in ^ ag_global)
  flux_out <- sum(C_f_global * abun_out ^ ag_global)
  DFPct_actually_happened <- flux_out/flux_in
  rmse <- (DFPct_actually_happened - DFpct_target)^2
  rmse
}

optFun <- function(abun_in, DFpct){
  opt <- optimize(shuffle_tune, c(0, 2), abun_in = abun_in, DFpct_target = DFpct)
  opt$minimum
}