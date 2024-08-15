TestExtinctionPoint_ModelB <- function(seed, params, setup, PLOT=TRUE){
  
  set.seed(seed)
  
  Density_m <- as.data.frame(matrix(NA, ncol = 8, nrow = L)) 
  colnames(Density_m)<- c('Flies', 'Wasps', 'r', 'temp', 'PropFliesSurv', 'w_att_temp_WOTRANS','w_att_tempWithTRANS','WaspRiskFactor')
  Density_m$Flies[1] <- setup$Start_flies
  Density_m$Wasps[1] <- setup$Start_wasps
  
  
  
  if(setup$AutoCorrelated == FALSE){
    Env_vec <- c( rep(setup$Av_temp, 
                      setup$burnin), 
                  seq( from = setup$Av_temp,
                       to = setup$MaxTemp, 
                       length.out = setup$CC_time )) + rnorm(L, 
                                                             mean = 0,
                                                             sd = setup$temp_sd)
  }else{
    
    
    Env_vec <- c( rep(setup$Av_temp, 
                      setup$burnin), 
                  seq( from = setup$Av_temp,
                       to = setup$MaxTemp, 
                       length.out = setup$CC_time )) + Make_AC_data(N = L,
                                                                    AC = setup$autocorrelation,
                                                                    sd = setup$temp_sd)
  }
  
  
  
  Env_vec[ Env_vec>(params$ctmax -0.1) ] <- (params$ctmax-0.1)
  
  ### Fly temperature dependence
  
  for( t in 2: L){
    
    Prev_Flies <- Density_m[t-1,1]
    Prev_Wasps <- Density_m[t-1,2]
    temp = Env_vec[t]
    G0_temp = Env_vec[t-1]
    
    #calculating growth rate of population based on parameters and environmental temp
    #with transgenerational effects and without    
    
    if(setup$transgenerational_fly){
      r <- EPC_trans(G1_temp = temp,
                     G0_temp = G0_temp,
                     cTmax = params$ctmax) 
    }else{
      r <- params$rmax* ifelse(temp < params$topt,
                               exp(-((temp -  params$topt)/
                                       (2 * params$fly_a))^2),
                               1 - ((temp -  params$topt)/
                                      ( params$topt -  params$ctmax))^2)
    }
    
    ### Wasp temperature dependence
    ## Wasp model fitted from data gives us proportion of wasps attacked (ranging from 0-1)
    ## Most direct way to build that is just to scale the 'wasp_a' rate parameter by this
    ## not perfect, but should capture response dynamics
    
    w_att_temp_WITHTRANS = params$wasp_aMODB * predict(waspglm2_poly,
                                                   newdata = data.frame(G0 = G0_temp,
                                                                        G1 = temp) ,
                                                   type='response' )
    
    w_att_tempWITHOUTTRANS = params$wasp_aMODB * predict(waspglm2_noG0,
                                                     newdata = data.frame(G1 = temp) ,
                                                     type='response' ) 
    #with transgenerational effects and without  
    
    if(setup$transgenerational_wasp){
      w_att_temp = w_att_temp_WITHTRANS
    }else{
      w_att_temp = w_att_tempWITHOUTTRANS
    }
    
    # Calculating projected fly population at next time step
    #based on growth rate, previous fly density and carrying capacity
    Next_Flies_potential <- r * (1/params$fly_fail)*Prev_Flies * (1 - (Prev_Flies / params$fly_K )) 
    
    
    # PropFliesSurv = How many flies survive
    
    w_k<- params$wasp_k
    WaspRiskFactor = (1+((w_att_temp*Prev_Wasps)/
                           (w_k*(1 + w_att_temp*params$wasp_h*Next_Flies_potential)) ) )^w_k    
    ## from page 50 (Hassell)
    
    PropFliesSurv = 1/WaspRiskFactor
    
    Next_Flies <- Next_Flies_potential*(PropFliesSurv)
    Next_Wasps <- Next_Flies_potential*(1-PropFliesSurv) +10
    
    Density_m$PropFliesSurv[t] <- PropFliesSurv
    Density_m$r[t] <- r
    Density_m$temp[t] <- temp
    Density_m$Flies[t] <- Next_Flies  
    Density_m$Wasps[t] <- Next_Wasps 
    Density_m$w_att_temp_WOTRANS[t] <- w_att_tempWITHOUTTRANS
    Density_m$w_att_tempWithTRANS[t] <- w_att_temp_WITHTRANS
    Density_m$WaspRiskFactor[t] <- WaspRiskFactor
    
  }
    if(PLOT){
    
    data.frame(Generation = 1:L, 
               Temperature = Env_vec, 
               Flies = Density_m$Flies,
               Wasps = Density_m$Wasps) ->DF
    
    
    EXTINCT<- which.max(Density_m$Flies[201:nrow(Density_m)] < setup$threshold)+200
    
    DF %>%
      ggplot( aes( x = Generation, y = Temperature))+
      geom_line(linewidth = 0.05)+
      theme_bw()+
      geom_hline(yintercept = 30, linetype = 'dashed') -> TempPlot
    
    DF %>%
      ggplot( aes( x = Generation, y = Flies))+
      geom_line(linewidth = 0.05)+
      theme_bw()+
      scale_y_continuous(limits = c( 0, 32000))+
      geom_vline(xintercept = EXTINCT, col = 'red') -> FlyPlot
    
    DF %>%
      ggplot( aes( x = Generation, y = Wasps))+
      geom_line(linewidth = 0.05)+
      scale_y_continuous(limits = c( 0, 32000))+
      theme_bw() -> WaspPlot
    
    plot_grid(TempPlot, FlyPlot,WaspPlot, ncol = 3,
              labels = c('a) Temperature Trajectory', 
                         'b) Fly Dynamics', 
                         'c) Wasp Dynamics'),hjust = -0.1,
              label_fontface = 'plain',
              scale = 0.9) -> XXX
    
    return(XXX)
  }
  
  
#  print(Density_m[1:50,])
  if( all( Density_m$Flies[201:nrow(Density_m)] > setup$threshold)){return(NA)}
  return(which.max(Density_m$Flies[201:nrow(Density_m)] < setup$threshold)+200)
}