mPBPK_model <- function(){
  ini({
    
    
    # Defining initial guesses for log-transformed estimated model parameters
    l_jtb <- log(2.9792)
    l_jpt <- log(0.1523)
    l_jpt_NSG <- log(0.001)
    l_jps <- log(1)
    l_jpi <- log(10)
    l_kcl0_B <- log(0.3213)
    l_Tmax <- log(1)
    l_T50 <- c(log(10), log(20.9980), log(30))
    l_g <- c(log(1), log(3.1925), Inf)
    l_Tmax_BL6 <- log(1)
    l_T50_BL6 <- c(log(1), log(5.9952), log(10))
    l_g_BL6 <- c(log(1), log(3.8853), Inf)
    l_kcl0_S<- log(1)
    l_kcl0_I <- log(10)
    l_kcl0_P <- log(0.0292)
    l_V0_B <- log(0.9621)
    
    
    # Setting inter-individual variability initial guesses
    # All parameters listed below, fixed effect parameter IIVs are commented out
    eta.jtb ~ 0.1
    eta.jpt ~ 0.1
    eta.jps ~ 0.1
    #eta.jpi ~ 0.1
    #eta.kcl0_B ~ .1
    eta.Tmax ~ .1
    eta.T50 ~ .1
    #eta.g ~ .1
    eta.Tmax_BL6 ~ .1
    eta.T50_BL6 ~ .1
    #eta.g_BL6 ~ .1
    eta.kcl0_S ~ .1
    eta.kcl0_I ~ .1
    eta.kcl0_P ~ .1
    #eta.V0_B ~ .1
    
    
    # Setting error model parameters
    prop.err <- c(0,0.3,1)
    add.err <- c(0, 1, Inf)
    
  })
  
  model({
    
    # Setting estimated model parameters based on fixed and random effect values
    jtb <- exp(l_jtb + eta.jtb)
    jpt <- exp(l_jpt + eta.jpt)
    jpt_NSG <- exp(l_jpt_NSG)
    jps <- exp(l_jps + eta.jps)
    jpi <- exp(l_jpi)
    # Setting clearance parameters
    kcl0_B <- exp(l_kcl0_B)
    Tmax <- exp(l_Tmax + eta.Tmax)
    T50 <- exp(l_T50 + eta.T50)
    g <- exp(l_g)
    Tmax_BL6 <- exp(l_Tmax_BL6 + eta.Tmax_BL6)
    T50_BL6 <- exp(l_T50_BL6 + eta.T50_BL6)
    g_BL6 <- exp(l_g_BL6)
    kcl0_S<- exp(l_kcl0_S + eta.kcl0_S)
    kcl0_I <- exp(l_kcl0_I + eta.kcl0_I)
    kcl0_P <- exp(l_kcl0_P + eta.kcl0_P)
    # Compartment volumes
    V0_B <- exp(l_V0_B)
    
    
    # Setting constant model parameters here
    jbl <- 2.4
    jlp <- 2.4
    V_T <- 0.05
    V_S <- .1
    V_I <- .1
    V_P <- 0.05
    V_L <- 0.05
    median_weight <- 19
    # Parameters that define sign of Tmax value for eac hstrain
    Tmax_S <- -1
    Tmax_S_BL6 <- 1
    
    # Calculating all clearance rate constants based on weight.
    kcl_B = kcl0_B*(WEIGHT/median_weight)^-0.15
    kcl_S = kcl0_S*(WEIGHT/median_weight)^-0.15
    kcl_I = kcl0_I*(WEIGHT/median_weight)^-0.15
    kcl_P = kcl0_P*(WEIGHT/median_weight)^-0.15
    
    # Calculating volume of distribution based on weight
    V_B = V0_B*(WEIGHT/median_weight)
    
    # Defining concentrations
    C_B = B/V_B
    C_S = S/V_S
    C_I = I/V_I
    C_T = TI/V_T
    C_P = P/V_P
    C_L = L/V_L

    
    
    # Calculating differential equations
    d/dt(B) = jbl*C_L - jtb*C_B - 
      (RAG2+BCD)*exp((Tmax_S*Tmax*time^g)/(time^g + T50^g))*B*kcl_B - 
      BL6*exp((Tmax_S_BL6*Tmax_BL6*time^g_BL6)/(time^g_BL6 + T50_BL6^g_BL6))*B*kcl_B - 
      NSG*B*kcl_B
    
    d/dt(S) = -jps*C_S - S*kcl_S
    
    d/dt(I) = -jpi*C_I - I*kcl_I
    
    d/dt(TI) = jtb*C_B - kcl_P*TI - 
      (jpt_NSG^NSG)*jpt*C_T
    
    d/dt(P) = jpi*C_I + jps*C_S - jlp*C_P - kcl_P*P + 
      (jpt_NSG^NSG)*jpt*C_T
    
    d/dt(L) = jlp*C_P - jbl*C_L - kcl_P*L
    
    
    B_obs = C_B
    # Defining proportional error 
    B_obs ~ prop(prop.err) + add(add.err)
  })
  
  
  
}