# Defining Model Equations
ode <- "
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

    # Defining sign of Tmax value based on strain
    Tmax_S = -1
    Tmax_S_BL6 = 1
    
    
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
"

# Defining model parameters
parameters <- c(jtb = 4.3075, jpt = 0.32,jpt_NSG = 1, jps = 11.8661,
                jpi = 42.0405, kcl0_B = 0.4395, jbl = 2.4, jlp = 2.4,
                Tmax = 1.3443, T50 = 21.2117, g = 3.2069, Tmax_BL6 = 1.8742, T50_BL6 = 7.6761, 
                g_BL6 = 3.1862, kcl0_S = 0.2644, 
                kcl0_I = 0.2813, kcl0_P = 0.1831, V0_B = 1.1878, V_T = 0.05, V_S = 0.1,
                V_I = 0.1, V_P = 0.05, V_L = 0.05, RAG2 = 1, BCD = 0, NSG = 0,BL6 = 0, 
                WEIGHT = 19,median_weight = 19)
# Defining Initial Condition
y0 <- c(0,0,0,0,0,0)


