#clearing the workspace
rm(list=ls())
graphics.off()
options(show.error.locations = TRUE)

# Calling required libraries
library(nlmixr2)
library(rxode2)
library(readxl)
library(data.table)
library(ggplot2)
library(xpose)
library(xpose.nlmixr2)
library(writexl)

# Loading workspace image of estimation results
filename = "estimation_results"
load(filename)

# Printing summary of estimation results
print(nlme_results)


# Defining model in ODE form for plotting with rxode2 package
source('mPBPK_model_ODE.R')
ode = ode
parameters = parameters
y0 = y0
model = RxODE(model = ode, modName = "model")



# Setting up groups of mice for plotting
IDs = unique(nlme_results[["ID"]])
groups = list((1:6),(7:11),(12:16), (17:22),(23:27),(28:30),(31:39),(40:43),(44:46),(47:53),(54:58),(59:63))
group_cmpts = c("B","I","S","B","S","I","B","S","I","B","I","S")
group_names = c("RAG2_IV","RAG2_IP","RAG2_SC","BL6_IV","BL6_SC",
                "BL6_IP","NSG_IV","NSG_SC","NSG_IP","BCD_IV","BCD_IP",
                "BCD_SC")

# Setting up simulation times for each group
tspan = list(seq(from = 0, to = 100, length.out = 1000),
             seq(from = 0, to = 100, length.out = 1000),
             seq(from = 0, to = 100, length.out = 1000),
             seq(from = 0, to = 30, length.out = 1000),
             seq(from = 0, to = 30, length.out = 1000),
             seq(from = 0, to = 30, length.out = 1000),
             seq(from = 0, to = 15, length.out = 1000),
             seq(from = 0, to = 15, length.out = 1000),
             seq(from = 0, to = 15, length.out = 1000),
             seq(from = 0, to = 100, length.out = 1000),
             seq(from = 0, to = 100, length.out = 1000),
             seq(from = 0, to = 100, length.out = 1000)
)

# Defining estimated parameters
fitted_parameters = c("jtb","jpt","jpt_NSG","jps","jpi",
                      "kcl0_B","Tmax","T50","g","Tmax_BL6","T50_BL6",
                      "g_BL6","kcl0_S","kcl0_I","kcl0_P","V0_B","WEIGHT",
                      "RAG2","BCD","BL6","NSG")

# Plotting predictions v. experimental data
# Iterating through each group, plotting expeirmental data and simulated timecourse with ODE model for each mouse
count = 1
for (G in groups){
  group_indices = G
  fig <- ggplot()
  for (ind in group_indices){
    # Extracting information on individual for simulation & plotting
    ind_results = nlme_results[nlme_results[["ID"]] == ind,]
    ind_DV = ind_results[["DV"]]
    ind_t = ind_results[["TIME"]]
    exp_data = data.frame(DV = ind_DV, time = ind_t)
    # Setting up dosing information for simulation
    ev <- eventTable()
    ev$add.dosing(dose = 5.6e4, nbr.doses = 1,cmt = group_cmpts[count])
    ev$add.sampling(tspan[[count]])
    # Setting model parameters based on individual estiamted values
    for (param in fitted_parameters){
      parameters[[param]] = ind_results[[param]][1]
    }
    # Running ODE model
    ind_pred = model$run(parameters,ev,y0)
    
    # Adding individuals data & predictions to group plot
    fig <- fig + theme_classic() + geom_line(data = ind_pred, aes(x = time, y = B_obs)) +
      geom_hline(yintercept = 2.5, size = 1,linetype = "dashed")+
      geom_point(data = exp_data,aes(x = time,y = DV),color = 'red') + 
      theme(text = element_text(size = 25))
  }
  # Formatting plot once all data is added
  fig <- fig  + xlab('Time (days)') + ylab('Serum [Ab] (ng/mL)') +
    scale_y_log10(labels = function(x) format(x,scientific = TRUE), limits = c(1e-1, 1e5))
  print(fig)
  count = count + 1
}

 
## Plotting EBEs for each strain
# Adding EBE information to NLME results with addCwres()
nlme_results = addCwres(nlme_results)
# Extracting information on individual etas
phi_SE = nlme_results$phiSE
eta_info = nlme_results$shrink
eta_names = colnames(phi_SE)
eta_names = eta_names[2:length(eta_names)]
eta_values = nlme_results$eta
IDs = phi_SE[["ID"]]
# Setting colors for plotting
strain_colors = list(RAG2 = "green3",BCD = "purple4",BL6 = "royalblue3",NSG = "red")


# Iterating through each eta and plotting each individuals EBE
for (eta in eta_names){
  fig <- ggplot() + theme_classic()
  # Adding individual EBEs
  for (ind in IDs){
    eta_ind = eta_values[[eta]][strtoi(ind)]
    eta_ind_SE = phi_SE[[eta]][strtoi(ind)]
    ind_results = nlme_results[nlme_results[["ID"]] == ind,]
    RAG2 = ind_results[["RAG2"]][1]
    BCD = ind_results[["BCD"]][1]
    NSG = ind_results[["NSG"]][1]
    BL6 = ind_results[["BL6"]][1]
    # Checking which strain to plot
    if (RAG2 == 1){
      c = strain_colors[["RAG2"]]
    } else if (BCD == 1){
      c = strain_colors[["BCD"]]
    } else if (NSG == 1){
      c = strain_colors[["NSG"]]
    } else if (BL6 == 1){
      c = strain_colors[["BL6"]]
    }
    
    x <- seq(eta_ind - 4*eta_ind_SE, eta_ind + 4*eta_ind_SE,length = 1000)
    y <- dnorm(x, mean = eta_ind, sd = eta_ind_SE)
    dat = data.frame(x = x, y = y)
    fig <- fig + geom_line(data = dat, aes(x = x, y = y),color = c, size = 0.5)
    fig <- fig + theme(text = element_text(size = 25))
  }
  # Adding labels to plot
  fig <- fig + xlab("Eta") + ylab("Density") + ggtitle(eta)
  print(fig)
}



# using xpose package to create model assessment plots
xp_object <- xpose_data_nlmixr2(nlme_results)
dv_vs_pred(xp_object) + theme_classic() + scale_y_log10() + scale_x_log10()
dv_vs_ipred(xp_object) + theme_classic() + scale_y_log10() + scale_x_log10()

