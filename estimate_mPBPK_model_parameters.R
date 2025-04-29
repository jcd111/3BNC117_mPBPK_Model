#clearing the workspace
rm(list=ls())
graphics.off()
options(show.error.locations = TRUE)


# Calling all required libraries
library(nlmixr2)
library(rxode2)
library(readxl)
library(data.table)
library(ggplot2)
library(xpose)
library(xpose.nlmixr2)
library(writexl)

# Defining Model
source('mPBPK_model.R')
model = mPBPK_model()

# Loading Data
dat <- read_excel("population_PK_data.xlsx")

# Running parameter estimtaion with nlmixr2() function
nlme_results <- nlmixr2(model, dat, est = "focei", control = foceiControl(covMethod = ""))
# Saving results as an excel file and a image of workspace (for plotting in plot_estimation_results)
write_xlsx(nlme_results,path = "estimation_results.xlsx")
save.image(file = "estimation_results")



