## ---- echo = FALSE, results = 'hide', message = F------------------------
library(Rnets)
#load(file = '..\\data\\NARMS_EC_DATA.rda')


## ------------------------------------------------------------------------

#Define the set of antimicrobials to include in the Rnet
ABX_LIST <- c('AMC', 'AXO', 'TIO', 'CIP', 'TET', 'STR', 'GEN', 'CHL')

#Estimate the Rnet
EC_all_Rnet <- Rnets::Rnet(Data = NARMS_EC_DATA, L1 =  0.25, V_set = ABX_LIST)
                
#View Results
summary(EC_all_Rnet)

## ------------------------------------------------------------------------
EC_all_L1Selection <- L1Selection(
            Data = NARMS_EC_DATA, 
            L1_set = seq(0.05, 0.50, 0.05),
            n_b = 1500,
            V_set = ABX_LIST
            )

round(EC_all_L1Selection@StARS_D, 4)

