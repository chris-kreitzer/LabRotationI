## Running the BradleyTerry Model
## 01/12/2020

source('~/Documents/MSKCC/05_IMPACT40K/CancerEvolution/Scripts/BradleyTerryUpdated.R')
x = BradleyTerryUpdate(mat = CRC_BT, game.cutoff = 8, 
                       plotting_refinement = T, 
                       COSMIC_cancer_type_abbreviation = c('colorectal', 'colon cancer', 'CRC'))
