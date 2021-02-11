## Lab Rotation I; Ulrich Technau
## Temporal order of point mutations in colorectal cancer
## Analysis start date: 02/01/2021
## 

## sample retrieval from cBIO: 
## working with primary colorectal cancer (adenocarcinomas), which are MSI_Stable; 
## working with 1 sample per patient:

## loading data (coming from cBIO portal):
CRC_sample_cohort = read.csv('~/Documents/ESB_Master/Lab Rotation I/CRC_cohort.tsv', sep = '\t')
Annotated_MAF_cbio = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/msk_impact_facets_annotated.ccf.maf', sep = '\t')


## some basic clinical description of the data
## Primary Tumor Site
sort(table(CRC_sample_cohort$Primary.Tumor.Site))
summary(CRC_sample_cohort$Mutation.Count)
CRC_sample_cohort$Sample.ID[which.max(CRC_sample_cohort$Mutation.Count)]

## exclude P-0035103-T01-IM6
CRC_sample_cohort = CRC_sample_cohort[!CRC_sample_cohort$Sample.ID %in% 'P-0035103-T01-IM6', ]

## look into CRC maf
CRC_maf = Annotated_MAF_cbio[Annotated_MAF_cbio$Tumor_Sample_Barcode %in% CRC_sample_cohort$Sample.ID,, drop = F]

## mutational burden
average.mutation.crc = as.data.frame(table(CRC_maf$Tumor_Sample_Barcode))
average.mutation.density = density(average.mutation.crc$Freq)

par(mgp = c(0.5, 0, 0), mar = c(3, 3, 2, 2))
plot(average.mutation.density$x,
     average.mutation.density$y, 
     xlab = 'Mutational load / sample', 
     ylab = 'Density',
     main = 'Primary Colorectal Adenocarcinoma\n Mutational Load distribution',
     ylim = c(0.0051, max(average.mutation.density$y)), 
     xaxs = 'i', 
     lwd = 2,
     type = 'l',
     axes = F)

axis(side = 1,
     at = c(0, 5, 10, 40),
     labels = c('0', '5', '10', '40'),
     tick = T,
     tck = -0.001)

qu.range = quantile(average.mutation.crc$Freq, probs = c(0.9, 0.99))
poly_range = average.mutation.density$x > qu.range[[1]] 
# poly_range2 = average.mutation.density$x > qu.range[[2]] 

polygon(c(qu.range[[1]], average.mutation.density$x[poly_range], max(average.mutation.crc$Freq)),
        c(0, average.mutation.density$y[poly_range], 0),                  
        col = "grey85")    

# polygon(c(qu.range[[2]], average.mutation.density$x[poly_range2], max(average.mutation.crc$Freq)),
#         c(0, average.mutation.density$y[poly_range2], 0),                  
#         col = "red")    

box(lty = 'solid', lwd = 0.8)


## prepare the data for BT modelling
CRC_BT = data.frame(Alteration = CRC_maf$Hugo_Symbol,
                    SAMPLE_ID = CRC_maf$Tumor_Sample_Barcode,
                    ccf = CRC_maf$ccf_expected_copies)
CRC_BT$is_clonal = ifelse(CRC_BT$ccf >= 0.8, 1, 0)
CRC_BT = CRC_BT[!is.na(CRC_BT$is_clonal), ]

# cbio(unique(CRC_maf$Tumor_Sample_Barcode))
