## Running an updated version of the Bradley Terry Model:
## The initial function was written by Bastien Nyguen:
## 
## Here, I included several updated steps (connectivity, game-off parameter, dead-ends, plotting refinements, etc.)
## Have a look into LabMeeting 12/15/2020 to see the updates:
## The documentation follows
## updates made on 12/03/2020 and 01/20/2021
## author: chris kreitzer


## packages required:
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(BradleyTerryScalable))
suppressPackageStartupMessages(library(Matrix))

# data-in
# 'mat: 4 column data frame with ('Alteration' (anything), 'SAMPLE_ID', 'ccf', 'is_clonal')
# 'ccf' either for mutations from IMPACT maf; or of SCNA via cf.em / purity
# game.cutoff: how many independent games should have been played (in independent patients); default = 2
# plotting_refinement: COSMIC and Bailey (cancer gene census) genes are loaded and only those genes in the
# BT model output which have any involvement in specific cancer type will be considered.
# For the plotting_refinement: COSMIC and BAILEY cancer type abbreviations need to be delivered (search within database)
# e.g., KIRC and KIRP for clear cell Renal carcinoma, etc.
# COSMIC and Bailey datasets are stored locally on my computer, hence need to be loaded

BradleyTerryUpdate = function(mat, 
                              game.cutoff = 2,
                              plotting_refinement = NULL,
                              COSMIC_cancer_type_abbreviation = NULL,
                              Bailey_cancer_type_abbreviation = NULL){
  
  print('specific cancer type abb. if you use plotting refinement. Needed for gene selection in COSMIC')
  
  if(!all(c('SAMPLE_ID', 'Alteration', 'ccf', 'is_clonal') %in% colnames(mat))){
    stop('The input matrix needs four columns; SAMPLE_ID, Alteration, ccf, is_clonal')
    
  }
  
  if(ncol(mat) > 4){
    mat = mat[, colnames(mat) %in% c('SAMPLE_ID', 'Alteration', 'ccf', 'is_clonal')]
  }
  
  # load the required package 
  require('BradleyTerryScalable')
  require(Matrix)
  
  sample.size = length(unique(mat$SAMPLE_ID))  
  output_list = vector(mode = 'list')
  uniq_sampleID = table(mat$SAMPLE_ID)
  uniq_sampleID = names(uniq_sampleID[uniq_sampleID > 1])
    
  sample_kept = vector()
  input_matrix = vector(mode = 'list')
    
  # make some modification to the input matrix
  for(i in 1:length(uniq_sampleID)){
    tmp = mat[mat$SAMPLE_ID %in% uniq_sampleID[i], ]
    tmp = tmp[order(tmp$ccf, decreasing = T), ]
    tmp = tmp[!duplicated(tmp$Alteration), ]
      
    if(nrow(tmp) > 1) {                                                
      sample_kept[i] = T
      tmp = tmp[sample(nrow(tmp)), ]
      tmp$ccf[tmp$ccf > .9] = 1
          
      input_mat = data.frame(t(combn(tmp$Alteration, m = 2)))
      input_mat$X1 = as.character(input_mat$X1)
      input_mat$X2 = as.character(input_mat$X2)
          
      player1 = tmp$ccf[match(input_mat$X1, tmp$Alteration)]
      player2 = tmp$ccf[match(input_mat$X2, tmp$Alteration)]
          
      input_mat$outcome = NA
      input_mat$outcome[which(player1 > player2)] = 'W1'
      input_mat$outcome[which(player2 > player1)] = 'W2'
      input_mat$outcome[which(player2 == player1)] = 'D'
      input_matrix[[i]] = input_mat
      } else {
        sample_kept[i] = F
      }
    }
    
    # merge the whole matrix
    sample_kept = uniq_sampleID[sample_kept]
    input_matrix = do.call(rbind, input_matrix)
    colnames(input_matrix)[1:2] = c('player1', 'player2')
    
    input_matrix_4col = codes_to_counts(input_matrix, c('W1', 'W2', 'D'))
    if(length(sample_kept) < 20) warning('Number of samples too low (< 20)')
    
    input_matrix_4col = input_matrix_4col[!is.na(input_matrix_4col$item1wins), ]
    
    
    ###################################
    ## BradleyTerry modeling
    input_mat_btdata = btdata(input_matrix_4col)
    
    
    ###################################
    ## extract 'dead-end' components (all games lost or won)
    dead_end_matrix = as.data.frame(as.matrix(input_mat_btdata$wins))
    
    games.overview = rowSums(dead_end_matrix)
    dead_end_loser = names(games.overview)[which(games.overview == 0)]
    
    dead_end_components = as.character(unlist(input_mat_btdata$components[lapply(input_mat_btdata$components, length) < 2]))
    
    single_end_out = data.frame()
    for(i in 1:length(dead_end_components)){
      sub = input_matrix_4col[input_matrix_4col$player1 == dead_end_components[i] | 
                                input_matrix_4col$player2 == dead_end_components[i], ]
      wins1 = sum(ifelse(sub$player1 == dead_end_components[i] & sub$item1wins == 1, TRUE, FALSE))
      wins2 = sum(ifelse(sub$player2 == dead_end_components[i] & sub$item2wins == 1, TRUE, FALSE))
      loss1 = sum(ifelse(sub$player1 == dead_end_components[i] & sub$item1wins == 0, TRUE, FALSE))
      loss2 = sum(ifelse(sub$player2 == dead_end_components[i] & sub$item2wins == 0, TRUE, FALSE))
      dead_end_summary = data.frame(Gene = dead_end_components[i],
                                    n.games = nrow(sub),
                                    wins = sum(wins1, wins2),
                                    lost = sum(loss1, loss2))
      
      single_end_out = rbind(single_end_out, dead_end_summary)
      
    }
    
    
    ## update data; exclude 'dead-end' components
    input_matrix_4col = input_matrix_4col[!input_matrix_4col$player1 %in% unique(single_end_out$Gene) &
                                            !input_matrix_4col$player2 %in% unique(single_end_out$Gene), ]
    
    ## make model fit and select gene comparisions
    input_mat_btdata = btdata(input_matrix_4col, return_graph = F)
    BT.updated.model = btfit(input_mat_btdata, a = 1) # MLE estimation

    
    ###################################
    ## look into comparisions (matches)
    if(!'full_dataset' %in% names(BT.updated.model$N)){
      max_size = names(which.max(lapply(BT.updated.model$N, function(x) length(x@i))))
      max_size = as.numeric(max_size)
      comparison.matrix = as.data.frame(as.matrix(BT.updated.model$N[[max_size]]))
    } else {
      comparison.matrix = as.data.frame(as.matrix(BT.updated.model$N$full_dataset))
    }
    
    
    comparison.out = data.frame()
    game.cutoff = game.cutoff
    
    for(i in 1:ncol(comparison.matrix)){
      
      if(any(comparison.matrix[, i] > game.cutoff)){
        goi = which(comparison.matrix[, i] > game.cutoff)
        goi.names = row.names(comparison.matrix)[goi]
        out = data.frame(player1 = colnames(comparison.matrix)[i],
                         player2 = goi.names,
                         n.games = comparison.matrix[goi.names, i])
        out = unique.data.frame(out)
        comparison.out = rbind(comparison.out, out)
      }
    }
    
    cols = c(1, 2)
    
    for(i in 1:nrow(comparison.out)){
      comparison.out[i, cols] = sort(comparison.out[i, cols])
    }
    
    comparison.out = comparison.out[!duplicated(comparison.out), ]

    
    ## select only those games which are above games_cutoff
    comparison.out$merge = paste(comparison.out$player1, comparison.out$player2, sep = '.')
    comparison.out$merge2 = paste(comparison.out$player2, comparison.out$player1, sep = '.')
    input_matrix_4col$merge = paste(input_matrix_4col$player1, input_matrix_4col$player2, sep = '.')
    
    input_matrix_4col = input_matrix_4col[input_matrix_4col$merge %in% comparison.out$merge | 
                                            input_matrix_4col$merge %in% comparison.out$merge2,, drop = F]
    
    input_matrix_4col$merge = NULL
    comparison.out$merge = NULL
    comparison.out$merge2 = NULL
    
    ## extract tournament participants (genes)
    remaining.genes = unique(as.vector(as.matrix(input_matrix_4col[, c(1, 2)])))
    remaining.genes = as.character(remaining.genes)
    
    
    ###################################
    ## add pseudo (reference) group: GeneX
    baseline.model = data.frame()
    
    for(i in 1:length(remaining.genes)){
      frame.out = data.frame(player1 = c(remaining.genes[i], 'GeneX'),
                             player2 = c('GeneX', remaining.genes[i]),
                             item1wins = c(1, 1),
                             item2wins = c(0, 0))
      baseline.model = rbind(baseline.model, frame.out)
    }
    
    input_matrix_4col = rbind(input_matrix_4col, baseline.model)
    
    
    ## BT MLE model fit 
    input_mat_btdata = btdata(input_matrix_4col, return_graph = F)
    BT.updated.model.MLE = btfit(input_mat_btdata, a = 1) # MLE estimation
    BT.updated.model.MAP = btfit(input_mat_btdata, a = 1.1) # MAP estimation
    
    output.capture = utils::capture.output(summary(input_mat_btdata))
    print(paste0('BT model input: ', output.capture[3], ' | MLE is used'))
    print(paste0(length(remaining.genes), ' genes are considered in the MLE BT model'))
    
    ## MLE and MAP output (summary)
    BT.model.MLE = summary(BT.updated.model.MLE, SE = T, ref = 'GeneX')
    BT.model.MAP = summary(BT.updated.model.MAP, SE = T, ref = 'GeneX')
    
    ## calculate incidence of mutations
    incidence.out = data.frame()
    freq.summary = as.data.frame(xtabs(~mat$SAMPLE_ID + mat$Alteration))
    colnames(freq.summary) = c('SAMPLE_ID', 'Alteration', 'Freq')
    
    for(i in 1:length(remaining.genes)){
      incidence = data.frame(Gene = remaining.genes[i],
                             freq = (length(freq.summary$Alteration[which(freq.summary$Alteration == remaining.genes[i] & freq.summary$Freq > 0)]) / sample.size) * 100)
      incidence.out = rbind(incidence.out, incidence)
    }
    
    
    
    ###################################
    # concentrate on tissue specific genes
    if(!is.null(plotting_refinement)){
      # mut frequencies must be over 4%
      data.plot = BT.model.MLE$item_summary
      data.plot$item = gsub('\\s+', replacement = '', x = data.plot$item)
      data.plot = data.plot[data.plot$item %in% incidence.out[which(incidence.out$freq > 4), 'Gene'], ]
      incidence.out = incidence.out[incidence.out$freq > 4, ]
      
      # just include genes with cancer-context
      # COSMIC; TIER 1 genes related to cancer-selected
      cancer.genes = gsub(pattern = '\\s+', replacement = '', remaining.genes)
      cancer.genes = cancer.genes[!grepl(pattern = '-', x = cancer.genes) & !grepl(pattern = '\\+', x = cancer.genes)]
      
      Cosmic_Cancer = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/cancer_gene_census_cosmic.csv', sep = ',')
      groups = as.character(COSMIC_cancer_type_abbreviation)
      
      cosmic.out = data.frame()
      for(i in 1:length(groups)){
        cosmic.sub = Cosmic_Cancer[grepl(pattern = groups[i], x = Cosmic_Cancer$Tumour.Types.Somatic.) & Cosmic_Cancer$Tier == 1, ]
        cosmic.out = rbind(cosmic.out, cosmic.sub)
      }
      
      # Bailey_comprehensive cancer gene list
      Bailey_consens = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Bailey_consens_cancergenes.txt', sep = '\t')
      Bailey_consens$Rescue.Notes = NULL
      
      group.bailey = as.character(Bailey_cancer_type_abbreviation)
      Bailey_selected = Bailey_consens$Gene[which(Bailey_consens$dLBC %in% group.bailey)]
      
      # comprehensive cancer-type specific genes
      cancer.genes.merged = c(unique(cosmic.out$Gene.Symbol), Bailey_selected)
      cancer.genes.merged = cancer.genes.merged[!duplicated(cancer.genes.merged)]
      
      ## prepare plotting output
      data.plot.scna = data.plot[grepl(pattern = '-', x = data.plot$item) | grepl(pattern = '\\+', x = data.plot$item), ]
      data.plot.mut = data.plot[!data.plot$item %in% data.plot.scna$item, ]
      data.plot.mut = data.plot.mut[data.plot.mut$item %in% cancer.genes.merged, ]
      
      data.plot = rbind(data.plot.mut, data.plot.scna)
      
    } else {
        data.plot = data.frame(output = 'empty')
      }
      
    
    # prepare output
    output_list[[1]] = input_matrix_4col
    output_list[[2]] = single_end_out
    output_list[[3]] = dead_end_loser
    output_list[[4]] = comparison.out
    output_list[[5]] = remaining.genes
    output_list[[6]] = BT.model.MLE
    output_list[[7]] = BT.model.MAP
    output_list[[8]] = incidence.out
    output_list[[9]] = data.plot
    output_list[[10]] = sample_kept
    
    names(output_list) = c('BT_input_matrix',
                           'dead_end_components',
                           'all_games_lost',
                           'gene_pairs_which_overcome_gamesThreshold',
                           'genes_considered_modelling',
                           'BT.MLE_output',
                           'BT.MAP.output',
                           'mut.frequencies',
                           'BT.model.plot',
                           'samples_kept')
    
    return(output_list)
    
}


# test = BradleyTerryUpdate(mat = NET.maf, game.cutoff = 3)

