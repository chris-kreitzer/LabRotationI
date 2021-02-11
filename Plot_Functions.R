## Plot BradleyTerry Output
## 
## plotting items are derived from BradleyTerryUpdate function; see script
## 05/01/2021


## plot in comes from the BradleyTerryUpdate function:
library(ggplot2)
library(cowplot)

plot.in = x$BT.model.plot
helper.lines.end = length(unique(plot.in$item)) + 1
plot.in$category = ifelse(plot.in$estimate > 0, 'clonal', 'subclonal')

cut_line = plot.in[plot.in$estimate >=0, ]
cut_line = cut_line$item[which.min(cut_line$estimate)]
cut_point = which(plot.in$item == cut_line)
cut_point = (length(unique(plot.in$item)) - cut_point) + 0.5

cancer_type = 'Colorectal Adenocarcinoma'
cancer_type_sample_size = length(x$samples_kept)


## base point range plot; BT model (strength) estimates;
## the relative ranking (mutations which precedes another)
BT.MLE.model = ggplot(plot.in,
                      aes(x = reorder(item, estimate), 
                          y = estimate, 
                          ymin = estimate - SE, 
                          ymax = estimate + SE)) +
  
  geom_pointrange(size = 0.5) +
  
  geom_linerange() +
  
  coord_flip() +
  
  geom_vline(xintercept = seq(0.5, helper.lines.end, 1), 
             size = 0.2, 
             linetype = 'dashed', 
             color = 'grey') +
  
  geom_vline(xintercept = cut_point, 
             color = '#8B2D76', 
             size = 1.2, 
             linetype = 'dashed') +
  
  scale_y_continuous(expand = c(0, 0)) +
  
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 12, 
                                   face = 'bold',
                                   color = 'black'),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(size = 0.6),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  
  labs(x = '', y = '', title = paste0(cancer_type, ' (n=', cancer_type_sample_size, ')'))

BT.MLE.model


## Mutational frequencies of kept mutations in BT modeling approach
plot.in = merge(plot.in, x$mut.frequencies, by.x = 'item', by.y = 'Gene', all.x = T)
plot.in$bar_color = ifelse(grepl(pattern = '-', x = plot.in$item), '#00008b',
                           ifelse(grepl(pattern = '\\+', x = plot.in$item), '#F5201E', 'grey15'))

upper.limit = max(ceiling(plot.in$freq))

mut.freq = ggplot(plot.in, 
                  aes(x = reorder(item, estimate), 
                      y = freq)) +
  
  geom_bar(stat = 'identity', fill = plot.in$bar_color) +
  
  coord_flip() +
  
  scale_y_continuous(position = 'right', 
                     expand = c(0, 0), 
                     breaks = pretty(x = plot.in$freq, n = 5, min.n = 0),
                     limits = c(0, upper.limit + 3)) +
  
  geom_hline(yintercept = seq(0, 80, 20), size = 0.3, linetype = 'solid', col = 'grey75') +
  
  geom_text(data = plot.in, aes(label = round(freq)), hjust = 1, col = 'white', size = 3) +
  
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.line.x.top = element_line(size = 0.6),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  labs(y = 'Mutational Frequency (%)')



## plot the driver mutation distribution among retained genes;
## x$BT.model.plot ~ whereas x = output of BradleyTerryUpdate function
plot.in = x$BT.model.plot
plot.in$item = as.factor(plot.in$item)
x.sort = reorder(plot.in$item, plot.in$estimate)
driver.distribution$region = factor(driver.distribution$region, levels = levels(x.sort))
driver.distribution = driver.distribution[!is.na(driver.distribution$region), ]
driver.distribution$Var1 = factor(driver.distribution$Var1, 
                                  levels = c('nonDriver', 'Driver'))


# make plot
driver.contribution = ggplot(driver.distribution, 
                             aes(x = region, 
                                 y = Freq, 
                                 fill = Var1)) +
  
  geom_bar(position = 'fill', stat = 'identity') +
  
  scale_y_continuous(position = 'right', 
                     expand = c(0,0), 
                     breaks = c(0, 0.5, 1)) +
  
  scale_fill_manual(values = c('Driver' = '#066832',
                               'nonDriver'= '#98CC9C'),
                    name = 'Category') +
  
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.line.x.top = element_line(size = 0.6),
        aspect.ratio = 4,
        legend.position = 'bottom',
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  
  labs(y = 'Relative contribution') +
  
  coord_flip()


## combine all three plots into one grid
plot_grid(BT.MLE.model, 
          mut.freq, 
          driver.contribution, 
          ncol = 3, 
          align = 'h', 
          axis = 'tblr', 
          rel_widths = c(4, 1, 1))

