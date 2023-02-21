library(tidyverse)
library(tidygraph)
library(ggdag)
library(ggraph)
library(cowplot)
library(patchwork)

dagified <- dagify('SnakeQUANT' ~ 'fastq',
                   'counts' ~ 'SnakeQUANT',
                   'SnakePOP' ~ 'CellType',
                   'SnakePOP' ~ 'counts',
                   'LatentDims' ~ 'SnakePOP',
                   'cluster' ~ 'LatentDims',
                   'umap' ~ 'LatentDims' ,
                   'ML' ~ 'CellType',
                   'projectedCellType' ~ 'ML',
                   'ML' ~ 'LatentDims',
                   'SnakeSCEIAD' ~ 'counts',
                   'SnakeSCEIAD' ~ 'projectedCellType',
                   'SnakeSCEIAD' ~  'CellType',
                   'SnakeSCEIAD' ~   'cluster',
                   'SnakeSCEIAD' ~   'umap',
                   'SnakeSCEIAD' ~ 'SampleInfo',
                   'SnakeQUANT' ~ 'SampleInfo',
                   'SnakePOP' ~ 'SampleInfo',
                   'DiffTesting' ~ 'SnakeSCEIAD',
                   'scEiaD' ~ 'SnakeSCEIAD',
                   'scEiaD' ~ 'DiffTesting',
                   'plae' ~ 'scEiaD')

tidy_dagitty(dagified)
set.seed(51345)
workflow_plot <- tidy_dagitty(dagified) %>% as_tbl_graph() %>% 
  mutate(type = case_when(name %in% c('fastq', 'CellType', 'SampleInfo') ~ 'Input',
                          name %in% c('SnakeQUANT', 'counts') ~ 'SQ',
                          name %in% c('umap', 'projectedCellType','cluster', 'counts', 'LatentDims', 'SnakePOP', 'ML') ~ 'SP',
                          name %in% c('scEiaD','DiffTesting', 'SnakeSCEIAD') ~ 'SS'),
         name = case_when(name == 'counts' ~ 'ambient cleaned counts',
                          name == 'CellType' ~ 'Published\nCell Types',
                          name == 'projectedCellType' ~ 'Learned\nCell Types',
                          name == 'DiffTesting' ~ 'Diff\nTesting',
                          name == 'scVIrefquery' ~ 'Optimal\nscVI Model',
                          name == 'LatentDims' ~ 'scVI Batch Corrected\nLatent Dims',
                          name == 'SampleInfo' ~ 'Sample\nInfo',
                          name == 'plae' ~ 'plae.nei.nih.gov',
                          TRUE ~ name)) %>% 
  ggraph(layout = 'auto') + 
  geom_edge_bend2(arrow = arrow(length = unit(3, 'mm')), 
                 end_cap = circle(5, 'mm'), start_cap = circle(2, 'mm')) + 
  #geom_node_point(size = 10) +
  geom_node_label(aes(label = name, color = type), alpha = 1) +
  scale_color_manual(values = c('Red','Blue','Purple','DarkGreen'), na.value = 'Black') +
  coord_cartesian(xlim=c(-1,6), ylim=c(0,12)) +
  theme_nothing()

workflow_plot  + theme(plot.margin = margin(0,0,0,0 , "cm"))
