ct_mat <- meta_filter %>% 
  filter(CellType %in% meta_filter$CellType_predict) %>% 
  filter(!is.na(CellType), !is.na(CellType_predict)) %>%
  select(CellType, CellType_predict) %>%
  group_by(CellType, CellType_predict) %>%
  summarise(Count = n()) %>% 
  mutate(Ratio = Count / sum(Count)) %>% 
  select(-Count) %>% 
  #filter(!CellType_predict %in% c('Astrocytes', 'Artery',  'Choriocapillaris', 'Tabula Muris','B-Cell','T-Cell', 'Smooth Muscle Cell', 'Unlabelled'), !CellType %in% c('Astrocytes','Artery','Choriocapillaris','Tabula Muris','B-Cell','T-Cell', 'Smooth Muscle Cell')) %>% 
  #filter(!CellType_predict %in% c('Astrocytes', 'Artery',  'Choriocapillaris', 'Tabula Muris','Smooth Muscle Cell', 'Unlabelled'), !CellType %in% c('Astrocytes','Artery','Choriocapillaris', 'Smooth Muscle Cell')) %>% 
  pivot_wider(values_from = Ratio, names_from = CellType_predict) 
ct_mat[is.na(ct_mat)] <- 0
#ct_mat <- data.frame(ct_mat)


rowN <- ct_mat$CellType
ct_mat <- ct_mat[,-1] 
ct_mat <- ct_mat[, colnames(ct_mat) %>% sort()]
ct_mat <- ct_mat %>% as.matrix()
row.names(ct_mat) <- rowN
ct_mat <- ct_mat[row.names(ct_mat) %>% rev(), ]

ct_confusion <- ComplexHeatmap::Heatmap(ct_mat, cluster_rows = FALSE, cluster_columns = FALSE, col=viridis::viridis(20), name =  'Recall', row_title = 'Published Cell Type', column_title = 'Predicted Cell Type')


ct_table <- meta_filter %>% 
  filter(CellType %in% meta_filter$CellType_predict) %>% 
  filter(!is.na(CellType), !is.na(CellType_predict)) %>%
  select(organism, CellType, CellType_predict) %>%
  group_by(organism, CellType, CellType_predict) %>%
  summarise(Count = n()) %>% 
  mutate(Ratio = Count / sum(Count))# %>% 
  #select(-Count) %>% filter(Ratio > 0.2)

ct_accuracy <- meta_filter %>% filter(Organ == 'Eye', !is.na(CellType), CellType == CellType_predict) %>% nrow() / (meta_filter %>% filter(Organ == 'Eye', !is.na(CellType), CellType == CellType_predict) %>% nrow() +
                                                                                                           meta_filter %>% filter(Organ == 'Eye', !is.na(CellType), CellType != CellType_predict) %>% nrow())
