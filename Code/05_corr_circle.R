library(tidyverse)
library(openxlsx)
library(circlize)

folder_path = "/path/to/pseudo/bulk/data"
patient_info = readRDS("../data/patient_info.rds") %>%
  filter(Treatment == "Untreated" & Tumor.Grade %in% c("IV","III","II"))
cluster_info = readRDS("../data/cluster_info.rds")


## generate cell type pairs
get_ct_pair = function(cell_popu_a, cell_popu_b){
  cell_type_a = cluster_info %>% filter(cell_popu == cell_popu_a) %>% pull(cell_type)
  cell_type_b = cluster_info %>% filter(cell_popu == cell_popu_b) %>% pull(cell_type)
  expand.grid(a = cell_type_a, b = cell_type_b) %>% 
    data.frame()
}

## load pseudobulk into concatenated matrix
get_norm_pb = function(cell_popu){
  cell_popu_alt = cell_popu
  pt = patient_info %>% filter(`cell_popu` == cell_popu_alt) %>% pull(Patient)
  pb = lapply(pt,function(i){
    d = readRDS(paste0(folder_path,cell_popu,"/pseudo_bulk/",i,".rds"))
    cn = colnames(d)
    colnames(d) = paste0(i,':',cn)
    d
  })
  pb = do.call(cbind,pb)
  pb
}
pb_myeloid=  get_norm_pb("myeloid")
pb_tumor=  get_norm_pb("tumor")
pb_t=  get_norm_pb("t")

get_LR_corr = function(cell_popu_a, cell_popu_b, cell_type_a, cell_type_b,
                       ligand, receptor){
  print(paste0(ligand,"-",receptor))
  pb_a = get(paste0("pb_",cell_popu_a))
  pb_b = get(paste0("pb_",cell_popu_b))
  ct_a = gsub('(.*):(.*)','\\2',colnames(pb_a))
  ct_b = gsub('(.*):(.*)','\\2',colnames(pb_b))
  pt_a = gsub('(.*):(.*)','\\1',colnames(pb_a))
  pt_b = gsub('(.*):(.*)','\\1',colnames(pb_b))
  
  idx_a= which(ct_a == cell_type_a);idx_b = which(ct_b == cell_type_b)
  pb_a = pb_a[,idx_a]; pb_b = pb_b[,idx_b]
  pt_a = pt_a[idx_a]; pt_b = pt_b[idx_b]
  pt = intersect(pt_a,pt_b)
  ptlist = patient_info %>% filter(cell_popu == cell_popu_a)
  grade = ptlist$Tumor.Grade[match(pt, ptlist$Patient)]
  
  ligand_expr = pb_a[ligand,match(pt,pt_a)]
  recep_expr = pb_b[receptor,match(pt,pt_b)]
  corr = cor(ligand_expr,recep_expr, method = "spearman", use= "complete")
  print(corr)
  print("---------------")
  return(corr)
}

LR_list = read.xlsx("/path/to/LR/list.xlsx")
corr = sapply(1:nrow(LR_list),function(i){
  get_LR_corr(LR_list$cell_popu_a[i], LR_list$cell_popu_b[i],
              LR_list$cell_type_a[i], LR_list$cell_type_b[i],
              LR_list$ligand[i], LR_list$receptor[i])
})
LR_list$corr = corr
LR_list = LR_list %>%
  mutate(class_b = case_when(
    cell_type_b == 'TREG' ~ 'TREG',
    cell_type_b %in% c("CD8 RM-GZMK","CD8 RM-XCL1","CD8 EX") ~ 'CD8 T',
    cell_type_b %in% c(" CD8 EFF","CD8 IFNG") ~ 'CD8 T',
    T ~ cell_type_b
  ))

#############################################################

circlize_plot = function(LR_list, cell_popu, highlight = T){
  circos.clear()
  dat = LR_list %>%
    filter(cell_popu_b == cell_popu | cell_popu_a == cell_popu) %>%
    arrange(class_b,cell_type_a) %>%
    dplyr::rename(from = ligand, to = receptor, value = corr)
  
  set.seed(123)
  link_arrow_col = 'red'
  
  
  ## redefine the colors
  if (cell_popu == 't'){
      LR_name = c(dat %>% arrange(cell_type_a) %>% pull(from),dat %>% arrange(class_b) %>% pull(to))
      grid.col = c(rep("#756bb1",nrow(dat)),rep("#43a2ca",nrow(dat)))
      names(grid.col) = LR_name
  }
  
   if (cell_popu == 'tumor'){
      LR_name = c( dat %>% filter(cell_type_a == 'E-MDSC') %>% pull(from), dat %>% filter(cell_type_b == 'E-MDSC') %>% pull(to),
                   dat %>% filter(cell_type_a == 'T4') %>% pull(from), dat %>% filter(cell_type_b == 'T4') %>% pull(to))
      grid.col = c(rep("#756bb1",nrow(dat)),rep("#e34a33",nrow(dat)))
      names(grid.col) = LR_name
      circos.par(gap.after = c(rep(1,17),5,rep(1,17),5))
  }

 
  chordDiagram(dat,
               order = LR_name,
               grid.col = grid.col,
               reduce = 0,
               directional = 1,  transparency = 1,
               link.arr.width = 2, 
               direction.type = "arrows",
               link.arr.col = link_arrow_col, link.arr.length = 1,
               link.arr.lwd = dat$value/max(dat$value)*20,
               annotationTrack = "grid", #c("grid"),#,"name"),
               annotationTrackHeight = c(0.01, 0.01),
               preAllocateTracks = list(list(
                 track.height = uh(5, "mm"),
                 track.margin = c(uh(90, "mm"),uh(0, "mm"))
               )))
  
  if (highlight){
     circos.track(track.index = 2, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", 
                niceFacing = TRUE, adj = c(-0.2, 0.5), cex = 3)
  }, bg.border = NA)
  }
 
  
  if (highlight){
    if (cell_popu == 't'){
          highlight.sector(dat$from[which(dat$cell_type_a=='E-MDSC')], track.index = 1, col = "#756bb1",text.vjust = "20mm",  
                     text = "E-MDSC", cex = 4, text.col = "black", niceFacing = TRUE)
          highlight.sector(dat$to[which(dat$class_b=='CD8 T')], track.index = 1, col = "#1c9099",text.vjust = "20mm",  
                     text = "CD8 T", cex = 4, text.col = "black", niceFacing = TRUE) 
    
    }
    
    if (cell_popu == 'tumor'){
          highlight.sector(c(dat$from[which(dat$cell_type_a=='E-MDSC')],
                             dat$to[which(dat$cell_type_b=='E-MDSC')]), track.index = 1, col = "#756bb1",text.vjust = "15mm",  
                     text = "E-MDSC", cex = 4, text.col = "black", niceFacing = TRUE)
          highlight.sector(c(dat$from[which(dat$cell_type_a=='T4')],
                             dat$to[which(dat$cell_type_b=='T4')]), track.index = 1, col = "#e34a33",text.vjust = "15mm",  
                     text = "T4", cex = 4, text.col = "black", niceFacing = TRUE)
                     
    }
  
  }
  
  if (!highlight){
    if (cell_popu == 't'){
          highlight.sector(dat$from[which(dat$cell_type_a=='E-MDSC')], track.index = 1, col = "#756bb1",text.vjust = "20mm",  
                     text = "", cex = 4, text.col = "black", niceFacing = TRUE)
          highlight.sector(dat$to[which(dat$class_b=='CD8 T')], track.index = 1, col = "#1c9099",text.vjust = "20mm",  
                     text = "", cex = 4, text.col = "black", niceFacing = TRUE) 
    }
    
    if (cell_popu == 'tumor'){
         highlight.sector(c(dat$from[which(dat$cell_type_a=='E-MDSC')],
                             dat$to[which(dat$cell_type_b=='E-MDSC')]), track.index = 1, col = "#756bb1",text.vjust = "15mm",  
                     text = "", cex = 4, text.col = "black", niceFacing = TRUE)
          highlight.sector(c(dat$from[which(dat$cell_type_a=='T4')],
                             dat$to[which(dat$cell_type_b=='T4')]), track.index = 1, col = "#e34a33",text.vjust = "15mm",  
                     text = "", cex = 4, text.col = "black", niceFacing = TRUE)
    }
    
  
  }
}

