library(tidyverse)
library(openxlsx)
library(limma)
folder_path = "/path/to/pseudo/bulk/data"
patient_info = readRDS("../data/patient_info.rds") %>%
  filter(Treatment == "Untreated" & Tumor.Grade %in% c("IV","III","II"))
cluster_info = readRDS("../data/cluster_info.rds")

myeloid_clu = cluster_info %>% filter(cell_popu == 'myeloid') %>% pull(cell_type)
tumor_clu = cluster_info %>% filter(cell_popu == 'tumor') %>% pull(cell_type)
t_clu = cluster_info %>% filter(cell_popu == 't') %>% pull(cell_type)


LR_list = read.xlsx("/path/to/LR/list",startRow = 2)

## generate cell type pairs
get_ct_pair = function(cell_popu_a, cell_popu_b){
  cell_type_a = cluster_info %>% filter(cell_popu == cell_popu_a) %>% pull(cell_type)
  cell_type_b = cluster_info %>% filter(cell_popu == cell_popu_b) %>% pull(cell_type)
  expand.grid(a = cell_type_a, b = cell_type_b) %>% 
    data.frame()
}

## load pseudobulk into concatenated matrix
read_pb = function(cell_popu){
  cell_popu_alt = cell_popu
  pt = patient_info %>% filter(`cell_popu` == cell_popu_alt) %>% pull(Patient)
  pb = lapply(pt,function(i){
    d = readRDS(paste0(folder_path,cell_popu,"/pseudo_bulk/",i,".rds"))
    colnames(d) = paste0(i,':',colnames(d))
    d
  })
  pb = do.call(cbind,pb)
  pb
}
pb_myeloid=  read_pb("myeloid")
pb_t = read_pb("t")
pb_tumor = read_pb("tumor")

## apply limma to identify DE genes
read_DEG = function(cell_popu, cell_type, label = 'Ligand', logFC_cutoff = 0, logFC_quantile = 0.025){
  cell_popu_alt = case_when(
    cell_popu =='myeloid' ~ 'Myeloid',
    cell_popu =='tumor' ~ 'Tumor',
    cell_popu == 't' ~ 'Lymph'
  )

  DE_result = readRDS(paste0("/path/to/DEG/result",
                               substr(cell_popu_alt,1,1),"/",gsub(' ','',cell_type),".rds"))

  cutneg = - logFC_cutoff; cutpos = logFC_cutoff;
  print(cutneg)
  print(cutpos)
  DE_result = DE_result %>%
    mutate(gene = rownames(.)) %>%
    mutate(!!(label):= gene) %>%
    mutate(signif = case_when(
      FDR < 0.05 & Foldchange > cutpos ~ "up",
      FDR < 0.05 & Foldchange < cutneg ~ "down",
      T ~ "no"
    ))
  
                  
  DEG = LR_list %>%
    left_join(DE_result %>% select(!!(label), Foldchange, FDR, signif)) %>%
    dplyr::rename(logFC = Foldchange, adj.P.Val = FDR)



}

## cell cell interaction
get_crosstalk = function(cell_popu_a, cell_popu_b, cell_type_a, cell_type_b, logFC_cutoff = 0.5, logFC_quantile = 0.025){
  DEL = read_DEG(cell_popu_a, cell_type_a, 'Ligand', logFC_cutoff,logFC_quantile)
  DER = read_DEG(cell_popu_b, cell_type_b, 'Receptor',logFC_cutoff, logFC_quantile)
  DELR = DEL %>%
    left_join(DER, by = c("Ligand","Receptor")) %>%
    select(Ligand, Receptor, logFC.x, logFC.y, signif.x,signif.y, adj.P.Val.x, adj.P.Val.y) %>%
    dplyr::rename(logFC.L = logFC.x, logFC.R = logFC.y, 
                  signif.L = signif.x, signif.R = signif.y,
                  FDR.L = adj.P.Val.x, FDR.R = adj.P.Val.y) %>%
    mutate(score = logFC.L * logFC.R) %>%
    na.omit()
  DELR_signif = DELR %>%
    filter(!(signif.L =='no' & signif.R == 'no')) %>%
    mutate(signif = ((signif.L == 'up') & (signif.R == 'up')) | 
                    ((signif.L == 'down') & (signif.R == 'down'))) %>%
    filter(signif) %>%
    arrange(desc(score))
  print(paste(nrow(DELR_signif),"signif results found!"))
  print(head(DELR_signif))
  print("----------------------------------------------")
  return(DELR_signif)
}


######
cell_popu_a = c("myeloid","myeloid","tumor","t","t")
cell_popu_b = c("tumor","t","myeloid","myeloid","t")

for(k in  1:length(cell_popu_a)){
  ct_pair = get_ct_pair(cell_popu_a[k],cell_popu_b[k])
  DELR = lapply(1:nrow(ct_pair), function(i){
    print(paste0(i,": ", ct_pair$a[i],":",ct_pair$b[i]))
    get_crosstalk(cell_popu_a[k],cell_popu_b[k],ct_pair$a[i],ct_pair$b[i], logFC_cutoff= 0)
  })
  names(DELR) = paste0(ct_pair$a,':', ct_pair$b)
  saveRDS(DELR, paste0("../data/FC_signif/",cell_popu_a[k],"_",cell_popu_b[k],".rds"))
}



