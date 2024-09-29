library(tidyverse)
library(openxlsx)
library(limma)
folder_path = "/path/to/pseudo/bulk/data"
patient_info = readRDS("../data/patient_info.rds") %>%
  filter(Treatment == "Untreated" & Tumor.Grade %in% c("IV","III","II"))
cluster_info = readRDS("../data/cluster_info.rds")
cluster_info$cell_type = as.character(cluster_info$cell_type)

myeloid_clu = cluster_info %>% filter(cell_popu == 'myeloid') %>% pull(cell_type)
tumor_clu = cluster_info %>% filter(cell_popu == 'tumor') %>% pull(cell_type)
t_clu = cluster_info %>% filter(cell_popu == 't') %>% pull(cell_type)

t_rm = c("CYC1","CYC2")
myeloid_rm = c("NK1","NK2","B cells", "Plasma")
myeloid_clu = setdiff(myeloid_clu, myeloid_rm)
t_clu = setdiff(t_clu,t_rm)
tumor_clu = c(paste0("T",1:10))

#
## generate cell type pairs
get_ct_pair = function(cell_popu_a, cell_popu_b){
  cell_type_a = get(paste0(cell_popu_a,"_clu"))#cluster_info %>% filter(cell_popu == cell_popu_a) %>% pull(cell_type)
  cell_type_b = get(paste0(cell_popu_b,"_clu")) #cluster_info %>% filter(cell_popu == cell_popu_b) %>% pull(cell_type)
  expand.grid(a = cell_type_a, b = cell_type_b) %>% 
    data.frame()
}

get_norm_pb = function(cell_popu, method ='norm'){
  cell_popu_alt = cell_popu
  pt = patient_info %>% filter(`cell_popu` == cell_popu_alt) %>% pull(Patient)
  pb = lapply(pt,function(i){
    d = readRDS(paste0(folder_path,cell_popu,"/pseudo_bulk/",i,".rds"))
    cn = colnames(d)
    if (method == 'norm'){
      d = t(apply(d,1,function(x){(x-min(x))/(max(x)-min(x))}))
      d[is.nan(d)] = 0
    }
    colnames(d) = paste0(i,':',cn)
    d
  })
  pb = do.call(cbind,pb)
  pb
}
pb_myeloid_norm=  get_norm_pb("myeloid")
pb_tumor_norm=  get_norm_pb("tumor")
pb_t_norm =  get_norm_pb("t")

pb_myeloid_orig = get_norm_pb("myeloid", method = 'orig')
pb_tumor_orig = get_norm_pb("tumor", method = 'orig')
pb_t_orig = get_norm_pb("t", method = 'orig')

get_norm_score = function(cell_popu_a, cell_popu_b, cell_type_a, cell_type_b, type = 'all', method ='norm'){
  pb_a = get(paste0("pb_",cell_popu_a,'_',method))
  pb_b = get(paste0("pb_",cell_popu_b,'_',method))
  ct_a = gsub('(.*):(.*)','\\2',colnames(pb_a))
  ct_b = gsub('(.*):(.*)','\\2',colnames(pb_b))
  pt_a = gsub('(.*):(.*)','\\1',colnames(pb_a))
  pt_b = gsub('(.*):(.*)','\\1',colnames(pb_b))
  
  idx_a= which(ct_a == cell_type_a);idx_b = which(ct_b == cell_type_b)
  pt_a = pt_a[idx_a]; pt_b = pt_b[idx_b]
  pt = intersect(pt_a,pt_b)
  ptlist = patient_info %>% filter(cell_popu == cell_popu_a)
  grade = ptlist$Tumor.Grade[match(pt, ptlist$Patient)]
  
  pb_a = pb_a[,idx_a]; pb_b = pb_b[,idx_b]
  print(str(pb_a))
  print(str(pb_b))
  
  if (type == 'all'){
    LR_list = read.xlsx("/path/to/LR/list",startRow = 2) 
    LR_list = LR_list %>%
      mutate(LR = paste0(Ligand,"-", Receptor))
  }
  if (type == 'immune'){
    LR_list = read.table("./path/to/immune/LR/list",header = T)
    LR_list = LR_list %>%
      mutate(LR = paste0(Ligand,"-", Receptor))
  }
  
  if(nrow(LR_list)==0)
    return(rep(NA,length(pt)))
  
  score = t(sapply(1:nrow(LR_list), function(i){
    ligand = LR_list$Ligand[i]
    receptor = LR_list$Receptor[i]
    if (!ligand %in% rownames(pb_a) | !receptor %in% rownames(pb_b))
      return(rep(NA,length(pt)))
    ligand_expr = pb_a[ligand,match(pt,pt_a)]
    recep_expr = pb_b[receptor,match(pt,pt_b)]
    score = ligand_expr * recep_expr
  }))
  colnames(score) = pt
  rownames(score) = LR_list$LR
  score
}

get_signif_LR = function(cell_popu_a, cell_popu_b, cell_type_a,cell_type_b, type = 'immune', method = 'norm', enrich_cutoff = 0.2){
  s = try(get_norm_score(cell_popu_a,cell_popu_b, cell_type_a,cell_type_b, type, method))
  if (class(s)=='try-error') 
    return(data.frame(LR = character(),unadj.p = numeric(),fdr = numeric(),enrich = numeric()))
  s = na.omit(s)
  if (length(s)==0)
    return(data.frame(LR = character(),unadj.p = numeric(),fdr = numeric(),enrich = numeric()))
 
  
  pt = colnames(s)
  ptlist = patient_info %>% filter(cell_popu == cell_popu_a)
  grade = ptlist$Tumor.Grade[match(pt, ptlist$Patient)]
  gradecoef = ifelse(grade == 'IV',1,-1)
  enrich_score = apply(s,1,function(x) mean(x[grade=='IV'])-mean(x[grade!='IV']))
 
  unadj_pval = sapply(1:nrow(s), function(i){
    sh = s[i,which(grade=='IV')]; sl = s[i,which(grade!='IV')]
    if (length(sh)==0 | length(sl)==0)
      return(NA)
    wt = wilcox.test(sl,sh)
    wt$p.value
  })
  fdr = p.adjust(unadj_pval,method = "fdr")
  result = data.frame(LR = rownames(s), unadj.p = unadj_pval,fdr = fdr, enrich = enrich_score) %>% 
    arrange(desc(enrich))
 
  signif_result = result %>% filter(fdr <0.05 & enrich > enrich_cutoff)
  print(paste(nrow(signif_result),"signif results found!"))
  print(head(signif_result,20))
  print("----------------------------------------------")
  return(result)
}



cell_popu_a = c("myeloid","myeloid","tumor","t") #,"t")
cell_popu_b = c("tumor","t","myeloid","myeloid") #,"t")
type = c("immune","all") 
method = c('norm','orig')
opt = expand.grid(method = method, type = type)

for (t in 1:nrow(opt)){
  print(opt[t,])
  for(k in  1:length(cell_popu_a)){
  ct_pair = get_ct_pair(cell_popu_a[k],cell_popu_b[k])
  signif.LR = lapply(1:nrow(ct_pair), function(i){
    print(paste0(i,": ", ct_pair$a[i],":",ct_pair$b[i]))
    get_signif_LR(cell_popu_a[k],cell_popu_b[k], ct_pair$a[i], ct_pair$b[i],
                  type = opt$type[t], method = opt$method[t])
  })
  names(signif.LR) = paste0(ct_pair$a,':', ct_pair$b)
  dir.create(paste0("../data/",opt$method[t]))
  dir.create(paste0("../data/",opt$method[t],"/",opt$type[t]))
  saveRDS(signif.LR, paste0("../data/",opt$method[t],"/",opt$type[t],"/",cell_popu_a[k],"_",cell_popu_b[k],".rds"))
}
}