library(tidyverse)
library(openxlsx)

get_summary = function(cell_popu_a,cell_popu_b, type = 'all', method = 'norm', cutoff = NULL, direction = 'pos'){
  if (type == 'all' & method == 'norm')
    result = readRDS(paste0("/path/to/normalized/score",cell_popu_a,"_",cell_popu_b,".rds"))
  if (type == 'immune' & method == 'norm')
    result = readRDS(paste0("/path/to/immune/normalized/score",cell_popu_a,"_",cell_popu_b,".rds"))  
  if (type == 'all' & method == 'orig')
    result = readRDS(paste0("/path/to/raw/score",cell_popu_a,"_",cell_popu_b,".rds"))
  if (type == 'immune' & method == 'orig')
    result = readRDS(paste0("/path/to/immune/raw/score",cell_popu_a,"_",cell_popu_b,".rds")) 
    
  result = readRDS(paste0("../data/",method,"/",type,"/",
                          cell_popu_a, "_", cell_popu_b, ".rds"))
  
  cell_type_a = gsub('(.*):(.*)','\\1',names(result))
  cell_type_b=  gsub('(.*):(.*)','\\2',names(result))
  
  ct_pair = names(result)

  result = lapply(1:length(result),function(i){
    result[[i]] %>% 
      mutate(signif = case_when(
        direction == 'pos' ~ (fdr<0.05 & enrich > cutoff),
        direction == 'neg' ~ (fdr<0.05 & enrich < -cutoff),
      )) %>%
      mutate(ct_pair = names(result)[i]) %>%
      mutate(a = cell_type_a[i],
             b = cell_type_b[i]) %>%
      mutate(Ligand = gsub('(.*)-(.*)','\\1',LR),
             Receptor = gsub('(.*)-(.*)','\\2',LR))
  })
  names(result)= ct_pair
  
  signif_num =sapply(1:length(result),function(i){
    nrow(result[[i]] %>% filter(signif))
  })
  names(signif_num) = names(result)
  
  signif_result = lapply(1:length(result),function(i){
    result[[i]] %>% filter(signif)
  })
  signif_result = signif_result[order(signif_num)]
  signif_result = do.call(rbind, signif_result)
  if(nrow(signif_result) == 0){
    return(NULL)
  }
  LR_count = table(signif_result$LR)
  
  la = lapply(unique(cell_type_a), function(ct){
    idx = which(cell_type_a == ct)
    la = lapply(idx, function(i){
      result[[i]] %>% filter(signif) %>% pull(Ligand)
    })
    la = do.call(c,la)
    if (length(la)>=1) la = unique(la)
    la
  })
  names(la) = unique(cell_type_a)
  
  rb = lapply(unique(cell_type_b), function(ct){
    idx = which(cell_type_b == ct)
    rb = lapply(idx, function(i){
      result[[i]] %>% filter(signif) %>% pull(Receptor)
    })
    rb = do.call(c,rb)
    if (length(rb)>=1) rb = unique(rb)
    rb
  })
  names(rb) = unique(cell_type_b)
  num_a = sapply(la,length)
  num_b = sapply(rb,length)
  
  LR_count = as.data.frame(LR_count); colnames(LR_count) = c("LR","count")
  LR_num = data.frame(cell_type_pair = names(signif_num),count = signif_num)
  L_num = data.frame(sender_cell = names(num_a),count = num_a)
  R_num = data.frame(receiver_cell = names(num_b),count = num_b)
  ans = list(signif_result = signif_result,
             LR_count = LR_count %>% arrange(desc(count)),
             LR_num = LR_num %>% arrange(desc(count)),
             L_num = L_num %>% arrange(desc(count)),
             R_num = R_num %>% arrange(desc(count)))
  ans
  
}
cell_popu_a = c("myeloid","myeloid", "tumor","t") 
cell_popu_b = c("tumor","t", "myeloid","myeloid")

opt = expand.grid(
  method = c('norm'),
  type = c('all','immune'), 
  direction = c('pos','neg'),
  cutoff = c(0.2,0.3,0.4,0.5)) %>%
  rbind(expand.grid(
  method = c('orig'),
  type = c('all','immune'), 
  direction = c('pos','neg'),
  cutoff = c(5,10,15,20))
  ) %>%
  arrange(method, type, direction)
                    


for(t in 1:nrow(opt)){
  print(opt[t,])
  for(i in 1:length(cell_popu_a)){
  signif_num = get_summary(cell_popu_a[i],cell_popu_b[i], type = opt$type[t], method = opt$method[t], cutoff = opt$cutoff[t], direction = opt$direction[t] )
  write.xlsx(signif_num, file = paste0("../data/",opt$method[t],"/",opt$type[t],"/",opt$method[t],"_",opt$type[t],"_", cell_popu_a[i],"_",cell_popu_b[i],"_", opt$cutoff[t],"_", opt$direction[t],"_summary.xlsx"))
  saveRDS(signif_num, paste0("../data/",opt$method[t],"/",opt$type[t],"/",opt$method[t],"_",opt$type[t],"_", cell_popu_a[i],"_",cell_popu_b[i],"_", opt$cutoff[t],"_", opt$direction[t],"_summary.rds"))
}
}
