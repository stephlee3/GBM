library(tidyverse)
library(openxlsx)

get_summary = function(cell_popu_a,cell_popu_b, type = "all"){
  if (type == "immune"){
    result = readRDS(paste0("../data/FC_signif_immune/",cell_popu_a,"_",cell_popu_b,".rds"))
  } else 
    result = readRDS(paste0("../data/FC_signif/",cell_popu_a,"_",cell_popu_b,".rds"))
   
  cell_type_a = gsub('(.*):(.*)','\\1',names(result))
  cell_type_b=  gsub('(.*):(.*)','\\2',names(result))
  
  signif_num =sapply(1:length(result),function(i){
    nrow(result[[i]] %>% filter(signif))
  })
  names(signif_num) = names(result)
  print(str(signif_num))
  
  signif_result = lapply(1:length(result),function(i){
    result[[i]] %>% filter(signif) %>%
    mutate(LR = paste0(Ligand,"-",Receptor)) %>%
      mutate(ct_pair = names(result)[i]) %>%
      mutate(a = cell_type_a[i],
             b = cell_type_b[i])
  })
  signif_result = signif_result[order(signif_num)]
  signif_result = do.call(rbind, signif_result)
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
cell_popu_a = c("myeloid","myeloid","tumor","t","t")
cell_popu_b = c("tumor","t","myeloid","myeloid","t")
for(i in 1:length(cell_popu_a)){
  signif_num = get_summary(cell_popu_a[i],cell_popu_b[i])
  write.xlsx(signif_num, file = paste0("../data/FC_signif/",cell_popu_a[i],"_",cell_popu_b[i],"_summary.xlsx"))
  saveRDS(signif_num, paste0("../data/FC_signif/",cell_popu_a[i],"_",cell_popu_b[i],"_summary.rds"))
}



