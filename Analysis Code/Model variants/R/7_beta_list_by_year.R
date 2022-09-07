# save beta_by_year matrices of every year in list

library(tosca)
library(ldaPrototype)
library(rollinglda)

meta <- readRDS("./input data/meta.RDS") # meta dataframe
years_unique <- sort(unique(meta$PY))

sim_mean_names_by_year <- readRDS("./Model variants/RDS/sim_mean_names_by_year.RDS")
top_candidates <- names(head(sort(sim_mean_names_by_year, decreasing = TRUE), 5))



for (j in 1:length(top_candidates)){
  
  # docs_all from candidate model
  docs_all <- readRDS(paste0("./Model variants/", top_candidates[j],"/docs_all.RDS"))
  
  # complete candidate model
  lda <- readRDS(paste0("./Model variants/", top_candidates[j],"/lda.RDS"))
  
  K = getK(lda)
  eta = getEta(lda) # = 1/K
  assignments = getAssignments(lda)
  vocab = colnames(getTopics(lda))
  years = meta$PY[match(names(docs_all), meta$DFK)]
  
  
  topics_chunks = lapply(sort(unique(years)), function(x){
    tmp = table(factor(unlist(assignments[years == x])+1, levels = 1:K), 
                factor(unlist(lapply(docs_all[years == x], function(y) y[1,]))+1, levels = seq_len(length(vocab))))
    tmp = matrix(as.integer(tmp), nrow = K)
    colnames(tmp) = vocab
    tmp
  })
  
  topics_chunks_all <- topics_chunks
  
  beta_list_by_year <- list()
  for (i in 1:length(sort(unique(years)))){
    beta_list_by_year[[i]] <- (topics_chunks_all[[i]] + eta)/(rowSums(topics_chunks_all[[i]]) + ncol(topics_chunks_all[[i]]) * eta)
    
  }
  
  saveRDS(beta_list_by_year, file = paste0("./Model variants/", top_candidates[j], "/beta_list_by_year.RDS"))
  
  }

  



