### create topic_evo with beta by year

# rebuild keywords helper function
rebuild_keywords <- function(x){
  x <- gsub("RRR1", "(", x, ignore.case = TRUE)
  x <- gsub("RRR2", ")", x, ignore.case = TRUE)
  x <- gsub("RRR3", "/", x, ignore.case = TRUE)
  x <- gsub("RRR4_", "'", x, ignore.case = TRUE)
  x <- gsub("_as_", " as ", x, ignore.case = TRUE)
  x <- gsub("_", " ", x)
  return(x)
}


make_topic_evo_by_year <- function(lda, docs_all, meta){
  
  library(rollinglda)
  library(tosca)
  
  # see buzz_search.R
  assignments = getAssignments(lda)
  
  vocab = colnames(getTopics(lda))
  years = meta$PY[match(names(docs_all), meta$DFK)]
  K = getK(lda)
  
  topics_chunks = lapply(sort(unique(years)), function(x){
    tmp = table(factor(unlist(assignments[years == x])+1, levels = 1:K), 
                factor(unlist(lapply(docs_all[years == x], function(y) y[1,]))+1, levels = seq_len(length(vocab))))
    tmp = matrix(as.integer(tmp), nrow = K)
    colnames(tmp) = vocab
    tmp
  })
  
  topwords_chunks = lapply(topics_chunks, topWords, numWords = 10)
  topwords_chunks = lapply(1:K, function(k) sapply(seq_along(topwords_chunks), function(t) topwords_chunks[[t]][,k]))
  topwords_chunks = lapply(topwords_chunks, function(x){colnames(x) = sort(unique(years)); x})
  
  # rebuild keywords
  topwords_chunks <- lapply(topwords_chunks, rebuild_keywords)
  
  return(topwords_chunks)
  
}

meta <- readRDS("./input data/meta.RDS")
years_unique <- sort(unique(meta$PY))

sim_mean_names_by_year <- readRDS("./Model variants/RDS/sim_mean_names_by_year.RDS")
top_candidates <- names(head(sort(sim_mean_names_by_year, decreasing = TRUE), 5))



for (j in 1:length(top_candidates)){
  
  # docs_all from candidate model
  docs_all <- readRDS(paste0("./Model variants/", top_candidates[j],"/docs_all.RDS"))
  
  # complete candidate model
  lda <- readRDS(paste0("./Model variants/", top_candidates[j],"/lda.RDS"))
  
  topic_evo_by_year <- make_topic_evo_by_year(lda, docs_all, meta)
  
  saveRDS(topic_evo_by_year, file = paste0("./Model variants/", top_candidates[j],"/topic_evo_by_year.RDS"))
  
}

