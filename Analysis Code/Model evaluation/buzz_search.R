library(rollinglda)
library(tosca)

meta = readRDS("./input data/meta.RDS")

roll = readRDS("./Demo model 3 keywords/RDS_rolling/lda.RDS")
assignments = getAssignments(roll)
docs = readRDS("./Demo model 3 keywords/RDS_rolling/docs_all.RDS")
vocab = colnames(getTopics(roll))
years = meta$PY[match(names(docs), meta$DFK)]
K = getK(roll)

topics_chunks = lapply(sort(unique(years)), function(x){
  tmp = table(factor(unlist(assignments[years == x])+1, levels = 1:K), 
              factor(unlist(lapply(docs[years == x], function(y) y[1,]))+1, levels = seq_len(length(vocab))))
  tmp = matrix(as.integer(tmp), nrow = K)
  colnames(tmp) = vocab
  tmp
})

topwords_chunks = lapply(topics_chunks, topWords, numWords = 50)
topwords_chunks = lapply(1:K, function(k) sapply(seq_along(topwords_chunks), function(t) topwords_chunks[[t]][,k]))
topwords_chunks = lapply(topwords_chunks, function(x){colnames(x) = sort(unique(years)); x})


terms_climate = c("temperature", "temperature_effects", "weather", "warming",
                  vocab[which(grepl("clima", vocab))])
terms_open = c("open", "openly", "reproducibility",
               vocab[which(grepl("replication", vocab))])

search_buzz = function(terms){
  buzzwords = list()
  for(i in seq_along(topwords_chunks)){
    x = topwords_chunks[[i]]
    buzz = unlist(lapply(sort(unique(years)), function(y){
      cand = x[1:10, as.character(y)]
      buzz = cand[cand %in% terms]
      pos = match(buzz, cand)
      if(length(buzz) > 0) paste0(y, "|", pos, ":", buzz)
    }))
    if(length(buzz) > 0) buzzwords[[as.character(i)]] = buzz
  }
  buzzwords
}

search_buzz(terms_open)

