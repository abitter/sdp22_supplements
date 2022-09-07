# get top terms from beta
get_top_terms <- function(x){
  terms <- names(head(sort(x, decreasing = TRUE), 10))
  return(paste(terms, collapse = ", "))
}

beta <- t(getEstimators(getLDA(res1))$phi)
topterms <- unlist(apply(beta, 2, get_top_terms))

head(topterms)



## Topic evolution without importance weighting ----

topics_new <- list()
for (i in 2:length(rolling_years)){
  temp <- readRDS(paste0("./Demo model 3/RDS_rolling/beta_", i, ".RDS"))
  topics_new[[i]] <- unlist(apply(t(temp), 2, function(x){names(head(sort(x, decreasing = TRUE), 20))}))
}

# function
get_evo <- function(t){
  topic_evo <- cbind(topics1[,t], topics_new[[2]][,t]) # topics_new starts with [[2]], see for loop below
  for (j in 3:length(topics_new)){
    topic_evo <- cbind(topic_evo, topics_new[[j]][,t])
  }
  colnames(topic_evo) <- rolling_years # rolling_years includes last year of init data
  return(topic_evo)
}

# list of topic evolutions
topic_evo <- list()
for (t in 1:K){
  topic_evo[[t]] <- get_evo(t)
}

topic_evo[[228]]
