# Compute cosine similarities
# between each variant's 2020 topics and the 2020 ldaprototype model


# Get files ---------------------------------------------------------------


# see ldaPrototype_2020.R
beta_2020 <- readRDS("./Candidates 2020/beta_2020.RDS")


# files from "3_get variant beta.R"
variant_names <- readRDS("./Model comparisons/variant_names.RDS")
beta_2020_candidates <- readRDS("./Model comparisons/beta_2020_candidates.RDS")
beta_2020_candidates_by_year <- readRDS("./Model comparisons/beta_2020_candidates_by_year.RDS")






# Cosine similarities using RollingLDA beta -------------------------------


# for each candidate model, save the vector of max cosine similarity with the proto model
# Each vector has length = number of topics in the proto model
sim_vectors <- list()
sim_mean <- list()
sim_vectors_by_year <- list()
sim_mean_by_year <- list()




# rolling beta similarities -----------------------------------------------



for (j in 1:length(beta_2020_candidates)){
  
  ## rolling beta ##
  
  
  
  ## prepare beta matrices ----
  
  beta_proto <- beta_2020
  beta_candidate <- beta_2020_candidates[[j]]
  
  
  # topics from different models have different vocabularies
  # for computing cosine similarity, match vocabularies in beta matrices
  vocab_proto <- colnames(beta_2020)
  vocab_candidate <- colnames(beta_candidate)
  vocab_both <- unique(c(vocab_proto, vocab_candidate))
  
  # table(vocab_both %in% vocab_proto)
  # table(vocab_both %in% vocab_candidate)
  
  
  # add cols for terms that are not included and set their probabilities to 0
  
  # proto
  temp_mat_proto <- matrix(0, nrow = nrow(beta_proto),
                           ncol = length(vocab_both) - length(vocab_proto))
  colnames(temp_mat_proto) <- vocab_both[!(vocab_both %in% vocab_proto)]
  beta_proto <- cbind(beta_proto, temp_mat_proto)
  
  # candidate
  temp_mat_candidate <- matrix(0, nrow = nrow(beta_candidate),
                               ncol = length(vocab_both) - length(vocab_candidate))
  colnames(temp_mat_candidate) <- vocab_both[!(vocab_both %in% vocab_candidate)]
  beta_candidate <- cbind(beta_candidate, temp_mat_candidate)
  
  
  # sort cols
  beta_proto <- beta_proto[, order(colnames(beta_proto))]
  beta_candidate <- beta_candidate[, order(colnames(beta_candidate))]
  
  
  
  ## cosine similarity ----
  
  max_sim <- list() # saves the max cosine similarity for each proto topic i
  
  for (i in 1:nrow(beta_proto)){
    cos_sim <- apply(beta_candidate, 1, function(x){
      as.numeric(lsa::cosine(beta_proto[i,], x))})
    max_sim[[i]] <- max(cos_sim)
    names(max_sim[[i]]) <- which.max(cos_sim)
  }
  
  sim_vectors[[j]] <- unlist(max_sim)
  sim_mean[[j]] <- mean(unlist(max_sim))
  
}

sim_mean_names <- unlist(sim_mean)
names(sim_mean_names) <- variant_names

saveRDS(sim_vectors, file = "./Model comparisons/sim_vectors.RDS")
saveRDS(sim_mean_names, file = "./Model comparisons/sim_mean_names.RDS")



# beta by year similarities -----------------------------------------------

for (j in 1:length(beta_2020_candidates_by_year)){
  
  
  ## prepare beta matrices ----
  
  beta_proto <- beta_2020
  beta_candidate <- beta_2020_candidates_by_year[[j]]
  
  
  # topics from different models have different vocabularies
  # for computing cosine similarity, match vocabularies in beta matrices
  vocab_proto <- colnames(beta_2020)
  vocab_candidate <- colnames(beta_candidate)
  vocab_both <- unique(c(vocab_proto, vocab_candidate))
  
  # table(vocab_both %in% vocab_proto)
  # table(vocab_both %in% vocab_candidate)
  
  
  # add cols for terms that are not included and set their probabilities to 0
  
  # proto
  temp_mat_proto <- matrix(0, nrow = nrow(beta_proto),
                           ncol = length(vocab_both) - length(vocab_proto))
  colnames(temp_mat_proto) <- vocab_both[!(vocab_both %in% vocab_proto)]
  beta_proto <- cbind(beta_proto, temp_mat_proto)
  
  # candidate
  temp_mat_candidate <- matrix(0, nrow = nrow(beta_candidate),
                               ncol = length(vocab_both) - length(vocab_candidate))
  colnames(temp_mat_candidate) <- vocab_both[!(vocab_both %in% vocab_candidate)]
  beta_candidate <- cbind(beta_candidate, temp_mat_candidate)
  
  
  # sort cols
  beta_proto <- beta_proto[, order(colnames(beta_proto))]
  beta_candidate <- beta_candidate[, order(colnames(beta_candidate))]
  
  
  
  ## cosine similarity ----
  
  max_sim <- list() # saves the max cosine similarity for each proto topic i
  
  for (i in 1:nrow(beta_proto)){
    cos_sim <- apply(beta_candidate, 1, function(x){
      as.numeric(lsa::cosine(beta_proto[i,], x))})
    max_sim[[i]] <- max(cos_sim)
    names(max_sim[[i]]) <- which.max(cos_sim)
  }
  
  sim_vectors_by_year[[j]] <- unlist(max_sim)
  sim_mean_by_year[[j]] <- mean(unlist(max_sim))
  
}

sim_mean_names_by_year <- unlist(sim_mean_by_year)
names(sim_mean_names_by_year) <- variant_names

saveRDS(sim_vectors_by_year, file = "./Model comparisons/sim_vectors_by_year.RDS")
saveRDS(sim_mean_names_by_year, file = "./Model comparisons/sim_mean_names_by_year.RDS")




# Final model -------------------------------------------------------------

# get model with max mean similarity
which.max(sim_mean_names)
which.max(sim_mean_names_by_year)

# top n models
head(sort(sim_mean_names, decreasing = TRUE), 10)
head(sort(sim_mean_names_by_year, decreasing = TRUE), 10)
