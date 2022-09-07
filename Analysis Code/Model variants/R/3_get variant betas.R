# get 2020 betas of each variant


library(rollinglda)
library(tosca)

meta <- readRDS("./input data/meta.RDS")

start_roll_seq <- seq(1990, 2015, 5)
K_seq <- seq(200, 500, 50)

variant_names <- "start"
beta_2020_candidates <- list()
beta_2020_candidates_by_year <- list()

for (start_roll in start_roll_seq){ # start of RollingLDA
  
  for (K in K_seq){ # number of topics
    
    variant_names <- c(variant_names, paste0("variant_", K, "_", start_roll))
    
    rolling_years_temp <- readRDS(paste0("./variants/variant_", K, "_", start_roll, "/RDS_init/rolling_years.RDS"))
    index_of_2020 <- which(rolling_years_temp == 2020)
    
    ## rolling beta
    index_of_list <- length(beta_2020_candidates)+1
    beta_2020_candidates[[index_of_list]] <- readRDS(paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling/beta_", index_of_2020, ".RDS"))
    
    
    ## beta by year
    roll <- readRDS(paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling/lda_", index_of_2020, ".RDS"))
    docs <- readRDS(paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling/docs_all.RDS"))
    
    assignments = getAssignments(roll)
    vocab = colnames(getTopics(roll))
    years = meta$PY[match(names(docs), meta$DFK)]
    
    # topics for 2020 (from buzz_search.R)
    topics_chunks = lapply(2020, function(x){
      tmp = table(factor(unlist(assignments[years == x])+1, levels = 1:K), 
                  factor(unlist(lapply(docs[years == x], function(y) y[1,]))+1, levels = seq_len(length(vocab))))
      tmp = matrix(as.integer(tmp), nrow = K)
      colnames(tmp) = vocab
      tmp
    })
    
    topics_chunks = topics_chunks[[1]]
    eta = getEta(roll) # = 1/K
    phi = (topics_chunks + eta)/(rowSums(topics_chunks) + ncol(topics_chunks) * eta)
    
    index_of_list <- length(beta_2020_candidates_by_year)+1
    beta_2020_candidates_by_year[[index_of_list]] <- phi
  }
}

variant_names <- variant_names[-1]

saveRDS(variant_names, file = "variant_names.RDS")
saveRDS(beta_2020_candidates, file = "beta_2020_candidates.RDS")
saveRDS(beta_2020_candidates_by_year, file = "beta_2020_candidates_by_year.RDS")