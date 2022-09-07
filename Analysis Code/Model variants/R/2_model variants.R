# Script 2: Compute model variants



# Data --------------------------------------------------------------------

meta <- readRDS("./input data/meta.RDS")
text <- readRDS("./input data/text.RDS")




# Settings ----------------------------------------------------------------

start_roll_seq <- seq(1990, 2015, 5)
K_seq <- seq(200, 500, 50)

alpha <- 0.0001
num.iterations <- 500
n_proto <- 25
min_docfreq <- 25
vocab.limit <- 10

years <- sort(unique(meta$PY))
current_year <- max(years)




# Libraries ---------------------------------------------------------------

library(quanteda)
library(tosca)
library(ldaPrototype)
library(ldaGibbs) # devtools::install_github("JonasRieger/ldaGibbs")
library(rcrossref)




# Stopwords and Lemma -----------------------------------------------------

stopwords_christ <- read.table("./ESM/Stopwords_Christ.csv")
stopwords_christ <- stopwords_christ$V1
stopwords_christ <- stopwords_christ[!is.na(stopwords_christ)]
stopwords_bittermann <- read.csv("./ESM/Stopwords_Bittermann.csv",
                                 header = FALSE)
stopwords_bittermann <- stopwords_bittermann$V1

stopwords_psyndex <- c("abstract", "copyright", "zpid", "released", "publisher", 
                       "translated", "deepl", "springer", "thieme","author", "original")



# Lemmatization

# https://github.com/michmech/lemmatization-lists/
lemma_en <- readr::read_delim("./ESM/lemmatization-en.txt", "\t", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE)

# baseform is in col 1
# switch cols
lemma_en <- lemma_en[,c(2,1)]
colnames(lemma_en) <- c("variant", "base")

# lower case
lemma_en <- as.data.frame(lemma_en) # due to some tidyverse tbls problem
lemma_en[,1] <- tolower(lemma_en[,1])
lemma_en[,2] <- tolower(lemma_en[,2])

# don't change "data" to "datum"
lemma_en <- lemma_en[-which(lemma_en[,1] == "data"),]

# don't change "sentencing" to "sentence"
lemma_en <- lemma_en[-which(lemma_en[,1] == "sentencing"),]

# remove duplicates
lemma_en <- lemma_en[!duplicated(lemma_en[,1]),]
rownames(lemma_en) <- NULL
rownames(lemma_en) <- 1:nrow(lemma_en)




# Beginn model variants double loop ---------------------------------------

dir.create("variants")

start <- Sys.time()

for (start_roll in start_roll_seq){ # start of RollingLDA

  for (K in K_seq){ # number of topics
    
    dir.create(paste0("./variants/variant_", K, "_", start_roll))
    dir.create(paste0("./variants/variant_", K, "_", start_roll, "/RDS_init"))
    dir.create(paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling"))
    dir.create(paste0("./variants/variant_", K, "_", start_roll, "/RDS_psychtopics"))
    
    
    
    # Initialization Data -----------------------------------------------------
    
    
    # use 1980-(start_roll-1 as initialization data
    # all subsequent years used for RollingLDAs
    # (start_roll-1 included, as for first rolling year the previous year is memory in the for loop
    rolling_years <- (start_roll-1):max(years)
    saveRDS(rolling_years, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_init/rolling_years.RDS"))
    
    
    keep <- which(meta$PY <= (start_roll-1))
    init_text <- text[keep]
    init_meta <- meta[keep,]
    
    
    
    # Preprocessing -----------------------------------------------------------
    
    
    # DFM
    corpus <- corpus(init_text)
    DFM <- corpus %>% 
      tokens() %>% 
      tokens_tolower() %>% 
      tokens_remove(pattern = c(stopwords("en"), stopwords("de"),
                                stopwords_christ, stopwords_bittermann, stopwords_psyndex)) %>%
      tokens(remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, remove_separators = TRUE) %>% 
      tokens_replace(., lemma_en$variant, lemma_en$base) %>%
      dfm() 
    
    DFM <- dfm_trim(DFM, min_docfreq = min_docfreq, min_termfreq = vocab.limit+1)
    
    
    # DTM for coherence calculation
    DTM <- convert(DFM, to = "topicmodels", omit_empty = TRUE)
    saveRDS(DTM, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_init/DTM.RDS"))
    
    # docs and vocab
    temp <- convert(DFM, to = "lda", omit_empty = TRUE)
    docs <- temp$documents
    vocab <- temp$vocab
    rm(temp)
    
    # set all freqs to 1, as expected by LDArep of ldaPrototype
    # https://github.com/JonasRieger/ldaPrototype/issues/10
    docs <- lapply(docs, function(x) rbind(rep(x[1,], x[2,]), 1L))
    
    # get indices of remaining docs
    # names(docs) matches rownames in init_text and init_meta when removing the "text" part.
    kept_docs <- as.numeric(gsub("text", "", names(docs)))
    
    # save meta of remaining docs
    meta_kept <- init_meta[kept_docs,] # will be master meta object
    
    
    
    
    # ldaPrototype ------------------------------------------------------------
    
    # https://github.com/JonasRieger/ldaPrototype
    
    start <- Sys.time()
    runif(1) # https://github.com/JonasRieger/ldaPrototype/issues/11
    res1 <- LDAPrototype(docs = docs, vocabLDA = vocab, K = K, alpha = alpha, n = n_proto, seeds = 1:n_proto, 
                         pm.backend = "multicore", ncpus = 8)
    proto1 <- getPrototype(res1) #= getLDA(res1)
    end <- Sys.time()
    end - start
    
    # save
    saveRDS(vocab, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_init/vocab.RDS"))
    saveRDS(docs, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_init/docs.RDS"))
    saveRDS(meta_kept, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_init/meta_kept.RDS"))
    saveRDS(res1, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_init/res1.RDS"))
    saveRDS(proto1, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_init/proto1.RDS"))
    
    # top terms
    topics1 <- tosca::topWords(getTopics(getLDA(res1)), 10)
    
    # theta
    theta1 <- t(getEstimators(getPrototype(res1))$theta)
    
    
    # topic assignments
    # NOTE: different indexing: assignment 0 is for topic 1
    assignments <- getAssignments(proto1)
    
    
    # init for RollingLDA
    init <- res1
    lda <- getLDA(init)
    assignments <- getAssignments(lda)
    docs_all <- docs
    names(docs_all) <- meta_kept$DFK # DFK is document identifier in PSYNDEX
    
    
    
    
    # RollingLDA --------------------------------------------------------------
    
    # https://github.com/JonasRieger/upi/blob/main/lda.R
    
    # 1. use last year (rolling_years[i-1]) as seed for modeling the topic assignments to docs from new year (rolling_years[i])
    # "The vocabulary is extended by words that occur more than five times in the new [...] articles 
    # and that were not included in the vocabulary before"
    
    
    topics_new <- list() # to save top terms of each new year
    
    for (i in 2:length(rolling_years)){ # starting with 2, as 1 is last year in init, needed for memory
      
      ## new data
      keep <- which(meta$PY == rolling_years[i])
      rolling_text <- text[keep]
      rolling_meta <- meta[keep,]
      
      # preprocessing
      DFM_rolling <- rolling_text %>%
        corpus() %>% 
        tokens() %>% 
        tokens_tolower() %>% 
        tokens_remove(pattern = c(stopwords("en"), stopwords("de"),
                                  stopwords_christ, stopwords_bittermann, stopwords_psyndex)) %>%
        tokens(remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, remove_separators = TRUE) %>% 
        tokens_replace(., lemma_en$variant, lemma_en$base) %>%
        dfm() %>%
        dfm_trim(min_docfreq = min_docfreq)
      
      # implementation of vocab.limit: keep words, that occur more often than vocab.limit-times, 
      # or that are already included in vocab
      DFM_rolling <- dfm_select(DFM_rolling,
                                colnames(DFM_rolling)[colSums(DFM_rolling) > vocab.limit | 
                                                        colnames(DFM_rolling) %in% vocab])
      
      # DTM for coherence calculation
      DTM_rolling <- convert(DFM_rolling, to = "topicmodels", omit_empty = TRUE)
      
      
      # convert
      temp <- convert(DFM_rolling, to = "lda", omit_empty = TRUE)
      docs_rolling <- temp$documents
      vocab_rolling <- temp$vocab
      rm(temp)
      
      # add new terms
      ind <- !(vocab_rolling %in% vocab)
      if(any(ind)){
        vocab <- c(vocab, vocab_rolling[ind])
      }
      
      # match vocab indices
      mtch <- match(vocab_rolling, vocab)-1
      # transform to integer (as expected by ldaGibbs)
      mtch <- as.integer(mtch)
      docs_rolling <- lapply(docs_rolling, function(x) rbind(mtch[x[1,]+1], x[2,]))
      
      # set all freqs to 1, as expected by LDArep of ldaPrototype
      # https://github.com/JonasRieger/ldaPrototype/issues/10
      docs_rolling <- lapply(docs_rolling, function(x) rbind(rep(x[1,], x[2,]), 1L))
      
      # get indices of remaining docs
      # names(docs_rolling) matches rownames in rolling_text and rolling_meta when removing the "text" part.
      kept_docs_rolling <- as.numeric(gsub("text", "", names(docs_rolling)))
      
      # meta
      rolling_meta <- rolling_meta[kept_docs_rolling,] # remove docs dropped
      meta_kept <- rbind(meta_kept, rolling_meta) # add to meta object with kept docs only
      
      # use DFK as names for docs
      names(docs_rolling) <- rolling_meta$DFK
      
      # add new docs and new doc names
      n.docs <- length(docs_all)
      docs_all[(n.docs+1):(n.docs+length(docs_rolling))] <- docs_rolling # add new docs to docs_all
      names(docs_all)[(n.docs+1):length(docs_all)] <- names(docs_rolling) # add new doc names
      
      # "The topic assignments of the new articles are initialized randomly" 
      assignments[(n.docs+1):length(docs_all)] <- lapply(
        sapply(docs_all[(n.docs+1):length(docs_all)], ncol),
        function(n) as.integer(sample(K, n, replace = TRUE)-1))
      
      
      ## modeling data ----
      
      # modeling data: last year (memory) + new year
      
      # get DFKs (PSYNDEX document identifiers) in memory data
      memory_docs <- meta_kept[meta_kept$PY == rolling_years[i-1], 1] # col 1 is DFK
      
      # which of the names in docs_all occur in this specified time window?
      ind_chunk <- names(docs_all) %in% c(memory_docs, names(docs_rolling))
      # table(ind_chunk)
      # ind_chunk == TRUE are all docs from memory + new year = modeling data
      
      # number of docs in "memory"
      # n.init = sum(obj$meta$date >= cache[i] & obj$meta$date < cuts[i]) # last three quarters
      n.init <- length(memory_docs)
      
      # docs_old <- docs # save old docs object
      docs <- docs_all[ind_chunk]
      
      topics <- matrix(as.integer(
        table(unlist(assignments[ind_chunk]) + 1, # Topics
              factor(unlist(lapply(docs, function(x) x[1,])) + 1, seq_len(length(vocab)), vocab)) # Vocab
      ), nrow = K)
      colnames(topics) <- vocab
      
      initial <- list(
        assignments = assignments[ind_chunk],
        topics = topics,
        topic_sums = matrix(as.integer(rowSums(topics))))
      
      
      ## LDA ----
      
      # "to model Latent Dirichlet Allocations with a subset of articles assignments fixed"
      
      # "[...] the Gibbs sampler iterates over each of the new articles again 200 times, 
      # while the topic assignments of all articles acting as initializing memory remain constant"
      res <- ldaGibbs::ldaGibbs(docs, K, vocab, num.iterations, alpha, eta = 1/K, initial = initial, n.init = n.init)
      
      # replace random assignments with new assignments from res
      assignments[(n.docs+1):length(docs_all)] <- res$assignments[(n.init+1):length(docs)]
      
      topics <- matrix(as.integer(
        table(unlist(assignments) + 1, # Topics
              factor(unlist(lapply(docs_all, function(x) x[1,])) + 1, seq_len(length(vocab)), vocab)) # Vocab
      ), nrow = K)
      colnames(topics) <- vocab
      document_sums <- cbind(getDocument_sums(lda), res$document_sums[,(n.init+1):length(docs)])
      
      lda <- LDA(assignments = assignments, topics = topics, document_sums = document_sums,
                 param = list(K = K, alpha = alpha, eta = 1/K, num.iterations = num.iterations))
      
      topics_new[[i]] <- tosca::topWords(getTopics(lda), 10) # importance weighting
      # no importance weighting: top terms by beta
      #topics_new[[i]] <- tosca::topWords(getTopics(lda), 10, byScore = FALSE)
      
      
      ## save lda, beta, DTM ----
      
      saveRDS(lda, 
              file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling/lda_", i, ".RDS"))
      
      beta <- getEstimators(lda)$phi
      saveRDS(beta,
              file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling/beta_", i, ".RDS"))
      
      saveRDS(DTM_rolling,
              file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling/DTM_", i, ".RDS"))
      
    }

    
    ## save objects ----
    
    # meta data of kept docs
    saveRDS(meta_kept, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling/meta_kept.RDS"))

    # list of top terms
    saveRDS(topics_new, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling/topics_new.RDS"))
    
    # docs
    saveRDS(docs_all, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling/docs_all.RDS"))
    
    

    
    
    # PsychTopics RDS ---------------------------------------------------------
    
    # objects needed:
    # topics_new
    # lda
    # meta_kept
    # years
    
    
    ## topic_evo ----
    
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
    
    
    # rebuild keywords
    rebuild_keywords <- function(x){
      x <- gsub("RRR1", "(", x, ignore.case = TRUE)
      x <- gsub("RRR2", ")", x, ignore.case = TRUE)
      x <- gsub("RRR3", "/", x, ignore.case = TRUE)
      x <- gsub("RRR4_", "'", x, ignore.case = TRUE)
      x <- gsub("_as_", " as ", x, ignore.case = TRUE)
      x <- gsub("_", " ", x)
      return(x)
    }
    
    topic_evo <- lapply(topic_evo, rebuild_keywords)
    
    saveRDS(topic_evo, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_psychtopics/topic_evo.RDS"))
    
    
    ## topic evo for psychtopics 2, by Zauad Shahreer
    
    make_topic_evo_concatenated = function(topic_evo_rds, directory = "./") {
      
      make_topic_evo_string = function(each) {
        
        each = as.data.frame(each)
        years = names(each)
        
        get_all_strings = function(year) {
          strings = glue::glue_collapse(each[[year]], sep = ", ")
          glue::glue("{year}: {strings}")
        }
        
        all_strings = sapply(years, get_all_strings)
        glue::glue_collapse(all_strings, sep = "\n")
        
      }
      
      topic_evo_concatenated = sapply(topic_evo, make_topic_evo_string)

      saveRDS(topic_evo_concatenated, 
              file = paste0("./variants/variant_", K, "_", start_roll, 
                            "/RDS_psychtopics/topic_evo_concatenated.RDS"))
    }
    
    make_topic_evo_concatenated(topic_evo)
    
    
    
    ## years ----
    
    saveRDS(years, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_psychtopics/years.RDS"))
    
    
    
    ## beta by year ----
    
    # get vocabulary size of last lda
    n_vocab <- ncol(lda$topics)
    terms_vocab <- colnames(lda$topics)
    
    # expand first beta to max vocab size
    beta1 <- getEstimators(getPrototype(res1))$phi
    temp_mat <- matrix(0, nrow = nrow(beta1), ncol = n_vocab - ncol(beta1))
    beta1 <- cbind(beta1, temp_mat)
    colnames(beta1) <- terms_vocab
    
    # initiate beta_list
    beta_list <- list()
    for (t in 1:K){
      beta_list[[t]] <- beta1[t,]
    }
    
    # add betas from other years
    for (i in 2:length(rolling_years)){ # 1 is init model
      
      # read data
      temp_beta <- readRDS(paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling/beta_", i, ".RDS"))
      
      # expand
      temp_mat <- matrix(0, nrow = nrow(temp_beta), ncol = n_vocab - ncol(temp_beta))
      temp_beta <- cbind(temp_beta, temp_mat)
      colnames(temp_beta) <- terms_vocab
      
      for (t in 1:K){
        beta_list[[t]] <- rbind(beta_list[[t]], temp_beta[t,])
      }
      
    }
    
    saveRDS(beta_list, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_rolling/beta_list.RDS"))
    
    # mean beta probabilities
    beta_means <- lapply(beta_list, colMeans)
    
    
    
    ## booster ----
    
    get_beta_booster <- function(x){
      betas <- unname(head(sort(x, decreasing = TRUE), 10))
      
      factor_list <- list()
      for (i in 1:(length(betas)-1)){
        factor_list[[i]] <- betas[i]/betas[10] # booster factors
      }
      return(round(unlist(factor_list), 2))
    }
    
    booster_list <- lapply(beta_means, get_beta_booster)
    
    booster <- as.data.frame(booster_list)
    rownames(booster) <- NULL
    colnames(booster) <- 1:K
    
    saveRDS(booster, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_psychtopics/booster.RDS"))
    
    
    
    ## topic evo search link ----
    
    createLink_evo <- function(val, boost) {
      list <- list()
      
      # use colMeans of booster for evo terms
      # Future update: compute separate booster objects for each year
      booster_means <- round(rowMeans(booster), 2)
      
      for (i in 1:length(val)){
        list[[i]] <- unlist(strsplit(val[i], ", ", fixed = TRUE))
        for (j in 1:9){
          list[[i]][j] <- paste0('"', list[[i]][j], '"%5E', booster_means[j]) # add boost factors for first 9 terms
        }
        list[[i]][10] <- paste0('"', list[[i]][10], '"') # Term 10 is reference, so no boosting
        list[[i]] <- paste0(list[[i]], collapse="+OR+")
        list[[i]] <- gsub("'", "%27", list[[i]])
      }
      val <- unlist(list)
      paste0("<a href='https://pubpsych.zpid.de/pubpsych/Search.action?q=%28%28", 
             val,"%29%29+DB%3DPSYNDEX&stats=TOP' target='_blank' class='btn btn-primary'>Search PSYNDEX</a>")
    }
    
    
    topic_evo_search <- list()
    for (i in 1: length(topic_evo)){
      val <- apply(topic_evo[[i]], 2, paste, collapse = ", ")
      link <- createLink_evo(val, booster)
      topic_evo_search[[i]] <- rbind(topic_evo[[i]], link)
    }
    
    saveRDS(topic_evo_search, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_psychtopics/topic_evo_search.RDS"))
    
    
    
    
    ## topic ----
    
    ## Top Terms
    # get top 10 terms across all years using mean beta
    get_top_terms <- function(x){
      terms <- names(head(sort(x, decreasing = TRUE), 10))
      return(paste(terms, collapse = ", "))
    }
    
    TopTerms <- unlist(lapply(beta_means, get_top_terms))
    TopTerms <- rebuild_keywords(TopTerms)
    
    
    ## n_docs
    
    # get theta from last lda (all thetas included in last lda)
    theta <- t(getEstimators(lda)$theta)
    
    # number of docs with theta > .5 per topic
    n_docs <- apply(theta, 2, function(x){unname(table(x > 0.5)[2])})
    
    
    ## Labels
    
    # for testing purposes, use first two top terms as Labels
    Label <- sapply(TopTerms, function(x) {paste(strsplit(x, ", ")[[1]][1:2], collapse = " ")})
    Label <- unname(Label)
    # add ID for disambiguation
    for (i in 1:K){
      Label[i] <- paste0("T", i, ": ", Label[i])
    }
    
    
    
    ## Empirical
    get_emp <- function(x){
      indices <- which(x > .5)
      emp_table <- table(meta_kept$EMP[indices])
      unname(round((emp_table[2] / sum(emp_table)*100), 2))
    }
    Empirical <- apply(theta, 2, get_emp)
    
    
    ## Journals
    # some JT fields contain ISSN, not journal titles
    
    # replacement table
    indices <- which(!duplicated(meta_kept$ISSN))
    issn_jt <- meta_kept[indices,6:5]
    issn_jt <- issn_jt[issn_jt$ISSN != "",]
    
    # use rcrossref to query issn not included in issn_jt
    issn_query <- function(x){
      if(!grepl("\\D", substr(x, 1, 4))){
        res <- rcrossref::cr_journals(issn = x)$data$title
        if (!is.null(res)){
          return(res)
        } else {return(x)}
      } else {return(x)}
    }
    
    get_journals <- function(x){
      indices <- which(x > .5)
      temp <-  meta_kept$JT[indices]
      temp <- gsub(",", "", temp)
      temp <- gsub("\\.", "", temp)
      temp <- gsub("&", " ", temp)
      temp <- gsub(" ", "_", temp)
      
      
      # replace ISSN with journal title
      temp <- unname(unlist(tokens_replace(tokens(temp), issn_jt$ISSN, issn_jt$JT)))
      #temp <- unname(sapply(temp, issn_query))
      
      # get top 3 journals
      temp <- gsub(",", "", temp)
      temp <- gsub("\\.", "", temp)
      temp <- gsub("&", " ", temp)
      temp <- gsub(" ", "_", temp)
      temp <- names(topfeatures(dfm(tokens(temp), tolower = FALSE), 3))
      temp <- gsub("_", " ", temp)
      
      return(paste(temp, collapse = " | "))
    }
    
    Journals <- apply(theta, 2, get_journals)
    
    
    ## data frame
    topic <- data.frame("ID" = 1:K, Label, TopTerms, n_docs, Empirical, Journals)
    saveRDS(topic, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_psychtopics/topic.RDS"))
    
    
    
    ## n_docs_year ----
    
    # make factor in order to include years with freq = 0 in tables
    PY_factor <- as.factor(meta_kept$PY)
    
    # first topic
    n_docs_year <- unname(table(PY_factor[which(theta[,1] > .5)]))
    
    # remaining topics
    for (i in 2:K){
      temp <- unname(table(PY_factor[which(theta[,i] > .5)]))
      n_docs_year <- cbind(n_docs_year, temp)
    }
    
    colnames(n_docs_year) <- Label
    rownames(n_docs_year) <- years
    
    
    saveRDS(n_docs_year, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_psychtopics/n_docs_year.RDS"))
    
    
    
    ## empirical_year ----
    
    
    empirical_list <- list()
    
    for (i in 1:K){
      
      # subset of meta_kept for topic i
      ind <- which(theta[,i] > .5)
      
      # within subset, share of empirical studies by year
      empirical_year <- by(meta_kept[ind,], PY_factor[ind], function(x){
        # table(x$EMP)
        temp <- table(x$EMP)
        val <- unname(temp[names(temp) == "1"])
        
        if (length(val) == 0){
          return(0)
        } else {
          res <- val/(sum(temp))
          return(res)
        }
        
      })
      
      empirical_list[[i]] <- as.vector(empirical_year)
      
    }
    
    empirical_year <- as.matrix(as.data.frame(empirical_list))
    empirical_year[is.na(empirical_year)] <- 0
    empirical_year <- round(empirical_year*100, 2)
    
    colnames(empirical_year) <- Label
    rownames(empirical_year) <- years
    
    
    saveRDS(empirical_year, 
            file = paste0("./variants/variant_", K, "_", start_roll, "/RDS_psychtopics/empirical_year.RDS"))
    
    
  }
}

end <- Sys.time()
end-start
