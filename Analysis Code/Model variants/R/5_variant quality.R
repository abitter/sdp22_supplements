# Inspect topic quality of candidate models for each year separately


library(tosca)
library(ldaPrototype)
library(tm)
library(quanteda)



# get files ---------------------------------------------------------------

meta <- readRDS("./input data/meta.RDS") # meta dataframe

years <- meta$PY
# years <- meta$PY[match(names(docs_all), meta$DFK)]
years_unique <- sort(unique(years))


# get candidate models with highest similarity to 2020 ldaprototype
sim_mean_names_by_year <- readRDS("./Model variants/RDS/sim_mean_names_by_year.RDS")
top_candidates <- names(head(sort(sim_mean_names_by_year, decreasing = TRUE), 5))




# dtm helper function -----------------------------------------------------

create_dtm = function(lda, docs_all, years, y = unique(years), to_dfm = FALSE){ # y is filter for years
  docs = docs_all[years %in% y]
  years = years[years %in% y]
  docs = lapply(docs, function(x) table(x[1,]))
  
  i = rep(seq_along(docs), lengths(docs))
  j = as.integer(unlist(lapply(docs, names))) + 1L
  v = unname(unlist(lapply(docs, as.numeric)))
  vocab = colnames(getTopics(lda))[seq_len(max(j))]
  dtm = slam::simple_triplet_matrix(i, j, v,
                                    dimnames = list(Docs = paste0("text", seq_along(docs)), 
                                                    Terms = vocab))
  dtm = as.DocumentTermMatrix(dtm, weighting = weightTf)
  
  # remove empty cols (i.e., terms not occurring in given year)
  dfm <- quanteda::as.dfm(dtm)
  dfm <- dfm[,colSums(dfm) != 0]
  
  # return dfm or dtm?
  if (to_dfm == TRUE){
    return(dfm)
  }
  else {return(quanteda::convert(dfm, to = "tm"))}
  
}

# a = create_dtm(lda, docs_all, years, 2000) # DTM for year 2000
# b = create_dtm(lda, docs_all, years, 2000:2005) # DTM for years 2000, 2001, 2002, 2003, 2004, 2005




# topic coherence by year -------------------------------------------------

# Mimno, D., Wallach, H. M., Talley, E., Leenders, M., & McCallum, A. (2011, July). 
# "Optimizing semantic coherence in topic models." In Proceedings of the Conference on Empirical Methods in 
# Natural Language Processing (pp. 262-272). Association for Computational Linguistics. Chicago

# modified from tmca_coherence() function (A. Niekler & G. Wiedemann)

# topicCoherence_by_year = function(lda, docs_all, years, y, N = 10){
coherence_by_year <- function(model, DTM, years, y, N = 10) {

  K <- ldaPrototype::getK(model)
  
  # topic matrix
  assignments = getAssignments(model)
  vocab = colnames(getTopics(model))
  
  tmp = table(factor(unlist(assignments[years %in% y])+1, levels = 1:K), 
              factor(unlist(lapply(docs_all[years %in% y], function(x) x[1,]))+1,
                     levels = seq_len(length(vocab))))
  topics_by_year = matrix(as.integer(tmp), nrow = K)
  colnames(topics_by_year) = vocab
  
  DTM = create_dtm(lda, docs_all, years, y)
  
  # Ensure matrix or Matrix-format (convert if slam)
  require(Matrix)
  require(slam)
  if (is.simple_triplet_matrix(DTM)) {
    DTM <- sparseMatrix(i=DTM$i, j=DTM$j, x=DTM$v, dims=c(DTM$nrow, DTM$ncol), dimnames = dimnames(DTM))
  }
  
  topNtermsPerTopic <- tosca::topWords(topics_by_year, N)
  allTopicModelTerms <- unique(as.vector(topNtermsPerTopic))
  
  DTMBIN <- DTM > 0
  
  documentFrequency <- colSums(DTMBIN)
  names(documentFrequency) <- colnames(DTMBIN)
  
  topNtermsPerTopic <- tosca::topWords(topics_by_year, N)
  allTopicModelTerms <- unique(as.vector(topNtermsPerTopic))
  
  DTMBIN <- DTMBIN[, allTopicModelTerms]
  DTMBINCooc <- t(DTMBIN) %*% DTMBIN
  DTMBINCooc <- t((DTMBINCooc + 1) / colSums(DTMBIN))
  DTMBINCooc <- log(DTMBINCooc)
  DTMBINCooc <- as.matrix(DTMBINCooc)
  
  coherence <- rep(0, K)
  pb <- txtProgressBar(max = K)
  for (topicIdx in 1:K) {
    setTxtProgressBar(pb, topicIdx)
    topWordsOfTopic <- topNtermsPerTopic[,topicIdx]
    
    coherence[topicIdx] <- 0
    for (m in 2:length(topWordsOfTopic)) {
      for (l in 1:(m-1)) {
        mTerm <- as.character(topWordsOfTopic[m])
        lTerm <- as.character(topWordsOfTopic[l])
        coherence[topicIdx] <- coherence[topicIdx] + DTMBINCooc[mTerm, lTerm]
      }
    }
  }
  close(pb)
  
  return(coherence)
}



# topic exclusivity by year -----------------------------------------------

# using LDAvis relevance score with lambda = 0 for putting emphasis on exclusivity
# Sievert, C., Shirley, K.E.: LDAvis: 
# A method for visualizing and interpreting topics. 
# In: Proceedings of the workshop on interactive language learning, visualization, and interfaces. 
# pp. 63–70 (2014).


# original codes by A. Niekler, G. Wiedemann, and A. Fischer

exclusivity_by_year <- function(lda, dtm, lambda = 0, num.words = 10){
  
  dfm <- quanteda::as.dfm(dtm)
  
  if (num.words == 0) 
    num.words = dim(dfm)[2]
  
  pwt <- as.matrix(getEstimators(lda)$phi)
  # keep only terms of selected year(s) y
  keep <- colnames(pwt) %in% colnames(dfm)
  pwt <- pwt[,keep]
  
  pw <- colSums(dfm)/sum(dfm)
  res <- apply(pwt, 1, function(x, num.words, pw, lambda) {
    x <- lambda * log(x) + (1 - lambda) * log(x/pw)
    return((sort(x, decreasing = TRUE)[1:num.words]))
  }, num.words, pw, lambda)
  return(colMeans(res))
}



# get quality -------------------------------------------------------------

candidates_coh <- list()
candidates_exc <- list()


for (j in 1:length(top_candidates)){
  
  # docs_all from candidate model
  docs_all <- readRDS(paste0("./Model variants/", top_candidates[j],"/docs_all.RDS"))
  
  # complete candidate model
  lda <- readRDS(paste0("./Model variants/", top_candidates[j],"/lda.RDS"))
  
  
  coh <- list()
  exc <- list()
  
  for (i in 1:(length(years_unique)-1)){
    dtm_i <- create_dtm(lda, docs_all, years, y = years_unique[i])
    coh[[i]]<- mean(coherence_by_year(lda, dtm_i, years, y = years_unique[i]))
    exc[[i]] <- mean(exclusivity_by_year(lda, dtm_i))
  }
  
  candidates_coh[[j]]<- unlist(coh)
  candidates_exc[[j]]<- unlist(exc)
}

saveRDS(candidates_coh, file = "./Model variants/RDS/candidates_coh.RDS")
saveRDS(candidates_exc, file = "./Model variants/RDS/candidates_exc.RDS")



# plots

cols <- c("blue", "green", "orange", "red", "black")

png("./Model variants/coherence_by_year.png")
plot(candidates_coh[[1]], type = "l", col = cols[1],
     main = "Semantic Coherence", xaxt = "n", ylab = "", xlab = "Year",
     ylim = c(min(unlist(candidates_coh)), max(unlist(candidates_coh))))
axis(1, at = 1:length(years_unique), labels = min(years_unique):max(years_unique))
for (i in 2:length(top_candidates)){
  lines(candidates_coh[[i]], col = cols[i])
}
legend("bottomright", legend = top_candidates, col = cols, lwd = 2)
dev.off()


png("./Model variants/exclusivity_by_year.png")
plot(candidates_exc[[1]], type = "l", col = cols[1],
     main = "Exclusivity", xaxt = "n", ylab = "", xlab = "Year",
     ylim = c(min(unlist(candidates_exc)), max(unlist(candidates_exc))))
axis(1, at = 1:length(years_unique), labels = min(years_unique):max(years_unique))
for (i in 2:length(top_candidates)){
  lines(candidates_exc[[i]], col = cols[i])
}
legend("bottomleft", legend = top_candidates, col = cols, lwd = 2)
dev.off()



## table with mean similarity,  mean coherence, and mean exclusivity ----

candidate_metrics <- data.frame(variant = top_candidates,
                                similarity = as.numeric(head(sort(sim_mean_names_by_year, decreasing = TRUE), 5)),
                                coherence = unlist(lapply(candidates_coh, mean)),
                                exclusivity = unlist(lapply(candidates_exc, mean)))
candidate_metrics$scaled_mean <- rowMeans(cbind(scale(candidate_metrics$similarity),
                                      scale(candidate_metrics$coherence),
                                      scale(candidate_metrics$exclusivity)))


write.csv(candidate_metrics, file = "./Model variants/candidate_metrics.csv", row.names = FALSE)


## Topic Coherence across all years, each topic separately ----

# modified from tmca_coherence() function (A. Niekler & G. Wiedemann)
topicCoherence <- function(model, DTM, N = 10) {
  
  # model <- ldaPrototype::getLDA(model)
  
  # Ensure matrix or Matrix-format (convert if slam)
  require(Matrix)
  require(slam)
  if (is.simple_triplet_matrix(DTM)) {
    DTM <- sparseMatrix(i=DTM$i, j=DTM$j, x=DTM$v, dims=c(DTM$nrow, DTM$ncol), dimnames = dimnames(DTM))
  }
  
  K <- ldaPrototype::getK(model)
  
  DTMBIN <- DTM > 0
  
  documentFrequency <- colSums(DTMBIN)
  names(documentFrequency) <- colnames(DTMBIN)
  
  topNtermsPerTopic <- tosca::topWords(getTopics(model), 10)
  allTopicModelTerms <- unique(as.vector(topNtermsPerTopic))
  
  DTMBIN <- DTMBIN[, allTopicModelTerms]
  DTMBINCooc <- t(DTMBIN) %*% DTMBIN
  DTMBINCooc <- t((DTMBINCooc + 1) / colSums(DTMBIN))
  DTMBINCooc <- log(DTMBINCooc)
  DTMBINCooc <- as.matrix(DTMBINCooc)
  
  coherence <- rep(0, K)
  pb <- txtProgressBar(max = K)
  for (topicIdx in 1:K) {
    setTxtProgressBar(pb, topicIdx)
    topWordsOfTopic <- topNtermsPerTopic[,topicIdx]
    
    coherence[topicIdx] <- 0
    for (m in 2:length(topWordsOfTopic)) {
      for (l in 1:(m-1)) {
        mTerm <- as.character(topWordsOfTopic[m])
        lTerm <- as.character(topWordsOfTopic[l])
        coherence[topicIdx] <- coherence[topicIdx] + DTMBINCooc[mTerm, lTerm]
      }
    }
  }
  close(pb)
  
  return(coherence)
}


### apply to a specific model ----

lda <- readRDS("./Model variants/variant_200_2010/lda.RDS")
docs_all <- readRDS("./Model variants/variant_200_2010/docs_all.RDS")

DTM <- create_dtm(lda, docs_all, years, 1980:2021)
coh <- topicCoherence(lda, DTM)

# add to topic table
topic <- readRDS("./Model variants/variant_200_2010/topic.RDS")
topic$coherence <- coh
