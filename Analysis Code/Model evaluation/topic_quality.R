### Topic Quality


# The following only works for the init models!
# FOR COMPUTING QUALITY FOR EACH TOPIC ACROSS ALL YEARS
# SEE ./Model variants/R/5_variant quality


# Semantic coherence ------------------------------------------------------



# Mimno, D., Wallach, H. M., Talley, E., Leenders, M., & McCallum, A. (2011, July). 
# "Optimizing semantic coherence in topic models." In Proceedings of the Conference on Empirical Methods in 
# Natural Language Processing (pp. 262-272). Association for Computational Linguistics. Chicago

# coh <- lapply(candidates_2020, function(x) {tosca::topicCoherence(getLDA(x), docs)}) # this takes a while

# Error in wordid[j] %in% x[1, x[2, ] == i] : 
#   unused arguments (function () 
#   {
#   }, integer(0))


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





# Exclusivity -------------------------------------------------------------



## exclusivity

# using LDAvis relevance score with lambda = 0 for putting emphasis on exclusivity
# Sievert, C., Shirley, K.E.: LDAvis: 
# A method for visualizing and interpreting topics. 
# In: Proceedings of the workshop on interactive language learning, visualization, and interfaces. 
# pp. 63–70 (2014).


# original codes by A. Niekler, G. Wiedemann, and A. Fischer

exclusivity <- function(lda, dtm, lambda = 0, num.words = 10){
  
  dfm <- quanteda::as.dfm(dtm)
  
  if (num.words == 0) 
    num.words = dim(dfm)[2]
  
  pwt <- as.matrix(getEstimators(getLDA(lda))$phi)
  pw <- colSums(dfm)/sum(dfm)
  res <- apply(pwt, 1, function(x, num.words, pw, lambda) {
    x <- lambda * log(x) + (1 - lambda) * log(x/pw)
    return((sort(x, decreasing = TRUE)[1:num.words]))
  }, num.words, pw, lambda)
  return(colMeans(res))
}





# Inspect Quality ---------------------------------------------------------

library(ldaPrototype)

DFM <- quanteda::as.dfm(DTM)

coh <- topicCoherence(res1, DTM)
exc <- exclusivity(res1, DFM)

plot(coh, exc, col = "white")
K <- ldaPrototype::getLDA(res1)$param$K
text(coh, exc, labels = 1:K, cex = 0.7)

topics1 <- tosca::topWords(getTopics(getLDA(res1)), 10, byScore = FALSE)
topic_evo <- readRDS("./Demo model 3 keywords/RDS_psychtopics/topic_evo.RDS")

topic_evo[[7]]
