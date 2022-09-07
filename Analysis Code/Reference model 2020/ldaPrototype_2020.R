###
# Reference model 2020
###


# Libraries ---------------------------------------------------------------

library(quanteda)
library(tosca)
library(ldaPrototype)



# Load Data ---------------------------------------------------------------

meta <- readRDS("./input data/meta.RDS")
text <- readRDS("./input data/text.RDS")

keep <- which(meta$PY == 2020)
text_2020 <- text[keep]
meta_2020 <- meta[keep,]


# Preprocessing -----------------------------------------------------------

# stopwords for scientific abstracts

# Christ et al. (2019): http://dx.doi.org/10.23668/psycharchives.2613
# Bittermann & Klos (2019): http://dx.doi.org/10.23668/psycharchives.2499
# German stopwords by Bittermann & Klos were translated using DeepL (https://www.deepl.com/translator)

stopwords_christ <- read.table("./stopwords/Stopwords_Christ.csv")
stopwords_christ <- stopwords_christ$V1
stopwords_christ <- stopwords_christ[!is.na(stopwords_christ)]
stopwords_bittermann <- read.csv("./stopwords/Stopwords_Bittermann.csv",
                                 header = FALSE)
stopwords_bittermann <- stopwords_bittermann$V1

stopwords_psyndex <- c("abstract", "copyright", "zpid", "released", "publisher", 
                       "translated", "deepl", "springer", "thieme", "author", "original")

# Lemmatization
lemma_en <- readRDS("./input data/lemma_en.RDS")

corpus <- corpus(text_2020)

DFM <- corpus %>% 
  tokens() %>% 
  tokens_tolower() %>% 
  tokens_remove(pattern = c(stopwords("en"), stopwords("de"),
                            stopwords_christ, stopwords_bittermann, stopwords_psyndex)) %>%
  tokens(remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, remove_separators = TRUE) %>%
  tokens_replace(., lemma_en$variant, lemma_en$base) %>%
  dfm() 

topfeatures(DFM)
dim(DFM)


DFM <- dfm_trim(DFM, min_docfreq = 15)
dim(DFM)
saveRDS(DFM, file = "./Candidates 2020/DFM.RDS")

# omit empty doc in DFM
empty <- as.numeric(which(rowSums(DFM) == 0))
DFM_omit_empty <- DFM[-empty,]
saveRDS(DFM, file = "./Candidates 2020/DFM.RDS")
saveRDS(DFM_omit_empty, file = "./Candidates 2020/DFM_omit_empty.RDS")


# DTM for coherence function
DTM <- convert(DFM, to = "topicmodels", omit_empty = TRUE)
saveRDS(DTM, file = "./Candidates 2020/DTM.RDS")


temp <- convert(DFM, to = "lda", omit_empty = TRUE)
docs <- temp$documents
vocab <- temp$vocab
rm(temp)

# set all freqs to 1, as expected by LDArep of ldaPrototype
# https://github.com/JonasRieger/ldaPrototype/issues/10
docs <- lapply(docs, function(x) rbind(rep(x[1,], x[2,]), 1L))

# get indices of remaining docs
# names(docs) matches rownames in text_2020 and meta_2020 when removing the "text" part.
kept_docs <- as.numeric(gsub("text", "", names(docs)))

# save meta of remaining docs
meta_kept <- meta_2020[kept_docs,] # will be master meta object



# Optimal K ---------------------------------------------------------------

## Settings ----

K_range <- seq(150, 300, 25)
alpha <- 0.0001
num.iterations <- 500
n_proto <- 25


## Run Candidate Models ----

candidates_2020 <- list()
start <- Sys.time()
for (i in 1:length(K_range)){
  runif(1)
  candidates_2020[[i]] <- LDAPrototype(docs = docs, vocabLDA = vocab, 
                                  K = K_range[i], alpha = alpha, n = n_proto, 
                                  seeds = 1:n_proto, 
                                  pm.backend = "multicore", ncpus = 8)
}
end <- Sys.time()
end - start

saveRDS(candidates_2020, file = "./RDS_keywords/candidates_2020_keywords.RDS")




## evaluation metrics ----

## semantic coherence

# Mimno, D., Wallach, H. M., Talley, E., Leenders, M., & McCallum, A. (2011, July). 
# "Optimizing semantic coherence in topic models." In Proceedings of the Conference on Empirical Methods in 
# Natural Language Processing (pp. 262-272). Association for Computational Linguistics. Chicago


# modified from tmca_coherence() function (A. Niekler & G. Wiedemann)
topicCoherence <- function(model, DTM, N = 10) {
  
  model <- ldaPrototype::getLDA(model)
  
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

coh <- lapply(candidates_2020, topicCoherence, DTM)
coh_mean <- unlist(lapply(coh, mean))



## exclusivity

# using LDAvis relevance score with lambda = 0 for putting emphasis on exclusivity
# Sievert, C., Shirley, K.E.: LDAvis: 
# A method for visualizing and interpreting topics. 
# In: Proceedings of the workshop on interactive language learning, visualization, and interfaces. 
# pp. 63–70 (2014).


# original codes by A. Niekler, G. Wiedemann, and A. Fischer

exclusivity <- function(lda, dfm, lambda = 0, num.words = 0){
  library(ldaPrototype)
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

exc <- lapply(candidates_2020, exclusivity, DFM_omit_empty, lambda = 0, num.words = 10)
exc_mean <- unlist(lapply(exc, mean))



# plot of scaled scores, including mean of semantic coherence & exclusivity
semcoh <- scale(coh_mean)
exclus <- scale(exc_mean)
semexc <- rowMeans(cbind(semcoh, exclus)) # mean of scaled sem & exc

plot(semcoh, type = "l", col = "blue", xaxt = "n", lwd = 2, ylim = c(-3, 3), 
     main = "Topic Quality", ylab = "Scaled Score", xlab = "Number of Topics")
lines(exclus, type = "l", col = "red", lwd = 2)
lines(semexc, type = "l", col = "black", lwd = 3)
abline(h = max(semexc), v = which.max(semexc), col = "gray", lwd = 2, lty = "dashed")
axis(1, at = 1:length(K_range), labels = K_range)
legend("bottomright", c("Semantic Coherence", "Exclusivity (LDAvis lambda = 0)", 
                        "Mean Coherence & Exclusivity"), 
       col = c("blue", "red", "black"), lty = "solid", lwd = 2)




# Top Terms ---------------------------------------------------------------

topterms_list <- list()
for (i in 1:length(candidates_2020)){
  temp <- tosca::topWords(getTopics(getLDA(candidates_2020[[i]])), 10)
  topterms_list[[i]] <- apply(temp, 2, paste, collapse = ", ")
}




# Topic Prevalence --------------------------------------------------------

theta <- getEstimators(getLDA(candidates_2020[[t]]))$theta
prevalence <- rowMeans(theta)

df <- data.frame(topterms_list[[t]], prevalence)




# Best model of 2020 ------------------------------------------------------

# Model K = 250

proto_2020 <- candidates_2020[[which(K_range == 250)]]
lda_2020 <- getLDA(proto_2020)
beta_2020 <- getEstimators(lda_2020)$phi
theta <- getEstimators(lda_2020)$theta

topterms <- tosca::topWords(getTopics(lda_2020), 10)
topterms <- apply(topterms, 2, paste, collapse = ", ")

prevalence <- rowMeans(theta)

topic <- data.frame(ID = 1:(lda_2020$param$K), topterms, prevalence)

saveRDS(proto_2020, file = "./Candidates 2020/proto_2020.RDS")
saveRDS(lda_2020, file = "./Candidates 2020/lda_2020.RDS")
saveRDS(beta_2020, file = "./Candidates 2020/beta_2020.RDS")

saveRDS(topic, file = "./Candidates 2020/topic_table_2020.RDS")
write.csv(topic, file = "./Candidates 2020/topic_table_2020.csv",row.names = FALSE)



