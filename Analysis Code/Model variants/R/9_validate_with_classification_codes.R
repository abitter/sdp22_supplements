# validate results using classification codes (SH)

library(quanteda)
library(ldaPrototype)

meta <- readRDS("./input data/meta.RDS")



# classifications ---------------------------------------------------------

dfm_sh <- dfm(tokens(meta$SH, remove_punct = TRUE))

# top 10 classification codes (to compare with models most prevalent topics)
# https://www.apa.org/pubs/databases/training/class-codes
top_sh <- topfeatures(dfm_sh, 10)
top_50_sh <- topfeatures(dfm_sh, 50)

# get labels
# Code labels can be retrieved from http://www.apa.org/pubs/databases/training/class-codes.aspx
apa_codes <- readr::read_delim("./Model evaluation/apa_codes.csv", 
                               ";", escape_double = FALSE, trim_ws = TRUE)

names(top_sh) <- stringi::stri_replace_all_regex(names(top_sh),
                                                      pattern = apa_codes$SH2,
                                                      replacement = apa_codes$SHfull,
                                                      vectorize = FALSE)

# 3310 Psychotherapy & Psychotherapeutic Counseling
# 2340 Cognitive Processes
# 3120 Personality Traits & Processes
# 3290 Physical & Somatoform & Psychogenic Disorders
# 2840 Psychosocial & Personality Development
# 2520 Neuropsychology & Neurology
# 3215 Neuroses & Anxiety Disorders
# 3530 Curriculum & Programs & Teaching Methods
# 3315 Psychoanalytic Therapy
# 3620 Personnel Management & Selection & Training




# topics ------------------------------------------------------------------

# classifications
classifications <- meta$SH

# candidate models (top 5 similarity with 2020)
candidates <- c("variant_200_2010", "variant_200_2005", "variant_200_1995", "variant_300_2010", "variant_200_2015")

# remove docs that were dropped in estimation processes
docs_all_cand <- list()
for (i in 1:length(candidates)){
  docs_all_cand[[i]] <- readRDS(paste0("./Model variants/", candidates[i],"/docs_all.RDS"))
}

classifications_cand <- list()
for (i in 1:length(candidates)){
  keep <- meta$DFK %in% names(docs_all_cand[[i]])
  classifications_cand[[i]] <- classifications[keep]
}

rm(classifications)


# get models' overall most prevalent topics
top_cand <- list()
for (i in 1:length(candidates)){
  topic <- readRDS(paste0("./Model variants/", candidates[i],"/topic.RDS"))
  top_cand[[i]] <- head(topic[order(-topic$n_docs), 1], 25)
}


# get theta
theta_cand <- list()
for (i in 1:length(candidates)){
  lda <- readRDS(paste0("./Model variants/", candidates[i],"/lda.RDS"))
  theta_cand[[i]] <- t(getEstimators(lda)$theta)
}


# validate with classification codes
source("./Model evaluation/validate_with_classifications.R")
validate_cand <- list()
for (i in 1:length(candidates)){
  validate_cand[[i]] <- validate_TM(theta_cand[[i]], classifications_cand[[i]], keep = top_cand[[i]],
                              threshold = 0.01, scale_limit = 0.2, codes = apa_codes)
}


# plots
for (i in 1:length(candidates)){
  png(paste0("./Model variants/validate_", candidates[i], ".png"), 
      width = 1200, height = 1200, res = 100)
  validate_cand[[i]][[3]]
  dev.off()
}



# compare classification frequency rankings

# get all relative frequencies of topics by classification code (sh)
topic_sh_cand <- list()
for (i in 1:length(candidates)){
  topic_sh_cand[[i]] <- validate_TM(theta_cand[[i]], classifications_cand[[i]], 
                                    class_threshold = 0, codes = apa_codes)[[1]]
}


# get all classification frequencies original metadata
sh_freqs <- as.data.frame(colSums(dfm_sh))
names(sh_freqs) <- "frequency"
sh_freqs$classification <- rownames(sh_freqs)
sh_freqs$classification <- stringi::stri_replace_all_regex(sh_freqs$classification,
                                                 pattern = apa_codes$SH2,
                                                 replacement = apa_codes$SHfull,
                                                 vectorize = FALSE)


# mean shares of classifications across all topics
topic_sh_shares_cand <- list()
for (i in 1:length(candidates)){
  temp <- as.data.frame(rowMeans(topic_sh_cand[[i]]))
  names(temp) <- paste0("share_", candidates[i])
  temp$classification <- rownames(temp)
  topic_sh_shares_cand[[i]] <- temp
}


# merge
merge <- dplyr::left_join(sh_freqs, topic_sh_shares_cand[[1]], by = "classification")

for (i in 2:length(candidates)){
  merge <- dplyr::left_join(merge, topic_sh_shares_cand[[i]], by = "classification")
}

# correlations of shares with actual metadata frequencies
cors <- list()
for (i in 1:length(candidates)){
  cors[[i]] <- cor(merge[,(i+2)], merge$frequency, use = "complete.obs")
}
cors <- unlist(cors)
names(cors) <- candidates


# save
# saveRDS(validate_cand, file = "./Model variants/validate_cand.RDS")
saveRDS(merge, file = "./Model variants/merge.RDS")
saveRDS(cors, file = "./Model variants/cors.RDS")


# add to candidate_metrics_table
candidate_metrics <- read.csv("././Model variants/candidate_metrics.csv")
candidate_metrics$correlations <- cors

candidate_metrics_scaled <- candidate_metrics[,-c(1,6)]
candidate_metrics_scaled[,2:5] <- apply(candidate_metrics_scaled[,2:5], 2, scale)
candidate_metrics_scaled$mean_scaled <- rowMeans(candidate_metrics_scaled[,2:5])

candidate_metrics <- candidate_metrics[,-c(1,6)]
candidate_metrics$mean_scaled <- candidate_metrics_scaled$mean_scaled

write.csv(candidate_metrics, file = "./Model variants/candidate_metrics.csv")

