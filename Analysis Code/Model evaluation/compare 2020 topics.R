# compare topics of 2020 only model with variant's 2020 topics

# get data ----

topic_2020 <- readRDS("./Candidates 2020/topic_table_2020.RDS")

# variant topic table
topic <- readRDS("./Model variants/variant_200_2010/topic.RDS")

# topics with memory effect
topic_evo <- readRDS("./Model variants/variant_200_2010/topic_evo.RDS")

# year-specific topics
topic_evo_by_year <- readRDS("./Model variants/variant_200_2010/topic_evo_by_year.RDS")



# get 2020 topics from topic_evo_by_year ----

# get index
ind_2020 <- which(colnames(topic_evo_by_year[[1]]) == "2020")

# extract 2020 topics
extract_2020 <- lapply(topic_evo_by_year, function(x){
  res <- x[,ind_2020]
  res <- paste(res, collapse = ", ")
  })

topic$Terms_2020 <- unlist(extract_2020)


# rename ----

topics_2020_reference <- topic_2020
topics_2020_variant <- topic[,c(7,4)]

names(topics_2020_reference)[2] <- "global_terms"
names(topics_2020_variant)[1] <- "evolution_terms"
topics_2020_variant$ID <- rownames(topics_2020_variant)
topics_2020_variant <- topics_2020_variant[,c(3,1,2)]



# inspect topics ----

View(topics_2020_reference)
View(topics_2020_variant)


# inspect similarities ----

variant_names <- readRDS("./Model variants/RDS/variant_names.RDS")
ind <- which(variant_names %in% "variant_200_2010")

sim_vectors_by_year <- readRDS("./Model variants/RDS/sim_vectors_by_year.RDS")

# for all 250 topics of 2020 model, show which topic of variant has the highest cosine similarity
sims <- unlist(sim_vectors_by_year[ind])

sims_df <- data.frame("ID_reference" = 1:250,
                      "cosine_sim" = sims,
                      "ID_variant" = names(sims),
                      "prevalence_reference" = topic_2020$prevalence)

# similarity less than 50 %
table(sims < .5)

ind <- which(sims < .5)
mean(sims_df$prevalence[ind])
mean(sims_df$prevalence)

# 45 not matched
temp <- sims_df[ind,]

# 11 not matched with above average prevalence
temp <- temp[temp$prevalence_reference > .004,]
ind_ref <- temp$ID_reference

# inspect these 11 topics
topics_2020_reference$global_terms[ind_ref]


# correlation similarity and prevalence
cor.test(sims_df$cosine_sim, sims_df$prevalence_reference)


topics_2020_reference$global_terms[236]
topics_2020_variant$evolution_terms[65]





# save ----

saveRDS(sims_df, file = "./Model variants/sims_df.RDS")
saveRDS(topics_2020_reference, file = "./Model variants/topics_2020_reference.RDS")
saveRDS(topics_2020_variant, file = "./Model variants/topics_2020_variant.RDS")

write.csv(sims_df, file = "./Model variants/sims_df.csv", row.names = FALSE)
write.csv(topics_2020_reference, file = "./Model variants/topics_2020_reference.csv", row.names = FALSE)
write.csv(topics_2020_variant, file = "./Model variants/topics_2020_variant.csv", row.names = FALSE)
