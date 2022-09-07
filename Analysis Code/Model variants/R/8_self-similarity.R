###
# topic self-similarity
###



## libraries ----

library(lsa)
library(ldaPrototype)



## load objects ----

meta <- readRDS("./input data/meta.RDS")

# all years used for RollingLDA (after init model)
unique_years <- sort(unique(meta$PY))


sim_mean_names_by_year <- readRDS("./Model variants/RDS/sim_mean_names_by_year.RDS")
top_candidates <- names(head(sort(sim_mean_names_by_year, decreasing = TRUE), 5))



for (j in 1:length(top_candidates)){
  
  # complete candidate model
  lda <- readRDS(paste0("./Model variants/", top_candidates[j],"/lda.RDS"))
  
  # list of year's beta_by_year probabilities (buzz_search.R) by topic
  beta_list <- readRDS(paste0("./Model variants/", top_candidates[j],"/beta_list_by_year.RDS"))
  
  K <- getK(lda)
  
  # topic table (for labels)
  topic <- readRDS(paste0("./Model variants/", top_candidates[j],"/topic.RDS"))
  
  # topic evolution (top 10 terms over unique_years)
  topic_evo <- readRDS(paste0("./Model variants/", top_candidates[j],"/topic_evo_by_year.RDS"))
  
  
  
  ## compute self-similarities across years ----
  
  sim_first <- list()
  sim_last <- list()
  sim_prev <- list()
  
  for (t in 1:K){
    
    # with first year
    res <- list()
    for (i in 1:length(unique_years)){
      res[[i]] <- as.numeric(lsa::cosine(beta_list[[1]][t,], beta_list[[i]][t,]))
    }
    sim_first[[t]] <- unlist(res)
    
    # with last year
    res <- list()
    for (i in 1:length(unique_years)){
      res[[i]] <- as.numeric(lsa::cosine(beta_list[[length(unique_years)]][t,], beta_list[[i]][t,]))
    }
    sim_last[[t]] <- unlist(res)
    
    # with previous year
    res <- list()
    res[[1]] <- NA
    for (i in 2:length(unique_years)){
      res[[i]] <- as.numeric(lsa::cosine(beta_list[[i]][t,], beta_list[[i-1]][t,]))
    }
    sim_prev[[t]] <- unlist(res)
    
  }
  
  saveRDS(sim_first, file = paste0("./Model variants/", top_candidates[j], "/sim_first.RDS"))
  saveRDS(sim_last, file = paste0("./Model variants/", top_candidates[j], "/sim_last.RDS"))
  saveRDS(sim_prev, file = paste0("./Model variants/", top_candidates[j], "/sim_prev.RDS"))
  
  }
  





## plot -----

sim_plot <- function(t, legend = TRUE) {
  plot(sim_prev[[t]], type = "l", ylab = "Cosine Similarity", xlab = "Year",
       xaxt = "n", main = topic$Label[t], ylim = c(0,1))
  axis(1, at = 1:length(unique_years), labels = unique_years)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey90")
  grid(NULL, NULL, lty = 7, col = "white") 
  lines(sim_prev[[t]], type = "l", lwd = 2, col = "black")
  lines(sim_first[[t]], type = "l", lwd = 2, col = "blue")
  lines(sim_last[[t]], type = "l", lwd = 2, col = "orange")
  if (legend == TRUE){
    legend(1, 0.2, legend = c("previous year ", "first year", "last year"), lty = 1, lwd = 2,
           col = c("black", "blue", "orange"), cex = 0.8)
  }
}

par(mfrow=c(5,5))

for (i in 1:K){
  sim_plot(i, legend = FALSE)
}

dev.off()



sim_plot(1) # check  example
topic_evo[[1]] # see topic evolution
