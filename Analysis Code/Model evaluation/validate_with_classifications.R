
validate_TM <- function(theta, classifications, codes, threshold = 0.1, class_threshold = 0.1,
                        keep = 1:k, scale_limit = 1, sort_class = TRUE) {
  
  # theta is theta matrix, with topics in cols
  # classifications is a character vector of length = nrow(theta)
  
  library(quanteda)
  
  k <- ncol(theta)
  theta_sh <- as.data.frame(theta)
  theta_sh$SH <- classifications
  theta_sh <- theta_sh[theta_sh$SH != "",] # SH for Subject Heading
  
  
  topic_i_list <- list()
  
  for (i in 1:k) {
    topic_i_list[[i]] <- theta_sh[theta_sh[,i] > 0.5,] 
    if(nrow(topic_i_list[[i]]) == 0){                  
      topic_i_list[[i]][1,] <- c(rep(0, k), "none")
    }                  
  }                                                    
  
  
  ### frequencies of main classifications by topic ----
  
  topic_i_split <- list() # split multiple classifications
  topic_i_table <- list() # frequency table
  
  for (i in 1:k) {
    
    topic_i_split[[i]] <- unlist(strsplit(topic_i_list[[i]][,length(topic_i_list[[i]])], ","))
    
    topic_i_table[[i]] <- as.data.frame(table(topic_i_split[[i]]))
    names(topic_i_table[[i]]) <- c("SH", i)
  }
  
  
  ### combine ----
  
  # how many different classifications?
  classifications_names <- colnames(dfm(tokens(corpus(classifications))))
  
  topic_sh <- as.data.frame(classifications_names)
  names(topic_sh)[1] <- "SH"
  
  
  for (i in 1:k){
    topic_sh <- dplyr::left_join(topic_sh, topic_i_table[[i]], by = "SH")
  }
  rownames(topic_sh) <- topic_sh$SH
  topic_sh <- topic_sh[,-1] # drop SH column
  
  
  # NA to Zeros
  for (i in 1:k) {
    topic_sh[,i] <- ifelse(is.na(topic_sh[,i]), 0, topic_sh[,i])
    topic_sh[,i] <- ifelse(topic_sh[,i] != 0,
                           topic_sh[,i]/sum(topic_sh[,i]),
                           0)
  }
  
  
  # keep only classifications with relative frequency > class_threshold
  cond <- apply(topic_sh, 1, max) > class_threshold
  topic_sh <- topic_sh[cond,]
  
  
  if (sort_class == TRUE){
    # sort by classifications code
    topic_sh <- topic_sh[order(rownames(topic_sh)),]
  } else {
    # sort by most prominent classifications
    topic_sh <- topic_sh[order(-rowMeans(topic_sh)),]
  }
  
  
  # replace APA codes with labels
  rownames(topic_sh) <- stringi::stri_replace_all_regex(rownames(topic_sh) ,
                                                     pattern = codes$SH2,
                                                     replacement = codes$SHfull,
                                                     vectorize = FALSE)

  
  ### count number of black cells ----
  
  cells <- list()
  for (i in 1:k){
    cells[[i]] <- length(which(topic_sh[,i] > threshold)) # count how many values are above threshold
  }
  cells <- unlist(cells)
  
  #cells <- ifelse(unlist(cells) > 1, 2, unlist(cells))
  cells <- table(cells)/k # relative frequencies
  
  
  
  ### plot ----
  
  plot <- lattice::levelplot(t(topic_sh[,keep]),
                             main = paste0("Relative frequencies of database classifications by topic (k = ", k, ")"),
                             cex.main = 0.5,
                             col.regions = gray (1000:0/1000), 
                             xlab = "Topic",
                             ylab = NULL,
                             at = seq(0, scale_limit, threshold),
                             pretty = TRUE)
  
  
  
  ### results ----
  
  results <- list()
  
  results[[1]] <- topic_sh
  results[[2]] <- cells
  results[[3]] <- plot
  
  return(results)
  
}
