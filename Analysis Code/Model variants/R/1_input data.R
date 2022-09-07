# Script 1: Create input data



# Libraries ---------------------------------------------------------------

library(cld2)
library(readr)
library(quanteda)




# Read and prepare data ---------------------------------------------------

psyndex <- readRDS("~/RollingLDA/psyndex.RDS")


## language checks ----

# check language of titles and abstracts. Use English only.

fill <- rep(NA, nrow(psyndex))
languages <- data.frame("TI" = fill, "TIU" = fill, "TIUE" = fill, "ABH" = fill, "ABN" = fill)

languages$TI <- cld2::detect_language(psyndex$TI)
languages$TIU <- cld2::detect_language(psyndex$TIU)
languages$TIUE <- cld2::detect_language(psyndex$TIUE)
languages$ABH <- cld2::detect_language(psyndex$ABH)
languages$ABN <- cld2::detect_language(psyndex$ABN)

languages[is.na(languages)] <- ""

saveRDS(languages, file = "./languages.RDS")

apply(languages, 2, table)


## text ----

# standardized keywords, called "Controlled Terms" (CTs) in psyndex

# transform
transform_keywords <- function(x){
  x <- tolower(x)
  x <- gsub(" ", "_", x)
  x <- gsub("\\(", "RRR1", x)
  x <- gsub("\\)", "RRR2", x)
  x <- gsub("\\/", "RRR3", x)
  x <- gsub("'_", "RRR4_", x) # word must not end with '
  x <- gsub("@", "_as_", x) # some issue in PSYNDEX
  x <- gsub(",_", " ", x)
  return(x)
}

ct <- transform_keywords(psyndex$CT)

# add all English text and CTs to a vector

text <- character(nrow(psyndex))

for (i in 1:length(text)){
  text[i] <- ifelse(languages$TI[i] == "en",
                    paste(text[i], psyndex$TI[i], collapse = " "),
                    text[i])
  text[i] <- ifelse(languages$TIU[i] == "en",
                    paste(text[i], psyndex$TIU[i], collapse = " "),
                    text[i])
  text[i] <- ifelse(languages$TIUE[i] == "en",
                    paste(text[i], psyndex$TIUE[i], collapse = " "),
                    text[i])
  text[i] <- ifelse(languages$ABH[i] == "en",
                    paste(text[i], psyndex$ABH[i], collapse = " "),
                    text[i])
  text[i] <- ifelse(languages$ABN[i] == "en",
                    paste(text[i], psyndex$ABN[i], collapse = " "),
                    text[i])
  text[i] <- paste(text[i], ct[i], collapse = " ")
}



# check for empty text
table(text == "")
remove <- which(text == "")

# remove empty cases from text and psyndex
text <- text[-remove]
psyndex <- psyndex[-remove,]



## meta ----

# create object with metadata
names(psyndex)
meta <- psyndex[,c(2,5,3,4,13,14)]

# create new variable EMP (empirical research) from CM
meta$EMP <- NA
cond_emp <- grepl("cm10", meta$CM) | grepl("cm131", meta$CM)

meta$EMP <- ifelse(cond_emp, 1, meta$EMP)
meta$EMP <- ifelse(!cond_emp, 0, meta$EMP)
meta$EMP <- ifelse(meta$CM == "cm", NA, meta$EMP) # for records without method information

table(meta$EMP)


## sort meta and text by year ----

temp <- cbind(meta, text)
temp <- temp[order(temp$PY),]

meta <- temp[,1:7]
text <- temp[,8]
rm(temp)


## save RDS ----

saveRDS(text, file = "./input data/text.RDS")
saveRDS(meta, file = "./input data/meta.RDS")
saveRDS(psyndex, file = "./input data/psyndex.RDS")