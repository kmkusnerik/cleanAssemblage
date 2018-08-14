#' cleanAssemblage: A cleaner function for taxa abundance data
#'
#' A function to clean imported taxa abundance datasets to a workable
#' datafile. The function combines samples of live assemblages from the
#' same site but which are different samples. It also removes all samples
#' below a minimum sample size.
#'
#' @param df A dataframe of imported taxa abundances. Should be arranged such
#' that rows are samples, columns at the start should describe sample
#' sample metadata while later columns after that should contain abundance of
#' various taxa
#' @param minN The minimum sample size desired
#' @param n.meta.col An integer of how many columns in your df are metadata
#' describing the samples (ie. not the counts of taxa abundance)
#'
#' @details The function combines taxa counts for live assemblage samples
#' coming from the same collection site (and thus from a single community).
#'  Finally, it discards all samples remaining below a minimum sie threshold.
#'
#' @return A dataset cleaned for live-dead analysis
#'
#' @note At the moment the function only works if the sample metadata is stored in
#' the first set of columns (1-n.meta.col) but can contain as many columns as needed.
#' The function also reduces the final cleaned dataset down to a number of 'important
#' metadata' currently listed in the first line of the function. Currently set to
#' return Location, Latitude, Longitude, Type, and System. All others will not be in
#' cleaned dataset. Ability to alter  which metadata is desired in the
#' cleaned dataset to be implemented later ...when I have the time.
#'
#' @author Kristopher Kusnerik\ email{kmkusnerik@@ufl.edu} based on code from
#' Michal Kowalewski (University of Florida)
#'

cleanAssemblage <- function (df, minN, n.meta.col){
  meta <- which(colnames(df) == 'Location' | colnames(df)=='Latitude'| colnames(df)=='Longitude'| colnames(df)=='Type'| colnames(df)=='System')
  site <- which(colnames(df) == 'Location')

  meta.rep <- rep(1:n.meta.col)
  junk.meta <- which(meta.rep %in% meta==FALSE)

  df.meta <- df[ ,meta]
  live.s <- which(df$Type == 'Live')
  live.d <- df[live.s, ]

  live.pool <-
    na.omit(sapply(live.d[, -c(1:n.meta.col)], function (x)
      tapply(x, live.d[, site], sum)))

  live.fct <- na.omit(sapply(live.d[, meta],
                             function (x)
                               tapply(x, live.d[, site], function (x)
                                 as.character(x[1]))))

  live.fin <- data.frame(live.fct, live.pool)
  dff <- rbind(df[-live.s, -junk.meta], live.fin)
  dff$Latitude <- as.numeric(dff$Latitude)
  dff$Longitude <- as.numeric(dff$Longitude)
  dff$Location <- droplevels(dff$Location)
  dff$Type <- droplevels(dff$Type)
  dff$System <- droplevels(dff$System)

  dff <- dff[which(rowSums(dff[, -(1:5)]) >= minN), ]
  return (dff)
}
