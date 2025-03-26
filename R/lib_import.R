
#########################################################################
#     GcDuo - R package for GCxGC processing and PARAFAC analysis
#     Copyright (C) 2023 Maria Llambrich & Lluc Semente
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################

##########################
#' getMSP
#'
#' Read library MSP file
#'
#' @param file library msp file
#' @return list of library entries
#' @export
#'

getMSP <- function(file){
  li <- NULL

  #Get library name
  l_id <- unlist(stringr::str_split(msp_file, "/"))
  library_id <- stringr::str_split(l_id[length(l_id)], "\\.")[[1]][1]

  msp <- readLines(msp_file)
  # remove empty lines
  msp <- msp[msp != '']
  #number of compounds
  ncomp <- grep('^NAME:', msp, ignore.case = TRUE)

  # Divide by entry
  splitFactorTmp <- rep(1:length(ncomp), diff(c(ncomp, length(msp) + 1)))

  li <- split(msp,f = splitFactorTmp)

  # Transform to table

  getmsp <- function(x){
    namet <- x[grep('^NAME:',x, ignore.case=TRUE)]
    name <- gsub('^NAME: ','',namet, ignore.case=TRUE)
    exactmasst <- x[grep('^EXACT.MASS:|^EXACTMASS:',x, ignore.case=TRUE)]
    exactmass <- gsub('EXACT.MASS: |^EXACTMASS: ','', exactmasst, ignore.case=TRUE)
    casnum <- x[grep('^CASNO:',x, ignore.case=TRUE)]
    cas <- gsub('CASNO: ','', casnum, ignore.case=TRUE)
    formt <- x[grep('^FORMULA: ',x, ignore.case=TRUE)]
    formula <- gsub('^FORMULA: ','',formt,ignore.case = TRUE)
    inchit <- x[grep('^INCHIKEY: |^INCHI: ',x, ignore.case=TRUE)]
    inchikey <- gsub('^INCHIKEY: |^INCHI: ','',inchit,ignore.case = TRUE)
    ids <- x[grep('ID: ',x,ignore.case = TRUE)]
    id <- gsub('ID: ','',ids,ignore.case=TRUE)
    comm <- x[grep('^COMMENT: ',x,ignore.case = TRUE)]
    rts <- x[grep('ri. |retention.index. |retentionindex. ',tolower(x),ignore.case = TRUE)][1]

    # Number of m/z fragments
    npt <- x[grep('^Num Peaks: ',x, ignore.case=TRUE)]
    np <- gsub('^Num Peaks: ','',npt,ignore.case = TRUE)

    if(as.numeric(np) > 0){
      # matrix of masses and intensities
      massIntIndx <- which(grepl('\\s[0-9]', x) & !grepl(': ', x))
      massesInts <- unlist(strsplit(x[massIntIndx], '; '))
      massesInts <- unlist(stringr::str_squish(massesInts))
      massesInts <- strsplit(massesInts, ' ')
      mz <- unlist(lapply(massesInts, '[[', 1))
      ins <- unlist(lapply(massesInts, '[[', 2))
      # if any NAs remove from indx
      spectra <- cbind.data.frame(mz=mz,ins=ins)
      return(list(name=name,exactmass=exactmass, cas = cas, formula=formula, inchikey=inchikey,
                  db.id = id, lib.id = library_id, ri = rts, np = np,spectra=spectra))
    }

  }

  li <- lapply(li, getmsp)
  return(li)
}

##########################
#' spectramatrix
#'
#' Converts library entries to matrix
#'
#' @param libs library list obtained with getMSP
#' @return matrix of spectra
#' @export
#'

spectramatrix <- function(libs){
  require(dplyr)
  seq_mz <- seq(30,600)
  ref_matrix <- matrix(nrow = length(libs), ncol = length(seq_mz))
  ref_names <- NULL
  for (j in 1:length(libs)) {
    print(j)
    spec_ref <- tibble(libs[[j]]$spectra) |>
      mutate(ins = gsub(";", "", ins))
    spec_ref_ok <- spec_ref %>% mutate(
      mz = round(as.numeric(as.character(mz)),0),
      ins = as.numeric(as.character(ins))) %>%
      group_by(mz) %>%
      summarise(mz = mean(mz),
                ins = sum(ins)) %>%
      filter(mz >= min(seq_mz) & mz <= max(seq_mz)) %>%
      arrange(mz)
    if (nrow(spec_ref_ok) < length(seq_mz) & nrow(spec_ref_ok) != 0) {
      spec_ref_ok <- spec_ref_ok %>%
        add_row(mz = seq_mz[!(seq_mz %in% .$mz)], ins = 0) %>%
        arrange(mz)
    }

    ref_names <- c(ref_names, libs[[j]]$name)
    if (nrow(spec_ref_ok) == length(seq_mz)) {

      spec_ref_ok$rel <- spec_ref_ok$ins/max(spec_ref_ok$ins)

      ref_matrix[j,] <- spec_ref_ok$rel
    }
  }
  row.names(ref_matrix) <- as.vector(ref_names)
  ref_matrix <- ref_matrix[complete.cases(ref_matrix),]
  colnames(ref_matrix) <- seq_mz

  return(ref_matrix)
}
