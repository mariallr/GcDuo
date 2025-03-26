# GcDuo: software for GCxGC-MS analysis
[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/) 

## Description

`GcDUO` package contains a collection of functions for computing processing of GCxGC-MS data covering: peak picking, deconvolution, alignemnt, and annotation. 

## Installation

#### Beta/Github release:

Installation using R package devtools:

```r
install.packages("devtools")
devtools::install_github("mariallr/GcDUO")
```

## Usage

## Read data: 

For reading the data the user needs to indicate: 
- Modulation time: it corresponds to the second column time
- m/z fragments range acquired: the range of m/z fragments acquired in by the mass spec. 

```r
GcDuoObject <- GcDuo::readFolderCDF(folderPath = path, 
                                    modulationTime = 2.5, #the second dimension time (seconds)
                                    mzRange = c(35, 440)) #range of m/z fragments acquired
```

We obtain the S3 object with the following information:
- `$data`: the 4D array with  the files intensities. The order is the following: samples x m/z fragments x retention 1 x m retentio.n- 2.
- `$time1d`: the time values cor respondig to the 1st chromatographic separation.
- `$time2d`: the time values correspondig to the 2nd chromatographic separation.
- `$mz`: the m/z fragments in the array
- `$files`: the CDF files contained in the object

Now we can do a visual check of the data.

## Plot raw data

We can visualize the data as a contour plot. This visualization is only possible sample by sample, indicate in the argument `sampleNum` the file position to plot.


```{r}
visualizeChrom3D(GcDuoObject, sampleNum = 5)
```

To compare between all samples we can perform a chromatogram TIC visualization. Indicate in the argument `samples` which samples you want to compare. 

```{r}
visualizeChrom2D(GcDuoObject, samples = 1:5)
```

# Peak picking 

This function detect the peaks contained in the array. The user needs to modify the signalNoiseRatios at convenience. 
- signalNoiseRatio1: relative to first dimension
- signalNoiseRation2: relative to second dimension
- win_size: the number of modulation in each window, if you observe big peaks increase it. However the minimum value recomended is 3 modulations. 

```{r}
DuoResults <- GcDuo::procesData(GcDuoObject, 
                             signalNoiseRatio1 = 50, #the threshold for intensity in RT1
                             signalNoiseRatio2 = 50,  #the threshold for intensity in RT2
                             win_size = 3)
```

The object contains: 

- `$ind_peak`: table of the features detected sample by sample
- `$spec_ind`: spectra of the features detected sample by sample
- `$cor_spectra`: consensus spectra of common features between samples

# Annotation

Annotation depends on the library used, for that the first step is to import your library to the R environment. 

Functions only execute:
```{r}
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

    # Nist
    rt_ <- unlist(regmatches(tolower(rts), gregexpr(
      'ri.*[0-9]*|retention.index.*[0-9]*|retentionindex.*[0-9]*', tolower(rts))))
    rtt <- gsub('ri.|retention.index. |retentionindex. ','',rt_, ignore.case=TRUE)

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
                  db.id = id, lib.id = library_id, ri = rtt, ri2 = rts, np = np,spectra=spectra))
    }

  }

  li <- lapply(li, getmsp)
  return(li)
}

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
```

Import your msp:

```{r}
libs <- getMSP(msp_file)

spectramatrix_libs <- spectramatrix(libs)

libs_dic <- NULL
for(i in 1:length(libs)) {
  libs_dic <- rbind(libs_dic, data.frame("comp" = libs[[i]]$name,
                                         "exactmass" = ifelse(length(libs[[i]]$exactmass) == 0,
                                                              NA, libs[[i]]$exactmass),
                                         "formula" = ifelse(length(libs[[i]]$formula) == 0, 
                                                            NA, libs[[i]]$formula),
                                         "inchikey" = ifelse(length(libs[[i]]$inchikey) == 0, 
                                                             NA, 
                                                             libs[[i]]$inchikey), 
                                          "ri" = ifelse(length(libs[[i]]$ri) == 0, 
                                                             NA, 
                                                             libs[[i]]$ri)))
}

lib <- list(data = libs_dic, spectra_matrix = spectramatrix_libs)

rm(libs_dic)
rm(spectramatrix_libs)
rm(libs)
```

Now we annotate the features:

- match_cutoff: minimum value to consider good match
- RI: if kovats retention index (RI) is used indicate here the RT-RI alkanes
- RIrange: the accepted variability for RI

```{r}
DuoResults <- GcDuo::duoID(DuoResults, 
                            lib, 
                            match_cutoff = 0.7, 
                            RI = NULL, 
                            RIrange = 30)
```

- `$id_peak` : annotations for each feature. 

# Data pulishment + relative quantification

To obtain the areas and intensities for each feature and each sample: 

- win_width: an average width for the peaks. 

```{r}
finalGcDuo <- GcDuo::ConsProcessing(GcDuoObject, 
                                    DuoResults, 
                                    lib, 
                                    win_width = 20)


```

We obtain the finalGcDuo object with: 

- `$peaks`: summary of the peaks detected
- `$area`: area for each peak and sample
- `$intensity`: intensity for each peak and sample
- `$spectra`: spectra for each peak and sample

## Explore the results

```{r}
# For the identified peaks
View(finalGcDuo$peaks)

#For the area of the peaks
View(finalGcDuo$area)
```

Plot some peak that you are interested: 

- peakid: the id of the peak obtained in `finalGcDuo$peaks`
- type: "tic" for total ion chromatogram or "eic" for extracted ion chromatogram. In the second option you need to indicate the m/z to visualize. 

```{r}
PlotPeak(GcDuoObject, finalGcDuo, peakid = "24-1", type = "eic", win_width = 30)

PlotPeak(GcDuoObject, finalGCDuo, peakid = "24-1", type = "tic", win_width = 10)
```

Compare spectra with the library: 

```{r}
SpectraView(GcDuoObject, finalGCDuo, peakid = "22-1", lib_comp = lib,
                             mzRange = c(30,600))
```

Thanks for using GcDuo :)

Please fill an issue if you have any question or problem
