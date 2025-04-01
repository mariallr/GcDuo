# GcDuo: software for GCxGC-MS analysis
[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/) 

## Description

`GcDUO` package contains a collection of functions for computing processing of GCxGC-MS data covering: peak picking, deconvolution, alignemnt, and annotation. 

## Installation

Required package 'EBImage`

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

BiocManager::install("EBImage")


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
