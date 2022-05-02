library(MALDIquant)
library(MALDIquantForeign)
tol=0.5
snr=4

spectra <- transformIntensity(CalDigest,method="sqrt")
spectra <- smoothIntensity(spectra, method="SavitzkyGolay",halfWindowSize=10)
spectra <- removeBaseline(spectra, method="SNIP",iterations=100)
spectra <- calibrateIntensity(spectra, method="TIC")
spectra <- alignSpectra(spectra,halfWindowSize=20,SNR=snr,tolerance=tol,warpingMethod="lowess")
peaks <- detectPeaks(spectra, method="MAD",halfWindowSize=20, SNR=snr)
peaks <- binPeaks(peaks, tolerance=0.1)
peaks <- filterPeaks(peaks, minFrequency=0.25)

featureMatrix <- intensityMatrix(peaks, spectra)


head(featureMatrix[, 1:3])


write.csv(featureMatrix, "TOl1_.csv")


rm(featureMatrix, peaks, spectra)


