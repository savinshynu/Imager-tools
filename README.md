# Imager-tools for Orville Wideband Imager at LWA-SV
The LWA-SV is the second LWA radio telescope station located at the Sevilleta Wildlife Refuge. The broadband imager at LWA-SV can produce all-sky images of the sky every 5 seconds in realtime (20 MHz bandwidth) and this constitute roughly ~ 1 TB of all-sky images stored in a LWA specific OIMS format (.oims files). This format is similar to the HDF5 architecture where the datasets are multidimensional array along with the metadata of the observation. The multidimensional array consists of [nTimes, nFreqs, nStokes, 128 Xpixels, 128 Ypixels] where nTimes: number of time samples, NFreqs: number of frequency channel, nStokes: number of Stokes parameters, Xpixel: number of pixels in the X axis and Ypixel: number of pixels in the Y axis.  Usually each one hour of all-sky images are stacked into a single .oims file and it is 5 dimensional array with a shape of [720, 196, 4, 128, 128] in the float 64 format.

The all-sky images from LWA-SV can be used to monitor energetic explosions and dyanamic events (transients) happening in our universe including our own Earth's atmosphere.

This repository contain python scripts to

- Read image dataset (.oims files). 
- Conduct basic analysis of LWA-SV all-sky images. 
- Transient search pipeline ingesting each .oims files conducting data preprocessing, continuous image subtraction, threshold based detection and listing the candidates in a text file with their locations and time of detections.
- Script to filter candidates in multiple stages:
  	* Filter nearby VLA Low Frequency Sky (VLSS) sources
	* Filter low SNR scintillation events
	* Filter airplanes
	* Filter narrowband radio frequency interference (RFI) events
- Tools to understand candidates arising from scintillation of radio sources and calibrate the data.
- Tools to plot the light curves, broadband spectra, images and movies of interesting candidates.

The ***Jupyter notebook [link](Orville_image_processing_pipeline.ipynb) demonstrates how the processing of images can result in the discovery of interesting transient sources in the sky.*** The notebooks provides an introduction on how use LWA images to find interesting transient sources or the energetic explosions occuring in our atmosphere and universe.

More description of the data products, transient search pipeline and results after the analysis of 2400 hours (48 TBs) of data are given in this published paper [link](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021JA029296). Please cite this work if your are using this code for your research work.