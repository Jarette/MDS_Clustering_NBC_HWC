# Identify karst groundwater flow paths using cave water drips
### Jarette Greene, Dr Khasif Mahmud
### Abstract: 
Cave drip water response to surface meteorological conditions is complex due to the heterogeneous nature of water movement in the karst unsaturated zone. Therefore, we aim to understand infiltration water hydrology in the Cretaceous karst formation of south-central Texas. This project utilizes the spatial survey of 20 automated cave drip loggers installed in Natural Bridge Caverns. Drip loggers were set up in approximate transects throughout the cavern from higher to lower elevations and have been monitored since May 2023. We analyzed drip time series for one hydrological year from all 20 sites taking advantage of an ensemble of drip loggers to extract common properties by clustering. A combination of multidimensional scaling and clustering by k means was utilized to classify similar drip types based on time series analysis. We used the coefficient of variation (COV) across different sampling frequencies (15 minutes to monthly) to determine the optimum sampling frequency that minimizes sampling artifacts while maximizing the capture of natural variability. The clustering reveals four unique drip regimes for better characterization of the limestone environment where flow occurs via fractures and storage. We found COV decreases with lower sampling frequency for all four cluster groups. We concluded that both daily and weekly sampling frequencies give the minimum COV, which does not change significantly with a finer sampling frequency and therefore we used one day as the optimum sampling frequency. By analyzing COV and cluster groupings, we discovered that age of limestone formation affects infiltration water flow pathways in heterogeneous karst formation.

### Files

|   #   | File             | Description                                        |
| :---: | ---------------- | -------------------------------------------------- |
|   1   | Main.cpp         | Main driver of my project that launches game.      |
|   2   | HelperClass.cpp  | Helper class that holds movement functions         |
|   3   | TextureClass.cpp | Helper class that assists with textures and images |

### Description:

In this project a MATLAB Script that was created by Dr Khasif Mahmud that was created to perform similar analysis on Golgotha Cave Australia in 2018. This 2018 script was then adapted and updated by Researcher Jarette Greene. This updated Script uses drip time series data collected from Natural Bridge Cavern, Texas, and Harrie Wood Cavern, Australia, and using multidimensional and the k-means algortithm to cluster the various drip time series. This clustering then allowed us to obtain a better understanding of water infiltration into these kasrt features and also gain a better understanding of the limestone enviroment (eg fracturing, pooling).

### Running Script

  - Change path to location where input data is being stored on users PC.
      'addpath(%insert your own path here)'
    
  - Begin running code using standard matlab commands or run button
      - user will then be propmted to pick which data they would like to be analyzed  (nbc - Natural Bridge Cavern, hwc - Harrie Wood Cave)
      - the user will then be prompted to pick the number of clusters (2 - 4) then the sampling frequency they would like to use ( daily, weekly, monthly)
   
  - After entering responses for all prompts the script then begins to produce results in the form of a variety of graphs and xlsx files containing statistical data

  ### NOTE:
  - For the prompt this script the only error checking/ text normalization performed is to convert all inputs into lower case, if user entered missed spelled or responses not considered the script will error or produce incorrect results

### Sample Results
  - Example of output graphs:

<img src="https://github.com/Jarette/MDS_Clustering_NBC_HWC/blob/main/images/Cluster%201.png" width="600">                       <img src="https://github.com/Jarette/MDS_Clustering_NBC_HWC/blob/main/images/Cluster%202.png" width="600"> 

  
  

