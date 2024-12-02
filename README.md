# Identify karst groundwater flow paths using cave water drips
### Jarette Greene, Dr Khasif Mahmud
### Abstract: 
Cave drip water response to surface meteorological conditions is complex due to the heterogeneous nature of water movement in the karst unsaturated zone. Therefore, we aim to understand infiltration water hydrology in the Cretaceous karst formation of south-central Texas. This project utilizes the spatial survey of 20 automated cave drip loggers installed in Natural Bridge Caverns. Drip loggers were set up in approximate transects throughout the cavern from higher to lower elevations and have been monitored since May 2023. We analyzed drip time series for one hydrological year from all 20 sites taking advantage of an ensemble of drip loggers to extract common properties by clustering. A combination of multidimensional scaling and clustering by k means was utilized to classify similar drip types based on time series analysis. We used the coefficient of variation (COV) across different sampling frequencies (15 minutes to monthly) to determine the optimum sampling frequency that minimizes sampling artifacts while maximizing the capture of natural variability. The clustering reveals four unique drip regimes for better characterization of the limestone environment where flow occurs via fractures and storage. We found COV decreases with lower sampling frequency for all four cluster groups. We concluded that both daily and weekly sampling frequencies give the minimum COV, which does not change significantly with a finer sampling frequency and therefore we used one day as the optimum sampling frequency. By analyzing COV and cluster groupings, we discovered that age of limestone formation affects infiltration water flow pathways in heterogeneous karst formation.

### Files

|   #   | File             | Description                                        |
| :---: | ---------------- | -------------------------------------------------- |
|   1   | [Drip_analysis_NBC_HWC](https://github.com/Jarette/MDS_Clustering_NBC_HWC/blob/main/Drip_analysis_NBC_HWC.m) | Matlab script to perform cluster analysis.      |
|   2   | [Data](https://github.com/Jarette/MDS_Clustering_NBC_HWC/tree/main/Data)  | Folder containing neccessary input files.       |

### Description:

In this project, a MATLAB Script was created by Dr Khasif Mahmud to perform a similar analysis on Golgotha Cave Australia in 2018. This 2018 script was then adapted and updated by Researcher Jarette Greene. This updated Script uses drip time series data collected from Natural Bridge Cavern, Texas, and Harrie Wood Cavern, Australia, and uses multidimensional and the k-means algorithm to cluster the various drip time series. This clustering then allowed us to obtain a better understanding of water infiltration into these karst features and also gain a better understanding of the limestone environment ( fracturing, pooling).

### Running Script

  - Change path to location where input data is being stored on users PC.
    
  - Begin running code using standard matlab commands or run button.
      - User will then be propmted to pick which data they would like to be analyzed  (nbc - Natural Bridge Cavern, hwc - Harrie Wood Cave).
      - The user will then be prompted to pick the number of clusters (2,3,4).
      - Then the sampling frequency they would like to use ( daily, weekly, monthly).
   
  - After entering responses for all prompts the script then begins to produce results in the form of a variety of graphs and xlsx files containing statistical data.

  ### NOTE:
  - For the prompt in this script, the only error checking/text normalization performed is converting all inputs into lowercase. If the user were to misspell or enter responses outside those considered by the script, the script would error or produce incorrect results.

### Sample Results
  - Example of output graphs:

<img src="https://github.com/Jarette/MDS_Clustering_NBC_HWC/blob/main/images/Cluster%201.png" width="500">                       <img src="https://github.com/Jarette/MDS_Clustering_NBC_HWC/blob/main/images/Cluster%202.png" width="500"> 
                                                                  <img src="https://github.com/Jarette/MDS_Clustering_NBC_HWC/blob/main/images/COV%20vs%20SAMP%20FREQ%20(Correclation%20Coeff).png" width="700">


## Researchers: 

###  Dr Khasif Mahmud
<img src="https://static.wixstatic.com/media/f14632_9e17a82d5cf84c42ac7825c9a9e123dd~mv2.jpg/v1/fill/w_348,h_488,al_c,q_80,usm_0.66_1.00_0.01,enc_avif,quality_auto/f14632_9e17a82d5cf84c42ac7825c9a9e123dd~mv2.jpg" width="500">

### Position:
Mentor

### Bio: 
[Link to Bio](https://kashifmahmud.wixsite.com/kashif-mahmud/about-me)


|---------------------------------------------------------------------------------------------------------------------------------|

### Jarette Greene
<img src="https://user-images.githubusercontent.com/111944626/186936341-5004ae15-9c21-4b9d-9a7d-d8aa302df928.jpeg" width="500">

### Position:
Research Assistant

### Email: 
jarettegreene09@gmail.com

### Bio:
[Link to Bio](https://kashifmahmud.wixsite.com/kashif-mahmud/jarette-greene)

