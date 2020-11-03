### This folder contains the R script to demonstrate how to use modis data to detect the forest fire impact signals. 
#### Set up a conda environment for the R
- https://docs.anaconda.com/anaconda/user-guide/tasks/using-r-language/
- using anancoda to install R and its packages in the atmosphere 
- install R packages is always difficulty, especially for sf R package in the ubutuo environment  
  this page contains useful information of how to install R spatial packages: https://www.r-bloggers.com/2020/03/installing-spatial-r-packages-on-ubuntu/

#### Use jupyter lab to download modis LAI maps and clip the LAI based on the watershed boundary map (American River watershed), and create the time series of LAI (4 days) per subbasin level. 
- jupyter lab --no-browser --allow-root --ip=0.0.0.0 --NotebookApp.token='' --NotebookApp.password='' --notebook-dir='/scratch/Foss2020Team6project/'
- 128.196.142.11:8888
- use "downloadmodis.ipynb" 
