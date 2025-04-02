# Satellite altimetry reveals intensifying global river water level variability
*****
## Introduction
This repository contains all necessary code for producing datasets and reproducing results for the manuscript "Satellite altimetry reveals intensifying global river water level variability" (in review). If you want to replicate our results from start to finish, the whole process boils down to two key steps:
+ ##### Step Ⅰ: Generate global virtual stations (VSs) and retrieve water levels from altimetry data.
+ ##### Step Ⅱ: Implement postprocessing on global river water levels, including fluctuations, change rates, extreme stage years, and seasonality.
Step I is computationally intensive, and we strongly recommend using parallel computing. Even with 20 cores running simultaneously, it will take about two months to complete the calculations. We anticipate that most users/readers will be more interested in Step II, which can be easily carried out on personal computers. To make it more convenient to reproduce and visualize our results, we've provided the global river water level dataset that was generated after completing Step I. A detailed user guide explaining our dataset's attributes is available via [this link](https://doi.org/10.5281/zenodo.14671453). This documentation has been uploaded along with our dataset and the processed results obtained after Step II.

*Next, we'll explain how to use our code (in **Windows**) and describe its various functions in detail.*
***
## Preparation
We've organized our code into four folders: `VS generation`, `Altimetry calculation`, `Post processing`, and `Dataset validation`. The first two folders correspond to Step I, while the last two are for Step II.
***
The code for Step I was written in **MATLAB** and requires version **2021b or later**, along with the **[Statistics and Machine Learning Toolbox](https://ww2.mathworks.cn/products/statistics.html)**, **[Financial Toolbox](https://ww2.mathworks.cn/products/finance.html)**, and **[Parallel Computing Toolbox](https://ww2.mathworks.cn/products/parallel-computing.html)**.
Note that we use parallel processing, and you'll need to download the [Progress monitor (progress bar) that works with parfor](https://ww2.mathworks.cn/matlabcentral/fileexchange/32101-progress-monitor-progress-bar-that-works-with-parfor) app from MATLAB Central File Exchange to monitor the calculation progress.
***
The code for Step II was writen in **Python** and requires version **3.10.0 or later**. You can run the code line by line in Jupyter Notebook after installing the packages mentioned in the `requirements.txt`. We recommend using a **[conda](https://www.anaconda.com/download/)** command to install these packages to ensure compatibility with your current Python version. Here's how to proceed:
```
conda install --yes --file requirements.txt
```
or if your Python is not equipped with **conda**:
```
pip install -r requirements.txt
```
If you encounter any issues while installing *[geemap](https://geemap.org/)* or *[arcpy](https://pro.arcgis.com/zh-cn/pro-app/latest/arcpy/get-started/installing-arcpy.htm)*, you'll need to visit their official documentation by clicking the hyperlinks on their names. Note that you'll still need to <u>purchase **[ArcGIS Pro](https://pro.arcgis.com/en/pro-app/index-geonet-allcontent.html)** to obtain a valid license</u> for using *arcpy* in a Python environment. Our basic calculations do not require these two packages. If you want to avoid these complications, you can apply the above process to `requirements_base.txt`.
When other packages encounter installation issues, you can opt to install them manually, which should resolve the respective problems.

#### Recommended approach
We also provide a Docker image for users, who are unfamiliar with Python package, to directly replicate our Python environment, which is available on [Docker Hub](https://www.docker.com/products/docker-hub/). This image includes our code and prepared datasets for Step II. To use it, first install [Docker Desktop](https://www.docker.com/) on your PC, then pull our image using the following command in the terminal:
```
docker pull fangchq/global_rivers:v202501
```
Then build a container euipped with our Python environment:
```
docker run -it --name my_river_test -p 8888:8888 fangchq/global_rivers:v202501
```
You'll go into the terminal of this new container, and activate our virtual environment:
```
exec bash
conda activate river_cal
```
Then you can start Jupyter Lab with a pre-mapped port:
```
jupyter lab --ip 0.0.0.0 --port 8888 --no-browser --allow-root
```
Copy the output link to your browser, and you'll be able to run our code successfully and reproduce our results.
***
## Running code
With all environments successfully set up, you can run our code to reproduce our results. Next, we'll describe all the code in this repository by explaining its functions, input files, and output files one by one. You can run our interactive code line by line to obtain the output for each section. Of course, each code line is thoroughly commended to help users understand and execute.
The input files are referenced by relative paths in the code. If any issues occur with the paths, you can manually revise them according to the following descriptions.
***
#### VS generation
##### *VS_creat_S3 - global.mlx*:
+ Function: This code is used to generate Sentinel-3 VSs across global rivers.
+ Input files: 
  (1) Continental and national boundaries (shapefile) in `./VS generation/countries`;
  (2) SWORD nodes (NetCDF) from [this web](http://gaia.geosci.unc.edu/SWORD/#:~:text=SWORD%20was%20developed%20by%20the%20SWOT%20Science%20Team%20to%20serve); 
  (3) Sentinel-3A orbit buffer (shapefile) in `./VS generation/S3A_orbit_buffer`. 
  (4) Sentinel-3B orbit buffer (shapefile) in `./VS generation/S3B_orbit_buffer`. 
+ Output files: Sentinel-3 VSs for globe and each continent (shapefile).
***
#### Altimetry calculation
##### *S3_parfor_global_river_water_level.mlx*:
+ Function: This code uses **parallel computing** to retrieve global river water levels from Sentinel-3 altimetry data with [our improved algorithm](https://github.com/Fangchq/An-improved-waveform-retracking-method/tree/master) and constructs a dataset across global VSs.
+ Input files:
  (1) Sentinel-3 altimetry data (NetCDF) from [Copernicus Data Space Ecosystem](https://browser.dataspace.copernicus.eu/?zoom=7&lat=48.33434&lng=39.08386&themeId=DEFAULT-THEME&visualizationUrl=U2FsdGVkX1%2B3eWB92jhCC7e6IMp4Dvdl7O%2FDkuWhw%2FNf5vBHTa8FKQFxbQTtILNzMwKqkRegWoMObKk1%2FEibLnrG54SFeAgLzYan%2Bm2KLry8mRvLjgPP4JuWEnDJHdX1&datasetId=S2_L2A_CDAS);
  (2) Sentinel-3A orbit (shapefile) in `./Altimetry calculation/Sentinel3 orbit/Sentinel-3-Absolute-Ground-Tracks/S3A_shp`; 
  (3) Sentinel-3B orbit (shapefile) in `./Altimetry calculation/Sentinel3 orbit/Sentinel-3-Absolute-Ground-Tracks/S3B_shp`; 
  (3) Sentinel-3 VS buffers (shapefile) from the output of ***VS_creat_S3 - global.mlx***;
  (4) Wet seasons (NetCDF) at `./Altimetry calculation/wet_season.nc`. 
+ Output files: Global and continental river water level datasets (JSON).

###### Note that the global river water level dataset is stored in `Result datasets/Global water level dataset (raw)`, uploaded to [Zenodo](https://doi.org/10.5281/zenodo.14671453), and serves as the input for the following processes.
***
#### Post processing
##### *Dataset_analyze_globe_VS_scale.ipynb*:
+ Function: This code employs a stringent quality examination on the global river water level dataset generated in **Altimetry calculation**. It then establishes a well-filtered dataset used for analyzing fluctuations and change rates at the VS scale.
+ Input files: 
  (1) Initial water level dataset (JSON) in `Result datasets/Global water level dataset (raw)`, which is uploaded to [Zenodo](https://doi.org/10.5281/zenodo.14671453);
  (2) Continental boundaries (shapefile) in `./Post processing/Continental_boundaries`; 
  (3) Wet seasons (NetCDF) at `./Altimetry calculation/wet_season.nc`.
+ Output files (uploaded to [Zenodo](https://doi.org/10.5281/zenodo.14671453)): Processed global and continental river water level datasets (JSON), stored in `Result datasets/Global water level datasets (processed)`. The description of their attributes can be found in Table 1 of [our technical documentation](https://doi.org/10.5281/zenodo.14671453).
##### *Dataset_analyze_globe_basin_scale.ipynb*:
+ Function: This code is used to analyze water level fluctuations, change rates, extreme stage years, and seasonal variations at the basin scale, based on the dataset from ***Dataset_analyze_globe_VS_scale.ipynb***.
+ Input files: 
  (1) Processed global river water level dataset (JSON) in `Result datasets/Global water level datasets (processed)`, which is uploaded to [Zenodo](https://doi.org/10.5281/zenodo.14671453);
  (2) Continental boundaries (shapefile) in `./Post processing/Continental_boundaries`; 
  (3) Wet seasons (NetCDF) at `./Altimetry calculation/wet_season.nc`;
  (4) Major hydrological basins (shapefile) in `Result datasets/Major_hydrological_basins`.
+ Output files (uploaded to [Zenodo](https://doi.org/10.5281/zenodo.14671453)): 
  (1) Basin-level fluctuations and change rates (MAT file), stored at `Result datasets/Basin-level MAT files/basin_data_fluc.mat`;
  (2) Basin-level seasonality and seasonal variations (MAT file), stored at `Result datasets/Basin-level MAT files/basin_season_data.mat`;
  (3) Basin-level extreme stage years (MAT file), stored at `Result datasets/Basin-level MAT files/basin_extreme_year_75per.mat`.
  The attributes of these outputs are explained in Table 3 of [our technical documentation](https://doi.org/10.5281/zenodo.14671453).

###### Note that:
 1. The basin-level analysis in ***Dataset_analyze_globe_basin_scale.ipynb*** uses parallel computing, with the default number of cores set to six. You can adjust this value according to your PC's configuration.
 2. *mpfun.py*, *mpfun_season.py*, and *mp_extreme_year2.py* are function files used in ***Dataset_analyze_globe_basin_scale.ipynb*** to support basin-level analysis.

***
#### Dataset validation
##### *Validate_with_insitu.ipynb*:
+ Function: This code is used to validate our global river water level dataset against in-situ measurements in the United States.
+ Input files: 
  (1) Processed global river water level dataset (JSON) in `Result datasets/Global water level datasets (processed)`, which is uploaded to [Zenodo](https://doi.org/10.5281/zenodo.14671453);
  (2) In-situ water levels (JSON) provided by [USGS](https://waterdata.usgs.gov/nwis/current/?type=dailystage&group_key=huc_cd&search_site_no_station_nm=06349700&site_no_name_select=siteno), stored at `Result datasets/Validation dataset/USGS_gauge_stage_update_to_20240901.json` (updated to 2024-9-1).
  (3) Table for matched VSs and gauging stations (excel file), stored at `./Dataset validation/Matching table.xlsx`.
+ Output files (uploaded to [Zenodo](https://doi.org/10.5281/zenodo.14671453)): River water level dataset with validation results (JSON), stored in `Result datasets/Validation dataset`. The description of its attributes can be found in Table 2 of [our technical documentation](https://doi.org/10.5281/zenodo.14671453).
##### *Dataset_Monte_experiment.ipynb*:
+ Function: This code is used to quantify VS-level and basin-level uncertainties in fluctuations and change rates through Monte Carlo simulations. The uncertainties are caused by altimetric errors, derived from in-situ validation.
+ Input files: 
  (1) Processed global river water level dataset (JSON) in `Result datasets/Global water level datasets (processed)`, which is uploaded to [Zenodo](https://doi.org/10.5281/zenodo.14671453);
  (2) Continental boundaries (shapefile) in `./Post processing/Continental_boundaries`; 
  (3) Major hydrological basins (shapefile) in `Result datasets/Major_hydrological_basins`.

###### Note that: 
*Monte_Carlo_VS.py* and *Monte_Carlo_basin.py* are function files used in ***Dataset_Monte_experiment.ipynb*** to support Monte Carlo simulations. Parallel computing is also applied to improve calculation efficiency. You can adjust the number of cores according to your PC's configuration.
