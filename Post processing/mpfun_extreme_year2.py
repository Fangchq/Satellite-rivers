import multiprocessing as mp
from tqdm.notebook import tqdm
import warnings
import numpy as np
from matplotlib.dates import date2num, num2date, DateFormatter, AutoDateLocator
import datetime
import geopandas as gpd
from statsmodels.robust.robust_linear_model import RLM
import statsmodels.api as sm
from matplotlib.ticker import MultipleLocator
import pandas as pd
from shapely.geometry import Point, Polygon, shape
from matplotlib import rcParams
from scipy.stats import gaussian_kde, linregress
from math import ceil, sqrt
import cartopy.feature as cfeature
from statsmodels.tsa.seasonal import STL
import pymannkendall as mk
from shapely.ops import unary_union
import xarray as xr
import matplotlib.colors as mcolors
import matplotlib.cm as cm
# Function to convert Excel serial date number to datetime 
def excel_num_to_date(excel_num):
    excel_epoch = datetime.datetime(1899, 12, 30)
    delta = datetime.timedelta(days=excel_num)
    return excel_epoch + delta
#Function to initialize global variables for worker processes
def init_worker(shared_data):
    global global_vars
    global_vars = shared_data
# Function to check if points are inside a polygon
def inpolygon(lons, lats, basin_polygon):
    points = [Point(lon, lat) for lon,lat in zip(lons,lats)]
    return [basin_polygon.contains(point) for point in points]
# Function to detemine the start and end months of wet seasons for each VS
def search_start_and_end_month(lon,lat,wet_season):
    lat_index=abs(wet_season['lat']-lat).argmin()
    lon_index=abs(wet_season['lon']-lon).argmin()
    start_month=wet_season['start_month'][lat_index,lon_index]
    end_month=wet_season['end_month'][lat_index,lon_index] 
    median_month=wet_season['median_month'][lat_index,lon_index]
    return start_month,end_month,median_month
# Function to identify the low and high stage years for each basin in the whole year analysis
def calculate_whole_year(date_collect,wse_collect):
    years=np.arange(2016,2025)
    dates_lim=[datetime.datetime(2016, 9, 1),datetime.datetime(2024, 11, 1)]
    dates_num=date2num(dates_lim)
    date_gaps=np.arange(dates_num[0],dates_num[1],30)
    date_records=[]
    stage_50=[]
    seasonal_year=[]
    for date_i in range(len(date_gaps)-1):
        date_record=(date_gaps[date_i]+date_gaps[date_i+1])/2
        date_mask= (date_collect>=date_gaps[date_i]) & (date_collect<date_gaps[date_i+1])
        if sum(date_mask)>2:
            date_records.append(date_record)
            wse_selected= wse_collect[date_mask]
            stage_50.append(np.median(wse_selected))
    stage_50_np=np.array(stage_50)
    #date_records_date_np=np.array(num2date(date_records))
    # ddetemine the thresholds of extreme low and high stages
    upper_wse=np.percentile(stage_50_np,75)
    upper_percents=[]
    lower_wse=np.percentile(stage_50_np,25)
    lower_percents=[]
    # calculate the percentages of extreme low and high stages within a year
    for year in years[::-1]:
        index= (date_records>= date2num(datetime.datetime(year,1,1))) & (date_records< date2num(datetime.datetime(year+1,1,1)))
        # ensure at least 7 monthly water levels
        if sum(index)>=7:
            seasonal_wse=stage_50_np[index]
            seasonal_year.append(year)
            upper_percents.append(sum(seasonal_wse>=upper_wse)/len(seasonal_wse))
            lower_percents.append(sum(seasonal_wse<=lower_wse)/len(seasonal_wse))
    # ensure at least 7 years of valid percentages
    if len(upper_percents)>=7:
        max_indices=np.argmax(upper_percents)
        Basin_collect_high_year_percent_most_whole_year=seasonal_year[max_indices]
        min_indices=np.argmax(lower_percents)
        Basin_collect_low_year_percent_most_whole_year=seasonal_year[min_indices]
    
        return Basin_collect_high_year_percent_most_whole_year, Basin_collect_low_year_percent_most_whole_year
    return -1,-1
    
# Function to identify the low and high stage years for each basin in the wet season analysis
def calculate_wet_season(date_collect,wse_collect,start_month,end_month):
    years=np.arange(2016,2025)
    dates_lim=[datetime.datetime(2016, 9, 1),datetime.datetime(2024, 11, 1)]
    dates_num=date2num(dates_lim)
    date_gaps=np.arange(dates_num[0],dates_num[1],30)
    date_records=[]
    stage_50=[]
    seasonal_year=[]
    summer_wse=[]
    summer_dates=[]
    # construct basin-level stage time series
    for date_i in range(len(date_gaps)-1):
        date_record=(date_gaps[date_i]+date_gaps[date_i+1])/2
        date_mask= (date_collect>=date_gaps[date_i]) & (date_collect<date_gaps[date_i+1])
        # ensure more than 2 valid measurements in the range window
        if sum(date_mask)>2:
            date_records.append(date_record)
            wse_selected= wse_collect[date_mask]
            stage_50.append(np.median(wse_selected))
    stage_50_np=np.array(stage_50)
    # extract wet season water levels
    for date_i,date in enumerate(date_records):
        if start_month>8:
            if (num2date(date).month>=start_month) | (num2date(date).month<=end_month):
                summer_wse.append(stage_50[date_i])
                summer_dates.append(date)
        else:
            if (num2date(date).month>=start_month) & (num2date(date).month<=end_month):
                summer_wse.append(stage_50[date_i])
                summer_dates.append(date)
    upper_wse=np.percentile(summer_wse,75)
    upper_percents=[]
    lower_wse=np.percentile(summer_wse,25)
    lower_percents=[]
    summer_years=[num2date(date).year for date in summer_dates]
    summer_year_np=np.array(summer_years)
    summer_wse_np=np.array(summer_wse)
    
    for year in years[::-1]:
        index= summer_year_np==year
        # ensure at least 3 monthly water levels in the wet season
        if sum(index)>=3:
            seasonal_wse=summer_wse_np[index]
            seasonal_year.append(year)
            upper_percents.append(sum(seasonal_wse>=upper_wse)/len(seasonal_wse))
            lower_percents.append(sum(seasonal_wse<=lower_wse)/len(seasonal_wse))
    
    if len(upper_percents)>=7:
        max_indices=np.argmax(upper_percents)# the year with the highest percentage of high stages
        Basin_collect_high_year_percent_most_wet_season=seasonal_year[max_indices]
        min_indices=np.argmin(upper_percents)# the year with the lowest percentage of high stages
        Basin_collect_low_year_percent_most_wet_season=seasonal_year[min_indices]
    
        return Basin_collect_high_year_percent_most_wet_season, Basin_collect_low_year_percent_most_wet_season
    return -1,-1

# Function to process a single basin
def process_single_basin(basin_ID):
    try:
        # Use global variables
        Basin_IDs = global_vars['Basin_IDs']
        china_basins = global_vars['china_basins']
        new_dataset2 = global_vars['new_dataset2']
        wet_season= global_vars['wet_season']
        lons=[item['center_lon'] for item in new_dataset2]
        lats=[item['center_lat'] for item in new_dataset2]

        Basin_collect_low_year_percent_most_whole_year=[]
        Basin_collect_high_year_percent_most_whole_year=[]
        Basin_collect_low_year_percent_most_wet_season=[]
        Basin_collect_high_year_percent_most_wet_season=[]
        time_basin_collect = []
        wse_basin_collect = []
        year_index=np.arange(2016,2025,1)
        basin_name = None
        count=0
        start_month_basin=[]
        end_month_basin=[]
        median_month_basin=[]

        for search_index, search_item in enumerate(Basin_IDs):
            if search_item == basin_ID:
                china_basin = china_basins.iloc[search_index]['geometry']
                in_basin = inpolygon(lons, lats, china_basin)
                # collect all stage measurements for VSs within the basin
                if sum(in_basin)>0:
                    time_basin = [new_dataset2[in_basin_i]['monitor_time'] 
                                for in_basin_i, in_basin_logical in enumerate(in_basin) if (in_basin_logical) & (new_dataset2[in_basin_i]['fluctuation2']<25)]
                    wse_basin = [new_dataset2[in_basin_i]['water_level_anomalies'] 
                               for in_basin_i, in_basin_logical in enumerate(in_basin) if (in_basin_logical) & (new_dataset2[in_basin_i]['fluctuation2']<25)]
                    time_basin_collect.append(np.concatenate(time_basin))
                    wse_basin_collect.append(np.concatenate(wse_basin))
                    if len(time_basin)==0:
                        continue

                    count+=len(time_basin)
                    
                    for in_basin_i, in_basin_logical in enumerate(in_basin):
                        if (in_basin_logical) & (new_dataset2[in_basin_i]['fluctuation2']<25):
                            start_month,end_month,median_month=search_start_and_end_month(lats[in_basin_i],lons[in_basin_i],wet_season)
                            start_month_basin.append(start_month)
                            end_month_basin.append(end_month)
                            median_month_basin.append(median_month)
                        
                    
        basin_name = china_basins[china_basins['MAJ_BAS']==basin_ID].iloc[0]['MAJ_NAME']
        # Ensure the basin has more than 35 valid VSs 
        if count > 35:
            month_basin_median=np.median(np.array(median_month_basin))
            if month_basin_median<3:
                month_basin_start=month_basin_median+12-2
                month_basin_end=month_basin_median+2
            elif month_basin_median>10:
                month_basin_start=month_basin_median-2
                month_basin_end=month_basin_median+2-12
            else:
                month_basin_start=month_basin_median-2
                month_basin_end=month_basin_median+2
                
            time_basin_collect = np.concatenate(time_basin_collect)
            wse_basin_collect = np.concatenate(wse_basin_collect)
            
            time_basin_collect_dates = [excel_num_to_date(date) for date in time_basin_collect]
            time_basin_collect_date_np = np.array(date2num(time_basin_collect_dates))
            # identify extreme stage years
            Basin_collect_high_year_percent_most_whole_year, Basin_collect_low_year_percent_most_whole_year=calculate_whole_year(time_basin_collect_date_np,wse_basin_collect)
            Basin_collect_high_year_percent_most_wet_season, Basin_collect_low_year_percent_most_wet_season=calculate_wet_season(time_basin_collect_date_np,wse_basin_collect,month_basin_start,month_basin_end)
            # -1 represents no valid outputs, as the basin did not pass the quality examination
            if (Basin_collect_high_year_percent_most_whole_year>0) & (Basin_collect_high_year_percent_most_wet_season>0):
                return {
                    'basin_name': basin_name,
                    'low_year_basin_whole_year': Basin_collect_low_year_percent_most_whole_year,
                    'high_year_basin_whole_year': Basin_collect_high_year_percent_most_whole_year,
                    'low_year_basin_wet_season': Basin_collect_low_year_percent_most_wet_season,
                    'high_year_basin_wet_season': Basin_collect_high_year_percent_most_wet_season,
                }
            elif Basin_collect_high_year_percent_most_whole_year>0:
                return {
                    'basin_name': basin_name,
                    'low_year_basin_whole_year': Basin_collect_low_year_percent_most_whole_year,
                    'high_year_basin_whole_year': Basin_collect_high_year_percent_most_whole_year,
                    'low_year_basin_wet_season': -1,
                    'high_year_basin_wet_season': -1,
                }
            elif Basin_collect_high_year_percent_most_wet_season>0:
                return {
                    'basin_name': basin_name,
                    'low_year_basin_whole_year': -1,
                    'high_year_basin_whole_year': -1,
                    'low_year_basin_wet_season': Basin_collect_low_year_percent_most_wet_season,
                    'high_year_basin_wet_season': Basin_collect_high_year_percent_most_wet_season,
                }
            else:
                return None       
        #return None
    except Exception as e:
        print(f"Error processing basin {basin_ID}: {str(e)}")
        return None