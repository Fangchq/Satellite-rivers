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
# Function to calculate the slope and fluctuations based on the developed basin-level time series
def calculate_slope(date_collect,wse_collect):
    dates_lim=[datetime.datetime(2016, 9, 1),datetime.datetime(2024, 11, 1)]
    dates_num=date2num(dates_lim)
    date_gaps=np.arange(dates_num[0],dates_num[1],30)
    date_records=[]
    stage_50=[]
    # construct basin-level stage time series by capturing the median level within the monthly sliding window
    for date_i in range(len(date_gaps)-1):
        date_record=(date_gaps[date_i]+date_gaps[date_i+1])/2
        date_mask= (date_collect>=date_gaps[date_i]) & (date_collect<date_gaps[date_i+1])
        if sum(date_mask)>2:
            date_records.append(date_record)
            wse_selected= wse_collect[date_mask]
            stage_50.append(np.median(wse_selected))
    # Use Robust Linear Model (RLM) to calculate slope
    model=RLM(stage_50,sm.add_constant(date_records),M=sm.robust.norms.TukeyBiweight())
    slope_median=model.fit().params[1]*365
    fluctuation1=max(stage_50)-min(stage_50)

    dates_lim=[datetime.datetime(2019, 11, 1),datetime.datetime(2024, 11, 1)]
    dates_num=date2num(dates_lim)
    date_gaps=np.arange(dates_num[0],dates_num[1],30)
    date_records=[]
    stage_50=[]
    for date_i in range(len(date_gaps)-1):
        date_record=(date_gaps[date_i]+date_gaps[date_i+1])/2
        date_mask= (date_collect>=date_gaps[date_i]) & (date_collect<date_gaps[date_i+1])
        if sum(date_mask)>2:
            date_records.append(date_record)
            wse_selected= wse_collect[date_mask]
            stage_50.append(np.median(wse_selected))
    
    model=RLM(stage_50,sm.add_constant(date_records),M=sm.robust.norms.TukeyBiweight())
    slope_median2=model.fit().params[1]*365
    fluctuation2=max(stage_50)-min(stage_50)
    return slope_median,slope_median2,fluctuation1,fluctuation2

# Function to process a single basin
def process_single_basin(basin_ID):
    try:
        # Use global variables
        Basin_IDs = global_vars['Basin_IDs']
        china_basins = global_vars['china_basins']
        lons = global_vars['lons']
        lats = global_vars['lats']
        fluctuation_np = global_vars['fluctuation_np']
        change_rate_np = global_vars['change_rate_np']
        new_dataset2 = global_vars['new_dataset2']
        
        fluctuation_basin = []
        change_rate_basin = []
        time_basin_collect = []
        wse_basin_collect = []
        basin_name = None

        for search_index, search_item in enumerate(Basin_IDs):
            if search_item == basin_ID:
                china_basin = china_basins.iloc[search_index]['geometry']
                in_basin = inpolygon(lons, lats, china_basin)
                # collect all stage measurements for VSs within the basin
                if sum(in_basin) > 0:
                    if sum(in_basin) == 1:
                        fluctuation_basin.extend([fluctuation_np[in_basin][0]])
                        change_rate_basin.extend([change_rate_np[in_basin][0]])
                    else:
                        fluctuation_basin.extend(fluctuation_np[in_basin])
                        change_rate_basin.extend(change_rate_np[in_basin])

                    in_fluctuation_index=np.array(fluctuation_basin)<25
                    
                    time_basin = [new_dataset2[in_basin_i]['monitor_time'] 
                                for in_basin_i, in_basin_logical in enumerate(in_basin) if (in_basin_logical) & (new_dataset2[in_basin_i]['fluctuation2']<25)]
                    wse_basin = [new_dataset2[in_basin_i]['water_level_anomalies'] 
                               for in_basin_i, in_basin_logical in enumerate(in_basin) if (in_basin_logical) & (new_dataset2[in_basin_i]['fluctuation2']<25)]
                    if len(time_basin)==0:
                        continue
                    time_basin_collect.append(np.concatenate(time_basin))
                    wse_basin_collect.append(np.concatenate(wse_basin))
                    basin_name = china_basins[china_basins['MAJ_BAS']==basin_ID].iloc[0]['MAJ_NAME']
        # Ensure the basin has more than 35 valid VSs 
        if sum(in_fluctuation_index) > 35:
            time_basin_collect = np.concatenate(time_basin_collect)
            wse_basin_collect = np.concatenate(wse_basin_collect)
            
            time_basin_collect_dates = [excel_num_to_date(date) for date in time_basin_collect]
            time_basin_collect_date_np = np.array(date2num(time_basin_collect_dates))
            # Calculate slope and fluctuation for the collected data
            slope_median1, slope_median2, fluctuation1, fluctuation2 = calculate_slope(
                time_basin_collect_date_np, wse_basin_collect)

            return {
                'basin_name': basin_name,
                'fluctuation_median': np.median(fluctuation_basin),
                'change_rate_median': np.median(change_rate_basin),
                'fluctuation': fluctuation_basin,
                'change_rate': change_rate_basin,
                'fluctuation2': fluctuation1,
                'change_rate2': slope_median1
            }
        return None
    except Exception as e:
        print(f"Error processing basin {basin_ID}: {str(e)}")
        return None