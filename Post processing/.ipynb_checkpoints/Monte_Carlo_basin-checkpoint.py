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

def excel_num_to_date(excel_num):
    excel_epoch = datetime.datetime(1899, 12, 30)  # Excel的纪元开始
    delta = datetime.timedelta(days=excel_num)
    return excel_epoch + delta

def init_worker(shared_data):
    """初始化工作进程的全局变量"""
    global global_vars
    global_vars = shared_data

def inpolygon(lons, lats, basin_polygon):
    points = [Point(lon, lat) for lon,lat in zip(lons,lats)]
    return [basin_polygon.contains(point) for point in points]

def calculate_slope(date_collect,wse_collect):
    dates_lim=[datetime.datetime(2016, 9, 1),datetime.datetime(2024, 11, 1)]
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

# 第二个单元格：定义处理函数
def process_single_basin(basin_ID):
    """处理单个盆地的函数"""
    try:
        # 使用全局变量
        Basin_IDs = global_vars['Basin_IDs']
        meas_error = global_vars['meas_error']
        china_basins = global_vars['china_basins']
        new_dataset2 = global_vars['new_dataset']
        lons=[item['center_lon'] for item in new_dataset2]
        lats=[item['center_lat'] for item in new_dataset2]
        fluctuation_np=np.array([item['fluctuation2'] for item in new_dataset2])
        change_rate_np=np.array([item['change_rate'] for item in new_dataset2])
        
        fluctuation_basin = []
        change_rate_basin = []
        time_basin_collect = []
        wse_basin_collect = []
        basin_name = None
        count=0
        
        for search_index, search_item in enumerate(Basin_IDs):
            if search_item == basin_ID:
                china_basin = china_basins.iloc[search_index]['geometry']
                in_basin = inpolygon(lons, lats, china_basin)
                
                if sum(in_basin) > 0:
                    if sum(in_basin) == 1:
                        fluctuation_basin.extend([fluctuation_np[in_basin][0]])
                        change_rate_basin.extend([change_rate_np[in_basin][0]])
                    else:
                        fluctuation_basin.extend(fluctuation_np[in_basin])
                        change_rate_basin.extend(change_rate_np[in_basin])

                    in_dataset=[item 
                                for item_i,item in enumerate(new_dataset2) if (in_basin[item_i]) & (item['fluctuation2']<25)]
                    time_basin = [item['monitor_time'] 
                                for item in in_dataset]
                    wse_basin = [item['water_level_anomalies'] 
                               for item in in_dataset]
                    if len(time_basin)==0:
                        continue

                    count+=len(time_basin)
                    time_basin_collect.append(np.concatenate(time_basin))
                    wse_basin_collect.append(np.concatenate(wse_basin))
                    basin_name = china_basins[china_basins['MAJ_BAS']==basin_ID].iloc[0]['MAJ_NAME']

        if count > 35:
            time_basin_collect = np.concatenate(time_basin_collect)
            wse_basin_collect = np.concatenate(wse_basin_collect)
            fluc_test=[]
            change_rate_test=[]

            time_basin_collect_dates = [excel_num_to_date(date) for date in time_basin_collect]
            time_basin_collect_date_np = np.array(date2num(time_basin_collect_dates))
            slope_median_raw, _, fluctuation_raw, _ = calculate_slope(
                time_basin_collect_date_np, wse_basin_collect)
            
            for sim_index in range(100):
                random_errors=np.random.normal(0,meas_error,size=wse_basin_collect.shape)
                wse_basin_collect_sim=wse_basin_collect+random_errors
            
                time_basin_collect_dates = [excel_num_to_date(date) for date in time_basin_collect]
                time_basin_collect_date_np = np.array(date2num(time_basin_collect_dates))
    
                slope_median1, slope_median2, fluctuation1, fluctuation2 = calculate_slope(
                    time_basin_collect_date_np, wse_basin_collect_sim)
                fluc_test.append(fluctuation1)
                change_rate_test.append(slope_median1)

            return {
                'basin_name': basin_name,
                'fluctuation_sim': fluc_test,
                'change_rate_sim': change_rate_test,
                'fluctuation_raw':fluctuation_raw,
                'change_rate_raw':slope_median_raw,
            }
        return None
    except Exception as e:
        print(f"Error processing basin {basin_ID}: {str(e)}")
        return None