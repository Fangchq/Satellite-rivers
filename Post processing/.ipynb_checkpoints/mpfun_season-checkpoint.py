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

def AE(wse_year):
    measuremnet_number=len(wse_year)
    wse_sum=sum(wse_year)
    ae=0
    for item in wse_year:
        if not item==0:
            ae=ae-item/wse_sum*(np.log2(item/wse_sum))
        else:
            ae=ae-0
    return ae/np.log2(measuremnet_number)*100

def calculate_slope(date_collect,wse_collect):
    years=np.arange(2016,2025)
    dates_lim=[datetime.datetime(2016, 9, 1),datetime.datetime(2024, 11, 1)]
    dates_num=date2num(dates_lim)
    date_gaps=np.arange(dates_num[0],dates_num[1],30)
    date_records=[]
    stage_50=[]
    AE_record=[]
    seasonal_max_stage=[]
    seasonal_min_stage=[]
    seasonal_fluctuation=[]
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
    for year in years:
        index= (date_records>= date2num(datetime.datetime(year,1,1))) & (date_records< date2num(datetime.datetime(year+1,1,1)))
        if sum(index)>=8:
            seasonal_wse=stage_50_np[index]
            seasonal_year.append(year)
            seasonal_max_stage.append(np.max(seasonal_wse))
            seasonal_min_stage.append(np.min(seasonal_wse))
            seasonal_fluctuation.append(np.max(seasonal_wse)-np.min(seasonal_wse))
            seasonal_wse_adjusted=np.array(seasonal_wse)-np.min(seasonal_wse)
            AE_record.append(AE(seasonal_wse_adjusted))
    
    model=RLM(AE_record,sm.add_constant(seasonal_year),M=sm.robust.norms.TukeyBiweight())
    slope_AE=model.fit().params[1]
    model=RLM(seasonal_max_stage,sm.add_constant(seasonal_year),M=sm.robust.norms.TukeyBiweight())
    slope_max_stage=model.fit().params[1]
    model=RLM(seasonal_min_stage,sm.add_constant(seasonal_year),M=sm.robust.norms.TukeyBiweight())
    slope_min_stage=model.fit().params[1]
    model=RLM(seasonal_fluctuation,sm.add_constant(seasonal_year),M=sm.robust.norms.TukeyBiweight())
    slope_fluctuation=model.fit().params[1]
    
    return np.median(AE_record),slope_AE,slope_max_stage,slope_min_stage,slope_fluctuation

# 第二个单元格：定义处理函数
def process_single_basin(basin_ID):
    """处理单个盆地的函数"""
    try:
        # 使用全局变量
        Basin_IDs = global_vars['Basin_IDs']
        china_basins = global_vars['china_basins']
        lons = global_vars['lons']
        lats = global_vars['lats']
        new_dataset2 = global_vars['new_dataset2']
        
        time_basin_collect = []
        wse_basin_collect = []
        basin_name = None
        count=0

        for search_index,search_item in enumerate(Basin_IDs):
            if search_item==basin_ID:
                china_basin=china_basins.iloc[search_index]['geometry']
                in_basin=inpolygon(lons,lats,china_basin)
                if sum(in_basin)>0:
                    time_basin = [new_dataset2[in_basin_i]['monitor_time'] 
                                for in_basin_i, in_basin_logical in enumerate(in_basin) if (in_basin_logical) & (new_dataset2[in_basin_i]['fluctuation2']<25)]
                    wse_basin = [new_dataset2[in_basin_i]['water_level_anomalies'] 
                               for in_basin_i, in_basin_logical in enumerate(in_basin) if (in_basin_logical) & (new_dataset2[in_basin_i]['fluctuation2']<25)]
                    if len(time_basin)==0:
                        continue

                    count+=len(time_basin)
                    time_basin_collect.append(np.concatenate(time_basin))
                    wse_basin_collect.append(np.concatenate(wse_basin))
                    
        
        basin_name = china_basins[china_basins['MAJ_BAS']==basin_ID].iloc[0]['MAJ_NAME']

        if count > 35:
            time_basin_collect = np.concatenate(time_basin_collect)
            wse_basin_collect = np.concatenate(wse_basin_collect)
            
            time_basin_collect_dates = [excel_num_to_date(date) for date in time_basin_collect]
            time_basin_collect_date_np = np.array(date2num(time_basin_collect_dates))

            AE_record,slope_AE,slope_max_stage,slope_min_stage,slope_fluctuation = calculate_slope(
                time_basin_collect_date_np, wse_basin_collect)

            return {
                'basin_name': basin_name,
                'AE_median': AE_record,
                'AE_slope': slope_AE,
                'year_max_slope': slope_max_stage,
                'year_min_slope': slope_min_stage,
                'year_fluctuation_slope': slope_fluctuation,
            }
        return None
        
    except Exception as e:
        print(f"Error processing basin {basin_ID}: {str(e)}")
        return None