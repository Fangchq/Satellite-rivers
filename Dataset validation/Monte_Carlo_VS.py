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

def init_worker(shared_data):
    """Initialize global variables for worker processes"""
    global global_vars
    global_vars = shared_data

def simulation(sim_num):
    """Perform Monte Carlo simulation on water level data
   
   Parameters:
   - sim_num: Simulation iteration number
   
   Performs:
   1. Adds random measurement errors to water levels
   2. Recalculates fluctuation and change rates
   3. Computes differences from original values
   
   Returns:
   - Dictionary with simulation results
   """
    new_dataset = global_vars['new_dataset']
    meas_error = global_vars['meas_error']
    fluctuation_np = np.array([item['fluctuation2'] for item in new_dataset])
    change_rate_np = np.array([item['change_rate'] for item in new_dataset])
    # Process each station
    for i, item in enumerate(new_dataset):
        # Get original measurements
        monitor_time=item['monitor_time']
        monitor_wse=item['monitor_wse']
        monitor_time_np=np.array(monitor_time)
        monitor_wse_np=np.array(monitor_wse)
        # Add random measurement errors
        random_errors=np.random.normal(0,meas_error,size=monitor_wse_np.shape)
        simulated_wse=monitor_wse_np+random_errors
        # Calculate simulated metrics
        fluctuation_simulated=np.max(simulated_wse)-np.min(simulated_wse)
        new_dataset[i]['fluctuation_simulated']=fluctuation_simulated
        model=RLM(simulated_wse,sm.add_constant(monitor_time_np),M=sm.robust.norms.TukeyBiweight())
        results=model.fit()
        new_dataset[i]['change_rate_simulated'] = results.params[1] * 365#m/year
        new_dataset[i]['p_value_for_change_simulated']=results.pvalues[1]
        # Calculate differences from original values
        new_dataset[i]['fluctuation_diff_simulated']=item['fluctuation2']-fluctuation_simulated
        new_dataset[i]['change_diff_simulated']=item['change_rate']-results.params[1] * 365
    # Collect results
    flu_diff=[item['fluctuation_diff_simulated'] for item in new_dataset]
    change_rate_diff=[item['change_diff_simulated'] for item in new_dataset]
    return {
        'simulate_index': sim_num,
        'flu_diff':flu_diff,
        'change_rate_diff': change_rate_diff,
    }