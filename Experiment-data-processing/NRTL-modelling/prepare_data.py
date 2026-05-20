import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = pd.read_csv('data/nrtl_data.csv')
T = data['T'].to_numpy()
x = data[['x1', 'x2', 'x3']].to_numpy()
y = data[['y1', 'y2', 'y3']].to_numpy()
gamma = data[['gamma1', 'gamma2', 'gamma3']].to_numpy()


np.savez('data/nrtl_prepared.npz', T=T, x=x, y=y, gamma=gamma)
print("Сохранено в data/nrtl_prepared.npz")
