from scipy.io import loadmat
import numpy as np
m = loadmat(r'D:\mathlab\github_stsi\STSI\Grid_Location_Parameters\Example_Patient1\LFD_Fxd_Vertices.mat')
print("keys:", [k for k in m.keys() if not k.startswith('__')])
tri = np.asarray(m['currytri'])
loc = np.asarray(m['curryloc'])
print("currytri shape:", tri.shape, tri.dtype)
print("curryloc shape:", loc.shape, loc.dtype)