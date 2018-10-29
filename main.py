import numpy as np
import scipy.io as sio

np.random.seed(1)
S = sio.loadmat('Lymph.mat')

lambda_var = 0.5
params_epsilon = 1e-2
mode = 'trace'
msg_urgency = 2
max_iter = 20

