import numpy as np
import scipy.io as sio
import Quic_wrapper
import ctypes

np.random.seed(1)
S = sio.loadmat('Lymph.mat')

lambda_var = 0.5
params_epsilon = 1e-2
mode = 'trace'
msg_urgency = 2
max_iter = 20

sData = S["data"]


def convertToCharArray(str):
    return ctypes.cast(ctypes.c_char_p(str.encode('utf-8')), ctypes.c_void_p)

def convertToInt(number):
    return ctypes.cast(ctypes.c_char_p(str(number).encode('utf-8')), ctypes.c_void_p)


rhs = [convertToCharArray(mode),
       ctypes.cast(sData.ctypes.data_as(ctypes.c_char_p), ctypes.c_void_p),
       convertToInt(lambda_var), convertToInt(params_epsilon),
       convertToInt(msg_urgency), convertToInt(max_iter),
       convertToInt(sData.shape[0]), convertToInt(sData.shape[1])]
lhs1 = convertToInt(0)
lhs2 =  convertToInt(0)
lhs3 = convertToInt(0)
lhs4 = convertToInt(0)
lhs5 = convertToInt(0)
lhs6 = convertToInt(0)
lhs7 = convertToInt(0)
lhs = [lhs1, lhs2, lhs3, lhs4, lhs5, lhs6, lhs7]

rhsArr = (ctypes.c_void_p * len(rhs))(*rhs)
lhsArr = (ctypes.c_void_p * len(lhs))(*lhs)
res = Quic_wrapper.main(len(lhs),
                        lhsArr,
                        len(rhs),
                        rhsArr)
print(res)
print(lhs[0])
