import ctypes

quic = ctypes.CDLL('./CPP_Code/QUIC_utils.o')
quic.main.argtypes = (ctypes.POINTER(ctypes.c_char), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int))


def main(fx_info, trace, maxIter, method):
    return quic.main(fx_info, trace, maxIter, method)

