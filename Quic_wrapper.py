import ctypes

quic = ctypes.CDLL('./CPP_Code/QUIC.o')
quic.main.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_void_p), ctypes.c_int, ctypes.POINTER(ctypes.c_void_p))


def main(nlhs, plhs, nrhs, prhs):
    return quic.main(nlhs, plhs, nrhs, prhs)

