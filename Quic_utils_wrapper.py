from ctypes import cdll
lib = cdll.LoadLibrary('./CPP_Code/QUIC_utils.so')

print(lib)
