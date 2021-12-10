from ctypes import *

dejavuc = CDLL("./libdejavu-api.so")   

print(dejavuc)
print(dejavuc.random_paths)  
