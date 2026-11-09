"""Optional compiled implementation of the unchanged forward derivative."""
import ctypes,os
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
_path=ROOT/'bin/libtmn_rhs.dylib'
_lib=ctypes.CDLL(str(_path)) if _path.exists() else None
enabled=_lib is not None and os.environ.get('TMN_FORCE_PYTHON_RHS')!='1'
if _lib is not None:
    _function=_lib.tmn_forward_rhs
    _function.argtypes=[ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)]+[ctypes.c_double]*10+[ctypes.c_int]
    _function.restype=None
def forward_rhs(state,low,N,b,mu,fc,fd,q,up,down,eta_dot,killed):
    if not enabled:return None
    n=len(state)-int(killed);out=np.empty_like(state)
    _function(n,state.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
              low,N,b,mu,fc,fd,q,up,down,eta_dot,int(killed))
    return out
