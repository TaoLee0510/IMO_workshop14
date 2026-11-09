"""Release completed ODE solver cycles without changing numerical operations."""
import gc
import ctypes
import sys
import model

# Return unused allocator pages to macOS after the solver's Python cycles have
# been collected. This does not touch live arrays or numerical operations.
_pressure_relief=None
if sys.platform=='darwin':
    try:
        _libc=ctypes.CDLL('/usr/lib/libSystem.B.dylib')
        _pressure_relief=_libc.malloc_zone_pressure_relief
        _pressure_relief.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
        _pressure_relief.restype=ctypes.c_size_t
    except (OSError,AttributeError):pass

def release():
    gc.collect()
    if _pressure_relief is not None:_pressure_relief(None,0)
def install():
    if getattr(model,'_tmn_memory_control',False):return
    model._tmn_memory_control=True
    original_rescale=model.ForwardFilter._rescale
    original_patient=model.patient_loglik
    counter=[0]
    def rescale(self):
        answer=original_rescale(self);counter[0]+=1
        if counter[0]%8==0:gc.collect()
        return answer
    def patient(*args,**kwargs):
        try:return original_patient(*args,**kwargs)
        finally:release()
    model.ForwardFilter._rescale=rescale
    model.patient_loglik=patient
    gc.set_threshold(200,5,5)
install()
