"""Risk-only acceleration with a certified absolute CDF remainder bound.

The fit always retains the full live survival/density calculation. Here a
forecast can stop once remaining event-free mass is below 1e-14 of accumulated
event mass: all later CDFs differ by at most that bound. No small risk is floored.
"""
import math,numpy as np
from model import ForwardFilter,Backgrounds,tilt_for_horizon,rates
class RiskFilter(ForwardFilter):
    risk_early_stop=False
    risk_stop_time=None
    def propagate(self,start,end):
        if not self.risk_early_stop:return super().propagate(start,end)
        if self.risk_stop_time is not None:return
        t=float(start)
        while t<end:
            birth,death=rates(np.array([float(self.cap)]),self.qfun(t),t,self.alpha,self.theta,self.p)
            decay=max(0.,math.sqrt(float(birth[0]))-math.sqrt(float(death[0])))**2
            step=min(120.,10/max(decay,1e-12));stop=min(float(end),t+step)
            super().propagate(t,stop);t=stop
            if self.physical_logmass()-self.log_cdf()<math.log(1e-14):
                self.risk_stop_time=t;break
def risk_curve(alpha,theta,p,founder_fraction,baseline_day=0.,targets=(365,730,1095,1460),refinement=2):
    if p.D_init>=p.K:return {float(t):1. for t in targets}
    if p.mu==0 and p.D_init==0:return {float(t):0. for t in targets}
    row={'baseline_vaf':founder_fraction/2,'baseline_day':baseline_day,'observations':[],'outcome':{'kind':'unknown'}}
    bg=Backgrounds([row],alpha,theta,p,refinement);result={}
    for target in targets:
        f=RiskFilter(p,alpha,theta,bg.for_patient(0),killed=True,cap=p.K-1,eta=tilt_for_horizon(alpha,theta,p,target,baseline_day),refinement=refinement)
        if baseline_day<0:f.propagate(baseline_day,0.);f.condition_event_free()
        f.risk_early_stop=True;f.propagate(0.,target)
        total=f.physical_logmass(include_absorbed=True)
        if abs(total)>2e-5:raise ArithmeticError(f'Risk mass drift {total}')
        lp=f.log_cdf()-total;result[float(target)]=math.exp(lp) if lp>-745 else 0.
    return result
