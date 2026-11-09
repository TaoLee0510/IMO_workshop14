"""Fixed-area gamma model and exponentially tilted shared-state forward filter.

Inference retains the explicit deterministic q=C/(N-D) approximation. Unknown
outcomes use an unrestricted (adaptively truncated, reflecting) D chain, never
the clinical killed chain. Event probabilities use exponential tilting so tiny
early-event flux is not replaced by a probability floor.
"""
from dataclasses import dataclass,replace
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit,logit,gammaln,logsumexp,gammainc


@dataclass(frozen=True)
class Parameters:
    N:int=100000
    b:float=1/28
    s_D:float=.02
    mu:float=.01
    area_days:float=1.
    D_init:int=0
    threshold_fraction:float=.01
    sigma_vaf:float=.25
    @property
    def K(self):return math.floor(self.threshold_fraction*self.N)+1


def fitness(t,alpha,theta,p):
    if alpha<=1 or theta<=0:raise ValueError('alpha>1 and theta>0 required')
    if t<0:return 1.,1.
    extra=0. if t==0 else p.area_days*math.exp((alpha-1)*math.log(t)-t/theta-gammaln(alpha)-alpha*math.log(theta))
    return 1+extra,1+extra+p.s_D


def peak_extra(alpha,theta,p):
    return fitness((alpha-1)*theta,alpha,theta,p)[0]-1


def tilt_for_horizon(alpha,theta,p,end,baseline_day):
    # A constant-rate branching approximation chooses a numerical tilt only.
    # The actual finite-N, time-dependent rates remain unchanged in the filter.
    duration=max(1.,end-baseline_day)
    area_to_end=p.area_days*float(gammainc(alpha,max(0.,end)/theta))
    effective_sd=(p.s_D*max(0.,end)+area_to_end)/duration
    birth=p.b*(1+effective_sd);death=p.b
    if effective_sd<1e-10:q=birth*duration/(1+birth*duration)
    else:
        decay=math.exp(-(birth-death)*duration)
        q=birth*(1-decay)/(birth-death*decay)
    return 1.01*(-math.log(max(q,1e-12)))


class Backgrounds:
    def __init__(self,cohort,alpha,theta,p,refinement=1):
        self.cohort=cohort;self.p=p
        q0=[]
        for row in cohort:
            founder=2*row['baseline_vaf'];y=p.D_init/p.N
            if not 0<=y<=founder<1 or p.D_init>=p.K:raise ValueError('Incompatible baseline/D_init')
            q=(founder-y)/(1-y)
            q0.append(q)
        self.initial=np.array(q0)
        self.zero=expit(logit(self.initial)+p.b*p.mu*np.array([r['baseline_day'] for r in cohort]))
        end=max([1460.]+[r['outcome'].get('day',0.) for r in cohort])
        pulse_end=min(end,(alpha-1)*theta+25*theta)
        fine_step=min(1/refinement,max(.005,theta/(8*refinement)))
        fine=np.arange(0,pulse_end,fine_step)
        coarse=np.arange(pulse_end,end,2/refinement)
        self.grid=np.unique(np.r_[fine,coarse,0.,end])
        def rhs(t,q):
            fm,_=fitness(t,alpha,theta,p)
            return p.b*q*(1-q)*((1-p.mu)*fm-1)/(1+(fm-1)*q)
        boundaries=sorted(set([0.,end]+[t for t in ((alpha-1)*theta,pulse_end) if 0<t<end]))
        pieces=[];state=self.zero
        for start,stop in zip(boundaries[:-1],boundaries[1:]):
            maxstep=min(10.,max(.01,theta/3)) if start<pulse_end else 30.
            sol=solve_ivp(rhs,(start,stop),state,method='DOP853',dense_output=True,rtol=2e-10,atol=1e-12,max_step=maxstep)
            if not sol.success:raise ArithmeticError(sol.message)
            pieces.append((start,stop,sol.sol));state=sol.y[:,-1]
        self.values=np.empty((len(cohort),len(self.grid)))
        for start,stop,sol in pieces:
            where=(self.grid>=start)&(self.grid<=stop)
            self.values[:,where]=sol(self.grid[where])
    def for_patient(self,index):
        base=self.cohort[index]['baseline_day'];q=self.initial[index];p=self.p
        def fun(t):
            if t<0:return float(expit(logit(q)-p.b*p.mu*(t-base)))
            return float(np.interp(t,self.grid,self.values[index]))
        return fun


def rates(d,q,t,alpha,theta,p):
    fm,fd=fitness(t,alpha,theta,p)
    remaining=p.N-d;c=remaining*q
    denom=remaining+(fm-1)*c+fd*d
    pd=(fd*d+p.mu*fm*c)/denom
    return p.b*remaining*pd,p.b*d*(1-pd)


class ForwardFilter:
    def __init__(self,p,alpha,theta,qfun,*,killed,cap,eta=0.,refinement=1,window=False,window_z=12.):
        self.p=p;self.alpha=alpha;self.theta=theta;self.qfun=qfun;self.killed=killed
        self.cap=p.K-1 if killed else min(p.N,int(cap))
        self.d=np.arange(self.cap+1,dtype=float)
        # A time-dependent branching tilt avoids amplifying unresolved tiny
        # early-time tails, which a fixed endpoint tilt can do in stiff solvers.
        self.dynamic_tilt=eta>0
        # Large initial D does not need a rare-seed tilt. Limit its initial
        # log weight so adaptive absolute tolerances cannot erase physical mass.
        # Immigration also amplifies a tilted distribution. Bound its initial
        # log moment when mu*C is large, while retaining strong tilt for rare seeds.
        immigration_shape=p.mu*p.N*max(0.,float(qfun(0.)))
        initial_eta=min(math.log(100.),30/(max(1,p.D_init)+immigration_shape)) if self.dynamic_tilt else 0.
        self.u=math.exp(-initial_eta)
        self.eta=-math.log(self.u);self.refinement=refinement;self.logscale=self.eta*p.D_init
        # Killed filters have a separate absorbing-flux state.
        self.state=np.zeros(self.cap+1+int(killed));self.state[p.D_init]=1.
        self.absorbed_logmass=-math.inf
        self.window=bool(window and not killed and not self.dynamic_tilt)
        self.window_z=float(window_z);self.window_max_size=len(self.d);self.window_edge_probability=0.
        self.negative_roundoff=0.;self.max_upper_probability=0.;self.nfev=0
    def _rescale(self):
        maximum=float(self.state.max())
        if maximum<=0 or not math.isfinite(maximum):raise ArithmeticError('Invalid tilted probability scale')
        minimum=float(self.state.min())
        if minimum < -max(1e-10*maximum,1e-25):raise ArithmeticError('Material negative forward probability')
        self.negative_roundoff=min(self.negative_roundoff,minimum/maximum)
        self.state=np.maximum(self.state,0.)/maximum
        self.logscale+=math.log(maximum)
    def physical_logmass(self,include_absorbed=False):
        size=len(self.d)
        values=self.state[:size]
        positive=values>0
        if not positive.any():return -math.inf
        indices=self.d
        live=self.logscale+float(logsumexp(np.log(values[positive])-self.eta*indices[positive]))
        return float(np.logaddexp(live,self.absorbed_logmass)) if include_absorbed and self.killed else live
    def condition_event_free(self):
        self.logscale-=self.physical_logmass()
        if self.killed:self.state[-1]=0.;self.absorbed_logmass=-math.inf
    def _move_window(self,start,end):
        if not self.window:return
        p=self.p;values=self.state;w=values/values.sum()
        m=float(w@self.d);v=float(w@((self.d-m)**2))
        def rhs(t,state):
            mean=float(np.clip(state[0],0,p.N));var=max(0.,float(state[1]));q=self.qfun(t)
            fm,fd=fitness(t,self.alpha,self.theta,p)
            den=p.N*(1+(fm-1)*q)+mean*((fd-1)-(fm-1)*q)
            num=p.mu*fm*q*p.N+mean*(fd-p.mu*fm*q)
            pd=num/den;dp=((fd-p.mu*fm*q)*den-num*((fd-1)-(fm-1)*q))/den**2
            birth=p.b*(p.N-mean)*pd;death=p.b*mean*(1-pd)
            return [birth-death,2*p.b*(p.N*dp-1)*var+birth+death]
        sol=solve_ivp(rhs,(start,end),[m,v],rtol=1e-5,atol=1e-6,max_step=max(.1,(end-start)/3))
        if not sol.success:raise ArithmeticError(sol.message)
        # Retain the current non-Gaussian tails, and pad the proposed future
        # support. Moments choose bounds only; probabilities follow the chain.
        tail=1e-35
        left=int(np.searchsorted(np.cumsum(w),tail))
        right=len(w)-1-int(np.searchsorted(np.cumsum(w[::-1]),tail))
        low=min(float(self.d[left]),float(np.min(sol.y[0]-self.window_z*np.sqrt(np.maximum(sol.y[1],0)))))-128
        high=max(float(self.d[right]),float(np.max(sol.y[0]+self.window_z*np.sqrt(np.maximum(sol.y[1],0)))))+256
        low=max(0,int(math.floor(low/128))*128);high=min(self.cap,max(low+512,int(math.ceil(high/128))*128))
        if high-low<512:low=max(0,high-512)
        newd=np.arange(low,high+1,dtype=float);newstate=np.zeros(len(newd))
        a=max(low,int(self.d[0]));b=min(high,int(self.d[-1]))
        if a<=b:newstate[a-low:b-low+1]=values[a-int(self.d[0]):b-int(self.d[0])+1]
        self.d=newd;self.state=newstate;self.window_max_size=max(self.window_max_size,len(newd))
    def propagate(self,start,end):
        if end<=start:return
        p=self.p;cap=self.cap;killed=self.killed
        chunk=60.
        if self.killed:
            bb,dd=rates(np.array([float(cap)]),self.qfun(start),start,self.alpha,self.theta,p)
            decay=max(0.,math.sqrt(float(bb[0]))-math.sqrt(float(dd[0])))**2
            chunk=min(chunk,5/max(decay,1e-12))
        if self.window:chunk=min(chunk,20.)
        proposed=sorted(set([end]+[t for t in [0.,(self.alpha-1)*self.theta]+list(np.arange(start+chunk,end,chunk)) if start+1e-8<t<end-1e-8]))
        boundaries=[]
        for stop in proposed:
            if boundaries and stop-boundaries[-1]<1e-8:boundaries[-1]=stop
            else:boundaries.append(stop)
        t=start
        for stop in boundaries:
            self._move_window(t,stop);d=self.d;size=len(d)
            if self.dynamic_tilt:
                def tilt_rhs(time,state):
                    fm,fd=fitness(time,self.alpha,self.theta,p)
                    birth_per_cell=p.b*fd/(1+(fm-1)*self.qfun(time))
                    u=float(state[0])
                    return [(1-u)*(birth_per_cell-p.b*u)]
                tilt_step=min(10.,max(.01,self.theta/3)) if 0<=t<(self.alpha-1)*self.theta+25*self.theta else 30.
                usol=solve_ivp(tilt_rhs,(t,stop),[self.u],method='DOP853',dense_output=True,
                    rtol=2e-12/self.refinement,atol=1e-14/self.refinement,max_step=tilt_step)
                if not usol.success:raise ArithmeticError(usol.message)
            def rhs(time,state):
                if self.dynamic_tilt:
                    u=float(usol.sol(time)[0])
                    upweight=1/u;downweight=u
                    eta_dot=-tilt_rhs(time,[u])[0]/u
                else:upweight=downweight=1.;eta_dot=0.
                birth,death=rates(d,self.qfun(time),time,self.alpha,self.theta,p)
                if not killed:
                    birth[-1]=0.;death[0]=0. # numerical reflection, checked against larger/full domains
                out=np.zeros_like(state)
                out[:size]=(-(birth+death)+eta_dot*d)*state[:size]
                out[1:size]+=upweight*birth[:-1]*state[:size-1]
                out[:size-1]+=downweight*death[1:]*state[1:size]
                if killed:out[-1]=upweight*birth[-1]*state[cap]+eta_dot*(cap+1)*state[-1]
                return out
            maxstep=min(10.,max(.01,self.theta/3)) if 0<=t<(self.alpha-1)*self.theta+25*self.theta else 30.
            sol=solve_ivp(rhs,(t,stop),self.state,method='LSODA',lband=1,uband=1,
                          rtol=2e-7/self.refinement,atol=1e-22/self.refinement,max_step=maxstep,t_eval=[stop])
            if not sol.success:raise ArithmeticError(sol.message)
            self.state=sol.y[:,-1];self.nfev+=sol.nfev
            if self.dynamic_tilt:self.u=float(usol.y[0,-1]);self.eta=-math.log(self.u)
            if killed:
                if self.state[-1]>0:
                    increment=self.logscale+math.log(self.state[-1])-self.eta*(cap+1)
                    self.absorbed_logmass=float(np.logaddexp(self.absorbed_logmass,increment))
                self.state[-1]=0.
            self._rescale();t=stop
            if not killed:
                lm=self.physical_logmass()
                probs=self.state/self.state.sum()
                if self.d[-1]==cap:self.max_upper_probability=max(self.max_upper_probability,float(probs[-1]))
                if self.window:
                    edge=max(float(probs[0]) if self.d[0]>0 else 0.,float(probs[-1]) if self.d[-1]<cap else 0.)
                    self.window_edge_probability=max(self.window_edge_probability,edge)
    def observation(self,vaf,time):
        marked=self.qfun(time)+(1-self.qfun(time))*self.d/self.p.N
        valid=(marked>0)&(marked<1)
        logs=np.full_like(marked,-math.inf)
        logs[valid]=-.5*((logit(2*vaf)-logit(marked[valid]))/self.p.sigma_vaf)**2-math.log(self.p.sigma_vaf*math.sqrt(2*math.pi))
        shift=float(np.max(logs))
        self.state[:len(self.d)]*=np.exp(logs-shift)
        if self.killed:
            self.state[-1]=0.;self.absorbed_logmass=-math.inf
        self.logscale+=shift;self._rescale()
    def predicted_vaf(self,time):
        values=self.state[:len(self.d)];positive=values>0
        logs=np.full(len(values),-math.inf)
        logs[positive]=np.log(values[positive])-self.eta*self.d[positive]
        weights=np.exp(logs-logsumexp(logs))
        marked=self.qfun(time)+(1-self.qfun(time))*self.d/self.p.N
        return float(np.sum(weights*marked)/2)
    def event_logdensity(self,time):
        birth,_=rates(self.d,self.qfun(time),time,self.alpha,self.theta,self.p)
        value=self.state[self.cap]
        if value<=0 or birth[-1]<=0:return -math.inf
        return self.logscale+math.log(value)+math.log(birth[-1])-self.eta*self.cap
    def log_cdf(self):
        if not self.killed:raise ValueError('CDF requires a killed chain')
        return self.absorbed_logmass


def suggested_cap(row,index,backgrounds,alpha,theta,p,end,mincap=512):
    """Conservative moment-based proposal, verified by larger-cap likelihoods.

    This proposal only sets a numerical reflecting boundary; it is not an event
    barrier or a distributional approximation used in the likelihood itself.
    """
    if mincap>=p.N:return p.N
    qfun=backgrounds.for_patient(index)
    def rhs(t,state):
        m=float(np.clip(state[0],0,p.N));v=max(0.,float(state[1]));q=qfun(t)
        fm,fd=fitness(t,alpha,theta,p)
        denominator=p.N*(1+(fm-1)*q)+m*((fd-1)-(fm-1)*q)
        numerator=p.mu*fm*q*p.N+m*(fd-p.mu*fm*q)
        pd=numerator/denominator
        derivative=((fd-p.mu*fm*q)*denominator-numerator*((fd-1)-(fm-1)*q))/denominator**2
        birth=p.b*(p.N-m)*pd;death=p.b*m*(1-pd)
        return [birth-death,2*p.b*(p.N*derivative-1)*v+birth+death]
    start=float(row['baseline_day'])
    if end<=start:return min(p.N,max(int(mincap),p.D_init+512))
    sol=solve_ivp(rhs,(start,end),[p.D_init,0.],rtol=2e-6,atol=1e-6,
                  max_step=min(10.,max(.1,theta)),dense_output=False)
    if not sol.success:raise ArithmeticError(sol.message)
    upper=float(np.max(sol.y[0]+12*np.sqrt(np.maximum(sol.y[1],0))))+256
    proposed=max(512,p.D_init+512,int(math.ceil(upper/256))*256)
    return min(p.N,max(int(mincap),int(math.ceil(proposed*(mincap/512)/256))*256))


def patient_loglik(row,index,backgrounds,alpha,theta,p,cap=512,refinement=1,force_unknown=False):
    outcome={'kind':'unknown'} if force_unknown else row['outcome']
    killed=outcome['kind']!='unknown'
    end=outcome.get('day',max([0.]+[r['day'] for r in row['observations']]))
    if killed and p.mu==0 and p.D_init==0 and outcome['kind']=='exact':return -math.inf,{}
    if not killed and not row['observations']:
        return 0.,{'nfev':0,'cap':0,'eta':0.,'max_upper_probability':0.,'negative_roundoff':0.,'VAF_predictions':[]}
    actual_cap=cap if killed else suggested_cap(row,index,backgrounds,alpha,theta,p,end,cap)
    duration=max(1.,end-row['baseline_day'])
    # Sharp observation likelihoods can otherwise amplify absolute-error noise
    # in the far D tail even without an observed clinical outcome.
    eta=tilt_for_horizon(alpha,theta,p,end,row['baseline_day']) if killed or p.sigma_vaf<.15 else 0.
    filt=ForwardFilter(p,alpha,theta,backgrounds.for_patient(index),killed=killed,cap=actual_cap,eta=eta,refinement=refinement,
                       window=not killed and cap<p.N,window_z=12*math.sqrt(cap/512))
    start=row['baseline_day']
    if start<0:
        filt.propagate(start,0.)
        if killed:
            # Condition on being initially event-free at CAR infusion.
            filt.condition_event_free()
    start=0.;predictions=[]
    for observation in row['observations']:
        time=observation['day']
        if killed and (time>end or outcome['kind']=='exact' and time>=end):raise ValueError('Post-outcome VAF in clinical process')
        filt.propagate(start,time)
        predictions.append({'day':time,'observed_vaf':observation['vaf'],'predicted_vaf':filt.predicted_vaf(time)})
        filt.observation(observation['vaf'],time);start=time
    filt.propagate(start,end)
    ll=filt.event_logdensity(end) if outcome['kind']=='exact' else filt.physical_logmass()
    return ll,{'nfev':filt.nfev,'cap':filt.cap,'eta':eta,'max_upper_probability':filt.max_upper_probability,
               'window_max_size':filt.window_max_size,'window_edge_probability':filt.window_edge_probability,
               'negative_roundoff':filt.negative_roundoff,'VAF_predictions':predictions}


def cohort_loglik(cohort,alpha,theta,p,cap=512,refinement=1,decompose=False):
    backgrounds=Backgrounds(cohort,alpha,theta,p,refinement=refinement)
    details=[];total=0.
    for i,row in enumerate(cohort):
        ll,diag=patient_loglik(row,i,backgrounds,alpha,theta,p,cap,refinement)
        if not math.isfinite(ll):return -math.inf,[]
        entry={'patient_id':row['patient_id'],'gene':row['gene'],'outcome':row['outcome']['kind'],'joint_loglik':ll,**diag}
        if decompose:
            marginal=ll if row['outcome']['kind']=='unknown' else patient_loglik(row,i,backgrounds,alpha,theta,p,cap,refinement,True)[0]
            entry['VAF_loglik']=marginal;entry['clinical_given_VAF_loglik']=ll-marginal
        details.append(entry);total+=ll
    return total,details


def risk_curve(alpha,theta,p,founder_fraction,baseline_day=0.,targets=(365,730,1095,1460),refinement=2):
    if p.D_init>=p.K:return {float(t):1. for t in targets}
    if p.mu==0 and p.D_init==0:return {float(t):0. for t in targets}
    row={'baseline_vaf':founder_fraction/2,'baseline_day':baseline_day,'observations':[],'outcome':{'kind':'unknown'}}
    bg=Backgrounds([row],alpha,theta,p,refinement)
    result={}
    for target in targets:
        eta=tilt_for_horizon(alpha,theta,p,target,baseline_day)
        f=ForwardFilter(p,alpha,theta,bg.for_patient(0),killed=True,cap=p.K-1,eta=eta,refinement=refinement)
        if baseline_day<0:
            f.propagate(baseline_day,0.);f.condition_event_free()
        f.propagate(0.,target)
        logcdf=f.log_cdf();total=f.physical_logmass(include_absorbed=True)
        if abs(total)>2e-5:raise ArithmeticError(f'Risk mass drift {total}')
        result[float(target)]=math.exp(logcdf-total) if logcdf-total>-745 else 0.
        if not -1e-10<=result[float(target)]<=1+1e-8:raise ArithmeticError('Risk outside [0,1]')
    return result
