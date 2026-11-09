"""Gamma Moran mechanisms and a shared-state VAF/event forward likelihood.

The production-sized filter is HYBRID: M/(N-D) is deterministic, D is a
stochastic birth/death chain. No claim of an exact two-dimensional Moran fit.
Backward Euler keeps the forward probabilities positive. Refine the timestep
before using likelihoods or risks for substantive inference.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded
from scipy.special import expit, logit, gammaln


@dataclass(frozen=True)
class Parameters:
    N: int = 100000
    b: float = 1/7
    amplitude: float = 0.1
    gamma_normalization: str = 'peak'
    mu: float = 1e-6
    r_D: float = 0.003
    D_init: int = 0
    threshold_fraction: float = 0.01
    sigma_vaf: float = 0.25
    diagnosis_delay_days: float = 0

    @property
    def K(self):
        return math.floor(self.threshold_fraction*self.N)+1

    def check(self):
        assert self.N >= 2 and self.b > 0 and self.amplitude >= 0
        assert 0 <= self.mu <= 1 and self.r_D >= 0
        assert 0 < self.threshold_fraction < 1 and self.sigma_vaf > 0
        assert 0 <= self.D_init < self.K
        assert self.gamma_normalization in ('peak', 'area')
        if self.diagnosis_delay_days != 0:
            raise NotImplementedError('Nonzero diagnosis delay needs a separate observation model.')


def fitness(t, alpha, theta, p):
    if alpha <= 1 or theta <= 0:
        raise ValueError('Require gamma shape alpha>1 and scale theta>0.')
    times = np.asarray(t, dtype=float)
    h = np.zeros_like(times)
    mask = times > 0
    tt = times[mask]
    if p.gamma_normalization == 'peak':
        mode = (alpha-1)*theta
        h[mask] = np.exp((alpha-1)*np.log(tt/mode)-(tt-mode)/theta)
    else:
        h[mask] = np.exp((alpha-1)*np.log(tt)-tt/theta-gammaln(alpha)-alpha*np.log(theta))
    fm = 1+p.amplitude*h
    fd = fm + np.where(times >= 0, p.r_D/p.b, 0.)
    if times.ndim == 0:
        return float(fm), float(fd)
    return fm, fd


def background(t0, end, founder_fraction, alpha, theta, p):
    if not 0 <= p.D_init/p.N <= founder_fraction < 1:
        raise ValueError('Initial D exceeds founder burden or baseline mapping is invalid.')
    q0 = (founder_fraction-p.D_init/p.N)/(1-p.D_init/p.N)
    def rhs(t, state):
        q = state[0]
        fm, _ = fitness(t, alpha, theta, p)
        return [p.b*q*(1-q)*((1-p.mu)*fm-1)/(1+q*(fm-1))]
    # Splitting stops adaptive integration from stepping across a narrow pulse.
    boundaries = sorted(set([t0, end]+[t for t in (0., (alpha-1)*theta, (alpha-1)*theta+20*theta) if t0 < t < end]))
    segments = []
    q = q0
    for start, stop in zip(boundaries[:-1], boundaries[1:]):
        step = min(10., max(0.02, theta/3)) if 0 <= start < (alpha-1)*theta+20*theta else 20.
        sol = solve_ivp(rhs, (start,stop), [q], method='DOP853', dense_output=True,
                        max_step=step, rtol=2e-10, atol=1e-12)
        if not sol.success:
            raise ArithmeticError(sol.message)
        segments.append((start,stop,sol.sol))
        q = float(sol.y[0,-1])
    def evaluate(t):
        if t == t0:
            return q0
        for start,stop,sol in segments:
            if start-1e-10 <= t <= stop+1e-10:
                return float(sol(t)[0])
        raise ValueError('Time outside deterministic background interval.')
    return evaluate


def rates(d, q, t, alpha, theta, p):
    fm, fd = fitness(t, alpha, theta, p)
    remaining = p.N-d
    m = remaining*q
    total_fitness = remaining+(fm-1)*m+fd*d
    pd = (fd*d+p.mu*fm*m)/total_fitness
    return p.b*remaining*pd, p.b*d*(1-pd)


def propagate(prob, start, stop, qfun, alpha, theta, p, *, killed, max_step=1.):
    """Implicit tridiagonal forward evolution; return distribution and lost flux."""
    prob = np.asarray(prob, dtype=float).copy()
    d = np.arange(len(prob), dtype=float)
    if not killed and len(prob) != p.N+1:
        raise ValueError('Unknown-outcome physical process must span D=0..N.')
    if killed and len(prob) != p.K:
        raise ValueError('Killed process must span D=0..K-1.')
    if stop < start:
        raise ValueError('Backward propagation is not supported.')
    loss = 0.
    initial_mass = prob.sum()
    band = np.zeros((3,len(prob)))
    boundaries = sorted(set([stop]+[t for t in (0., (alpha-1)*theta) if start < t < stop]))
    t = start
    for end in boundaries:
        while t < end-1e-12:
            # Resolve the gamma pulse as well as the global BE timestep.
            h = min(max_step, end-t)
            if 0 <= t <= (alpha-1)*theta+15*theta:
                h = min(h, theta/8)
            middle = t+h/2
            birth, death = rates(d, qfun(middle), middle, alpha, theta, p)
            band[1] = 1+h*(birth+death)
            band[0,1:] = -h*death[1:]
            band[2,:-1] = -h*birth[:-1]
            prob = solve_banded((1,1), band, prob, check_finite=False)
            if killed:
                loss += h*birth[-1]*prob[-1]
            t += h
    if prob.min() < -1e-12 or abs(prob.sum()+loss-initial_mass) > 1e-8:
        raise ArithmeticError('Forward probability positivity/mass check failed.')
    return prob, float(loss)


def initial_at_zero(t0, founder_fraction, qfun, alpha, theta, p, max_step):
    if t0 > 0:
        raise ValueError('Post-CAR baseline needs an explicit left-entry model.')
    prob = np.zeros(p.K)
    prob[p.D_init] = 1.
    if t0 < 0:
        prob, _ = propagate(prob,t0,0.,qfun,alpha,theta,p,killed=True,max_step=max_step)
        mass = prob.sum()
        if mass <= 0:
            raise ArithmeticError('No event-free mass at infusion.')
        prob /= mass
    return prob


def observation_log_weights(vaf, q, d, p):
    if not 0 < vaf < 0.5:
        raise ValueError('Observation outside positive diploid-heterozygous mapping.')
    marked = q+(1-q)*d/p.N
    mu_z = np.full_like(marked, np.nan)
    inside = (marked > 0)&(marked < 1)
    mu_z[inside] = logit(marked[inside])
    logw = np.full_like(marked, -np.inf)
    logw[inside] = -.5*((logit(2*vaf)-mu_z[inside])/p.sigma_vaf)**2-math.log(p.sigma_vaf*math.sqrt(2*math.pi))
    return logw


def joint_loglik(*, alpha, theta, p, baseline_day, baseline_vaf, observations, outcome, max_step=1.):
    """One patient's same-state approximate joint likelihood, conditional on baseline.

    outcome: {'kind': 'unknown'}, {'kind':'censor','day':C},
    {'kind':'exact','day':T}, or {'kind':'interval','left':L,'right':R}.
    Interval VAF must precede L. Only positive prediagnosis observations supported.
    """
    p.check()
    kind = outcome['kind']
    if kind not in ('unknown','censor','exact','interval'):
        raise ValueError('Unknown outcome encoding.')
    obs = sorted(observations)
    last = max([0.]+[t for t,v in obs])
    end = last if kind == 'unknown' else float(outcome['right'] if kind == 'interval' else outcome['day'])
    if last > end or any(t < 0 for t,v in obs):
        raise ValueError('Only post-infusion, pre-outcome observations supported.')
    if kind == 'exact' and any(t >= end for t,v in obs):
        raise ValueError('Exact-event VAF must be strictly prediagnosis.')
    if kind == 'interval' and last > outcome['left']:
        raise ValueError('Interval-event filter currently requires VAF at or before interval start.')
    qfun = background(baseline_day,max(end,1e-8),2*baseline_vaf,alpha,theta,p)
    prob = initial_at_zero(baseline_day,2*baseline_vaf,qfun,alpha,theta,p,max_step)
    killed = kind != 'unknown'
    if not killed:
        prob = np.pad(prob,(0,p.N+1-p.K))
    log_likelihood = 0.
    t = 0.
    d = np.arange(len(prob),dtype=float)
    for day,vaf in obs:
        prob,_ = propagate(prob,t,day,qfun,alpha,theta,p,killed=killed,max_step=max_step)
        logw = observation_log_weights(vaf,qfun(day),d,p)
        finite = logw[np.isfinite(logw)]
        if not len(finite):
            return -math.inf
        shift = finite.max()
        prob *= np.exp(logw-shift)
        z = prob.sum()
        if z <= 0:
            return -math.inf
        log_likelihood += math.log(z)+shift
        prob /= z
        t = day
    if kind == 'interval':
        prob,_ = propagate(prob,t,outcome['left'],qfun,alpha,theta,p,killed=True,max_step=max_step)
        prob,event_probability = propagate(prob,outcome['left'],end,qfun,alpha,theta,p,killed=True,max_step=max_step)
        contribution = event_probability
    else:
        prob,_ = propagate(prob,t,end,qfun,alpha,theta,p,killed=killed,max_step=max_step)
        if kind == 'exact':
            birth,_ = rates(d,qfun(end),end,alpha,theta,p)
            contribution = birth[-1]*prob[-1]
        else:
            contribution = prob.sum()
    return log_likelihood+math.log(contribution) if contribution > 0 else -math.inf


def first_passage(*, alpha, theta, p, baseline_vaf, baseline_day=0., targets=(365,730,1095,1460), max_step=1.):
    p.check()
    targets = sorted(set(float(t) for t in targets))
    if p.mu == 0 and p.D_init == 0:
        return {t:0. for t in targets}
    qfun = background(baseline_day,max(targets),2*baseline_vaf,alpha,theta,p)
    prob = initial_at_zero(baseline_day,2*baseline_vaf,qfun,alpha,theta,p,max_step)
    answer = {}
    lost,t = 0.,0.
    for end in targets:
        prob,loss = propagate(prob,t,end,qfun,alpha,theta,p,killed=True,max_step=max_step)
        lost += loss
        answer[end] = lost
        t = end
    return answer


def first_passage_adaptive(*, alpha, theta, p, baseline_vaf, baseline_day=0., targets=(365,730,1095,1460), refinement=1.):
    """Banded LSODA with explicit absorbing flux; resolves small early risks."""
    p.check()
    targets=sorted(set(float(t) for t in targets))
    if p.mu==0 and p.D_init==0:
        return {t:0. for t in targets}
    qfun=background(baseline_day,max(targets),2*baseline_vaf,alpha,theta,p)
    initial=initial_at_zero(baseline_day,2*baseline_vaf,qfun,alpha,theta,p,.25)
    prob=np.r_[initial,0.]
    d=np.arange(p.K,dtype=float)
    t=0.; answer={}
    boundaries=sorted(set(targets+[v for v in ((alpha-1)*theta,(alpha-1)*theta+20*theta) if 0<v<targets[-1]]))
    for end in boundaries:
        def rhs(time,state):
            birth,death=rates(d,qfun(time),time,alpha,theta,p)
            out=np.zeros(p.K+1)
            out[:p.K]=-(birth+death)*state[:p.K]
            out[1:p.K]+=birth[:-1]*state[:p.K-1]
            out[:p.K-1]+=death[1:]*state[1:p.K]
            out[-1]=birth[-1]*state[p.K-1]
            return out
        step=min(10.,max(.02,theta/4)) if t<(alpha-1)*theta+20*theta else 30.
        sol=solve_ivp(rhs,(t,end),prob,method='LSODA',lband=1,uband=1,
            rtol=2e-8/refinement,atol=1e-13/refinement,max_step=step)
        if not sol.success:
            raise ArithmeticError(sol.message)
        if abs(sol.y.sum(axis=0)-1).max()>1e-8 or sol.y.min() < -1e-9:
            raise ArithmeticError('Adaptive forward-chain mass or positivity error.')
        prob=sol.y[:,-1]; t=end
        if end in targets:
            if not 0<=prob[-1]<=1:
                raise ArithmeticError('Adaptive absorption probability outside [0,1].')
            answer[end]=float(prob[-1])
    return answer


def meanfield_predictions(alpha, theta, p, baselines, requested):
    """Coupled mean-field VAF diagnostic only; not the forward joint filter."""
    count = len(baselines)
    init = np.zeros((2,count))
    for i,(day,vaf) in enumerate(baselines):
        founder = 2*vaf
        y = p.D_init/p.N
        if day > 0 or not 0 <= y <= founder < 1:
            raise ValueError('Incompatible baseline or D_init.')
        x = (founder-y)*math.exp(p.b*p.mu*day)
        init[:,i] = x,founder-x
    def rhs(t,state):
        x,y = state.reshape(2,count)
        fm,fd = fitness(t,alpha,theta,p)
        fbar = 1-x-y+fm*x+fd*y
        return np.concatenate([p.b*((1-p.mu)*fm*x/fbar-x), p.b*((fd*y+p.mu*fm*x)/fbar-y)])
    times = sorted(set(float(t) for i,t in requested))
    if any(t < 0 for t in times):
        raise ValueError('Diagnostic currently requires post-CAR follow-up.')
    end = max(times)
    boundaries = sorted(set([0.,end]+[t for t in ((alpha-1)*theta,(alpha-1)*theta+20*theta) if 0 < t < end]))
    segments=[]
    state=init.ravel()
    for start,stop in zip(boundaries[:-1],boundaries[1:]):
        step=min(20.,max(.05,theta/2)) if start<(alpha-1)*theta+20*theta else 30.
        sol=solve_ivp(rhs,(start,stop),state,method='DOP853',dense_output=True,max_step=step,rtol=2e-8,atol=1e-11)
        if not sol.success:
            raise ArithmeticError(sol.message)
        state=sol.y[:,-1]; segments.append((start,stop,sol.sol))
    out=[]
    for i,t in requested:
        state=init.ravel() if t==0 else next(sol(t) for start,stop,sol in segments if start-1e-10<=t<=stop+1e-10)
        marked=state[i]+state[count+i]
        if not -1e-8<=marked<=1+1e-8:
            raise ArithmeticError('Mean-field marked fraction outside [0,1].')
        out.append(marked/2)
    return np.array(out)
