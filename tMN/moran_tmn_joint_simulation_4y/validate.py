"""Independent numerical and stochastic checks for the confirmed fixed-area study."""
import json
import math
import subprocess
from dataclasses import replace
from pathlib import Path
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from scipy.special import logit,logsumexp
from model import Parameters,fitness,ForwardFilter,rates,risk_curve,cohort_loglik
from simulations import simulate_job

ROOT=Path(__file__).resolve().parent


def validate(config,cohort):
    p=Parameters(**config['reference']);checks={}
    assert fitness(-1,3,7,p)==(1.,1.)
    checks['pre_CAR_all_fitness_equal_one']=True
    integrals=[]
    for a in config['simulation']['alpha']:
        for theta in config['simulation']['theta_days']:
            value=quad(lambda t:fitness(t,a,theta,p)[0]-1,0,np.inf,epsabs=1e-10)[0]
            assert abs(value-p.area_days)<1e-8
            integrals.append({'alpha':a,'theta_days':theta,'integral_days':value,'peak_extra_fitness':fitness((a-1)*theta,a,theta,p)[0]-1})
    checks['fixed_area_integrals']=integrals
    zero=replace(p,mu=0)
    ll,_=cohort_loglik(cohort,3,7,zero)
    assert ll==-math.inf
    checks['observed_events_force_zero_likelihood_when_no_D_can_arise']=True
    # Independent sparse matrix exponential for a genuinely tiny event density
    # in the same finite-N, stochastic-D generator; no probability floor.
    pp=replace(p,mu=0,D_init=1)
    t=127.;eta=.95*math.log1p(1/(pp.b*t));q=.005
    f=ForwardFilter(pp,3,7,lambda time:q,killed=True,cap=pp.K-1,eta=eta,refinement=5)
    # area=0 gives a constant generator, making the matrix exponential independent.
    pp=replace(pp,area_days=0)
    f=ForwardFilter(pp,3,7,lambda time:q,killed=True,cap=pp.K-1,eta=eta,refinement=5)
    d=np.arange(pp.K,dtype=float);birth,death=rates(d,q,0,3,7,pp)
    diagonal=np.r_[-(birth+death),0.]
    lower=np.r_[np.exp(eta)*birth[:-1],np.exp(eta)*birth[-1]]
    upper=np.r_[np.exp(-eta)*death[1:],0.]
    generator=diags([lower,diagonal,upper],[-1,0,1],format='csr')
    initial=np.zeros(pp.K+1);initial[1]=1.
    exact=expm_multiply(generator*t,initial)
    exact_log_density=eta+math.log(exact[pp.K-1])+math.log(birth[-1])-eta*(pp.K-1)
    f.propagate(0,t);got=f.event_logdensity(t)
    error=abs(got-exact_log_density)
    assert error<.02,(got,exact_log_density,error)
    checks['tiny_event_density_independent_matrix_check']={'day':t,'matrix_log_density':exact_log_density,'filter_log_density':got,'absolute_log_error':error}
    seeded=[]
    for initial_D in (10,100):
        pp_seeded=replace(pp,D_init=initial_D)
        initial=np.zeros(pp.K+1);initial[initial_D]=1.
        exact_seeded=expm_multiply(generator*t,initial)
        exact_ll=eta*initial_D+math.log(exact_seeded[pp.K-1])+math.log(birth[-1])-eta*(pp.K-1)
        filt=ForwardFilter(pp_seeded,3,7,lambda time:q,killed=True,cap=pp.K-1,eta=eta,refinement=5)
        filt.propagate(0,t);value=filt.event_logdensity(t)
        assert abs(value-exact_ll)<.02,(initial_D,value,exact_ll)
        seeded.append({'D_init':initial_D,'matrix_log_density':exact_ll,'filter_log_density':value,'absolute_log_error':abs(value-exact_ll)})
    checks['preexisting_D_independent_density_checks']=seeded
    # Independent rare VAF-observation likelihood on an unrestricted chain.
    # This catches absolute tail noise amplified by a narrow measurement model.
    obs_p=replace(p,area_days=0,D_init=0,sigma_vaf=.1)
    obs_cap=2048;obs_q=.04;obs_time=185.;constant_eta=.12
    obs_d=np.arange(obs_cap+1,dtype=float)
    ub,down=rates(obs_d,obs_q,0,3,7,obs_p);ub[-1]=0.
    A=diags([np.exp(constant_eta)*ub[:-1],-(ub+down),np.exp(-constant_eta)*down[1:]],[-1,0,1],format='csr')
    vec=np.zeros(obs_cap+1);vec[0]=1.;exact_vec=expm_multiply(A*obs_time,vec)
    marked=obs_q+(1-obs_q)*obs_d/obs_p.N
    logweights=-.5*((logit(.38)-logit(marked))/obs_p.sigma_vaf)**2-math.log(obs_p.sigma_vaf*math.sqrt(2*math.pi))
    positive=exact_vec>0
    exact_vaf_ll=float(logsumexp(np.log(exact_vec[positive])-constant_eta*obs_d[positive]+logweights[positive]))
    obs_filter=ForwardFilter(obs_p,3,7,lambda time:obs_q,killed=False,cap=obs_cap,eta=.1,refinement=5)
    obs_filter.propagate(0,obs_time)
    positive=obs_filter.state>0
    actual_vaf_ll=obs_filter.logscale+float(logsumexp(np.log(obs_filter.state[positive])-obs_filter.eta*obs_d[positive]+logweights[positive]))
    assert abs(actual_vaf_ll-exact_vaf_ll)<.02,(actual_vaf_ll,exact_vaf_ll)
    checks['narrow_VAF_likelihood_independent_matrix_check']={'sigma_vaf':.1,'matrix_loglik':exact_vaf_ll,'filter_loglik':actual_vaf_ll,'absolute_log_error':abs(actual_vaf_ll-exact_vaf_ll)}
    # Small-N full stochastic W/C/D forward equation vs q-background approximation.
    small=replace(p,N=80,mu=.002,threshold_fraction=.05)
    states=[(c,d) for d in range(small.K) for c in range(small.N-d+1)]
    lookup={state:i for i,state in enumerate(states)};size=len(states)
    cc=np.array([x[0] for x in states]);dd=np.array([x[1] for x in states]);ww=small.N-cc-dd
    destinations=[(1,0),(0,1),(-1,0),(-1,1),(0,-1),(1,-1)]
    dest_indices=[np.array([lookup.get((int(c+dc),int(d+ddelta)),size) for c,d in states]) for dc,ddelta in destinations]
    def rhs(time,prob):
        fc,fd=fitness(time,3,5,small);F=ww+fc*cc+fd*dd
        pw=ww/F;pc=(1-small.mu)*fc*cc/F;pd=(fd*dd+small.mu*fc*cc)/F
        transitions=[small.b*ww*pc,small.b*ww*pd,small.b*cc*pw,small.b*cc*pd,small.b*dd*pw,small.b*dd*pc]
        out=np.zeros(size+1)
        for rate,dest in zip(transitions,dest_indices):
            flow=rate*prob[:size];out[:size]-=flow;np.add.at(out,dest,flow)
        return out
    initial=np.zeros(size+1);initial[lookup[(16,0)]]=1
    exact=solve_ivp(rhs,(0,30),initial,t_eval=[10,20,30],rtol=1e-8,atol=1e-12,max_step=.25)
    assert exact.success
    hybrid=risk_curve(3,5,small,.2,targets=[10,20,30],refinement=5)
    difference=max(abs(exact.y[-1,j]-hybrid[t]) for j,t in enumerate([10,20,30]))
    assert difference<.005,difference
    checks['small_N_full_Moran_comparison']={'N':80,'exact':exact.y[-1].tolist(),'hybrid':hybrid,'max_absolute_difference':difference}
    # Native exact simulator: neutrality, population conservation and all-zero D.
    simulation=dict(config['simulation'],horizon_days=30,replicates=200,initial_C_fraction=.2)
    neutral=replace(p,N=200,area_days=0,s_D=0,mu=0)
    meta=simulate_job(('validation_neutral',neutral.__dict__,3.,7.,0.,simulation))
    data=np.load(ROOT/'results/trajectories'/f"{meta['name']}.npz")['counts']
    assert not data[:,:,1].any()
    expected=40.;variance=40*(200-40)*(1-math.exp(-2*neutral.b*30/200))
    mean=float(data[:,-1,0].mean());se=math.sqrt(variance/200)
    assert abs(mean-expected)<4*se,(mean,expected,se)
    checks['exact_Moran_neutral_mean']={'expected_C':expected,'observed_C_mean':mean,'Monte_Carlo_SE':se,'zero_D_invariant':True}
    # Check source outcomes enter exactly once and VAF-only unknown is not censoring.
    checks['cohort']={'patients':len(cohort),'events':sum(r['outcome']['kind']=='exact' for r in cohort),'censor':sum(r['outcome']['kind']=='censor' for r in cohort),'unknown':sum(r['outcome']['kind']=='unknown' for r in cohort)}
    assert checks['cohort']=={'patients':14,'events':2,'censor':1,'unknown':11}
    checks['status']='passed'
    checks['scope']='Initial invariant, independent tiny-flux, small-N and neutral-simulation checks; not clinical validation.'
    (ROOT/'results/validation.json').write_text(json.dumps(checks,indent=2))
    print('VALIDATION PASSED',flush=True)
    return checks
