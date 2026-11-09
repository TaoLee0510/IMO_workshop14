"""English figures at frozen fitted parameters; no outcome-conditioned risk forecasts."""
import math,json,csv,argparse
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from scipy.stats import norm,gamma
from scipy.integrate import quad
from model import fitness
BLUE='#1765b0';RED='#c82d36';GRAY='#687481';GOLD='#c58a2f'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
def read_csv(name):
    rows=list(csv.DictReader((ROOT/'results'/name).open()))
    for r in rows:
        for k,v in list(r.items()):
            if v=='':r[k]=None
            elif k not in ['patient_id','gene','sample_type','outcome','kind','conditioning']:
                try:r[k]=float(v)
                except ValueError:pass
    return rows
def save(fig,name):
    path=ROOT/'figures'/name;path.parent.mkdir(exist_ok=True,parents=True)
    fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight');fig.savefig(path.with_suffix('.png'),dpi=180,bbox_inches='tight');plt.close(fig)
def event_line(ax,row,years=False):
    out=row['outcome']
    if out['kind']!='unknown':
        x=out['day']/(365 if years else 1)
        ax.axvline(x,color=RED if out['kind']=='exact' else GRAY,ls=':',lw=1.1,label='Observed tMN' if out['kind']=='exact' else 'Right censor')
def vaf_panel(ax,row,rec,early=False,legend=False):
    c=rec['curves'];x=np.array([p['day'] for p in c]);mean=np.array([p['mean_VAF'] for p in c])*100
    ax.fill_between(x,[p['lower_VAF']*100 for p in c],[p['upper_VAF']*100 for p in c],color=BLUE,alpha=.14,label='95% observation prediction interval')
    ax.plot(x,mean,color=BLUE,lw=1.8,label='Sequential predictive mean')
    ax.scatter([row['baseline_day']],[row['baseline_vaf']*100],s=40,marker='D',facecolors='white',edgecolors='black',zorder=4,label='Baseline')
    if rec['predictions']:
        ax.scatter([p['day'] for p in rec['predictions']],[p['observed_VAF']*100 for p in rec['predictions']],color='black',s=22,zorder=4,label='Observed follow-up')
        ax.scatter([p['day'] for p in rec['predictions']],[p['mean_VAF']*100 for p in rec['predictions']],facecolors='white',edgecolors=BLUE,s=20,zorder=4,label='One-step predictive mean')
    else:ax.text(.03,.9,'Baseline only; no follow-up VAF',transform=ax.transAxes,fontsize=8)
    ax.axvline(0,color=GRAY,lw=.8,alpha=.6);event_line(ax,row)
    end=min(1460,row['outcome'].get('day',1460)) if row['outcome']['kind']=='exact' else 1460
    ax.set(xlim=(min(-5,row['baseline_day']-2),min(90,end) if early else end),xlabel='Days after CAR infusion',ylabel='VAF (%)')
    ax.set_ylim(bottom=0,top=min(50,max(max(p['upper_VAF'] for p in c)*100,row['baseline_vaf']*100,*[p['observed_VAF']*100 for p in rec['predictions']])*1.12))
    ax.grid(axis='y',alpha=.15)
    if legend:ax.legend(fontsize=7,loc='best')
def trajectory_panel(ax,row,counts,index,early=False,legend=False):
    days=np.arange(counts.shape[1]);sl=days<=90 if early else days<=1460
    color=BLUE if index==1 else RED;values=counts[:,:,index]/P.N*100
    ax.plot(days[sl]/365,values[:,sl].T,color=color,alpha=.035,lw=.45,rasterized=True)
    ax.plot(days[sl]/365,values[:,sl].mean(axis=0),color=color,lw=2,label='Mean of all 200')
    ax.set(xlim=(0,(90 if early else 1460)/365),xlabel='Years after CAR infusion',ylabel=f"{'Single mutant C' if index==1 else 'Double hit D'} (%)",ylim=(0,None))
    ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,4));event_line(ax,row,years=True);ax.grid(axis='y',alpha=.15)
    if np.max(values)==0:ax.text(.5,.6,'All 200 trajectories remain at zero',ha='center',transform=ax.transAxes,fontsize=8)
    if legend:ax.legend(handles=[Line2D([0],[0],color=color,alpha=.25,lw=.7,label='Individual path'),Line2D([0],[0],color=color,lw=2,label='Mean of all 200')],fontsize=8)
def risk_panel(ax,row,points,log=False,legend=False):
    x=[p['day']/365 for p in points];y=[p['probability'] for p in points]
    if log:
        ax.plot(x,[p['log10_probability'] for p in points],color=RED,lw=1.7)
        ax.set_ylabel('log10 cumulative tMN probability')
    else:
        ax.plot([0]+x,[0]+y,color=RED,lw=1.7,label='Baseline-initialized model risk')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,3));ax.set_ylabel('Cumulative tMN probability')
        ax.set_ylim(bottom=0)
    p365=next(p for p in points if p['day']==365)
    ax.scatter([1],[p365['log10_probability'] if log else p365['probability']],color=RED,s=25,zorder=3)
    event_line(ax,row,years=True);ax.axvline(1,color=GRAY,lw=.7,alpha=.5)
    ax.set(xlim=(0,4),xlabel='Years after CAR infusion');ax.grid(axis='y',alpha=.15)
    if legend:ax.legend(fontsize=7)

def gamma_figure():
    t=np.linspace(0,60,1600);g=gamma.pdf(t,A,scale=THETA)
    fig,axes=plt.subplots(1,2,figsize=(11,3.6),layout='constrained')
    axes[0].fill_between(t,g,color=BLUE,alpha=.18);axes[0].plot(t,g,color=BLUE,lw=2)
    axes[0].axvline((A-1)*THETA,color=GRAY,ls=':',label=f'Mode {(A-1)*THETA:.3f} days')
    axes[0].set(xlabel='Days after CAR infusion',ylabel='Gamma density (per day)',title=f'alpha = {A:.5f}, theta = {THETA:g} day')
    axes[0].text(.97,.88,'Integral on [0, infinity) = 1',ha='right',transform=axes[0].transAxes,fontsize=9);axes[0].legend(fontsize=8)
    tt=np.linspace(-10,60,1800);fc=np.array([fitness(t,A,THETA,P)[0] for t in tt]);fd=np.array([fitness(t,A,THETA,P)[1] for t in tt])
    axes[1].plot(tt,np.ones_like(tt),color=GRAY,label='Wild type W');axes[1].plot(tt,fc,color=BLUE,label='Single mutant C');axes[1].plot(tt,fd,color=RED,label='Double hit D')
    axes[1].axvline(0,color=GRAY,ls=':',lw=.8);axes[1].set(xlabel='Days after CAR infusion',ylabel='Relative fitness',title='Fitness functions used in this fit');axes[1].legend(fontsize=8)
    save(fig,'gamma_and_fitness')
    integral=quad(lambda u:math.exp((A-1)*math.log(u)-u-math.lgamma(A)) if u>0 else 0,0,np.inf)[0]
    assert abs(integral-1)<1e-9
    write_json(ROOT/'results/gamma_check.json',{'alpha':A,'theta_days':THETA,'integral_0_infinity':integral,'peak_day':(A-1)*THETA,'peak_density_per_day':float(gamma.pdf((A-1)*THETA,A,scale=THETA))})

def diagnostics():
    rows=read_csv('VAF_diagnostics.csv');patients={r['patient_id']:i for i,r in enumerate(COHORT)}
    palette=plt.cm.tab20(np.linspace(0,1,14));colors=[palette[patients[r['patient_id']]] for r in rows]
    fig,ax=plt.subplots(2,2,figsize=(11,9),layout='constrained')
    obs=np.array([r['observed_VAF'] for r in rows])*100;pred=np.array([r['mean_VAF'] for r in rows])*100
    ax[0,0].scatter(obs,pred,color=BLUE,s=25);ax[0,0].plot([0,50],[0,50],ls='--',color=GRAY)
    ax[0,0].set(xlabel='Observed VAF (%)',ylabel='One-step predicted VAF (%)',title='Observed and predicted VAF',xlim=(0,50),ylim=(0,50))
    residual=np.array([r['predictive_residual'] for r in rows]);theory=norm.ppf((np.arange(len(rows))+.5)/len(rows))
    ax[0,1].scatter(theory,np.sort(residual),color=BLUE,s=22);bounds=[theory.min(),theory.max()]
    ax[0,1].plot(bounds,bounds,ls='--',color=GRAY);ax[0,1].set(xlim=(theory.min()-.2,theory.max()+.2),xlabel='Standard normal quantile',ylabel='Predictive quantile residual',title='Predictive residual Q-Q plot')
    ax[1,0].scatter([r['day'] for r in rows],residual,color=BLUE,s=25);ax[1,0].axhline(0,color=GRAY,lw=1);ax[1,0].axhspan(-1.96,1.96,color=GRAY,alpha=.08)
    ax[1,0].set(xlabel='Days after CAR infusion',ylabel='Predictive quantile residual',title='Residuals over time')
    detail=FIT['patient_details'];y=np.arange(14);v=np.array([-r['VAF_loglik'] for r in detail]);c=np.array([-r['clinical_given_VAF_loglik'] for r in detail])
    ax[1,1].barh(y,v,color=BLUE,label='VAF NLL');ax[1,1].barh(y,c,left=v,color=RED,label='Clinical given VAF NLL')
    ax[1,1].set_yticks(y,[r['patient_id'].replace('Patient_','P') for r in detail],fontsize=8);ax[1,1].invert_yaxis();ax[1,1].set(xlabel='Negative log likelihood',title='Patient contributions');ax[1,1].legend(fontsize=7)
    save(fig,'fit_diagnostics')
    fig,ax=plt.subplots(figsize=(11,4),layout='constrained')
    values=[[r['predictive_residual'] for r in rows if r['patient_id']==p['patient_id']] for p in COHORT]
    for i,val in enumerate(values):
        if val:ax.scatter(np.full(len(val),i),val,color=palette[i],s=25)
    ax.set_xticks(np.arange(14),[r['patient_id'].replace('Patient_','P') for r in COHORT],rotation=45)
    ax.axhline(0,color=GRAY,lw=1);ax.axhspan(-1.96,1.96,color=GRAY,alpha=.08);ax.set(ylabel='Predictive quantile residual',title='VAF residuals by patient; baseline-only patients have no residuals')
    save(fig,'residuals_by_patient')
    write_json(ROOT/'results/diagnostic_summary.json',{'followup_observations':len(rows),'VAF_RMSE_percentage_points':float(np.sqrt(np.mean((obs-pred)**2))),
      'predictive_interval_coverage_fraction':sum(r['lower_VAF']<=r['observed_VAF']<=r['upper_VAF'] for r in rows)/len(rows),
      'PIT_rounded_to_endpoint_count':sum(r['PIT'] in [0,1] for r in rows),
      'quantile_residual_note':'Normal scores use stable log-CDF/log-survival inversion with no tail probability clipping; rounded PIT and stable log tails are retained'})

def likelihood_figures():
    rows=read_csv('likelihood_surface.csv');aa=sorted(set(r['alpha'] for r in rows));tt=sorted(set(r['theta_days'] for r in rows))
    z=np.array([[next(r['delta_nll'] for r in rows if r['alpha']==a and r['theta_days']==t) for a in aa] for t in tt])
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    mesh=ax.pcolormesh(aa,tt,z,shading='auto',cmap='viridis',vmin=0,vmax=min(6,z.max()));fig.colorbar(mesh,ax=ax,label='Delta joint NLL from frozen fit')
    ax.contour(aa,tt,z,levels=[.25,.5,1,2,4],colors='white',linewidths=.6)
    ax.scatter([A],[THETA],marker='*',s=100,color=RED,edgecolor='white',zorder=4,clip_on=False)
    ax.set(yscale='log',xlabel='Gamma shape alpha',ylabel='Gamma scale theta (days)',title='Joint likelihood surface; theta lower bound = 1 day');save(fig,'likelihood_surface')
    rows=read_csv('likelihood_profiles.csv');fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
    for ax,kind,label in zip(axes,['alpha_profile','theta_profile'],['Alpha (shape)','Theta (days)']):
        points=sorted([r for r in rows if r['kind']==kind],key=lambda r:r['fixed_value'])
        ax.plot([r['fixed_value'] for r in points],[r['delta_nll'] for r in points],marker='o',color=BLUE,ms=4)
        if kind=='theta_profile':ax.set_xscale('log')
        ax.set(xlabel=label,ylabel='Delta joint NLL',title='Nuisance parameter re-optimized');ax.axhline(0,color=GRAY,lw=.7)
    save(fig,'likelihood_profiles')
    allpoints=rows+read_csv('likelihood_surface.csv');best=min(allpoints,key=lambda r:r['joint_nll'])
    write_json(ROOT/'results/likelihood_diagnostic_summary.json',{'frozen_joint_nll':FIT['joint_nll'],'best_diagnostic_point':best,'better_than_frozen_by_more_than_0.01':best['joint_nll']<FIT['joint_nll']-.01,'profile_interpretation':'Conditional searches at boundaries; no chi-square confidence interval asserted'})

def heatmap():
    rows=sorted(read_csv('patient_risks.csv'),key=lambda r:r['baseline_founder_fraction'])
    values=np.array([r['year1_probability'] for r in rows]);assert np.all(values>0)
    fig,ax=plt.subplots(figsize=(9,8.3),layout='constrained')
    image=ax.imshow(values[:,None],aspect='auto',cmap='YlOrRd',norm=LogNorm(vmin=values.min(),vmax=values.max()))
    labels=[]
    for r in rows:
        status=f"tMN day {r['outcome_day']:.0f}" if r['outcome']=='exact' else ('censored day 1460' if r['outcome']=='censor' else 'outcome unknown')
        labels.append(r['patient_id'].replace('Patient_','P')+'  '+r['gene']+f"  VAF={100*r['baseline_VAF']:.3f}%\nC0={100*r['baseline_founder_fraction']:.3f}%, baseline day {r['baseline_day']:.0f}; "+status)
    ax.set_yticks(np.arange(14),labels,fontsize=8);ax.set_xticks([0],['tMN probability by day 365'])
    for i,v in enumerate(values):
        tone=image.norm(v);ax.text(0,i,f'{v:.3e}',ha='center',va='center',color='white' if tone>.6 else '#202020',fontsize=10)
    ax.set_title('One-year patient-specific tMN model probability\nOwn baseline state; frozen fitted parameters',fontsize=11,pad=12)
    fig.colorbar(image,ax=ax,label='Probability (logarithmic color scale)',shrink=.8,pad=.03);save(fig,'one_year_tMN_heatmap')

def patient_figures():
    vaf_all={};risk_all={}
    for row in COHORT:
        pid=row['patient_id'];vaf=json.loads((ROOT/'results/patients'/f'{pid}_VAF.json').read_text());risk=json.loads((ROOT/'results/patients'/f'{pid}_risks.json').read_text())['points']
        counts=np.load(ROOT/'results/trajectories'/f'{pid}.npz')['counts'];vaf_all[pid]=vaf;risk_all[pid]=risk
        title=f"{pid.replace('_',' ')}   {row['gene']}   {row['sample_type']}"
        fig,axes=plt.subplots(1,2,figsize=(11,3.6),layout='constrained');vaf_panel(axes[0],row,vaf,legend=True);vaf_panel(axes[1],row,vaf,early=True)
        axes[0].set_title('Full prediagnosis follow-up');axes[1].set_title('First 90 days');fig.suptitle(title);save(fig,f'patients/{pid}_VAF')
        fig,axes=plt.subplots(2,2,figsize=(10,6.5),layout='constrained')
        for j,index in enumerate([1,2]):
            trajectory_panel(axes[0,j],row,counts,index,legend=True);trajectory_panel(axes[1,j],row,counts,index,early=True)
            axes[0,j].set_title('Four-year model continuation');axes[1,j].set_title('First 90 days')
        init=initial(row);fig.suptitle(title+f"\nBaseline day {row['baseline_day']:g}: W={init['initial_W']:,}, C={init['initial_C']:,}, D={init['initial_D']}; 200 paths",fontsize=11)
        save(fig,f'patients/{pid}_trajectories')
        fig,axes=plt.subplots(1,2,figsize=(10,3.5),layout='constrained');risk_panel(axes[0],row,risk,legend=True);risk_panel(axes[1],row,risk,log=True);fig.suptitle(title);save(fig,f'patients/{pid}_tMN')
        fig,axes=plt.subplots(3,2,figsize=(10,11.2),layout='constrained')
        vaf_panel(axes[0,0],row,vaf);vaf_panel(axes[0,1],row,vaf,early=True)
        axes[0,0].set_title('VAF: observations and prediction');axes[0,1].set_title('VAF: first 90 days')
        trajectory_panel(axes[1,0],row,counts,1);trajectory_panel(axes[1,1],row,counts,2)
        axes[1,0].set_title('Single mutant C: all 200 paths');axes[1,1].set_title('Double hit D: all 200 paths')
        risk_panel(axes[2,0],row,risk);risk_panel(axes[2,1],row,risk,log=True)
        axes[2,0].set_title('Baseline-initialized tMN risk');axes[2,1].set_title('tMN risk on log10 scale')
        fig.suptitle(title+f"\nInitial W={init['initial_W']:,}, C={init['initial_C']:,}, D=0 at day {row['baseline_day']:g}",fontsize=12)
        save(fig,f'patients/{pid}_dashboard');print('PLOT',pid,flush=True)
    for kind in ['VAF','C','D','tMN']:
        fig,axes=plt.subplots(4,4,figsize=(13,10),layout='constrained')
        for ax,row in zip(axes.flat,COHORT):
            pid=row['patient_id']
            if kind=='VAF':vaf_panel(ax,row,vaf_all[pid])
            elif kind=='tMN':risk_panel(ax,row,risk_all[pid],log=True)
            else:trajectory_panel(ax,row,np.load(ROOT/'results/trajectories'/f'{pid}.npz')['counts'],1 if kind=='C' else 2)
            ax.set_title(pid.replace('_',' '),fontsize=9);ax.tick_params(labelsize=7);ax.xaxis.label.set_size(7);ax.yaxis.label.set_size(7)
            for text in ax.texts:text.set_fontsize(6)
        for ax in axes.flat[14:]:ax.set_visible(False)
        fig.suptitle({'VAF':'Patient VAF fits','C':'Patient single-mutant trajectories','D':'Patient double-hit trajectories','tMN':'Patient tMN risks'}[kind])
        save(fig,f'overview_{kind}')
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['core','likelihood','all'],default='all');args=parser.parse_args()
    verify_inputs()
    if args.stage in ['core','all']:gamma_figure();diagnostics();heatmap();patient_figures()
    if args.stage in ['likelihood','all']:likelihood_figures()
    print(f'Figures {args.stage} complete',flush=True)
if __name__=='__main__':main()
