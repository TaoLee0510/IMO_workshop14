// Exact inhomogeneous three-type Moran simulation with state-specific thinning.
// Type-preserving replacement events are analytically skipped, not approximated.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

static uint64_t splitmix(uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
int main(int argc,char**argv){
 try{
  if(argc!=15)throw std::runtime_error("prefix N b sd mu area alpha theta C0 D0 horizon replicates seed threshold_fraction");
  const std::string prefix=argv[1];const int N=std::stoi(argv[2]);
  const double b=std::stod(argv[3]),sd=std::stod(argv[4]),mu=std::stod(argv[5]),area=std::stod(argv[6]);
  const double alpha=std::stod(argv[7]),theta=std::stod(argv[8]);
  const int c0=std::stoi(argv[9]),d0=std::stoi(argv[10]),days=std::stoi(argv[11]),reps=std::stoi(argv[12]);
  const uint64_t master=std::stoull(argv[13]);const int K=int(std::floor(std::stod(argv[14])*N))+1;
  if(N<2||c0<0||d0<0||c0+d0>N||b<=0||sd<0||mu<0||mu>1||alpha<=1||theta<=0||area<0)throw std::runtime_error("invalid parameters");
  const double lgamma=std::lgamma(alpha),logtheta=std::log(theta),mode=(alpha-1)*theta;
  const double max_extra=area==0?0:area*std::exp((alpha-1)*std::log(mode)-mode/theta-lgamma-alpha*logtheta);
  const double fmax=1+max_extra;
  auto gamma_extra=[&](double time){return time<=0||area==0?0.:area*std::exp((alpha-1)*std::log(time)-time/theta-lgamma-alpha*logtheta);};
  double cutoff=std::max(mode,theta);
  while(gamma_extra(cutoff)>1e-17)cutoff=2*cutoff+theta;
  double low=mode,high=cutoff;
  for(int i=0;i<60;++i){double middle=(low+high)/2;if(gamma_extra(middle)>1e-17)low=middle;else high=middle;}
  cutoff=high;
  std::ofstream raw(prefix+".bin",std::ios::binary),summary(prefix+".csv");
  if(!raw||!summary)throw std::runtime_error("cannot open output");
  summary<<"replicate,seed,first_hit_day,ever_D,final_C,final_D,proposals,type_changes\n";
  const auto begun=std::chrono::steady_clock::now();
  uint64_t proposals_all=0,changes_all=0;
  for(int rep=0;rep<reps;++rep){
   const uint64_t seed=splitmix(master+uint64_t(rep));std::mt19937_64 gen(seed);
   auto uniform=[&](){return (double((gen()>>11))+0.5)*0x1.0p-53;};
   int c=c0,d=d0,nextday=0;double t=0,first=d>=K?0:-1;bool ever=d>0;
   uint64_t proposals=0,changes=0;std::vector<int32_t> series((days+1)*2);
   auto save=[&](){series[2*nextday]=c;series[2*nextday+1]=d;++nextday;};save();
   while(t<days){
    const int w=N-c-d;
    // Gamma is decreasing after its mode: its current value then bounds
    // every future proposal time. This avoids a loose peak bound for years.
    const double future_fmax=t<mode?fmax:(t>=cutoff?1.:1+gamma_extra(t));
    const double aw=w,ac=(1-mu)*future_fmax*c,ad=(future_fmax+sd)*d+mu*future_fmax*c;
    const double bound=b*(aw*(N-w)+ac*(N-c)+ad*(N-d))/(N+sd*d);
    if(bound<=0){while(nextday<=days)save();break;}
    const double proposed=t-std::log(uniform())/bound;
    while(nextday<=days&&nextday<=proposed)save();
    if(proposed>days)break;t=proposed;++proposals;
    const double extra=t>=cutoff?0.:gamma_extra(t);
    const double fc=1+extra,fd=fc+sd,F=w+fc*c+fd*d;
    const double pw=w/F,pc=(1-mu)*fc*c/F,pd=(fd*d+mu*fc*c)/F;
    const double r[6]={b*w*pc,b*w*pd,b*c*pw,b*c*pd,b*d*pw,b*d*pc};
    const double total=r[0]+r[1]+r[2]+r[3]+r[4]+r[5];
    if(total>bound*(1+1e-10))throw std::runtime_error("invalid thinning bound");
    double u=uniform()*bound;int transition=0;
    for(;transition<6&&u>=r[transition];++transition)u-=r[transition];
    if(transition==6)continue;
    switch(transition){case 0:++c;break;case 1:++d;break;case 2:--c;break;case 3:--c;++d;break;case 4:--d;break;case 5:++c;--d;break;}
    ++changes;if(c<0||d<0||c+d>N)throw std::runtime_error("population conservation failure");
    ever=ever||d>0;if(first<0&&d>=K)first=t;
   }
   while(nextday<=days)save();
   raw.write(reinterpret_cast<const char*>(series.data()),series.size()*sizeof(int32_t));
   summary<<rep<<','<<seed<<',';if(first>=0)summary<<std::setprecision(17)<<first;
   summary<<','<<int(ever)<<','<<c<<','<<d<<','<<proposals<<','<<changes<<'\n';
   proposals_all+=proposals;changes_all+=changes;
  }
  std::cout<<"{\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begun).count()
           <<",\"proposals\":"<<proposals_all<<",\"type_changes\":"<<changes_all<<",\"replicates\":"<<reps<<"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
 return 0;
}
