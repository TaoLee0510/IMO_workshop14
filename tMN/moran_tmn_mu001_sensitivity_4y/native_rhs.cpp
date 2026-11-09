// The same tridiagonal forward derivative as model.py, with no approximation.
// Compile without FMA contraction to preserve the NumPy operation order.
#include <vector>
#include <cstddef>
extern "C" void tmn_forward_rhs(int n,const double*y,double*out,double low,
  double N,double b,double mu,double fc,double fd,double q,double up,double down,double eta_dot,int killed){
  thread_local std::vector<double> birth,death;birth.resize(n);death.resize(n);
  for(int i=0;i<n;++i){
    const double d=low+i,remaining=N-d,c=remaining*q;
    const double denominator=remaining+(fc-1)*c+fd*d;
    const double pd=(fd*d+mu*fc*c)/denominator;
    birth[i]=b*remaining*pd;death[i]=b*d*(1-pd);
  }
  if(!killed){birth[n-1]=0;death[0]=0;}
  for(int i=0;i<n;++i)out[i]=(-(birth[i]+death[i])+eta_dot*(low+i))*y[i];
  for(int i=0;i<n-1;++i)out[i+1]+=(up*birth[i])*y[i];
  for(int i=0;i<n-1;++i)out[i]+=(down*death[i+1])*y[i+1];
  if(killed)out[n]=(up*birth[n-1])*y[n-1]+(eta_dot*(low+n))*y[n];
}
