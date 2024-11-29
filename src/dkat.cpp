#include <RcppArmadillo.h>
#include "linreg.h"
#include "svd.h"
#include "crossprod.h"
#include <Rcpp.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
int sign(double x) {
  return (x > 0) - (x < 0);
}
// [[Rcpp::export(.dkat_statistic)]]
long double dkat_statistic(arma::mat Ux,arma::mat dx, arma::mat Uy, arma::mat dy) {
  arma::mat P = Ux.t() * Uy;
  long double Fstar = arma::trace(dx * P * dy * P.t());
  return(Fstar);
}

// [[Rcpp::export(.dkat_speed_cpp)]]
long double dkat(arma::mat Ux,arma::mat dx, arma::mat Uy, arma::mat dy) {
  long int n = Uy.n_rows;

  arma::mat W = Ux * dx * Ux.t(); // / trace(K*K);
  arma::mat A = Uy * dy * Uy.t(); // / trace(L*L);

  long double Fstar = arma::trace(A*W);
  long double T = arma::accu(dy);
  double Ts = arma::accu(dx);
  long double meankrv = T * Ts / (n-1);
  double T2 = arma::accu(pow(dy.diag(), 2));//arma::trace(A * A);
  long double S2 = arma::accu(arma::pow(A.diag(), 2));
  long double T2s =  arma::accu(pow(dx.diag(), 2));
  double S2s = arma::accu(arma::pow(W.diag(), 2));

  long double temp1 = 2*((n-1)*T2-pow(T,2))*((n-1)*T2s-pow(Ts,2))/pow((n-1),2)/(n+1)/(n-2);
  long double temp21 = n*(n+1)*S2- (n-1)*(pow(T,2)+2*T2);
  long double temp22 = n*(n+1)*S2s- (n-1)*(pow(Ts,2)+2*T2s);
  // long int temp23 = (n+1) * n * (n-1) * (n-2) * (n-3);
  long double temp2=log(abs(temp21)) + log(abs(temp22)) - log(n+1.0) - log(n) - log(n-1) -log(n-2) -log(n-3);
  temp2 = exp(temp2);
  temp2 = sign(temp21)*sign(temp22) * temp2;
  // long double temp2 = temp21*temp22/temp23;
  long double variancekrv=temp1+temp2;

  arma::mat A2 = A * A;
  long double T3=arma::accu(pow(dy.diag(), 3));
  long double S3=arma::accu(pow(A.diag(),3));
  long double U = arma::accu(pow(A,3));
  long double R= arma::sum(A.diag().t() * A2.diag()); // fix up?
  long double B= arma::sum(A.diag().t() * A * A.diag());

  arma::mat W2 = W * W;
  long double T3s=arma::accu(pow(dx.diag(), 3)); long double S3s=arma::accu(pow(W.diag(),3));
  long double Us = arma::accu(pow(W,3));
  long double Rs= arma::sum(W.diag().t() * W2.diag()); // fix up?
  long double Bs= arma::sum(W.diag().t() * W * A.diag());

  long double t1=pow(n,2)*(n+1)*(pow(n,2)+15*n-4)*S3*S3s;
  long double t2=4*(pow(n,4)-8*pow(n,3)+19*pow(n,2)-4*n-16)*U*Us;
  long  double t3=24*(pow(n,2)-n-4) * (U*Bs+B*Us);
  long  double t4=6*(pow(n,4)-8*pow(n,3)+21*pow(n,2)-6*n-24)*B*Bs;

  long double t5=12*(pow(n,4)-pow(n,3)-8*pow(n,2)+36*n-48)*R*Rs;
  long double t6=12*(pow(n,3)-2*pow(n,2)+9*n-12)*(T*S2*Rs+R*Ts*S2s);
  long double t7=3*(pow(n,4)-4*pow(n,3)-2*pow(n,2)+9*n-12)*T*Ts*S2*S2s;

  long double t81=(pow(n, 3)-3*pow(n,2)-2*n+8)*(R*Us+U*Rs);
  long double t82=(pow(n,3)-2*pow(n,2)-3*n+12)*(R*Bs+B*Rs);
  long double t8=24*(t81+t82);
  long double t9=12*(pow(n,2)-n+4)*(T*S2*Us+U*Ts*S2s);
  long double t10=6*(2*pow(n,3)-7*pow(n,2)-3*n+12)*(T*S2*Bs+B*Ts*S2s);
  long  double t11=-2*n*(n-1)*(pow(n,2)-n+4)*((2*U+3*B)*S3s+(2*Us+3*Bs)*S3);

  long double t12=-3*n*pow(n-1,2)*(n+4)*((T*S2+4*R)*S3s+(Ts*S2s+4*Rs)*S3);
  long double  t13=2*n*(n-1)*(n-2)*((pow(T,3)+6*T*T2+8*T3)*S3s+(pow(Ts,3)+6*Ts*T2s+8*T3s)*S3);

  long double t14=pow(T,3)*((pow(n,3)-9*pow(n,2)+23*n-14)*pow(Ts,3)+6*(n-4)*Ts*T2s+8*T3s);
  long double t15=6*T*T2*((n-4)*pow(Ts,3)+(pow(n,3)-9*pow(n,2)+24*n-14)*Ts*T2s+4*(n-3)*T3s);
  long double t16=8*T3*(pow(Ts,3)+3*(n-3)*Ts*T2s+(pow(n,3)-9*pow(n,2)+26*n-22)*T3s);
  long double t17=-16*(pow(T,3)*Us+U*pow(Ts,3))-6*(T*T2*Us+U*Ts*T2s)*(2*pow(n,2)-10*n+16);
  long double t18=-8*(T3*Us+U*T3s)*(3*pow(n,2)-15*n+16)-(pow(T,3)*Bs+B*pow(Ts,3))*(6*pow(n,2)-30*n+24);
  long double t19=-6*(T*T2*Bs+B*Ts*T2s)*(4*pow(n,2)-20*n+24)-8*(T3*Bs+B*T3s)*(3*pow(n,2)-15*n+24);
  long double t201=24*(pow(T,3)*Rs+R*pow(Ts,3))+6*(T*T2*Rs+R*Ts*T2s)*(2*pow(n,2)-10*n+24);
  long double t202=8*(T3*Rs+R*T3s)*(3*pow(n,2)-15*n+24)+(3*pow(n,2)-15*n+6)*(pow(T,3)*Ts*S2s+T*S2*pow(Ts,3));
  long double t203=6*(T*T2*Ts*S2s+Ts*T2s*T*S2)*(pow(n,2)-5*n+6)+48*(T3*Ts*S2s+T3s*T*S2);
  long  double t20=-(n-2)*(t201+t202+t203);
  long double temp31=t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20;

  // long long int temp32=n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5);
  long double mom3=log(temp31)- log(n) - log(n-1) - log(n-2) -log(n-3) -log(n-4) - log(n-5); //temp32;
  mom3 = exp(mom3);
  long double skewnesskrv= (mom3-3*meankrv*variancekrv-pow(meankrv,3))/pow(variancekrv,1.5);
  long double shape = 4/pow(skewnesskrv,2);
  long double location = meankrv - 2 * sqrt(variancekrv) / skewnesskrv;
  long double scale=sqrt(variancekrv)*skewnesskrv/2;
  Environment pkg = Environment::namespace_env("PearsonDS");
  Function f = pkg["ppearsonIII"];
  SEXP p = f(Fstar, Named("shape")=shape, _["location"]=location, _["scale"] = scale, _["lower.tail"] = false);
  NumericVector pp(p);
  return(pp[0]);
}
