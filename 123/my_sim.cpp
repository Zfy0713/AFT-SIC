#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
//#include <chrono>
#include <random>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

int n,k;
double **xx,*yy,*vv,*cumw,*cs,*f,*psi,*W,*B;
double *a, *temp;
int *delta;
#define SQR(x) ((x)*(x))
typedef struct
{
  int index;
  double v;
  double y;
  double B;
  int delta;
}
data_object;


void    swap(double *x,double *y);
double  criterion(int m, double alpha[]);
void    sort_alpha(int m, int n, double **xx, double alpha[], double vv[], double yy[],double B[], int delta[]);
void    gen_data(int m, int n, double alpha[], double **xx, double yy[], double vv[], int delta[], int seed);
void    convexmin(int n, double cumw[], double cs[], double y[]);
double  best_nearby (int m, double delta[], double point[], double prevbest,
                     double f(int m, double alpha[]), int *funevals);
int     hooke(int m, double startpt[], double endpt[], double rho,
              double eps, int itermax, double f(int m, double alpha[]));
int     CompareTime(const void *a, const void *b);
void    statistics(int n, int m, double **estimate, double *mean, double **covar);
double Normal(double z);
double normalCDF(double value);
void Normalize(int m, double alpha[]);

// [[Rcpp::export]]

List Compute_SSE(){
  int             i,j,m,iter,iter2,NumIt,seed, nIterations,EMmaxiter,sgn;
  double          *alpha,*alpha0,rho;
  double          **covar,*mean,**estimate, *alpha_init, *alpha_new;
  double          eps,sum;
  clock_t         StartTime, StopTime;
  double          Time_simulation;
  seed=159;
  
  n=1000;
  NumIt=20;
  m=5;
  
  EMmaxiter = 20;
  nIterations=1000;
  
  iter2=1000;
  eps=1.0e-10;
  
  
  xx = new double *[n];
  for (i=0;i<n;i++)
    xx[i] = new double [m];
  
  estimate = new double *[NumIt];
  for (i=0;i<NumIt;i++)
    estimate[i] = new double [m];
  
  covar = new double *[m];
  for (i=0;i<m;i++)
    covar[i] = new double [m];
  
  mean = new double[m];
  
  yy= new double[n];
  vv= new double[n];
  alpha= new double[m];
  alpha0= new double[m];
  alpha_init= new double[m];
  alpha_new = new double[m];
  f= new double[m];
  psi  = new double[n];
  delta = new int[n];
  cumw = new double[n+1];
  cs = new double[n+1];
  B = new double[n];
  W = new double[n];
  a = new double[n];
  temp = new double[n];
  
  
  std::mt19937 generator(seed+55);
  std::uniform_real_distribution<> dis(-2,2);
  std::normal_distribution<double> dis_normal(0.0,1);
  
  for (i=0;i<m;i++){
    alpha[i]=alpha_init[i]=1/sqrt(m);
    alpha0[i] = dis(generator);
  }
  Normalize(m, alpha0);
  
  StartTime = clock();
  Rcout << "True value:" ;
  for(i=0;i<m;i++){
    Rcout<<alpha0[i]<<" ";
  }
  Rcout<<endl;

  // Rcout << "     Iteration  " << "  alpha1  " << "      alpha2  "<< "       alpha3  " << std::endl << std::endl;
  
  
  for (iter=0;iter<NumIt;iter++)
  {
    
    gen_data(m, n, alpha0, xx, yy, vv, delta, seed+iter);
    std::mt19937 generator(seed+iter);
    for (i=0;i<m;i++){
      // alpha_init[i]=1/sqrt(m);
      // alpha_init[i] = alpha0[i] + dis_normal(generator);
      alpha_init[i] = dis(generator);
    }
    Normalize(m, alpha_init);
    Rcout <<setw(10)<<iter+1<<"th rep "<<"Init_Value:";
    for(i=0;i<m;i++){
      Rcout<< setprecision(6)<< alpha_init[i]<<" ";
    }
    Rcout << endl;
    rho=0.5;
    //-----------add:EM steps------------------
    double myeps = 1.0;
    double sum;
    for(i=0;i<n;i++){
      psi[i] = yy[i]; //initial value for unknown link function.
    }
    
    for(i=0;i<n;i++)
      B[i] = yy[i];
    
    k=1;
    for(k=1; k<EMmaxiter; k++){
      
      for(i=0;i<m;i++)
        alpha[i] = 0;
      
      iter2=hooke(m,alpha_init,alpha,rho,eps,nIterations,&criterion);
      
      Normalize(m, alpha);
      
      sum=0;
      for(i=0;i<m;i++){
        sum += SQR(alpha[i] - alpha_init[i]);
      }
      myeps = sum / m;
      
      for(i=0;i<m;i++){
        alpha_init[i] = alpha[i];
      }
      
      
      Rcout <<setw(10)<< k<<"th "<<"Iter: ";
      for (int i = 0; i < m; i++){
        Rcout << setprecision(6)<<setw(15)<< alpha_init[i] <<" ";
      }
      Rcout << " myeps:"<<myeps<<endl;

      if(myeps < 1e-8)
        break;
    }
    //-------------------------------------
    
    sum=0;
    for (i=0;i<m;i++)
      sum += SQR(alpha[i]);
    
    sgn = (alpha[0]>0) - (alpha[0]<0);
    for (j=0;j<m;j++)
      estimate[iter][j]= sgn * alpha[j]/sqrt(sum);
    
    Rcout <<endl;
    Rcout  << setw(10) << iter+1 <<"th replication:";
    for(i=0;i<m;i++){
      Rcout<< setprecision(6) <<  setw(15) << alpha[i];
    }
    Rcout << endl;
  }
  
  StopTime  = clock();
  Time_simulation   = (double)(StopTime - StartTime)/(double)CLOCKS_PER_SEC;
  
  Rcout << std::endl << std::endl;
  Rcout << "The computations took    " << setprecision(10) << Time_simulation << "   seconds"  << std::endl;
  
  statistics(NumIt,m,estimate,mean,covar);
  
  NumericMatrix out1 = NumericMatrix(NumIt,m);
  
  for (i=0;i<NumIt;i++)
  {
    for (j=0;j<m;j++)
      out1(i,j) = estimate[i][j];
  }
  
  NumericVector out2 = NumericVector(m);
  
  for (i=0;i<m;i++)
    out2(i)=mean[i];
  
  
  NumericMatrix out3 = NumericMatrix(m,m);
  
  for (i=0;i<m;i++)
  {
    for (j=0;j<m;j++)
      out3(i,j)=n*covar[i][j];
  }
  
  
  ofstream file0_("estimate.txt");
  
  if (file0_.is_open())
  {
    for (i=0;i<NumIt;i++)
    {
      for (j=0;j<m;j++)
        file0_ << setprecision(11) << setw(20) << estimate[i][j];
      file0_ << "\n";
    }
    file0_.close();
  }
  
  
  Rcout << "Making output list" << std::endl;
  
  // make the list for the output
  
  List out = List::create(Rcpp::Named("SSE")=out1,Rcpp::Named("means")=out2,Rcpp::Named("covariance_matrix")=out3);
  
  
  Rcout << "Freeing memory" << std::endl;
  
  // free memory
  
  delete[] mean, delete[] yy, delete[] vv, delete[] alpha0, delete[] alpha,
  delete[] alpha_init, delete[] f, delete[] psi, delete[] cumw, delete[] cs;
  delete[] W, delete[] B;
  for (i = 0;i < n;i++) delete[] xx[i];
  delete[] xx;
  
  for (i = 0;i < NumIt;i++) delete[] estimate[i];
  delete[] estimate;
  
  for (i = 0;i < m;i++) delete[] covar[i];
  delete[] covar;
  
  return out;
  
}

double Normal(double z){  return exp((-1)*z*z/2)/sqrt(2*PI);}

double normalCDF(double value){ return 0.5 * erfc(-value * M_SQRT1_2); }

void Normalize(int m, double alpha[]){
  int i,j;
  double sum=0;
  for(i=0;i<m;i++)
    sum += SQR(alpha[i]);
  for(i=0;i<m;i++)
    alpha[i] /= sqrt(sum);
}

void gen_data(int m, int n, double alpha[], double **xx, double yy[], double vv[], int delta[], int seed)
{
  int i,j;
  
  std::mt19937 generator(seed);
  std::uniform_real_distribution<> dis(0,1);
  std::normal_distribution<double> dis_normal(0.0,1.0);
  std::bernoulli_distribution b(0.6);
  
  for (i=0;i<n;i++)
  {
    for (j=0;j<m;j++)
      xx[i][j] = 1+dis(generator);
    
    yy[i]=0;
    
    for (j=0;j<m;j++)
      yy[i] += alpha[j]*xx[i][j];
    
    yy[i] = pow(yy[i],3)+dis_normal(generator);
    // yy[i] = std::sin(yy[i]) + dis_normal(generator);
    // yy[i] = pow(yy[i],2)+dis_normal(generator);
    delta[i] = b(generator);
    
  }
  
  sort_alpha(m,n,xx,alpha,vv,yy,B,delta);
}

void sort_alpha(int m, int n, double **xx, double alpha[], double vv[], double yy[], double B[], int delta[])
{
  int i,j,*ind;
  double **xx_new;
  data_object *obs;
  
  obs= new data_object[n];
  ind= new int[n];
  
  xx_new = new double *[n];
  for (i=0;i<n;i++)
    xx_new[i] = new double [m];
  
  for (i=0;i<n;i++)
  {
    vv[i]=0;
    for (j=0;j<m;j++)
      vv[i] += alpha[j]*xx[i][j];
  }
  
  for (i=0;i<n;i++)
  {
    obs[i].index=i;
    obs[i].v=vv[i];
    obs[i].y=yy[i];
    obs[i].delta = delta[i];
    obs[i].B = B[i];
  }
  
  qsort(obs,n,sizeof(data_object),CompareTime);
  
  for (i=0;i<n;i++)
    ind[i]=obs[i].index;
  
  
  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      xx_new[i][j]=xx[ind[i]][j];
  
  for (i=0;i<n;i++)
  {
    for (j=0;j<m;j++)
      xx[i][j]=xx_new[i][j];
    vv[i]=obs[i].v;
    yy[i]=obs[i].y;
    delta[i] = obs[i].delta;
    B[i] = obs[i].B;
  }
  
  delete[] obs;
  
  delete[] ind;
  for (i=0;i<n;i++)
    delete[] xx_new[i];
  delete[] xx_new;
}

void swap(double *x,double *y)
{
  double temp;
  temp=*x;
  *x=*y;
  *y=temp;
}

double criterion(int m, double alpha[])
{
  int i,j;
  double lambda, sum;
  
  //E-step---------------------------
  for (i = 0; i < n; i++)
    a[i] = B[i] - psi[i];
  
  for(i=0;i<n;i++){
    temp[i] = Normal(a[i]) / max((1-normalCDF(a[i])), 1e-8);
  }
  
  for(i=0;i<n;i++){
    if(delta[i]==1)
      W[i] = yy[i];
    if(delta[i]==0)
      W[i] = psi[i]+temp[i];
  }
  
  for(i=0;i<n;i++){
    yy[i] = W[i];
  }
  //---------------------------------
  
  sort_alpha(m,n,xx,alpha,vv,yy,B,delta);
  
  for (i=1;i<=n;i++)
  {
    cumw[i]=i*1.0;
    cs[i]=cs[i-1]+yy[i-1];
  }
  
  convexmin(n,cumw,cs,psi);
  
  for (j=0;j<m;j++)
    f[j]=0;
  
  for (j=0;j<m;j++)
  {
    for (i=0;i<n;i++)
      f[j] += xx[i][j]*(psi[i]-yy[i]);
  }
  
  lambda=0;
  
  for (i=0;i<m;i++)
    lambda += alpha[i]*f[i];
  
  for (i=0;i<m;i++)
    f[i] -= lambda*alpha[i];
  
  for (i=0;i<m;i++)
    f[i]/=n;
  
  sum=0;
  
  for (i=0;i<m;i++)
    sum += SQR(f[i]);
  
  return sqrt(sum);
  
}

void statistics(int n, int m, double **estimate, double *mean, double **covar)
{
  int i,j,k;
  double sum;
  
  for (j=0;j<m;j++)
  {
    sum=0;
    for (i=0;i<n;i++)
      sum += estimate[i][j];
    mean[j]=sum/n;
  }
  
  for (j=0;j<m;j++)
  {
    for (k=0;k<m;k++)
    {
      sum=0;
      for (i=0;i<n;i++)
        sum += (estimate[i][j]-mean[j])*(estimate[i][k]-mean[k]);
      covar[j][k] = sum/n;
    }
  }
}

void convexmin(int n, double cumw[], double cs[], double y[])
{
  int i, j, m;
  
  y[0] = cs[1]/cumw[1];
  for (i=2;i<=n;i++)
  {
    y[i-1] = (cs[i]-cs[i-1])/(cumw[i]-cumw[i-1]);
    if (y[i-2]>y[i-1])
    {
      j = i;
      while (y[j-2] > y[i-1] && j>1)
      {
        j--;
        if (j>1)
          y[i-1] = (cs[i]-cs[j-1])/(cumw[i]-cumw[j-1]);
        else
          y[i-1] = cs[i]/cumw[i];
        for (m=j;m<i;m++) y[m-1] = y[i-1];
      }
    }
  }
}

double best_nearby (int m, double delta[], double point[], double prevbest,
                    double f(int m, double alpha[]), int *funevals)
{
  double ftmp,minf,*z;
  int i;
  
  z = new double[m];
  
  minf = prevbest;
  
  for ( i = 0; i < m; i++ )
    z[i] = point[i];
  
  for ( i = 0; i < m; i++ )
  {
    z[i] = point[i] + delta[i];
    
    ftmp = f(m,z);
    *funevals = *funevals + 1;
    
    if ( ftmp < minf )
      minf = ftmp;
    else
    {
      delta[i] = - delta[i];
      z[i] = point[i] + delta[i];
      ftmp = f(m,z);
      *funevals = *funevals + 1;
      
      if ( ftmp < minf )
        minf = ftmp;
      else
        z[i] = point[i];
    }
  }
  
  for ( i = 0; i < m; i++ )
    point[i] = z[i];
  
  delete [] z;
  
  return minf;
}

int hooke(int m, double startpt[], double endpt[], double rho, double eps,
          int itermax, double f(int m, double alpha[]))
{
  double *delta,fbefore;
  int i,iters,keep,funevals,count;
  double newf,*newx,steplength,tmp;
  bool verbose = false;
  double *xbefore;
  
  delta = new double[m];
  newx = new double[m];
  xbefore = new double[m];
  
  for ( i = 0; i < m; i++ )
    xbefore[i] = newx[i] = startpt[i];
  
  for ( i = 0; i < m; i++ )
  {
    if ( startpt[i] == 0.0 )
      delta[i] = rho;
    else
      delta[i] = rho*fabs(startpt[i]);
  }
  
  funevals = 0;
  steplength = rho;
  iters = 0;
  
  
  fbefore = f(m,newx);
  funevals = funevals + 1;
  newf = fbefore;
  
  while ( iters < itermax && eps < steplength )
  {
    iters = iters + 1;
    
    if (verbose)
    {
      cout << "\n";
      cout << "  FUNEVALS, = " << funevals
           << "  F(X) = " << fbefore << "\n";
      
      for ( i = 0; i < m; i++ )
      {
        cout << "  " << i + 1
             << "  " << xbefore[i] << "\n";
      }
    }
    //
    //  Find best new alpha, one coordinate at a time.
    //
    for ( i = 0; i < m; i++ )
      newx[i] = xbefore[i];
    
    
    
    newf = best_nearby(m,delta,newx,fbefore,f,&funevals);
    //
    //  If we made some improvements, pursue that direction.
    //
    keep = 1;
    count=0;
    
    while (newf<fbefore && keep == 1 && count<=100)
    {
      count++;
      for ( i = 0; i < m; i++ )
      {
        //
        //  Arrange the sign of DELTA.
        //
        if ( newx[i] <= xbefore[i] )
          delta[i] = - fabs(delta[i]);
        else
          delta[i] = fabs(delta[i]);
        //
        //  Now, move further in this direction.
        //
        tmp = xbefore[i];
        xbefore[i] = newx[i];
        newx[i] = newx[i] + newx[i] - tmp;
      }
      
      fbefore = newf;
      
      newf = best_nearby(m,delta,newx,fbefore,f,&funevals);
      //
      //  If the further (optimistic) move was bad...
      //
      if (fbefore <= newf)
        break;
      //
      //  Make sure that the differences between the new and the old points
      //  are due to actual displacements; beware of roundoff errors that
      //  might cause NEWF < FBEFORE.
      //
      keep = 0;
      
      for ( i = 0; i < m; i++ )
      {
        if ( 0.5 * fabs(delta[i]) < fabs(newx[i]-xbefore[i]))
        {
          keep = 1;
          break;
        }
      }
    }
    
    if (eps <= steplength && fbefore <= newf)
    {
      steplength = steplength * rho;
      for ( i = 0; i < m; i++ )
        delta[i] = delta[i] * rho;
    }
    
  }
  
  for ( i = 0; i < m; i++ )
    endpt[i] = xbefore[i];
  
  delete [] delta;
  delete [] newx;
  delete [] xbefore;
  
  return iters;
}

int CompareTime(const void *a, const void *b)
{
  if ((*(data_object *) a).v < (*(data_object *) b).v)
    return -1;
  if ((*(data_object *) a).v > (*(data_object *) b).v)
    return 1;
  return 0;
}

/***R
Compute_SSE()
*/

