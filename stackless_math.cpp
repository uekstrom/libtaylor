#include <cmath>
#include <cassert>
#include <iostream>

using namespace std;

/* pow and variants */

/* x^p */
template<typename T>
void taylor_pow_forwards(T *out, const T &x0, int ndeg, double p)
{
  out[0] = pow(x0,p);
  for (int k=0;k<ndeg;k++)
    out[k+1] = ((p-k)/(k+1))*out[k]/x0;
}

/*
  Forward recursion seems better but this one can accept x0=0.
  What to do?
 */
template<typename T>
void taylor_pow_backwards(T *out, const T &x0, int ndeg, double p)
{
  int kmax = ndeg;
  if (p == int(p) and p < ndeg)
    {
      kmax = p;
      for (int k=kmax+1;k<=ndeg;k++)
	out[k] = 0;
    }
  double fac = 1;
  for (int k=1;k<=kmax;k++)
    fac *= (p-k+1)/k;
  out[kmax] = fac*pow(x0,p-kmax);
  for (int k=kmax-1;k>=0;k--)
    out[k] = ((k+1)/(p-k))*x0*out[k+1];
}

template<typename T>
void taylor_pow(T *out, const T &x0, int ndeg, double p)
{
  if (p == int(p) and p >= 0) // backwards recursion
    {
      int kmax = p;
      if (kmax > ndeg)
	kmax = ndeg;
      for (int k=kmax+1;k<=ndeg;k++)
	out[k] = 0;
      double fac = 1;
      for (int k=1;k<=kmax;k++)
	fac *= (p-k+1)/k;
      out[kmax] = fac*pow(x0,p-kmax);
      for (int k=kmax-1;k>=0;k--)
	out[k] = (k+1)/(p-k)*x0*out[k+1];
    }
  else // forward recursion
    {
      out[0] = pow(x0,p);
      for (int k=0;k<ndeg;k++)
	out[k+1] = ((p-k)/(k+1))*out[k]/x0;
    }
}

template<typename T>
void taylor_sqrt(T *out, const T &x0, int ndeg)
{
  out[0] = sqrt(x0);
  for (int k=0;k<ndeg;k++)
    out[k+1] = ((0.5-k)/(k+1))*out[k]/x0;
}

template<typename T>
void taylor_cbrt(T *out, const T &x0, int ndeg)
{
  out[0] = cbrt(x0);
  if (ndeg<1)
    return;
  for (int k=0;k<ndeg;k++)
    out[k+1] = ((1.0/3.0-k)/(k+1))*out[k]/x0;
}

#if 0
/* (a+bx^2)^p, the derivative of many functions */
template<typename T>
void taylor_pow_q(T *out, const T &x0, T tmp[2], int ndeg, double a, double b, double p)
{
  tmp[0] = a + b*x0*x0;
  out[0] = pow(tmp[0],p);
  if (ndeg<1)
    return;
  tmp[1] = x0/tmp[0];
  out[1] = 2*b*p*tmp[1]*out[0];
  if (ndeg<2)
    return;
  for (int k=1;k<ndeg;k++)
    out[k+1] = (2*(p-k)*b/(k+1))*tmp[1]*out[k] + ((2*p-k+1)*b/(k+1))*out[k-1]/tmp[0];
}
#endif

/* integral (a+bx^2)^p with unset constant term */
template<typename T>
void taylor_integral_pow_q(T *out, const T &x0, T tmp[2], int ndeg, double a, double b, double p)
{
  if (ndeg<=0)
    return;
  tmp[0] = a + b*x0*x0;
  out[1] = pow(tmp[0],p);
  tmp[1] = x0/tmp[0];
  if (ndeg<2)
    return;
  out[2] = b*p*tmp[1]*out[1];
  if (ndeg<3)
    return;
  for (int k=1;k<ndeg-1;k++)
    out[k+2] = (2*(p-k)*b/(k+2))*tmp[1]*out[k+1] + (k*(2*p-k+1)*b/((k+2)*(k+1)))*out[k]/tmp[0];
}


// specialization of ..pow_q for p = -0.5
template<typename T>
void taylor_integral_invsqrt_q(T *out, const T &x0, T tmp[2], int ndeg, double a, double b)
{
  if (ndeg<=0)
    return;
  tmp[0] = a + b*x0*x0;
  tmp[1] = 1/tmp[0];
  out[1] = sqrt(tmp[1]);
  tmp[1] *= x0;
  if (ndeg<2)
    return;
  out[2] = -0.5*b*tmp[1]*out[1];
  if (ndeg<3)
    return;
  for (int k=1;k<ndeg-1;k++)
    out[k+2] = (2*(-0.5-k)*b/(k+2))*tmp[1]*out[k+1] + (k*(-k)*b/((k+2)*(k+1)))*out[k]/tmp[0];
}

// specialization of ..pow_q for p = -1
template<typename T>
void taylor_integral_inv_q(T *out, const T &x0, T tmp[2], int ndeg, double a, double b)
{
  // p = -1
  if (ndeg<=0)
    return;
  tmp[0] = a + b*x0*x0;
  out[1] = 1/tmp[0];
  tmp[1] = x0/tmp[0];
  if (ndeg<2)
    return;
  out[2] = -b*tmp[1]*out[1];
  if (ndeg<3)
    return;
  for (int k=1;k<ndeg-1;k++)
    out[k+2] = (2*(-1-k)*b/(k+2))*tmp[1]*out[k+1] + (k*(-k-1)*b/((k+2)*(k+1)))*out[k]/tmp[0];
}


/* Inverse trig functions */

template<typename T>
void taylor_atan(T *out, const T &x0, T tmp[2], int ndeg)
{
  taylor_integral_inv_q(out,x0,tmp,ndeg,1.0,1.0);
  out[0] = atan(x0);
}

template<typename T>
void taylor_asin(T *out, const T &x0, T tmp[2], int ndeg)
{
  taylor_integral_invsqrt_q(out,x0,tmp,ndeg,1.0,-1.0);
  out[0] = asin(x0);
}

template<typename T>
void taylor_acos(T *out, const T &x0, T tmp[2], int ndeg)
{
  taylor_integral_invsqrt_q(out,x0,tmp,ndeg,1.0,-1.0);
  out[0] = acos(x0);
  for (int k=1;k<=ndeg;k++) // TODO: Build this functionality into the integral_pow_q function
    out[k] *= -1;
}

template<typename T>
void taylor_asinh(T *out, const T &x0, T tmp[2], int ndeg)
{
  taylor_integral_invsqrt_q(out,x0,tmp,ndeg,1.0,1.0);
  out[0] = asinh(x0);
}

template<typename T>
void taylor_acosh(T *out, const T &x0, T tmp[2], int ndeg)
{
  taylor_integral_invsqrt_q(out,x0,tmp,ndeg,-1.0,1.0);
  out[0] = acosh(x0);
}

template<typename T>
void taylor_atanh(T *out, const T &x0, T tmp[2], int ndeg)
{
  taylor_integral_inv_q(out,x0,tmp,ndeg,1.0,-1.0);
  out[0] = atanh(x0);
}

/* Trig functions */

template<typename T>
void taylor_sin(T *out, const T &x0, int ndeg)
{
  out[0] = sin(x0);
  if (ndeg > 0)
    out[1] = cos(x0);
  double fac = 1.0;
  for (int k=2;k<=ndeg;k++)
    {
      fac /= k;
      out[k] = fac*(1 - 2*((k/2) % 2))*out[k % 2];
    }
}

template<typename T>
void taylor_cos(T *out, const T &x0, int ndeg)
{
  out[0] = cos(x0);
  if (ndeg > 0)
    out[1] = -sin(x0);
  double fac = 1.0;
  for (int k=2;k<=ndeg;k++)
    {
      fac /= k;
      out[k] = fac*(1 - 2*((k/2) % 2))*out[k % 2];
    }
}

template<typename T>
void taylor_tan(T *out, const T &x0, int ndeg)
{
  out[0] = tan(x0);
  for (int k=0;k<ndeg;k++)
    {
      out[k+1] = 1/(k+1);
      for (int i=0;i<=k;i++)
	out[k+1] += out[i]*out[k-i]/(k+1);
    }
}

/* exp, log and erf */

template<typename T>
void taylor_exp(T *out, const T &x0, int ndeg)
{
  out[0] = exp(x0);
  for (int k=0;k<ndeg;k++)
    out[k+1] = out[k]/(k+1);
}

template<typename T>
void taylor_expm1(T *out, const T &x0, int ndeg)
{
  out[0] = expm1(x0);
#ifdef TAYLOR_FAST
  if (ndeg>0)
    out[1] = out[0]+1;
#else
  if (ndeg>0)
    out[1] = exp(x0);
#endif
  for (int k=1;k<ndeg;k++)
    out[k+1] = out[k]/(k+1);
}

template<typename T>
void taylor_log(T *out, const T &x0, int ndeg)
{
  out[0] = log(x0);
  if (ndeg < 1)
    return;
  out[1] = 1/x0;
  if (ndeg < 2)
    return;
  for (int k=0;k<ndeg-1;k++)
    out[k+2] = (-(k+1)/double(k+2))*out[k+1]/x0;
}

template<typename T>
void taylor_erf(T *out, const T &x0, T tmp[1], int ndeg)
{
  out[0] = erf(x0);
  if (ndeg < 1)
    return;
  tmp[0] = -x0*x0;
  out[1] = 2/sqrt(M_PI)*exp(tmp[0]);
  if (ndeg < 2)
    return;
  for (int k=0;k<ndeg-1;k++)
    out[k+2] = (-2.0/(k+2))*out[k+1]*x0 + (k*-2.0/((k+2)*(k+1)))*out[k];
}

template<typename T>
void taylor_erfc(T *out, const T &x0, T tmp[1], int ndeg)
{
  out[0] = erfc(x0);
  if (ndeg < 1)
    return;
  tmp[0] = -x0*x0;
  out[1] = -2/sqrt(M_PI)*exp(tmp[0]);
  if (ndeg < 2)
    return;
  for (int k=0;k<ndeg-1;k++)
    out[k+2] = (-2.0/(k+2))*out[k+1]*x0 + (k*-2.0/((k+2)*(k+1)))*out[k];
}

/* exp(p*x**2) */
template<typename T>
void taylor_gauss(T *out, const T &x0, T tmp[1], int ndeg, double p)
{
  tmp[0] = p*x0*x0;
  out[0] = exp(tmp[0]);
  if (ndeg < 1)
    return;
  out[1] = 2*p*x0*out[0];
  if (ndeg < 2)
    return;
  for (int k=1;k<ndeg;k++)
    out[k+1] = 2*p/(k+1)*out[k]*x0 + 2*p/(k+1)*out[k-1];
}



#if 0
// This is for vectorization use. Must be checked for ndeg=0!
/*
  Upward-downward recursive method for f(x) = asinh(sqrt(x))/sqrt(x)
  This is the vectorizable version, a branching version should be
  made separately.


  The function satisfies the recurrence relation
  2x f^(n+1) = h^(n) - (2n+1)f^(n), where h(x) = 1/sqrt(1+x).
  This is upwards unstable for small x. To use downward recursion
  one can compute the derivative f^(n) directly using the theory of 
  Hypergeometric functions. We have that f(x) = 2F1(1/2,1/2;3/2;-x).
  This Taylor series converges more slowly for higher derivatives,
  and should not be used. 
  Using 15.1.7 and 15.2.6 in Abramovitz and Stegun we can write
  f^(n) = (-1)^n*(1/2)_n^2/(3/2)_n*(1+x)^(1/2-n) 2F1(1,1;3/2+n;-x),
  where (.)_n is the rising Pochhammer symbol.
  This hypergeometric function converges more rapidly for larger
  n, and is suitable for use. 
  
  The strategy is to use downward recursion for small x,
  and upward for large. The value f^(0) can be accurately 
  computed using the original definition of the function for 
  all positive x, so always use that.

  If branches are acceptable some small savings can be had by testing
  x at the top of the function and skip the taylor part if not needed.
 */
template<typename T>
void taylor_asinh_sqrtx_div_sqrtx_interleaved(T *out, const T &x0, T tmp[2], int ndeg)
{
#ifdef TAYLOR_FAST
  // These constants are appropriate for double precision.
  // Roughly xc = exp(-log(E)/nterms), where E is the relative error wanted.
  const double xc = 0.01; const int nterms = 6; // Less accurate but good up to order 3.
#else
  const double xc = 0.1; const int nterms = 12; // More accurate and to higher orders.
  // For very high accuracy use xc=0.5 and nterms=40.
#endif
  /* 
     assert x0>0.
     TODO: Taylor coeffs instead of derivatives!
   */
#ifdef TAYLOR_ASINH_NO_LOG
  tmp[0] = sqrt(x0); // here we can use x0+1e-40 or something to allow for x0=0.
  out[0] = asinh(tmp[0]); // possibly more accurate than log(sqrt(x)+sqrt(1+x)) for small x.
  out[0] /= tmp[0];
  if (ndeg<1)
    return;
  tmp[0] = 1 + x0;
  tmp[1] = sqrt(tmp[0]);   // accumulate fn into tmp[1]
  tmp[0] = 1/tmp[0];
#else
  tmp[0] = 1 + x0;
  tmp[1] = sqrt(tmp[0]);   // accumulate fn into tmp[1]
  tmp[0] = sqrt(x0); // here we can use x0+1e-40 or something to allow for x0=0. 
  out[0] = tmp[0] + tmp[1];
  out[0] = log(out[0]);
  out[0] /= tmp[0];
  if (ndeg<1)
    return;
  tmp[0] = 1/(1+x0);
#endif
  if (ndeg > 1)
    {
      out[1] = -0.5*tmp[0]/tmp[1]; //h1. Now tmp[1] is free
      tmp[1] = 0.5/tmp[1];
      tmp[1] -= 0.5*out[0];
      tmp[1] /= x0; // tmp[1] contains f1 here
      for (int n=1;n<ndeg-1;n++)
	{
	  tmp[1] = (out[n]-(2*n+1)*tmp[1])/(2*(n+1)*x0);
	  out[n+1] = -(n+0.5)/(n+1.0)*tmp[0]*out[n];  // Now hn is in out[1]..out[ndeg-1]
	}
      out[ndeg] = (out[ndeg-1] - (2*ndeg-1)*tmp[1])/(2*ndeg*x0);
    }
  else
    {
      tmp[1] = 0.5/tmp[1];
      out[1] = (tmp[1]-0.5*out[0])/x0; // put f1 in out[1]
    }
  // Coefficient conversion until here
  // Taylor part. This takes about 30% of the cpu time with -O2 and 12 terms, 15% with 6 terms.
  // If we can branch we can just skip the series when x0<xc.
  if (ndeg>1)
    tmp[0] = (-(ndeg-0.5)/(ndeg*(2*ndeg+1.0)))*out[ndeg-1];
  else
    tmp[0] = (-1.0/3.0)*tmp[1];
  tmp[1] = tmp[0];
  for (int k=1; k < nterms; k++ )
    {
      tmp[0] *= (-k/(0.5+ndeg+k))*x0;
      tmp[1] += tmp[0];
    }
  // Pick best fn
  if (x0 < xc)
    out[ndeg] = tmp[1];
  for (int n=ndeg-1;n>0;n--)
    out[n] = (out[n] - 2*x0*(n+1)*out[n+1])/(2*n+1);
  if (x0 < xc)
    {
      out[0] = (1/sqrt(1+x0) - 2*x0*out[1]);
      // Use the Taylor recursion for out[0] for small x to allow for the fast expression asinh(sqrt(x)) = log(sqrt(x)+sqrt(1+x))
    }
}
#endif


/*
  Upward-downward recursive method for f(x) = asinh(sqrt(x))/sqrt(x)
  
  The function satisfies the recurrence relation
  2x f^(n+1) = h^(n) - (2n+1)f^(n), where h(x) = 1/sqrt(1+x).
  This is upwards unstable for small x. To use downward recursion
  one can compute the derivative f^(n) directly using the theory of 
  Hypergeometric functions. We have that f(x) = 2F1(1/2,1/2;3/2;-x).
  This Taylor series converges more slowly for higher derivatives,
  and should not be used. 
  Using 15.1.7 and 15.2.6 in Abramovitz and Stegun we can write
  f^(n) = (-1)^n*(1/2)_n^2/(3/2)_n*(1+x)^(1/2-n) 2F1(1,1;3/2+n;-x),
  where (.)_n is the rising Pochhammer symbol.
  This hypergeometric function converges more rapidly for larger
  n, and is suitable for use. 
  
  The strategy is to use downward recursion for small x,
  and upward for large. The value f^(0) can be accurately 
  computed using the original definition of the function for 
  all positive x, but the asinh function is typically much
  slower than log. Therefore use Taylor expansion for f^(0)
  for small x, and use asinh(sqrt(x)) = log(sqrt(x)+sqrt(1+x)) 
  for large x.

 */
#define TAYLOR_FAST
template<typename T>
void taylor_asinh_sqrtx_div_sqrtx(T *out, const T &x0, T tmp[2], int ndeg)
{
#ifdef TAYLOR_FAST
  // These constants are appropriate for double precision.
  // Roughly xc = exp(-log(E)/nterms), where E is the relative error wanted.
  const double xc = 0.01; const int nterms = 6; // Less accurate but good up to order 3.
#else
  const double xc = 0.1; const int nterms = 14; // More accurate and to higher orders.
  // For very high accuracy use xc=0.5 and nterms=40.
#endif
  if (x0 <= xc) // At order 2 with nterms = 14 this branch takes 30% longer 
                // than the upwards branch, smaller relative difference for high orders.
    { // Downward recursion from Taylor. Put the 1/sqrt(1+x) coefficients in out for future use.
      out[0] = 1/sqrt(1+x0);
      for (int n=0;n<ndeg-1;n++)
	out[n+1] = -(n+0.5)/(n+1.0)/(1+x0)*out[n];
      double t;  // Do the Taylor sum
      if (ndeg>0)
	t = (0.5-ndeg)/(ndeg*(2*ndeg+1.0))*out[ndeg-1];
      else
	t = 1/out[0];
      out[ndeg] = t;
      for (int k=1; k < nterms; k++ )
	{
	  t *= (-k/(0.5+ndeg+k))*x0;
	  out[ndeg] += t;
	}      
      for (int n=ndeg-1;n>=0;n--) // Recurse down
	out[n] = (out[n] - 2*x0*(n+1)*out[n+1])/(2*n+1);
    }
  else
    {
      // Upwards recursion.
      double s = sqrt(x0);
      double t = sqrt(1+x0);
      out[0] = log(s+t)/s; //accurate because x0 is not too small. Unclear if asinh() is faster than this.
      if (ndeg<1)
	return;
      double hn = 1/t;
      out[1] = (hn-out[0])/(2*x0);
      for (int n=1;n<ndeg;n++)
	{
	  hn *= -(n-0.5)/(n*(1+x0));
	  out[n+1] = (hn - (2*n+1)*out[n])/(2*(n+1)*x0);
	}
    }
}

template<typename T>
void taylor_sqrtx_asinh_sqrtx(T *out, const T &x0, T tmp[2], int ndeg)
{
  taylor_asinh_sqrtx_div_sqrtx(out,x0,tmp,ndeg);
  for (int k=ndeg;k>0;k--)
    out[k] = x0*out[k] + out[k-1];
  out[0] = x0*out[0];
}

int main()
{
  double ref_009[] = {0.08870191426902674,
		      0.9717015552162798,
		      -0.1483848805546426,
		      0.06109901810484648,
		      -0.03330295625104985,
		      0.02076536927391453};

  int ndeg = 5;
  double c[ndeg+1];
  double tmp[2];
  double x = (double)(0.09);
  cout << "ndeg " << ndeg << endl;
  taylor_sqrtx_asinh_sqrtx(c,x,tmp,ndeg);
  cout.precision(17);
  for (int k=0;k<=ndeg;k++)
    cout << k << " " << c[k] << " " << c[k] - ref_009[k] << endl;

  return 0;
}
