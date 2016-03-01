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
  if (ndeg>0)
    out[1] = out[0]+1; // TODO: would it help to use exp() here?
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


// At x0=0 the Taylor expansion multiplied by 2*k-1 (term wise) is equal to the taylor exp. of
// 2x/sqrt(1+x). 
template<typename T>
void taylor_sqrtx_asinh_sqrtx(T *out, const T &x0, T tmp[2], int ndeg)
{
  tmp[0] = sqrt(x0);
  out[0] = asinh(tmp[0]);
  out[0] *= tmp[0];
  if (ndeg < 1)
    return;
  out[1] = 2*p*x0*out[0];
  if (ndeg < 2)
    return;

}


int main()
{
  int ndeg = 5;
  long double c[ndeg+1];
  long double tmp[2];
  taylor_atanh(c,(long double)1e-6,tmp,ndeg);
  cout.precision(17);
  for (int k=0;k<=ndeg;k++)
    cout << k << " " << c[k] << endl;
  return 0;
}
