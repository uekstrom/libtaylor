#ifndef TAYLOR_H
#define TAYLOR_H

#include <cmath>
#include "polymul.h"

template<class T, int Nvar, int Ndeg>
class taylor : public polynomial<T, Nvar, Ndeg>
{
public:
  taylor(void) {}
  taylor(T c0) : polynomial<T, Nvar,Ndeg>(c0) {} 
  // Set the constant term to c0 and the first order term of
  // variable var to var_value.
  taylor(T c0, int var, T var_value = 1) : polynomial<T, Nvar,Ndeg>(c0) 
  {
    assert(var>=0);
    assert(var<Nvar);
    if (Ndeg>0)
      this->c[var+1] = var_value;
  }
  template<class S>
  taylor<T,Nvar,Ndeg> &operator=(const S& c0)
  {
    polynomial<T,Nvar,Ndeg>::operator=(c0);
    return *this;
  }
  // Multiply each term with the correct factorials to get 
  // derivatives, i.e. x^2y^3 is multiplied by 2!3!
  void deriv_facs(void)
  {
    polydfac(*this);
  }
  void integrate(void)
  {
    assert(Nvar == 1 && "Implement this..");
    for (int i=Ndeg;i>0;i--)
      this->c[i] = this->c[i-1]/i;
    this->c[0] = 0.0;
  }
  // Set x := alpha*x, y = :alpha*y etc
  void stretch(T alpha)
  {
    T an = alpha;
    int k = 1;
    for (int i=1;i<=Ndeg;i++)
      {
	for (;k<polylen(Nvar,i);k++)
	  this->c[k] *= an;
	an *= alpha;
      }
  }
  void invert_parity(void)
  {
    this->stretch(-1); // TODO: optimized version of stretch(-1)
  }
  // Calculate a shifted version of this taylor expansion of one
  // variable.  out(x) ~= this(x+dx). Of course the shifted taylor
  // series is no longer exact around the new "zero" point. For this
  // reason NdegOut may be choosen to be smaller than Ndeg.
  // TODO: implement for more variables, and without runtime polylen.
  template<int NdegOut>
  void shift(taylor<T,1,NdegOut> &out, const T dx[Nvar]) const
  {
      assert(Nvar == 1);
      assert(NdegOut <= Ndeg);
      T dxpow[Ndeg+1];
      out = 0;
      dxpow[0] = 1;
      for (int i=1;i<=Ndeg;i++)
	  dxpow[i] = dx[0]*dxpow[i-1];
      for (int i=0;i<=NdegOut;i++)
	  for (int j=i;j<=Ndeg;j++)
	      out[i] += polylen(j-i,i)*dxpow[j-i]*(*this)[j];
  }
  taylor<T, Nvar,Ndeg> operator-(void) const
    {
      taylor<T, Nvar,Ndeg> res = *this;
      for (int i=0;i<res.size;i++)
	res[i] *= -1;
      return res;
    }
  void operator-=(const taylor<T, Nvar,Ndeg>& t)
    {
      for (int i=0;i<this->size;i++)
	this->c[i] -= t.c[i];
    }
  void operator+=(const taylor<T, Nvar,Ndeg>& t)
    {
      for (int i=0;i<this->size;i++)
	this->c[i] += t.c[i];
    }
  void operator*=(const T& scale)
    {
      for (int i=0;i<this->size;i++)
	this->c[i] *= scale;
    }
  template<int Ndeg2>
  void operator*=(const taylor<T, Nvar,Ndeg2>& t)
    {
      taylormul(*this,t);
    }
  /* Put sum_i coeff[i]*(this - this[0])^i in res,
     used when evaluating analytical functions of this */
  template<int Nres>
  void compose(taylor<T, Nvar,Nres>& res, const taylor<T,1,Nres> &coeff) const
    {
      assert(Nres >= Ndeg);
      taylor<T, Nvar,Ndeg> tmp = *this;
      tmp[0] = T(0);
      res = 0;
      for (int i=Nres;i>0;i--)
	{
	  res[0] += coeff[i];
	  res *= tmp;
	}
      res[0] += coeff[0];
    }
  T dot(const taylor<T, Nvar, Ndeg> &t) const
  {
    T sum = 0;
    for (int i=0;i<this->size;i++)
      sum += this->c[i]*t.c[i];
    return sum;
  }
};

// <> comparisons are taken to mean comparing the constant
// coefficient. This makes the transition from numbers to
// taylor objects easier.
template<class S, class T, int Nvar, int Ndeg>
bool operator<(const S &x, const taylor<T, Nvar,Ndeg> &t)
{
  return x < t[0];
}

template<class S, class T, int Nvar, int Ndeg>
bool operator>(const S &x, const taylor<T, Nvar,Ndeg> &t)
{
  return x > t[0];
}

template<class S, class T, int Nvar, int Ndeg>
bool operator<(const taylor<T, Nvar,Ndeg> &t, const S &x)
{
  return t[0] < x;
}

template<class S, class T, int Nvar, int Ndeg>
bool operator>(const taylor<T, Nvar,Ndeg> &t, const S &x)
{
  return t[0] > x;
}

template<class T, int Nvar, int Ndeg, class S>
taylor<T, Nvar, Ndeg> operator*(const S& x, const taylor<T, Nvar,Ndeg>& t)
{
  taylor<T, Nvar,Ndeg> tmp;
  for (int i=0;i<tmp.size;i++)
    tmp[i] = x*t[i];
  return tmp;
}

template<class T, int Nvar, int Ndeg, class S>
taylor<T, Nvar, Ndeg> operator*(const taylor<T, Nvar,Ndeg>& t, const S& x)
{
  taylor<T, Nvar,Ndeg> tmp;
  for (int i=0;i<tmp.size;i++)
    tmp[i] = x*t[i];
  return tmp;
}

template<class T, int Nvar, int Ndeg>
taylor<T, Nvar, Ndeg> operator*(const taylor<T, Nvar,Ndeg>& t1, 
				const taylor<T, Nvar,Ndeg>& t2)
{
  taylor<T, Nvar,Ndeg> tmp;
  taylormul(tmp,t1,t2);
  return tmp;
}

template<class T, int Nvar, int Ndeg, class S>
taylor<T, Nvar, Ndeg> operator+(const S& x, const taylor<T, Nvar,Ndeg>& t)
{
  taylor<T, Nvar,Ndeg> tmp = t;
  tmp[0] += x;
  return tmp;
}

template<class T, int Nvar, int Ndeg, class S>
taylor<T, Nvar, Ndeg> operator+(const taylor<T, Nvar,Ndeg>& t, const S& x)
{
  taylor<T, Nvar,Ndeg> tmp = t;
  tmp[0] += x;
  return tmp;
}

template<class T, int Nvar, int Ndeg>
taylor<T, Nvar, Ndeg> operator+(const taylor<T, Nvar,Ndeg>& t1, const taylor<T, Nvar,Ndeg>& t2)
{
  taylor<T, Nvar,Ndeg> tmp;
  for (int i=0;i<tmp.size;i++)
    tmp[i] = t1[i]+t2[i];
  return tmp;
}

template<class T, int Nvar, int Ndeg, class S>
taylor<T, Nvar, Ndeg> operator-(const S& x, const taylor<T, Nvar,Ndeg>& t)
{
  taylor<T, Nvar,Ndeg> tmp = -t;
  tmp[0] += x;
  return tmp;
}

template<class T, int Nvar, int Ndeg, class S>
taylor<T, Nvar, Ndeg> operator-(const taylor<T, Nvar,Ndeg>& t, const S &x)
{
  taylor<T, Nvar,Ndeg> tmp = t;
  tmp[0] -= x;
  return tmp;
}

template<class T, int Nvar, int Ndeg>
taylor<T, Nvar, Ndeg> operator-(const taylor<T, Nvar,Ndeg>& t1, 
				const taylor<T, Nvar,Ndeg>& t2)
{
  taylor<T, Nvar,Ndeg> tmp;
  for (int i=0;i<tmp.size;i++)
    tmp[i] = t1[i]-t2[i];
  return tmp;
}

#include "taylor_math.h"

#endif
