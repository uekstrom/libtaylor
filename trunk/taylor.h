#ifndef TAYLOR_H
#define TAYLOR_H

#include <cmath>
#include "polymul.h"

template<class T, int Nvar, int Ndeg>
class taylor : public polynomial<T, Nvar, Ndeg>
{
public:
  typedef taylor<T,Nvar,Ndeg> type;
  using polynomial<T, Nvar, Ndeg>::size;
  using polynomial<T, Nvar, Ndeg>::c;

  taylor(void) {}
  // This does not work when T is a int. One ugly solution is to
  // specialize taylor for int T's, but perhaps not so useful?
  template<typename S>
  taylor(const S &c0) : polynomial<T, Nvar,Ndeg>(T(c0)) {} 
  //  taylor(double c0) : polynomial<T, Nvar,Ndeg>(T(c0)) {} 
  //  taylor(const T &c0) : polynomial<T, Nvar,Ndeg>(c0) {} 
  // Set the constant term to c0 and the first order term of
  // variable var to var_value.
  template<typename S>
  taylor(const S &c0, int var) : polynomial<T, Nvar,Ndeg>(T(c0)) 
  {
    assert(var>=0);
    assert(var<Nvar);
    if (Ndeg>0)
      c[var+1] = 1;
  }
  template<typename S, typename U>
  taylor(const S &c0, int var, const U &var_value) : polynomial<T, Nvar,Ndeg>(T(c0)) 
  {
    assert(var>=0);
    assert(var<Nvar);
    if (Ndeg>0)
      c[var+1] = var_value;
  }
  template<typename S>
  taylor<T,Nvar,Ndeg> &operator=(const S& c0)
  {
    polynomial<T,Nvar,Ndeg>::operator=(T(c0));
    return *this;
  }
  template<int N>
  taylor<T,Nvar-1,N> &pick_order(void)
  {
    assert(N<=Ndeg);
    return *reinterpret_cast<taylor<T,Nvar-1,N> *>
      (c+polymul_internal::polylen<Nvar,N-1>::len);
  }
  template<int N>
  const taylor<T,Nvar-1,N> &pick_order(void) const
  {
    assert(N<=Ndeg);
    return *reinterpret_cast<const taylor<T,Nvar-1,N> *>
      (c+polymul_internal::polylen<Nvar,N-1>::len);
  }
  template<int N>
  void pick_order(taylor<T,Nvar-1,N> &p) const
  {
    p = *reinterpret_cast<const taylor<T,Nvar-1,N> *>
      (c+polymul_internal::polylen<Nvar,N-1>::len);
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
      c[i] = c[i-1]/i;
    c[0] = 0;
  }
  // Set x := alpha*x, y = :alpha*y etc
  void stretch(T alpha)
  {
    T an = alpha;
    int k = 1;
    for (int i=1;i<=Ndeg;i++)
      {
	for (;k<polylen(Nvar,i);k++)
	  c[k] *= an;
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
      for (int i=0;i<size;i++)
	c[i] -= t.c[i];
    }
  void operator+=(const taylor<T, Nvar,Ndeg>& t)
    {
      for (int i=0;i<size;i++)
	c[i] += t.c[i];
    }
  template<class S>
  void operator-=(const S& x)
    {
      c[0] -= x;
    }
  template<class S>
  void operator+=(const S& x)
    {
      c[0] += x;
    }
  template<class S>
  void operator*=(const S& scale)
    {
      for (int i=0;i<size;i++)
	c[i] *= scale;
    }
  template<class S>
  void operator/=(const S& scale)
    {
      for (int i=0;i<size;i++)
	c[i] /= scale;
    }
  void operator/=(const taylor<T,Nvar,Ndeg>& t)
    {
      taylor<T,Nvar,Ndeg> tinv = 1/t;
      *this*=tinv;
    }
  void operator*=(int scale)
    {
      for (int i=0;i<size;i++)
	c[i] *= scale;
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
      tmp[0] = 0;
      taylorcompose0(res,tmp,coeff.c);
    }
  T dot(const taylor<T, Nvar, Ndeg> &t) const
  {
    T sum = 0;
    for (int i=0;i<size;i++)
      sum += c[i]*t.c[i];
    return sum;
  }
  const T &operator[](int i) const
  {
    assert(i>=0);
    assert(i<this->size);
    return c[i];
  }
  T &operator[](int i)
  {
    assert(i>=0);
    assert(i<this->size);
    return c[i];
  }
};

// Define a taylor polynomial with taylor polynomial coefficients,
// i.e. a tensoring of polynomial spaces.
template<int Ndim, class T, int Nvar, int Ndeg>
struct tensored_taylor
{
  typedef taylor<typename 
    tensored_taylor<Ndim-1,T,Nvar,Ndeg>::type,Nvar,Ndeg> type;
};

template<class T, int Nvar, int Ndeg>
struct tensored_taylor<1,T,Nvar,Ndeg>
{
  typedef taylor<T,Nvar,Ndeg> type;
};

template<class T, int Nvar, int Ndeg>
void as_taylor(taylor<T,Nvar,Ndeg> *&ptr, T *data)
{
  ptr = reinterpret_cast<taylor<T,Nvar,Ndeg> *>(data);
}

template<class T, int Nvar, int Ndeg>
void as_taylor(const taylor<T,Nvar,Ndeg> *&ptr, const T *data)
{
  ptr = reinterpret_cast<const taylor<T,Nvar,Ndeg> *>(data);
}

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

template<class S, class T, int Nvar, int Ndeg>
bool operator!=(const taylor<T, Nvar,Ndeg> &t, const S &x)
{
  return t[0] != x;
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

#ifndef TAYLOR_NOSTDCXX

#include <iostream>

template<class num, int Nvar, int Ndeg>
std::ostream &operator<<(std::ostream& stream, const taylor<num,Nvar,Ndeg> &t)
{
  stream << "{" << t[0];
  for (int i=1;i<t.size;i++)
    stream << ", " << t[i];
  stream << "}";
  return stream;
} 

#endif

#endif
