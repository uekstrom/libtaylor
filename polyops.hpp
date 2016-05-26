// Copyright (c) Ulf Ekstrom 2016 uekstrom@gmail.com

#ifndef POLYOPS_HPP
#define POLYOPS_HPP
#include <cassert>

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

/* nchoosek(nvar+ndeg,ndeg) */
int polylen(int nvar, int ndeg)
{
  int i, val = 1;
  if (ndeg < 0) // this simplifies some recursions
    return 0;
  for (i=1;i<=ndeg;i++)
    val = (val*(nvar+i))/i;
  return val;
}

// TODO: It is possible to optimize the use of polylen below. This is only worth it if
// one does not vectorize the coefficients. In that case it may be better to unroll the
// functions completely.


// For all of the functions below the c[] and a[] array may alias
// and the outcome is still as expected. The b[] array may not alias c[].


// Multiply polynomials with coefficients a and b, and add the result to c.
template<typename T>
void polymul_add(int nvar,
		 T c[],
		 const T a[], int adeg,
		 const T b[], int bdeg)
{
  int i,j;
  if (nvar == 0)
    {
      c[0] += a[0]*b[0];
    }
  else
    {
      for (i=adeg+bdeg;i>=0;i--)
	{
	  int u = MIN(bdeg,i);
	  int l = MAX(0,i-adeg);
	  for (j=l;j<=u;j++)
	    polymul_add(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-j-1), i-j, b + polylen(nvar,j-1), j);
	}
    }
}

// Multiply polynomials with coefficients a and b, and add the result to c.
// Output is truncated at order cdeg
template<typename T>
void polymul_add_trunc(int nvar,
		       T c[], int cdeg,
		 const T a[], int adeg,
		 const T b[], int bdeg)
{
  int i,j;
  assert(cdeg <= adeg+bdeg);
  if (nvar == 0)
    {
      c[0] += a[0]*b[0];
    }
  else
    {
      for (i=cdeg;i>=0;i--)
	{
	  int u = MIN(bdeg,i);
	  int l = MAX(0,i-adeg);
	  for (j=l;j<=u;j++)
	    polymul_add(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-j-1), i-j, b + polylen(nvar,j-1), j);
	}
    }
}

// Multiply polynomials with coefficients a and b, setting c to the result.
template<typename T>
void polymul_set(int nvar,
		   T c[],
		   const T a[], int adeg,
		   const T b[], int bdeg)
{
  int i,j;
  if (nvar == 0)
    {
      c[0] = a[0]*b[0];
    }
  else
    {
      for (i=adeg+bdeg;i>=0;i--)
	{
	  int u = MIN(bdeg,i);
	  int l = MAX(0,i-adeg);
	  assert(l<=u);
	  polymul_set(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-l-1), i-l, b + polylen(nvar,l-1), l);
	  for (j=l+1;j<=u;j++)
	    polymul_add(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-j-1), i-j, b + polylen(nvar,j-1), j);
	}
    }
}

// Multiply polynomials with coefficients a and b, setting c to the result. Only
// coefficients up to degree cdeg are affected.
template<typename T>
void polymul_set_trunc(int nvar,
		       T c[], int cdeg,
		       const T a[], int adeg,
		       const T b[], int bdeg)
{
  int i,j;
  assert(cdeg <= adeg+bdeg);
  if (nvar == 0)
    {
      c[0] = a[0]*b[0];
    }
  else
    {
      for (i=cdeg;i>=0;i--)
	{
	  int u = MIN(bdeg,i);
	  int l = MAX(0,i-adeg);
	  assert(l<=u);
	  polymul_set  (nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-l-1), i-l, b + polylen(nvar,l-1), l);
	  for (j=l+1;j<=u;j++)
	    polymul_add(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-j-1), i-j, b + polylen(nvar,j-1), j);
	}
    }
}

// Multiply polynomials with coefficients a and b, setting c to the result. Only
// coefficients up to degree cdeg are affected.
// b[0] = 0 is assumed, and c[0] is undefined at exit
template<typename T>
void polymul_set_trunc_noconstb(int nvar,
				T c[], int cdeg,
				const T a[], int adeg,
				const T b[], int bdeg)
{
  int i,j;
  assert(cdeg <= adeg+bdeg);
  if (nvar == 0)
    {
      
    }
  else
    {
      for (i=cdeg;i>0;i--)
	{
	  int u = MIN(bdeg,i);
	  int l = MAX(1,i-adeg);
	  if (l<=u)
	    polymul_set  (nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-l-1), i-l, b + polylen(nvar,l-1), l);
	  for (j=l+1;j<=u;j++)
	    polymul_add(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-j-1), i-j, b + polylen(nvar,j-1), j);
	}
    }
}


/*
  Compute y = sum_i=0,ndeg c[i]*(x-x0)^i = c[0] + (x-x0)*(c[1] + (x-x0)*(..))
 */
template<typename T, typename S>
void taylor_compose(int nvar, int ndeg,
		   T y[], const T x[], const S c[])
{
  y[0] = c[ndeg];
  for (int i=ndeg-1;i>=0;i--)
    {
      polymul_set_trunc_noconstb(nvar,y,ndeg-i,y,ndeg-i-1,x,ndeg-i);
      y[0] = c[i];
    }
}


/* adeg = bdeg = ndeg, and truncate output at ndeg. */
template<typename T>
void taylormul_add(int nvar, int ndeg,
		   T c[],
		   const T a[],
		   const T b[])
{
  int i,j;
  if (nvar == 0)
    {
      c[0] += a[0]*b[0];
    }
  else
    {
      for (i=ndeg;i>=0;i--)
	{
	  for (j=0;j<=i;j++)
	    polymul_add(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-j-1), i-j, b + polylen(nvar,j-1), j);
	}
    }
}

template<typename T>
void taylormul_set(int nvar, int ndeg,
		   T c[],
		   const T a[],
		   const T b[])
{
  int i,j;
  if (nvar == 0)
    {
      c[0] = a[0]*b[0];
    }
  else
    {
      for (i=ndeg;i>=0;i--)
	{
	  polymul_set(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-1), i, b, 0);
	  for (j=1;j<=i;j++)
	    polymul_add(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-j-1), i-j, b + polylen(nvar,j-1), j);
	}
    }
}



#endif
