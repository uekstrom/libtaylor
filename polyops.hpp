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

// For all of the functions below the c[] and a[] array may alias
// and the outcome is still as expected.

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
	  polymul_set(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-l-1), i-l, b + polylen(nvar,l-1), l);
	  for (j=l+1;j<=u;j++)
	    polymul_add(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-j-1), i-j, b + polylen(nvar,j-1), j);
	}
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



// c += a*b with b[0] implicitly set to 0.
template<typename T>
void polymul_add_trunc_noconstb(int nvar,
				T c[], int cdeg,
				const T a[], int adeg,
				const T b[], int bdeg)
{
  int i,j;
  assert(cdeg <= adeg+bdeg);
  if (nvar == 0)
    {
      return;
    }
  else
    {
      for (i=cdeg;i>=0;i--)
	{
	  int u = MIN(bdeg,i);
	  int l = MAX(0,i-adeg);
	  if (l == 0)
	    polymul_add_trunc_noconstb(nvar-1, c + polylen(nvar,i-1), i, a + polylen(nvar,i-1), i, b, 0);
	  else
	    polymul_add(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-l-1), i-l, b + polylen(nvar,l-1), l);
	  for (j=l+1;j<=u;j++)
	    polymul_add(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-j-1), i-j, b + polylen(nvar,j-1), j);
	}
    }
}

// c = a*b with b[0] implicitly set to 0.
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
      c[0] = 0;
    }
  else
    {
      for (i=cdeg;i>=0;i--)
	{
	  int u = MIN(bdeg,i);
	  int l = MAX(0,i-adeg);
	  if (l == 0)
	    polymul_set_trunc_noconstb(nvar-1, c + polylen(nvar,i-1), i, a + polylen(nvar,i-1), i, b, 0);
	  else
	    polymul_set(nvar-1, c + polylen(nvar,i-1), a + polylen(nvar,i-l-1), i-l, b + polylen(nvar,l-1), l);
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
  int plen = polylen(nvar,ndeg);
  y[0] = c[ndeg];
  for (int i=1;i<plen;i++)
    y[i] = 0;
  for (int i=ndeg-1;i>=0;i--)
    {
      polymul_set_trunc_noconstb(nvar,y,ndeg-i,y,ndeg-i-1,x,ndeg-i);
      y[0] += c[i];
    }
}


#endif
