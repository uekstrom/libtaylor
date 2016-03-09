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

#include <stdio.h>

namespace polymul{
#include "polymul.h"
};

int main()
{
  polymul::polynomial<double,4,3> a,b,ct;
  polymul::polynomial<double,4,6> c,c2,aa;
  int i;
  for (int i=0;i<polylen(4,6);i++)
    {
      c[i] = 0;
      c2[i] = 0;
      aa[i] = 0;
    }
  for (int i=0;i<polylen(4,3);i++)
    {
      aa[i] = 1+i;
      a[i] = 1+i;
      b[i] = 2+i;
    }
  polymul::polymul(c2,a,b);
  polymul::taylormul(ct,a,b);
  taylormul_set(4,3,c.c,a.c,b.c);
  for (i=0;i<polylen(4,3);i++)
    {
      if (c[i] != ct[i])
	printf("%i %lf\n",i, c[i]-ct[i]);
    }
  return 0;
}
