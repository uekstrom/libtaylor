#include<cstdio>
#include "polyops.hpp"
#include "polymul.hpp"


int main()
{
  polymul::polynomial<double,3,3> a,aa,b;
  polymul::polynomial<double,1,3> c;
  int i;
  for (int i=0;i<polylen(1,3);i++)
    {
      c[i] = 1.1+i;
    }
  for (int i=0;i<polylen(3,3);i++)
    {
      a[i] = 0;
      b[i] = 2.13+i;
    }
  b[0] = 0;
  polymul::taylorcompose0<double,3,3,3>(a,b,c.c);
  taylor_compose(3,3,aa.c,b.c,c.c);
  //  taylormul_set(4,3,c.c,a.c,b.c);
  for (i=0;i<polylen(3,3);i++)
    {
      //      if (a[i] != aa[i])
      printf("%i %lf %lf %lf\n",i, a[i], aa[i],a[i] - aa[i]);
    }
  return 0;
} 

int main2()
{
  polymul::polynomial<double,3,2> resnew,resold,a,b;
  int i;
  for (int i=0;i<polylen(3,2);i++)
    {
      a[i] = 2.13+i;
      b[i] = 2*i+1.1;
      resnew[i] = 17;
      resold[i] = 17;
    }
  polymul_set_trunc_noconstb(3, resnew.c, 2, a.c, 2, b.c, 2);
  polymul::taylormul(resold,a,b);

  for (i=0;i<polylen(3,2);i++)
    {
      //      if (a[i] != aa[i])
      printf("%i %lf %lf %lf\n",i, resnew[i], resold[i], resnew[i] - resold[i]);
    }
  return 0;
} 

#if 0
int mainapa()
{
  polymul::polynomial<double,4,3> a,b,ct;
  polymul::polynomial<double,4,3> c,c2,aa;
  int i;
  for (int i=0;i<polylen(4,3);i++)
    {
      c[i] = i;
      c2[i] = i;
      aa[i] = 0;
    }
  for (int i=0;i<polylen(4,3);i++)
    {
      aa[i] = 1+i;
      a[i] = 1+i;
      b[i] = 2+i;
    }
  //polymul_add(4,c.c,a.c,3,b.c,3);
  a[0] = 0;
  polymul_add_trunc_noconstb(4,c.c,3,c.c,3,a.c,3);
  polymul::taylormul(c2,a); // c2 = c2*a
  
  for (i=0;i<polylen(4,3);i++)
    {
      if (c[i] != c2[i]-i)
	printf("%i %lf\n",i, c[i]-c2[i]-i);
    }
  return 0;
}

#endif
