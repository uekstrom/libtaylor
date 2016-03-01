<TeXmacs|1.99.1>

<style|generic>

<\body>
  Long story short

  <\equation>
    f<around*|(|x|)>=<frac|asinh<sqrt|x>|<sqrt|x>>=<rsub|2>F<rsub|1><around*|(|<frac|1|2>,<frac|1|2>;<frac|3|2>,-x|)>=<sqrt|1+x>F<around*|(|1,1;<frac|3|2>;-x|)>
  </equation>

  (15.1.7 in A&S). Using 15.2.6

  <\equation>
    f<rsup|<around*|(|n|)>><around*|(|x|)>=<around*|(|-1|)><rsup|n><frac|<around*|(|1/2|)><rsub|n><rsup|2>|<around*|(|3/2|)><rsub|n>><around*|(|1+x|)><rsup|1/2-n>F<around*|(|1,1;<frac|3|2>+n;-x|)>
  </equation>

  The taylor coefficients are

  <\equation>
    F<around*|(|1,1;<frac|3|2>+n;-x|)>=<big|sum><rsub|k\<geqslant\>0><frac|<around*|(|-1|)><rsup|k>k!|<around*|(|3/2+n|)><rsub|k>>x<rsup|k>
  </equation>

  The prefactor satisfies

  <\equation>
    C<rsub|n>=<around*|(|-1|)><rsup|n><frac|<around*|(|1/2|)><rsub|n><rsup|2>|<around*|(|3/2|)><rsub|n>><around*|(|1+x|)><rsup|1/2-n>=-<frac|<around*|(|n-1/2|)><rsup|2>|<around*|(|n+1/2|)><around*|(|1+x|)>>C<rsub|n-1>
  </equation>

  \;

  For <math|0\<leqslant\>x\<leqslant\>1> we can compute the hypergeometric
  function at highest order <math|n> by accelerated Taylor series. Then apply

  <\equation>
    2x f<rsup|<around*|(|n+1|)>>=<around*|(|<frac|1|<sqrt|1+x>>|)><rsup|<around*|(|n|)>>-<around*|(|2n+1|)>f<rsup|<around*|(|n|)>><around*|(|x|)>
  </equation>

  <\equation>
    <around*|(|<frac|1|<sqrt|1+x>>|)><rsup|<around*|(|n|)>>=<around*|(|-1|)><rsup|n><around*|(|1/2|)><rsub|n><around*|(|1+x|)><rsup|-1/2-n>=h<rsup|<around*|(|n|)>>
  </equation>

  Verified ^^

  With the recurrence

  <\equation>
    h<rsup|<around*|(|n+1|)>>=-<frac|1|2>\<cdot\><frac|2n+1|1+x>h<rsup|<around*|(|n|)>>
  </equation>

  <\equation>
    f<rsup|<around*|(|n|)>><around*|(|x|)>=<frac|h<rsup|<around*|(|n|)>>-2x
    f<rsup|<around*|(|n+1|)>>|2n+1>
  </equation>

  Use <math|h<rsup|<around*|(|n|)>>> in the hypergeometric formula and
  simplifying pochhammer symbols gives

  <\equation>
    f<rsup|<around*|(|n|)>><around*|(|x|)>=-<frac|1|2>\<cdot\><frac|2n-1|2n+1>h<rsup|<around*|(|n-1|)>>F<around*|(|1,1;<frac|3|2>+n;-x|)>
  </equation>

  verified ^^

  Strategy: Calculate <math|h<rsup|<around*|(|n|)>>> up to <math|n-1> using
  upward recursion. Multiply highest element by hypergeometric function.
  Correct using downward recursion.

  \;

  \;

  For the asymptotic region we can use

  <\equation>
    F<around*|(|1,1;<frac|3|2>+n;-x|)>=<frac|<around*|(|1+x|)><rsup|-1>|\<Gamma\><around*|(|n+1/2|)>><big|sum><rsub|k\<geqslant\>0><frac|<around*|(|1/2+n|)><rsub|k>|k!><around*|(|1+x|)><rsup|-k>\<times\><around*|(|log<around*|(|1+x|)>+\<psi\><around*|(|k+1|)>-\<psi\><around*|(|1/2+n+k|)>|)>
  </equation>

  <\equation>
    \<psi\><around*|(|x+1|)>=\<psi\><around*|(|x|)>+<frac|1|x>;\<psi\><around*|(|1|)>=-\<gamma\>;\<psi\><around*|(|1/2|)>=-\<gamma\>-2log
    2
  </equation>

  <\equation>
    \<psi\><around*|(|1+k|)>-\<psi\><around*|(|1/2+n+k|)>=\<psi\><around*|(|1|)>-\<psi\><around*|(|1/2+n|)>+<big|sum><rsub|j=1><rsup|k><around*|(|<frac|1|j><rsup|>-<frac|1|n+j-1/2>|)>
  </equation>

  Defining <math|A<rsub|n>=><math|log<around*|(|1+x|)>>+<math|\<psi\><around*|(|1|)>-\<psi\><around*|(|1/2+n|)>>
  we have

  <\equation>
    F<around*|(|1,1;<frac|3|2>+n;-x|)>=<frac|<around*|(|1+x|)><rsup|-1>|\<Gamma\><around*|(|n+1/2|)>><big|sum><rsub|k=0><rsup|\<infty\>><frac|<around*|(|1/2+n|)><rsub|k>|k!><around*|(|1+x|)><rsup|-k>\<times\><around*|[|A<rsub|n>+<big|sum><rsub|j=1><rsup|k><around*|(|<frac|1|j><rsup|>-<frac|1|n+j-1/2>|)>|]>
  </equation>

  <\equation>
    f<rsup|<around*|(|n|)>><around*|(|x|)>=C<rsub|n><around*|(|1+x|)><rsup|-n-1/2><big|sum><rsub|k=0><rsup|\<infty\>><frac|<around*|(|1/2+n|)><rsub|k>|k!><around*|(|1+x|)><rsup|-k>\<times\><around*|[|A<rsub|n>+<big|sum><rsub|j=1><rsup|k><around*|(|<frac|1|j><rsup|>-<frac|1|j+n-1/2>|)>|]>
  </equation>

  \;

  <\equation>
    <big|sum><rsub|k=0><rsup|\<infty\>><frac|<around*|(|1/2+n|)><rsub|k>|k!><around*|(|1+x|)><rsup|-k>=<around*|(|<frac|x+1|x>|)><rsup|n+1/2><rsup|>\<Rightarrow\>
  </equation>

  <\equation>
    f<rsup|<around*|(|n|)>><around*|(|x|)>=C<rsub|n><frac|log<around*|(|1+x|)>+\<psi\><around*|(|1|)>-\<psi\><around*|(|1/2+n|)>|x<rsup|n+1/2>>+<big|sum>\<ldots\>
  </equation>

  \;

  We can combine this with

  <\equation>
    F<around*|(|1,1;<frac|3|2>;-x|)>=F<around*|(|1/2,1/2;3/2;-4x<around*|(|1+x|)>|)>
  </equation>

  \;

  (verified^^)

  Follows that

  <\equation>
    asinh<around*|(|<sqrt|x>|)>=<frac|asinh<around*|(|<sqrt|4x<around*|(|1+x|)>>|)>|2>=
  </equation>

  <\equation>
    <frac|1|2><around*|[|log<around*|(|2<sqrt|4x<around*|(|1+x|)>>|)>+<big|sum><rsub|n=1><rsup|\<infty\>><around*|(|-1|)><rsup|n-1><frac|<around*|(|2n-1|)>!!|2n<around*|(|2n|)>!!>\<cdot\><frac|1|<around*|(|4
    x<around*|(|1+x|)>|)><rsup|n>>|]>
  </equation>

  \;

  <\equation>
    M<around*|(|x|)>=x<around*|(|1+x|)>
  </equation>

  <\equation>
    log M<around*|(|x|)>=log x+log<around*|(|1+x|)>
  </equation>

  <\equation>
    log M<around*|(|M<around*|(|x|)>|)>=log M<around*|(|1+M|)>=log M+log
    <around*|(|1+M|)>
  </equation>

  \;

  \;

  \;

  And asymptotic expansion (should be quadratically convergent?)

  <\equation>
    -4x<around*|(|1+x|)>=-z\<Rightarrow\>x=<frac|<sqrt|1+z>-1|2>
  </equation>

  \;

  Old stuff:

  By 15.3.22 in A&S we have

  <\equation>
    F<rsub|><around*|(|<frac|1|2>,<frac|1|2>;<frac|3|2>,-x|)>=F<around*|(|1,1;<frac|3|2>;<frac|1|2>-<frac|1|2><sqrt|1+x>|)>=
  </equation>

  <\equation>
    2<around*|(|1+<sqrt|1+x>|)><rsup|-1>F<around*|(|1,<frac|1|2>;<frac|3|2>;<frac|<sqrt|1+x>-1|<sqrt|1+x>+1>|)>
  </equation>

  But also

  <\equation>
    x F<rsub|><around*|(|<frac|1|2>,<frac|1|2>;<frac|3|2>,-x|)>+<frac|3|2><sqrt|1+x>-<frac|3|2><around*|(|1+x|)>F<around*|(|<frac|1|2>,<frac|1|2>;<frac|1|2>;-x|)>=0
  </equation>

  \;

  We also have\ 

  <\equation>
    F<around*|(|1,<frac|1|2>;<frac|3|2>;x|)>=<frac|atanh<around*|(|<sqrt|x>|)>|<sqrt|x>>
  </equation>

  --

  \;

  <\equation>
    f<rsup|<around*|(|n|)>>=<frac|h<rsup|<around*|(|n|)>>-2x
    f<rsup|<around*|(|n+1|)>>|2n+1>
  </equation>

  <\equation>
    <around*|[|h<rsup|<around*|(|n|)>>-<around*|(|2n+1|)>f<rsup|<around*|(|n|)>>|]>=2x
    f<rsup|<around*|(|n+1|)>>
  </equation>

  <\equation>
    F<rsup|<around*|(|n|)>>=x<rsup|n>f<rsup|<around*|(|n|)>>
  </equation>

  <\equation>
    h<rsup|<around*|(|n|)>>-x<rsup|-n><around*|(|2n+1|)>F<rsup|<around*|(|n|)>>=2x
    x<rsup|-n-1>F<rsup|<around*|(|n+1|)>>\<Leftrightarrow\>
  </equation>

  <\equation>
    <around*|(|x<rsup|n>h<rsup|<around*|(|n|)>>-<around*|(|2n+1|)>F<rsup|<around*|(|n|)>>|)>/2=
    F<rsup|<around*|(|n+1|)>>
  </equation>

  \;

  <\equation>
    x<rsup|-n>F<rsup|<around*|(|n|)>>=<frac|h<rsup|<around*|(|n|)>>-2x
    x<rsup|-n-1> F<rsup|<around*|(|n+1|)>>|2n+1>
  </equation>

  <\equation>
    x<rsup|n>h<rsup|<around*|(|n|)>>-<around*|(|2n+1|)>F<rsup|<around*|(|n|)>>=
    F<rsup|<around*|(|n+1|)>>
  </equation>

  \;

  <\equation>
    h<rsup|<around*|(|n+1|)>>=-<around*|(|1+x|)><rsup|-1><around*|(|n+1/2|)>h<rsup|<around*|(|n|)>>
  </equation>

  <\equation>
    x<rsup|n+1>h<rsup|<around*|(|n+1|)>>=-<frac|x|1+x><around*|(|n+1/2|)><around*|[|x<rsup|n>h<rsup|<around*|(|n|)>>|]>
  </equation>

  \;

  <\equation>
    f<rsup|<around*|(|n|)>>=<frac|h<rsup|<around*|(|n|)>>|2n+1>-<frac|2x|<around*|(|2n+1|)><around*|(|2n+3|)>><around*|(|h<rsup|<around*|(|n+1|)>>-2x
    f<rsup|<around*|(|n+2|)>>|)>=
  </equation>

  <\equation>
    <frac|h<rsup|<around*|(|n|)>>|2n+1>+<frac|2x|<around*|(|2n+1|)><around*|(|2n+3|)>><around*|(|<around*|(|n+1/2|)><frac|h<rsup|<around*|(|n|)>>|1+x>+2x
    f<rsup|<around*|(|n+2|)>>|)>=
  </equation>

  <\equation>
    <frac|h<rsup|<around*|(|n|)>>|2n+1>+<frac|x|2n+3>\<cdot\><frac|h<rsup|<around*|(|n|)>>|1+x>+<frac|4x<rsup|2>|<around*|(|2n+1|)><around*|(|2n+3|)>>
    f<rsup|<around*|(|n+2|)>>=
  </equation>

  <\equation>
    h<rsup|<around*|(|n|)>><around*|(|x|)><big|sum><rsub|k\<geqslant\>0><rsup|M-1><frac|1|2n+2k+1><around*|(|<frac|x|1+x>|)><rsup|k>+<frac|<around*|(|2<rsup|>x|)><rsup|M>|<around*|(|2n+1|)><around*|(|2n+3|)>\<cdots\><around*|(|2n+2M+1|)>>f<rsup|<around*|(|n+M|)>>=f<rsup|<around*|(|n|)>>
  </equation>

  --

  <\equation>
    <frac|<around*|(|a|)><rsub|k>|<around*|(|a+1|)><rsub|k>>=<frac|a<around*|(|a+1|)>\<cdots\><around*|(|a+k-1|)>|<around*|(|a+1|)>\<cdots\><around*|(|a+2|)>\<cdots\><around*|(|a+k|)>>=<frac|a|a+k>
  </equation>

  <\equation>
    <frac|1|2n+2k+1>=<frac|1|2>\<cdot\><frac|1|n+1/2+k>=<frac|1|2n+1>\<cdot\><frac|<around*|(|n+1/2|)><rsub|k>|<around*|(|n+3/2|)><rsub|k>>
  </equation>

  <\equation>
    <around*|(|1|)><rsub|k>=k!
  </equation>

  --

  <\equation>
    h<rsup|<around*|(|n|)>><around*|(|x|)><big|sum><rsub|k\<geqslant\>0><rsup|\<infty\>><frac|1|2n+2k+1><around*|(|<frac|x|1+x>|)><rsup|k>=
    <frac|h<rsup|<around*|(|n|)>><rsub|><around*|(|x|)>|2n+1>
    <rsub|2>F<rsub|1><around*|(|1,n+1/2;n+3/2;<frac|x|1+x>|)>=f<rsup|<around*|(|n|)>><around*|(|x|)>
  </equation>

  <\equation>
    b-c=-1
  </equation>

  <\equation>
    \ <rsub|2>F<rsub|1><around*|(|1,n+1/2;n+3/2;<frac|x|1+x>|)>=
    <around*|(|1+x|)><rsub|2>F<rsub|1><around*|(|1,1;n+3/2;-x|)>
  </equation>

  Verified^^

  \;

  <\equation>
    <rsub|2>F<rsub|1><around*|(|1,1;n+3/2;-x|)>=<big|sum><rsub|k\<geqslant\>0><frac|<around*|(|1|)><rsub|k><around*|(|1|)><rsub|k>|<around*|(|n+3/2|)><rsub|k>k!><around*|(|-x|)><rsup|k>
  </equation>

  \;

  Using the Euler transformation one can even write

  <\equation>
    <around*|(|1+x|)><rsup|-n+1/2>F<rsub|><around*|(|1,1;n+3/2;-x|)>=F<around*|(|n+1/2,n+1/2,n+3/2,-x|)>
  </equation>

  and from there write <math|f<rsup|<around*|(|n|)>>> as a scaled
  hypergeometric function (related to the deriative of
  <math|F<around*|(|1/2,1/2,3/2|)>>).

  \;

  \;

  The rhs expression seems slightly favourable numerically. Sumalt
  acceleration is not successful. Odd/even terms must be handled separately.\ 

  \;

  --

  <\equation>
    F<rsub|><around*|(|1,1;n+3/2;-x|)>=<frac|<around*|(|1+x|)><rsup|-1>|\<Gamma\><around*|(|1|)>\<Gamma\><around*|(|n+1/2|)>><big|sum><rsub|j=0><rsup|\<infty\>><frac|<around*|(|1|)><rsub|j><around*|(|n+1/2|)><rsub|j>|j!j!><around*|(|1+x|)><rsup|-j>\<times\>
  </equation>

  <\equation>
    <around*|[|log<around*|(|1+x|)>+\<psi\><around*|(|j+1|)>-\<psi\><around*|(|n+1/2-j|)>|]>
  </equation>

  --

  <\equation>
    <sqrt|1+x><sqrt|1-x>=<sqrt|1-x<rsup|2>>
  </equation>

  --

  <\equation>
    \ <rsub|>F<rsub|><around*|(|1,n+1/2;n+3/2;<frac|x|1+x>|)>=C<rsub|n><big|int><rsub|0><rsup|1>t<rsup|n-1/2><around*|(|1-t<frac|x|1+x>|)><rsup|-1>d
    t
  </equation>

  <\equation>
    F<around*|(|1,1,n+3/2,-z|)>=A<rsub|n><big|int><rsub|0><rsup|1><around*|(|1-t|)><rsup|n-1/2><around*|(|1+t
    z|)><rsup|-1>d t
  </equation>

  <\equation>
    t=1-u<rsup|2>
  </equation>

  <\equation>
    dt=-2u d u
  </equation>

  <\equation>
    I=2<big|int><rsub|0><rsup|1>u<rsup|2n><around*|(|1+<around*|(|1-u<rsup|2>|)>z|)><rsup|-1>d
    u
  </equation>

  <\equation>
    <around*|(|1-u<rsup|2>|)>z=s
  </equation>

  <\equation>
    u=<sqrt|1-<frac|s|z>>
  </equation>

  <\equation>
    u<rsup|2>=1-<frac|s|z>
  </equation>

  <\equation>
    du =-<frac|1|2z<sqrt|1-<frac|s|z>>>
  </equation>

  <\equation>
    I=<frac|1|z><big|int><rsub|0><rsup|z><around*|(|1-<frac|s|z>|)><rsup|n-1/2><frac|1|1+s>
    ds
  </equation>

  --

  <\equation>
    -log<around*|(|c<rsub|k>z<rsup|k>|)>=-log<around*|(|c<rsub|k>|)>-k
    log<around*|(|z|)>
  </equation>
</body>