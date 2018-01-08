# JGMM

We present a derivation of Gaussian [Mixture Model](https://en.wikipedia.org/wiki/Mixture_model) in the [J programming language](http://jsoftware.com/).
The derivation will follow the work of [Leonid Perlovsky](https://en.wikipedia.org/wiki/Leonid_Perlovsky) on [Dynamic Logic](http://www.springer.com/fr/book/9783642228292).

## Some basic utilities

```
load 'plot trig numeric'

NB.utils......................................................................

NB. if y is an empty array, return an empty list, otherwise, return y
as_z_=: ]`(($0)"_)@.(0:e.$)
B1_z_=: <"1               NB. rank-1 box
CP_z_=:{@(,&<)            NB. cartesian product
notin_z_=: [ {~ [: I.@:-. [ e."0 1 ] NB. elements of x not in y

id_z_=: =@i.              NB. identity matrix of size y
MP_z_=: +/ . *            NB. matrix product
det_z_=: -/ . *           NB. determinant
QF_z_=: ] MP"1 [ MP"2 1 ] NB. quadratic form
```

## Multivariate normal distribution

We provide some code to draw data from a multivariate normal distribution.
From the uniform random number generator of `J`, we build a 1-dimensional normal random number generator. The verb `randn` comes from the addon `math/mt/rand.ijs`.

We follow a ["widely used method"](https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution) (Wikipedia dixit), to draw random values from a multivariate normal distribution (see verb `randmultin`). The J/Lapack interface is used to compute the Cholesky decomposition of the covariance matrix.

Finally, we compute the coordinates of the 2-sigma concentration ellipse for a 2d normal distribution (see verb `cellipse`). We only compute the points for the positive quadrant and we then proceed by symmetry.

An ellipse is a set of points in a plane whose distances from two fixed points (viz. foci) add to a constant. We first consider the equation of an ellipse in standard form (i.e., centred on the origin and aligned with the axis of the orthonormal basis).

![ellipse](media/ellipse-1.png)


$$foci: (\pm c,0) \;with\; c>0$$
\\[
\begin{equation*}
\begin{split}
\sqrt{(x+c)^2+y^2}+\sqrt{(x-c)^2+y^2} & = 2a \\
\frac{x^2}{a^2}+\frac{y^2}{b^2} & = 1 \;with\; b=\sqrt{a^2-c^2} \\
y & = \sqrt{\left( 1 - \frac{x^2}{a^2} \right)b^2}
\end{split}
\end{equation*}
\\]

The verb `el` implements this last equation to compute the ordinate of an ellipse in standard form.

The standard form equation of an ellipse can be written in matrix notation:
\\[1=\frac{x^2}{a^2}+\frac{y^2}{b^2} = 
\begin{bmatrix}x & y\end{bmatrix}
\begin{bmatrix}
\frac{1}{a^2} & 0 \\
0 & \frac{1}{b^2}
\end{bmatrix}
\begin{bmatrix}x \\ y\end{bmatrix} =
X^T \Lambda \, X\\]

To model a generic 2d ellipse, we can start from a standard form one to which we apply a linear transformation combining translation, scaling and rotation. For convenience of notation, we will note the scaling factor $\sqrt{k}$, the translation vector $M$ and the rotation matrix $R$:

\\[
\begin{equation*}
\begin{split}
\tilde{X} & = M + \sqrt{k} RX \\
X & = \sqrt{k}^{-1} R^T (\tilde{X}-M)
\end{split}
\end{equation*}
\\]

Thus, applying the linear transformation to the equation of the ellipse, we have:

\\[
\begin{equation*}
\begin{split}
X^T \Lambda \, X & = 1 \\
\sqrt{k}^{-2} (\tilde{X}-M)^T R \Lambda R^T (\tilde{X}-M) & = 1 \\
(\tilde{X}-M)^T R \Lambda R^T (\tilde{X}-M) & = k \;(Eq.1)\\
\end{split}
\end{equation*}
\\]

Let us introduce the probability density function of a 2d normal distribution with covariance $C$ and mean $M$. We use the notation $D$ (like "Difference") for $X-M$.
$$G(X) = \frac{1}{2\pi \sqrt{det C}} \; exp \left[ -\frac{1}{2} D^T C^{-1} D \right]$$

We recognise in the equation for the set of points of equal density the one of an ellipse:

\\[
\begin{equation*}
\begin{split}
& \;\;\;\;\; \left\{ X: G(X) = k_1 \right\} \\
&= \left\{ X: D^T C^{-1} D = k \right\} \; with \; k=-2ln\left(2\pi k_1 \sqrt{det C}\right) \;(Eq.2)
\end{split}
\end{equation*}
\\]

Matching $(Eq.1)$ and $(Eq.2)$, we have: $C^{-1} \leftrightarrow R \Lambda R^T$, or equivalently: $C \leftrightarrow R \Lambda^{-1} R^T \; (Eq.3)$ (given that the inverse of an orthogonal matrix, here the rotation matrix, is its transpose).
Therefore, to draw the concentration ellipse, we compute the eigen decomposition of the inverse of the covariance matrix. Then, we compute the coordinates of the standard form ellipse ($\Lambda$ gives $1/a^2$ and $1/b^2$) and we transform it with a scaling factor $\sqrt{k}$, a rotation $R$ and a translation $M$.

Given a vector $v$, the projection of the data $X$ on $v$ is $v^T X$. The variance of the projected data is $v^T C v$ (see this very nice blog post about a [geometric interpretation of the covariance matrix](http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/)). The $v$ maximising the covariance of the projected data is the largest eigenvector of $C$ (see this [straightforward derivation](https://en.wikipedia.org/wiki/Rayleigh_quotient#Formulation_using_Lagrange_multipliers) based on the method of Lagrange multipliers).

Therefore, citing the [aforementioned article](http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/): 
> [...]the largest eigenvector of the covariance matrix always points into the direction of the largest variance of the data, and the magnitude of this vector equals the corresponding eigenvalue. The second largest eigenvector is always orthogonal to the largest eigenvector, and points into the direction of the second largest spread of the data.

From $(Eq.3)$, we have that $a^2$ is the magnitude of the vector that points into the direction of the largest variance which is also, by definition, $\sigma_x^2$, the square of the standard deviation. Thus, we understand why we are speaking of the $\sqrt{k}\sigma$-concentration ellipse. With the code below (see verb `cellipse`), we draw the $2\sigma$-concentration ellipse.

```
NB.RandN......................................................................

NB. draw values from a multivariate normal distribution
NB. with mean vector M and covariance matrix C
NB. find coordinates of the 2-sigma concentration ellipse
coclass 'RandN'

randu=: (?@$&0) :((p.~ -~/\)~ $:) NB. Uniform distribution U(a,b) with support (a,b)
rande=: -@^.@randu : (* $:) NB. Exponential distribution E(μ) with mean μ
randn=: (($,) +.@:(%:@(2&rande) r. 0 2p1&randu)@>.@-:@(*/)) : (p. $:) NB. Normal distribution N(μ,σ^2) of real numbers

require 'math/lapack'
require 'math/lapack/potrf'
chol=: potrf_jlapack_
require 'math/lapack/geev'
eig=: geev_jlapack_

create=: monad define
('M';'C')=:y
update''
)

destroy=: codestroy

update=: monad define
IC=: %. C             NB. inverse of the covariance matrix
CC=: chol C           NB. cholesky decomposition of the covariance matrix
t=. eig IC
L=: (>1{t) * id #C    NB. eigen values on a diag matrix
R=: >0{t              NB. rotation matrix (eigen vectors) IC -: R mp L mp |:R
)

NB. https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
NB. y random vectors of a #M-dimensional multivariate normal distribution
randmultin=: 3 : 'M +"1 CC MP"(_ 1) randn y , #M'

el=: dyad define NB. ordinate of an ellipse in std form
'a b'=. x
%: (*:b) * 1 - (*:y) % *:a
)

cellipse=: monad define NB. concentration ellipse at (%:k)-sigma
'u v'=. L MP 1#~#C             NB. std form ellipse (u -: %*:a) *. (v -: %*:b)
                               NB. 1 -: (u**:x) + v**:y
k=. 4
'a b'=. %: % u,v               NB. ellipse with axis a and b
elab=. (a,b)&el
n=. 10                         NB. half the number of points in the positive quadrant
absc0x=. }. steps 0 , a , +:n  NB. abscisses of the points in the positive quadrant
ordi0x=. elab absc0x           NB. ordinate of the points in the positive quadrant
absc=. (|. - absc0x) , 0 , absc0x
orditop=. (|. ordi0x) , (elab 0) , ordi0x
ordiall=. orditop , - |. orditop
coord=. (absc , |. absc) ,"0 ordiall
ncoord=. M +"1 (%:k) * R MP"(2 1) coord NB. transform the std form ellipse into the new basis defined by M and IC
<"1 |: ncoord
)

cocurrent 'base'
```

## Gaussian Mixture Model

```
NB.GMM.......................................................................

NB. make a multivariate normal distribution with mean 1{y and covariance 2{y
NB. draw 0{y objects from this distribution
mkobj=: dyad define
('n',":x)=: conew 'RandN'
n=. ". 'n',":x
create__n }.y
('x',":x)=: randmultin__n >{.y
)

mkdataset=: monad define
isomode=. i.#y
isomode mkobj"0 1 y
X=: ". }: , ([: ,&',' [: 'x'&, ":)"(0) isomode
perm=: ?~ #X
X=: perm { X NB. random permutation of the dataset
truth=: perm { isomode #~ ". }: , ([: ,&'),' [: '(#x'&, ":)"(0) isomode
)

initdataset=: monad define
mode=: ~.,>class
iclass=: (mode&notin)each class NB. modes not in a class ("Inverse" of class) 
K=: #mode
dim=:{:$X
)

NB. y is a boxed list of the number of examples known a priori for each class
initteacher=: monad define
rsel=. [ {~ ] ? #@[ NB. random selection of y elements of x
teacher=: (truth&([: I. [: +./ ="1 0)each class) rsel each y
)

dataset1=: 3 : 0 
mkdataset (50;(0 2);2 2$0.5 0 0 0.5),:(1000;(_0.5 2);2 2$3 0.6 0.6 1)
nbclass=: 2
trueclass=: class=: (,0);(,1)
initdataset''
initteacher (10;200)
)

dataset2=: 3 : 0
mkdataset (25;( 3  9);2 2$0.5 0 0 1),(20;(10  6);2 2$0.5 0 0 1),(5;(17 16);2 2$0.5 0 0 1),(20;( 3 12);2 2$1 0 0 0.5),(20;(12  6);2 2$1 0 0 0.5),:(10;(17 13);2 2$1 0 0 0.5)
nbclass=: 2
trueclass=: class=: (0 1 2);(3 4 5)
initdataset''
initteacher (25;25)
)

NB. Compute Initial Means à la kmeans++
CIM=: (],rndcenter)^:(]`(<:@:[)`(,:@:seed@:]))
seed=: {~ ?@#
NB. generate x rnd integers in i.#y with probability proportional to list of weights y
wghtprob=: 1&$: :((% {:)@:(+/\)@:] I. [ ?@$ 0:)"0 1
dist=: +/&.:*:@:-"1
rndcenter=: [ {~ [: wghtprob [: <./ dist/~

NB. Update F, the Fuzzy to crisp association between data and models,
NB. given the prior knowledge of a teacher
UF=: monad define
if. #,>teacher do.
toprob=. [: (%"1+/)@as {
merge=. (<@;)"1 @: |:
tocp=. 1 : '([: ((,@:u) ; (<@,)) CP)"1' NB. apply u to the cartesian product of x and y
m0=. merge > class  (toprob&F) tocp each teacher
m1=. merge > iclass (0:"0)     tocp each teacher
F=: F ((>@{.@])`(>@{:@])`([))} merge m0,:m1
end.
)

init=: monad define
t=:0 NB. time
M=: K CIM X
C0=: (+/%#) */~"1 (] -"1 +/%#) X
('detmin';'detmax')=: (1e_4&* ; 1e4&*) det C0
CSensor=: C0%25 NB. covariance corresponding to the sensor precision
C=: K#,:C0
F=: K#,:(#X)#%K NB. initial Fuzzyness
UF''
NB.F=: (="1 0 /:~@~.) truth NB. perfect teacher
sinned=: $0
conv=:0 NB. convergence, boolean
)

iter=: monad define
N=: +/"1 F
R=: N % #X
D=: (K#,:X) -"1 M
MLC=: N %~ +/"3 F * */~"1 D NB. max likelihood estimate of the covariances
PC=: C
C=: CSensor sinned} MLC
detC=: det C
sin=: I. (<&detmin +. >&detmax) detC
sinned=: sinned , sin
MSinned=: sinned { M
C=: CSensor sin} C
invC=: %.C
detC=: (det ` (detC"_) @. (0=#sin)) C
exp=: ^ --:1 * invC QF D
PDF=: exp * ((o.2)^--:dim) * %%:detC
F=: (%"1 +/) R*PDF
UF''
PM=: M
M=: MSinned sinned} N %~ +/"2 (K#,:X) *"(3 2) F
tconv=: [: *./@, [ > |@-/@] NB. test convergence
conv=: (1e_1 tconv M,:PM) *. 1e_2 tconv C,:PC
t=: >:t
)

run=: monad define
while. (-.conv) *. t<100 do. iter'' end.
CE=: - +/^:2 > ([: (* ^.) [: +/ {&F) each class NB. classification entropy
PRAUC''
)

metarun0=: monad define
init''
run''
('LM';'LC';'LCE';'LAUC')=: (LM,M);(LC,C);(LCE,CE);(LAUC,AUC)
)

metarun1=: monad define
init''
run''
('LM';'LC';'LCE';'LAUC')=: (,: M);(,: C);(,: CE);(,: AUC)
metarun0^:9''
NB. keep the results of the run with minimum classification entropy
('M';'C';'CE';'AUC')=: ((i. <./) LCE)&{ each LM;LC;LCE;LAUC
)

metarun2=: monad define
mkclass=: ([: i. each ;/) +each ([: ;/ 0: , [: }: +/\)
classmask=: 0: = ; @: (i.each) @ ;
modeinc=: nbmode=: nbclass # 1
PCE=: _ [ CE=: 1e6
('BM';'BC';'BCE';'BAUC';'BClass';'BMode')=: 6$a:
nbiter=: 0
while. (0 < +/modeinc) *. (PCE>CE) *. nbiter<10 do.
  class=: mkclass nbmode
  initdataset''
  PCE=: CE
  metarun1''
  ('BM';'BC';'BCE';'BAUC';'BClass';'BMode')=: (BM,<M);(BC,<C);(BCE,<CE);(BAUC,<AUC);(BClass,<class);<(BMode,<mode)
  modeinc=: > *./each (classmask nbmode) (<;.1) N>3 NB. ISO classes whose # of modes can increase
  nbmode=: nbmode + modeinc
  nbiter=: >: nbiter
end.
('M';'C';'CE';'AUC';'class';'mode')=: >@(((i. <./) >BCE)&{) each BM;BC;BCE;BAUC;BClass;<BMode
)

PRAUC=: monad define NB. area under the precision-recall curve
LR=:(%"1 +/)PDF NB. likelihood ratios
mask=: 0 (;teacher)} 1 #~ #truth NB. to remove the data points known to the teacher
LR=: mask #"1 LR
class0=: +./ (mask#truth) ="1 0 >{.trueclass
TP=: 3 : '+/(+./ (>0{class) { LR>y) *. class0'
FP=: 3 : '+/(+./ (>0{class) { LR>y) *. -.class0'
FN=: 3 : '(+/class0)-TP y'
precision=: TP % TP + FP
recall=: TP % TP + FN
clean=. (([: ~. {."1) |:@,: {."1 >.//. {:"1) NB. for equal values of recall keep the greater precision
NB. data for the parametric precision-recall curve. For recall 0 the precision is 1.
AUCData=: clean 0 1 ,~ }: (recall,precision)"0 steps 0 1 100
AUC=: +/ 2 (|@-/ (-:@*/@[ + {.@[ * ]) {:@{.)\ AUCData NB. area under the PR curve
NB.'stick,line' plot <"1|:AUCData
)

end0=: monad define
n=. ". 'n',":y
destroy__n''
)

end=: monad define
end0"0 mode
)

class_color=:'blue';'red'
class_style=:'markersize 0.1';'markersize 0.1'
estimate_color=:'green';'yellow'

draw_class=: monad define
color=: >y{class_color
style=: >y{class_style
draw_mode"0 >y{trueclass
)

draw_mode=: monad define
obj=. ". 'n',":y
dat=. ". 'x',":y
pd 'color ',color
pd 'type line ; pensize 3'
pd cellipse__obj''
pd 'type marker ; markers circle'
pd style
pd <"1 |: dat
)

draw_estimate=: monad define
n=. conew 'RandN'
create__n (y{M);(y{C)
pd 'color green'
pd 'type line ; pensize 3'
pd cellipse__n''
destroy__n''
)

draw=: monad define
pd 'reset'
draw_class"0 i.#class
draw_estimate"0 i.#mode
pd 'show'
)
```


