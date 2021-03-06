0 : 0
Copyright 2018 Pierre-Edouard Portier
peportier.me

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
)


load 'plot trig numeric'

NB.utils......................................................................

NB. if y is an empty array, return an empty list, otherwise, return y
as_z_=: ]`(($0)"_)@.(0:e.$)
B1_z_=: <"1               NB. rank-1 box
CP_z_=:{@(,&<)            NB. cartesian product

id_z_=: =@i.              NB. identity matrix of size y
MP_z_=: +/ . *            NB. matrix product
det_z_=: -/ . *           NB. determinant
QF_z_=: ] MP"1 [ MP"2 1 ] NB. quadratic form


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
mode=: ;class
iclass=: (mode&-.)each class NB. modes not in a class ("Inverse" of class) 
K=: #mode
dim=:{:$X
)

NB. y is a boxed list of the number of examples known a priori for each class
initteacher=: monad define
rsel=. [ {~ ] ? #@[ NB. random selection of y elements of x
teacher=: (truth&([: I. [: +./ ="1 0)each class) rsel each y
)

dataset1=: 3 : 0 
mkdataset (50;(0 2);2 2$0.5 0 0 0.5),:(50;(_2 4);2 2$1 0.5 0.5 1)
nbclass=: 2
trueclass=: class=: (,0);(,1)
initdataset''
initteacher (0;0)
)

dataset2=: 3 : 0 
mkdataset (50;(0 2);2 2$0.5 0 0 0.5),:(1000;(_0.5 2);2 2$3 0.6 0.6 1)
nbclass=: 2
trueclass=: class=: (,0);(,1)
initdataset''
initteacher (10;0)
)

dataset3=: 3 : 0
mkdataset (25;( 3  9);2 2$0.5 0 0 1),(20;(10  6);2 2$0.5 0 0 1),(5;(17 16);2 2$0.5 0 0 1),(20;( 3 12);2 2$1 0 0 0.5),(20;(12  6);2 2$1 0 0 0.5),:(10;(17 13);2 2$1 0 0 0.5)
nbclass=: 2
trueclass=: class=: (0 1 2);(3 4 5)
initdataset''
initteacher (25;25)
)

end0=: monad define
n=. ". 'n',":y
destroy__n''
)

end=: monad define
end0"0 ;trueclass
)


class_color=:'blue';'red'
class_style=:'markersize 0.1';'markersize 0.1'

plot_class=: monad define
color=: >y{class_color
style=: >y{class_style
plot_mode"0 >y{trueclass
)

plot_mode=: monad define
obj=. ". 'n',":y
dat=. ". 'x',":y
pd 'color ',color
pd 'type line ; pensize 3'
pd cellipse__obj''
pd 'type marker ; markers circle'
pd style
pd <"1 |: dat
)

plot_dataset=: monad define
pd 'reset'
plot_class"0 i.#trueclass
pd 'show'
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


NB. Compute Initial Means à la kmeans++
CIM=: (],rndcenter)^:(]`(<:@:[)`(,:@:seed@:]))
seed=: {~ ?@#
NB. generate x rnd integers in i.#y with probability proportional to list of weights y
wghtprob=: 1&$: :((% {:)@:(+/\)@:] I. [ ?@$ 0:)"0 1
dist=: +/&.:*:@:-"1
rndcenter=: [ {~ [: wghtprob [: <./ dist/~

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


run=: monad define
while. (-.conv) *. t<100 do. iter'' end.
CE=: - +/^:2 > ([: (* ^.) [: +/ {&F) each class NB. classification entropy
PRAUC''
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


plot_estimate=: monad define
n=. conew 'RandN'
create__n (y{M);(y{C)
pd 'color green'
pd 'type line ; pensize 3'
pd cellipse__n''
destroy__n''
)

plot_all=: monad define
pd 'reset'
plot_class"0 i.#trueclass
plot_estimate"0 i.#mode
pd 'show'
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

