load 'plot trig numeric'
load '~home/documents/shared/j/mixtures/utils.ijs'
load '~home/documents/shared/j/mixtures/randn.ijs'

NB. make a multivariate normal distribution with mean 1{y and covariance 2{y
NB. draw 0{y objects from this distribution
mkobj=: dyad define
('n',":x)=: conew 'RandN'
n=. ". 'n',":x
create__n }.y
('x',":x)=: randmultin__n >{.y
)

mkdataset=: monad define
iomodes=. i.#y NB. indexes of the modes
iomodes mkobj"0 1 y
xs=: ". }: , ([: ,&',' [: 'x'&, ":)"(0) iomodes
perm=: ?~ #xs
xs=: perm { xs NB. the all dataset in a random permutation
truth=: perm { iomodes #~ ". }: , ([: ,&'),' [: '(#x'&, ":)"(0) iomodes
)

initdataset=: monad define
mode=: ~.,>class
iclass=: (mode&notin)e class NB. modes not in a class ("Inverse" of class) 
k=: #mode
dim=:{:$xs
)

dataset1=: 3 : 0 
mkdataset (50;(0 2);2 2$0.5 0 0 0.5),:(1000;(_0.5 2);2 2$3 0.6 0.6 1)
class=: (,0);(,1)
initdataset''
)

dataset2=: 3 : 0
mkdataset (25;( 3  9);2 2$0.5 0 0 1),(20;(10  6);2 2$0.5 0 0 1),(5;(17 16);2 2$0.5 0 0 1),(20;( 3 12);2 2$1 0 0 0.5),(20;(12  6);2 2$1 0 0 0.5),:(10;(17 13);2 2$1 0 0 0.5)
class=: (0 1 2);(3 4 5)
initdataset''
)

NB. Compute Initial Means à la kmeans++
cim=: (],rndcenter)^:(]`(<:@:[)`(,:@:seed@:]))
seed=: {~ ?@#
rndcenter=: [ {~ [: wghtprob [: <./ dst/~

NB. Update F, the Fuzzy to crisp association between data and models,
NB. given the prior knowledge of a teacher
uf=: monad define
if. #,>teacher do.
toprob=: [: (%"1+/)@as {
merge=: (<@;)"1 @: |:
NB. apply u to the cartesian product of x and y,
tocp=: 1 : '([: ((,@:u) ; (<@,)) CP)"1'
m0=: merge > class  (toprob&F) tocp e teacher
m1=: merge > iclass (0:"0)     tocp e teacher
f=: f ((>@{.@])`(>@{:@])`([))} merge m0,:m1
end.
)

NB. update N, the Number of objects associated with each mode,
NB. given the prior knowledge of some of the ratios
un=: monad define
bn=. ({"0 1)&n e class           NB. boxed # of objects (one box per class)
mr=. (%+/)e bn                   NB. ratio of the modes within a class
nc=. > ([: +/ ] {"0 1 n"_) e class NB. # of objects per class
rnc=. ar * +/n                   NB. rectified # of objects per class
irc=. I.*rnc                     NB. IO rectified classes
icc=. (i.#class) notin irc       NB. IO classes used for compensation
drm=. ; (irc{class) ([: |: ,:)e (irc{mr) *e irc{rnc-nc NB. deltas for the rectified modes
cost=. - +/ {:"1 drm             NB. cost of the rectifications, to be compensated with the remaining modes
ccr=. (%+/) ; +/e icc { bn       NB. ratios of the classes used for compensation
dcm=. ; (icc{CLS) ([: |: ,:)e (icc { mr) *e cost * ccr NB. deltas for the compensating modes
bion=: B1@{."1                   NB. extract boxed indexes of n from dcm or drm
n=: (drm,dcm) (({:"1@[ + bion@[ { ]) ` (bion@[) ` ])} n
)

init=: monad define
t=:0 NB. time
rsel=. [ {~ ] ? #@[ NB. random selection of y elements of x
teacher=: (truth&([: I. [: +./ ="1 0)e class) rsel e (10;200)
m=: k cim xs
c0=: (+/%#) */~"1 (] -"1 +/%#) xs
c=: k#,:c0
f=: k#,:(#xs)#%k NB. initial Fuzzyness
uf''
NB.f=: (="1 0 /:~@~.) truth NB. perfect teacher
ar=: 0.05 0 NB. a priori knowledge of the class ratio
conv=:0 NB. convergence, boolean
draw''
)

iter=: monad define
d=: (k#,:xs) -"1 m
n=: +/"1 f
un''
r=: n % #xs
mlc=: n %~ +/"3 f * */~"1 d NB. maximum likelihood estimate of the covariances
pc=: c
c=: mlc
detc=: det c
sin=: I. detc < 1e_1 NB. indicator for the covariances that became SINgular
c=: c0 sin} c
invc=: %.c
detc=: (det ` (detc"_) @. (0=#sin)) c
exp=: ^ --:1 * invc QF d
pdf=: exp * ((o.2)^--:dim) * %%:detc
f=: (%"1 +/) r*pdf
uf''
pm=: m
m=: n %~ +/"2 (k#,:xs) *"(3 2) f
conv=: *./> (c;m) ([: *./@:,@:<&1e_2@:| -)e (pc;pm)
t=: >:t
)

run=: monad define
while. (-.conv) *. t<100 do. iter'' end.
draw''
prauc''
)

occ=: monad define NB. operating characteristics curve
lr=:(%"1 +/)pdf NB. likelihood ratios
pt=. 3 : '((-/ % (#x1)"_) , {: % (#x0)"_) (#,+/) <&(#x0) I. ({.lr)>y'
plot <"1|:pt"(0) steps 0 1 50
)

prauc=: monad define NB. area under the precision-recall curve
lr=:(%"1 +/)pdf NB. likelihood ratios
class0=: +./ truth ="1 0 >{.class
tp=: 3 : '+/({.lr>y) *. class0'
fp=: 3 : '+/({.lr>y) *. -.class0'
fn=: 3 : '(+/class0)-tp y'
pr=: tp % tp + fp
rc=: tp % tp + fn
data=. 0 1 ,~ }: (rc,pr)"0 steps 0 1 100 NB. data for the parametric precision-recall curve. For recall 0 the precision is 1.
+/ 2 (|@-/ (-:@*/@[ + {.@[ * ]) {:@{.)\ data
NB.'stick,line' plot <"1|:(rc,pr)"0 steps 0 1 10
)

end0=: monad define
n=: ". 'n',":y
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
draw_mode"0 >y{class
)

draw_mode=: monad define
obj=: ". 'n',":y
dat=: ". 'x',":y
pd 'color ',color
pd 'type line ; pensize 3'
pd cellipse__obj''
pd 'type marker ; markers circle'
pd style
pd <"1 |: dat
)

draw_estimate=: monad define
nt=. conew 'RandN'
create__nt (y{m);(y{c)
pd 'color green'
pd 'type line ; pensize 3'
pd cellipse__nt''
destroy__nt''
)

draw=: monad define
pd 'reset'
draw_class"0 i.#class
draw_estimate"0 i.#mode
pd 'show'
)