load 'plot trig numeric'
load '~home/documents/shared/j/mixtures/utils.ijs'
load '~home/documents/shared/j/mixtures/randn.ijs'

DATASET=: monad define
N0=: conew 'RandN'
create__N0 (2 0);2 2$0.5 0 0 0.5
X0=: 2 randmultin__N0 50
N1=: conew 'RandN'
create__N1 (1 0.5);2 2$3 0.6 0.6 1
X1=: 2 randmultin__N1 1000
X=: X0 , X1
)

CIM=: (],rndcenter)^:(]`(<:@:[)`(,:@:seed@:])) NB. compute initial means à la kmeans++
seed=: {~ ?@#
rndcenter=: [ {~ [: wghtprob [: <./ dst/~

UF=: monad define
if. #,>T do.
toprob=: [: (%"1+/)@as {
merge=: (<@;)"1 @: |:
applyToCP=: 1 : '([: ((,@:u) ; (<@,)) CP)"1'
m0=: merge > CLS  (toprob&F) applyToCP e T
m1=: merge > ICLS (0:"0)     applyToCP e T
F=: F ((>@{.@])`(>@{:@])`([))} merge m0,:m1
end.
)

UN=: monad define
bn=. ({"0 1)&N e CLS             NB. boxed # of objects (one box per class)
mr=. (%+/)e bn                   NB. ratio of the modes within a class
nc=. > ([: +/ ] {"0 1 N"_) e CLS NB. # of objects per class
rnc=. AR * +/N                   NB. rectified # of objects per class
irc=. I.*rnc                     NB. IO rectified classes
icc=. (i.#CLS) notin irc         NB. IO classes used for compensation
drm=. ; (irc{CLS) ([: |: ,:)e (irc{mr) *e irc{rnc-nc NB. deltas for the rectified modes
cost=. - +/ {:"1 drm             NB. cost of the rectifications, to be compensated with the remaining modes
ccr=. (%+/) ; +/e icc { bn       NB. ratios of the classes used for compensation
dcm=. ; (icc{CLS) ([: |: ,:)e (icc { mr) *e cost * ccr NB. deltas for the compensating modes
bion=: B1@{."1                   NB. extract boxed indexes of n from dcm or drm
N=: (drm,dcm) (({:"1@[ + bion@[ { ]) ` (bion@[) ` ])} N
)

INIT=: monad define
CLS=: (,0);,1 NB. classes
MOD=: ~.,>CLS NB. modes
ICLS=: (MOD&notin)e CLS NB. modes not in a class ("Inverse" of class) 
K=: 2
d=:{:$X
t=:0 NB. time
T=: (i.5) ; (#X0)+i.100 NB. teacher
M=: K CIM X
C0=: (+/%#) */~"1 (] -"1 +/%#) X
C=: K#,:C0
F=: K#,:(#X)#%K NB. initial Fuzzyness (Fuzzy to crisp association between data and models)
UF''
NB.F=: (((#X0)#1) , (#X1)#0) ,: ((#X0)#0) , (#X1)#1 NB. perfect teacher
AR=: 0.05 0 NB. a priori knowledge of the ratio
CNV=:0
draw''
)

I=: monad define
D=: (K#,:X) -"1 M
N=: +/"1 F
UN''
R=: N % #X
MLC=: N %~ +/"3 F * */~"1 D
PC=: C
C=: MLC
DC=: det C
SIN=: I. DC < 1e_1
C=: C0 SIN} C
IC=: %.C
DC=: (det ` (DC"_) @. (0=#SIN)) C
EXP=: ^ --:1 * IC QF D
PDF=: EXP * ((o.2)^--:d) * %%:DC
F=: (%"1 +/) R*PDF
UF''
PM=: M
M=: N %~ +/"2 (K#,:X) *"(3 2) F
CNV=: *./> (C;M) ([: *./@:,@:<&1e_2@:| -)e (PC;PM)
t=: >:t
)

RUN=: monad define
while. (-.CNV) *. t<100 do. I'' end.
draw''
)

OC=: monad define NB. operating characteristics
LR=:(%"1 +/)PDF NB. likelihood ratios
OCPT=: 3 : '((-/ % (#X1)"_) , {: % (#X0)"_) (#,+/) <&(#X0) I. ({.LR)>y'
plot <"1|:OCPT"(0) steps 0 1 50
)

END=: monad define
destroy__N1''
destroy__N2''
)

class_color=:'blue';'red'
class_style=:'markersize 0.1';'markersize 0.1'
estimate_color=:'green';'yellow'

draw_class=: monad define
color=. >y{class_color
obj=: ". 'N',":y
dat=: ". 'X',":y
pd 'color ',color
pd 'type line ; pensize 3'
pd cellipse__obj''
pd 'type marker ; markers circle'
pd >y{class_style
pd <"1 |: dat
)

draw_estimate=: monad define
color=. >y{estimate_color
NT=. conew 'RandN'
create__NT (y{M);(y{C)
pd 'color ',color
pd 'type line ; pensize 3'
pd cellipse__NT''
destroy__NT''
)

draw=: monad define
pd 'reset'
draw_class"0 i.K
draw_estimate"0 i.K
pd 'show'
)