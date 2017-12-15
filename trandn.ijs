load 'plot trig numeric'
load '~home/documents/shared/j/mixtures/math.ijs'
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
if. #T do.
NB. each cell of CU contains Coordinates of F that must be Updated
NB. the first element of a cell is the one for which the teacher gives probability TP
NB. the other K-1 elements will share equally 1-TP
CU=. ({: ,~"0 {. , {. ihole K"_)"1 T
NB. each cell of NP contains the new probabilities
NB. for the elements of the corresponding cell of IM
NP0=. [ * F {~ [: B1 ({.@])
NP1=. (1-[ % (<:K)"_) * F {~ [: B1 }.@]
NP=. (% +/)"1 TP (NP0,NP1)"0 2 CU
F=: (,NP) (,B1 CU)} F
end.
)

INIT=: monad define
class=: 2 1 $ i.2
K=: 2
d=:{:$X
t=:0
M=: K CIM X
C0=: (+/%#) */~"1 (] -"1 +/%#) X
C=: K#,:C0
F=: K#,:(#X)#%K NB. initial fuzzyness
T=: (0 ,"0 i.0) , 1 ,"(0) (#X0)+i.0 NB. teacher
TP=:(#T)#1 NB. teacher confidence probability
UF''
NB.F=: (((#X0)#1) , (#X1)#0) ,: ((#X0)#0) , (#X1)#1 NB. perfect teacher
AR=: 0.05 NB. approximate apriori knowledge of the ratio
ARTOL=: 0.1 NB. tolerance on AR
NR=: 2 2 $ 0,#X
NR=: /:~"1 |: (AR ([ (pstv@-,+) *) ARTOL) ([: (-/ ,~ {:) ] , *)"0 #X NB. range of possible values for N
UN=: (putin"0 1)&NR NB. update N to put it in the range NR
CNV=:0
draw''
)

I=: monad define
D=: (K#,:X) -"1 M
N=: UN +/"1 F
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