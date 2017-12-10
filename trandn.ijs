load 'plot trig numeric'
load '~home/documents/shared/j/mixtures/math.ijs'
load '~home/documents/shared/j/mixtures/randn.ijs'

DATASET=: monad define
N0=: conew 'RandN'
create__N0 (4 _1);2 2$0.5 0 0 0.5
X0=: 2 randmultin__N0 100
N1=: conew 'RandN'
create__N1 (1 0.5);2 2$3 0.6 0.6 1
X1=: 2 randmultin__N1 100000
X=: X0 , X1
)

UF=: ] ` ( ((K*#T) $ 1,(<:K)#0) & ( (,|:(i.K) CP T)} ) ) @. (0<#T)

CIM=: (],rndcenter)^:(]`(<:@:[)`(,:@:seed@:])) NB. compute initial means à la kmeans++
seed=: {~ ?@#
rndcenter=: [ {~ [: wghtprob [: <./ dst/~

INIT=: monad define
K=: 2
d=:{:$X
t=:0
M=: K CIM X
C0=: (+/%#) */~"1 (] -"1 +/%#) X
C=: K#,:C0
F=: K#,:(#X)#%K NB. initial fuzzyness
T=: i.6 NB. teacher
F=: UF F
NB.F=: (((#X0)#0) , (#X1)#1) ,: ((#X0)#1) , (#X1)#0 NB. perfect teacher
N=: +/"1 F
N=: (#X0) , (#X1) NB. perfect knowledge of the ratio
L=: __
CNV=:0
draw''
)

I=: monad define
D=: (K#,:X) -"1 M
NB.N=: +/"1 F
R=: N % #X
MLC=: N %~ +/"3 F * */~"1 D
PC=: C
C=: MLC
SIN=: I. (AME"1^:2@:| < 1e_1"_) C
C=: C0 SIN} C
IC=: %.C
DC=: det C
EXP=: ^ --:1 * IC QF D
PDF=: EXP * ((o.2)^--:d) * %%:DC
F=: (%"1 +/) R*PDF
F=: UF F
PM=: M
M=: N %~ +/"2 (K#,:X) *"(3 2) F
L=: +/^.+/ R*PDF
CNV=: *./> (C;M) ([: *./@:,@:<&1e_1@:| -)e (PC;PM)
t=: >:t
)

RUN=: monad define
while. (-.CNV) *. t<100 do. I'' end.
draw''
)

END=: monad define
destroy__N1''
destroy__N2''
)

class_color=:'blue';'red'
class_style=:'markersize 0.3';'markersize 0.15'
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