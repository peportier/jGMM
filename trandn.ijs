load 'plot trig numeric'
load '~home/documents/shared/j/mixtures/math.ijs'
load '~home/documents/shared/j/mixtures/randn.ijs'

N1=: conew 'RandN'
create__N1 (1 1);2 2$3 0.6 0.6 1

N2=: conew 'RandN'
create__N2 (2 _1);2 2$1 0 0 0.5

X1=: randn2d__N1 10
X2=: randn2d__N2 100
X=: X1 , X2

QF=: ] mp"1 [ mp"2 1 ] NB. quadratic form
AME=: +/%#
e=:each

INIT=: monad define
K=: 2
d=:{:$X
t=:0
dt=:0.001
c=: 0.1
M=: (AME + [: ?@>:@>. 5 %~ >.@>./ - <.@<./) e K#< X
C0=: (+/%#) */~"1 (] -"1 +/%#) X
C=: K#< C0
F=: (];]) (#X)#%K NB. initial fuzzyness
NB.N=: +/ e F
N=: 10;100
L=: __
CNV=:0
draw''
)

I=: monad define
D=: X&(-"1) e M
NB.N=: +/ e F
R=: N %e #X
MLC=: N (+/@] % [)e F ([ * */~"1@])e D
PC=: C
NB.C=: MLC +e C0 *e ^-c*t
C=: MLC
SIN=: I. > (AME^:2@:| < 1e_3"_)e C NB. index of the cov matrixes becoming too small
C=: (<C0) SIN} C
IC=: %.e C
DC=: det e C
EXP=: (^ @ *&(-%2))e IC QF e D
PDF=: EXP *e (*&((o.2)^--:d) @ %@%:)e DC
F=: <"1 (%"1 +/"2) > R *e PDF NB. (+/ -: {:@$ # 1:)> F
DM=: +/e F *e D mp e IC
PM=: M
M=: M (+ *&dt)e DM
t=: >:t
L=: +/^.+/"2 >PDF
CNVC=: >(1e_3 <~ [: AME^:2 |) e C -e PC
CNVP=: >(1e_3 <~ [: AME |) e M -e PM
CNV=: *./ CNVC , CNVP
)

RUN=: monad define
while. (-.CNV) *. t<1000 do. I'' end.
draw''
)

draw_data=: monad define
pd 'reset'
pd 'color blue'
pd 'type line'
pd cellipse__N1''
pd 'type marker ; markers circle ; markersize 0.2'
pd <"1 |: X1

pd 'color red'
pd 'type line'
pd cellipse__N2''
pd 'type marker ; markers circle ; markersize 0.2'
pd <"1 |: X2
pd 'show'
)

draw_estimate=: monad define
NT1=: conew 'RandN'
NT2=: conew 'RandN'
create__NT1 (>0{M);(>0{C)
create__NT2 (>1{M);(>1{C)
pd 'color green'
pd 'type line'
pd cellipse__NT1''
pd 'color yellow'
pd cellipse__NT2''
pd 'show'
destroy__NT1''
destroy__NT2''
)

draw=: monad define
draw_data''
draw_estimate''
)

NB. 2-sigma concentration ellipse
NB. in=: ({.@] < [)*.([ < {:@]) NB. x in (y1 y2)
NB. g=: 3 : 'in&(4 (-,+) 0.1) (y-M__N2) mp IC__N2 mp (y-M__N2)'
NB. X=: Y=: steps _5 5 100
NB. Z=: X (g@:,"(0))/ Y
NB. 'type surface' plot X;Y;Z

END=: monad define
destroy__N1''
destroy__N2''
)