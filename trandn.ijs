load 'plot trig numeric'
load '~home/documents/shared/j/mixtures/math.ijs'
load '~home/documents/shared/j/mixtures/randn.ijs'

DATASET=: monad define
N0=: conew 'RandN'
create__N0 (1 0.5);2 2$3 0.6 0.6 1
N1=: conew 'RandN'
create__N1 (2 1);2 2$0.5 0 0 0.5
X0=: 2 randmultin__N0 100
X1=: 2 randmultin__N1 10
X=: X0 , X1
)

INIT=: monad define
K=: 2
d=:{:$X
t=:0
M=: (AME + [: ?@>:@>. 5 %~ >.@>./ - <.@<./)"2 K#,:X
C0=: (+/%#) */~"1 (] -"1 +/%#) X
C=: K#,:C0
F=: K#,:(#X)#%K NB. initial fuzzyness
NB.N=: +/"1 F
N=: 10 100
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
SIN=: I. (AME"1^:2@:| < 1e_3"_) C
C=: C0 SIN} C
IC=: %.C
DC=: det C
EXP=: ^ --:1 * IC QF D
PDF=: EXP * ((o.2)^--:d) * %%:DC
F=: (%"1 +/) R*PDF
PM=: M
M=: N %~ +/"2 (K#,:X) *"(3 2) F
L=: +/^.+/ R*PDF
CNV=: *./> (C;M) ([: *./@:,@:<&1e_3 -)e (PC;PM)
t=: >:t
)

RUN=: monad define
while. (-.CNV) *. t<1000 do. I'' end.
draw''
)

END=: monad define
destroy__N1''
destroy__N2''
)

class_colors=:'blue';'red'
estimate_colors=:'green';'yellow'

draw_class=: monad define
color=. >y{class_colors
obj=: ". 'N',":y
dat=: ". 'X',":y
pd 'color ',color
pd 'type line'
pd cellipse__obj''
pd 'type marker ; markers circle ; markersize 0.2'
pd <"1 |: dat
)

draw_estimate=: monad define
color=. >y{estimate_colors
NT=. conew 'RandN'
create__NT (y{M);(y{C)
pd 'color ',color
pd 'type line'
pd cellipse__NT''
destroy__NT''
)

draw=: monad define
pd 'reset'
draw_class"0 i.K
draw_estimate"0 i.K
pd 'show'
)