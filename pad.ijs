NB. MM=:M -e PM
NB. CC=: C (-:@+)e PC
NB. ICC=: %.e CC
NB. B1=: *&(%8)e ICC QF e MM
NB. B2=: (*&(%2) @ ^.)e (det e CC) %e (% @ %: @ det)e PC mp e C
NB. B=: B1 +e B2 NB. Bhattacharyya distance between two successive iterations

NB. 2-sigma concentration ellipse
NB. in=: ({.@] < [)*.([ < {:@]) NB. x in (y1 y2)
NB. g=: 3 : 'in&(4 (-,+) 0.1) (y-M__N2) mp IC__N2 mp (y-M__N2)'
NB. X=: Y=: steps _5 5 100
NB. Z=: X (g@:,"(0))/ Y
NB. 'type surface' plot X;Y;Z

NB.CIM=: AME + [: (] * [: ? # # 0:) 5 %~ >./ - <./
NB.M=: CIM"2 K#,:X

NB.SIN=: I. (AME"1^:2@:| < 1e_1"_) C

NB.L=: +/^.+/ R*PDF

NB.UF=: 3 : 'if. #T do. ((K*#T) $ 1,(<:K)#0) (,|:(i.K) CP T)} y else. y end.'

NB. very unbalanced dataset
DATASET=: monad define
N0=: conew 'RandN'
create__N0 (4 _1);2 2$0.5 0 0 0.5
X0=: 2 randmultin__N0 100
N1=: conew 'RandN'
create__N1 (1 0.5);2 2$3 0.6 0.6 1
X1=: 2 randmultin__N1 100000
X=: X0 , X1
)

F=: F , {. F
F=: 0.1 0.2 0.3 * F

c=: (0 1);(2 3 4);,5   NB. classes
n=:  1 1   3 3 2   2   NB. # of objects per mode
ar=: 0.2   0       0.1 NB. a priori ratio, 0 when no a priori

((1 2);(2 3 2);,2) -: bn=: ({"0 1)&n e c    NB. boxed # of objects (one box per class)
mr=: (%+/)e bn                              NB. ratio of the modes within a class
3   7 2   -: nc=: > ([: +/ ] {"0 1 n"_) e c NB. # of objects per class
2.4 0 1.2 -: rnc=: ar * +/n                 NB. rectified # of objects per class
0     2   -: irc=: I.*rnc                   NB. IO rectified classes
1         -: icc=: (i.#c) notin irc         NB. IO classes used for compensation
drm=: ; (irc{c) ([: |: ,:)e (irc{mr) *e irc{rnc-nc NB. deltas for the rectified modes
cost=: - +/ {:"1 drm NB. cost of the rectifications, to be compensated with the remaining modes
ccr=: (%+/) ; +/e icc { bn NB. ratios of the classes used for compensation
dcm=: ; (icc{c) ([: |: ,:)e (icc { mr) *e cost * ccr NB. deltas for the compensating modes
bion=: B1@{."1 NB. extract boxed indexes of n from dcm or drm
n=: (drm,dcm) (({:"1@[ + bion@[ { ]) ` (bion@[) ` ])} n

LR=:(%"1 +/)PDF NB. likelihood ratios
+/({.LR>0.8) *. TRUTH
TP=: 3 : '+/({.LR>y) *. TRUTH'
FP=: 3 : '+/({.LR>y) *. -.TRUTH'
FN=: 3 : '(+/TRUTH)-TP y'
PR=: TP % TP + FP
RC=: TP % TP + FN