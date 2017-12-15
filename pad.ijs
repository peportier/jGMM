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
