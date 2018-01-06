NB. 2-sigma concentration ellipse
NB. in=: ({.@] < [)*.([ < {:@]) NB. x in (y1 y2)
NB. g=: 3 : 'in&(4 (-,+) 0.1) (y-M__N2) mp IC__N2 mp (y-M__N2)'
NB. X=: Y=: steps _5 5 100
NB. Z=: X (g@:,"(0))/ Y
NB. 'type surface' plot X;Y;Z

NB.L=: +/^.+/ R*PDF