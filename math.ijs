mp=: +/ . *
det=: -/ . *
id=: =@i.
QF=: ] mp"1 [ mp"2 1 ] NB. quadratic form
AME=: +/%#
CP=:{@(,&<) NB. cartesian product
wghtprob=: 1&$: :((% {:)@:(+/\)@:] I. [ ?@$ 0:)"0 1 NB. rnd gen from lst of weights
dst=: +/&.:*:@:-"1
in=: ({.@] < [)*.([ < {:@]) NB. x in (y1 y2)
putin=: ({:@] ` ({.@]) @. ([<({.@]))) ` [ @. in NB. e.g. 5 -: 3 putin 5 15
pstv=: 0: ` ] @. (>&0)
B1=: <"1
ihole=: ({. , ] {.~ [: - #@] - [ + 1:) i. NB. 0 1 3 4 -: 2 igap 5
e=:each