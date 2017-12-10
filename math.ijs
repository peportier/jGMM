mp=: +/ . *
det=: -/ . *
id=: =@i.
QF=: ] mp"1 [ mp"2 1 ] NB. quadratic form
AME=: +/%#
CP=:{@(,&<) NB. cartesian product
wghtprob=: 1&$: :((% {:)@:(+/\)@:] I. [ ?@$ 0:)"0 1 NB. rnd gen from lst of weights
dst=: +/&.:*:@:-"1
e=:each