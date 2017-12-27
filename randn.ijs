NB. draw values from a multivariate normal distribution
NB. with mean vector M and covariance matrix C
NB. find coordinates of the 2-sigma concentration ellipse
coclass 'RandN'

load '~home/documents/shared/j/mixtures/utils.ijs'

randu=: (?@$&0) :((p.~ -~/\)~ $:)
rande=: -@^.@randu : (* $:)
randn=: (($,) +.@:(%:@(2&rande) r. 0 2p1&randu)@>.@-:@(*/)) : (p. $:)

require 'math/lapack'
require 'math/lapack/potrf'
chol=: potrf_jlapack_
require 'math/lapack/geev'
eig=: geev_jlapack_

create=: monad define
  ('M';'C')=:y
  update''
)

destroy =: codestroy

update=: monad define
  IC=: %. C
  CC=: chol C
  t=. eig IC
  L=: (>1{t) * id #C    NB. eigen values on a diag matrix
  R=: >0{t              NB. rotation matrix (eigen vectors)
                        NB. IC -: R mp L mp |:R
)

NB. https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
NB. y random vectors of a #M-dimensional multivariate normal distribution
randmultin=: 3 : 'M +"1 CC mp"(_ 1) randn y , #M'

el=: dyad define NB. ordinate of an ellipse in std form
  'a b'=. x
  %: (*:b) * 1 - (*:y) % *:a
)

cellipse=: monad define NB. concentration ellipse at (%:k)-sigma
  'u v'=. L mp 1#~#C             NB. std form ellipse (u -: %*:a) *. (v -: %*:b)
                                 NB. 1 -: (u**:x) + v**:y
  k=. 4
  'a b'=. %: % u,v               NB. ellipse with axis a and b
  elab=. (a,b)&el
  n=. 10                         NB. half the number of points in the positive quadrant
  absc0x=. }. steps 0 , a , 2*n  NB. abscisses of the points in the positive quadrant
  ordi0x=. elab absc0x           NB. ordinate of the points in the positive quadrant
  absc=. (|. - absc0x) , 0 , absc0x
  orditop=. (|. ordi0x) , (elab 0) , ordi0x
  ordiall=. orditop , - |. orditop
  coord=. (absc , |. absc) ,"0 ordiall
  ncoord=. M +"1 (%:k) * R mp"(2 1) coord NB. transform the std form ellipse into the new basis defined by M and IC
  <"1 |: ncoord
)

cocurrent 'base'