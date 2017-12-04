NB. ---------------------------------------------------------
NB. coclass 'RandN'
NB. draw values from a multivariate normal distribution
NB. with mean vector M and covariance matrix C
NB. find coordinates of the 2-sigma concentration ellipse

NB. ---------------------------------------------------------
NB. randu
NB.
NB. Description:
NB.   Uniform distribution U(a,b) with support (a,b)
NB.
NB. Syntax:
NB.   S=. [supp] randu sh
NB. where
NB.   sh   - r-vector of non-negative integers, shape of S
NB.   supp - optional 2-vector (a,b), a<b, support of
NB.          distribution, open interval, default is (0 1)
NB.   S    - sh-array of values s ~ U(a,b)
NB.   r    ≥ 0, the rank of S
NB.
NB. Formula:
NB.   s ← a + (b-a)*u
NB. where
NB.   u ~ U(0,1)
NB.   s ~ U(a,b)

NB. ---------------------------------------------------------
NB. rande
NB.
NB. Description:
NB.   Exponential distribution E(μ) with mean μ
NB.
NB. Syntax:
NB.   S=. [μ] rande sh
NB. where
NB.   sh - r-vector of non-negative integers, shape of S
NB.   μ  > 0, optional mean, default is 1
NB.   S  - sh-array of values s ~ E(μ)
NB.   r  ≥ 0, the rank of S
NB.
NB. Formula:
NB.   s ← -μ*log(1-u)
NB. where
NB.   u ~ U(0,1)
NB.   s ~ E(μ)
NB.
NB. Notes:
NB. - randu is used instead of (1-randu) since they are
NB.   equivalent under following assumptions:
NB.   - randu's support is open interval (0,1)
NB.   - randu's values aren't used outside elsewhere

NB. ---------------------------------------------------------
NB. randn
NB.
NB. Description:
NB.   Normal distribution N(μ,σ^2) of real numbers
NB.   with mean μ and variance σ^2
NB.
NB. Syntax:
NB.   S=. [par] randnf sh
NB. where
NB.   sh  - r-vector of non-negative integers, shape of S
NB.   par - optional 2-vector (μ,σ), σ>0, parameters of
NB.         distribution, default is (0 1)
NB.   S   - sh-array of values s ~ N(μ,σ^2)
NB.   r   ≥ 0, the rank of S
NB.
NB. Formula:
NB.   nf1 ← sqrt(e)*cos(u)         NB. see...
NB.   nf2 ← sqrt(e)*sin(u)         NB. ...[1]
NB.   n3  ← σ*nxx + μ
NB. where
NB.   u            ~ U(0,2*π)
NB.   e            ~ E(2)
NB.   nf1,nf2,nc,n ~ N(0,1)
NB.   n3           ~ N(μ,σ^2)
NB.
NB. Notes:
NB. - randn requests only ⌈Π{sh[:]}/2⌉ numbers from RNG
NB.
NB. References:
NB. [1] G. E. P. Box, Mervin E. Muller. A Note on the
NB.     Generation of Random Normal Deviates. The Annals of
NB.     Mathematical Statistics, 1958, Vol. 29, No. 2, pp.
NB.     610-611.

NB. ---------------------------------------------------------
NB. randn2d
NB.
NB. https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution