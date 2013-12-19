from __future__ import division
from math import *


half = 1/2

third = 1/3

twothirds = 2/3

quarter = 1/4


root_pi = sqrt(pi)



def inverse(x_):

    return 1/float(x_)

def squared(x_):

    return x_**2

def cubed(x_):

    return x_**3

def cubed_root(x_):

    return x_**third


def sgn(x_):

    return copysign(1,x_)

def heaviside_step(x_):

    return half*(1+sgn(x_))


def binomial(n_, k_):

    return factorial(n_)/(factorial(k_)*factorial(n_-k_))

def double_binomial(n_, k1_, k2_):

    return binomial(n_,k1_)*binomial(n_,k2_)

def multinomial(*klist):

    n = 0
    k_ = 1
    for k in klist:
        k_ = k_ * factorial(k)
        n = n + k

    return factorial(n)/k_

def factorial_approx(n_):

    fact = sqrt(pi)*(n_/e)**n_
    fact *= (((8*n_ + 4)*n_ + 1)*n_ + 1/30.)**(1./6.)
    if isinstance(n_, int):
        fact = int(fact)
    return fact
