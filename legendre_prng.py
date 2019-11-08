import numpy as np
from typing import *
from sympy.ntheory import isprime, factorint, generate
from bitarray import bitarray
from random import randint
from time import time
from dateutil.relativedelta import relativedelta as rd
from utils import cyclic_bitarray
fmt = '{0.minutes} m, {0.seconds} s'

# a nice example prime: 492876863
def legendre(a: int, p: int) -> int:
    """
    computes the legendre symbol of a mod p where p is an odd prime,
    algorithm taken from https://martin-thoma.com/how-to-calculate-the-legendre-symbol/
    """
    # try to speed this up using jacobi symbol or fast expenentiation 
    # then try to implemented in C using gmp / boost (it contains legendre and jacobi symbols)
    assert isprime(p) and p > 2, "ERROR: p has to be an odd prime"
    if a > p or a < 0:
        return legendre(a % p, p)
    elif a in [0, 1]:
        return a
    elif a == 2:
        if p % 8 in [-1, 1]:
            return 1
        else:
            return -1
    elif a == p - 1:
        if p % 4 == 1:
            return 1
        else:
            return -1
    elif not isprime(a):
        facs = factorint(a)
        return np.prod(np.array([legendre(f,p) ** facs[f] for f in facs.keys()]))
    else:
        if p % 4 == 1 or a % 4 == 1:
            return legendre(p,a)
        else:
            return - legendre(p,a)

def prng(key: int, seed: int, p: int, length: int) -> bitarray:
    # assert that n is not too large (10 * log2^10(p) ?!?)
    # is it better to return a bitarray?
    return bitarray([legendre(seed + key + i, p) == 1 for i in range(length)]) 

def bf_v0(p: int, hint=None, rand=False)-> int:
    """
    naive bruteforce:
    given the pseudorandom stream f_k(0), ... ,f_k(t) as hint
    returns the key thanks to a bruteforce search
    """
    if not hint:
        randkey = randint(10**3,10**5) if rand else 75162
        print('random' if rand else 'hardcoded' + ' key: ' + str(randkey))
        hint = prng(randkey, 0, p, 10**4)
    k = 0
    t = time()
    while True:
        if hint == prng(k, 0, p, len(hint)):
            print("execution time: " + str(time() - t))
            return k
        k += 1

def bf_v1(p: int, hint=None, par=2, rand=False)-> int:
    """
    improvement:
    avoid calculating the same legendre symbol too many times
    """
    if not hint:
        randkey = randint(0, p) if rand else 75162
        print('random ' if rand else 'hardcoded ')
        print('key: ' + str(randkey))
        hint = prng(randkey, 0, p, 10**4)
    k = 0
    t = time()
    while True:
        stream = prng(k, 0, p, len(hint))
        if hint == stream:
            print("execution time: " + str(time() - t))
            return k
        candidates = []
        for i in range(int(len(hint) / par)):
            if hint[i:] == stream[:-i]:
                candidates.append(k - i)
            if hint[:-i] == stream[i:]:
                candidates.append(k + i)
        for key in candidates:
            if hint == prng(key, 0, p, len(hint)):
                print("execution time: " + str(time() - t))
                return key
        k += 2 * int(len(hint) / par)
        print(k)

def bf_v2(p: int, hint=None, rand=False)-> int:
    """
    improvement: x10 speedup over v1
    not all symbols need to be calculated
    """
    if not hint:
        randkey = randint(0, p) if rand else 9075162
        print('key: ' + str(randkey))
        t = time()
        hint = prng(randkey, 0, p, 10**4)
        print("hint derivation time: " + fmt.format(rd(seconds=time()-t)))
    blocklen = len(hint)
    counter_sym = 0
    counter_rep = 0
    calc_syms = set()
    t = time()
    c = 0
    candidates = set(range(blocklen))
    while True:
        prev = c
        if candidates == set():
            c = prev + blocklen
        else:
            c = min(candidates)
        candidates.update(range(blocklen + prev, blocklen + c))
        i = 0
        for i in range(blocklen):
            next_index = c + blocklen - i
            if next_index in calc_syms:
                counter_rep += 1
                continue
            calc_syms.add(next_index)
            counter_sym += 1
            symbol = legendre(next_index, p) == 1
            # workaround a "set changed size during iteration" RuntimeError
            for k in candidates.copy():
                rel_index = k - c + i
                if rel_index <= blocklen and symbol != hint[-rel_index]:
                    candidates.remove(k)
            if c not in candidates:
                break
        if c in candidates:
            print("bruteforce execution time: " + fmt.format(rd(seconds=time() - t)))
            print("total number of symbol calulations = " + str(counter_sym))
            print("skipped redundant symbol calulations = " + str(counter_rep))
            assert randkey == c or not randkey, "wrong result, please debug"
            return c

# speedup? by trying to empty the whole set all the times (going back by 2*log b for example)


if __name__ == "__main__":
    # target prime size ~40 bits, first prize at ~80 bits
    bf_v(492876863)


