import numpy as np
from typing import *
from sympy.ntheory import isprime, factorint
from sympy import randprime, legendre_symbol, jacobi_symbol
from bitarray import bitarray
from random import randint
from time import time
from dateutil.relativedelta import relativedelta as rd
from utils import cyclic_bitarray
fmt = '{0.minutes} m, {0.seconds} s'

def legendre(a: int, p: int) -> int:
    # assert isprime(p) and p > 2, "ERROR: p has to be an odd prime"
    #return _legendre(a, p) == 1
    return jacobi_symbol(a,p) == 1

def _legendre(a: int, p: int) -> int:
    """
    computes the legendre symbol of a mod p where p is an odd prime,
    algorithm taken from https://martin-thoma.com/how-to-calculate-the-legendre-symbol/
    """
    # try to speed this up using jacobi symbol or fast expenentiation 
    # then try to implemented in C using gmp / boost (it contains legendre and jacobi symbols)
    if a > p or a < 0:
        return _legendre(a % p, p)
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
        return np.prod(np.array([_legendre(f,p) ** facs[f] for f in facs.keys()]))
    else:
        if p % 4 == 1 or a % 4 == 1:
            return _legendre(p,a)
        else:
            return - _legendre(p,a)

def prng(key: int, seed: int, p: int, length: int) -> bitarray:
    # assert that n is not too large (10 * log2^10(p) ?!?)
    result = bitarray()
    for i in range(length):
        result.append(legendre(seed + key + i, p))
    return result

def bf_v0(p: int, hint=None, rand=False)-> int:
    """
    naive bruteforce:
    given the pseudorandom stream f_k(0), ... ,f_k(t) as hint
    returns the key thanks to a bruteforce search
    """
    k = 0
    t = time()
    while True:
        if hint == prng(k, 0, p, len(hint)):
            print("execution time: " + str(time() - t))
            return k
        k += 1

def bf_v1(p: int, hint, par=0.1)-> int:
    """
    improvement:
    avoid calculating the same legendre symbol too many times
    """
    k = 0
    t = time()
    print("bruteforcing with version 1...")
    while True:
        stream = prng(k, 0, p, len(hint))
        if hint == stream:
            print("execution time: " + str(time() - t))
            return k % p
        candidates = []
        for i in range(int(len(hint) / par)):
            if hint[i:] == stream[:-i]:
                candidates.append(k - i)
            if hint[:-i] == stream[i:]:
                candidates.append(k + i)
        for key in candidates:
            if hint == prng(key, 0, p, len(hint)):
                print("execution time: " + str(time() - t))
                return key % p
        k += 2 * int(len(hint) / par)

def bf_v2(p: int, hint)-> int:
    """
    improvement: x10 speedup over v1
    not all symbols need to be calculated
    """
    blocklen = len(hint)
    counter_sym = 0
    counter_rep = 0
    calc_syms = cyclic_bitarray(blocklen, False)
    t = time()
    c = 0
    candidates = set(range(blocklen))
    print("bruteforcing with version 2...")
    while True:
        prev = c
        if candidates == set():
            c = prev + blocklen
        else:
            c = min(candidates)
        candidates.update(range(blocklen + prev, blocklen + c))
        calc_syms.shift(c - prev)
        for i in range(blocklen):
            next_index = c + blocklen - i
            if calc_syms.get(next_index - c - 1):
                counter_rep += 1
                continue
            calc_syms.set(next_index - c - 1, True)
            counter_sym += 1
            symbol = legendre(next_index, p)
            # workaround a "set changed size during iteration" RuntimeError
            for k in candidates.copy():
                rel_index = k - c + i
                if rel_index > blocklen:
                    break
                if symbol != hint[-rel_index]:
                    candidates.remove(k)
            if c not in candidates:
                break
        if c in candidates:
            print("total number of symbol calulations = " + str(counter_sym))
            print("skipped redundant symbol calulations = " + str(counter_rep))
            return c

def bf_v3(p: int, hint)-> int:
    """
    potential improvement: optimize internal data structures and memory access
    as of now waay worse
    """
    blocklen = len(hint)
    counter_sym = 0
    counter_rep = 0
    calc_syms = cyclic_bitarray(blocklen, False)
    candidates = cyclic_bitarray(blocklen, True)
    t = time()
    c = 0
    print("bruteforcing with version 3...")
    while True:
        offset = candidates.first(True)
        # print("offset " + str(offset))
        c += offset
        # print(c)
        candidates.shift(offset)
        calc_syms.shift(offset)
        for i in range(blocklen):
            to_calc_offset = blocklen - i - 1
            if calc_syms.get(to_calc_offset):
                counter_rep += 1
                continue
            calc_syms.set(to_calc_offset, True)
            counter_sym += 1
            symbol = legendre_symbol(c + to_calc_offset, p)
            for k in range(blocklen - i):
                if not candidates.get(k):
                    continue
                rel_index = k + i + 1
                if symbol != hint[-rel_index]:
                    candidates.set(k, False)
            if not candidates.get(0):
                break
        if candidates.get(0):
            print("total number of symbol calulations = " + str(counter_sym))
            print("skipped redundant symbol calulations = " + str(counter_rep))
            return c

def bruteforce(security_bits, key=None, version=None, stream_length=10**5):
    """
    bruteforce a legendre prng stream generated with a prime p of given security bits
    if key is not specified it's generated randomly in {0, ..., p-1}
    uses the latest version of the bf function if none specified
    """
    p = randprime(2**security_bits, 2**(security_bits + 1))
    if not key:
        key = randint(0, p)
    print('key: ' + str(key) + ", p: " + str(p))
    t = time()
    hint = prng(key, 0, p, stream_length)
    print("hint derivation time: " + fmt.format(rd(seconds=time()-t)))
    t = time()
    result = bf_v2(p, hint)
    print("bruteforce execution time: " + fmt.format(rd(seconds=time() - t)))
    print("result: " + str(result))
    assert key == result, "wrong result, please debug"
    return result

if __name__ == "__main__":
    # target prime size ~40 bits, first prize at ~80 bits
    bruteforce(23, stream_length=10**5)


