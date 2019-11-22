import numpy as np
from typing import *
from sympy.ntheory import isprime, factorint
from sympy import randprime, legendre_symbol, jacobi_symbol
from bitarray import bitarray
import gmpy2
from random import randint
from time import time
from itertools import repeat
from dateutil.relativedelta import relativedelta as rd
from utils import cyclic_bitarray
import argparse
fmt = '{0.minutes} m, {0.seconds} s'

def legendre(a: int, p: int) -> int:
    # assert isprime(p) and p > 2, "ERROR: p has to be an odd prime"
    # return legendre_symbol(a, p) == 1
    return gmpy2.legendre(a,p) == 1

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

def bf_v2(p: int, hint, start_from=0, confidence_bits=100)-> int:
    """
    improvement: x10 speedup over v1
    not all symbols need to be calculated
    """
    blocklen = len(hint)
    counter_sym = 0
    counter_rep = 0
    calc_syms = cyclic_bitarray(blocklen, False)
    t = time()
    c = start_from
    candidates = set(range(c, c + blocklen))
    print("bruteforcing with version 2...")
    med_len = 0
    avg_inner_time = 0
    while True:
        prev = c
        if candidates == set():
            c = prev + blocklen
        else:
            c = min(candidates)
        if time() - t > 20:
            t = time()
            print(c)
        candidates.update(range(blocklen + prev, blocklen + c))
        for i in range(confidence_bits):
            next_index = c + blocklen - i
            counter_sym += 1
            med_len = med_len + len(candidates)
            symbol = legendre(next_index, p)
            t1 = time()
            # workaround a "set changed size during iteration" RuntimeError
            for k in sorted(list(candidates)):
                rel_index = k - c + i
                if rel_index == 0:
                    continue
                if rel_index > blocklen:
                    break 
                if symbol != hint[-rel_index]:
                    candidates.remove(k)
            avg_inner_time += time() - t1
            if c not in candidates:
                break
        if c in candidates:
            print("total number of symbol calulations = " + str(counter_sym))
            print("average length of candidates: " + str(med_len / counter_sym))
            print("average inner cycle time: " + str(avg_inner_time / counter_sym))
            return c

# bitwise operations on int are significantly faster than on bitarrays: shall we switch?
def bf_v3(p: int, hint, start_from=0, confidence_bits=100)-> int:
    """
    potential improvement: optimize internal data structures and memory access
    as of now ~2x speedup over v2 for proper values of the hint length, approx 1000
    """
    blocklen = len(hint)
    counter_sym = 0
    candidates = cyclic_bitarray(blocklen, True)
    t = time()
    c = start_from
    inner_time = 0
    shifts_time = 0
    first_time = 0
    print("bruteforcing with version 3...")
    while True:
        t1 = time()
        offset = candidates.first(True)
        first_time += time() - t1
        # print("offset " + str(offset))
        c += offset
        if time() - t > 20:
            t = time()
            print(c)
        candidates.shift(offset)
        for i in range(confidence_bits):
            to_calc_offset = blocklen - i - 1
            counter_sym += 1
            symbol = legendre(c + to_calc_offset, p)
            t1 = time()
            # mod represents a vector containing the values of hint[j] == symbol
            # ordered compatibly with candidates
            if symbol:
                mod = hint[-(i + 1): -(blocklen + 1): -1]
            else:
                mod = ~hint[-(i + 1): -(blocklen + 1): -1] 
            inner_time += time() - t1
            t2 = time()
            candidates.bitwise_and(mod, blocklen - i)
            shifts_time += time() - t2
            if not candidates.get(0):
                 break
        if candidates.get(0):
            print("total number of symbol calulations = " + str(counter_sym))
            print("inner cycle summed time = " + fmt.format(rd(seconds=inner_time)))
            print("large internal cycle summed time = " + fmt.format(rd(seconds=shifts_time)))
            return c

# TODO: pass bf version as parameter
def bruteforce(security_bits, stream_length, keyspace_bits=None, key=None, p=None):
    """ 
    bruteforce a legendre prng stream generated with a prime p of given security bits
    if key is not specified it's generated randomly in {0, ..., p-1}
    uses the latest version of the bf function if none specified
    """
    if not p:
        p = randprime(2**security_bits, 2**(security_bits + 1))
    if not key:
        key = randint(p // 2, p)
    print('key: ' + str(key) + ", p: " + str(p))
    start_from = 0
    if keyspace_bits:
        start_from = max(0, key - 2**keyspace_bits)
        print('start_from: ' + str(start_from))
    t = time()
    hint = prng(key, 0, p, stream_length)
    print("hint derivation time: " + fmt.format(rd(seconds=time()-t)))
    t = time()
    result = bf_v3(p, hint, start_from)
    assert key == result, "wrong result, please debug"
    print("bruteforce v3 execution time: " + fmt.format(rd(seconds=time() - t)))
    t = time()
    result = bf_v2(p, hint, start_from)
    assert key == result, "wrong result, please debug"
    print("bruteforce v2 execution time: " + fmt.format(rd(seconds=time() - t)))
    return result

if __name__ == "__main__":
    # target prime size ~40 bits, first prize at ~80 bits
    parser = argparse.ArgumentParser(description='Bruteforce the legendre prng.')
    parser.add_argument('security_bits', type=int,
            help='bitlength of the prime number and the key')
    parser.add_argument('stream_length', type=int,
            help='bitlength of the stream passed as hint for the bruteforce')
    parser.add_argument('keyspace_bits', type=int,
            help='bitlength of the keyspace to search from')
    args = parser.parse_args()
    #could try to build up a test case to try and use
    bruteforce(args.security_bits, args.stream_length, args.keyspace_bits)


