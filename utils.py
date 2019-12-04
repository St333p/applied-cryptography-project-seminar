from bitarray import bitarray
from typing import *
from time import time
from itertools import repeat

setup = "from utils import cyclic_bitarray, bitarray\nfrom random import randint\na = bitarray([randint(0,1) for _ in range(1000000)])\nb = bitarray([randint(0,1) for _ in range(1000000)])\nc = cyclic_bitarray(1000000, random=True)"

a = 10

# consider switching away from bitarrays, alternatives: intbitset
class cyclic_bitarray:
    """
    this offers a x1000 speedup in shifting sequence compared to the 
    naive same operation on bitarrays
    """
    def __init__(self, length: int, default: bool=False, random=False):
        self.length = length
        self.default = default
        if random:
            self.array = bitarray(length)
        else:
            self.array = bitarray([self.default] * self.length)
        self.index = 0

    def get(self, index: int) -> bool:
        return self.array[(self.index + index) % self.length]
    
    def set(self, a: int, replacement):
        # check consistency of inputs?
        if a < 0 or a >= self.length:
            raise Exception("Out of bounds")
        target = (a + self.index) % self.length
        if replacement in (0,1):
            self.array[target] = replacement
            return
        l = len(replacement)
        if a + l > self.length:
            raise Exception("Substring out of bounds")
        if target + l <= self.length:
            self.array[target: target + l] = bitarray(replacement)
        else:
            exceeding = (target + l) % self.length 
            self.array[target:] = bitarray(replacement[:-exceeding])
            self.array[:exceeding] = bitarray(replacement[-exceeding:])
        assert len(self.array) == self.length, "something went wrong"

    def shift(self, steps: int):
        assert steps <= self.length, "cannot shift more than bitarray length"
        self.set(0, [self.default]*steps)
        self.index = (self.index + steps) % self.length
    
    def bitwise_and_v2(self, stream: bitarray, l: int):
        # TODO: this is the speed bottleneck, improve it!
        t = time()
        target = self.index
        # l = len(stream)
        #print(time() - t)
        t = time()
        if target + l <= self.length:
            self.array[target: target + l] &= stream
        else:
            exceeding = (target + l) % self.length
            self.array[target:] &= stream[:-exceeding]
            self.array[:exceeding] &= stream[-exceeding:]
       # print(time() - t)
    
    def bitwise_and(self, stream: bitarray, l: int):
        # TODO: this is the speed bottleneck, improve it!
        t = time()
        target = self.index
        # l = len(stream)
        #print(time() - t)
        t = time()
        if target + l <= self.length:
            s = bitarray(repeat(False,target)) + stream + bitarray(repeat(False, self.length - target - l))
        else:
            exceeding = (target + l) % self.length
            s = stream[-exceeding:] + bitarray(repeat(False, self.length - l)) + stream[:-exceeding]
        self.array &= s

    def first(self, value: bool)-> int:
        '''
        return the index of the first occurence of value in the array
        '''
        try:
            return self.array.index(value, self.index) - self.index
        except ValueError:
            try:
                return self.array.index(value, 0, self.index) - self.index + self.length
            except ValueError:
                return self.length
        # this is an efficient implementation of the following desired syntax
        #for i in range(self.length):
        #    if self.get(i) == value:
        #        return i
        #return self.length 
    
    def any(self, a, b):
        return self.array[self.index + a: self.index + b].any() 

    def __str__(self):
        return 'cyclic_' + str(self.array[self.index:] + self.array[:self.index])

