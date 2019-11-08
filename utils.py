from bitarray import bitarray
from typing import *

class cyclic_bitarray:
    
    def __init__(self, length: int, default: bool=False):
        self.length = length
        self.default = default
        self.array = bitarray([self.default] * self.length)
        self.index = 0

    def set(self, a: int, replacement):
        # check consistency of inputs?
        if a < 0 or a >= self.length:
            raise Exception("Out of bounds")
        target = (a + self.index) % self.length
        if replacement in (0,1):
            self.array[target] = replacement
            return
        l = len(replacement)
        if target + l < self.length:
            self.array[target: target + l] = bitarray(replacement)
        else:
            exceeding = (target + l) % self.length
            self.array[target:] = bitarray(replacement[:-exceeding])
            self.array[:exceeding] = bitarray(replacement[-exceeding:])
        assert len(self.array) == self.length, "something went wrong"

    def shift(self, steps: int):
        self.set(self.index, [self.default]*steps)
        self.index += steps
    
    def __str__(self):
        return 'cyclic_' + str(self.array)

