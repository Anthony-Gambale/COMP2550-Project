
"""
Testing indexing functions classically
"""

from operator import index
from pickletools import uint1
from qiskit.quantum_info.operators import Operator
import numpy as np
from math import log, ceil, factorial


log2 = lambda x: ceil(log(x,2))


def lehmer_encode(permutation):
    lehmer_code = []
    for i in range(len(permutation)):
        next_entry = 0
        for j in range(i, len(permutation)):
            if permutation[j] < permutation[i]:
                next_entry += 1
        lehmer_code.append(next_entry)
    return lehmer_code


def invert(permutation):
    inverse = [-1 for _ in range(len(permutation))]
    for idx,val in enumerate(permutation):
        inverse[val] = idx
    return inverse


def compose_permutations(permutations):
    def compose_helper(permutation1, permutation2):
        # apply permutation 2 to permutation 1's list
        assert len(permutation1) == len(permutation2), "Permutations must have same domain"
        rtn = [0 for _ in range(len(permutation1))]
        for idx, value in enumerate(permutation2):
            entry = permutation1[value] # grab entry
            rtn[idx] = entry            # move it to where it needs to be
        return rtn
    rtn = permutations[0]
    for perm in permutations[1:]: rtn = compose_helper(rtn, perm)
    return rtn


def rank_perm(permutation):
    n = len(permutation)
    lehmer_sequence = lehmer_encode(permutation)
    return sum([(n-1-idx)*lehmer_sequence[idx] for idx in range(n-1)])


def integer_to_binary(n):
    rtn = []
    while n > 0:
        rtn.append(n&1)
        n = n >> 1
    return rtn # lsb is index 0, msb is index len-1


def binary_to_integer(bits):
    rtn = 0
    acc = 1
    for val in bits:
        rtn += val * acc
        acc *= 2
    return rtn


def permutation_to_integer_repr(permutation):
    register_size = log2(len(permutation))
    rtn = []
    for value in permutation:
        next_bitstring = integer_to_binary(value)
        while len(next_bitstring) < register_size:
            next_bitstring.append(0)
        rtn = rtn + next_bitstring
    return binary_to_integer(rtn)


def rank_perm_as_integer(n, bitstring_length, perm_size):
    bitstring = integer_to_binary(n)
    while len(bitstring) < bitstring_length:
        bitstring.append(0) # pad with 0s in the most significant bits
    chunks = [bitstring[i:i+log2(perm_size)] for i in range(0,bitstring_length,log2(perm_size))]
    permutation = [binary_to_integer(chunk) for chunk in chunks]
    assert len(permutation) == perm_size, "permutation length is incorrect"
    # check that the permutation is injective (valid)
    for i in range(len(permutation)):
        if permutation[i] >= perm_size: return -1 # this is not a permutation!
        for j in range(i+1,len(permutation)):
            if permutation[i] == permutation[j]: return -1 # this is not a permutation!
    return rank_perm(permutation)


def indexing_unitary(n):
    nqubits = n * log2(n)
    rtn_matrix = np.zeros((2**nqubits,2**nqubits))
    current_invalids_row = factorial(n) # the number after the first n! indices will be where the non-permutations get indexed
    for permutation in range(0, len(rtn_matrix)):
        index_value = rank_perm_as_integer(permutation, nqubits, n)
        if index_value != -1:
            rtn_matrix[index_value,permutation] = 1
        else:
            rtn_matrix[current_invalids_row,permutation] = 1
            current_invalids_row += 1
    return Operator(rtn_matrix)


def unindexing_unitary(n):
    return indexing_unitary(n).transpose()


def tests():
    assert binary_to_integer(integer_to_binary(5)) == 5, "binary to and from integer is broken"
    assert lehmer_encode([1,5,0,6,3,4,2]) == [1, 4, 0, 3, 1, 1, 0], "lehmer encoding is broken"
    assert compose_permutations([invert([2,0,3,1]), [2,0,3,1], [0,1,2,3], invert([2,3,1,0]), [2,3,1,0]]) == [0,1,2,3], "permutation composition or inversion is broken"
    rp = rank_perm
    assert [rp([0,1,2]),rp([0,2,1]),rp([1,0,2]),rp([1,2,0]),rp([2,0,1]),rp([2,1,0])] == [0, 1, 2, 3, 4, 5], "permutation ranking is broken"
    rpi = rank_perm_as_integer
    assert [rpi(2, 2, 2), rpi(1, 2, 2)] == [0, 1], "permutation ranking as integers is broken"
    assert indexing_unitary(3).is_unitary(), "indexing unitary is not unitary"
    assert unindexing_unitary(3).is_unitary(), "unindexing unitary is not unitary"
    assert permutation_to_integer_repr([0,1,2]) == 36, "permutation_to_integer_repr is broken"
    print('tests pass')

if __name__ == "__main__":
    tests()