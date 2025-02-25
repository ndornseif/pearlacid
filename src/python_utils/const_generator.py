"""
Copyright 2025 N. Dornseif

Dual-licensed under Apache 2.0 and MIT terms.
Generates random constants.
"""

import random


def random_constant(bits: int, count: int = 1) -> str:
    constants = list()
    for i in range(count):
        constants.append(hex(random.getrandbits(bits)))
    return str(constants)


print(f"8 bits: {random_constant(8, 4)}")
print(f"16 bits: {random_constant(16, 4)}")
print(f"32 bits: {random_constant(32, 4)}")
print(f"64 bits: {random_constant(64, 4)}")
print(f"128 bits: {random_constant(128, 4)}")
