"""
Copyright 2025 N. Dornseif

Dual-licensed under Apache 2.0 and MIT terms.
Generates matrix test cases,
for finding the rank of a binary matrix over GF(2).
"""

from pathlib import Path
import random

import numpy as np
import sympy
from sympy.polys.domains import GF
from sympy.polys.matrices import DomainMatrix

# All matrices are square for this usecase.
MATRIX_SIZE = 32
RESULT_FILE = "test_matrices.txt"

RANDOM_REDUCED_TEST_CASES = [256, 128, 128, 64, 64, 64, 64]
SPARSE_REPLACED_TEST_CASES = [1, 64, 64, 64, 64]
DENSE_REPLACED_TEST_CASES = [1, 64, 64, 64, 64]


def get_gf2_matrix_rank(matrix: np.ndarray) -> int:
    """Get the rank of a binary matrix over GF(2)."""
    matrix_gf2 = DomainMatrix.from_Matrix(sympy.Matrix(matrix.tolist())).convert_to(GF(2))
    return matrix_gf2.rank()


def generate_random_binary_matrix() -> np.ndarray:
    """Generate a random binary matrix of size MATRIX_SIZE x MATRIX_SIZE."""
    return np.random.randint(2, size=(MATRIX_SIZE, MATRIX_SIZE))


def generate_random_matrix_row() -> np.ndarray:
    """Generate a MATRIX_SIZEd row filled with random binary data."""
    return np.random.randint(2, size=MATRIX_SIZE)


def generate_full_binary_matrix() -> np.ndarray:
    """Generate a binary matrix of size MATRIX_SIZE x MATRIX_SIZE where all elements are one."""
    return np.full((MATRIX_SIZE, MATRIX_SIZE), 1, dtype=int)


def generate_empty_binary_matrix() -> np.ndarray:
    """Generate a binary matrix of size MATRIX_SIZE x MATRIX_SIZE where all elements are zero."""
    return np.zeros((MATRIX_SIZE, MATRIX_SIZE), dtype=int)


def randomized_rank_reduction(input_matrix: np.ndarray, steps: int) -> np.ndarray:
    """Reduce the average rank of a matrix by replacing 'steps' rows
    with randomly selected different rows."""
    matrix = input_matrix.copy()
    for i in range(min(steps, len(matrix))):
        matrix[i] = random.choice(matrix)
    return matrix


def insert_randomized_rank(input_matrix: np.ndarray, steps: int) -> np.ndarray:
    """Randomly choose a row 'steps' times and replace it with a randomized row."""
    matrix = input_matrix.copy()
    for i in range(steps):
        matrix[random.randint(0, len(matrix) - 1)] = generate_random_matrix_row()
    return matrix


def generate_sparse_test_matrix(replacement_steps: int) -> np.ndarray:
    matrix = generate_empty_binary_matrix()
    matrix = insert_randomized_rank(matrix, replacement_steps)
    return matrix


def generate_dense_test_matrix(replacement_steps: int) -> np.ndarray:
    matrix = generate_full_binary_matrix()
    matrix = insert_randomized_rank(matrix, replacement_steps)
    return matrix


def generate_intemediate_test_matrix(replacement_steps: int) -> np.ndarray:
    matrix = generate_random_binary_matrix()
    matrix = randomized_rank_reduction(matrix, replacement_steps)
    return matrix


def repr_matrix_as_intlist(matrix: np.ndarray) -> list[int]:
    """Represent each row in a binary matrix as an integer."""
    return [int("".join(map(str, row)), 2) for row in matrix]


def format_test_matrix(matrix: np.ndarray) -> str:
    """Format a test matrix as integer list and determine its rank over GF(2).
    Prepare both for inclusion in rust code."""
    int_rows = repr_matrix_as_intlist(matrix)
    rank = get_gf2_matrix_rank(matrix)
    return f"TestMatrix {{matrix:{int_rows},rank:{rank}}},\n"


with open(Path(__file__).with_name(RESULT_FILE),"w") as fd:
    # Identity matrix
    fd.write(format_test_matrix(np.eye(MATRIX_SIZE, dtype=int)))

    for reduction, count in enumerate(RANDOM_REDUCED_TEST_CASES):
        print(f"Generating {count} random matrices with {reduction} reduction steps.")
        for i in range(count):
            fd.write(format_test_matrix(generate_intemediate_test_matrix(reduction)))

    for replacement, count in enumerate(SPARSE_REPLACED_TEST_CASES):
        print(f"Generating {count} sparse matrices with {replacement} replacement steps.")
        for i in range(count):
            fd.write(format_test_matrix(generate_sparse_test_matrix(replacement)))

    for replacement, count in enumerate(DENSE_REPLACED_TEST_CASES):
        print(f"Generating {count} dense matrices with {replacement} replacement steps.")
        for i in range(count):
            fd.write(format_test_matrix(generate_dense_test_matrix(replacement)))
