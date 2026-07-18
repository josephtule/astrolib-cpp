# Copyright 2025-2026 Joseph Tu Le
# SPDX-License-Identifier: Apache-2.0

import re
import sympy as sp

x, y, z = sp.symbols("x y z", real=True)

J2c, J3c, J4c, J5c, J6c = sp.symbols("J2 J3 J4 J5 J6", real=True)
mu, R = sp.symbols("mu R", real=True)

r_mag_s, r_mag2_s = sp.symbols("r_mag r_mag2", positive=True, real=True)

vars = sp.Matrix([x, y, z])

r_mag2 = x**2 + y**2 + z**2
r_mag = sp.sqrt(r_mag2)


def radius_subs(expr):
    expr = sp.factor(sp.simplify(expr))

    subs_map = {
        r_mag2: r_mag2_s,
        r_mag2 ** sp.Rational(1, 2): r_mag_s,
        r_mag2 ** sp.Rational(3, 2): r_mag_s**3,
        r_mag2 ** sp.Rational(5, 2): r_mag_s**5,
        r_mag2 ** sp.Rational(7, 2): r_mag_s**7,
        r_mag2 ** sp.Rational(9, 2): r_mag_s**9,
        r_mag2 ** sp.Rational(11, 2): r_mag_s**11,
        r_mag2 ** sp.Rational(13, 2): r_mag_s**13,
        r_mag2 ** sp.Rational(15, 2): r_mag_s**15,
        r_mag2 ** sp.Rational(17, 2): r_mag_s**17,
    }

    expr = expr.subs(subs_map)
    expr = sp.factor(sp.simplify(expr))

    return expr


def clean_matrix(M):
    return M.applyfunc(radius_subs)


def make_jacobian(a):
    return clean_matrix(a.jacobian(vars))


def cpp_code(expr):
    code = sp.ccode(sp.simplify(expr))

    code = code.replace("pow(", "std::pow(")

    # replace Jn when it's its own entity/character block
    code = re.sub(r"\bJ2\b", "J(2)", code)
    code = re.sub(r"\bJ3\b", "J(3)", code)
    code = re.sub(r"\bJ4\b", "J(4)", code)
    code = re.sub(r"\bJ5\b", "J(5)", code)
    code = re.sub(r"\bJ6\b", "J(6)", code)

    return code


def print_case_body(degree, A, indent="        "):
    name = f"j{degree}"

    print(f"{indent}mat3d G_{name};")

    for i in range(3):
        for j in range(3):
            expr = cpp_code(A[i, j])
            print(f"{indent}G_{name}({i}, {j}) = {expr};")

    print(f"{indent}G.block<3, 3>(3, 0) += G_{name};")


aJ2 = (
    -sp.Rational(3, 2)
    * J2c
    * mu
    * R**2
    / r_mag**5
    * sp.Matrix(
        [
            x * (1 - 5 * z**2 / r_mag2),
            y * (1 - 5 * z**2 / r_mag2),
            z * (3 - 5 * z**2 / r_mag2),
        ]
    )
)

aJ3 = (
    -sp.Rational(1, 2)
    * J3c
    * mu
    * R**3
    / r_mag**7
    * sp.Matrix(
        [
            5 * x * z * (3 - 7 * z**2 / r_mag2),
            5 * y * z * (3 - 7 * z**2 / r_mag2),
            -(3 * r_mag2 - 30 * z**2 + 35 * z**4 / r_mag2),
        ]
    )
)

aJ4 = (
    sp.Rational(5, 8)
    * J4c
    * mu
    * R**4
    / r_mag**7
    * sp.Matrix(
        [
            3 * x * (1 - 14 * z**2 / r_mag2 + 21 * z**4 / r_mag2**2),
            3 * y * (1 - 14 * z**2 / r_mag2 + 21 * z**4 / r_mag2**2),
            z * (15 - 70 * z**2 / r_mag2 + 63 * z**4 / r_mag2**2),
        ]
    )
)

aJ5 = (
    sp.Rational(3, 8)
    * J5c
    * mu
    * R**5
    / r_mag**9
    * sp.Matrix(
        [
            7 * x * (5 * z - 30 * z**3 / r_mag2 + 33 * z**5 / r_mag2**2),
            7 * y * (5 * z - 30 * z**3 / r_mag2 + 33 * z**5 / r_mag2**2),
            -(5 * r_mag2 - 105 * z**2 + 315 * z**4 / r_mag2 - 231 * z**6 / r_mag2**2),
        ]
    )
)

aJ6 = (
    -sp.Rational(7, 16)
    * J6c
    * mu
    * R**6
    / r_mag**9
    * sp.Matrix(
        [
            x
            * (
                5
                - 135 * z**2 / r_mag2
                + 495 * z**4 / r_mag2**2
                - 429 * z**6 / r_mag2**3
            ),
            y
            * (
                5
                - 135 * z**2 / r_mag2
                + 495 * z**4 / r_mag2**2
                - 429 * z**6 / r_mag2**3
            ),
            z
            * (
                35
                - 315 * z**2 / r_mag2
                + 693 * z**4 / r_mag2**2
                - 429 * z**6 / r_mag2**3
            ),
        ]
    )
)

jacobians = {
    2: make_jacobian(aJ2),
    3: make_jacobian(aJ3),
    4: make_jacobian(aJ4),
    5: make_jacobian(aJ5),
    6: make_jacobian(aJ6),
}


for degree, A in jacobians.items():
    is_symmetric = sp.simplify(A - A.T) == sp.zeros(3, 3)

    if not is_symmetric:
        raise RuntimeError(f"J{degree} Jacobian is not symmetric")

for degree in [6, 5, 4, 3, 2]:
    print_case_body(degree, jacobians[degree])
    print("\n")
