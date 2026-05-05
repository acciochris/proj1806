import math


def mmul(A, B):
    return [
        [A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
        [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]],
    ]


def dot(v1, v2):
    return sum([a * b for a, b in zip(v1, v2)])


def transp(A):
    return [[A[0][0], A[1][0]], [A[0][1], A[1][1]]]


def norm(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1])


def scale(v, c):
    if isinstance(v, (int, float)):
        return v * c
    return [scale(x, c) for x in v]


def add(v1, v2):
    if isinstance(v1, (int, float)):
        return v1 + v2
    return [add(x, y) for x, y in zip(v1, v2)]


def det(A):
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def inv(A):
    return scale([[A[1][1], -A[0][1]], [-A[1][0], A[0][0]]], 1 / det(A))


def pmat(A):
    for row in A:
        for c in row:
            print(f"{c:.8f}", end=" ")
        print()


def peig1(A):
    print(f"${A[0][0]:.8f}$", end=", ")


def peig2(A):
    print(f"${A[1][1]:.8f}$", end=", ")


def qr_decomp(A):
    a1, a2 = transp(A)
    q1 = scale(a1, 1 / norm(a1))
    v2 = add(a2, scale(q1, -dot(q1, a2)))
    q2 = scale(v2, 1 / norm(v2))
    q = transp([q1, q2])
    r = mmul(inv(q), A)
    return q, r


def qr(A, k=10):
    for i in range(k):
        # print("Iteration", i)
        peig2(A)
        # print()
        q, r = qr_decomp(A)
        A = mmul(r, q)


if __name__ == "__main__":
    qr([[2, 1], [1, 2]], k=7)
