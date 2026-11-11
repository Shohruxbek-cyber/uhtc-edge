"""Exact thermoelastic reference solutions for V11 (validation/V11_plan.md), constant E, nu and alpha, stress-free at T_ref.

Hollow sphere a < r < b, spherically symmetric temperature T(r) (Timoshenko and Goodier, Theory of Elasticity, 3rd ed., Section 152):
    u(r) = (1 + nu)/(1 - nu) alpha / r^2 int_a^r (T - T_ref) s^2 ds + C1 r + C2 / r^2,
    sigma_r = E / ((1 + nu)(1 - 2 nu)) [(1 - nu) u' + 2 nu u / r - (1 + nu) alpha (T - T_ref)],
    sigma_t = E / ((1 + nu)(1 - 2 nu)) [nu u' + u / r - (1 + nu) alpha (T - T_ref)]     (hoop),
with C1, C2 from sigma_r(a) = sigma_r(b) = 0.
Thin free disc 0 < r < b, plane stress, radial temperature T(r) (same book, Section 150):
    u(r) = (1 + nu) alpha / r int_0^r (T - T_ref) s ds + C1 r,
    sigma_r = E / (1 - nu^2) [u' + nu u / r - (1 + nu) alpha (T - T_ref)],
    sigma_t = E / (1 - nu^2) [u / r + nu u' - (1 + nu) alpha (T - T_ref)],
with C1 from sigma_r(b) = 0. Integrals by adaptive quadrature; selfcheck() verifies equilibrium and the boundary conditions.
"""
import numpy as np
from scipy.integrate import quad


def sphere(T, a, b, E, nu, alpha, T_ref=0.0):
    k = (1 + nu) / (1 - nu) * alpha
    I = lambda r: quad(lambda s: (T(s) - T_ref) * s * s, a, r, epsabs=0, epsrel=1e-13, limit=200)[0]  # noqa: E731
    lam = E / ((1 + nu) * (1 - 2 * nu))

    def parts(r, C1, C2):
        u = k * I(r) / r ** 2 + C1 * r + C2 / r ** 2
        du = k * ((T(r) - T_ref) - 2 * I(r) / r ** 3) + C1 - 2 * C2 / r ** 3
        th = alpha * (T(r) - T_ref)
        sr = lam * ((1 - nu) * du + 2 * nu * u / r - (1 + nu) * th)
        st = lam * (nu * du + u / r - (1 + nu) * th)
        return sr, st
    # sigma_r is linear in (C1, C2): solve the 2 x 2 system for sigma_r(a) = sigma_r(b) = 0
    f0 = np.array([parts(a, 0, 0)[0], parts(b, 0, 0)[0]])
    J = np.array([[parts(a, 1, 0)[0] - f0[0], parts(a, 0, 1)[0] - f0[0]],
                  [parts(b, 1, 0)[0] - f0[1], parts(b, 0, 1)[0] - f0[1]]])
    C1, C2 = np.linalg.solve(J, -f0)
    return lambda r: parts(r, C1, C2)


def disc(T, b, E, nu, alpha, T_ref=0.0):
    I = lambda r: quad(lambda s: (T(s) - T_ref) * s, 0.0, r, epsabs=0, epsrel=1e-13, limit=200)[0]  # noqa: E731
    lam = E / (1 - nu ** 2)

    def parts(r, C1):
        u = (1 + nu) * alpha * I(r) / r + C1 * r
        du = (1 + nu) * alpha * ((T(r) - T_ref) - I(r) / r ** 2) + C1
        th = alpha * (T(r) - T_ref)
        return lam * (du + nu * u / r - (1 + nu) * th), lam * (u / r + nu * du - (1 + nu) * th)
    f0 = parts(b, 0)[0]
    C1 = -f0 / (parts(b, 1)[0] - f0)
    return lambda r: parts(r, C1)


def selfcheck():
    """Equilibrium (central differences) and boundary conditions of both solutions; returns the largest relative residual."""
    E, nu, al = 400e9, 0.14, 6e-6
    a, b = 0.5e-3, 1e-3
    T = lambda r: 1000.0 * (1 / r - 1 / b) / (1 / a - 1 / b)  # noqa: E731  steady conduction, T(a) = 1000, T(b) = 0
    s = sphere(T, a, b, E, nu, al)
    scale = E * al * 1000.0
    res = [abs(s(a)[0]) / scale, abs(s(b)[0]) / scale]
    for r in np.linspace(0.55e-3, 0.95e-3, 5):
        h = 1e-7
        dsr = (s(r + h)[0] - s(r - h)[0]) / (2 * h)
        sr, st = s(r)
        res.append(abs(dsr + 2 * (sr - st) / r) * r / scale)  # d sigma_r/dr + 2 (sigma_r - sigma_t)/r = 0
    Td = lambda r: 1000.0 * (r / b) ** 2  # noqa: E731
    d = disc(Td, b, E, nu, al)
    res.append(abs(d(b)[0]) / scale)
    for r in np.linspace(0.1e-3, 0.9e-3, 5):
        h = 1e-7
        dsr = (d(r + h)[0] - d(r - h)[0]) / (2 * h)
        sr, st = d(r)
        res.append(abs(dsr + (sr - st) / r) * r / scale)  # d sigma_r/dr + (sigma_r - sigma_t)/r = 0
    return max(res)


if __name__ == "__main__":
    print("largest relative residual of the exact solutions:", selfcheck())
    E, nu, al = 400e9, 0.14, 6e-6
    a, b = 0.5e-3, 1e-3
    s = sphere(lambda r: 1000.0 * (1 / r - 1 / b) / (1 / a - 1 / b), a, b, E, nu, al)
    print("sphere hoop stress at a and b (MPa):", s(a)[1] / 1e6, s(b)[1] / 1e6)
    d = disc(lambda r: 1000.0 * (r / b) ** 2, b, E, nu, al)
    print("disc hoop stress at 0+ and b (MPa):", d(1e-9)[1] / 1e6, d(b)[1] / 1e6, " closed form E alpha dT/4, -E alpha dT/2:",
          E * al * 1000 / 4 / 1e6, -E * al * 1000 / 2 / 1e6)
    # sphere: closed form of Timoshenko and Goodier for steady conduction, hoop stress at r (T_a at a, 0 at b)
    k = al * E * 1000.0 / (1 - nu) * a * b / (b ** 3 - a ** 3)
    for r in (a, b):
        print("sphere hoop closed form at", r, (k * (a + b - (b * b + a * b + a * a) / (2 * r) - a * a * b * b / (2 * r ** 3))) / 1e6)
