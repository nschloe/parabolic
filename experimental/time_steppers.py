# -*- coding: utf-8 -*-
#
import numpy


class Heun(object):
    '''
    Heun's method for :math:`u' = F(u)`.
    https://en.wikipedia.org/wiki/Heun's_method
    '''
    order = 2.0

    def __init__(self, problem):
        self.problem = problem

        # alpha = 0.5
        # alpha = 2.0 / 3.0
        alpha = 1.0

        self.tableau = {
            'A': [[0.0,   0.0], [alpha, 0.0]],
            'b': [1.0 - 1.0 / (2 * alpha), 1.0 / (2 * alpha)],
            'c': [0.0, alpha],
            }
        return

    def step(self, u0, t, dt):
        return _runge_kutta_step(self.problem, self.tableau, u0, t, dt)


# def rk4_step(
#         V,
#         F,
#         u0,
#         t, dt,
#         sympy_dirichlet_bcs=[],
#         tol=1.0e-10,
#         verbose=True
#         ):
#     '''Classical RK4.
#     '''
#     c = [0.0, 0.5, 0.5, 1.0]
#     A = [[0.0, 0.0, 0.0, 0.0],
#          [0.5, 0.0, 0.0, 0.0],
#          [0.0, 0.5, 0.0, 0.0],
#          [0.0, 0.0, 1.0, 0.0]]
#     b = [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0]
#
#     return runge_kutta_step(
#             A, b, c,
#             V, F, u0, t, dt,
#             sympy_dirichlet_bcs=sympy_dirichlet_bcs,
#             tol=tol,
#             verbose=verbose
#             )
#
#
# def rkf_step(
#         V,
#         F,
#         u0,
#         t, dt,
#         sympy_dirichlet_bcs=[],
#         tol=1.0e-10,
#         verbose=True
#         ):
#     '''Runge--Kutta--Fehlberg method.
#     '''
#     c = [0.0, 0.25, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5]
#     A = [[0.0,         0.0,         0.0,         0.0,         0.0,    0.0],
#          [0.25,        0.0,         0.0,         0.0,         0.0,    0.0],
#          [3./32,       9./32,       0.0,         0.0,         0.0,    0.0],
#          [1932./2197, -7200./2197,  7296./2197,  0.0,         0.0,    0.0],
#          [439./216,   -8.,          3680./513,  -845./4104,   0.0,    0.0],
#          [-8./27,      2.,         -3544./2565,  1859./4104, -11./40, 0.0]]
#     # b = [25./216, 0.0, 1408./2565, 2197./4104, -1./5,  0.0] # 4th order
#     # 5th order
#     b = [16./135, 0.0, 6656./12825,  28561./56430, -9./50, 2./55]
#
#     return runge_kutta_step(
#             A, b, c,
#             V, F, u0, t, dt,
#             sympy_dirichlet_bcs=sympy_dirichlet_bcs,
#             tol=tol,
#             verbose=verbose
#             )

def _runge_kutta_step(
        problem, tableau, u0, t, dt
        ):
    A = numpy.array(tableau['A'])
    b = tableau['b']
    c = tableau['c']

    # Make sure that the scheme is strictly lower-triangular.
    s = len(tableau['b'])
    # Can't handle implicit methods yet.
    assert numpy.all(abs(A[numpy.triu_indices(s)]) < 1.0e-15)

    # # For the boundary values, see
    # #
    # #   Intermediate Boundary Conditions for Runge-Kutta Time Integration of
    # #   Initial-Boundary Value Problems,
    # #   D. Pathria,
    # #   <http://www.math.uh.edu/~hjm/june1995/p00379-p00388.pdf>.
    # #
    # tt = sympy.symbols('t')
    # BCS = []
    # # Get boundary conditions and their derivatives.
    # for k in range(2):
    #     BCS.append([])
    #     for boundary, expr in sympy_dirichlet_bcs:
    #         # Form k-th derivative.
    #         DexprDt = sympy.diff(expr, tt, k)
    #         # TODO set degree of expression
    #         BCS[-1].append(
    #             DirichletBC(
    #                 V,
    #                 Expression(sympy.printing.ccode(DexprDt), t=t + dt),
    #                 boundary
    #                 )
    #             )

    # Compute the stage values.
    k = [u0.copy() for i in range(s)]
    for i in range(s):
        U = u0.copy()
        for j in range(i):
            if A[i][j] != 0.0:
                U.vector()[:] += dt * A[i][j] * k[j].vector()

        L = problem.eval_alpha_M_beta_F(0.0, 1.0, U, t + c[i]*dt)
        # TODO boundary conditions!
        # for g in BCS[1]:
        #     g.t = t + c[i] * dt
        k[i].assign(problem.solve_alpha_M_beta_F(1.0, 0.0, L, t + c[i]*dt))

    # Put it all together.
    U = u0.copy()
    for i in range(s):
        U.vector()[:] += dt * b[i] * k[i].vector()

    # TODO boundary conditions
    # for g in BCS[0]:
    #     g.t = t + dt
    theta = problem.solve_alpha_M_beta_F(1.0, 0.0, U, t+dt)
    return theta
