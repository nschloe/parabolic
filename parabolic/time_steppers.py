# -*- coding: utf-8 -*-
#
'''
Routines for one time-stepping of the general equation

.. math::
    \\frac{du}{dt} = F(u).

'''


class Dummy():
    '''
    Dummy method for :math:`u' = F(u)`.
    '''
    name = 'Dummy'
    order = 0.0

    def __init__(self, problem):
        self.problem = problem
        return

    def step(self, u0, t, dt):
        # (u1 - u0) / dt = 0
        # u1 = u0
        return u0


class ExplicitEuler():
    '''
    Explicit Euler method for :math:`u' = F(u)`.
    '''
    name = 'Explicit Euler'
    order = 1.0

    def __init__(self, problem):
        self.problem = problem
        return

    def step(self, u0, t, dt):
        # (u{k+1} - u{k}) / dt = F(u{k}, t)
        # u{k+1} = u{k} + dt * F(u{k}, t)
        b = self.problem.eval_alpha_M_beta_F(1.0, dt, u0, t)
        u1 = self.problem.solve_alpha_M_beta_F(1.0, 0.0, b, t+dt)
        return u1


class ImplicitEuler():
    '''
    Implicit Euler method for :math:`u' = F(u)`.
    '''
    name = 'Implicit Euler'
    order = 1.0

    def __init__(self, problem):
        self.problem = problem
        return

    def step(self, u0, t, dt):
        # (u{k+1} - u{k}) / dt = F(u{k+1}, t+dt)
        # u{k+1} - dt * F(u{k}, t+dt) = u{k}
        b = self.problem.eval_alpha_M_beta_F(1.0, 0.0, u0, t)
        u1 = self.problem.solve_alpha_M_beta_F(1.0, -dt, b, t+dt)
        return u1


class Trapezoidal():
    '''
    Trapezoidal method for :math:`u' = F(u)`. (Known as Crank-Nicolson in the
    ODE context.)
    '''
    name = 'Trapezoidal'
    order = 2.0

    def __init__(self, problem):
        self.problem = problem
        return

    def step(self, u0, t, dt):
        # (u{k+1} - u{k}) / dt = 1/2 * (F(u{k+1}, t+dt) + F(u{k}, t))
        # u{k+1} - dt/2 * F(u{k+1}, t+dt) = u{k} + dt/2 * F(u{k}, t)
        b = self.problem.eval_alpha_M_beta_F(1.0, 0.5*dt, u0, t)
        u1 = self.problem.solve_alpha_M_beta_F(1.0, -0.5*dt, b, t+dt)
        return u1


# def heun_step(
#         V,
#         F,
#         u0,
#         t, dt,
#         sympy_dirichlet_bcs,
#         tol=1.0e-10,
#         verbose=True
#         ):
#     # Heun & variants.
#     # alpha = 0.5
#     alpha = 2.0 / 3.0
#     # alpha = 1.0
#     c = [0.0, alpha]
#     A = [[0.0,   0.0],
#          [alpha, 0.0]]
#     b = [1.0 - 1.0 / (2 * alpha), 1.0 / (2 * alpha)]
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
#
#
# def runge_kutta_step(
#         A, b, c,
#         V,
#         F,
#         u0,
#         t, dt,
#         sympy_dirichlet_bcs=[],
#         tol=1.0e-10,
#         verbose=True
#         ):
#     # Make sure that the scheme is strictly upper-triangular.
#     s = len(b)
#     A = numpy.array(A)
#     # Can't handle fully implicit methods.
#     assert numpy.all(abs(A[numpy.triu_indices(s)]) < 1.0e-15)
#
#     u = TrialFunction(V)
#     v = TestFunction(V)
#
#     solver_params = {
#             'linear_solver': 'iterative',
#             'symmetric': True,
#             'preconditioner': 'hypre_amg',
#             'krylov_solver': {
#                 'relative_tolerance': tol,
#                 'absolute_tolerance': 0.0,
#                 'maximum_iterations': 100,
#                 'monitor_convergence': verbose
#                 }
#             }
#
#     # For the boundary values, see
#     #
#     #   Intermediate Boundary Conditions for Runge-Kutta Time Integration of
#     #   Initial-Boundary Value Problems,
#     #   D. Pathria,
#     #   <http://www.math.uh.edu/~hjm/june1995/p00379-p00388.pdf>.
#     #
#     tt = sympy.symbols('t')
#     BCS = []
#     # Get boundary conditions and their derivatives.
#     for k in range(2):
#         BCS.append([])
#         for boundary, expr in sympy_dirichlet_bcs:
#             # Form k-th derivative.
#             DexprDt = sympy.diff(expr, tt, k)
#             # TODO set degree of expression
#             BCS[-1].append(
#                 DirichletBC(
#                     V,
#                     Expression(sympy.printing.ccode(DexprDt), t=t + dt),
#                     boundary
#                     )
#                 )
#
#     # Use the Constant() syntax to avoid compiling separate expressions for
#     # different values of dt.
#     ddt = Constant(dt)
#
#     # Compute the stage values.
#     k = []
#     for i in range(s):
#         U = u0
#         for j in range(i):
#             U += ddt * A[i][j] * k[j]
#         L = F(t + c[i] * dt, U, v)
#         k.append(Function(V))
#         # Using this bc is somewhat random.
#         # TODO come up with something better here.
#         for g in BCS[1]:
#             g.t = t + c[i] * dt
#         solve(
#             u * v * dx == L, k[i],
#             bcs=BCS[1],
#             solver_parameters=solver_params
#             )
#         # plot(k[-1])
#         # interactive()
#
#     # Put it all together.
#     U = u0
#     for i in range(s):
#         U += b[i] * k[i]
#     theta = Function(V)
#     for g in BCS[0]:
#         g.t = t + dt
#     solve(
#         u * v * dx == (u0 + ddt * U) * v * dx,
#         theta,
#         bcs=BCS[0],
#         solver_parameters=solver_params
#         )
#
#     return theta
