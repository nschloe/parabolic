# -*- coding: utf-8 -*-
#
import parabolic

# pylint: disable=import-error
from dolfin import (
    FunctionSpace, DirichletBC, Function, grad, dx, dot, UnitSquareMesh,
    TrialFunction, TestFunction, assemble, Constant, XDMFFile, KrylovSolver,
    solve
    )


def test_heat_equation_fenics():
    # Define problem
    class Heat(object):
        '''
        u' = \\Delta u + f
        '''
        def __init__(self, V):
            self.V = V
            u = TrialFunction(V)
            v = TestFunction(V)
            self.M = assemble(u * v * dx)
            self.A = assemble(-dot(grad(u), grad(v)) * dx)
            self.b = assemble(1.0 * v * dx)
            self.bcs = DirichletBC(self.V, 0.0, 'on_boundary')
            return

        # pylint: disable=unused-argument
        def eval_alpha_M_beta_F(self, alpha, beta, u, t):
            # Evaluate  alpha * M * u + beta * F(u, t).
            uvec = u.vector()
            return alpha * (self.M * uvec) + beta * (self.A * uvec + self.b)

        def solve_alpha_M_beta_F(self, alpha, beta, b, t):
            # Solve  alpha * M * u + beta * F(u, t) = b  for u.
            A = alpha * self.M + beta * self.A

            rhs = b - beta * self.b
            self.bcs.apply(A, rhs)

            solver = KrylovSolver('gmres', 'ilu')
            solver.parameters['relative_tolerance'] = 1.0e-13
            solver.parameters['absolute_tolerance'] = 0.0
            solver.parameters['maximum_iterations'] = 100
            solver.parameters['monitor_convergence'] = True
            solver.set_operator(A)

            u = Function(self.V)
            solver.solve(u.vector(), rhs)
            return u

    # create initial guess
    mesh = UnitSquareMesh(20, 20, 'crossed')
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    u0 = Function(V)
    solve(u*v*dx == Constant(0.0)*v*dx, u0)

    u1 = Function(V)
    u1.assign(u0)

    # create time stepper
    # stepper = parabolic.Dummy(Heat(V))
    # stepper = parabolic.ExplicitEuler(Heat(V))
    stepper = parabolic.ImplicitEuler(Heat(V))
    # stepper = parabolic.Trapezoidal(Heat(V))

    # step
    t = 0.0
    dt = 1.0e-3
    with XDMFFile('heat.xdmf') as xf:
        xf.write(u1, t)
        for k in range(10):
            u1.assign(stepper.step(u0, t, dt))
            u0.assign(u1)
            t += dt
            xf.write(u1, t)
    return


if __name__ == '__main__':
    test_heat_equation_fenics()
