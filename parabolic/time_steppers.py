# -*- coding: utf-8 -*-
#
'''
Routines for one time-stepping of the general equation

.. math::
    \\frac{du}{dt} = F(u).

'''


class ExplicitEuler(object):
    '''
    Explicit Euler method for :math:`u' = F(u)`.
    '''
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


class ImplicitEuler(object):
    '''
    Implicit Euler method for :math:`u' = F(u)`.
    '''
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


class Trapezoidal(object):
    '''
    Trapezoidal method for :math:`u' = F(u)`. (Known as Crank-Nicolson if
    combined with a second-order discretization in time, or used in an ODE
    context.)
    '''
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
