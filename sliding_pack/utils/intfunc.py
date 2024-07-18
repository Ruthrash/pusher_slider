# Author: Joao Moura
# Contact: jpousad@ed.ac.uk
# Date: 15/12/2020
# -------------------------------------------------------------------
# Description:
# 
# Integration functions based on scipy integrate and symbolic casadi
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Import libraries
# -------------------------------------------------------------------
import numpy as np
from scipy import integrate
import casadi as cs

# -------------------------------------------------------------------
# python integration lambda functions (1D and 2D)
# -------------------------------------------------------------------
square_np = lambda sq_side: integrate.dblquad(lambda x, y: np.sqrt((x**2)
    + (y**2)), - sq_side/2, sq_side/2, -sq_side/2, sq_side/2)[0]
quad_np = lambda sq_side: integrate.quad(lambda var: var**2,
    - sq_side/2, sq_side/2)[0]
# -------------------------------------------------------------------
# casadi auxiliary variables
# -------------------------------------------------------------------
# Fixed step Runge-Kutta 4 integrator
M = 4  # RK4 steps per interval
N = 4  # number of control intervals
sLenght = cs.SX.sym('sLenght')
xLenght = cs.SX.sym('xLenght')
yLenght = cs.SX.sym('yLenght')
x = cs.SX.sym('x')
y = cs.SX.sym('y')
DX = xLenght/(N*M)
DY = yLenght/(N*M)
# -------------------------------------------------------------------
# 1D casadi integration of g
# cost function
g = cs.Function('g_ext', [x], [DX, (x**2)*DX], {"allow_free":True})
Q = 0  # initialize cost
xx = -xLenght/2  # initialize initial cond
for n in range(N):
    for m in range(M):
        k1, k1_q = g(xx)
        k2, k2_q = g(xx + k1/2)
        k3, k3_q = g(xx + k2/2)
        k4, k4_q = g(xx + k3)
        Q += (k1_q + 2*k2_q + 2*k3_q + k4_q)/6
        xx += (k1 + 2*k2 + 2*k3 + k4)/6
quad_cs = cs.Function('quad_cs', [xLenght], [Q])
# -------------------------------------------------------------------
# 2D casadi integration of g
g = cs.Function('h_ext', [x, y], [DX, DY, (cs.sqrt((x**2)+(y**2)))*DX*DY], {"allow_free":True})
Q = 0  # initialize cost
yy = -yLenght/2  # initialize initial cond
for ny in range(N):
    for my in range(M):
        xx = -xLenght/2
        for nx in range(N):
            for mx in range(M):
                k1_x, k1_y, k1_q = g(xx, yy)
                k2_x, k2_y, k2_q = g(xx + k1_x/2, yy + k1_y/2)
                k3_x, k3_y, k3_q = g(xx + k2_x/2, yy + k2_y/2)
                k4_x, k4_y, k4_q = g(xx + k3_x, yy + k3_y)
                Q += (k1_q + 2*k2_q + 2*k3_q + k4_q)/6
                xx += (k1_x + 2*k2_x + 2*k3_x + k4_x)/6
        yy += (k1_y + 2*k2_y + 2*k3_y + k4_y)/6
rect_cs = cs.Function('rect_cs', [xLenght, yLenght], [Q])
square_cs = cs.Function('square_cs', [sLenght], [rect_cs(sLenght, sLenght)])
