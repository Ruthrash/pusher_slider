## Author: Joao Moura
## Contact: jpousad@ed.ac.uk
## Date: 15/12/2020
## -------------------------------------------------------------------
## Description:
## 
## Functions for outputting different nominal trajectories
## -------------------------------------------------------------------

## -------------------------------------------------------------------
## Import libraries
## -------------------------------------------------------------------
import numpy as np
import casadi as cs

## Generate Nominal Trajectory (line)
def generate_traj_line(x_f, y_f, N, N_MPC):
    x_nom = np.linspace(0.0, x_f, N)
    y_nom = np.linspace(0.0, y_f, N)
    # return x_nom, y_nom
    return np.concatenate((x_nom, x_f+x_nom[1:N_MPC+1]), axis=0), np.concatenate((y_nom, y_f+y_nom[1:N_MPC+1]), axis=0)
def generate_traj_circle(theta_i, theta_f, radious, N, N_MPC):
    s = np.linspace(theta_i, theta_f, N)
    x_nom = radious*np.cos(s)
    y_nom = radious*np.sin(s)
    # initial position at the origin
    x_nom -= x_nom[0]
    y_nom -= y_nom[0]
    # return x_nom, y_nom
    return np.concatenate((x_nom, x_nom[1:N_MPC+1]), axis=0), np.concatenate((y_nom, y_nom[1:N_MPC+1]), axis=0)
def generate_traj_ellipse(theta_i, theta_f, radious_x, radious_y, N, N_MPC):
    s = np.linspace(theta_i, theta_f, N)
    x_nom = radious_x*np.cos(s)
    y_nom = radious_y*np.sin(s)
    # initial position at the origin
    x_nom -= x_nom[0]
    y_nom -= y_nom[0]
    # return x_nom, y_nom
    return np.concatenate((x_nom, x_nom[1:N_MPC+1]), axis=0), np.concatenate((y_nom, y_nom[1:N_MPC+1]), axis=0)
def generate_traj_eight(side_lenght, N, N_MPC):
    s = np.linspace(0.0, 2*np.pi, N)
    x_nom = side_lenght*np.sin(s)
    y_nom = side_lenght*np.sin(s)*np.cos(s)
    # return x_nom, y_nom
    return np.concatenate((x_nom, x_nom[1:N_MPC+1]), axis=0), np.concatenate((y_nom, y_nom[1:N_MPC+1]), axis=0)
def compute_nomState_from_nomTraj(x_data, y_data, dt):
    # assign two first state trajectories
    x0_nom = x_data
    x1_nom = y_data
    # compute diff for plannar traj
    Dx0_nom = np.diff(x0_nom)
    Dx1_nom = np.diff(x1_nom)
    # compute traj angle 
    ND = len(Dx0_nom)
    x2_nom = np.empty(ND)
    theta = 0.0
    for i in range(ND):
        c, s = np.cos(theta), np.sin(theta)
        R = np.array(((c, s), (-s, c)))
        Dx_new = R.dot(np.array((Dx0_nom[i],Dx1_nom[i])))
        theta += np.arctan2(Dx_new[1], Dx_new[0])
        x2_nom[i] = theta
    x2_nom = np.append(x2_nom, x2_nom[-1])
    Dx2_nom = np.diff(x2_nom)
    # specify angle of the pusher relative to slider
    x3_nom = np.zeros(x0_nom.shape)
    Dx3_nom = np.diff(x3_nom)
    # stack state and derivative of state
    # x_nom = np.vstack((x0_nom, x1_nom, x2_nom, x3_nom))
    x_nom = cs.horzcat(x0_nom, x1_nom, x2_nom, x3_nom).T
    dx_nom = cs.horzcat(Dx0_nom, Dx1_nom, Dx2_nom, Dx3_nom).T/dt
    return x_nom, dx_nom
