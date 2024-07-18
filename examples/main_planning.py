# Author: Joao Moura
# Date: 21/08/2020
#  -------------------------------------------------------------------
# Description:
#  This script implements a non-linear program (NLP) model predictive controller (MPC)
#  for tracking a trajectory of a square slider object with a single
#  and sliding contact pusher.
#  -------------------------------------------------------------------

#  import libraries
#  -------------------------------------------------------------------
import sys
import yaml
import numpy as np
import pandas as pd
import casadi as cs
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#  -------------------------------------------------------------------
import sliding_pack
#  -------------------------------------------------------------------

# Get config files
#  -------------------------------------------------------------------
planning_config = sliding_pack.load_config('planning_config.yaml')
#  -------------------------------------------------------------------

# Set Problem constants
#  -------------------------------------------------------------------
T = 2.5  # time of the simulation is seconds
freq = 25  # number of increments per second
show_anim = True
save_to_file = True
#  -------------------------------------------------------------------
# Computing Problem constants
#  -------------------------------------------------------------------
dt = 1.0/freq  # sampling time
N = int(T*freq)  # total number of iterations
#  -------------------------------------------------------------------

# define system dynamics
#  -------------------------------------------------------------------
dyn = sliding_pack.dyn.Sys_sq_slider_quasi_static_ellip_lim_surf(
        planning_config['dynamics'],
        planning_config['TO']['contactMode']
)
#  -------------------------------------------------------------------

# Generate Nominal Trajectory
#  -------------------------------------------------------------------
X_goal = planning_config['TO']['X_goal']
print(X_goal)
x0_nom, x1_nom = sliding_pack.traj.generate_traj_line(X_goal[0], X_goal[1], N, 0)
# x0_nom, x1_nom = sliding_pack.traj.generate_traj_line(0.3, 0.4, N, 0)
# x0_nom, x1_nom = sliding_pack.traj.generate_traj_circle(-np.pi/2, 3*np.pi/2, 0.1, N, 0)
# x0_nom, x1_nom = sliding_pack.traj.generate_traj_eight(0.2, N, 0)
#  -------------------------------------------------------------------
# stack state and derivative of state
X_nom_val, _ = sliding_pack.traj.compute_nomState_from_nomTraj(x0_nom, x1_nom, dt)
#  ------------------------------------------------------------------

# Compute nominal actions for sticking contact
#  ------------------------------------------------------------------
optObj = sliding_pack.to.buildOptObj(
        dyn, N, planning_config['TO'], dt=dt, useGoalFlag=True)
# Set obstacles
#  ------------------------------------------------------------------
if optObj.numObs==0:
    obsCentre = None
    obsRadius = None
elif optObj.numObs==2:
    obsCentre = [[0.2, 0.2], [0.1, 0.5]]
    obsRadius = [0.05, 0.05]
elif optObj.numObs==3:
    obsCentre = [[0.2, 0.2], [0.0, 0.4], [0.3, 0.0]]
    obsRadius = [0.05, 0.05, 0.05]
#  ------------------------------------------------------------------
x_init = [0., 0., -20.*(np.pi/180.), 0.]
beta = [
    planning_config['dynamics']['xLenght'],
    planning_config['dynamics']['yLenght'],
    planning_config['dynamics']['pusherRadious']
]
# x_init = [0., 0., 340.*(np.pi/180.), 0.]
resultFlag, X_nom_val_opt, U_nom_val_opt, other_opt, _, t_opt = optObj.solveProblem(
        0, x_init, beta,
        X_warmStart=X_nom_val,
        obsCentre=obsCentre, obsRadius=obsRadius,
        X_goal_val=X_goal)
f_d = cs.Function('f_d', [dyn.x, dyn.u], [dyn.x + dyn.f(dyn.x, dyn.u, beta)*dt])
f_rollout = f_d.mapaccum(N-1)
print('comp time: ', t_opt)
p_new = cs.Function('p_new', [dyn.x], [dyn.p(dyn.x, beta)])
p_map = p_new.map(N)
X_pusher_opt = p_map(X_nom_val_opt)
#  ------------------------------------------------------------------


if save_to_file:
    #  Save data to file using pandas
    #  -------------------------------------------------------------------
    df_state = pd.DataFrame(
                    np.array(cs.vertcat(X_nom_val_opt,X_pusher_opt)).transpose(),
                    columns=['x_slider', 'y_slider', 'theta_slider', 'psi_pusher', 'x_pusher', 'y_pusher'])
    df_state.to_csv('planning_positive_angle_state.csv',
                    float_format='%.5f')
    df_action = pd.DataFrame(
                    np.array(U_nom_val_opt).transpose(),
                    columns=['u0', 'u1', 'u3', 'u3'])
    df_action.to_csv('planning_positive_angle_action.csv',
                     float_format='%.5f')
    #  -------------------------------------------------------------------

# Animation
#  -------------------------------------------------------------------
plt.rcParams['figure.dpi'] = 150
if show_anim:
    #  ---------------------------------------------------------------
    fig, ax = sliding_pack.plots.plot_nominal_traj(
                x0_nom, x1_nom, plot_title='')
    # add computed nominal trajectory
    X_nom_val_opt = np.array(X_nom_val_opt)
    ax.plot(X_nom_val_opt[0, :], X_nom_val_opt[1, :], color='blue',
            linewidth=2.0, linestyle='dashed')
    # add obstacles
    if optObj.numObs > 0:
        for i in range(len(obsCentre)):
            circle_i = plt.Circle(obsCentre[i], obsRadius[i], color='b')
            ax.add_patch(circle_i)
    # set window size
    fig.set_size_inches(8, 6, forward=True)
    # get slider and pusher patches
    dyn.set_patches(ax, X_nom_val_opt, beta)
    # call the animation
    ani = animation.FuncAnimation(
            fig,
            dyn.animate,
            fargs=(ax, X_nom_val_opt, beta),
            frames=N-1,
            interval=dt*1000,  # microseconds
            blit=True,
            repeat=False,
    )
    # to save animation, uncomment the line below:
    # ani.save('planning_with_obstacles.mp4', fps=25, extra_args=['-vcodec', 'libx264'])
#  -------------------------------------------------------------------

# # Plot Optimization Results
# #  -------------------------------------------------------------------
# fig, axs = plt.subplots(3, 4, sharex=True)
# fig.set_size_inches(10, 10, forward=True)
# t_Nx = np.linspace(0, T, N)
# t_Nu = np.linspace(0, T, N-1)
# ctrl_g = dyn.g_u.map(N-1)
# ctrl_g_val = ctrl_g(U_nom_val_opt, other_opt)
# #  -------------------------------------------------------------------
# # plot position
# for i in range(dyn.Nx):
#     axs[0, i].plot(t_Nx, X_nom_val[i, 0:N].T, color='red',
#                    linestyle='--', label='nom')
#     axs[0, i].plot(t_Nx, X_nom_val_opt[i, 0:N].T, color='blue',
#                    linestyle='--', label='plan')
#     handles, labels = axs[0, i].get_legend_handles_labels()
#     axs[0, i].legend(handles, labels)
#     axs[0, i].set_xlabel('time [s]')
#     axs[0, i].set_ylabel('x%d' % i)
#     axs[0, i].grid()
# #  -------------------------------------------------------------------
# # plot extra variables
# for i in range(dyn.Nz):
#     axs[1, 2].plot(t_Nu, other_opt[i, :].T, label='s%d' % i)
# handles, labels = axs[1, 2].get_legend_handles_labels()
# axs[1, 2].legend(handles, labels)
# axs[1, 2].set_xlabel('time [s]')
# axs[1, 2].set_ylabel('extra vars')
# axs[1, 2].grid()
# #  -------------------------------------------------------------------
# # plot constraints
# for i in range(dyn.Ng_u):
#     axs[1, 3].plot(t_Nu, ctrl_g_val[i, :].T, label='g%d' % i)
# handles, labels = axs[1, 3].get_legend_handles_labels()
# axs[1, 3].legend(handles, labels)
# axs[1, 3].set_xlabel('time [s]')
# axs[1, 3].set_ylabel('constraints')
# axs[1, 3].grid()
# #  -------------------------------------------------------------------
# # plot actions
# for i in range(dyn.Nu):
#     axs[2, i].plot(t_Nu, U_nom_val_opt[i, 0:N-1].T, color='blue',
#                    linestyle='--', label='plan')
#     handles, labels = axs[2, i].get_legend_handles_labels()
#     axs[2, i].legend(handles, labels)
#     axs[2, i].set_xlabel('time [s]')
#     axs[2, i].set_ylabel('u%d' % i)
#     axs[2, i].grid()
# #  -------------------------------------------------------------------

#  -------------------------------------------------------------------
plt.show()
#  -------------------------------------------------------------------
