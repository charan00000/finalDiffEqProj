#seal populations: https://www.gbif.org/occurrence/charts?taxon_key=5219369&year=2019
#polar bear populations: https://www.gbif.org/occurrence/charts?taxon_key=2433451&year=2019

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize

def lotka_volterra(X, t, ALPHA, BETA, GAMMA, DELTA, A, B, C):
    x, y = X
    dx_dt = x*(ALPHA-BETA*y)
    dy_dt = y*(-GAMMA+DELTA*x)-(B*y*y)/(C+y)+A*y
    return [dx_dt, dy_dt]

def objective(params):
    alpha, beta, gamma, delta, a, b, c = params
    sol = odeint(lotka_volterra, initial_populations, times, args=(alpha, beta, gamma, delta, a, b, c))
    x_pred, y_pred = sol.T
    return np.sum((prey_populations - x_pred)**2 + (predator_populations - y_pred)**2)

times = [0, 1, 2, 3, 4, 5, 6]  
prey_populations = [32, 45, 13, 9, 35, 43, 40]  
predator_populations = [23, 43, 34, 46, 12, 31, 10]  

initial_populations = [prey_populations[0], predator_populations[0]]

guess_params = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]  

res = minimize(objective, guess_params, method='L-BFGS-B', bounds = [(0, None)] * 7)

a_opt, b_opt, g_opt, d_opt, A_opt, B_opt, C_opt = res.x

print(f"Optimal parameters are alpha={a_opt}, beta={b_opt}, gamma={g_opt}, delta={d_opt}, A={A_opt}, B={B_opt}, C={C_opt}")