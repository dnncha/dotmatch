"""Numerical JACKS kernel adapted from Felicity Allen and Leopold Parts.

Copyright 2018 Felicity Allen and Leopold Parts. Distributed under the MIT
license in licenses/JACKS-MIT.txt. Source: JACKS commit
 dd5c4be5e83baa5ee79a589b0a8c3a2ac3f7d6ad, jacks/jacks/infer.py.
Changes: use NumPy instead of removed SciPy NumPy aliases; omit logging;
optionally report iteration and stopping diagnostics without changing updates.
The variational updates and stopping rule are unchanged.
"""

import numpy as np

def nandot(x1, x2):
    if len(x1.shape) == 1 and len(x2.shape) == 2:
        x1T = np.tile(x1, [x2.shape[1],1]).transpose()
        return np.nansum(np.multiply(x1T,x2), axis=0)
    elif len(x2.shape) == 1 and len(x1.shape) == 2:
        x2T = np.tile(x2, [x1.shape[0],1])
        return np.nansum(np.multiply(x1,x2T), axis=1)
    elif len(x1.shape) == 1 and len(x2.shape) == 1:
        return np.nansum(np.multiply(x1,x2))
    return None

def inferJACKSGene(data, data_err, ctrl, ctrl_err, n_iter, tol=0.1, mu0_x=1, var0_x=1.0, mu0_w=0.0, var0_w=1e4, tau_prior_strength=0.5, fixed_x=None, apply_w_hp = False, diagnostics=None):

    if not isinstance(n_iter, int) or isinstance(n_iter, bool) or n_iter < 1:
        raise ValueError("n_iter must be a positive integer")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be finite and positive")

    #Adjust estimated variances if needed
    data_err[np.isnan(data_err)] = 2.0 # very uncertain if a single replicate

    #The control can be specified once for each sample, or common across all cell lines
    if ctrl.shape != data.shape:
        #If only 1 control replicate, use mean variance from data across cell lines for that guide
        ctrl_err[np.isnan(ctrl_err)] = np.nanmean(data_err, axis=1).reshape(np.isnan(ctrl_err).shape)[np.isnan(ctrl_err)]
        y = (data.T - ctrl).T
        tau_pr_den = tau_prior_strength*1.0*((data_err**2).T + ctrl_err**2 + 1e-2).T
    else:
        #If only 1 control replicate, use data variances for ctrls as well
        ctrl_err[np.isnan(ctrl_err)] = data_err[np.isnan(ctrl_err)]
        y = data - ctrl
        tau_pr_den = tau_prior_strength*1.0*(data_err**2 + ctrl_err**2 + 1e-2)

    #Run the inference
    G,L = y.shape
    if fixed_x is None:
        x1 = mu0_x*np.ones(G)
        x2 = x1**2
    else:
        x1 = fixed_x['X1']
        x2 = fixed_x['X2']

    w1 = np.nanmedian(y, axis=0)

    tau = tau_prior_strength*1.0/tau_pr_den

    w2 = w1**2
    bound = lowerBound(x1,x2,w1,w2,y,tau)
    converged = False
    for i in range(n_iter):
        last_bound = bound
        if fixed_x is None: x1,x2 = updateX(w1,w2,tau,y,mu0_x,var0_x)
        if apply_w_hp and len(w1) > 1: mu0_w, var0_w = w1.mean(), w1.var()*3+1e-4  # hierarchical update on w (to encourage w's together - use with caution!)
        w1,w2 = updateW(x1,x2,tau,y,mu0_w,var0_w)
        tau = updateTau(x1, x2, w1, w2, y, tau_prior_strength, tau_pr_den)
        bound = lowerBound(x1,x2,w1,w2,y,tau)
        change = abs(last_bound - bound)
        if change < tol:
            converged = True
            break
    if diagnostics is not None:
        diagnostics.update(iterations=i + 1, converged=bool(converged),
                           termination_reason="bound_tolerance" if converged else "iteration_limit",
                           final_bound_change=float(change) if np.isfinite(change) else None)
    return y, tau, x1, x2, w1, w2

def updateX(w1, w2, tau, y, mu0_x, var0_x):
    x1 = (mu0_x/var0_x + nandot((y.T).T*tau,w1))/(nandot(tau,w2) + 1.0/var0_x)
    x2 = x1**2 + 1.0/(nandot(tau,w2)+1.0/var0_x)
    wadj = 0.5/len(x1)
    #Normalize by the median-emphasized mean of x
    x1m = x1.mean() + 2*wadj*np.nanmedian(x1) - wadj*x1.max() - wadj*x1.min()
    return x1/x1m, x2/x1m/x1m

def updateW(x1, x2, tau, y, mu0_ws, var0_w):
    w1 = (mu0_ws/var0_w + nandot(x1,(y.T).T*tau))/(nandot(x2,tau)+1.0/var0_w)
    w2 = w1**2 + 1.0/(nandot(x2,tau)+1.0/var0_w)
    return w1, w2

def updateTau(x1, x2, w1, w2, y, tau_prior_strength, tau_pr_den):
    b_star = y**2 - 2*y*(np.outer(x1,w1)) +np.outer(x2,w2)
    tau = (tau_prior_strength + 0.5)/(tau_pr_den + 0.5*b_star)
    return tau

def lowerBound(x1,x2,w1,w2,y,tau):
    xw = np.outer(x1,w1)
    return np.nansum(tau*(y**2 + np.outer(x2,w2) -2*xw*y))
