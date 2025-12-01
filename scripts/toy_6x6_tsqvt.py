"""
6x6 J-self-adjoint toy model for TSQVT
Reproducible script: eigenvalue flow, gap, M_mn, fidelity, perturbation sensitivity, adiabatic bound.

Dependencies:
  numpy, scipy, matplotlib

Run:
  python toy_6x6_tsqvt.py
"""

import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

# ---------------------------
# Reproducibility and helpers
# ---------------------------
#np.random.seed(42)
np.random.seed(123)          # reproducibilidad (opcional)
n_trials = 1000              # en test_perturbations llamarlo con n_trials=n_trials


def sort_eigvals_and_vecs(vals, vecs):
    """Sort eigenvalues (real) and corresponding eigenvectors by ascending eigenvalue."""
    idx = np.argsort(vals.real)
    return vals[idx], vecs[:, idx]

def j_inner(u, v, J):
    """Krein pairing <u|v>_K = u^* J v (complex conjugate transpose)"""
    return np.vdot(u, J.dot(v))

def os_projector(J):
    """Osterwalder-Schrader projector Pos = (I + J)/2"""
    n = J.shape[0]
    return 0.5 * (np.eye(n) + J)

# ---------------------------
# Model definition
# ---------------------------
def build_6x6_D(ρ, a_vals, g_of_ρ, W):
    """
    Build the 6x6 J-self-adjoint Dirac-twistor matrix:
      D(ρ) = [[ A(ρ),  g(ρ) W ],
              [ -g(ρ) W^†, -A(ρ) ]]
    where A(ρ) = diag(a1(ρ), a2(ρ), a3(ρ)).
    """
    A = np.diag(a_vals(ρ))
    g = g_of_ρ(ρ)
    top = np.hstack([A, g * W])
    bottom = np.hstack([-g * W.conj().T, -A])
    D = np.vstack([top, bottom])
    return D

# Example parameter functions
def a_vals_example(ρ):
    # smooth lifting of degeneracy: linear shift
    base = np.array([ -0.6, -0.2, 0.1 ])
    return base + 0.5 * ρ * np.array([0.1, 0.05, -0.02])

def g_of_p_example(ρ):
    # coupling grows smoothly from 0 to 1
    return 0.5 * (1 - np.cos(np.pi * ρ))  # smooth monotone from 0 to 1

# Random coupling W with controlled norm (recalibrated)
W0 = 1e-3 * (np.random.randn(3,3) + 1j * np.random.randn(3,3))


# Fundamental symmetry J = diag(I3, -I3)
J = np.diag([1,1,1,-1,-1,-1])

# ---------------------------
# Spectral flow and diagnostics
# ---------------------------
def spectral_flow(ρ_grid, a_vals, g_of_ρ, W):
    n = 6
    eigvals = np.zeros((len(ρ_grid), n), dtype=np.complex128)
    eigvecs = np.zeros((len(ρ_grid), n, n), dtype=np.complex128)
    for i, ρ in enumerate(ρ_grid):
        D = build_6x6_D(ρ, a_vals, g_of_ρ, W)
        vals, vecs = la.eig(D)
        vals, vecs = sort_eigvals_and_vecs(vals, vecs)
        eigvals[i, :] = vals
        eigvecs[i, :, :] = vecs
    return eigvals, eigvecs

# Compute M_mn via finite difference derivative of D(ρ)
def compute_M_matrix(ρ, a_vals, g_of_ρ, W, eps=1e-6):
    D_plus = build_6x6_D(ρ + eps, a_vals, g_of_ρ, W)
    D_minus = build_6x6_D(ρ - eps, a_vals, g_of_ρ, W)
    dD_dp = (D_plus - D_minus) / (2 * eps)
    vals, vecs = la.eig(build_6x6_D(ρ, a_vals, g_of_ρ, W))
    vals, vecs = sort_eigvals_and_vecs(vals, vecs)
    M = np.zeros((6,6), dtype=np.complex128)
    for m in range(6):
        for n in range(6):
            M[m,n] = np.vdot(vecs[:,m].conj(), dD_dp.dot(vecs[:,n]))
    return M, vals, vecs

# ---------------------------
# Fidelity and perturbation tests
# ---------------------------
def fidelity_between_states(psi, phi):
    # fidelity for pure states: |<psi|phi>|^2 normalized
    psi = psi / la.norm(psi)
    phi = phi / la.norm(phi)
    return np.abs(np.vdot(psi.conj(), phi))**2

def test_perturbations(ρ_target, a_vals, g_of_ρ, W, perturbation_scale=1e-3, n_trials=200):
    D0 = build_6x6_D(ρ_target, a_vals, g_of_ρ, W)
    vals0, vecs0 = la.eig(D0)
    vals0, vecs0 = sort_eigvals_and_vecs(vals0, vecs0)
    # choose a target eigenstate index (e.g., middle of cluster)
    target_idx = 2  # example
    psi0 = vecs0[:, target_idx]
    fidelities = []
    gaps = []
    for _ in range(n_trials):
        # small random Hermitian perturbation respecting J-self-adjoint structure approximately
        H = perturbation_scale * (np.random.randn(6,6) + 1j*np.random.randn(6,6))
        H = 0.5 * (H + H.conj().T)
        Dp = D0 + H
        vals_ρ, vecs_ρ = la.eig(Dp)
        vals_ρ, vecs_ρ = sort_eigvals_and_vecs(vals_ρ, vecs_ρ)
        psi_ρ = vecs_ρ[:, target_idx]
        fidelities.append(fidelity_between_states(psi0, psi_ρ))
        # compute minimal gap in perturbed spectrum
        sorted_reals = np.sort(vals_ρ.real)
        gaps.append(np.min(np.diff(sorted_reals)))
    return np.array(fidelities), np.array(gaps)
    
# Check J-self-adjointness (numerical tolerance)
def check_J_self_adjoint(D, J, tol=1e-10):
    lhs = D.conj().T.dot(J)
    rhs = J.dot(D)
    return np.linalg.norm(lhs - rhs) < tol

# Example usage
D_test = build_6x6_D(0.5, a_vals_example, g_of_p_example, W0)
assert check_J_self_adjoint(D_test, J, tol=1e-8), "Warning: D(ρ) not J-self-adjoint within tolerance"

def compute_dDdp_converged(ρ, a_vals, g_of_ρ, W, eps_list=[1e-6,1e-7,1e-8]):
    dD_list = []
    for eps in eps_list:
        Dp = build_6x6_D(ρ+eps, a_vals, g_of_ρ, W)
        Dm = build_6x6_D(ρ-eps, a_vals, g_of_ρ, W)
        dD = (Dp - Dm) / (2*eps)
        dD_list.append(dD)
    # compute pairwise norms to check convergence
    diffs = [np.linalg.norm(dD_list[i] - dD_list[i+1]) for i in range(len(dD_list)-1)]
    return dD_list[-1], diffs

# Example
dDdp, diffs = compute_dDdp_converged(0.5, a_vals_example, g_of_p_example, W0)
print("Derivative diffs:", diffs)

def compute_diagnostics(ρ_grid, a_vals, g_of_ρ, W):
    M_max = 0.0
    Delta_min = np.inf
    for ρ in ρ_grid:
        D = build_6x6_D(ρ, a_vals, g_of_ρ, W)
        vals, vecs = la.eig(D)
        vals_sorted = np.sort(vals.real)
        Delta_min = min(Delta_min, np.min(np.abs(np.diff(vals_sorted))))
        dDdp, _ = compute_dDdp_converged(ρ, a_vals, g_of_ρ, W)
        # compute M matrix
        _, vecs_sorted = sort_eigvals_and_vecs(vals, vecs)
        for m in range(6):
            for n in range(6):
                if m != n:
                    M_val = np.abs(np.vdot(vecs_sorted[:,m].conj(), dDdp.dot(vecs_sorted[:,n])))
                    M_max = max(M_max, M_val)
    return M_max, Delta_min

M_max, Delta_min = compute_diagnostics(np.linspace(0,1,41), a_vals_example, g_of_p_example, W0)
print("M_max, Delta_min:", M_max, Delta_min)


def estimate_dotp_max(M_max, Delta, I_size=6, epsilon_target=1e-6, hbar=1.0):
    factor = np.sqrt(epsilon_target / (2*(I_size-1)))
    dotp = (hbar * Delta**2 / M_max) * factor
    return dotp

dotp_max_est = estimate_dotp_max(M_max, Delta_min, I_size=6, epsilon_target=1e-6)
print("Estimated dotp_max for eps=1e-6:", dotp_max_est)


from scipy.integrate import solve_ivp

def p_protocol(t, T):
    # smooth cosine protocol from 0 to 1 on [0,T]
    return 0.5*(1 - np.cos(np.pi * t / T))

def H_t(t, T, a_vals, g_of_ρ, W):
    ρ = p_protocol(t, T)
    return build_6x6_D(ρ, a_vals, g_of_ρ, W)

def tdse_rhs(t, psi_flat, T, a_vals, g_of_ρ, W):
    psi = psi_flat.view(np.complex128)
    H = H_t(t, T, a_vals, g_of_ρ, W)
    dpsi = -1j * H.dot(psi)
    return dpsi.view(np.float64)  # real vector for solver

# initial state: pick instantaneous eigenvector at t0
T = 1e4  # choose T so that dotp ~ desired value
t_span = (0.0, T)
D0 = build_6x6_D(0.0, a_vals_example, g_of_p_example, W0)
vals0, vecs0 = la.eig(D0)
vals0, vecs0 = sort_eigvals_and_vecs(vals0, vecs0)
psi0 = vecs0[:,2]  # example initial eigenstate
psi0_flat = psi0.view(np.float64)

sol = solve_ivp(tdse_rhs, t_span, psi0_flat, args=(T, a_vals_example, g_of_p_example, W0),
                method='RK45', atol=1e-8, rtol=1e-7)
psi_final = sol.y[:, -1].view(np.complex128)
# project onto instantaneous eigenstate at final p=1
Df = build_6x6_D(1.0, a_vals_example, g_of_p_example, W0)
valf, vecf = la.eig(Df)
valf, vecf = sort_eigvals_and_vecs(valf, vecf)
target_idx = 2
F = fidelity_between_states(psi_final, vecf[:, target_idx])
print("TDSE final fidelity:", F)


# ---------------------------
# Adiabatic bound diagnostic
# ---------------------------
def adiabatic_diagnostic(ρ_grid, a_vals, g_of_ρ, W, dotp_max):
    # compute M_max and Delta along p_grid
    M_max = 0.0
    Delta_min = np.inf
    for ρ in ρ_grid:
        M, vals, vecs = compute_M_matrix(ρ, a_vals, g_of_ρ, W)
        # gap among cluster (all 6 levels)
        sorted_vals = np.sort(vals.real)
        deltas = np.abs(np.diff(sorted_vals))
        Delta_min = min(Delta_min, np.min(deltas))
        # off-diagonal maxima
        for m in range(6):
            for n in range(6):
                if m != n:
                    M_max = max(M_max, np.abs(M[m,n]))
    # adiabatic bound estimate for P_max (leading term)
    P_est = (M_max * dotp_max / (1.0 * (Delta_min**2)))**2
    return {'M_max': M_max, 'Delta_min': Delta_min, 'P_est': P_est}

# ---------------------------
# Main execution: produce figures
# ---------------------------
def main():
    # parameter grid
    ρ_grid = np.linspace(0.0, 1.0, 201)
    eigvals, eigvecs = spectral_flow(ρ_grid, a_vals_example, g_of_p_example, W0)

    # Figure 1: eigenvalue flow
    plt.figure(figsize=(6.5,4))
    for i in range(6):
        plt.plot(ρ_grid, eigvals[:, i].real, label=f'$\lambda_{i+1}$')
    plt.xlabel('ρ')
    plt.ylabel('Eigenvalues (real part)')
    plt.title('Eigenvalue flow $\lambda_i(ρ)$ (6x6 model)')
    plt.grid(True)
    plt.tight_layout()
    plt.legend(loc='upper left', fontsize='small')
    plt.savefig('../figures/fig_eigenflow.png', dpi=300)

    # Figure 2: minimal gap vs ρ
    gaps = np.min(np.diff(eigvals.real, axis=1), axis=1)
    plt.figure(figsize=(6.5,3))
    plt.plot(ρ_grid, gaps, '-k')
    plt.xlabel('ρ')
    plt.ylabel('Minimal adjacent gap')
    plt.title('Minimal adjacent gap vs ρ')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../figures/fig_gap.png', dpi=300)

    # Compute M_max and Delta_min diagnostics
    diag = adiabatic_diagnostic(np.linspace(0,1,41), a_vals_example, g_of_p_example, W0, dotp_max=1e-3)
    print("Adiabatic diagnostic:", diag)

    # Figure 3: perturbation fidelity histogram at ρ_target
    ρ_target = 1.0  # final measurement regime
    fidelities, pert_gaps = test_perturbations(ρ_target, a_vals_example, g_of_p_example, W0,
                                              perturbation_scale=1e-3, n_trials=500)
    plt.figure(figsize=(6.5,3))
    plt.hist(fidelities, bins=40, color='C0', alpha=0.8)
    plt.xlabel('Fidelity')
    plt.ylabel('Counts')
    plt.title(f'Fidelity distribution under perturbations (ρ={ρ_target})')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../figures/fig_fidelity_hist.png', dpi=300)

    # Figure 4: fidelity vs perturbation scale sweep
    scales = np.logspace(-5, -2, 8)
    mean_fids = []
    for s in scales:
        fids, _ = test_perturbations(ρ_target, a_vals_example, g_of_p_example, W0,
                                     perturbation_scale=s, n_trials=200)
        mean_fids.append(np.mean(fids))
    plt.figure(figsize=(6.5,3))
    plt.semilogx(scales, 1 - np.array(mean_fids), marker='o')
    plt.xlabel('Perturbation scale (norm)')
    plt.ylabel('Mean infidelity (1 - F)')
    plt.title('Mean infidelity vs perturbation scale')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../figures/fig_infidelity_vs_scale.png', dpi=300)

    plt.show()

if __name__ == '__main__':
    main()
