#!/usr/bin/env python3
"""
Lambda Suppression Figure for TSQVT Paper
==========================================

Generates publication-quality figure showing the cosmological constant
suppression mechanism using actual TSQVT numerical results.

Based on: test_matched_filter_cc.py validation results
Author: Generated for TSQVT paper
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib import gridspec

# Set publication-quality parameters
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11
plt.rcParams['font.family'] = 'serif'
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['mathtext.fontset'] = 'cm'

# Physical constants (from test results)
M_Pl_GeV = 1.221e19  # GeV
m_star = 1.0  # GeV
Lambda_obs_planck = 1.260e-123  # dimensionless, Planck units
Lambda_QFT_planck = 1.0  # dimensionless, Planck units

# Compute suppression factor
ratio_m_M = m_star / M_Pl_GeV  # 8.191e-20
ln_M_m = np.log(M_Pl_GeV / m_star)  # 43.95
rho2 = (ratio_m_M**6) / ln_M_m  # 6.871e-117
Lambda_eff = rho2 * Lambda_QFT_planck  # 6.871e-117

# Fine-tuning measures
Delta_before = Lambda_QFT_planck / Lambda_obs_planck  # 7.94e122
Delta_after = Lambda_eff / Lambda_obs_planck  # 5.45e6
suppression_orders = np.log10(Delta_before / Delta_after)  # 116.2

print("═"*70)
print("LAMBDA SUPPRESSION FIGURE - NUMERICAL VALUES")
print("═"*70)
print(f"M_Pl = {M_Pl_GeV:.3e} GeV")
print(f"m_* = {m_star} GeV")
print(f"Λ_QFT = {Lambda_QFT_planck} (Planck units)")
print(f"Λ_obs = {Lambda_obs_planck:.3e} (Planck units)")
print(f"⟨ρ²⟩ = {rho2:.3e}")
print(f"Λ_eff = {Lambda_eff:.3e} (Planck units)")
print(f"Suppression: {suppression_orders:.1f} orders of magnitude")
print(f"Δ_before = {Delta_before:.2e}")
print(f"Δ_after = {Delta_after:.2e}")
print("═"*70)

# Create figure with GridSpec for flexible layout
fig = plt.figure(figsize=(14, 9))
gs = gridspec.GridSpec(3, 3, figure=fig, height_ratios=[1, 0.8, 0.4],
                       width_ratios=[1, 1, 1], hspace=0.4, wspace=0.3)

# ═══════════════════════════════════════════════════════════
# PANEL 1: Naive QFT (top left)
# ═══════════════════════════════════════════════════════════
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 10)
ax1.axis('off')

# Title
ax1.text(5, 9.5, 'Naive QFT', ha='center', va='top', 
         fontsize=14, fontweight='bold', color='darkblue')

# Box around panel
fancy_box1 = FancyBboxPatch((0.5, 1), 9, 8, boxstyle="round,pad=0.1",
                            edgecolor='darkblue', facecolor='lightblue', 
                            alpha=0.2, linewidth=2)
ax1.add_patch(fancy_box1)

# Planck scale indicator
ax1.plot([1, 9], [8, 8], 'b-', linewidth=2)
ax1.text(5, 8.3, r'$M_{\rm Pl} \sim 10^{19}\,{\rm GeV}$', 
         ha='center', fontsize=10)

# Vacuum modes (many oscillations)
x_osc = np.linspace(1, 9, 200)
for y_center in [7, 6.5, 6, 5.5, 5, 4.5]:
    y_osc = y_center + 0.3 * np.sin(10 * x_osc)
    ax1.plot(x_osc, y_osc, 'r-', linewidth=1.5, alpha=0.7)

# Sum symbol
ax1.text(5, 4, r'$\sum_{\mathrm{all\ modes}}$', 
         ha='center', va='center', fontsize=20)

# Result box
result_box1 = FancyBboxPatch((2, 1.2), 6, 1.5, boxstyle="round,pad=0.1",
                             edgecolor='red', facecolor='lightyellow', 
                             linewidth=2)
ax1.add_patch(result_box1)
ax1.text(5, 2, r'$\Lambda_{\rm QFT} \sim \mathcal{O}(1)$', 
         ha='center', va='center', fontsize=13, fontweight='bold')
ax1.text(5, 1.5, '(Planck units)', ha='center', va='center', fontsize=9)

# ═══════════════════════════════════════════════════════════
# PANEL 2: TSQVT Filter (top center)
# ═══════════════════════════════════════════════════════════
ax2 = fig.add_subplot(gs[0, 1])
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 10)
ax2.axis('off')

# Title
ax2.text(5, 9.5, 'TSQVT Spectral Filter', ha='center', va='top',
         fontsize=14, fontweight='bold', color='purple')

# Box around panel
fancy_box2 = FancyBboxPatch((0.5, 1), 9, 8, boxstyle="round,pad=0.1",
                            edgecolor='purple', facecolor='lavender',
                            alpha=0.2, linewidth=2)
ax2.add_patch(fancy_box2)

# Filter funnel
funnel_x = [2, 8, 6.5, 3.5, 2]
funnel_y = [8, 8, 5, 5, 8]
ax2.fill(funnel_x, funnel_y, color='purple', alpha=0.3)
ax2.plot(funnel_x, funnel_y, 'purple', linewidth=2)

# Filtered modes (fewer, smaller oscillations)
x_filt = np.linspace(3.5, 6.5, 100)
for y_center, amp in zip([6.5, 6, 5.5], [0.2, 0.15, 0.1]):
    y_filt = y_center + amp * np.sin(8 * x_filt)
    ax2.plot(x_filt, y_filt, 'g-', linewidth=2, alpha=0.8)

# Suppression formula box
formula_box = FancyBboxPatch((1.5, 7), 7, 1.2, boxstyle="round,pad=0.1",
                             edgecolor='darkgreen', facecolor='lightgreen',
                             linewidth=1.5)
ax2.add_patch(formula_box)
ax2.text(5, 7.6, r'$\langle\rho^2\rangle = \frac{(m_\ast/M_{\rm Pl})^6}{\ln(M_{\rm Pl}/m_\ast)}$',
         ha='center', va='center', fontsize=10)

# m_* scale
ax2.text(5, 3.8, r'$m_\ast \sim 1\,{\rm GeV}$', ha='center',
         fontsize=11, style='italic', color='darkgreen')

# Result box
result_box2 = FancyBboxPatch((2, 1.2), 6, 1.5, boxstyle="round,pad=0.1",
                             edgecolor='purple', facecolor='lavender',
                             linewidth=2)
ax2.add_patch(result_box2)
ax2.text(5, 2.2, r'$\Lambda_{\rm eff} \sim 10^{-117}$', 
         ha='center', va='center', fontsize=13, fontweight='bold')
ax2.text(5, 1.5, f'suppressed by {suppression_orders:.0f} orders',
         ha='center', va='center', fontsize=9, fontweight='bold')

# ═══════════════════════════════════════════════════════════
# PANEL 3: Observable Lambda (top right)
# ═══════════════════════════════════════════════════════════
ax3 = fig.add_subplot(gs[0, 2])
ax3.set_xlim(0, 10)
ax3.set_ylim(0, 10)
ax3.axis('off')

# Title
ax3.text(5, 9.5, r'Observable $\Lambda$', ha='center', va='top',
         fontsize=14, fontweight='bold', color='darkgreen')

# Box around panel
fancy_box3 = FancyBboxPatch((0.5, 1), 9, 8, boxstyle="round,pad=0.1",
                            edgecolor='darkgreen', facecolor='lightgreen',
                            alpha=0.2, linewidth=2)
ax3.add_patch(fancy_box3)

# CMB/telescope icon
circle = plt.Circle((5, 6.5), 0.8, color='darkgreen', fill=False, linewidth=2)
ax3.add_patch(circle)
for angle in np.linspace(0, 2*np.pi, 8, endpoint=False):
    x_ray = 5 + 0.8 * np.cos(angle)
    y_ray = 6.5 + 0.8 * np.sin(angle)
    x_end = 5 + 1.2 * np.cos(angle)
    y_end = 6.5 + 1.2 * np.sin(angle)
    ax3.plot([x_ray, x_end], [y_ray, y_end], 'darkgreen', linewidth=2)

# Small residual oscillation
x_tiny = np.linspace(3, 7, 50)
y_tiny = 4.5 + 0.15 * np.sin(6 * x_tiny)
ax3.plot(x_tiny, y_tiny, 'g-', linewidth=3)

# Observed value box
result_box3 = FancyBboxPatch((2, 2.5), 6, 1.5, boxstyle="round,pad=0.1",
                             edgecolor='darkgreen', facecolor='palegreen',
                             linewidth=2)
ax3.add_patch(result_box3)
ax3.text(5, 3.5, r'$\Lambda_{\rm obs} \sim 10^{-123}$',
         ha='center', va='center', fontsize=13, fontweight='bold')
ax3.text(5, 2.9, '(dark energy)', ha='center', va='center', fontsize=9)

# Residual tuning indicator
residual_box = FancyBboxPatch((2.5, 1), 5, 1, boxstyle="round,pad=0.05",
                              edgecolor='orange', facecolor='lightyellow',
                              linewidth=1.5)
ax3.add_patch(residual_box)
ax3.text(5, 1.65, r'Residual: $\Delta_{\rm after} \sim 10^6$',
         ha='center', va='center', fontsize=9, color='darkorange')
ax3.text(5, 1.25, '(tractable hierarchy)', ha='center', va='center',
         fontsize=8, style='italic')

# ═══════════════════════════════════════════════════════════
# ARROWS BETWEEN PANELS
# ═══════════════════════════════════════════════════════════
# Create arrows in a separate subplot that spans all columns
ax_arrows = fig.add_subplot(gs[0, :])
ax_arrows.set_xlim(0, 30)
ax_arrows.set_ylim(0, 10)
ax_arrows.axis('off')

# Arrow 1: QFT → Filter
arrow1 = FancyArrowPatch((9.8, 5), (11, 5), arrowstyle='->', 
                         mutation_scale=30, linewidth=3, 
                         color='blue', alpha=0.7)
ax_arrows.add_patch(arrow1)
ax_arrows.text(10.4, 5.5, 'Spectral', ha='center', fontsize=9, color='blue')
ax_arrows.text(10.4, 4.5, 'Exclusion', ha='center', fontsize=9, color='blue')

# Arrow 2: Filter → Observable
arrow2 = FancyArrowPatch((19.2, 5), (20.4, 5), arrowstyle='->',
                         mutation_scale=30, linewidth=3,
                         color='purple', alpha=0.7)
ax_arrows.add_patch(arrow2)
ax_arrows.text(19.8, 5.5, 'Geometric', ha='center', fontsize=9, color='purple')
ax_arrows.text(19.8, 4.5, 'projection', ha='center', fontsize=9, color='purple')

# ═══════════════════════════════════════════════════════════
# MIDDLE ROW: Matched Filter Analogy Box
# ═══════════════════════════════════════════════════════════
ax_analogy = fig.add_subplot(gs[1, :])
ax_analogy.set_xlim(0, 10)
ax_analogy.set_ylim(0, 10)
ax_analogy.axis('off')

# Inset box
analogy_box = FancyBboxPatch((1, 2), 8, 6, boxstyle="round,pad=0.15",
                             edgecolor='black', facecolor='lightyellow',
                             linewidth=2)
ax_analogy.add_patch(analogy_box)

# Title
ax_analogy.text(5, 7.5, 'Matched Filter Analogy', ha='center', va='center',
                fontsize=13, fontweight='bold')

# LIGO example
ax_analogy.text(1.5, 6, 'LIGO:', ha='left', va='center',
                fontsize=11, fontweight='bold', color='darkblue')
ax_analogy.text(5, 5.5, r'Raw strain data $\longrightarrow$ GW signal',
                ha='center', va='center', fontsize=10)
ax_analogy.text(5, 4.9, 'Suppresses incoherent noise, extracts chirp template',
                ha='center', va='center', fontsize=9, style='italic')

# Separator line
ax_analogy.plot([1.5, 8.5], [4.3, 4.3], 'k--', linewidth=1, alpha=0.5)

# TSQVT parallel
ax_analogy.text(1.5, 3.7, 'TSQVT:', ha='left', va='center',
                fontsize=11, fontweight='bold', color='purple')
ax_analogy.text(5, 3.2, r'Vacuum $\Lambda_{\rm QFT} \longrightarrow$ Observable $\Lambda$',
                ha='center', va='center', fontsize=10)
ax_analogy.text(5, 2.6, 'Suppresses spectral noise, extracts geometric modes',
                ha='center', va='center', fontsize=9, style='italic')

# ═══════════════════════════════════════════════════════════
# BOTTOM ROW: Logarithmic Scale Bar
# ═══════════════════════════════════════════════════════════
ax_scale = fig.add_subplot(gs[2, :])
ax_scale.set_xlim(-4, 4)
ax_scale.set_ylim(-1, 1)
ax_scale.axis('off')

# Scale bar background
scale_bar = mpatches.Rectangle((-3.5, 0), 7, 0.3, 
                               facecolor='lightgray', edgecolor='black',
                               linewidth=1)
ax_scale.add_patch(scale_bar)

# Tick marks and labels
tick_positions = [-3.5, -1.75, 0, 1.75, 3.3]
tick_labels = [r'$10^{0}$', r'$10^{-60}$', r'$10^{-120}$', 
               r'$\Lambda_{\rm eff}$', r'$\Lambda_{\rm obs}$']
for pos, label in zip(tick_positions, tick_labels):
    ax_scale.plot([pos, pos], [0, 0.3], 'k-', linewidth=2)
    ax_scale.text(pos, -0.2, label, ha='center', fontsize=10)

# Suppression span (brace)
ax_scale.plot([-3.5, 1.75], [-0.4, -0.4], 'b-', linewidth=2)
ax_scale.plot([-3.5, -3.5], [-0.35, -0.45], 'b-', linewidth=2)
ax_scale.plot([1.75, 1.75], [-0.35, -0.45], 'b-', linewidth=2)
ax_scale.text(-0.875, -0.65, '116 orders suppression', ha='center',
              fontsize=10, color='blue', fontweight='bold')

# Residual span (brace)
ax_scale.plot([1.75, 3.3], [-0.7, -0.7], 'orange', linewidth=2)
ax_scale.plot([1.75, 1.75], [-0.65, -0.75], 'orange', linewidth=2)
ax_scale.plot([3.3, 3.3], [-0.65, -0.75], 'orange', linewidth=2)
ax_scale.text(2.525, -0.9, r'$\sim 10^6$ residual', ha='center',
              fontsize=9, color='darkorange', fontweight='bold')

# Overall title
fig.suptitle('TSQVT Suppression of the Cosmological Constant',
             fontsize=16, fontweight='bold', y=0.98)

# Save figure
plt.tight_layout(rect=[0, 0.02, 1, 0.96])
output_file = '../figures/lambda_suppression_figure_MATPLOTLIB.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\n✓ Figure saved to: {output_file}")

# Also save as PDF for publication
output_pdf = '../figures/lambda_suppression_figure_MATPLOTLIB.pdf'
plt.savefig(output_pdf, bbox_inches='tight')
print(f"✓ PDF version saved to: {output_pdf}")

plt.show()

print("\n" + "═"*70)
print("FIGURE GENERATION COMPLETE")
print("="*70)
print("Use PNG for preview, PDF for paper submission")
print("="*70)
