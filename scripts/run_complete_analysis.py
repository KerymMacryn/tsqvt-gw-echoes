"""
Complete TSQVT Echo Analysis Pipeline - SYNTHETIC DATA VERSION

This script runs the full analysis using SYNTHETIC data instead of downloading
from LIGO. This allows testing the complete pipeline without gwpy/h5py.

Usage:
    python run_complete_analysis_synthetic.py --event GW150914 --M_star 1.0 --xi 0.1

Author: Mohamed H.M. Makraini
Date: November 2025
"""

import argparse
import json
import os
import sys
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

# Add src to path
SRC_DIR = os.path.join(os.path.dirname(__file__), "..", "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

# Only import echo_model and matched_filter (no data_retrieval needed)
from echo_model import TSQVTEchoModel
from matched_filter import MatchedFilterSearcher


# Event catalog (same as in data_retrieval.py)
EVENT_CATALOG = {
    'GW150914': {
        'gps_time': 1126259462.4,
        'mass_1': 36.2,
        'mass_2': 29.1,
        'detectors': ['H1', 'L1']
    },
    'GW151226': {
        'gps_time': 1135136350.6,
        'mass_1': 14.2,
        'mass_2': 7.5,
        'detectors': ['H1', 'L1']
    },
    'GW170814': {
        'gps_time': 1186741861.5,
        'mass_1': 30.5,
        'mass_2': 25.3,
        'detectors': ['H1', 'L1', 'V1']
    },
}


def generate_synthetic_data(duration, sample_rate, M_total, include_echo=True, 
                           M_star=1.0, xi=0.1, gamma=0.05, snr=5.0):
    """
    Generate synthetic gravitational wave data with echoes.
    
    Parameters
    ----------
    duration : float
        Duration in seconds
    sample_rate : int
        Sample rate in Hz
    M_total : float
        Total mass in solar masses
    include_echo : bool
        Whether to include echoes
    M_star : float
        Spectral gap scale (GeV)
    xi : float
        Spectral coupling
    gamma : float
        Damping parameter
    snr : float
        Signal-to-noise ratio
        
    Returns
    -------
    times : np.ndarray
        Time array
    data : np.ndarray
        Synthetic strain data with noise
    clean_signal : np.ndarray
        Clean signal without noise (for reference)
    """
    n_samples = int(duration * sample_rate)
    times = np.arange(n_samples) / sample_rate
    
    # Generate signal
    model = TSQVTEchoModel(
        M_total=M_total,
        M_star=M_star,
        xi=xi,
        gamma=gamma
    )
    
    _, signal = model.generate_echo(
        t_merger=duration / 2,
        n_echoes=3 if include_echo else 0,
        sample_rate=sample_rate,
        duration=duration,
        include_ringdown=True
    )
    
    # Generate colored noise (LIGO-like)
    noise = np.random.randn(n_samples)
    noise_fft = np.fft.rfft(noise)
    freqs = np.fft.rfftfreq(n_samples, 1.0/sample_rate)
    
    # 1/f^2 PSD shape (simplified LIGO noise curve)
    psd = 1e-23 * (1 + (30.0/np.maximum(freqs, 1.0))**2)
    noise_fft *= np.sqrt(psd * sample_rate / 2)
    noise_colored = np.fft.irfft(noise_fft, n=n_samples)
    
    # Scale signal to achieve desired SNR
    signal_power = np.sqrt(np.mean(signal**2))
    noise_power = np.sqrt(np.mean(noise_colored**2))
    
    if signal_power > 0:
        signal_scaled = signal * (snr * noise_power / signal_power)
    else:
        signal_scaled = signal
    
    # Combined data
    data = signal_scaled + noise_colored
    
    return times, data, signal_scaled


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="TSQVT Gravitational Wave Echo Analysis (Synthetic Data)"
    )

    # Event selection
    parser.add_argument(
        "--event",
        type=str,
        default="GW150914",
        help="Event name (e.g., GW150914)",
    )

    # TSQVT parameters
    parser.add_argument(
        "--M_star",
        type=float,
        default=1.0,
        help="Spectral gap scale in GeV (default: 1.0)",
    )

    parser.add_argument(
        "--xi",
        type=float,
        default=0.1,
        help="Spectral coupling parameter (default: 0.1)",
    )

    parser.add_argument(
        "--gamma",
        type=float,
        default=0.05,
        help="Damping parameter (default: 0.05)",
    )

    parser.add_argument(
        "--beta",
        type=float,
        default=0.02,
        help="Frequency modulation (default: 0.02)",
    )

    # Analysis parameters
    parser.add_argument(
        "--n_echoes",
        type=int,
        default=3,
        help="Number of echoes to include (default: 3)",
    )

    parser.add_argument(
        "--snr_threshold",
        type=float,
        default=4.0,
        help="SNR detection threshold (default: 4.0)",
    )
    
    parser.add_argument(
        "--signal_snr",
        type=float,
        default=5.0,
        help="Injected signal SNR (default: 5.0)",
    )

    parser.add_argument(
        "--search_window",
        nargs=2,
        type=float,
        default=[0.1, 2.0],
        metavar=("T_MIN", "T_MAX"),
        help="Search window in seconds after merger (default: 0.1 2.0)",
    )

    # I/O
    parser.add_argument(
        "--output_dir",
        type=str,
        default="results_synthetic",
        help="Output base directory for results (default: 'results_synthetic')",
    )

    return parser.parse_args()


def main():
    """Main analysis pipeline."""
    args = parse_arguments()

    print("=" * 80)
    print("TSQVT GRAVITATIONAL WAVE ECHO ANALYSIS - SYNTHETIC DATA")
    print("=" * 80)
    print(f"Event: {args.event} (simulated)")
    print(f"Signal SNR: {args.signal_snr}")
    print(
        f"TSQVT Parameters: M_*={args.M_star} GeV, "
        f"ξ={args.xi}, γ={args.gamma}, β={args.beta}"
    )
    print("=" * 80)
    print("\nNOTE: Using synthetic data (no LIGO download required)")
    print("      For real LIGO data analysis, install gwpy and use run_complete_analysis.py")
    print("=" * 80)

    # Create output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(args.output_dir, f"{args.event}_synthetic_{timestamp}")
    os.makedirs(output_dir, exist_ok=True)

    print(f"\nOutput directory: {output_dir}\n")

    # ======================================================================
    # STEP 1: Generate synthetic data
    # ======================================================================
    print("\n" + "=" * 80)
    print("STEP 1: GENERATE SYNTHETIC DATA")
    print("=" * 80)

    # Get event info
    if args.event not in EVENT_CATALOG:
        print(f"ERROR: Unknown event {args.event}")
        print(f"Available events: {list(EVENT_CATALOG.keys())}")
        return

    event_info = EVENT_CATALOG[args.event]
    M_total = event_info["mass_1"] + event_info["mass_2"]
    q = event_info["mass_1"] / event_info["mass_2"]

    print(f"\nBlack hole parameters:")
    print(f"  Total mass: {M_total:.1f} M_sun")
    print(f"  Mass ratio: {q:.2f}")

    sample_rate = 4096
    duration = 32.0

    print(f"\nGenerating synthetic data...")
    print(f"  Duration: {duration} s")
    print(f"  Sample rate: {sample_rate} Hz")
    print(f"  Signal SNR: {args.signal_snr}")

    times_array, data_array, clean_signal = generate_synthetic_data(
        duration=duration,
        sample_rate=sample_rate,
        M_total=M_total,
        include_echo=True,
        M_star=args.M_star,
        xi=args.xi,
        gamma=args.gamma,
        snr=args.signal_snr
    )

    print(f"  Generated: {len(data_array)} samples")

    # Save synthetic data
    data_dir = os.path.join(output_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    np.save(os.path.join(data_dir, "times.npy"), times_array)
    np.save(os.path.join(data_dir, "strain_noisy.npy"), data_array)
    np.save(os.path.join(data_dir, "strain_clean.npy"), clean_signal)

    # ======================================================================
    # STEP 2: Generate TSQVT echo template
    # ======================================================================
    print("\n" + "=" * 80)
    print("STEP 2: GENERATE TSQVT ECHO TEMPLATE")
    print("=" * 80)

    # Initialize echo model
    echo_model = TSQVTEchoModel(
        M_total=M_total,
        M_star=args.M_star,
        xi=args.xi,
        gamma=args.gamma,
        beta=args.beta,
        q=q,
    )

    print("\nRingdown properties:")
    print(f"  QNM frequency: {echo_model.f_qnm:.1f} Hz")
    print(f"  QNM damping time: {echo_model.tau_qnm * 1000:.2f} ms")

    print("\nEcho predictions:")
    for n in range(1, args.n_echoes + 1):
        dt = echo_model.echo_delay(n)
        A = echo_model.echo_amplitude(n)
        df = echo_model.frequency_shift(n)
        print(f"  Echo {n}: Δt={dt * 1000:.1f} ms, A={A:.3f}, Δf/f={df:.4f}")

    # Generate template (echoes only, for matched filtering)
    print(f"\nGenerating template...")

    template_times, template_strain = echo_model.generate_echo(
        t_merger=duration / 2.0,
        n_echoes=args.n_echoes,
        sample_rate=sample_rate,
        duration=duration,
        include_ringdown=False,  # Only echoes for template
    )

    print(f"Template generated: {len(template_strain)} samples")

    # Save template
    template_dir = os.path.join(output_dir, "templates")
    os.makedirs(template_dir, exist_ok=True)
    np.save(os.path.join(template_dir, "template.npy"), template_strain)

    # ======================================================================
    # STEP 3: Matched filter search
    # ======================================================================
    print("\n" + "=" * 80)
    print("STEP 3: MATCHED FILTER SEARCH")
    print("=" * 80)

    searcher = MatchedFilterSearcher(
        sample_rate=sample_rate,
        fmin=30.0,
        fmax=500.0,
    )

    print("\nSearch parameters:")
    print(
        f"  Window: {args.search_window[0]:.2f} - "
        f"{args.search_window[1]:.2f} s after merger"
    )
    print(f"  SNR threshold: {args.snr_threshold}")

    results = searcher.search_for_echoes(
        data_array,
        template_strain,
        search_window=tuple(args.search_window),
        snr_threshold=args.snr_threshold,
        t_merger=duration / 2.0,
    )

    print("\nSearch results:")
    print(f"  Max SNR: {results['max_snr']:.2f}")
    print(f"  Candidates: {results['n_candidates']}")

    if results["n_candidates"] > 0:
        print("\n  Top candidates:")
        for i, (t, snr) in enumerate(
            zip(results["peaks"][:5], results["peak_snrs"][:5]), start=1
        ):
            dt_merger = t - duration / 2.0
            print(
                f"    {i}. t={t:.3f}s "
                f"(Δt={dt_merger * 1000:.1f}ms), SNR={snr:.2f}"
            )

    # Compute significance
    print("\nComputing statistical significance...")
    p_value, sigma = searcher.compute_significance(results["snr"])
    print(f"  p-value: {p_value:.4e}")
    print(f"  Significance: {sigma:.2f} σ")

    # Save results
    results_file = os.path.join(output_dir, "search_results.json")

    results_dict = {
        "event": args.event + "_synthetic",
        "data_type": "synthetic",
        "signal_snr": args.signal_snr,
        "timestamp": timestamp,
        "tsqvt_params": {
            "M_star": args.M_star,
            "xi": args.xi,
            "gamma": args.gamma,
            "beta": args.beta,
        },
        "max_snr": float(results["max_snr"]),
        "n_candidates": int(results["n_candidates"]),
        "p_value": float(p_value),
        "significance_sigma": float(sigma),
        "peaks": [float(t) for t in results["peaks"]],
        "peak_snrs": [float(s) for s in results["peak_snrs"]],
    }

    with open(results_file, "w", encoding="utf-8") as f:
        json.dump(results_dict, f, indent=2)

    print(f"\nResults saved to: {results_file}")

    # ======================================================================
    # STEP 4: Generate plots
    # ======================================================================
    print("\n" + "=" * 80)
    print("STEP 4: GENERATE PLOTS")
    print("=" * 80)

    figures_dir = os.path.join(output_dir, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    # Plot 1: Data (clean vs noisy)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))
    
    ax1.plot(times_array, clean_signal, 'b-', linewidth=1.0, alpha=0.8)
    ax1.axvline(duration / 2.0, color='r', linestyle='--', alpha=0.5, label='Merger')
    ax1.set_ylabel('Strain', fontsize=12)
    ax1.set_title(f'Clean Signal (SNR={args.signal_snr})', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(times_array, data_array, 'k-', linewidth=0.5, alpha=0.7)
    ax2.axvline(duration / 2.0, color='r', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Time (s)', fontsize=12)
    ax2.set_ylabel('Strain', fontsize=12)
    ax2.set_title('Noisy Data (with detector noise)', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "data.png"), dpi=300)
    print("  Saved: data.png")
    plt.close()

    # Plot 2: SNR timeseries
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.plot(results["times"], results["snr"], "k-", linewidth=0.5, alpha=0.7)
    ax.axhline(
        args.snr_threshold,
        color='r',
        linestyle="--",
        linewidth=2,
        label=f"Threshold ({args.snr_threshold})",
    )
    ax.axvline(
        duration / 2.0,
        color='g',
        linestyle="--",
        alpha=0.5,
        label="Merger",
    )

    # Mark expected echo times
    for n in range(1, args.n_echoes + 1):
        t_echo = duration / 2.0 + echo_model.echo_delay(n)
        ax.axvline(t_echo, color='orange', linestyle=':', alpha=0.7, linewidth=1.0)

    if results["n_candidates"] > 0:
        ax.plot(
            results["peaks"],
            results["peak_snrs"],
            "ro",
            markersize=8,
            label="Detected",
            zorder=5
        )

    ax.set_xlabel("Time (s)", fontsize=12)
    ax.set_ylabel("SNR", fontsize=12)
    ax.set_title(
        f"{args.event} (Synthetic): Matched Filter SNR",
        fontsize=14,
    )
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "snr_timeseries.png"), dpi=300)
    print("  Saved: snr_timeseries.png")
    plt.close()

    # Plot 3: Template
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(template_times, template_strain, "b-", linewidth=1.0)
    ax.axvline(duration / 2.0, color='r', linestyle='--', alpha=0.5, label='Merger')
    
    # Mark echo positions
    for n in range(1, args.n_echoes + 1):
        t_echo = duration / 2.0 + echo_model.echo_delay(n)
        ax.axvline(t_echo, color='orange', linestyle=':', alpha=0.7)
    
    ax.set_xlabel("Time (s)", fontsize=12)
    ax.set_ylabel("Strain", fontsize=12)
    ax.set_title(
        f"TSQVT Echo Template (M_*={args.M_star} GeV)",
        fontsize=14,
    )
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "template.png"), dpi=300)
    print("  Saved: template.png")
    plt.close()

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nAll results saved to: {output_dir}")
    print("\nSummary:")
    print(f"  Injected SNR: {args.signal_snr}")
    print(f"  Max detected SNR: {results['max_snr']:.2f}")
    print(f"  Candidates: {results['n_candidates']}")
    print(f"  Significance: {sigma:.2f} σ (p={p_value:.2e})")

    if sigma > 5.0:
        print(f"\n  *** CLEAR DETECTION (>5σ) ***")
    elif sigma > 3.0:
        print(f"\n  *** SIGNIFICANT DETECTION (>3σ) ***")
    elif sigma > 2.0:
        print(f"\n  *** MARGINAL DETECTION ({sigma:.1f}σ) ***")
    else:
        print("\n  No significant detection")
    
    print("\n" + "=" * 80)
    print("NOTE: This analysis used SYNTHETIC data")
    print("For real LIGO data:")
    print("  1. Install gwpy: pip install gwpy h5py")
    print("  2. Run: python run_complete_analysis.py --event GW150914")
    print("=" * 80)


if __name__ == "__main__":
    main()