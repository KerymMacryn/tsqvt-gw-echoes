"""
Matched Filter Search for Gravitational Wave Echoes

Implements matched filtering techniques to search for TSQVT-predicted echoes
in LIGO/Virgo data. Computes signal-to-noise ratio (SNR) timeseries and
identifies candidate echo detections.

Author: Mohamed H.M. Makraini
Date: November 2025
"""

import numpy as np
from scipy import signal, fft
from scipy.interpolate import interp1d
from typing import Tuple, Dict, List, Optional
import warnings


class MatchedFilterSearcher:
    """
    Perform matched filter search for gravitational wave echoes.
    
    The matched filter SNR is computed as:
        ρ(t) = <d(t)|h> / sqrt(<h|h>)
    
    where the inner product is:
        <a|b> = 4 Re ∫ ã*(f) b̃(f) / S_n(f) df
    
    Parameters
    ----------
    sample_rate : int, optional
        Sample rate in Hz (default: 4096)
    fmin : float, optional
        Minimum frequency for analysis (default: 30.0 Hz)
    fmax : float, optional
        Maximum frequency for analysis (default: 500.0 Hz)
    """
    
    def __init__(
        self,
        sample_rate: int = 4096,
        fmin: float = 30.0,
        fmax: float = 500.0
    ):
        self.sample_rate = sample_rate
        self.fmin = fmin
        self.fmax = fmax
        self.dt = 1.0 / sample_rate
    
    def compute_psd(
        self,
        data: np.ndarray,
        fft_length: int = 4096,
        overlap: float = 0.5
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimate power spectral density using Welch's method.
        
        Parameters
        ----------
        data : np.ndarray
            Time-domain data
        fft_length : int, optional
            FFT length (default: 4096)
        overlap : float, optional
            Overlap fraction (default: 0.5)
            
        Returns
        -------
        freqs : np.ndarray
            Frequency array
        psd : np.ndarray
            One-sided power spectral density
        """
        nperseg = fft_length
        noverlap = int(overlap * nperseg)
        
        freqs, psd = signal.welch(
            data,
            fs=self.sample_rate,
            nperseg=nperseg,
            noverlap=noverlap,
            scaling='density',
            average='median'
        )
        
        return freqs, psd
    
    def whiten_data(
        self,
        data: np.ndarray,
        psd: np.ndarray,
        freqs: np.ndarray
    ) -> np.ndarray:
        """
        Whiten data using given PSD.
        
        Parameters
        ----------
        data : np.ndarray
            Time-domain data
        psd : np.ndarray
            Power spectral density
        freqs : np.ndarray
            Frequency array for PSD
            
        Returns
        -------
        np.ndarray
            Whitened data
        """
        # FFT of data
        data_fft = fft.rfft(data)
        data_freqs = fft.rfftfreq(len(data), self.dt)
        
        # Interpolate PSD to data frequencies
        psd_interp_func = interp1d(
            freqs,
            psd,
            kind='linear',
            bounds_error=False,
            fill_value=(psd[0], psd[-1])
        )
        psd_interp = psd_interp_func(data_freqs)
        
        # Avoid division by zero
        psd_interp[psd_interp == 0] = np.inf
        
        # Whiten in frequency domain
        data_fft_white = data_fft / np.sqrt(psd_interp * self.sample_rate / 2)
        
        # Apply frequency mask
        mask = (data_freqs >= self.fmin) & (data_freqs <= self.fmax)
        data_fft_white[~mask] = 0
        
        # Transform back
        data_white = fft.irfft(data_fft_white, n=len(data))
        
        return data_white
    
    def matched_filter_snr(
        self,
        data: np.ndarray,
        template: np.ndarray,
        psd: np.ndarray,
        freqs: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute matched filter SNR timeseries.
        
        Parameters
        ----------
        data : np.ndarray
            Detector strain data (time domain)
        template : np.ndarray
            Echo template (time domain, same length as data)
        psd : np.ndarray
            Power spectral density
        freqs : np.ndarray
            Frequency array for PSD
            
        Returns
        -------
        times : np.ndarray
            Time array
        snr : np.ndarray
            SNR timeseries
        """
        # Ensure same length
        if len(data) != len(template):
            raise ValueError("Data and template must have same length")
        
        # FFT
        data_fft = fft.rfft(data)
        template_fft = fft.rfft(template)
        
        # Frequency array
        data_freqs = fft.rfftfreq(len(data), self.dt)
        
        # Interpolate PSD
        psd_interp_func = interp1d(
            freqs,
            psd,
            kind='linear',
            bounds_error=False,
            fill_value=(psd[0], psd[-1])
        )
        psd_interp = psd_interp_func(data_freqs)
        psd_interp[psd_interp == 0] = np.inf
        
        # Frequency mask
        mask = (data_freqs >= self.fmin) & (data_freqs <= self.fmax)
        
        # Compute inner product <d|h> in frequency domain
        integrand = np.conj(data_fft) * template_fft / psd_interp
        integrand[~mask] = 0  # Zero outside band
        
        # IFFT to get SNR timeseries
        snr_complex = 4.0 * fft.irfft(integrand, n=len(data))
        
        # Normalize by <h|h>
        hh_integrand = np.abs(template_fft)**2 / psd_interp
        hh_integrand[~mask] = 0
        
        df = data_freqs[1] - data_freqs[0]
        hh = 4.0 * np.sum(hh_integrand) * df
        
        if hh <= 0:
            warnings.warn("Template norm is zero or negative")
            snr = np.zeros(len(data))
        else:
            snr = np.abs(snr_complex) / np.sqrt(hh)
        
        times = np.arange(len(data)) * self.dt
        
        return times, snr
    
    def find_peaks(
        self,
        snr: np.ndarray,
        threshold: float = 4.0,
        min_distance: float = 0.05
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Find peaks in SNR timeseries above threshold.
        
        Parameters
        ----------
        snr : np.ndarray
            SNR timeseries
        threshold : float, optional
            SNR threshold (default: 4.0)
        min_distance : float, optional
            Minimum distance between peaks in seconds (default: 0.05)
            
        Returns
        -------
        peak_indices : np.ndarray
            Indices of peaks
        peak_snrs : np.ndarray
            SNR values at peaks
        """
        # Convert min_distance to samples
        min_samples = int(min_distance * self.sample_rate)
        
        # Find peaks
        peak_indices, properties = signal.find_peaks(
            snr,
            height=threshold,
            distance=min_samples
        )
        
        peak_snrs = properties['peak_heights']
        
        # Sort by SNR (descending)
        sort_idx = np.argsort(peak_snrs)[::-1]
        peak_indices = peak_indices[sort_idx]
        peak_snrs = peak_snrs[sort_idx]
        
        return peak_indices, peak_snrs
    
    def search_for_echoes(
        self,
        data: np.ndarray,
        template: np.ndarray,
        search_window: Tuple[float, float] = (0.1, 2.0),
        snr_threshold: float = 4.0,
        t_merger: float = 0.0
    ) -> Dict:
        """
        Complete echo search pipeline.
        
        Parameters
        ----------
        data : np.ndarray
            Detector data (preprocessed)
        template : np.ndarray
            Echo template
        search_window : tuple of float, optional
            Time window to search (relative to merger) in seconds
        snr_threshold : float, optional
            SNR detection threshold
        t_merger : float, optional
            Merger time in seconds (relative to data start)
            
        Returns
        -------
        dict
            Search results containing:
            - 'snr': SNR timeseries
            - 'times': Time array
            - 'peaks': Peak times
            - 'peak_snrs': Peak SNR values
            - 'max_snr': Maximum SNR
            - 'n_candidates': Number of candidates above threshold
        """
        # Estimate PSD from data
        print("Computing PSD...")
        freqs, psd = self.compute_psd(data)
        
        # Compute matched filter SNR
        print("Computing matched filter SNR...")
        times, snr = self.matched_filter_snr(data, template, psd, freqs)
        
        # Define search window mask
        window_start = t_merger + search_window[0]
        window_end = t_merger + search_window[1]
        window_mask = (times >= window_start) & (times <= window_end)
        
        # Find peaks in search window
        print("Finding peaks...")
        snr_windowed = snr.copy()
        snr_windowed[~window_mask] = 0  # Zero outside window
        
        peak_indices, peak_snrs = self.find_peaks(snr_windowed, threshold=snr_threshold)
        peak_times = times[peak_indices]
        
        results = {
            'snr': snr,
            'times': times,
            'peaks': peak_times,
            'peak_snrs': peak_snrs,
            'max_snr': np.max(snr[window_mask]) if np.any(window_mask) else 0.0,
            'n_candidates': len(peak_times),
            'search_window': search_window,
            'threshold': snr_threshold
        }
        
        print(f"Found {len(peak_times)} candidates above SNR={snr_threshold}")
        
        return results
    
    def compute_significance(
        self,
        snr: np.ndarray,
        background_duration: float = 1.0,
        n_trials: int = 1000
    ) -> Tuple[float, float]:
        """
        Estimate statistical significance using time-shifted background.
        
        Parameters
        ----------
        snr : np.ndarray
            SNR timeseries
        background_duration : float, optional
            Duration for background estimation in seconds
        n_trials : int, optional
            Number of time-shift trials
            
        Returns
        -------
        p_value : float
            False alarm probability
        sigma : float
            Significance in standard deviations
        """
        max_snr_signal = np.max(snr)
        
        # Estimate background by looking at off-peak regions
        # Simple approach: use first/last portions as background
        n_bg_samples = int(background_duration * self.sample_rate)
        background = np.concatenate([snr[:n_bg_samples], snr[-n_bg_samples:]])
        
        # Estimate false alarm rate
        n_above_threshold = np.sum(background >= max_snr_signal)
        p_value = n_above_threshold / len(background)
        
        # Convert to sigma
        from scipy.stats import norm
        if p_value > 0:
            sigma = norm.ppf(1 - p_value)
        else:
            sigma = norm.ppf(1 - 1/len(background))  # Upper limit
        
        return p_value, sigma


def matched_filter_search(
    data_file: str,
    template: np.ndarray,
    search_window: Tuple[float, float] = (0.1, 2.0),
    snr_threshold: float = 4.0,
    sample_rate: int = 4096
) -> Dict:
    """
    Convenience function for matched filter search.
    
    Parameters
    ----------
    data_file : str
        Path to preprocessed data file (HDF5 or numpy)
    template : np.ndarray
        Echo template
    search_window : tuple, optional
        Search window relative to merger (seconds)
    snr_threshold : float, optional
        Detection threshold
    sample_rate : int, optional
        Sample rate
        
    Returns
    -------
    dict
        Search results
    """
    # Load data
    if data_file.endswith('.npy'):
        data = np.load(data_file)
    elif data_file.endswith('.h5') or data_file.endswith('.hdf5'):
        import h5py
        with h5py.File(data_file, 'r') as f:
            # Assume strain is stored under 'strain' key
            data = f['strain'][:]
    else:
        raise ValueError("Unsupported file format. Use .npy or .h5")
    
    # Initialize searcher
    searcher = MatchedFilterSearcher(sample_rate=sample_rate)
    
    # Perform search
    results = searcher.search_for_echoes(
        data,
        template,
        search_window=search_window,
        snr_threshold=snr_threshold
    )
    
    return results


if __name__ == "__main__":
    # Example usage
    print("Matched Filter Search - Example")
    print("=" * 60)
    
    # Simulate some data
    sample_rate = 4096
    duration = 4.0
    t = np.arange(0, duration, 1.0 / sample_rate)
    
    # Simulated signal: damped sinusoid
    f0 = 250.0  # Hz
    tau = 0.1   # seconds
    signal_true = np.exp(-t / tau) * np.sin(2 * np.pi * f0 * t)
    signal_true[:int(0.5 * sample_rate)] = 0  # Zero before t=0.5s
    
    # Add noise
    noise = np.random.randn(len(t)) * 0.5
    data = signal_true + noise
    
    # Template (same as signal)
    template = signal_true.copy()
    
    # Initialize searcher
    searcher = MatchedFilterSearcher(sample_rate=sample_rate)
    
    # Compute PSD
    print("\nComputing PSD...")
    freqs, psd = searcher.compute_psd(data)
    print(f"PSD computed: {len(freqs)} frequency bins")
    
    # Matched filter
    print("\nComputing matched filter SNR...")
    times, snr = searcher.matched_filter_snr(data, template, psd, freqs)
    print(f"Max SNR: {np.max(snr):.2f}")
    
    # Find peaks
    print("\nFinding peaks...")
    peak_indices, peak_snrs = searcher.find_peaks(snr, threshold=5.0)
    print(f"Found {len(peak_indices)} peaks above threshold")
    
    for i, (idx, snr_val) in enumerate(zip(peak_indices[:3], peak_snrs[:3])):
        print(f"  Peak {i+1}: t={times[idx]:.3f}s, SNR={snr_val:.2f}")
