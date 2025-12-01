"""
TSQVT Echo Waveform Model

This module implements the gravitational wave echo model predicted by
Twistorial Spectral Quantum Vacuum Theory (TSQVT). Echoes arise from
transient oscillations in the geometric condensation parameter ρ(x,t)
near newly formed black hole horizons.

Author: Mohamed H.M. Makraini
Date: November 2025
"""

import numpy as np
from scipy import signal
from scipy.optimize import minimize
from typing import Tuple, Dict, Optional
import warnings


class TSQVTEchoModel:
    """
    Generate gravitational wave echoes from TSQVT predictions.
    
    The echo waveform is modeled as:
        h_echo(t) = Σ_n A_n(M, M_*) h_ringdown(t - Δt_n) exp(-γΔt_n/M)
    
    where:
        - Δt_n = 2M ln(M/M_Pl) + ξn(M/M_*) is the echo delay
        - A_n is the amplitude suppression factor
        - h_ringdown is the fundamental ringdown mode
        - M is the black hole mass
        - M_* is the spectral gap scale (TSQVT parameter)
        - ξ is the spectral coupling strength
        - γ is the damping parameter
    
    Parameters
    ----------
    M_total : float
        Total black hole mass in solar masses
    M_star : float, optional
        Spectral gap scale in GeV (default: 1.0)
    xi : float, optional
        Spectral coupling parameter (default: 0.1)
    gamma : float, optional
        Damping parameter (default: 0.05)
    beta : float, optional
        Frequency modulation parameter (default: 0.02)
    q : float, optional
        Mass ratio m1/m2 (default: 1.0)
    chi_eff : float, optional
        Effective spin parameter (default: 0.0)
    """
    
    # Physical constants
    G = 6.67430e-11  # m^3 kg^-1 s^-2
    c = 299792458.0  # m/s
    M_sun = 1.98847e30  # kg
    M_Pl_GeV = 1.22e19  # Planck mass in GeV
    
    def __init__(
        self,
        M_total: float,
        M_star: float = 1.0,
        xi: float = 0.1,
        gamma: float = 0.05,
        beta: float = 0.02,
        q: float = 1.0,
        chi_eff: float = 0.0
    ):
        self.M_total = M_total  # Solar masses
        self.M_star = M_star    # GeV
        self.xi = xi
        self.gamma = gamma
        self.beta = beta
        self.q = q
        self.chi_eff = chi_eff
        
        # Compute derived quantities
        self.M_total_kg = M_total * self.M_sun
        self.M_total_seconds = self._geometric_mass(M_total)
        
        # Compute ringdown parameters
        self._compute_ringdown_params()
        
    def _geometric_mass(self, M_solar: float) -> float:
        """Convert mass to geometric units (seconds)."""
        return self.G * M_solar * self.M_sun / self.c**3
    
    def _compute_ringdown_params(self):
        """
        Compute quasi-normal mode (QNM) frequencies and damping times.
        Uses fits from Berti et al. (2009) for Kerr black holes.
        """
        a = self.chi_eff  # Dimensionless spin
        
        # Fundamental (l=2, m=2, n=0) mode
        # Frequency fit
        f1 = 1.5251 - 1.1568*(1-a)**0.1292
        f2 = 0.7000 + 1.4187*(1-a)**(-0.4990)
        
        # Damping time fit  
        q1 = 0.7000 + 1.4187*(1-a)**(-0.4990)
        q2 = -0.0739*(1-a)**(-0.3350)
        
        # QNM frequency in Hz
        self.f_qnm = (f1 + f2) / (2 * np.pi * self.M_total_seconds)
        
        # Damping time in seconds
        self.tau_qnm = 2 * self.M_total_seconds / (q1 + q2)
        
    def echo_delay(self, n: int) -> float:
        """
        Compute the time delay for the n-th echo.
        
        TSQVT Prediction:
            Δt_n = 2M ln(M/M_Pl) + ξ n (M/M_*)
        
        The first term is the light-crossing time logarithmically enhanced
        by quantum geometry. The second term is the spectral gap contribution.
        
        Parameters
        ----------
        n : int
            Echo number (n=1 for first echo)
            
        Returns
        -------
        float
            Echo delay in seconds
        """
        # Base delay: light crossing enhanced by quantum geometry
        M_ratio = (self.M_total * self.M_sun * self.c**2 / 1e9) / self.M_Pl_GeV
        t_base = 2 * self.M_total_seconds * np.log(M_ratio)
        
        # Spectral gap contribution
        # Convert M_total to GeV: M [GeV] = M [kg] * c^2 / (1.602e-10 J/GeV)
        M_GeV = self.M_total * self.M_sun * self.c**2 / 1.602e-10
        t_spectral = n * self.xi * (M_GeV / self.M_star) * self.M_total_seconds
        
        return t_base + t_spectral
    
    def echo_amplitude(self, n: int) -> float:
        """
        Compute amplitude suppression for n-th echo.
        
        TSQVT Prediction:
            A_n = A_0 exp(-γ Δt_n / M)
        
        where γ parameterizes dissipation during ρ oscillations.
        
        Parameters
        ----------
        n : int
            Echo number
            
        Returns
        -------
        float
            Amplitude factor (dimensionless)
        """
        dt_n = self.echo_delay(n)
        return np.exp(-self.gamma * dt_n / self.M_total_seconds)
    
    def frequency_shift(self, n: int) -> float:
        """
        Compute frequency shift for n-th echo due to ρ dynamics.
        
        TSQVT Prediction:
            Δf_n / f_qnm = β (ρ_horizon(t_n) - ρ_∞)
        
        where ρ_horizon oscillates during transient phase.
        
        Parameters
        ----------
        n : int
            Echo number
            
        Returns
        -------
        float
            Fractional frequency shift
        """
        # Model ρ_horizon as damped oscillation
        t_n = self.echo_delay(n)
        rho_horizon = 1 - 0.3 * np.exp(-t_n / (10 * self.M_total_seconds)) * \
                      np.cos(2 * np.pi * self.f_qnm * t_n)
        rho_infinity = 1.0
        
        return self.beta * (rho_horizon - rho_infinity)
    
    def ringdown_waveform(
        self,
        times: np.ndarray,
        t_merger: float = 0.0,
        amplitude: float = 1.0
    ) -> np.ndarray:
        """
        Generate fundamental ringdown waveform.
        
        h(t) = A exp(-(t-t_merger)/τ) cos(2πf(t-t_merger)) for t > t_merger
        
        Parameters
        ----------
        times : np.ndarray
            Time array in seconds
        t_merger : float, optional
            Merger time in seconds
        amplitude : float, optional
            Overall amplitude
            
        Returns
        -------
        np.ndarray
            Ringdown strain
        """
        h = np.zeros_like(times)
        mask = times > t_merger
        
        t_rel = times[mask] - t_merger
        envelope = amplitude * np.exp(-t_rel / self.tau_qnm)
        oscillation = np.cos(2 * np.pi * self.f_qnm * t_rel)
        
        h[mask] = envelope * oscillation
        
        return h
    
    def generate_echo(
        self,
        t_merger: float = 0.0,
        n_echoes: int = 3,
        sample_rate: int = 4096,
        duration: float = 4.0,
        include_ringdown: bool = True
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate complete echo waveform including original ringdown and echoes.
        
        Parameters
        ----------
        t_merger : float, optional
            Merger time in seconds (default: 0.0)
        n_echoes : int, optional
            Number of echoes to include (default: 3)
        sample_rate : int, optional
            Sample rate in Hz (default: 4096)
        duration : float, optional
            Total duration in seconds (default: 4.0)
        include_ringdown : bool, optional
            Include original ringdown (default: True)
            
        Returns
        -------
        times : np.ndarray
            Time array
        strain : np.ndarray
            Total strain including ringdown + echoes
        """
        # Create time array
        times = np.arange(0, duration, 1.0 / sample_rate)
        strain = np.zeros_like(times)
        
        # Add original ringdown
        if include_ringdown:
            strain += self.ringdown_waveform(times, t_merger, amplitude=1.0)
        
        # Add echoes
        for n in range(1, n_echoes + 1):
            # Compute echo properties
            dt_n = self.echo_delay(n)
            A_n = self.echo_amplitude(n)
            df_n = self.frequency_shift(n)
            
            # Time-shifted and modulated echo
            t_echo = t_merger + dt_n
            
            if t_echo > duration:
                warnings.warn(f"Echo {n} at t={t_echo:.3f}s exceeds duration")
                break
            
            # Generate echo with frequency shift
            mask = times > t_echo
            t_rel = times[mask] - t_echo
            
            f_echo = self.f_qnm * (1 + df_n)
            envelope = A_n * np.exp(-t_rel / self.tau_qnm)
            oscillation = np.cos(2 * np.pi * f_echo * t_rel)
            
            strain[mask] += envelope * oscillation
        
        return times, strain
    
    def compute_snr(
        self,
        data: np.ndarray,
        template: np.ndarray,
        psd: np.ndarray,
        sample_rate: int = 4096
    ) -> float:
        """
        Compute matched-filter SNR between data and template.
        
        SNR = <d|h> / sqrt(<h|h>)
        
        where <a|b> = 4 Re ∫ ã*(f) b̃(f) / S_n(f) df
        
        Parameters
        ----------
        data : np.ndarray
            Detector data
        template : np.ndarray
            Echo template
        psd : np.ndarray
            Power spectral density
        sample_rate : int, optional
            Sample rate in Hz
            
        Returns
        -------
        float
            Signal-to-noise ratio
        """
        # Fourier transform
        data_fft = np.fft.rfft(data)
        template_fft = np.fft.rfft(template)
        
        # Frequency array
        freqs = np.fft.rfftfreq(len(data), 1.0 / sample_rate)
        
        # Avoid division by zero
        psd_interp = np.interp(freqs, np.linspace(0, sample_rate/2, len(psd)), psd)
        psd_interp[psd_interp == 0] = np.inf
        
        # Inner products
        df = freqs[1] - freqs[0]
        integrand_dh = 4 * np.real(np.conj(data_fft) * template_fft / psd_interp)
        integrand_hh = 4 * np.real(np.conj(template_fft) * template_fft / psd_interp)
        
        dh = np.sum(integrand_dh) * df
        hh = np.sum(integrand_hh) * df
        
        if hh <= 0:
            return 0.0
        
        return dh / np.sqrt(hh)
    
    def parameter_estimation(
        self,
        data: np.ndarray,
        times: np.ndarray,
        param_bounds: Dict[str, Tuple[float, float]],
        n_iterations: int = 1000
    ) -> Dict[str, float]:
        """
        Estimate TSQVT parameters from data using maximum likelihood.
        
        Parameters
        ----------
        data : np.ndarray
            Observed strain data
        times : np.ndarray
            Time array
        param_bounds : dict
            Parameter bounds: {'M_star': (min, max), 'xi': (min, max), ...}
        n_iterations : int, optional
            Number of optimization iterations
            
        Returns
        -------
        dict
            Best-fit parameters
        """
        def neg_log_likelihood(params):
            """Negative log-likelihood for optimization."""
            M_star, xi, gamma = params
            
            # Update model
            self.M_star = M_star
            self.xi = xi
            self.gamma = gamma
            
            # Generate template
            _, template = self.generate_echo(
                t_merger=0.0,
                n_echoes=3,
                sample_rate=int(1.0 / (times[1] - times[0])),
                duration=times[-1] - times[0]
            )
            
            # Residual
            residual = data - template
            
            # Simple Gaussian likelihood
            return 0.5 * np.sum(residual**2)
        
        # Initial guess (center of bounds)
        p0 = []
        bounds = []
        for key in ['M_star', 'xi', 'gamma']:
            if key in param_bounds:
                low, high = param_bounds[key]
                p0.append(0.5 * (low + high))
                bounds.append((low, high))
        
        # Optimize
        result = minimize(
            neg_log_likelihood,
            x0=p0,
            bounds=bounds,
            method='L-BFGS-B',
            options={'maxiter': n_iterations}
        )
        
        return {
            'M_star': result.x[0],
            'xi': result.x[1],
            'gamma': result.x[2],
            'neg_log_L': result.fun
        }


def add_detector_noise(
    signal: np.ndarray,
    psd: np.ndarray,
    sample_rate: int = 4096,
    seed: Optional[int] = None
) -> np.ndarray:
    """
    Add colored Gaussian noise with specified PSD to signal.
    
    Parameters
    ----------
    signal : np.ndarray
        Clean signal
    psd : np.ndarray
        One-sided power spectral density
    sample_rate : int, optional
        Sample rate in Hz
    seed : int, optional
        Random seed for reproducibility
        
    Returns
    -------
    np.ndarray
        Noisy signal
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Generate white noise
    noise_white = np.random.randn(len(signal))
    
    # Color the noise with PSD
    noise_fft = np.fft.rfft(noise_white)
    freqs = np.fft.rfftfreq(len(signal), 1.0 / sample_rate)
    
    # Interpolate PSD to FFT frequencies
    psd_interp = np.interp(freqs, np.linspace(0, sample_rate/2, len(psd)), psd)
    
    # Apply coloring
    noise_fft *= np.sqrt(psd_interp * sample_rate / 2)
    noise_colored = np.fft.irfft(noise_fft, n=len(signal))
    
    return signal + noise_colored


if __name__ == "__main__":
    # Example usage
    print("TSQVT Echo Model - Example")
    print("=" * 50)
    
    # Initialize model for GW150914-like system
    model = TSQVTEchoModel(
        M_total=65.0,   # Solar masses
        M_star=1.0,     # GeV
        xi=0.1,
        gamma=0.05,
        q=0.85,         # Mass ratio
        chi_eff=0.0
    )
    
    print(f"Black hole mass: {model.M_total} M_sun")
    print(f"QNM frequency: {model.f_qnm:.1f} Hz")
    print(f"QNM damping time: {model.tau_qnm*1000:.2f} ms")
    print()
    
    # Compute echo delays
    print("Echo delays:")
    for n in range(1, 4):
        dt = model.echo_delay(n)
        A = model.echo_amplitude(n)
        df = model.frequency_shift(n)
        print(f"  Echo {n}: Δt = {dt*1000:.1f} ms, A = {A:.3f}, Δf/f = {df:.4f}")
    print()
    
    # Generate waveform
    times, strain = model.generate_echo(n_echoes=3, duration=2.0)
    
    print(f"Generated waveform: {len(strain)} samples")
    print(f"Duration: {times[-1]:.2f} s")
    print(f"Sample rate: {1.0/(times[1]-times[0]):.0f} Hz")
