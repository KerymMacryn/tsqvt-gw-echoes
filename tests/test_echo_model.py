"""
Unit tests for TSQVT echo model

Run with: pytest tests/

Author: Mohamed H.M. Makraini
"""

import numpy as np
import pytest
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from echo_model import TSQVTEchoModel


class TestTSQVTEchoModel:
    """Test suite for TSQVTEchoModel class."""
    
    @pytest.fixture
    def model(self):
        """Create a standard model for testing."""
        return TSQVTEchoModel(
            M_total=65.0,
            M_star=1.0,
            xi=0.1,
            gamma=0.05,
            beta=0.02
        )
    
    def test_initialization(self, model):
        """Test model initialization."""
        assert model.M_total == 65.0
        assert model.M_star == 1.0
        assert model.xi == 0.1
        assert model.gamma == 0.05
        assert model.beta == 0.02
    
    def test_geometric_mass(self, model):
        """Test geometric mass conversion."""
        M_seconds = model._geometric_mass(65.0)
        # Should be ~0.1 seconds for 65 solar masses
        assert 0.05 < M_seconds < 0.15
    
    def test_echo_delay_positive(self, model):
        """Test that echo delays are positive."""
        for n in range(1, 5):
            dt = model.echo_delay(n)
            assert dt > 0, f"Echo delay {n} should be positive"
    
    def test_echo_delay_increasing(self, model):
        """Test that echo delays increase with n."""
        delays = [model.echo_delay(n) for n in range(1, 5)]
        assert all(delays[i] < delays[i+1] for i in range(len(delays)-1)), \
            "Echo delays should increase monotonically"
    
    def test_echo_amplitude_decreasing(self, model):
        """Test that echo amplitudes decrease with n."""
        amplitudes = [model.echo_amplitude(n) for n in range(1, 5)]
        assert all(amplitudes[i] > amplitudes[i+1] for i in range(len(amplitudes)-1)), \
            "Echo amplitudes should decrease monotonically"
    
    def test_echo_amplitude_bounds(self, model):
        """Test that echo amplitudes are in (0, 1]."""
        for n in range(1, 5):
            A = model.echo_amplitude(n)
            assert 0 < A <= 1, f"Echo amplitude {n} should be in (0, 1]"
    
    def test_frequency_shift_small(self, model):
        """Test that frequency shifts are small."""
        for n in range(1, 5):
            df = model.frequency_shift(n)
            assert abs(df) < 0.5, f"Frequency shift {n} should be << 1"
    
    def test_ringdown_waveform_shape(self, model):
        """Test ringdown waveform has correct properties."""
        times = np.linspace(0, 1, 1000)
        h = model.ringdown_waveform(times, t_merger=0.5)
        
        # Should be zero before merger
        assert np.all(h[times < 0.5] == 0)
        
        # Should decay after merger
        h_after = h[times > 0.5]
        envelope = np.abs(h_after)
        # Check exponential decay (envelope should decrease)
        assert envelope[0] > envelope[-1]
    
    def test_generate_echo_length(self, model):
        """Test that generated echo has correct length."""
        sample_rate = 4096
        duration = 2.0
        
        times, strain = model.generate_echo(
            t_merger=1.0,
            n_echoes=3,
            sample_rate=sample_rate,
            duration=duration
        )
        
        expected_length = int(duration * sample_rate)
        assert len(times) == expected_length
        assert len(strain) == expected_length
    
    def test_generate_echo_contains_echoes(self, model):
        """Test that generated waveform contains echo signals."""
        times, strain = model.generate_echo(
            t_merger=0.5,
            n_echoes=3,
            sample_rate=4096,
            duration=2.0,
            include_ringdown=False
        )
        
        # Should be zero before merger
        assert np.all(strain[times < 0.5] == 0)
        
        # Should have non-zero signal after first echo time
        t_echo_1 = 0.5 + model.echo_delay(1)
        assert np.any(strain[times > t_echo_1] != 0)
    
    def test_different_masses(self):
        """Test model works for different black hole masses."""
        masses = [10.0, 30.0, 65.0, 100.0]
        
        for M in masses:
            model = TSQVTEchoModel(M_total=M)
            dt1 = model.echo_delay(1)
            A1 = model.echo_amplitude(1)
            
            assert dt1 > 0
            assert 0 < A1 <= 1
    
    def test_parameter_variations(self):
        """Test model responds correctly to parameter variations."""
        base_model = TSQVTEchoModel(M_total=65.0, xi=0.1)
        base_delay = base_model.echo_delay(1)
        
        # Larger xi should give longer delay (for spectral contribution)
        large_xi_model = TSQVTEchoModel(M_total=65.0, xi=0.5)
        large_xi_delay = large_xi_model.echo_delay(1)
        
        assert large_xi_delay > base_delay, "Larger xi should give longer delay"
    
    def test_zero_echoes(self, model):
        """Test that n_echoes=0 gives only ringdown."""
        times, strain = model.generate_echo(
            t_merger=1.0,
            n_echoes=0,
            sample_rate=4096,
            duration=2.0,
            include_ringdown=True
        )
        
        # With no echoes, signal should decay quickly
        # (only ringdown, which has short damping time)
        t_late = 1.0 + 0.5  # 500ms after merger
        mask_late = times > t_late
        
        # Signal should be very small at late times
        assert np.max(np.abs(strain[mask_late])) < 0.1


def test_add_detector_noise():
    """Test noise generation function."""
    from echo_model import add_detector_noise
    
    # Create simple signal
    signal = np.sin(2 * np.pi * 100 * np.linspace(0, 1, 4096))
    
    # Simple PSD (white noise)
    psd = np.ones(2049) * 1e-22
    
    # Add noise
    noisy = add_detector_noise(signal, psd, sample_rate=4096, seed=42)
    
    assert len(noisy) == len(signal)
    assert not np.array_equal(noisy, signal), "Noise should modify signal"
    
    # Check reproducibility with same seed
    noisy2 = add_detector_noise(signal, psd, sample_rate=4096, seed=42)
    np.testing.assert_array_equal(noisy, noisy2)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
