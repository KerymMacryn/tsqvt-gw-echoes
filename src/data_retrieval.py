"""
LIGO Data Retrieval and Preprocessing

Download gravitational wave strain data from LIGO/Virgo Open Science Center
and perform standard preprocessing (whitening, bandpass filtering, etc.)

Author: Mohamed H.M. Makraini
Date: November 2025
"""

import numpy as np
import os
import json
from typing import List, Tuple, Dict, Optional, TYPE_CHECKING
import warnings

# Use TYPE_CHECKING to avoid runtime errors with optional dependencies
if TYPE_CHECKING:
    from gwpy.timeseries import TimeSeries
    from gwpy.signal import filter_design
else:
    TimeSeries = None
    filter_design = None

# Try importing gwpy
try:
    from gwpy.timeseries import TimeSeries
    from gwpy.signal import filter_design
    GWPY_AVAILABLE = True
except ImportError:
    GWPY_AVAILABLE = False
    warnings.warn("gwpy not available. Install with: pip install gwpy")

# Try importing h5py
try:
    import h5py
    H5PY_AVAILABLE = True
except ImportError:
    H5PY_AVAILABLE = False
    warnings.warn("h5py not available. Install with: pip install h5py")


class LIGODataRetriever:
    """
    Download and preprocess LIGO/Virgo gravitational wave data.
    
    Parameters
    ----------
    cache_dir : str, optional
        Directory for caching downloaded data (default: 'data/raw/')
    """
    
    # Known GWTC-3 events (subset for testing)
    GWTC3_EVENTS = {
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
        'GW170817': {  # Binary neutron star
            'gps_time': 1187008882.4,
            'mass_1': 1.46,
            'mass_2': 1.27,
            'detectors': ['H1', 'L1', 'V1']
        }
    }
    
    def __init__(self, cache_dir: str = 'data/raw/'):
        if not GWPY_AVAILABLE:
            raise ImportError(
                "gwpy is required for data retrieval. "
                "Install with: pip install gwpy\n"
                "Or install all dependencies: pip install -r requirements.txt"
            )
        
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        
        # Load event catalog
        self.event_catalog = self.GWTC3_EVENTS
    
    def download_event_data(
        self,
        event_name: str,
        detectors: Optional[List[str]] = None,
        segment_duration: int = 32,
        sample_rate: int = 4096
    ) -> Dict[str, 'TimeSeries']:
        """
        Download strain data for a specific GW event.
        
        Parameters
        ----------
        event_name : str
            Event name (e.g., 'GW150914')
        detectors : list of str, optional
            Detector names (e.g., ['H1', 'L1']). If None, uses all available.
        segment_duration : int, optional
            Duration of data segment in seconds (default: 32)
        sample_rate : int, optional
            Sample rate in Hz (default: 4096)
            
        Returns
        -------
        dict
            Dictionary mapping detector name to TimeSeries object
        """
        if event_name not in self.event_catalog:
            raise ValueError(f"Unknown event: {event_name}")
        
        event_info = self.event_catalog[event_name]
        gps_time = event_info['gps_time']
        
        if detectors is None:
            detectors = event_info['detectors']
        
        data = {}
        
        for det in detectors:
            print(f"Downloading {event_name} from {det}...")
            
            # Check cache
            cache_file = os.path.join(
                self.cache_dir,
                f"{event_name}_{det}_{segment_duration}s.hdf5"
            )
            
            if os.path.exists(cache_file):
                print(f"  Loading from cache: {cache_file}")
                strain = TimeSeries.read(cache_file, format='hdf5')
            else:
                try:
                    # Download from GWOSC
                    strain = TimeSeries.fetch_open_data(
                        det,
                        gps_time - segment_duration // 2,
                        gps_time + segment_duration // 2,
                        sample_rate=sample_rate,
                        cache=True
                    )
                    
                    # Save to cache
                    strain.write(cache_file, format='hdf5', overwrite=True)
                    print(f"  Saved to cache: {cache_file}")
                    
                except Exception as e:
                    print(f"  Error downloading {det}: {e}")
                    continue
            
            data[det] = strain
        
        return data
    
    def whiten_data(
        self,
        strain: 'TimeSeries',
        fft_duration: float = 4.0,
        overlap: float = 2.0
    ) -> 'TimeSeries':
        """
        Whiten strain data using Welch's method for PSD estimation.
        
        Parameters
        ----------
        strain : TimeSeries
            Raw strain data
        fft_duration : float, optional
            Duration of FFT segments for PSD estimation (default: 4.0 s)
        overlap : float, optional
            Overlap between segments (default: 2.0 s)
            
        Returns
        -------
        TimeSeries
            Whitened strain data
        """
        # Estimate PSD using Welch's method
        psd = strain.psd(
            fftlength=fft_duration,
            overlap=overlap,
            method='median'
        )
        
        # Whiten
        whitened = strain.whiten(
            fftlength=fft_duration,
            overlap=overlap,
            asd=psd.sqrt()
        )
        
        return whitened
    
    def bandpass_filter(
        self,
        strain: 'TimeSeries',
        flow: float = 30.0,
        fhigh: float = 500.0,
        order: int = 8
    ) -> 'TimeSeries':
        """
        Apply bandpass filter to strain data.
        
        Parameters
        ----------
        strain : TimeSeries
            Strain data
        flow : float, optional
            Low-frequency cutoff in Hz (default: 30.0)
        fhigh : float, optional
            High-frequency cutoff in Hz (default: 500.0)
        order : int, optional
            Filter order (default: 8)
            
        Returns
        -------
        TimeSeries
            Filtered strain data
        """
        # Design Butterworth bandpass filter
        bp = filter_design.bandpass(flow, fhigh, strain.sample_rate, type='butter', order=order)
        
        # Apply filter
        filtered = strain.filter(bp, filtfilt=True)
        
        return filtered
    
    def remove_glitches(
        self,
        strain: 'TimeSeries',
        threshold: float = 20.0
    ) -> 'TimeSeries':
        """
        Remove transient glitches using gating.
        
        Parameters
        ----------
        strain : TimeSeries
            Strain data
        threshold : float, optional
            SNR threshold for glitch identification (default: 20.0)
            
        Returns
        -------
        TimeSeries
            Glitch-removed strain data
        """
        # Identify high-SNR transients
        whitened = self.whiten_data(strain)
        
        # Find peaks above threshold
        glitch_times = whitened.times[np.abs(whitened.value) > threshold]
        
        if len(glitch_times) > 0:
            print(f"  Found {len(glitch_times)} potential glitches")
            
            # Gate out glitches (zero out ±0.5s around each glitch)
            cleaned = strain.copy()
            for t_glitch in glitch_times:
                mask = np.abs(strain.times.value - t_glitch.value) < 0.5
                # Apply Tukey window to smoothly zero out
                window = np.ones(np.sum(mask))
                window[:len(window)//10] = np.blackman(2 * len(window)//10)[:len(window)//10]
                window[-len(window)//10:] = np.blackman(2 * len(window)//10)[-len(window)//10:]
                cleaned.value[mask] *= (1 - window)
            
            return cleaned
        else:
            return strain
    
    def preprocess_pipeline(
        self,
        strain: 'TimeSeries',
        flow: float = 30.0,
        fhigh: float = 500.0,
        whiten: bool = True,
        remove_glitches: bool = True
    ) -> 'TimeSeries':
        """
        Complete preprocessing pipeline.
        
        Steps:
        1. Bandpass filter
        2. Remove glitches (optional)
        3. Whiten (optional)
        
        Parameters
        ----------
        strain : TimeSeries
            Raw strain data
        flow : float, optional
            Low-frequency cutoff (default: 30.0 Hz)
        fhigh : float, optional
            High-frequency cutoff (default: 500.0 Hz)
        whiten : bool, optional
            Apply whitening (default: True)
        remove_glitches : bool, optional
            Remove transient glitches (default: True)
            
        Returns
        -------
        TimeSeries
            Preprocessed strain data
        """
        print("Preprocessing pipeline:")
        
        # 1. Bandpass filter
        print("  1. Bandpass filtering...")
        processed = self.bandpass_filter(strain, flow, fhigh)
        
        # 2. Remove glitches
        if remove_glitches:
            print("  2. Removing glitches...")
            processed = self.remove_glitches(processed)
        
        # 3. Whiten
        if whiten:
            print("  3. Whitening...")
            processed = self.whiten_data(processed)
        
        print("  Done!")
        return processed
    
    def save_processed_data(
        self,
        data: Dict[str, 'TimeSeries'],
        output_dir: str = 'data/processed/',
        event_name: str = 'event'
    ):
        """
        Save preprocessed data to HDF5 files.
        
        Parameters
        ----------
        data : dict
            Dictionary mapping detector names to TimeSeries
        output_dir : str, optional
            Output directory
        event_name : str, optional
            Event name for filename
        """
        os.makedirs(output_dir, exist_ok=True)
        
        for det, strain in data.items():
            output_file = os.path.join(output_dir, f"{event_name}_{det}_processed.hdf5")
            strain.write(output_file, format='hdf5', overwrite=True)
            print(f"Saved {det} data to {output_file}")
    
    def load_processed_data(
        self,
        event_name: str,
        detector: str,
        data_dir: str = 'data/processed/'
    ) -> 'TimeSeries':
        """
        Load preprocessed data from HDF5.
        
        Parameters
        ----------
        event_name : str
            Event name
        detector : str
            Detector name (e.g., 'H1')
        data_dir : str, optional
            Data directory
            
        Returns
        -------
        TimeSeries
            Preprocessed strain data
        """
        filepath = os.path.join(data_dir, f"{event_name}_{detector}_processed.hdf5")
        
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        
        return TimeSeries.read(filepath, format='hdf5')
    
    def get_event_info(self, event_name: str) -> Dict:
        """Get event metadata."""
        if event_name not in self.event_catalog:
            raise ValueError(f"Unknown event: {event_name}")
        return self.event_catalog[event_name]


def download_all_gwtc3_events(
    output_dir: str = 'data/raw/',
    events: Optional[List[str]] = None
):
    """
    Download all GWTC-3 events (or subset).
    
    Parameters
    ----------
    output_dir : str, optional
        Output directory
    events : list of str, optional
        Specific events to download. If None, downloads all.
    """
    retriever = LIGODataRetriever(cache_dir=output_dir)
    
    if events is None:
        events = list(retriever.event_catalog.keys())
    
    print(f"Downloading {len(events)} events...")
    print("=" * 60)
    
    for event_name in events:
        print(f"\nEvent: {event_name}")
        print("-" * 60)
        
        try:
            data = retriever.download_event_data(event_name, segment_duration=32)
            print(f"✓ Successfully downloaded {event_name}")
        except Exception as e:
            print(f"✗ Error downloading {event_name}: {e}")
    
    print("\n" + "=" * 60)
    print("Download complete!")


if __name__ == "__main__":
    # Example usage
    print("LIGO Data Retrieval - Example")
    print("=" * 60)
    
    # Initialize retriever
    retriever = LIGODataRetriever()
    
    # Download GW150914
    print("\nDownloading GW150914...")
    data = retriever.download_event_data('GW150914', segment_duration=32)
    
    # Preprocess H1 data
    if 'H1' in data:
        print("\nPreprocessing H1 data...")
        strain_h1 = data['H1']
        processed_h1 = retriever.preprocess_pipeline(strain_h1)
        
        print(f"\nOriginal strain: {strain_h1}")
        print(f"Processed strain: {processed_h1}")
        
        # Save processed data
        retriever.save_processed_data(
            {'H1': processed_h1},
            output_dir='data/processed/',
            event_name='GW150914'
        )
