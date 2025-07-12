import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, optimize, interpolate
from scipy.ndimage import gaussian_filter1d
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

class TelescopeDataCleaner:
    """
    Advanced telescope data cleaning pipeline for exoplanet atmospheric detection
    """
    
    def __init__(self):
        self.wavelength_range = None
        self.reference_spectrum = None
        self.calibration_params = {}
        
    def remove_cosmic_rays(self, spectrum, threshold=5):
        """Remove cosmic ray hits using median filtering"""
        median_filtered = signal.medfilt(spectrum, kernel_size=5)
        deviation = np.abs(spectrum - median_filtered)
        mad = np.median(deviation)
        
        # Identify cosmic ray pixels
        cosmic_mask = deviation > threshold * mad
        
        # Replace cosmic rays with interpolated values
        clean_spectrum = spectrum.copy()
        if np.any(cosmic_mask):
            indices = np.arange(len(spectrum))
            clean_spectrum[cosmic_mask] = np.interp(
                indices[cosmic_mask], 
                indices[~cosmic_mask], 
                spectrum[~cosmic_mask]
            )
        
        return clean_spectrum, cosmic_mask
    
    def correct_instrumental_response(self, spectrum, wavelength):
        """Correct for instrumental response and detector non-linearity"""
        # Simulate instrumental response correction
        # In practice, this would use calibration lamp data
        
        # Wavelength-dependent efficiency correction
        efficiency = 1.0 - 0.1 * np.exp(-(wavelength - 1000)**2 / (2 * 200**2))
        corrected_spectrum = spectrum / efficiency
        
        # Non-linearity correction (quadratic)
        linearity_factor = 1.0 + 1e-6 * spectrum**2
        corrected_spectrum = corrected_spectrum / linearity_factor
        
        return corrected_spectrum
    
    def remove_telluric_absorption(self, spectrum, wavelength, telluric_model=None):
        """Remove Earth's atmospheric absorption lines"""
        if telluric_model is None:
            # Create synthetic telluric model
            telluric_lines = [760, 820, 940, 1130, 1400]  # nm
            telluric_model = np.ones_like(wavelength)
            
            for line_center in telluric_lines:
                if wavelength.min() <= line_center <= wavelength.max():
                    absorption = 0.3 * np.exp(-(wavelength - line_center)**2 / (2 * 5**2))
                    telluric_model -= absorption
        
        # Divide out telluric absorption
        corrected_spectrum = spectrum / np.maximum(telluric_model, 0.1)
        
        return corrected_spectrum, telluric_model
    
    def denoise_spectrum(self, spectrum, method='savgol', **kwargs):
        """Apply denoising using various methods"""
        if method == 'savgol':
            window_length = kwargs.get('window_length', 15)
            polyorder = kwargs.get('polyorder', 3)
            return signal.savgol_filter(spectrum, window_length, polyorder)
        
        elif method == 'gaussian':
            sigma = kwargs.get('sigma', 1.0)
            return gaussian_filter1d(spectrum, sigma)
        
        elif method == 'median':
            kernel_size = kwargs.get('kernel_size', 5)
            return signal.medfilt(spectrum, kernel_size)
        
        else:
            return spectrum

class AtmosphericSignatureDetector:
    """
    Detect and analyze atmospheric spectral signatures of exoplanets
    """
    
    def __init__(self):
        self.molecular_templates = self._load_molecular_templates()
        self.detection_threshold = 3.0  # sigma
        
    def _load_molecular_templates(self):
        """Load molecular absorption templates"""
        # Simplified molecular templates (wavelength in nm, relative absorption)
        templates = {
            'H2O': {
                'wavelengths': np.array([940, 1130, 1400, 1870]),
                'depths': np.array([0.02, 0.015, 0.025, 0.03]),
                'widths': np.array([10, 15, 20, 25])
            },
            'CO2': {
                'wavelengths': np.array([1430, 1600, 2010]),
                'depths': np.array([0.01, 0.008, 0.012]),
                'widths': np.array([8, 12, 15])
            },
            'CH4': {
                'wavelengths': np.array([1660, 2300]),
                'depths': np.array([0.005, 0.008]),
                'widths': np.array([20, 30])
            },
            'O2': {
                'wavelengths': np.array([760, 1270]),
                'depths': np.array([0.03, 0.01]),
                'widths': np.array([5, 8])
            },
            'O3': {
                'wavelengths': np.array([600, 1000]),
                'depths': np.array([0.02, 0.015]),
                'widths': np.array([50, 30])
            }
        }
        return templates
    
    def generate_molecular_template(self, molecule, wavelength_grid):
        """Generate molecular absorption template on given wavelength grid"""
        if molecule not in self.molecular_templates:
            return np.zeros_like(wavelength_grid)
        
        template_data = self.molecular_templates[molecule]
        template = np.ones_like(wavelength_grid)
        
        for center, depth, width in zip(
            template_data['wavelengths'],
            template_data['depths'],
            template_data['widths']
        ):
            if wavelength_grid.min() <= center <= wavelength_grid.max():
                absorption = depth * np.exp(-(wavelength_grid - center)**2 / (2 * width**2))
                template -= absorption
        
        return template
    
    def cross_correlate_molecules(self, spectrum, wavelength, molecules=None):
        """Cross-correlate spectrum with molecular templates"""
        if molecules is None:
            molecules = list(self.molecular_templates.keys())
        
        results = {}
        
        for molecule in molecules:
            template = self.generate_molecular_template(molecule, wavelength)
            
            # Normalize both spectrum and template
            spectrum_norm = (spectrum - np.mean(spectrum)) / np.std(spectrum)
            template_norm = (template - np.mean(template)) / np.std(template)
            
            # Cross-correlation
            correlation = signal.correlate(spectrum_norm, template_norm, mode='same')
            correlation = correlation / len(spectrum)
            
            # Find peak correlation
            max_corr = np.max(np.abs(correlation))
            max_idx = np.argmax(np.abs(correlation))
            
            # Calculate significance (simplified)
            noise_level = np.std(correlation[correlation < np.percentile(correlation, 90)])
            significance = max_corr / noise_level if noise_level > 0 else 0
            
            results[molecule] = {
                'correlation': correlation,
                'max_correlation': max_corr,
                'significance': significance,
                'detected': significance > self.detection_threshold
            }
        
        return results
    
    def fit_molecular_abundances(self, spectrum, wavelength, molecules):
        """Fit molecular abundances using least squares"""
        # Create design matrix with molecular templates
        templates = []
        for molecule in molecules:
            template = self.generate_molecular_template(molecule, wavelength)
            templates.append(1 - template)  # Convert to absorption
        
        templates = np.array(templates).T
        
        # Add continuum term
        continuum = np.ones((len(wavelength), 1))
        design_matrix = np.hstack([continuum, templates])
        
        # Solve least squares
        try:
            coefficients, residuals, rank, s = np.linalg.lstsq(
                design_matrix, spectrum, rcond=None
            )
            
            continuum_level = coefficients[0]
            abundances = coefficients[1:]
            
            # Calculate fit quality
            model_spectrum = design_matrix @ coefficients
            chi_squared = np.sum((spectrum - model_spectrum)**2)
            reduced_chi_squared = chi_squared / (len(spectrum) - len(coefficients))
            
            return {
                'abundances': dict(zip(molecules, abundances)),
                'continuum': continuum_level,
                'model_spectrum': model_spectrum,
                'chi_squared': reduced_chi_squared,
                'residuals': spectrum - model_spectrum
            }
            
        except np.linalg.LinAlgError:
            return None

class ExoplanetSpectralPipeline:
    """
    Complete pipeline for exoplanet atmospheric analysis
    """
    
    def __init__(self):
        self.cleaner = TelescopeDataCleaner()
        self.detector = AtmosphericSignatureDetector()
        self.results = {}
        
    def process_spectrum(self, raw_spectrum, wavelength, plot_results=True):
        """Complete processing pipeline"""
        print("Processing exoplanet spectrum...")
        
        # Step 1: Remove cosmic rays
        spectrum, cosmic_mask = self.cleaner.remove_cosmic_rays(raw_spectrum)
        print(f"Removed {np.sum(cosmic_mask)} cosmic ray hits")
        
        # Step 2: Correct instrumental response
        spectrum = self.cleaner.correct_instrumental_response(spectrum, wavelength)
        print("Applied instrumental response correction")
        
        # Step 3: Remove telluric absorption
        spectrum, telluric_model = self.cleaner.remove_telluric_absorption(spectrum, wavelength)
        print("Removed telluric absorption")
        
        # Step 4: Denoise
        spectrum = self.cleaner.denoise_spectrum(spectrum, method='savgol')
        print("Applied denoising")
        
        # Step 5: Detect molecular signatures
        molecules = ['H2O', 'CO2', 'CH4', 'O2', 'O3']
        correlations = self.detector.cross_correlate_molecules(spectrum, wavelength, molecules)
        
        # Step 6: Fit molecular abundances
        detected_molecules = [mol for mol, data in correlations.items() if data['detected']]
        if detected_molecules:
            abundance_fit = self.detector.fit_molecular_abundances(
                spectrum, wavelength, detected_molecules
            )
        else:
            abundance_fit = None
        
        # Store results
        self.results = {
            'processed_spectrum': spectrum,
            'wavelength': wavelength,
            'cosmic_rays': cosmic_mask,
            'telluric_model': telluric_model,
            'correlations': correlations,
            'abundance_fit': abundance_fit,
            'detected_molecules': detected_molecules
        }
        
        if plot_results:
            self.plot_results()
        
        return self.results
    
    def plot_results(self):
        """Plot analysis results"""
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
        # Plot 1: Raw vs processed spectrum
        ax1 = axes[0, 0]
        ax1.plot(self.results['wavelength'], self.results['processed_spectrum'], 
                'b-', alpha=0.8, label='Processed')
        ax1.set_xlabel('Wavelength (nm)')
        ax1.set_ylabel('Flux')
        ax1.set_title('Processed Spectrum')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Plot 2: Molecular detections
        ax2 = axes[0, 1]
        molecules = list(self.results['correlations'].keys())
        significances = [self.results['correlations'][mol]['significance'] for mol in molecules]
        colors = ['red' if sig > 3 else 'blue' for sig in significances]
        
        bars = ax2.bar(molecules, significances, color=colors, alpha=0.7)
        ax2.axhline(y=3, color='red', linestyle='--', alpha=0.7, label='Detection threshold')
        ax2.set_ylabel('Detection Significance (σ)')
        ax2.set_title('Molecular Detection Results')
        ax2.legend()
        ax2.tick_params(axis='x', rotation=45)
        
        # Plot 3: Best-fit model (if available)
        ax3 = axes[1, 0]
        if self.results['abundance_fit'] is not None:
            ax3.plot(self.results['wavelength'], self.results['processed_spectrum'], 
                    'k-', alpha=0.8, label='Data')
            ax3.plot(self.results['wavelength'], self.results['abundance_fit']['model_spectrum'], 
                    'r--', alpha=0.8, label='Model fit')
            ax3.set_xlabel('Wavelength (nm)')
            ax3.set_ylabel('Flux')
            ax3.set_title('Molecular Abundance Fit')
            ax3.legend()
        else:
            ax3.text(0.5, 0.5, 'No significant\nmolecular detections', 
                    ha='center', va='center', transform=ax3.transAxes)
            ax3.set_title('Molecular Abundance Fit')
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Residuals
        ax4 = axes[1, 1]
        if self.results['abundance_fit'] is not None:
            residuals = self.results['abundance_fit']['residuals']
            ax4.plot(self.results['wavelength'], residuals, 'g-', alpha=0.8)
            ax4.axhline(y=0, color='black', linestyle='-', alpha=0.5)
            ax4.set_xlabel('Wavelength (nm)')
            ax4.set_ylabel('Residuals')
            ax4.set_title(f"Fit Residuals (χ²ᵣ = {self.results['abundance_fit']['chi_squared']:.2f})")
        else:
            ax4.text(0.5, 0.5, 'No fit available', 
                    ha='center', va='center', transform=ax4.transAxes)
            ax4.set_title('Fit Residuals')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        # Print detection summary
        print("\n" + "="*50)
        print("ATMOSPHERIC ANALYSIS SUMMARY")
        print("="*50)
        
        detected = self.results['detected_molecules']
        if detected:
            print(f"Detected molecules: {', '.join(detected)}")
            
            if self.results['abundance_fit'] is not None:
                print("\nMolecular abundances (relative):")
                for mol, abundance in self.results['abundance_fit']['abundances'].items():
                    print(f"  {mol}: {abundance:.4f}")
                
                print(f"\nFit quality (χ²ᵣ): {self.results['abundance_fit']['chi_squared']:.3f}")
        else:
            print("No significant molecular detections found.")
        
        print("="*50)

# Demo function to generate synthetic exoplanet data
def generate_synthetic_exoplanet_data(wavelength_range=(500, 2500), n_points=2000):
    """Generate synthetic exoplanet transit spectrum with realistic noise"""
    wavelength = np.linspace(wavelength_range[0], wavelength_range[1], n_points)
    
    # Base stellar spectrum (simplified blackbody)
    T_star = 5800  # K
    h, c, k = 6.626e-34, 2.998e8, 1.381e-23
    wl_m = wavelength * 1e-9
    
    # Planck function
    flux = (2 * h * c**2 / wl_m**5) / (np.exp(h * c / (wl_m * k * T_star)) - 1)
    flux = flux / np.max(flux)
    
    # Add atmospheric absorption
    detector = AtmosphericSignatureDetector()
    
    # Water vapor
    h2o_template = detector.generate_molecular_template('H2O', wavelength)
    flux *= (1 - 0.003 * (1 - h2o_template))  # 0.3% absorption depth
    
    # Carbon dioxide  
    co2_template = detector.generate_molecular_template('CO2', wavelength)
    flux *= (1 - 0.002 * (1 - co2_template))  # 0.2% absorption depth
    
    # Add realistic noise
    # Photon noise
    photon_noise = np.random.poisson(flux * 10000) / 10000 - flux
    
    # Systematic noise
    systematic_noise = 0.0005 * np.sin(2 * np.pi * wavelength / 100)
    
    # Random noise
    random_noise = np.random.normal(0, 0.0002, len(wavelength))
    
    # Add cosmic rays
    n_cosmic = 10
    cosmic_indices = np.random.choice(len(wavelength), n_cosmic, replace=False)
    cosmic_strength = np.random.exponential(0.01, n_cosmic)
    
    noisy_flux = flux + photon_noise + systematic_noise + random_noise
    noisy_flux[cosmic_indices] += cosmic_strength
    
    return wavelength, noisy_flux

# Example usage
if __name__ == "__main__":
    print("Exoplanet Atmospheric Spectral Analysis Pipeline")
    print("=" * 50)
    
    # Generate synthetic data
    print("Generating synthetic exoplanet transit spectrum...")
    wavelength, raw_spectrum = generate_synthetic_exoplanet_data()
    
    # Initialize pipeline
    pipeline = ExoplanetSpectralPipeline()
    
    # Process the spectrum
    results = pipeline.process_spectrum(raw_spectrum, wavelength)
    
    print("\nPipeline processing complete!")
    print(f"Analysis covered {len(wavelength)} spectral points")
    print(f"Wavelength range: {wavelength.min():.1f} - {wavelength.max():.1f} nm")