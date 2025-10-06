"""
Plotting module for HP-POX model.

Generates all required plots in Virtuhcon style.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Dict, Any, List
import matplotlib.patches as mpatches

# Set style
plt.style.use('default')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['grid.alpha'] = 0.3


class Plotting:
    """Plotting utilities for HP-POX model."""
    
    def __init__(self):
        """Initialize plotting utilities."""
        self.colors = {
            'H2': '#FF6B6B',      # Red
            'CO': '#4ECDC4',      # Teal
            'CO2': '#45B7D1',     # Blue
            'CH4': '#96CEB4',     # Green
            'H2O': '#FFEAA7',     # Yellow
            'O2': '#DDA0DD',      # Plum
            'N2': '#98D8C8'       # Mint
        }
    
    def generate_axial_profiles(self, pfr_result: Dict[str, Any], output_dir: Path) -> None:
        """Generate 4-panel axial profiles plot."""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('HP-POX Axial Profiles', fontsize=16, fontweight='bold')
        
        z = pfr_result['z_m']
        T = pfr_result['temperature_K'] - 273.15  # Convert to °C
        P = pfr_result['pressure_Pa'] / 1e5  # Convert to bar
        u = pfr_result['velocity_m_s']
        residence_time = pfr_result['residence_time_s']
        
        # Panel 1: Temperature
        ax1 = axes[0, 0]
        ax1.plot(z, T, 'b-', linewidth=2, label='Temperature')
        ax1.axvline(x=pfr_result['ignition_length_peak_m'], color='r', linestyle='--', 
                   label=f'Peak dT/dz: {pfr_result["ignition_length_peak_m"]:.3f} m')
        ax1.axvline(x=pfr_result['ignition_length_threshold_m'], color='orange', linestyle=':', 
                   label=f'T threshold: {pfr_result["ignition_length_threshold_m"]:.3f} m')
        ax1.set_xlabel('Axial Position (m)')
        ax1.set_ylabel('Temperature (°C)')
        ax1.set_title('Temperature Profile')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Panel 2: Species
        ax2 = axes[0, 1]
        key_species = ['H2', 'CO', 'CO2', 'CH4', 'H2O', 'O2']
        
        for species in key_species:
            if species in pfr_result['species_names']:
                idx = pfr_result['species_names'].index(species)
                ax2.plot(z, pfr_result['mole_fractions'][idx, :] * 100, 
                        color=self.colors.get(species, '#000000'), 
                        linewidth=2, label=species)
        
        ax2.set_xlabel('Axial Position (m)')
        ax2.set_ylabel('Mole Fraction (%)')
        ax2.set_title('Species Profiles')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Panel 3: Pressure
        ax3 = axes[1, 0]
        ax3.plot(z, P, 'b-', linewidth=2, label='Pressure')
        ax3.set_xlabel('Axial Position (m)')
        ax3.set_ylabel('Pressure (bar)')
        ax3.set_title('Pressure Profile')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Panel 4: Residence time
        ax4 = axes[1, 1]
        ax4.plot(z, residence_time, 'b-', linewidth=2, label='Residence Time')
        ax4.set_xlabel('Axial Position (m)')
        ax4.set_ylabel('Residence Time (s)')
        ax4.set_title('Residence Time Profile')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'axial_profiles_plot.png', bbox_inches='tight')
        plt.close()
    
    def generate_ignition_markers(self, pfr_result: Dict[str, Any], output_dir: Path) -> None:
        """Generate temperature plot with ignition markers."""
        fig, ax = plt.subplots(figsize=(12, 8))
        
        z = pfr_result['z_m']
        T = pfr_result['temperature_K'] - 273.15  # Convert to °C
        
        # Plot temperature
        ax.plot(z, T, 'b-', linewidth=3, label='Temperature')
        
        # Add ignition markers
        ign_peak = pfr_result['ignition_length_peak_m']
        ign_thresh = pfr_result['ignition_length_threshold_m']
        
        ax.axvline(x=ign_peak, color='r', linestyle='--', linewidth=2,
                  label=f'Peak dT/dz: {ign_peak:.3f} m')
        ax.axvline(x=ign_thresh, color='orange', linestyle=':', linewidth=2,
                  label=f'T threshold: {ign_thresh:.3f} m')
        
        # Highlight ignition zone
        ax.axvspan(0, ign_peak, alpha=0.2, color='red', label='Ignition Zone')
        
        ax.set_xlabel('Axial Position (m)')
        ax.set_ylabel('Temperature (°C)')
        ax.set_title('Temperature Profile with Ignition Markers')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'ignition_length_markers.png', bbox_inches='tight')
        plt.close()
    
    def generate_outlet_composition(self, kpi_result: Dict[str, Any], output_dir: Path) -> None:
        """Generate outlet composition table plot."""
        # This is handled by the bar plot, but we can add a text table here if needed
        pass
    
    def generate_outlet_barplot(self, kpi_result: Dict[str, Any], output_dir: Path) -> None:
        """Generate outlet composition bar chart."""
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Get dry composition data
        dry_comp = kpi_result['dry_composition_pct']
        targets = kpi_result['targets']['dry_composition_pct']
        
        species = list(dry_comp.keys())
        values = list(dry_comp.values())
        target_values = [targets.get(s, 0.0) for s in species]
        
        # Create bars
        x = np.arange(len(species))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, values, width, label='Actual', 
                      color=[self.colors.get(s, '#000000') for s in species], alpha=0.8)
        bars2 = ax.bar(x + width/2, target_values, width, label='Target', 
                      color=[self.colors.get(s, '#000000') for s in species], alpha=0.4)
        
        # Add value labels
        for i, (bar, value) in enumerate(zip(bars1, values)):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{value:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        for i, (bar, value) in enumerate(zip(bars2, target_values)):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{value:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        # Calculate and display errors
        errors = kpi_result['errors_pp']
        error_text = "Errors (pp): " + ", ".join([f"{s}: {e:.1f}" for s, e in errors.items()])
        
        ax.set_xlabel('Species')
        ax.set_ylabel('Mole Fraction (%)')
        ax.set_title('Outlet Composition (Dry Basis)')
        ax.set_xticks(x)
        ax.set_xticklabels(species)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add error text
        ax.text(0.02, 0.98, error_text, transform=ax.transAxes, 
               verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(output_dir / 'outlet_barplot.png', bbox_inches='tight')
        plt.close()
