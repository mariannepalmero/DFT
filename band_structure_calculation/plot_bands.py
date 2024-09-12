import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import numpy as np

# Set global matplotlib parameters for consistent plot appearance
plt.rcParams["figure.dpi"]=150 # Set resolution for the figure (DPI)
plt.rcParams["figure.facecolor"]="white" # Set figure background color
plt.rcParams["figure.figsize"]=(8, 6) # Set figure size (width, height in inches)

# Fermi energy (in eV) - adjust this value based on your specific SCF calculation output
fermi = 4.3104

# Function to plot band structure without spin component specification
def plot_band_structure_without_spin(data_file, fermi_energy):
    """
    Plots the band structure without spin component specification.

    Parameters:
    - data_file: str, path to the input file containing band structure data
    - fermi_energy: float, the Fermi energy level for the material (eV)
    """
    # Load data from file
    data = np.loadtxt(data_file)
    k_points = np.unique(data[:, 0]) # Extract unique k-points
    bands = np.reshape(data[:, 1], (-1, len(k_points)))  # Reshape data into bands

    # Plot each band
    for band in range(len(bands)):
        plt.plot(k_points, (bands[band, :]) - fermi_energy, linewidth=1, alpha=0.5, color='k')

    plt.xlim(min(k_points), max(k_points))  # Set x-axis limits to k-point range

# Function to plot band structure with spin component specification
def plot_band_structure_with_spin(data_file_up, data_file_down, fermi_energy):
    """
    Plots the band structure with spin component specification for spin-up and spin-down.

    Parameters:
    - data_file_up: str, path to the input file containing spin-up band data
    - data_file_down: str, path to the input file containing spin-down band data
    - fermi_energy: float, the Fermi energy level for the material (eV)
    """
    # Load spin-up and spin-down data
    data_up = np.loadtxt(data_file_up)
    data_down = np.loadtxt(data_file_down)

    k_points = np.unique(data_up[:, 0])# Extract unique k-points
    bands_up = np.reshape(data_up[:, 1], (-1, len(k_points)))   # Reshape spin-up data into bands
    bands_down = np.reshape(data_down[:, 1], (-1, len(k_points)))  # Reshape spin-down data into bands

    # Plot each spin-up band in red and spin-down band in blue
    for band in range(len(bands_up)):
        plt.plot(k_points, (bands_up[band, :]) - fermi_energy, linewidth=1, alpha=0.5, color='red')
        plt.plot(k_points, (bands_down[band, :]) - fermi_energy, linewidth=1, alpha=0.5, color='blue')

    plt.xlim(min(k_points), max(k_points))  # Set x-axis limits to k-point range

# Main plotting configuration
def plot_band_structure(fermi_energy):
    """
    Configures and displays the plot for band structure with high symmetry points and labels.

    Parameters:
    - fermi_energy: float, the Fermi energy level (eV)
    """
    # Plot horizontal Fermi level line at E-Ef = 0
    plt.axhline(0, linestyle=(0, (5, 5)), linewidth=0.75, color='k', alpha=0.5)

    # Add vertical lines for high symmetry k-points (modify based on the specific path in your data)
    high_symmetry_kpoints = [0.5774, 0.9107, 1.5774]
    for point in high_symmetry_kpoints:
        plt.axvline(point, linewidth=0.75, color='k', alpha=0.5)

    # Customize x-ticks with k-point labels (change these based on your path)
    plt.xticks(ticks=[0, 0.5774, 0.9107, 1.5774], labels=['$\Gamma$', 'M', 'K', '$\Gamma$'])

    # Set y-label for energy relative to Fermi level
    plt.ylabel("E - E$_{F}$ (eV)")

    # Set y-axis range (adjust based on data)
    plt.ylim(-6, 6)

    # Save plot to PNG file
    plt.savefig("plot_surf-1.png", format='png')
    plt.show()

# Uncomment the appropriate function depending on whether spin is considered
''''
# Example without spin component (specify the correct data file):
# plot_band_structure_without_spin('surf-1.dat.gnu', fermi)

# Example with spin component (specify the correct data files for spin-up and spin-down):
# plot_band_structure_with_spin('surf-1_up.dat.gnu', 'surf-1_down.dat.gnu', fermi)
'''

# Finally, configure and display the plot
plot_band_structure(fermi)
