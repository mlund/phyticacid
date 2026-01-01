#!/usr/bin/env python3
"""
Simple visualization and analysis tools for MD simulation
Requires: matplotlib, numpy
"""

import sys
import os

def plot_xvg(filename, title=None, ylabel=None):
    """Plot XVG file from Gromacs"""
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("Error: matplotlib and numpy required for plotting")
        print("Install with: pip install matplotlib numpy")
        sys.exit(1)

    # Read XVG file
    x_data = []
    y_data = []

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                # Parse title and labels from XVG comments
                if '@    title' in line and title is None:
                    title = line.split('"')[1]
                if '@    yaxis  label' in line and ylabel is None:
                    ylabel = line.split('"')[1]
                continue

            parts = line.split()
            if len(parts) >= 2:
                try:
                    x_data.append(float(parts[0]))
                    y_data.append(float(parts[1]))
                except ValueError:
                    continue

    x_data = np.array(x_data)
    y_data = np.array(y_data)

    # Create plot
    plt.figure(figsize=(10, 6))
    plt.plot(x_data, y_data, linewidth=1.5)
    plt.xlabel('Time (ps)')
    plt.ylabel(ylabel or 'Value')
    plt.title(title or os.path.basename(filename))
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save plot
    output_file = filename.replace('.xvg', '.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Plot saved to {output_file}")

    # Show statistics
    print(f"\nStatistics for {os.path.basename(filename)}:")
    print(f"  Mean: {np.mean(y_data):.3f}")
    print(f"  Std Dev: {np.std(y_data):.3f}")
    print(f"  Min: {np.min(y_data):.3f}")
    print(f"  Max: {np.max(y_data):.3f}")
    print()

def plot_all_xvg_in_directory(directory='analysis'):
    """Plot all XVG files in a directory"""
    if not os.path.exists(directory):
        print(f"Error: Directory '{directory}' not found")
        return

    xvg_files = [f for f in os.listdir(directory) if f.endswith('.xvg')]

    if not xvg_files:
        print(f"No XVG files found in '{directory}'")
        return

    print(f"Found {len(xvg_files)} XVG files in '{directory}'")
    print()

    for xvg_file in xvg_files:
        filepath = os.path.join(directory, xvg_file)
        try:
            plot_xvg(filepath)
        except Exception as e:
            print(f"Error plotting {xvg_file}: {e}")

def create_analysis_summary():
    """Create a summary of analysis results"""
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("Error: matplotlib and numpy required")
        sys.exit(1)

    analysis_dir = 'analysis'
    if not os.path.exists(analysis_dir):
        print(f"Error: '{analysis_dir}' directory not found")
        print("Run analysis.sh first")
        return

    # Create multi-panel figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Phytic Acid MD Simulation - Summary', fontsize=16, fontweight='bold')

    plots = [
        ('rmsd.xvg', 'RMSD', 'RMSD (nm)'),
        ('gyrate.xvg', 'Radius of Gyration', 'Rg (nm)'),
        ('energy.xvg', 'Energy', 'Energy (kJ/mol)'),
        ('temperature.xvg', 'Temperature', 'Temperature (K)'),
        ('pressure.xvg', 'Pressure', 'Pressure (bar)'),
        ('density.xvg', 'Density', 'Density (kg/m³)')
    ]

    for idx, (filename, title, ylabel) in enumerate(plots):
        ax = axes[idx // 3, idx % 3]
        filepath = os.path.join(analysis_dir, filename)

        if os.path.exists(filepath):
            try:
                x_data = []
                y_data = []

                with open(filepath, 'r') as f:
                    for line in f:
                        if line.startswith('#') or line.startswith('@'):
                            continue
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                x_data.append(float(parts[0]))
                                y_data.append(float(parts[1]))
                            except ValueError:
                                continue

                if x_data and y_data:
                    ax.plot(x_data, y_data, linewidth=1)
                    ax.set_xlabel('Time (ps)', fontsize=9)
                    ax.set_ylabel(ylabel, fontsize=9)
                    ax.set_title(title, fontsize=10, fontweight='bold')
                    ax.grid(True, alpha=0.3)
                else:
                    ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                    ax.set_title(title, fontsize=10)
            except Exception as e:
                ax.text(0.5, 0.5, f'Error: {str(e)[:20]}', ha='center', va='center', fontsize=8)
                ax.set_title(title, fontsize=10)
        else:
            ax.text(0.5, 0.5, 'File not found', ha='center', va='center')
            ax.set_title(title, fontsize=10)

    plt.tight_layout()
    output_file = 'analysis_summary.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Summary plot saved to {output_file}")

def main():
    print("=" * 60)
    print("Phytic Acid MD Simulation - Visualization")
    print("=" * 60)
    print()

    if len(sys.argv) > 1:
        if sys.argv[1] == 'summary':
            create_analysis_summary()
        elif sys.argv[1] == 'all':
            plot_all_xvg_in_directory()
        elif os.path.exists(sys.argv[1]):
            plot_xvg(sys.argv[1])
        else:
            print(f"Error: File '{sys.argv[1]}' not found")
    else:
        print("Usage:")
        print("  python visualize.py <file.xvg>     - Plot single XVG file")
        print("  python visualize.py all            - Plot all XVG files in analysis/")
        print("  python visualize.py summary        - Create summary figure")
        print()
        print("Examples:")
        print("  python visualize.py analysis/rmsd.xvg")
        print("  python visualize.py all")
        print("  python visualize.py summary")

if __name__ == '__main__':
    main()
