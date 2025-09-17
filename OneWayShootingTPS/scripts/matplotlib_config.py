# matplotlib_config.py
import matplotlib.pyplot as plt

colors = {
    'red': '#A61C3C',
    'purple': '#691C77',
    'darkblue': '#4E5593',
    'lightblue': '#6DB1DE',
    'darkgreen': '#82CBAB',
    'lightgreen': '#ACDA67'
}

plt.rcParams.update({
    # 'font.family': 'sans-serif',
    # 'font.sans-serif': ['Arial', 'Helvetica', 'sans-serif'],
    'font.family': 'serif',  # Use a LaTeX-compatible serif font
    'font.serif': ['Computer Modern Roman'],  # Match LaTeX's default font
    'figure.figsize': (10, 6),
    'figure.dpi': 100,
    'text.usetex': True,
    'font.size': 18,
    'axes.labelsize': 20,
    'axes.titlesize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'lines.linewidth': 2,
    'lines.markersize': 8,
    'axes.grid': True,
    'grid.color': 'gray',
    'grid.linestyle': '--',
    'grid.linewidth': 0.5,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'xtick.major.width': 1,
    'ytick.major.width': 1,
    'xtick.minor.size': 4,
    'ytick.minor.size': 4,
    'xtick.minor.width': 0.8,
    'ytick.minor.width': 0.8,
    'axes.prop_cycle': plt.cycler('color', colors)
})
