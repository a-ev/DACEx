import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable

# Set larger font sizes globally
plt.rcParams.update({'font.size': 18})

# Load 2D data
data = np.loadtxt('data_out/comparison_data_2d.txt', skiprows=1)

# Extract columns - NEW FORMAT with approximation values
# Columns: x, y, true, cheb_val_4, taylor_val_4, cheb_error_4, taylor_error_4, cheb_val_6, taylor_val_6, cheb_error_6, taylor_error_6, cheb_val_8, taylor_val_8, cheb_error_8, taylor_error_8
x_vals = data[:, 0]
y_vals = data[:, 1]
true_vals = data[:, 2]
cheb_err_4 = data[:, 5]
taylor_err_4 = data[:, 6]
cheb_err_6 = data[:, 9]
taylor_err_6 = data[:, 10]
cheb_err_8 = data[:, 13]
taylor_err_8 = data[:, 14]

# Determine grid shape
n_unique_x = len(np.unique(x_vals))
n_unique_y = len(np.unique(y_vals))
print(f"Grid shape: {n_unique_x} x {n_unique_y}")

# Reshape data into 2D grids
x_grid = x_vals.reshape(n_unique_y, n_unique_x)
y_grid = y_vals.reshape(n_unique_y, n_unique_x)
true_grid = true_vals.reshape(n_unique_y, n_unique_x)
cheb_err_4_grid = cheb_err_4.reshape(n_unique_y, n_unique_x)
taylor_err_4_grid = taylor_err_4.reshape(n_unique_y, n_unique_x)
cheb_err_6_grid = cheb_err_6.reshape(n_unique_y, n_unique_x)
taylor_err_6_grid = taylor_err_6.reshape(n_unique_y, n_unique_x)
cheb_err_8_grid = cheb_err_8.reshape(n_unique_y, n_unique_x)
taylor_err_8_grid = taylor_err_8.reshape(n_unique_y, n_unique_x)

# Create figure with 2x3 subplots (3 orders, 2 methods)
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
# fig.suptitle('2D Gaussian Error Heatmaps: Chebyshev vs Taylor', fontsize=16, fontweight='bold')

orders = [4, 6, 8]
orders_data = [
    (cheb_err_4_grid, taylor_err_4_grid),
    (cheb_err_6_grid, taylor_err_6_grid),
    (cheb_err_8_grid, taylor_err_8_grid)
]

# Find global max across BOTH methods and ALL orders
vmax_global = max([grid[0].max() for grid in orders_data] + [grid[1].max() for grid in orders_data])

# Store all images for colorbars
images = []

for col, order in enumerate(orders):
    cheb_grid, taylor_grid = orders_data[col]
    
    # Chebyshev heatmap (top row)
    ax = axes[0, col]
    im = ax.contourf(x_grid, y_grid, cheb_grid, levels=20, cmap='RdYlBu_r', vmin=0, vmax=vmax_global)
    # Use logarithmically spaced contour levels for better distribution
    cheb_max = cheb_grid.max()
    # Restore 6 Chebyshev contour levels for detail
    cheb_levels = np.logspace(np.log10(cheb_max/100), np.log10(cheb_max), 6)
    contours = ax.contour(x_grid, y_grid, cheb_grid, levels=cheb_levels, colors='black', alpha=0.8, linewidths=1.2)
    # Only label every other contour, and skip near-zero levels to avoid '0.00' labels
    # For Chebyshev 8th order, use exponential notation for labels; for others, use '%.2f'
    label_levels = cheb_levels[::2]
    def smart_fmt(x):
        return f'{x:.1e}' if abs(x) < 0.01 else f'{x:.2f}'
    ax.clabel(contours, levels=label_levels, inline=True, fontsize=14, fmt=smart_fmt, colors='black')
    ax.set_xlabel('x', fontsize=20)
    ax.set_ylabel('y', fontsize=20)
    ax.set_title(f'Chebyshev Order {order}', fontsize=18)
    ax.set_aspect('equal')
    images.append(im)
    print(f"Chebyshev Order {order}: Max = {cheb_grid.max():.4f}")
    
    # Taylor heatmap (bottom row)
    ax = axes[1, col]
    im = ax.contourf(x_grid, y_grid, taylor_grid, levels=20, cmap='RdYlBu_r', vmin=0, vmax=vmax_global)
    # Use logarithmically spaced contour levels for better distribution
    taylor_max = taylor_grid.max()
    taylor_levels = np.logspace(np.log10(taylor_max/100), np.log10(taylor_max), 6)
    contours = ax.contour(x_grid, y_grid, taylor_grid, levels=taylor_levels, colors='black', alpha=0.8, linewidths=1.2)
    ax.clabel(contours, inline=True, fontsize=14, fmt='%.2f', colors='black')
    ax.set_xlabel('x', fontsize=20)
    ax.set_ylabel('y', fontsize=20)
    ax.set_title(f'Taylor Order {order}', fontsize=18)
    ax.set_aspect('equal')
    images.append(im)
    print(f"Taylor Order {order}: Max = {taylor_grid.max():.4f}")

# Adjust spacing to prevent overlap between subplot titles and axis labels
plt.subplots_adjust(hspace=0.35, wspace=0.3)

# Add a single colorbar for all plots using ScalarMappable
norm = mcolors.Normalize(vmin=0, vmax=vmax_global)
sm = ScalarMappable(cmap='RdYlBu_r', norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Error', fontsize=18)

plt.savefig('data_out/gaussian_2d_heatmaps.png', dpi=150, bbox_inches='tight')
print("\nPlot saved as data_out/gaussian_2d_heatmaps.png")
plt.show()

# Print summary statistics
print("\n" + "="*70)
print("2D FULL GRID ERROR SUMMARY")
print("="*70)

for order in [4, 6, 8]:
    if order == 4:
        cheb_g, taylor_g = cheb_err_4_grid, taylor_err_4_grid
    elif order == 6:
        cheb_g, taylor_g = cheb_err_6_grid, taylor_err_6_grid
    else:
        cheb_g, taylor_g = cheb_err_8_grid, taylor_err_8_grid
    
    print(f"\nOrder {order}:")
    print(f"  Chebyshev - Max: {cheb_g.max():.6f}, Mean: {cheb_g.mean():.6f}")
    print(f"  Taylor    - Max: {taylor_g.max():.6f}, Mean: {taylor_g.mean():.6f}")
    
    if cheb_g.max() < taylor_g.max():
        imp = (taylor_g.max() - cheb_g.max()) / taylor_g.max() * 100
        print(f"  Chebyshev better by {imp:.1f}%")
    else:
        imp = (cheb_g.max() - taylor_g.max()) / cheb_g.max() * 100
        print(f"  Taylor better by {imp:.1f}%")

print("="*70)

plt.show()

