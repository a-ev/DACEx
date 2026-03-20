import numpy as np
import matplotlib.pyplot as plt

# Match font sizing used in the 2D heatmap plot
plt.rcParams.update({'font.size': 18})

# Load 1D pure data - NEW FORMAT with approximation values
data_1d = np.loadtxt('data_out/comparison_data.txt', skiprows=1)
x_1d = data_1d[:, 0]
true_1d = data_1d[:, 1]
# Columns: x, true, cheb_val_4, taylor_val_4, cheb_err_4, taylor_err_4, cheb_val_6, taylor_val_6, cheb_err_6, taylor_err_6, cheb_val_8, taylor_val_8, cheb_err_8, taylor_err_8
cheb_val_1d_4 = data_1d[:, 2]
taylor_val_1d_4 = data_1d[:, 3]
cheb_err_1d_4 = data_1d[:, 4]
taylor_err_1d_4 = data_1d[:, 5]
cheb_val_1d_6 = data_1d[:, 6]
taylor_val_1d_6 = data_1d[:, 7]
cheb_err_1d_6 = data_1d[:, 8]
taylor_err_1d_6 = data_1d[:, 9]
cheb_val_1d_8 = data_1d[:, 10]
taylor_val_1d_8 = data_1d[:, 11]
cheb_err_1d_8 = data_1d[:, 12]
taylor_err_1d_8 = data_1d[:, 13]

# Create figure with 2 rows (approximations and errors) and 3 columns (orders)
fig, axes = plt.subplots(2, 3, figsize=(15, 8))
# No overall title per output request

orders = [4, 6, 8]
cheb_val_1d = [cheb_val_1d_4, cheb_val_1d_6, cheb_val_1d_8]
taylor_val_1d = [taylor_val_1d_4, taylor_val_1d_6, taylor_val_1d_8]
cheb_err_1d = [cheb_err_1d_4, cheb_err_1d_6, cheb_err_1d_8]
taylor_err_1d = [taylor_err_1d_4, taylor_err_1d_6, taylor_err_1d_8]

for idx, order in enumerate(orders):
    # Top row: Approximated function values
    ax_func = axes[0, idx]
    ax_func.plot(x_1d, true_1d, 'k-', linewidth=2.5, label='True function', alpha=0.8, zorder=1)
    ax_func.plot(x_1d, cheb_val_1d[idx], 'b-', linewidth=2, label='Chebyshev', alpha=0.7, zorder=3)
    ax_func.plot(x_1d, taylor_val_1d[idx], 'r-', linewidth=2, label='Taylor', alpha=0.7, zorder=2)

    ax_func.set_xlabel('x', fontsize=20)
    ax_func.set_ylabel('exp(-x^2)', fontsize=20)
    ax_func.set_title(f'Order {order}', fontsize=18)
    ax_func.grid(True, alpha=0.3)
    ax_func.legend(fontsize=12, loc='best')

    # Bottom row: Errors
    ax_err = axes[1, idx]
    ax_err.plot(x_1d, cheb_err_1d[idx], 'b-', linewidth=2, label='Chebyshev')
    ax_err.plot(x_1d, taylor_err_1d[idx], 'r-', linewidth=2, label='Taylor')

    ax_err.set_xlabel('x', fontsize=20)
    ax_err.set_ylabel('Absolute Error', fontsize=20)
    ax_err.set_yscale('log')
    ax_err.grid(True, alpha=0.3, which='both')
    ax_err.legend(fontsize=12)

plt.tight_layout()
plt.savefig('data_out/gaussian_1d_comparison.png', dpi=150, bbox_inches='tight')
print('\nPlot saved as data_out/gaussian_1d_comparison.png')
plt.show()

