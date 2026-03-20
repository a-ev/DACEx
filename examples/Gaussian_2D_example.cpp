#include "dacex/DACEx.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <filesystem>

using namespace DACE;

// 2D Gaussian function
double gaussian2D(double x, double y) {
    return std::exp(-x * x - y * y);
}

int main() {
    std::cout << "2D Gaussian Function Approximation Test\n";
    std::cout << "========================================\n";
    std::cout << "Comparing 4th, 6th, and 8th order Chebyshev and Taylor polynomials\n";
    std::cout << "Approximating 2D Gaussian: exp(-x^2 - y^2) on domain [-1.5, 1.5] x [-1.5, 1.5]\n\n";

    // Domain parameters
    double x_min = -1.5;
    double x_max = 1.5;
    double center = (x_min + x_max) / 2.0;
    double half_width = (x_max - x_min) / 2.0;

    // Create 2D evaluation grid (same for all orders)
    // Use odd number of points so that x=0 and y=0 are naturally included
    int grid_points = 101;  // 101x101 grid (odd number centers at 0)
    std::vector<std::vector<double>> eval_points;
    std::vector<double> true_vals;

    for (int i = 0; i < grid_points; ++i) {
        for (int j = 0; j < grid_points; ++j) {
            double x = x_min + (x_max - x_min) * i / (grid_points - 1);
            double y = x_min + (x_max - x_min) * j / (grid_points - 1);
            std::vector<double> point = {x, y};
            eval_points.push_back(point);
            true_vals.push_back(gaussian2D(x, y));
        }
    }

    // Test orders 4, 6, 8
    std::vector<int> orders = {4, 6, 8};
    std::cout << std::fixed << std::setprecision(6);

    // Storage for all results
    std::vector<std::vector<double>> all_cheb_errors(orders.size());
    std::vector<std::vector<double>> all_taylor_errors(orders.size());
    std::vector<std::vector<double>> all_cheb_vals(orders.size());
    std::vector<std::vector<double>> all_taylor_vals(orders.size());

    for (int order_idx = 0; order_idx < orders.size(); ++order_idx) {
        int order = orders[order_idx];

        // Initialize for 2 variables
        DA::init(order, 2);
        DACE::DA::setEps(1e-40);
        
        // Create Chebyshev variables
        Polynomial<Chebyshev> x_cheb = center + half_width*DA(1);
        Polynomial<Chebyshev> y_cheb = center + half_width*DA(2);
        Polynomial<Chebyshev> cheb_poly = exp(-x_cheb * x_cheb - y_cheb * y_cheb);

        // Create Taylor variables
        Polynomial<Taylor> x_taylor = center + DA(1);
        Polynomial<Taylor> y_taylor = center + DA(2);
        Polynomial<Taylor> taylor_poly = exp(-x_taylor * x_taylor - y_taylor * y_taylor);

        // Evaluate at all grid points
        std::vector<double> cheb_errors;
        std::vector<double> taylor_errors;
        std::vector<double> cheb_vals;
        std::vector<double> taylor_vals;
        double cheb_max_err = 0.0, cheb_mean_err = 0.0;
        double taylor_max_err = 0.0, taylor_mean_err = 0.0;

        // Chebyshev evaluation
        for (int i = 0; i < eval_points.size(); ++i) {
            std::vector<double> point = eval_points[i];
            double x = point[0];
            double y = point[1];

            // Chebyshev approximation
            std::vector<double> point_cheb(2);
            point_cheb[0] = x / half_width;
            point_cheb[1] = y / half_width;
            double cheb_val = cheb_poly.eval(point_cheb);
            double cheb_err = std::abs(cheb_val - true_vals[i]);
            cheb_errors.push_back(cheb_err);
            cheb_vals.push_back(cheb_val);
            cheb_max_err = std::max(cheb_max_err, cheb_err);
        }

        // Taylor evaluation
        for (int i = 0; i < eval_points.size(); ++i) {
            std::vector<double> point = eval_points[i];
            double x = point[0];
            double y = point[1];

            // Taylor approximation
            std::vector<double> point_taylor(2);
            point_taylor[0] = x;
            point_taylor[1] = y;
            double taylor_val = taylor_poly.eval(point_taylor);
            double taylor_err = std::abs(taylor_val - true_vals[i]);
            taylor_errors.push_back(taylor_err);
            taylor_vals.push_back(taylor_val);
            taylor_max_err = std::max(taylor_max_err, taylor_err);
        }

        // Compute means
        cheb_mean_err = std::accumulate(cheb_errors.begin(), cheb_errors.end(), 0.0) / cheb_errors.size();
        taylor_mean_err = std::accumulate(taylor_errors.begin(), taylor_errors.end(), 0.0) / taylor_errors.size();

        // Compute standard deviations
        double cheb_variance = 0.0, taylor_variance = 0.0;
        for (double err : cheb_errors) {
            cheb_variance += (err - cheb_mean_err) * (err - cheb_mean_err);
        }
        for (double err : taylor_errors) {
            taylor_variance += (err - taylor_mean_err) * (err - taylor_mean_err);
        }
        cheb_variance /= cheb_errors.size();
        taylor_variance /= taylor_errors.size();
        double cheb_std = std::sqrt(cheb_variance);
        double taylor_std = std::sqrt(taylor_variance);

        // Print statistics
        std::cout << "=== Order " << order << " ===\n";
        std::cout << std::setprecision(6);
        std::cout << "Chebyshev - Mean: " << cheb_mean_err << ", 1-sigma: " << cheb_std << ", 2-sigma: " << 2*cheb_std << ", 3-sigma: " << 3*cheb_std << ", Max: " << cheb_max_err << "\n";
        std::cout << "Taylor    - Mean: " << taylor_mean_err << ", 1-sigma: " << taylor_std << ", 2-sigma: " << 2*taylor_std << ", 3-sigma: " << 3*taylor_std << ", Max: " << taylor_max_err << "\n";

        double improvement;
        if (cheb_max_err < taylor_max_err) {
            improvement = (taylor_max_err - cheb_max_err) / taylor_max_err * 100.0;
            std::cout << "Chebyshev is " << improvement << "% better (max error)\n\n";
        } else {
            improvement = (cheb_max_err - taylor_max_err) / cheb_max_err * 100.0;
            std::cout << "Chebyshev is " << improvement << "% better (max error)\n\n";
        }

        // Store errors
        all_cheb_errors[order_idx] = cheb_errors;
        all_taylor_errors[order_idx] = taylor_errors;
        all_cheb_vals[order_idx] = cheb_vals;
        all_taylor_vals[order_idx] = taylor_vals;
    }

    // Write comparison data to file
    // Find project root by looking for CMakeLists.txt
    std::filesystem::path exe_path = std::filesystem::current_path();
    std::filesystem::path project_root = exe_path;
    while (project_root.has_parent_path() && !std::filesystem::exists(project_root / "CMakeLists.txt")) {
        project_root = project_root.parent_path();
    }
    std::filesystem::path data_out_dir = project_root / "data_out";
    std::filesystem::create_directories(data_out_dir);
    std::string output_filename = (data_out_dir / "comparison_data_2d.txt").string();
    std::ofstream outfile(output_filename);
    outfile << "x y true_gaussian ";
    for (int i = 0; i < orders.size(); ++i) {
        outfile << "cheb_val_" << orders[i] << " taylor_val_" << orders[i] << " "
                << "cheb_error_" << orders[i] << " taylor_error_" << orders[i] << " ";
    }
    outfile << "\n";

    for (int i = 0; i < eval_points.size(); ++i) {
        outfile << std::scientific << std::setprecision(6);
        outfile << eval_points[i][0] << " " << eval_points[i][1] << " " << true_vals[i] << " ";
        for (int order_idx = 0; order_idx < orders.size(); ++order_idx) {
            outfile << all_cheb_vals[order_idx][i] << " "
                    << all_taylor_vals[order_idx][i] << " "
                    << all_cheb_errors[order_idx][i] << " "
                    << all_taylor_errors[order_idx][i];
            if (order_idx < orders.size() - 1) outfile << " ";
        }
        outfile << "\n";
    }
    outfile.close();

    std::cout << "Data written to comparison_data_2d.txt\n";
    std::cout << "Use Python/matplotlib to create plots.\n";

    return 0;
}

