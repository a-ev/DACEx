#include "dacex/DACEx.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <filesystem>

using namespace DACE;

// Gaussian function
double gaussian(double x) {
    return std::exp(-x * x);
}

int main() {
    std::cout << "Gaussian Function Approximation Test\n";
    std::cout << "====================================\n";
    std::cout << "Comparing 4th, 6th, and 8th order Chebyshev and Taylor polynomials\n";
    std::cout << "Approximating Gaussian: exp(-x^2) on domain [-1.5, 1.5]\n\n";

    // Domain parameters - centered at origin
    double x_min = -1.5;
    double x_max = 1.5;
    double center = (x_min + x_max) / 2.0;
    double half_width = (x_max - x_min) / 2.0;

    // Create evaluation grid (same for all orders)
    int num_points = 100;
    std::vector<double> eval_points;
    std::vector<double> true_vals;

    for (int i = 0; i < num_points; ++i) {
        double x = x_min + (x_max - x_min) * i / (num_points - 1);
        eval_points.push_back(x);
        true_vals.push_back(gaussian(x));
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
        std::cout << "\n=== Order " << order << " ===\n";
        DA::init(order, 1);
        DACE::DA::setEps(1e-40);

        // Time Chebyshev variable creation
        auto cheb_var_start = std::chrono::high_resolution_clock::now();
        Polynomial<Chebyshev> x_cheb = center + half_width*DA(1);
        auto cheb_var_end = std::chrono::high_resolution_clock::now();
        
        // Time Chebyshev Gaussian computation
        auto cheb_gauss_start = std::chrono::high_resolution_clock::now();
        Polynomial<Chebyshev> cheb_approx = exp(-x_cheb * x_cheb);
        auto cheb_gauss_end = std::chrono::high_resolution_clock::now();

        // Time Taylor variable creation
        auto taylor_var_start = std::chrono::high_resolution_clock::now();
        Polynomial<Taylor> x_taylor = center + DA(1);
        auto taylor_var_end = std::chrono::high_resolution_clock::now();
        
        // Time Taylor Gaussian computation
        auto taylor_gauss_start = std::chrono::high_resolution_clock::now();
        Polynomial<Taylor> taylor_approx = exp(-x_taylor * x_taylor);
        auto taylor_gauss_end = std::chrono::high_resolution_clock::now();

        // Evaluate and compute errors
        std::vector<double> cheb_err_vals, taylor_err_vals;
        std::vector<double> cheb_vals, taylor_vals;
        
        // Time Chebyshev evaluation
        auto cheb_eval_start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num_points; ++i) {
            double x = eval_points[i];
            double true_val = true_vals[i];

            // Chebyshev approximation
            std::vector<double> point_cheb(1);
            point_cheb[0] = x / half_width;
            double cheb_val = cheb_approx.eval(point_cheb);
            double cheb_err = std::abs(cheb_val - true_val);
            cheb_err_vals.push_back(cheb_err);
            cheb_vals.push_back(cheb_val);
        }
        auto cheb_eval_end = std::chrono::high_resolution_clock::now();

        // Time Taylor evaluation
        auto taylor_eval_start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num_points; ++i) {
            double x = eval_points[i];
            double true_val = true_vals[i];

            // Taylor approximation
            std::vector<double> point_taylor(1);
            point_taylor[0] = x;
            double taylor_val = taylor_approx.eval(point_taylor);
            double taylor_err = std::abs(taylor_val - true_val);
            taylor_err_vals.push_back(taylor_err);
            taylor_vals.push_back(taylor_val);
        }
        auto taylor_eval_end = std::chrono::high_resolution_clock::now();

        // Statistics
        double cheb_max = *std::max_element(cheb_err_vals.begin(), cheb_err_vals.end());
        double cheb_mean = std::accumulate(cheb_err_vals.begin(), cheb_err_vals.end(), 0.0) / cheb_err_vals.size();
        double taylor_max = *std::max_element(taylor_err_vals.begin(), taylor_err_vals.end());
        double taylor_mean = std::accumulate(taylor_err_vals.begin(), taylor_err_vals.end(), 0.0) / taylor_err_vals.size();

        std::cout << "Chebyshev: Max Error = " << cheb_max << ", Mean Error = " << cheb_mean << "\n";
        std::cout << "Taylor:    Max Error = " << taylor_max << ", Mean Error = " << taylor_mean << "\n";
        double improvement = (taylor_max - cheb_max) / taylor_max * 100.0;
        std::cout << "Chebyshev is " << improvement << "% better (max error)\n";

        // Compute timing results
        auto cheb_var_time = std::chrono::duration_cast<std::chrono::microseconds>(cheb_var_end - cheb_var_start).count();
        auto cheb_gauss_time = std::chrono::duration_cast<std::chrono::microseconds>(cheb_gauss_end - cheb_gauss_start).count();
        auto taylor_var_time = std::chrono::duration_cast<std::chrono::microseconds>(taylor_var_end - taylor_var_start).count();
        auto taylor_gauss_time = std::chrono::duration_cast<std::chrono::microseconds>(taylor_gauss_end - taylor_gauss_start).count();
        auto cheb_eval_time = std::chrono::duration_cast<std::chrono::microseconds>(cheb_eval_end - cheb_eval_start).count();
        auto taylor_eval_time = std::chrono::duration_cast<std::chrono::microseconds>(taylor_eval_end - taylor_eval_start).count();

        std::cout << "\nTiming (microseconds):\n";
        std::cout << "  Chebyshev - Variable: " << cheb_var_time << " µs, Gaussian: " << cheb_gauss_time << " µs, Evaluation: " << cheb_eval_time << " µs\n";
        std::cout << "  Taylor    - Variable: " << taylor_var_time << " µs, Gaussian: " << taylor_gauss_time << " µs, Evaluation: " << taylor_eval_time << " µs\n";
        std::cout << "  Ratio (Cheb/Taylor) - Variable: " << std::fixed << std::setprecision(2) 
                  << (double)cheb_var_time/taylor_var_time << "x, Gaussian: " 
                  << (double)cheb_gauss_time/taylor_gauss_time << "x, Evaluation: " 
                  << (double)cheb_eval_time/taylor_eval_time << "x\n";

        // Store errors and values for plotting
        all_cheb_errors[order_idx] = cheb_err_vals;
        all_taylor_errors[order_idx] = taylor_err_vals;
        all_cheb_vals[order_idx] = cheb_vals;
        all_taylor_vals[order_idx] = taylor_vals;
    }

    // Write data to file for plotting
    // Find project root by looking for CMakeLists.txt
    std::filesystem::path exe_path = std::filesystem::current_path();
    std::filesystem::path project_root = exe_path;
    while (project_root.has_parent_path() && !std::filesystem::exists(project_root / "CMakeLists.txt")) {
        project_root = project_root.parent_path();
    }
    std::filesystem::path data_out_dir = project_root / "data_out";
    std::filesystem::create_directories(data_out_dir);
    std::string output_filename = (data_out_dir / "comparison_data.txt").string();
    std::ofstream outfile(output_filename);
    outfile << "x true_gaussian ";
    for (int i = 0; i < orders.size(); ++i) {
        outfile << "cheb_val_" << orders[i] << " taylor_val_" << orders[i] << " "
                << "cheb_error_" << orders[i] << " taylor_error_" << orders[i] << " ";
    }
    outfile << "\n";

    for (int i = 0; i < num_points; ++i) {
        outfile << std::fixed << std::setprecision(8)
                << eval_points[i] << " "
                << true_vals[i] << " ";
        for (int j = 0; j < orders.size(); ++j) {
            outfile << all_cheb_vals[j][i] << " "
                    << all_taylor_vals[j][i] << " "
                    << all_cheb_errors[j][i] << " "
                    << all_taylor_errors[j][i] << " ";
        }
        outfile << "\n";
    }
    outfile.close();

    std::cout << "\nData written to comparison_data.txt\n";
    std::cout << "Use Python/matplotlib to create plots.\n";

    return 0;
}

