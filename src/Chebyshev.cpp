// Copyright (c) 2026 Adam Evans (@a-ev)

#include "dacex/Chebyshev.h"
#include <unordered_map>
#include <functional>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct ChebyshevSplitCache {
    struct SplitResult {
        double sum_coeff;
        unsigned int sum_exp;
        double diff_coeff;
        unsigned int diff_exp;
    };
    
    std::map<std::pair<unsigned int, unsigned int>, SplitResult> cache;
    
    SplitResult computeAndCache(unsigned int n, unsigned int m) {
        auto key = std::make_pair(n, m);
        auto it = cache.find(key);
        if (it != cache.end()) {
            return it->second;
        }
        
        SplitResult result;
        result.sum_coeff = 0.5;
        result.sum_exp = n + m;
        result.diff_coeff = 0.5;
        result.diff_exp = (n > m) ? (n - m) : (m - n);
        
        cache[key] = result;
        return result;
    }
};

static ChebyshevSplitCache g_split_cache;



Polynomial<Chebyshev>::Polynomial(const Polynomial<Chebyshev> &p) : PolynomialBase() {
    m_da = p.m_da;
}

Polynomial<Chebyshev>::Polynomial(Polynomial<Chebyshev> &&p) noexcept : PolynomialBase() {
    m_da = std::move(p.m_da);
}



Polynomial<Chebyshev>& Polynomial<Chebyshev>::operator=(Polynomial<Chebyshev> &&p) noexcept {
    if (this != &p) {
        m_da = p.m_da;
    }
    return *this;
}

Polynomial<Chebyshev>& Polynomial<Chebyshev>::operator=(const Polynomial<Chebyshev> &p) {
    if (this != &p) {
        m_da = std::move(p.m_da);
    }
    return *this;
}

Polynomial<Chebyshev>& Polynomial<Chebyshev>::operator=(const double c) {
    m_da = c;
    return *this;
}



Polynomial<Chebyshev>& Polynomial<Chebyshev>::operator+=(const Polynomial<Chebyshev> &p) {
    *this = *this + p;
    return *this;
}

Polynomial<Chebyshev>& Polynomial<Chebyshev>::operator+=(const double c) {
    *this = *this + c;
    return *this;
}

Polynomial<Chebyshev>& Polynomial<Chebyshev>::operator-=(const Polynomial<Chebyshev> &p) {
    *this = *this - p;
    return *this;
}

Polynomial<Chebyshev>& Polynomial<Chebyshev>::operator-=(const double c) {
    *this = *this - c;
    return *this;
}

Polynomial<Chebyshev>& Polynomial<Chebyshev>::operator*=(const Polynomial<Chebyshev> &p) {
    *this = *this * p;  // Delegate to binary operator (Chebyshev multiplication algebra)
    return *this;
}

Polynomial<Chebyshev>& Polynomial<Chebyshev>::operator*=(const double c) {
    *this = *this * c;  // Delegate to binary operator (scalar multiplication)
    return *this;
}

Polynomial<Chebyshev>& Polynomial<Chebyshev>::operator/=(const Polynomial<Chebyshev> &p) {
    *this = *this / p;  // Delegate to binary operator (Chebyshev division algebra)
    return *this;
}

Polynomial<Chebyshev>& Polynomial<Chebyshev>::operator/=(const double c) {
    *this = *this / c;  // Delegate to binary operator (scalar division)
    return *this;
}

// ============================================================================
// Helper functions
// ============================================================================

// --- Helper: Term Accumulator for Chebyshev Indices ---
void addTermToDA(DACE::DA &da, const std::vector<unsigned int>& indices, double coeff) {
    // Use DACE's epsilon threshold instead of hardcoded value
    if (std::abs(coeff) < DACE::DA::getEps()) return;
    
    // In Chebyshev algebra, indices represent T_n1(x1) * T_n2(x2) * ...
    // We need to ADD this coefficient to the existing coefficient at these indices
    
    double existing_coeff = da.getCoefficient(indices);
    double new_coeff = existing_coeff + coeff;
    da.setCoefficient(indices, new_coeff);
}

static Polynomial<Chebyshev> evaluateOnInterval(const Polynomial<Chebyshev>& p,
                                                const std::function<double(double)>& func,
                                                double min_p,
                                                double max_p);

// --- Helper: Inverse (1/P) using Newton-Raphson ---
Polynomial<Chebyshev> inv(const Polynomial<Chebyshev> &p) {
    // Always use algebraic Newton-Raphson for inversion
    // Check for singularity
    double c0 = p.cons();
    if (std::abs(c0) < DACE::DA::getEps()) {
        throw std::runtime_error("Chebyshev Division: Divisor must have a non-zero constant part.");
    }

    // Normalize P: P_norm = p / c0 (so constant part = 1)
    Polynomial<Chebyshev> P_norm = p * (1.0 / c0);

    // Initial Guess: X = 1
    Polynomial<Chebyshev> X = 1.0; 

    // Newton-Raphson iterations: X = X * (2 - P_norm * X)
    // General formula: iterations = ceil(log2(order+1)) + 2
    // Quadratic convergence requires log2(order) iterations for order-n polynomial precision
    // +2 is a safety buffer to guarantee convergence in all cases
    unsigned int max_order = DACE::DA::getMaxOrder();
    int iterations = static_cast<int>(std::ceil(std::log2(max_order + 1))) + 2;

    for (int i = 0; i < iterations; ++i) {
        Polynomial<Chebyshev> prod = P_norm * X;
        Polynomial<Chebyshev> correction = 2.0 - prod;
        X = X * correction;
    }

    // Denormalize: return 1/c0 * X
    return X * (1.0 / c0);
}

// ============================================================================
// Basic Arithmetic Operations (Skeleton implementations)
// ============================================================================

Polynomial<Chebyshev> Polynomial<Chebyshev>::operator-() const {
    // Negation is basis-independent: just negate all coefficients
    return Polynomial<Chebyshev>(-m_da);
}

Polynomial<Chebyshev> operator+(const Polynomial<Chebyshev> &p1, const Polynomial<Chebyshev> &p2) {
    // Addition is basis-independent: coefficient-wise addition works for any basis
    // (p1 coefficients) + (p2 coefficients) = (sum of coefficients)
    return Polynomial<Chebyshev>(p1.m_da + p2.m_da);
}

Polynomial<Chebyshev> operator+(const Polynomial<Chebyshev> &p, const double c) {
    return Polynomial<Chebyshev>(p.m_da + c);
}

Polynomial<Chebyshev> operator+(const double c, const Polynomial<Chebyshev> &p) {
    return Polynomial<Chebyshev>(c + p.m_da);
}

Polynomial<Chebyshev> operator-(const Polynomial<Chebyshev> &p1, const Polynomial<Chebyshev> &p2) {
    return Polynomial<Chebyshev>(p1.m_da - p2.m_da);
}

Polynomial<Chebyshev> operator-(const Polynomial<Chebyshev> &p, const double c) {
    return Polynomial<Chebyshev>(p.m_da - c);
}

Polynomial<Chebyshev> operator-(const double c, const Polynomial<Chebyshev> &p) {
    return Polynomial<Chebyshev>(c - p.m_da);
}

Polynomial<Chebyshev> operator*(const Polynomial<Chebyshev> &p, const double c) {
    return Polynomial<Chebyshev>(p.m_da * c);
}

Polynomial<Chebyshev> operator*(const double c, const Polynomial<Chebyshev> &p) {
    return Polynomial<Chebyshev>(c * p.m_da);
}

Polynomial<Chebyshev> operator*(const Polynomial<Chebyshev> &p1, const Polynomial<Chebyshev> &p2) {
    const DACE::DA &da1 = p1.m_da;
    const DACE::DA &da2 = p2.m_da;

    unsigned int nvar = DACE::DA::getMaxVariables();
    unsigned int hard_max_order = DACE::DA::getMaxOrder();

    auto terms1 = da1.getMonomials();
    auto terms2 = da2.getMonomials();

    DACE::DA result_da = 0.0;
    std::vector<unsigned int> new_exponents(nvar);

    for (const auto& t1 : terms1) {
        for (const auto& t2 : terms2) {
            // Identify dimensions where both terms have nonzero indices
            std::vector<unsigned int> active_dims;
            for (unsigned int d = 0; d < nvar; ++d) {
                if (t1.m_jj[d] != 0 && t2.m_jj[d] != 0) {
                    active_dims.push_back(d);
                }
            }
            
            int num_active = static_cast<int>(active_dims.size());
            int num_branches = 1 << num_active;  // 2^(num_active)
            double base_coeff = t1.m_coeff * t2.m_coeff;
            
            for (int k = 0; k < num_branches; ++k) {
                double current_coeff = base_coeff;
                
                // First: set exponents for passive dimensions (at least one is zero)
                for (unsigned int d = 0; d < nvar; ++d) {
                    if (t1.m_jj[d] == 0 || t2.m_jj[d] == 0) {
                        new_exponents[d] = t1.m_jj[d] + t2.m_jj[d];
                    }
                }
                
                // Second: set exponents for active dimensions (both nonzero) based on branch k
                bool valid_term = true;
                for (int active_idx = 0; active_idx < num_active; ++active_idx) {
                    unsigned int d = active_dims[active_idx];
                    unsigned int n = t1.m_jj[d];
                    unsigned int m = t2.m_jj[d];
                    
                    current_coeff *= 0.5;  // Chebyshev product rule: T_n * T_m = 0.5*(T_{n+m} + T_{|n-m|})

                    if ((k >> active_idx) & 1) {
                        new_exponents[d] = (n > m) ? (n - m) : (m - n);
                    } else {
                        new_exponents[d] = n + m;
                    }
                }
                
                // Third: compute total order and check truncation
                unsigned int current_total_order = 0;
                for (unsigned int d = 0; d < nvar; ++d) {
                    current_total_order += new_exponents[d];
                }
                
                if (current_total_order <= hard_max_order) {
                    addTermToDA(result_da, new_exponents, current_coeff);
                }
            }
        }
    }

    return Polynomial<Chebyshev>(result_da);
}

Polynomial<Chebyshev> operator/(const Polynomial<Chebyshev> &p, const double c) {
    // Scalar division is basis-independent: divide all coefficients by scalar
    return Polynomial<Chebyshev>(p.m_da / c);
}

Polynomial<Chebyshev> operator/(const double c, const Polynomial<Chebyshev> &p) {
    return c * inv(p);
}

Polynomial<Chebyshev> operator/(const Polynomial<Chebyshev> &p1, const Polynomial<Chebyshev> &p2) {
    return p1 * inv(p2);
}

Polynomial<Chebyshev> Polynomial<Chebyshev>::deriv(const unsigned int var_idx) const {
    if (var_idx < 1 || var_idx > DACE::DA::getMaxVariables()) {
        throw std::runtime_error("Derivative: Variable index out of bounds");
    }

    unsigned int max_order = DACE::DA::getMaxOrder();
    unsigned int max_vars = DACE::DA::getMaxVariables();
    
    std::vector<Polynomial<Chebyshev>> T(max_order + 1);
    Polynomial<Chebyshev> x = DACE::DA(var_idx);
    
    T[0] = 1.0;
    if (max_order >= 1) {
        T[1] = x;
        for (unsigned int k = 2; k <= max_order; ++k) {
            T[k] = (T[k-1] * x * 2.0) - T[k-2];
        }
    }
    
    Polynomial<Chebyshev> result = 0.0;
    std::vector<DACE::Monomial> monomials = m_da.getMonomials();
    
    for (const auto& mono : monomials) {
        double coeff = mono.m_coeff;
        if (std::abs(coeff) < DACE::DA::getEps()) continue;
        
        std::vector<unsigned int> indices = mono.m_jj;
        unsigned int n = indices[var_idx - 1];
        
        if (n == 0) continue;
        
        for (int k = n - 1; k >= 0; k -= 2) {
            double deriv_coeff = coeff * (2.0 * n);
            
            // Halve coefficient for T_0 term
            if (k == 0) deriv_coeff *= 0.5;
            
            // Build the term: deriv_coeff * T_k(x_var) * Product of T_j(other vars)
            Polynomial<Chebyshev> term = T[k];  // T_k for the differentiated variable
            
            // Multiply by Chebyshev terms for all OTHER variables
            for (unsigned int var = 0; var < max_vars; ++var) {
                if (var == var_idx - 1) continue;  // Skip the differentiated variable
                
                unsigned int order = indices[var];
                if (order == 0) continue;  // T_0 = 1, no multiplication needed
                
                // Need to build T_order for this variable
                // Build using recurrence with that variable
                Polynomial<Chebyshev> var_x = DACE::DA(var + 1);
                Polynomial<Chebyshev> T_var_0 = 1.0;
                Polynomial<Chebyshev> T_var_1 = var_x;
                
                if (order == 0) {
                    term = term * T_var_0;
                } else if (order == 1) {
                    term = term * T_var_1;
                } else {
                    Polynomial<Chebyshev> T_var_prev = T_var_0;
                    Polynomial<Chebyshev> T_var_curr = T_var_1;
                    for (unsigned int i = 2; i <= order; ++i) {
                        Polynomial<Chebyshev> T_var_next = (T_var_curr * var_x * 2.0) - T_var_prev;
                        T_var_prev = T_var_curr;
                        T_var_curr = T_var_next;
                    }
                    term = term * T_var_curr;
                }
            }
            
            // Add this derivative term to result
            result = result + (term * deriv_coeff);
        }
    }

    return result;
}

// Mixed partial derivatives - sequential application of single variable derivatives
Polynomial<Chebyshev> Polynomial<Chebyshev>::deriv(const std::vector<unsigned int> ind) const {
    // Start with the current polynomial
    Polynomial<Chebyshev> result = *this;
    
    // Apply derivatives in sequence for each variable
    for (unsigned int var = 0; var < ind.size(); ++var) {
        unsigned int deriv_order = ind[var];
        
        // Apply derivative 'deriv_order' times with respect to variable (var+1)
        for (unsigned int i = 0; i < deriv_order; ++i) {
            result = result.deriv(var + 1);  // DACE variables are 1-indexed
        }
    }
    
    return result;
}

Polynomial<Chebyshev> Polynomial<Chebyshev>::integ(const unsigned int var_idx) const {
    // 1. Safety Check
    if (var_idx < 1 || var_idx > DACE::DA::getMaxVariables()) {
        throw std::runtime_error("Integral: Variable index out of bounds");
    }

    const DACE::DA& da_ref = m_da;
    unsigned int max_order = DACE::DA::getMaxOrder();

    // ---------------------------------------------------------
    // Step 1: Pre-calculate Chebyshev Basis T_k(x)
    // We need up to T_{max+1} for the integral formula.
    // ---------------------------------------------------------
    std::vector<Polynomial<Chebyshev>> T(max_order + 2);
    
    // Define the integration variable x
    Polynomial<Chebyshev> x = DACE::DA(var_idx);

    // Recurrence: T_0 = 1, T_1 = x, T_k = 2xT_{k-1} - T_{k-2}
    T[0] = 1.0;
    T[1] = x;
    
    for (unsigned int k = 2; k <= max_order + 1; ++k) {
        T[k] = (T[k-1] * x * 2.0) - T[k-2];
    }

    // ---------------------------------------------------------
    // Step 2: Accumulate the Integral
    // ---------------------------------------------------------
    Polynomial<Chebyshev> result = 0.0;

    for (unsigned int n = 0; n <= max_order; ++n) {
        // Get coefficient for T_n(x)
        std::vector<unsigned int> indices(DACE::DA::getMaxVariables(), 0);
        indices[var_idx - 1] = n;
        
        double c_n = da_ref.getCoefficient(indices);

        // Skip zeros
        if (std::abs(c_n) < DACE::DA::getEps()) continue;

        if (n == 0) {
            // Int(T_0) = T_1
            result = result + (T[1] * c_n);
        }
        else if (n == 1) {
            // Int(T_1) = 0.25*T_2 + 0.25*T_0
            result = result + (T[2] * (0.25 * c_n));
            result = result + (T[0] * (0.25 * c_n));
        }
        else {
            // Int(T_n) = 0.5 * [ T_{n+1}/(n+1) - T_{n-1}/(n-1) ]
            double term1_scale = 0.5 * c_n / (double)(n + 1);
            double term2_scale = 0.5 * c_n / (double)(n - 1);

            result = result + (T[n + 1] * term1_scale);
            result = result - (T[n - 1] * term2_scale);
        }
    }

    return result;
}

// Mixed integration - sequential application of single variable integration  
Polynomial<Chebyshev> Polynomial<Chebyshev>::integ(const std::vector<unsigned int> ind) const {
    // Start with the current polynomial
    Polynomial<Chebyshev> result = *this;
    
    // Apply integration in sequence for each variable
    for (unsigned int var = 0; var < ind.size(); ++var) {
        unsigned int integ_order = ind[var];
        
        // Apply integration 'integ_order' times with respect to variable (var+1)
        for (unsigned int i = 0; i < integ_order; ++i) {
            result = result.integ(var + 1);  // DACE variables are 1-indexed
        }
    }
    
    return result;
}

// ============================================================================
// I/O routines
// ============================================================================

// toString() function is now provided by PolynomialBase

std::string Polynomial<Chebyshev>::toString() const {
    // Custom toString that replaces "EXPONENTS" with "INDICES"
    std::stringstream buffer;

    buffer << m_da; // Use the built-in DACE print
    std::string output = buffer.str();
    
    size_t pos = output.find("EXPONENTS");
    if (pos != std::string::npos) {
        output.replace(pos, 9, "INDICES  "); // Pad with spaces to maintain alignment
    }

    return output;
}

// write() function is now provided by PolynomialBase

std::ostream& operator<<(std::ostream& os, const Polynomial<Chebyshev>& poly) {
    // Uses normal <<, but finds & replaces "EXPONENTS" with "INDICES"
    std::stringstream buffer;

    buffer << poly.m_da; // Use the built-in DACE print (which produces the table)
    std::string output = buffer.str(); // Extract the string
    
    size_t pos = output.find("EXPONENTS");
    if (pos != std::string::npos) {
        output.replace(pos, 9, "INDICES  "); // Pad with spaces to maintain alignment
    }

    os << output;

    return os;
}

std::istream& operator>> (std::istream &in, Polynomial<Chebyshev> &p) {
    in >> p.m_da;
    return in;
}



// ============================================================================
// Test routines
// ============================================================================
static double chebVal(unsigned int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    double t_prev = 1.0, t_curr = x, two_x = 2.0 * x;
    for (unsigned int k = 2; k <= n; ++k) {
        double t_next = two_x * t_curr - t_prev;
        t_prev = t_curr; t_curr = t_next;
    }
    return t_curr;
}

static Polynomial<Chebyshev> chebPoly(unsigned int n, const Polynomial<Chebyshev>& P) {
    if (n == 0) return 1.0;
    if (n == 1) return P;
    Polynomial<Chebyshev> t_prev = 1.0, t_curr = P, two_P = P * 2.0;
    for (unsigned int k = 2; k <= n; ++k) {
        Polynomial<Chebyshev> t_next = (two_P * t_curr) - t_prev;
        t_prev = t_curr; t_curr = t_next;
    }
    return t_curr;
}

// --- The Universal Worker (Template) ---
// This sits locally in the .cpp file. It handles the recursion logic once.
template <typename T>
T evalInternalLogic(const DACE::DA* da_ptr, const std::vector<T>& args) {
    if (args.size() != DACE::DA::getMaxVariables()) {
        throw std::runtime_error("Dimension Mismatch in eval");
    }

    // Recursive Lambda
    std::function<T(const DACE::DA*, std::vector<unsigned int>&, int, int)> recurse;
    
    recurse = [&](const DACE::DA* da, std::vector<unsigned int>& powers, int dim, int current_order) -> T {
        T sum(0.0);
        int max_vars = DACE::DA::getMaxVariables();
        int max_order = DACE::DA::getMaxOrder();

        if (dim == max_vars) {
            double coeff = da->getCoefficient(powers);
            if (std::abs(coeff) < DACE::DA::getEps()) return T(0.0);

            T term(coeff);
            for (int v = 0; v < max_vars; ++v) {
                if (powers[v] > 0) {
                    // Use "if constexpr" to choose the right helper at compile time
                    if constexpr (std::is_same_v<T, double>) {
                        term = term * chebVal(powers[v], args[v]);
                    } else {
                        term = term * chebPoly(powers[v], args[v]);
                    }
                }
            }
            return term;
        }

        for (int p = 0; p <= (max_order - current_order); ++p) {
            powers[dim] = p;
            sum = sum + recurse(da, powers, dim + 1, current_order + p);
        }
        return sum;
    };

    std::vector<unsigned int> powers(args.size(), 0);
    return recurse(da_ptr, powers, 0, 0);
}

// --- Class Member Implementations ---

// 1. Symbolic Composition (DACE Vector)
Polynomial<Chebyshev> Polynomial<Chebyshev>::eval(const DACE::AlgebraicVector<Polynomial<Chebyshev>>& args) const {
    
    // Convert DACE vector to std::vector
    std::vector<Polynomial<Chebyshev>> std_args(args.size());
    for(size_t i=0; i<args.size(); ++i) std_args[i] = args[i];

    return evalInternalLogic<Polynomial<Chebyshev>>(&m_da, std_args);
}

// 2. Numeric Evaluation (DACE Vector)
double Polynomial<Chebyshev>::eval(const DACE::AlgebraicVector<double>& args) const {

    std::vector<double> std_args(args.size());
    for(size_t i = 0; i < args.size(); ++i) {
        std_args[i] = args[i];
        // Validate that each value is within [-1, 1] for Chebyshev polynomials
        if (std_args[i] < -1.0 || std_args[i] > 1.0) {
            throw std::runtime_error("Chebyshev evaluation error: All values must be in range [-1, 1]. "
                                   "Value at index " + std::to_string(i) + " is " + std::to_string(std_args[i]));
        }
    }

    return evalInternalLogic<double>(&m_da, std_args);
}

// 3. Numeric Evaluation (Std Vector) -> THIS FIXES YOUR COMPILER ERROR
double Polynomial<Chebyshev>::eval(const std::vector<double>& args) const {
    // Validate that each value is within [-1, 1] for Chebyshev polynomials
    for(size_t i = 0; i < args.size(); ++i) {
        if (args[i] < -1.0 || args[i] > 1.0) {
            throw std::runtime_error("Chebyshev evaluation error: All values must be in range [-1, 1]. "
                                   "Value at index " + std::to_string(i) + " is " + std::to_string(args[i]));
        }
    }
    return evalInternalLogic<double>(&m_da, args);
}

// ============================================================================
// Norm and Estimation Routines
// ============================================================================

double Polynomial<Chebyshev>::abs() const {
    // Maximum absolute value of all coefficients
    return m_da.abs();
}

double Polynomial<Chebyshev>::norm(const unsigned int type) const {
    // Compute different types of norms:
    // type 0: Max norm (default)
    // type 1: Sum norm
    // type >1: Vector norm of given type
    return m_da.norm(type);
}

std::vector<double> Polynomial<Chebyshev>::orderNorm(const unsigned int var, const unsigned int type) const {
    // Extract order-sorted norms
    // var 0: Terms sorted by total order (default)
    // var >0: Terms sorted by exponent of variable var
    return m_da.orderNorm(var, type);
}

std::vector<double> Polynomial<Chebyshev>::estimNorm(const unsigned int var, const unsigned int type, const unsigned int nc) const {
    // Estimate order-sorted norms up to order nc
    return m_da.estimNorm(var, type, nc);
}

std::vector<double> Polynomial<Chebyshev>::estimNorm(std::vector<double> &err, const unsigned int var, const unsigned int type, const unsigned int nc) const {
    // Estimate order-sorted norms with error estimates
    return m_da.estimNorm(err, var, type, nc);
}

DACE::Interval Polynomial<Chebyshev>::bound() const {
    // Compute lower and upper bounds over [-1,1]^n
    return m_da.bound();
}

// ============================================================================
// Chebyshev Series Coefficient Computation
// ============================================================================

namespace ChebyshevSeries {

// Cache for precomputed coefficients
static std::unordered_map<std::string, std::vector<double>> coefficient_cache;

// Compute Chebyshev coefficients using FFT-based method
std::vector<double> computeChebyshevCoeffs(std::function<double(double)> func, unsigned int order) {
    std::vector<double> coeffs(order + 1, 0.0);
    // General formula: N = 8 * (order + 1)
    // Nyquist requires 2×order samples; we use 8× for robust 8× Nyquist oversampling for DCT accuracy
    unsigned int N = 8 * (order + 1);
    
    // Sample the function at Chebyshev points
    std::vector<double> samples(N);
    for (unsigned int k = 0; k < N; ++k) {
        double theta = M_PI * (k + 0.5) / N;
        double x = std::cos(theta);
        samples[k] = func(x);
    }
    
    // Compute Chebyshev coefficients using DCT-like transform
    for (unsigned int n = 0; n <= order; ++n) {
        double sum = 0.0;
        for (unsigned int k = 0; k < N; ++k) {
            double theta = M_PI * (k + 0.5) / N;
            sum += samples[k] * std::cos(n * theta);
        }
        coeffs[n] = (n == 0) ? sum / N : 2.0 * sum / N;
    }
    
    return coeffs;
}

std::vector<double> getSinCoeffs(unsigned int order) {
    std::string key = "sin_" + std::to_string(order);
    if (coefficient_cache.find(key) != coefficient_cache.end()) {
        return coefficient_cache[key];
    }
    
    auto coeffs = computeChebyshevCoeffs([](double x) { return std::sin(x); }, order);
    coefficient_cache[key] = coeffs;
    return coeffs;
}

std::vector<double> getCosCoeffs(unsigned int order) {
    std::string key = "cos_" + std::to_string(order);
    if (coefficient_cache.find(key) != coefficient_cache.end()) {
        return coefficient_cache[key];
    }
    
    auto coeffs = computeChebyshevCoeffs([](double x) { return std::cos(x); }, order);
    coefficient_cache[key] = coeffs;
    return coeffs;
}

std::vector<double> getExpCoeffs(unsigned int order) {
    std::string key = "exp_" + std::to_string(order);
    if (coefficient_cache.find(key) != coefficient_cache.end()) {
        return coefficient_cache[key];
    }
    
    auto coeffs = computeChebyshevCoeffs([](double x) { return std::exp(x); }, order);
    coefficient_cache[key] = coeffs;
    return coeffs;
}

std::vector<double> getLogCoeffs(unsigned int order) {
    std::string key = "log_" + std::to_string(order);
    if (coefficient_cache.find(key) != coefficient_cache.end()) {
        return coefficient_cache[key];
    }
    
    // For log, we work on [0.1, 4.0] mapped to [-1, 1]
    // Transform: x in [-1,1] -> t in [0.1,4.0]: t = 1.95*x + 2.05
    auto coeffs = computeChebyshevCoeffs([](double x) { 
        double t = 1.95 * x + 2.05;  // Map [-1,1] to [0.1,4.0]
        if (t <= 0) return -10.0;  // Avoid log(0)
        return std::log(t); 
    }, order);
    coefficient_cache[key] = coeffs;
    return coeffs;
}

std::vector<double> getSqrtCoeffs(unsigned int order) {
    std::string key = "sqrt_" + std::to_string(order);
    if (coefficient_cache.find(key) != coefficient_cache.end()) {
        return coefficient_cache[key];
    }
    
    // For sqrt, work on [0.1, 4.0] mapped to [-1, 1]  
    // Transform: x in [-1,1] -> t in [0.1,4.0]: t = 1.95*(x+1) + 0.1 = 1.95*x + 2.05
    auto coeffs = computeChebyshevCoeffs([](double x) { 
        double t = 1.95 * x + 2.05;  // Map [-1,1] to [0.1,4.0]
        if (t < 0) return 0.0;
        return std::sqrt(t); 
    }, order);
    coefficient_cache[key] = coeffs;
    return coeffs;
}

std::vector<double> getAsinCoeffs(unsigned int order) {
    std::string key = "asin_" + std::to_string(order);
    if (coefficient_cache.find(key) != coefficient_cache.end()) {
        return coefficient_cache[key];
    }
    
    auto coeffs = computeChebyshevCoeffs([](double x) -> double { 
        if (x < -1.0) return -M_PI/2.0;
        if (x > 1.0) return M_PI/2.0;
        return std::asin(x); 
    }, order);
    coefficient_cache[key] = coeffs;
    return coeffs;
}

std::vector<double> getAtanCoeffs(unsigned int order) {
    std::string key = "atan_" + std::to_string(order);
    if (coefficient_cache.find(key) != coefficient_cache.end()) {
        return coefficient_cache[key];
    }
    
    auto coeffs = computeChebyshevCoeffs([](double x) { return std::atan(x); }, order);
    coefficient_cache[key] = coeffs;
    return coeffs;
}

// Clenshaw algorithm for evaluating Chebyshev series
Polynomial<Chebyshev> evaluateSeries(const std::vector<double>& coeffs, 
                                      const Polynomial<Chebyshev>& p) {
    if (coeffs.empty()) {
        return Polynomial<Chebyshev>(0.0);
    }
    
    if (coeffs.size() == 1) {
        return Polynomial<Chebyshev>(coeffs[0]);
    }
    
    // Clenshaw algorithm: evaluate sum_{n=0}^N c_n T_n(p)
    // Using recurrence: b_{N+1} = b_{N+2} = 0
    //                   b_k = c_k + 2*p*b_{k+1} - b_{k+2}
    // Result: (b_0 - b_2) / 2 + c_0  OR  c_0 + p*b_1 - b_2
    
    int N = coeffs.size() - 1;
    Polynomial<Chebyshev> b_next2(0.0);  // b_{k+2}
    Polynomial<Chebyshev> b_next1(0.0);  // b_{k+1}
    
    // Backward iteration from k = N down to k = 1
    for (int k = N; k >= 1; --k) {
        Polynomial<Chebyshev> b_k = coeffs[k] + (p * b_next1) * 2.0 - b_next2;
        b_next2 = b_next1;
        b_next1 = b_k;
    }
    
    // Final step: result = c_0 + p*b_1 - b_2
    return coeffs[0] + p * b_next1 - b_next2;
}

} // namespace ChebyshevSeries

static Polynomial<Chebyshev> evaluateOnInterval(const Polynomial<Chebyshev>& p,
                                                const std::function<double(double)>& func,
                                                double min_p,
                                                double max_p) {
    double center = 0.5 * (min_p + max_p);
    double radius = 0.5 * (max_p - min_p);
    if (radius <= 1e-12) {
        return Polynomial<Chebyshev>(func(center));
    }
    unsigned int order = DACE::DA::getMaxOrder();
    auto coeffs = ChebyshevSeries::computeChebyshevCoeffs([&](double x) {
        return func(center + radius * x);
    }, order);
    return ChebyshevSeries::evaluateSeries(coeffs, (p - center) / radius);
}

// ============================================================================
// Transcendental Functions Implementation
// ============================================================================

Polynomial<Chebyshev> sin(const Polynomial<Chebyshev>& p) {
    DACE::Interval p_bounds = p.bound();
    return evaluateOnInterval(p, [](double x) { return std::sin(x); },
                              p_bounds.m_lb, p_bounds.m_ub);
}

Polynomial<Chebyshev> cos(const Polynomial<Chebyshev>& p) {
    DACE::Interval p_bounds = p.bound();
    return evaluateOnInterval(p, [](double x) { return std::cos(x); },
                              p_bounds.m_lb, p_bounds.m_ub);
}

Polynomial<Chebyshev> tan(const Polynomial<Chebyshev>& p) {
    return sin(p) / cos(p);
}

Polynomial<Chebyshev> asin(const Polynomial<Chebyshev>& p) {
    DACE::Interval p_bounds = p.bound();
    double min_p = p_bounds.m_lb;
    double max_p = p_bounds.m_ub;
    
    // asin domain is [-1, 1]
    const double eps = 1e-12;
    if (min_p < -1.0) min_p = -1.0;
    if (max_p > 1.0) max_p = 1.0;
    if (max_p < min_p + eps) max_p = min_p + eps;
    
    return evaluateOnInterval(p, [](double x) { return std::asin(x); }, min_p, max_p);
}

Polynomial<Chebyshev> acos(const Polynomial<Chebyshev>& p) {
    DACE::Interval p_bounds = p.bound();
    double min_p = p_bounds.m_lb;
    double max_p = p_bounds.m_ub;
    
    // acos domain is [-1, 1]
    const double eps = 1e-12;
    if (min_p < -1.0) min_p = -1.0;
    if (max_p > 1.0) max_p = 1.0;
    if (max_p < min_p + eps) max_p = min_p + eps;
    
    return evaluateOnInterval(p, [](double x) { return std::acos(x); }, min_p, max_p);
}

Polynomial<Chebyshev> atan(const Polynomial<Chebyshev>& p) {
    DACE::Interval p_bounds = p.bound();
    return evaluateOnInterval(p, [](double x) { return std::atan(x); },
                              p_bounds.m_lb, p_bounds.m_ub);
}

Polynomial<Chebyshev> atan2(const Polynomial<Chebyshev>& y, const Polynomial<Chebyshev>& x) {
    // Four-quadrant arctangent: atan2(y, x) returns angle in [-π, π]
    // Get constant parts to determine quadrant
    const double cx = x.cons();
    const double cy = y.cons();
    
    const double PI = 3.14159265358979323846;
    const double PI_2 = PI / 2.0;
    
    // Handle special case where both constants are zero
    if (cx == 0.0 && cy == 0.0) {
        return Polynomial<Chebyshev>(0.0);
    }
    
    // Choose algorithm based on which constant is larger (for numerical stability)
    if (std::abs(cy) > std::abs(cx)) {
        // Use atan(x/y) approach and adjust by ±π/2
        Polynomial<Chebyshev> ratio = x / y;
        Polynomial<Chebyshev> result = atan(ratio);
        
        if (cy < 0.0) {
            result = -PI_2 - result;  // -(π/2 + atan(x/y))
        } else {
            result = PI_2 - result;   // π/2 - atan(x/y)
        }
        return result;
    } else {
        // Use atan(y/x) approach
        Polynomial<Chebyshev> ratio = y / x;
        Polynomial<Chebyshev> result = atan(ratio);
        
        if (cx < 0.0) {
            if (cy > 0.0) {
                result = result + PI;  // atan(y/x) + π
            } else {
                result = result - PI;  // atan(y/x) - π
            }
        }
        return result;
    }
}

Polynomial<Chebyshev> exp(const Polynomial<Chebyshev>& p) {
    DACE::Interval p_bounds = p.bound();
    return evaluateOnInterval(p, [](double x) { return std::exp(x); },
                              p_bounds.m_lb, p_bounds.m_ub);
}

Polynomial<Chebyshev> log(const Polynomial<Chebyshev>& p) {
    DACE::Interval p_bounds = p.bound();
    double min_p = p_bounds.m_lb;
    double max_p = p_bounds.m_ub;
    
    // log domain is (0, inf)
    const double eps = 1e-12;
    if (min_p < eps) min_p = eps;
    if (max_p < min_p + eps) max_p = min_p + eps;
    
    return evaluateOnInterval(p, [](double x) { return std::log(x); }, min_p, max_p);
}

Polynomial<Chebyshev> sqrt(const Polynomial<Chebyshev>& p) {
    // Use algebraic Newton-Raphson iteration for sqrt
    // Iteration: x_{n+1} = 0.5 * (x_n + p/x_n)
    
    double c0 = p.cons();
    if (c0 < DACE::DA::getEps()) {
        throw std::runtime_error("Chebyshev sqrt: Argument must have positive constant part.");
    }
    
    // Initial guess: sqrt of constant part
    double sqrt_c0 = std::sqrt(c0);
    Polynomial<Chebyshev> X(sqrt_c0);
    
    // Normalize: P_norm = p / c0 (so constant part = 1)
    Polynomial<Chebyshev> P_norm = p * (1.0 / c0);
    
    // Newton-Raphson iterations: X = 0.5 * (X + P_norm/X)
    // General formula: iterations = ceil(log2(order+1)) + 2
    // Quadratic convergence requires log2(order) iterations for order-n polynomial precision
    // +2 is a safety buffer to guarantee convergence in all cases
    unsigned int max_order = DACE::DA::getMaxOrder();
    int iterations = static_cast<int>(std::ceil(std::log2(max_order + 1))) + 2;
    
    for (int i = 0; i < iterations; ++i) {
        Polynomial<Chebyshev> X_new = (X + P_norm / X) * 0.5;
        X = X_new;
    }
    
    // Denormalize: multiply by sqrt(c0)
    return X * sqrt_c0;
}

Polynomial<Chebyshev> pow(const Polynomial<Chebyshev>& p, double exponent) {
    // pow(p, n) = exp(n * log(p))
    return exp(exponent * log(p));
}

Polynomial<Chebyshev> pow(const Polynomial<Chebyshev>& base, const Polynomial<Chebyshev>& exponent) {
    // pow(a, b) = exp(b * log(a))
    return exp(exponent * log(base));
}

Polynomial<Chebyshev> sinh(const Polynomial<Chebyshev>& p) {
    // sinh(x) = (exp(x) - exp(-x)) / 2
    auto exp_p = exp(p);
    auto exp_neg_p = exp(-p);
    return (exp_p - exp_neg_p) / 2.0;
}

Polynomial<Chebyshev> cosh(const Polynomial<Chebyshev>& p) {
    // cosh(x) = (exp(x) + exp(-x)) / 2
    auto exp_p = exp(p);
    auto exp_neg_p = exp(-p);
    return (exp_p + exp_neg_p) / 2.0;
}

Polynomial<Chebyshev> tanh(const Polynomial<Chebyshev>& p) {
    return sinh(p) / cosh(p);
}

// Inverse hyperbolic functions
Polynomial<Chebyshev> asinh(const Polynomial<Chebyshev>& p) {
    // asinh(x) = log(x + sqrt(x^2 + 1))
    return log(p + sqrt(p*p + 1.0));
}

Polynomial<Chebyshev> acosh(const Polynomial<Chebyshev>& p) {
    // acosh(x) = log(x + sqrt(x^2 - 1))
    // Valid for x >= 1
    return log(p + sqrt(p*p - 1.0));
}

Polynomial<Chebyshev> atanh(const Polynomial<Chebyshev>& p) {
    // atanh(x) = 0.5 * log((1+x)/(1-x))
    // Valid for |x| < 1
    return 0.5 * log((1.0 + p) / (1.0 - p));
}

// Logarithmic variants
Polynomial<Chebyshev> log10(const Polynomial<Chebyshev>& p) {
    // log10(x) = log(x) / log(10)
    return log(p) / std::log(10.0);
}

Polynomial<Chebyshev> log2(const Polynomial<Chebyshev>& p) {
    // log2(x) = log(x) / log(2)
    return log(p) / std::log(2.0);
}

Polynomial<Chebyshev> logb(const Polynomial<Chebyshev>& p, double base) {
    // logb(x) = log(x) / log(base)
    return log(p) / std::log(base);
}

// Root functions
Polynomial<Chebyshev> root(const Polynomial<Chebyshev>& p, double n) {
    // root(x, n) = x^(1/n)
    return pow(p, 1.0 / n);
}

Polynomial<Chebyshev> cbrt(const Polynomial<Chebyshev>& p) {
    // cbrt(x) = x^(1/3)
    return pow(p, 1.0 / 3.0);
}

Polynomial<Chebyshev> icrt(const Polynomial<Chebyshev>& p) {
    // icrt(x) = x^(-1/3)
    return pow(p, -1.0 / 3.0);
}

Polynomial<Chebyshev> isrt(const Polynomial<Chebyshev>& p) {
    // isrt(x) = x^(-1/2) = 1/sqrt(x)
    return pow(p, -0.5);
}

// Utility functions
Polynomial<Chebyshev> sqr(const Polynomial<Chebyshev>& p) {
    // sqr(x) = x^2
    return p * p;
}

// Norm and analysis functions (free function versions)
double abs(const Polynomial<Chebyshev>& p) {
    return p.abs();
}

double norm(const Polynomial<Chebyshev>& p, unsigned int type) {
    return p.norm(type);
}

std::vector<double> orderNorm(const Polynomial<Chebyshev>& p, unsigned int var, unsigned int type) {
    return p.orderNorm(var, type);
}

std::vector<double> estimNorm(const Polynomial<Chebyshev>& p, unsigned int var, unsigned int type, unsigned int nc) {
    return p.estimNorm(var, type, nc);
}

std::vector<double> estimNorm(const Polynomial<Chebyshev>& p, std::vector<double> &err, unsigned int var, unsigned int type, unsigned int nc) {
    return p.estimNorm(err, var, type, nc);
}

DACE::Interval bound(const Polynomial<Chebyshev>& p) {
    return p.bound();
}

/********************************************************************************
*     AlgebraicVector<Polynomial<Chebyshev>> Inversion
*********************************************************************************/
namespace DACE {

// Helper function to invert a matrix (simplified version)
void matrix_inverse(std::vector<std::vector<double>>& A) {
    const size_t n = A.size();
    
    // Create augmented matrix [A | I]
    std::vector<std::vector<double>> aug(n, std::vector<double>(2*n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            aug[i][j] = A[i][j];
        }
        aug[i][i + n] = 1.0;
    }
    
    // Gauss-Jordan elimination
    for (size_t i = 0; i < n; ++i) {
        // Find pivot
        size_t maxRow = i;
        double maxVal = std::abs(aug[i][i]);
        for (size_t k = i + 1; k < n; ++k) {
            if (std::abs(aug[k][i]) > maxVal) {
                maxVal = std::abs(aug[k][i]);
                maxRow = k;
            }
        }
        
        // Swap rows
        if (maxRow != i) {
            std::swap(aug[i], aug[maxRow]);
        }
        
        // Check for singularity
        if (std::abs(aug[i][i]) < 1e-14) {
            throw std::runtime_error("AlgebraicVector<Chebyshev>::invert: Matrix is singular");
        }
        
        // Scale pivot row
        double pivot = aug[i][i];
        for (size_t j = 0; j < 2*n; ++j) {
            aug[i][j] /= pivot;
        }
        
        // Eliminate column
        for (size_t k = 0; k < n; ++k) {
            if (k != i) {
                double factor = aug[k][i];
                for (size_t j = 0; j < 2*n; ++j) {
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
    }
    
    // Extract inverse from augmented matrix
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A[i][j] = aug[i][j + n];
        }
    }
}

template<>
AlgebraicVector<Polynomial<Chebyshev>> AlgebraicVector<Polynomial<Chebyshev>>::invert() const {
    /*! Invert the Chebyshev polynomial map given by AlgebraicVector<Polynomial<Chebyshev>>.
       Uses native Chebyshev composition and algebra throughout.
       \return the inverted polynomial map in Chebyshev basis
       \throw std::runtime_error if dimension exceeds max variables or matrix is singular
       \note This implementation stays entirely in Chebyshev basis - never converts to Taylor
    */
    const unsigned int ord = Polynomial<Chebyshev>::getTO();
    const size_t nvar = this->size();

    if (nvar > Polynomial<Chebyshev>::getMaxVariables()) {
        throw std::runtime_error("AlgebraicVector<Chebyshev>::invert: dimension exceeds max variables");
    }

    // Create identity vector in Chebyshev basis
    AlgebraicVector<Polynomial<Chebyshev>> DDA(nvar);
    for (size_t i = 0; i < nvar; ++i) {
        DDA[i] = Polynomial<Chebyshev>(static_cast<int>(i + 1), 1.0);
    }

    // Extract constant part AC from each polynomial
    AlgebraicVector<double> AC(nvar);
    for (size_t i = 0; i < nvar; ++i) {
        AC[i] = (*this)[i].cons();
    }

    // Extract linear coefficients matrix from underlying DA objects
    // m_da stores coefficients in Chebyshev basis, but linear terms are basis-independent
    std::vector<std::vector<double>> AL(nvar, std::vector<double>(nvar, 0.0));
    for (size_t i = 0; i < nvar; ++i) {
        const DACE::DA* da_ptr = (*this)[i].getDA();
        for (size_t j = 0; j < nvar; ++j) {
            std::vector<unsigned int> jj(nvar, 0);
            jj[j] = 1;
            AL[i][j] = da_ptr->getCoefficient(jj);
        }
    }

    // Compute inverse of linear coefficient matrix
    std::vector<std::vector<double>> AI = AL;
    DACE::matrix_inverse(AI);

    // Build M = *this with constant removed
    AlgebraicVector<Polynomial<Chebyshev>> M(nvar);
    for (size_t i = 0; i < nvar; ++i) {
        M[i] = (*this)[i] - AC[i];
    }

    // Build AN = non-linear part of M (order >= 2)
    // We'll construct this by removing linear terms
    AlgebraicVector<Polynomial<Chebyshev>> AN(nvar);
    for (size_t i = 0; i < nvar; ++i) {
        // Start with M[i] and subtract linear part
        Polynomial<Chebyshev> linear_part(0.0);
        for (size_t j = 0; j < nvar; ++j) {
            if (AL[i][j] != 0.0) {
                linear_part += Polynomial<Chebyshev>(static_cast<int>(j + 1), AL[i][j]);
            }
        }
        AN[i] = M[i] - linear_part;
    }

    // Compute Linv = AI * DDA (linear map inversion)
    AlgebraicVector<Polynomial<Chebyshev>> Linv(nvar);
    for (size_t i = 0; i < nvar; ++i) {
        Linv[i] = 0.0;
        for (size_t j = 0; j < nvar; ++j) {
            Linv[i] += AI[i][j] * DDA[j];
        }
    }

    // Compute AI * AN for iterative composition
    AlgebraicVector<Polynomial<Chebyshev>> AIoAN(nvar);
    for (size_t i = 0; i < nvar; ++i) {
        AIoAN[i] = 0.0;
        for (size_t j = 0; j < nvar; ++j) {
            AIoAN[i] += AI[i][j] * AN[j];
        }
    }

    // Iteratively refine the inverse using Chebyshev composition
    AlgebraicVector<Polynomial<Chebyshev>> MI = Linv;
    for (unsigned int i = 1; i < ord; ++i) {
        Polynomial<Chebyshev>::setTO(i + 1);
        
        // Evaluate AIoAN at MI using Chebyshev composition
        AlgebraicVector<Polynomial<Chebyshev>> composed(nvar);
        for (size_t k = 0; k < nvar; ++k) {
            composed[k] = AIoAN[k].eval(MI);  // Native Chebyshev composition!
        }
        
        MI = Linv - composed;
    }

    // Final composition: MI evaluated at (DDA - AC)
    AlgebraicVector<Polynomial<Chebyshev>> shifted_identity(nvar);
    for (size_t i = 0; i < nvar; ++i) {
        shifted_identity[i] = DDA[i] - AC[i];
    }

    AlgebraicVector<Polynomial<Chebyshev>> result(nvar);
    for (size_t i = 0; i < nvar; ++i) {
        result[i] = MI[i].eval(shifted_identity);  // Native Chebyshev composition!
    }

    Polynomial<Chebyshev>::setTO(ord);  // Restore original truncation order
    return result;
}

} // namespace DACE


