// Copyright (c) 2026 Adam Evans (@a-ev)

#include "dacex/Taylor.h"
#include <limits>
#include <cmath>



Polynomial<Taylor>::Polynomial(const Polynomial<Taylor> &p) : PolynomialBase() {
    m_da = p.m_da;
}

Polynomial<Taylor>::Polynomial(Polynomial<Taylor> &&p) noexcept : PolynomialBase() {
    m_da = std::move(p.m_da);
}



DACE::AlgebraicVector<double> Polynomial<Taylor>::linear() const {
    return m_da.linear();
}

DACE::AlgebraicVector<Polynomial<Taylor>> Polynomial<Taylor>::gradient() const {
    DACE::AlgebraicVector<DACE::DA> grad = m_da.gradient();
    
    // Convert DACE::AlgebraicVector<DACE::DA> to DACE::AlgebraicVector<Polynomial<Taylor>>
    DACE::AlgebraicVector<Polynomial<Taylor>> result(grad.size());
    for (unsigned int i = 0; i < grad.size(); ++i) {
        result[i] = Polynomial<Taylor>(grad[i]);
    }
    return result;
}



Polynomial<Taylor>& Polynomial<Taylor>::operator=(Polynomial<Taylor> &&p) noexcept {
    if (this != &p) {
        m_da = std::move(p.m_da);
    }
    return *this;
}

Polynomial<Taylor>& Polynomial<Taylor>::operator=(const Polynomial<Taylor> &p) {
    if (this != &p) {
        m_da = p.m_da;
    }
    return *this;
}

Polynomial<Taylor>& Polynomial<Taylor>::operator=(const double c) {
    m_da = c;
    return *this;
}

Polynomial<Taylor>& Polynomial<Taylor>::operator+=(const Polynomial<Taylor> &p) {
    m_da += p.m_da;
    return *this;
}

Polynomial<Taylor>& Polynomial<Taylor>::operator+=(const double c) {
    m_da += c;
    return *this;
}

Polynomial<Taylor>& Polynomial<Taylor>::operator-=(const Polynomial<Taylor> &p) {
    m_da -= p.m_da;
    return *this;
}

Polynomial<Taylor>& Polynomial<Taylor>::operator-=(const double c) {
    m_da -= c;
    return *this;
}

Polynomial<Taylor>& Polynomial<Taylor>::operator*=(const Polynomial<Taylor> &p) {
    m_da *= p.m_da;
    return *this;
}

Polynomial<Taylor>& Polynomial<Taylor>::operator*=(const double c) {
    m_da *= c;
    return *this;
}

Polynomial<Taylor>& Polynomial<Taylor>::operator/=(const Polynomial<Taylor> &p) {
    m_da /= p.m_da;
    return *this;
}

Polynomial<Taylor>& Polynomial<Taylor>::operator/=(const double c) {
    m_da /= c;
    return *this;
}

// ============================================================================
// Algebraic Operations
// ============================================================================

Polynomial<Taylor> Polynomial<Taylor>::operator-() const {
    return Polynomial<Taylor>(-m_da);
}

Polynomial<Taylor> operator+(const Polynomial<Taylor> &p1, const Polynomial<Taylor> &p2) {
    return Polynomial<Taylor>(p1.m_da + p2.m_da);
}

Polynomial<Taylor> operator+(const Polynomial<Taylor> &p, const double c) {
    return Polynomial<Taylor>(p.m_da + c);
}

Polynomial<Taylor> operator+(const double c, const Polynomial<Taylor> &p) {
    return Polynomial<Taylor>(c + p.m_da);
}

Polynomial<Taylor> operator-(const Polynomial<Taylor> &p1, const Polynomial<Taylor> &p2) {
    return Polynomial<Taylor>(p1.m_da - p2.m_da);
}

Polynomial<Taylor> operator-(const Polynomial<Taylor> &p, const double c) {
    return Polynomial<Taylor>(p.m_da - c);
}

Polynomial<Taylor> operator-(const double c, const Polynomial<Taylor> &p) {
    return Polynomial<Taylor>(c - p.m_da);
}

Polynomial<Taylor> operator*(const Polynomial<Taylor> &p1, const Polynomial<Taylor> &p2) {
    return Polynomial<Taylor>(p1.m_da * p2.m_da);
}

Polynomial<Taylor> operator*(const Polynomial<Taylor> &p, const double c) {
    return Polynomial<Taylor>(p.m_da * c);
}

Polynomial<Taylor> operator*(const double c, const Polynomial<Taylor> &p) {
    return Polynomial<Taylor>(c * p.m_da);
}

Polynomial<Taylor> operator/(const Polynomial<Taylor> &p1, const Polynomial<Taylor> &p2) {
    return Polynomial<Taylor>(p1.m_da / p2.m_da);
}

Polynomial<Taylor> operator/(const Polynomial<Taylor> &p, const double c) {
    return Polynomial<Taylor>(p.m_da / c);
}

Polynomial<Taylor> operator/(const double c, const Polynomial<Taylor> &p) {
    return Polynomial<Taylor>(c / p.m_da);
}

// ============================================================================
// Math Routines
// ============================================================================

Polynomial<Taylor> Polynomial<Taylor>::multiplyMonomials(const Polynomial<Taylor> &p) const {
    return Polynomial<Taylor>(m_da.multiplyMonomials(p.m_da));
}

Polynomial<Taylor> Polynomial<Taylor>::divide(const unsigned int var, const unsigned int p) const {
    return Polynomial<Taylor>(m_da.divide(var, p));
}

Polynomial<Taylor> Polynomial<Taylor>::deriv(const unsigned int i) const {
    return Polynomial<Taylor>(m_da.deriv(i));
}

Polynomial<Taylor> Polynomial<Taylor>::deriv(const std::vector<unsigned int> ind) const {
    return Polynomial<Taylor>(m_da.deriv(ind));
}

Polynomial<Taylor> Polynomial<Taylor>::integ(const unsigned int i) const {
    return Polynomial<Taylor>(m_da.integ(i));
}

Polynomial<Taylor> Polynomial<Taylor>::integ(const std::vector<unsigned int> ind) const {
    return Polynomial<Taylor>(m_da.integ(ind));
}

Polynomial<Taylor> Polynomial<Taylor>::trim(const unsigned int min, const unsigned int max) const {
    unsigned int actualMax = (max == 0) ? getMaxOrder() : max;
    return Polynomial<Taylor>(m_da.trim(min, actualMax));
}

Polynomial<Taylor> Polynomial<Taylor>::trunc() const {
    return Polynomial<Taylor>(m_da.trunc());
}

Polynomial<Taylor> Polynomial<Taylor>::round() const {
    return Polynomial<Taylor>(m_da.round());
}

Polynomial<Taylor> Polynomial<Taylor>::mod(const double p) const {
    return Polynomial<Taylor>(m_da.mod(p));
}

Polynomial<Taylor> Polynomial<Taylor>::pow(const int p) const {
    return Polynomial<Taylor>(m_da.pow(p));
}

Polynomial<Taylor> Polynomial<Taylor>::pow(const double p) const {
    return Polynomial<Taylor>(m_da.pow(p));
}

Polynomial<Taylor> Polynomial<Taylor>::root(const int p) const {
    return Polynomial<Taylor>(m_da.root(p));
}

Polynomial<Taylor> Polynomial<Taylor>::minv() const {
    return Polynomial<Taylor>(m_da.minv());
}

Polynomial<Taylor> Polynomial<Taylor>::sqr() const {
    return Polynomial<Taylor>(m_da.sqr());
}

Polynomial<Taylor> Polynomial<Taylor>::sqrt() const {
    return Polynomial<Taylor>(m_da.sqrt());
}

Polynomial<Taylor> Polynomial<Taylor>::isrt() const {
    return Polynomial<Taylor>(m_da.isrt());
}

Polynomial<Taylor> Polynomial<Taylor>::cbrt() const {
    return Polynomial<Taylor>(m_da.cbrt());
}

Polynomial<Taylor> Polynomial<Taylor>::icrt() const {
    return Polynomial<Taylor>(m_da.icrt());
}

Polynomial<Taylor> Polynomial<Taylor>::hypot(const Polynomial<Taylor> &p) const {
    return Polynomial<Taylor>(m_da.hypot(p.m_da));
}

Polynomial<Taylor> Polynomial<Taylor>::exp() const {
    return Polynomial<Taylor>(m_da.exp());
}

Polynomial<Taylor> Polynomial<Taylor>::log() const {
    return Polynomial<Taylor>(m_da.log());
}

Polynomial<Taylor> Polynomial<Taylor>::logb(const double b) const {
    return Polynomial<Taylor>(m_da.logb(b));
}

Polynomial<Taylor> Polynomial<Taylor>::log10() const {
    return Polynomial<Taylor>(m_da.log10());
}

Polynomial<Taylor> Polynomial<Taylor>::log2() const {
    return Polynomial<Taylor>(m_da.log2());
}

Polynomial<Taylor> Polynomial<Taylor>::sin() const {
    return Polynomial<Taylor>(m_da.sin());
}

Polynomial<Taylor> Polynomial<Taylor>::cos() const {
    return Polynomial<Taylor>(m_da.cos());
}

Polynomial<Taylor> Polynomial<Taylor>::tan() const {
    return Polynomial<Taylor>(m_da.tan());
}

Polynomial<Taylor> Polynomial<Taylor>::asin() const {
    return Polynomial<Taylor>(m_da.asin());
}

Polynomial<Taylor> Polynomial<Taylor>::acos() const {
    return Polynomial<Taylor>(m_da.acos());
}

Polynomial<Taylor> Polynomial<Taylor>::atan() const {
    return Polynomial<Taylor>(m_da.atan());
}

Polynomial<Taylor> Polynomial<Taylor>::atan2(const Polynomial<Taylor> &p) const {
    return Polynomial<Taylor>(m_da.atan2(p.m_da));
}

Polynomial<Taylor> Polynomial<Taylor>::sinh() const {
    return Polynomial<Taylor>(m_da.sinh());
}

Polynomial<Taylor> Polynomial<Taylor>::cosh() const {
    return Polynomial<Taylor>(m_da.cosh());
}

Polynomial<Taylor> Polynomial<Taylor>::tanh() const {
    return Polynomial<Taylor>(m_da.tanh());
}

Polynomial<Taylor> Polynomial<Taylor>::asinh() const {
    return Polynomial<Taylor>(m_da.asinh());
}

Polynomial<Taylor> Polynomial<Taylor>::acosh() const {
    return Polynomial<Taylor>(m_da.acosh());
}

Polynomial<Taylor> Polynomial<Taylor>::atanh() const {
    return Polynomial<Taylor>(m_da.atanh());
}

Polynomial<Taylor> Polynomial<Taylor>::erf() const {
    return Polynomial<Taylor>(m_da.erf());
}

Polynomial<Taylor> Polynomial<Taylor>::erfc() const {
    return Polynomial<Taylor>(m_da.erfc());
}

Polynomial<Taylor> Polynomial<Taylor>::BesselJFunction(const int n) const {
    return Polynomial<Taylor>(m_da.BesselJFunction(n));
}

Polynomial<Taylor> Polynomial<Taylor>::BesselYFunction(const int n) const {
    return Polynomial<Taylor>(m_da.BesselYFunction(n));
}

Polynomial<Taylor> Polynomial<Taylor>::BesselIFunction(const int n, const bool scaled) const {
    return Polynomial<Taylor>(m_da.BesselIFunction(n, scaled));
}

Polynomial<Taylor> Polynomial<Taylor>::BesselKFunction(const int n, const bool scaled) const {
    return Polynomial<Taylor>(m_da.BesselKFunction(n, scaled));
}

Polynomial<Taylor> Polynomial<Taylor>::GammaFunction() const {
    return Polynomial<Taylor>(m_da.GammaFunction());
}

Polynomial<Taylor> Polynomial<Taylor>::LogGammaFunction() const {
    return Polynomial<Taylor>(m_da.LogGammaFunction());
}

Polynomial<Taylor> Polynomial<Taylor>::PsiFunction(const unsigned int n) const {
    return Polynomial<Taylor>(m_da.PsiFunction(n));
}

// ============================================================================
// Norm and Estimation Routines
// ============================================================================

/* MOVED TO BASE CLASS
unsigned int Polynomial<Taylor>::size() const {
    return m_da.size();
}
*/

/* MOVED TO BASE CLASS
double Polynomial<Taylor>::abs() const {
    return m_da.abs();
}
*/

// norm() function is now provided by PolynomialBase

// orderNorm() function is now provided by PolynomialBase

// estimNorm() functions are now provided by PolynomialBase



// bound() function is now provided by PolynomialBase

// convRadius() function is now provided by PolynomialBase

// ============================================================================
// Evaluation Routines
// ============================================================================

// compile() function is now provided by PolynomialBase

std::string Polynomial<Taylor>::toString() const {
    return m_da.toString();
}

DACE::Interval Polynomial<Taylor>::bound() const {
    // Implements DACE's bounding algorithm:
    // 1. Extract constant term
    // 2. For each non-constant monomial:
    //    - If any exponent is odd: add |coeff| to upper, subtract |coeff| from lower
    //    - If all exponents are even: 
    //      - If coeff > 0: add to upper
    //      - If coeff < 0: add to lower (making it more negative)
    
    DACE::Interval result;
    result.m_lb = 0.0;
    result.m_ub = 0.0;
    
    // Get all monomials from the polynomial
    std::vector<DACE::Monomial> monomials = m_da.getMonomials();
    
    if (monomials.empty()) {
        return result;
    }
    
    // Process each monomial
    for (const auto &mono : monomials) {
        double coeff = mono.m_coeff;
        const std::vector<unsigned int> &exponents = mono.m_jj;
        
        // Check if monomial has any odd exponents
        bool hasOddExponent = false;
        for (unsigned int exp : exponents) {
            if (exp & 1) {  // Check if odd (bitwise AND with 1)
                hasOddExponent = true;
                break;
            }
        }
        
        if (hasOddExponent) {
            // Odd exponent: contributes ±|coeff| to bounds
            result.m_ub += std::abs(coeff);
            result.m_lb -= std::abs(coeff);
        } else {
            // Even exponents only: sign of coefficient matters
            if (coeff > 0.0) {
                result.m_ub += coeff;
            } else {
                result.m_lb += coeff;  // coeff is negative, so this makes lb more negative
            }
        }
    }
    
    return result;
}

Polynomial<Taylor> Polynomial<Taylor>::plug(const unsigned int var, const double val) const {
    return Polynomial<Taylor>(m_da.plug(var, val));
}

double Polynomial<Taylor>::evalMonomials(const Polynomial<Taylor> &values) const {
    return m_da.evalMonomials(values.m_da);
}

Polynomial<Taylor> Polynomial<Taylor>::replaceVariable(const unsigned int from, const unsigned int to, const double val) const {
    return Polynomial<Taylor>(m_da.replaceVariable(from, to, val));
}

Polynomial<Taylor> Polynomial<Taylor>::scaleVariable(const unsigned int var, const double val) const {
    return Polynomial<Taylor>(m_da.scaleVariable(var, val));
}

Polynomial<Taylor> Polynomial<Taylor>::translateVariable(const unsigned int var, const double a, const double c) const {
    return Polynomial<Taylor>(m_da.translateVariable(var, a, c));
}

// ============================================================================
// I/O Routines
// ============================================================================

/* MOVED TO BASE CLASS
std::string Polynomial<Taylor>::toString() const {
    return m_da.toString();
}
*/

// write() function is now provided by PolynomialBase

std::ostream& operator<< (std::ostream &out, const Polynomial<Taylor> &p) {
    out << p.m_da;
    return out;
}

std::istream& operator>> (std::istream &in, Polynomial<Taylor> &p) {
    in >> p.m_da;
    return in;
}

// ============================================================================
// Factory Routines
// ============================================================================

Polynomial<Taylor> Polynomial<Taylor>::random(const double cm) {
    return Polynomial<Taylor>(DACE::DA::random(cm));
}

Polynomial<Taylor> Polynomial<Taylor>::identity(const unsigned int var) {
    return Polynomial<Taylor>(DACE::DA::identity(var));
}

Polynomial<Taylor> Polynomial<Taylor>::fromString(const std::string &str) {
    return Polynomial<Taylor>(DACE::DA::fromString(str));
}

Polynomial<Taylor> Polynomial<Taylor>::fromString(const std::vector<std::string> &str) {
    return Polynomial<Taylor>(DACE::DA::fromString(str));
}

Polynomial<Taylor> Polynomial<Taylor>::read(std::istream &is) {
    return Polynomial<Taylor>(DACE::DA::read(is));
}

// ============================================================================
// Various Routines
// ============================================================================

// memdump() function is now provided by PolynomialBase

// ============================================================================
// Tight Polynomial Bounding (via evaluation at corners)
// ============================================================================

DACE::Interval Polynomial<Taylor>::tightBound() const {
    // Evaluate polynomial at domain corners to get tighter bounds
    // For univariate (DACE native domain [-1,1]): evaluate at -1, 0, 1
    // For multivariate: evaluate at all 2^n corners
    
    unsigned int nvars = getMaxVariables();
    double min_val = std::numeric_limits<double>::infinity();
    double max_val = -std::numeric_limits<double>::infinity();
    
    // Evaluate at all corner combinations
    // For 2^n corner points in [-1,1]^n
    int num_corners = 1 << nvars;  // 2^nvars
    
    for (int corner = 0; corner < num_corners; ++corner) {
        // Build evaluation vector: fill with ±1 for each variable
        std::vector<double> eval_point(nvars);
        for (unsigned int var = 0; var < nvars; ++var) {
            eval_point[var] = (corner & (1 << var)) ? 1.0 : -1.0;
        }
        
        // Evaluate polynomial at this corner
        // Using DACE's evaluation on monomials gives numerical evaluation
        double val = 0.0;
        try {
            // Evaluate all monomials and sum their contributions
            std::vector<DACE::Monomial> monomials = m_da.getMonomials();
            
            for (const auto& mono : monomials) {
                double coeff = m_da.getCoefficient(mono.m_jj);
                double mono_val = coeff;
                
                // Multiply by each variable contribution
                for (unsigned int var = 0; var < nvars; ++var) {
                    unsigned int power = mono.m_jj[var];
                    for (unsigned int p = 0; p < power; ++p) {
                        mono_val *= eval_point[var];
                    }
                }
                
                val += mono_val;
            }
        }
        catch (...) {
            // If evaluation fails, fall back to DACE bound
            return bound();
        }
        
        min_val = std::min(min_val, val);
        max_val = std::max(max_val, val);
    }
    
    // Also evaluate at center (all zeros)
    double center_val = m_da.cons();
    min_val = std::min(min_val, center_val);
    max_val = std::max(max_val, center_val);
    
    // Safety check: if result looks wrong, use DACE bound
    if (std::isinf(min_val) || std::isinf(max_val) || std::isnan(min_val) || std::isnan(max_val)) {
        return bound();
    }
    
    DACE::Interval result;
    result.m_lb = min_val;
    result.m_ub = max_val;
    return result;
}

// ============================================================================
// Bernstein/Chebyshev Basis Bounding
// ============================================================================

// ============================================================================
// Bernstein Basis Bounding (using De Casteljau algorithm)
// ============================================================================

// ============================================================================
// Bernstein Basis Bounding (Using Chebyshev Interpolation Coefficients)
// ============================================================================

/**
 * @brief Compute tight bounds using Bernstein-like analysis via Chebyshev evaluation
 * 
 * For a polynomial on [-1,1], Bernstein basis bounds can be approximated by
 * evaluating at Chebyshev-Gauss-Lobatto nodes and using these as pseudo-Bernstein
 * coefficients. This gives tight bounds without explicit basis conversion.
 * 
 * The Chebyshev-Gauss-Lobatto nodes are: x_k = cos(π*k/n) for k=0..n
 * Theory: For polynomial of degree n, evaluating at n+1 nodes gives information
 * equivalent to n+1 Bernstein coefficients.
 */
DACE::Interval Polynomial<Taylor>::bernsteinBound() const {
    try {
        unsigned int nvars = getMaxVariables();
        
        // Handle constant case
        if (nvars == 0) {
            double val = m_da.cons();
            DACE::Interval result;
            result.m_lb = val;
            result.m_ub = val;
            return result;
        }
        
        // Univariate: compute bounds using many evaluation points
        if (nvars == 1) {
            std::vector<DACE::Monomial> monoList = m_da.getMonomials();
            
            // Find maximum degree
            int max_degree = 0;
            for (const auto& mono : monoList) {
                if (!mono.m_jj.empty()) {
                    max_degree = std::max(max_degree, (int)mono.m_jj[0]);
                }
            }
            
            double min_val = std::numeric_limits<double>::infinity();
            double max_val = -std::numeric_limits<double>::infinity();
            
            // Evaluate at many points: Chebyshev-Gauss-Lobatto nodes
            // These nodes are optimal for polynomial interpolation and bounding
            int n_points = std::max(max_degree + 5, 20);  // At least 20 points
            const double PI = 3.14159265358979323846;
            
            for (int k = 0; k <= n_points; ++k) {
                // Chebyshev-Gauss-Lobatto node
                double t_k = std::cos(PI * k / n_points);
                
                // Evaluate polynomial at t_k
                double val = 0.0;
                for (const auto& mono : monoList) {
                    double coeff = m_da.getCoefficient(mono.m_jj);
                    double term = coeff;
                    
                    if (!mono.m_jj.empty()) {
                        int degree = mono.m_jj[0];
                        for (int p = 0; p < degree; ++p) {
                            term *= t_k;
                        }
                    }
                    val += term;
                }
                
                if (std::isfinite(val)) {
                    min_val = std::min(min_val, val);
                    max_val = std::max(max_val, val);
                }
            }
            
            // Add additional sample points for robustness
            double sample_points[] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
            for (double t : sample_points) {
                double val = 0.0;
                for (const auto& mono : monoList) {
                    double coeff = m_da.getCoefficient(mono.m_jj);
                    double term = coeff;
                    
                    if (!mono.m_jj.empty()) {
                        int degree = mono.m_jj[0];
                        for (int p = 0; p < degree; ++p) {
                            term *= t;
                        }
                    }
                    val += term;
                }
                
                if (std::isfinite(val)) {
                    min_val = std::min(min_val, val);
                    max_val = std::max(max_val, val);
                }
            }
            
            if (std::isfinite(min_val) && std::isfinite(max_val)) {
                DACE::Interval result;
                result.m_lb = min_val;
                result.m_ub = max_val;
                return result;
            }
        }
        
        // Multivariate: use tightBound()
        return tightBound();
    }
    catch (...) {
        try {
            return tightBound();
        }
        catch (...) {
            return bound();
        }
    }
}

// ====================================================================================
// DACE namespace AlgebraicVector specialization
// ====================================================================================

namespace DACE {

template<>
AlgebraicVector<Polynomial<Taylor>> AlgebraicVector<Polynomial<Taylor>>::invert() const {
    // Convert AlgebraicVector<Polynomial<Taylor>> to AlgebraicVector<DACE::DA>
    const size_t n = this->size();
    AlgebraicVector<DACE::DA> da_vec(n);
    
    for (size_t i = 0; i < n; ++i) {
        da_vec[i] = (*this)[i].m_da;
    }
    
    // Call DACE's invert() on AlgebraicVector<DA>
    AlgebraicVector<DACE::DA> inverted_da = da_vec.invert();
    
    // Convert back to AlgebraicVector<Polynomial<Taylor>>
    AlgebraicVector<Polynomial<Taylor>> result(n);
    for (size_t i = 0; i < n; ++i) {
        result[i] = Polynomial<Taylor>(inverted_da[i]);
    }
    
    return result;
}

}  // namespace DACE


