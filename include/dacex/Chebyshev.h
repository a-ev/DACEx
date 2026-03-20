// Copyright (c) 2026 Adam Evans (@a-ev)

#ifndef CHEBYSHEV_H_
#define CHEBYSHEV_H_

#include "Polynomial.h"
#include <dace/Interval.h>

/** @brief Tag type for Chebyshev polynomial family specialization. */
struct Chebyshev {};

/** @brief Chebyshev-polynomial specialization of Polynomial. */
template<>
class Polynomial<Chebyshev> : public PolynomialBase {
    friend class DACE::AlgebraicVector<Polynomial<Chebyshev>>;

public:
    /**
     * @brief Construct a zero Chebyshev polynomial.
     */
    Polynomial() : PolynomialBase() {}
    /**
     * @brief Copy constructor.
     * @param p Source polynomial.
     */
    Polynomial(const Polynomial<Chebyshev> &p);
    /**
     * @brief Move constructor.
     * @param p Source polynomial.
     */
    Polynomial(Polynomial<Chebyshev> &&p) noexcept;

    /**
     * @brief Construct variable term.
     * @param i Variable index.
     * @param c Scaling coefficient.
     */
    explicit Polynomial(const int i, const double c = 1.0) : PolynomialBase(i, c) {}
    /**
     * @brief Construct variable term.
     * @param i Variable index.
     * @param c Scaling coefficient.
     */
    explicit Polynomial(const unsigned int i, const double c = 1.0) : PolynomialBase(i, c) {}
    /**
     * @brief Construct a constant polynomial.
     * @param c Constant value.
     */
    Polynomial(const double c) : PolynomialBase(c) {}
    /**
     * @brief Construct from a DACE polynomial object.
     * @param da Source DACE polynomial.
     */
    Polynomial(const DACE::DA &da) : PolynomialBase(da) {}

    /**
     * @brief Move assignment operator.
     * @param p Source polynomial.
     * @return Reference to this object.
     */
    Polynomial<Chebyshev>& operator=(Polynomial<Chebyshev> &&p) noexcept;
    /**
     * @brief Copy assignment operator.
     * @param p Source polynomial.
     * @return Reference to this object.
     */
    Polynomial<Chebyshev>& operator=(const Polynomial<Chebyshev> &p);
    /**
     * @brief Assign from scalar constant.
     * @param c Constant value.
     * @return Reference to this object.
     */
    Polynomial<Chebyshev>& operator=(const double c);

    /**
     * @brief In-place addition with another polynomial.
     * @param p Right-hand side polynomial.
     * @return Reference to this object.
     */
    Polynomial<Chebyshev>& operator+=(const Polynomial<Chebyshev> &p);
    /**
     * @brief In-place addition with scalar.
     * @param c Right-hand side scalar.
     * @return Reference to this object.
     */
    Polynomial<Chebyshev>& operator+=(const double c);
    /**
     * @brief In-place subtraction with another polynomial.
     * @param p Right-hand side polynomial.
     * @return Reference to this object.
     */
    Polynomial<Chebyshev>& operator-=(const Polynomial<Chebyshev> &p);
    /**
     * @brief In-place subtraction with scalar.
     * @param c Right-hand side scalar.
     * @return Reference to this object.
     */
    Polynomial<Chebyshev>& operator-=(const double c);
    /**
     * @brief In-place multiplication with another polynomial.
     * @param p Right-hand side polynomial.
     * @return Reference to this object.
     */
    Polynomial<Chebyshev>& operator*=(const Polynomial<Chebyshev> &p);
    /**
     * @brief In-place multiplication with scalar.
     * @param c Right-hand side scalar.
     * @return Reference to this object.
     */
    Polynomial<Chebyshev>& operator*=(const double c);
    /**
     * @brief In-place division by another polynomial.
     * @param p Right-hand side polynomial.
     * @return Reference to this object.
     */
    Polynomial<Chebyshev>& operator/=(const Polynomial<Chebyshev> &p);
    /**
     * @brief In-place division by scalar.
     * @param c Right-hand side scalar.
     * @return Reference to this object.
     */
    Polynomial<Chebyshev>& operator/=(const double c);

    /**
     * @brief Return the constant term of a polynomial.
     * @param p Input polynomial.
     * @return Constant coefficient.
     */
    friend double getConstant(const Polynomial<Chebyshev> &p);
    /**
     * @brief Compute multiplicative inverse of a polynomial.
     * @param p Input polynomial.
     * @return Inverse polynomial.
     */
    friend Polynomial<Chebyshev> inv(const Polynomial<Chebyshev> &p);

    /**
     * @brief Unary negation.
     * @return Negated polynomial.
     */
    Polynomial<Chebyshev> operator-() const;

    /** @brief Add two Chebyshev polynomials.
     * @param p1 Left polynomial.
     * @param p2 Right polynomial.
     * @return Sum polynomial.
     */
    friend Polynomial<Chebyshev> operator+(const Polynomial<Chebyshev> &p1, const Polynomial<Chebyshev> &p2);
    /** @brief Add scalar to polynomial.
     * @param p Polynomial term.
     * @param c Scalar term.
     * @return Sum polynomial.
     */
    friend Polynomial<Chebyshev> operator+(const Polynomial<Chebyshev> &p, const double c);
    /** @brief Add polynomial to scalar.
     * @param c Scalar term.
     * @param p Polynomial term.
     * @return Sum polynomial.
     */
    friend Polynomial<Chebyshev> operator+(const double c, const Polynomial<Chebyshev> &p);

    /** @brief Subtract two Chebyshev polynomials.
     * @param p1 Left polynomial.
     * @param p2 Right polynomial.
     * @return Difference polynomial.
     */
    friend Polynomial<Chebyshev> operator-(const Polynomial<Chebyshev> &p1, const Polynomial<Chebyshev> &p2);
    /** @brief Subtract scalar from polynomial.
     * @param p Polynomial term.
     * @param c Scalar term.
     * @return Difference polynomial.
     */
    friend Polynomial<Chebyshev> operator-(const Polynomial<Chebyshev> &p, const double c);
    /** @brief Subtract polynomial from scalar.
     * @param c Scalar term.
     * @param p Polynomial term.
     * @return Difference polynomial.
     */
    friend Polynomial<Chebyshev> operator-(const double c, const Polynomial<Chebyshev> &p);

    /** @brief Multiply two Chebyshev polynomials.
     * @param p1 Left polynomial.
     * @param p2 Right polynomial.
     * @return Product polynomial.
     */
    friend Polynomial<Chebyshev> operator*(const Polynomial<Chebyshev> &p1, const Polynomial<Chebyshev> &p2);
    /** @brief Multiply polynomial by scalar.
     * @param p Polynomial term.
     * @param c Scalar term.
     * @return Product polynomial.
     */
    friend Polynomial<Chebyshev> operator*(const Polynomial<Chebyshev> &p, const double c);
    /** @brief Multiply scalar by polynomial.
     * @param c Scalar term.
     * @param p Polynomial term.
     * @return Product polynomial.
     */
    friend Polynomial<Chebyshev> operator*(const double c, const Polynomial<Chebyshev> &p);

    /** @brief Divide one Chebyshev polynomial by another.
     * @param p1 Numerator polynomial.
     * @param p2 Denominator polynomial.
     * @return Quotient polynomial.
     */
    friend Polynomial<Chebyshev> operator/(const Polynomial<Chebyshev> &p1, const Polynomial<Chebyshev> &p2);
    /** @brief Divide polynomial by scalar.
     * @param p Numerator polynomial.
     * @param c Denominator scalar.
     * @return Quotient polynomial.
     */
    friend Polynomial<Chebyshev> operator/(const Polynomial<Chebyshev> &p, const double c);
    /** @brief Divide scalar by polynomial.
     * @param c Numerator scalar.
     * @param p Denominator polynomial.
     * @return Quotient polynomial.
     */
    friend Polynomial<Chebyshev> operator/(const double c, const Polynomial<Chebyshev> &p);

    /**
     * @brief Evaluate at polynomial-vector arguments.
     * @param args Algebraic-vector argument values.
     * @return Evaluated Chebyshev polynomial.
     */
    Polynomial<Chebyshev> eval(const DACE::AlgebraicVector<Polynomial<Chebyshev>>& args) const;
    /**
     * @brief Evaluate at numeric algebraic-vector arguments.
     * @param args Numeric argument values.
     * @return Evaluated scalar value.
     */
    double eval(const DACE::AlgebraicVector<double>& args) const;
    /**
     * @brief Evaluate at numeric vector arguments.
     * @param args Numeric argument values.
     * @return Evaluated scalar value.
     */
    double eval(const std::vector<double>& args) const;

    /**
     * @brief Compute a partial derivative.
     * @param var_idx 1-based variable index.
     * @return Derived polynomial.
     */
    Polynomial<Chebyshev> deriv(const unsigned int var_idx) const;
    /**
     * @brief Compute sequential mixed derivatives.
     * @param ind Derivative counts by variable.
     * @return Derived polynomial.
     */
    Polynomial<Chebyshev> deriv(const std::vector<unsigned int> ind) const;
    /**
     * @brief Compute an antiderivative.
     * @param i 1-based variable index.
     * @return Integrated polynomial.
     */
    Polynomial<Chebyshev> integ(const unsigned int i) const;
    /**
     * @brief Compute sequential mixed antiderivatives.
     * @param ind Integration counts by variable.
     * @return Integrated polynomial.
     */
    Polynomial<Chebyshev> integ(const std::vector<unsigned int> ind) const;
    
    /** @brief Maximum absolute coefficient value.
     * @return Maximum absolute coefficient.
     */
    double abs() const;
    /** @brief Norm value of selected type.
     * @param type Norm type selector.
     * @return Norm value.
     */
    double norm(const unsigned int type = 0) const;
    /** @brief Per-order norm sequence.
     * @param var Variable selector.
     * @param type Norm type selector.
     * @return Per-order norms.
     */
    std::vector<double> orderNorm(const unsigned int var = 0, const unsigned int type = 0) const;
    /** @brief Convergence norm estimate.
     * @param var Variable selector.
     * @param type Norm type selector.
     * @param nc Maximum order used in estimate.
     * @return Estimated norms.
     */
    std::vector<double> estimNorm(const unsigned int var = 0, const unsigned int type = 0, const unsigned int nc = DACE::DA::getMaxOrder()) const;
    /** @brief Convergence norm estimate with fit errors.
     * @param err Output fit-error values.
     * @param var Variable selector.
     * @param type Norm type selector.
     * @param nc Maximum order used in estimate.
     * @return Estimated norms.
     */
    std::vector<double> estimNorm(std::vector<double> &err, const unsigned int var = 0, const unsigned int type = 0, const unsigned int nc = DACE::DA::getMaxOrder()) const;
    /** @brief Interval range bound over the canonical Chebyshev domain.
     * @return Interval enclosure.
     */
    DACE::Interval bound() const;
    
    /** @brief Convert to formatted textual representation.
     * @return Human-readable polynomial string.
     */
    std::string toString() const;
    
    /** @brief Stream insertion operator.
     * @param out Output stream.
     * @param p Input polynomial.
     * @return Output stream reference.
     */
    friend std::ostream& operator<< (std::ostream &out, const Polynomial<Chebyshev> &p);
    /** @brief Stream extraction operator.
     * @param in Input stream.
     * @param p Output polynomial.
     * @return Input stream reference.
     */
    friend std::istream& operator>> (std::istream &in, Polynomial<Chebyshev> &p);
};

/** @brief Sine of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Sine polynomial.
 */
Polynomial<Chebyshev> sin(const Polynomial<Chebyshev>& p);
/** @brief Cosine of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Cosine polynomial.
 */
Polynomial<Chebyshev> cos(const Polynomial<Chebyshev>& p);
/** @brief Tangent of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Tangent polynomial.
 */
Polynomial<Chebyshev> tan(const Polynomial<Chebyshev>& p);

/** @brief Arcsine of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Arcsine polynomial.
 */
Polynomial<Chebyshev> asin(const Polynomial<Chebyshev>& p);
/** @brief Arccosine of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Arccosine polynomial.
 */
Polynomial<Chebyshev> acos(const Polynomial<Chebyshev>& p);
/** @brief Arctangent of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Arctangent polynomial.
 */
Polynomial<Chebyshev> atan(const Polynomial<Chebyshev>& p);
/** @brief Two-argument arctangent for Chebyshev polynomials.
 * @param y Y polynomial.
 * @param x X polynomial.
 * @return Angle polynomial.
 */
Polynomial<Chebyshev> atan2(const Polynomial<Chebyshev>& y, const Polynomial<Chebyshev>& x);

/** @brief Exponential of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Exponential polynomial.
 */
Polynomial<Chebyshev> exp(const Polynomial<Chebyshev>& p);
/** @brief Natural logarithm of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Logarithm polynomial.
 */
Polynomial<Chebyshev> log(const Polynomial<Chebyshev>& p);
/** @brief Square root of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Square-root polynomial.
 */
Polynomial<Chebyshev> sqrt(const Polynomial<Chebyshev>& p);
/** @brief Raise a polynomial to a scalar exponent.
 * @param p Base polynomial.
 * @param exponent Scalar exponent.
 * @return Power polynomial.
 */
Polynomial<Chebyshev> pow(const Polynomial<Chebyshev>& p, double exponent);
/** @brief Raise a polynomial to a polynomial exponent.
 * @param base Base polynomial.
 * @param exponent Exponent polynomial.
 * @return Power polynomial.
 */
Polynomial<Chebyshev> pow(const Polynomial<Chebyshev>& base, const Polynomial<Chebyshev>& exponent);

/** @brief Hyperbolic sine of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Hyperbolic-sine polynomial.
 */
Polynomial<Chebyshev> sinh(const Polynomial<Chebyshev>& p);
/** @brief Hyperbolic cosine of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Hyperbolic-cosine polynomial.
 */
Polynomial<Chebyshev> cosh(const Polynomial<Chebyshev>& p);
/** @brief Hyperbolic tangent of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Hyperbolic-tangent polynomial.
 */
Polynomial<Chebyshev> tanh(const Polynomial<Chebyshev>& p);

/** @brief Inverse hyperbolic sine of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Inverse-hyperbolic-sine polynomial.
 */
Polynomial<Chebyshev> asinh(const Polynomial<Chebyshev>& p);
/** @brief Inverse hyperbolic cosine of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Inverse-hyperbolic-cosine polynomial.
 */
Polynomial<Chebyshev> acosh(const Polynomial<Chebyshev>& p);
/** @brief Inverse hyperbolic tangent of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Inverse-hyperbolic-tangent polynomial.
 */
Polynomial<Chebyshev> atanh(const Polynomial<Chebyshev>& p);

/** @brief Base-10 logarithm of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Base-10 logarithm polynomial.
 */
Polynomial<Chebyshev> log10(const Polynomial<Chebyshev>& p);
/** @brief Base-2 logarithm of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Base-2 logarithm polynomial.
 */
Polynomial<Chebyshev> log2(const Polynomial<Chebyshev>& p);
/** @brief Logarithm with explicit base.
 * @param p Input polynomial.
 * @param base Logarithm base.
 * @return Logarithm polynomial.
 */
Polynomial<Chebyshev> logb(const Polynomial<Chebyshev>& p, double base);

/** @brief Principal nth root of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @param n Root order.
 * @return Root polynomial.
 */
Polynomial<Chebyshev> root(const Polynomial<Chebyshev>& p, double n);
/** @brief Cubic root of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Cubic-root polynomial.
 */
Polynomial<Chebyshev> cbrt(const Polynomial<Chebyshev>& p);
/** @brief Inverse cubic root of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Inverse-cubic-root polynomial.
 */
Polynomial<Chebyshev> icrt(const Polynomial<Chebyshev>& p);
/** @brief Inverse square root of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Inverse-square-root polynomial.
 */
Polynomial<Chebyshev> isrt(const Polynomial<Chebyshev>& p);

/** @brief Square of a Chebyshev polynomial.
 * @param p Input polynomial.
 * @return Squared polynomial.
 */
Polynomial<Chebyshev> sqr(const Polynomial<Chebyshev>& p);

/** @brief Maximum absolute coefficient value.
 * @param p Input polynomial.
 * @return Maximum absolute coefficient.
 */
double abs(const Polynomial<Chebyshev>& p);
/** @brief Norm value for selected type.
 * @param p Input polynomial.
 * @param type Norm type selector.
 * @return Norm value.
 */
double norm(const Polynomial<Chebyshev>& p, unsigned int type = 0);
/** @brief Per-order norm sequence.
 * @param p Input polynomial.
 * @param var Variable selector.
 * @param type Norm type selector.
 * @return Per-order norms.
 */
std::vector<double> orderNorm(const Polynomial<Chebyshev>& p, unsigned int var = 0, unsigned int type = 0);
/** @brief Convergence norm estimate.
 * @param p Input polynomial.
 * @param var Variable selector.
 * @param type Norm type selector.
 * @param nc Maximum order used in estimate.
 * @return Estimated norms.
 */
std::vector<double> estimNorm(const Polynomial<Chebyshev>& p, unsigned int var = 0, unsigned int type = 0, unsigned int nc = DACE::DA::getMaxOrder());
/** @brief Convergence estimate with fit errors.
 * @param p Input polynomial.
 * @param err Output fit-error values.
 * @param var Variable selector.
 * @param type Norm type selector.
 * @param nc Maximum order used in estimate.
 * @return Estimated norms.
 */
std::vector<double> estimNorm(const Polynomial<Chebyshev>& p, std::vector<double> &err, unsigned int var = 0, unsigned int type = 0, unsigned int nc = DACE::DA::getMaxOrder());
/** @brief Interval bound.
 * @param p Input polynomial.
 * @return Interval enclosure.
 */
DACE::Interval bound(const Polynomial<Chebyshev>& p);

/** @brief Helpers for building and evaluating Chebyshev series expansions. */
namespace ChebyshevSeries {
    /** @brief Get precomputed sine-series coefficients.
     * @param order Maximum order.
     * @return Coefficient vector.
     */
    std::vector<double> getSinCoeffs(unsigned int order);
    /** @brief Get precomputed cosine-series coefficients.
     * @param order Maximum order.
     * @return Coefficient vector.
     */
    std::vector<double> getCosCoeffs(unsigned int order);
    /** @brief Get precomputed exponential-series coefficients.
     * @param order Maximum order.
     * @return Coefficient vector.
     */
    std::vector<double> getExpCoeffs(unsigned int order);
    /** @brief Get precomputed logarithm-series coefficients.
     * @param order Maximum order.
     * @return Coefficient vector.
     */
    std::vector<double> getLogCoeffs(unsigned int order);
    /** @brief Get precomputed square-root-series coefficients.
     * @param order Maximum order.
     * @return Coefficient vector.
     */
    std::vector<double> getSqrtCoeffs(unsigned int order);
    /** @brief Get precomputed arcsine-series coefficients.
     * @param order Maximum order.
     * @return Coefficient vector.
     */
    std::vector<double> getAsinCoeffs(unsigned int order);
    /** @brief Get precomputed arctangent-series coefficients.
     * @param order Maximum order.
     * @return Coefficient vector.
     */
    std::vector<double> getAtanCoeffs(unsigned int order);
    
    /** @brief Evaluate a Chebyshev series with Clenshaw recurrence.
     * @param coeffs Series coefficients.
     * @param p Polynomial argument.
     * @return Evaluated polynomial.
     */
    Polynomial<Chebyshev> evaluateSeries(const std::vector<double>& coeffs, 
                                          const Polynomial<Chebyshev>& p);
}

namespace DACE {
    /** @brief Return the constant term of a Chebyshev polynomial.
     * @param p Input polynomial.
     * @return Constant coefficient.
     */
    inline double cons(const Polynomial<Chebyshev> &p) {
        return p.cons();
    }
    
    /** @brief Specialization of vector inversion for Chebyshev polynomial vectors. */
    template<>
    AlgebraicVector<Polynomial<Chebyshev>> AlgebraicVector<Polynomial<Chebyshev>>::invert() const;
}

#endif /* CHEBYSHEV_H_ */


