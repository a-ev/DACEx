// Copyright (c) 2026 Adam Evans (@a-ev)

#ifndef TAYLOR_H_
#define TAYLOR_H_

#include "Polynomial.h"

/** @brief Tag type for Taylor polynomial family specialization. */
struct Taylor {};

/** @brief Taylor-polynomial specialization of Polynomial. */
template<>
class Polynomial<Taylor> : public PolynomialBase {
    friend class DACE::AlgebraicVector<Polynomial<Taylor>>;

public:
    /**
     * @brief Bound polynomial using Bernstein basis enclosure over [-1,1]^n
     * Returns interval [min, max] of all Bernstein coefficients.
     */
    DACE::Interval bernsteinBound() const;
    /** @brief Construct a zero Taylor polynomial. */
    Polynomial() : PolynomialBase() {}
    /**
     * @brief Copy constructor.
     * @param p Source polynomial.
     */
    Polynomial(const Polynomial<Taylor> &p);
    /**
     * @brief Move constructor.
     * @param p Source polynomial.
     */
    Polynomial(Polynomial<Taylor> &&p) noexcept;
    
    /** @brief Construct variable term.
     * @param i Variable index.
     * @param c Scaling coefficient.
     */
    explicit Polynomial(const int i, const double c = 1.0) : PolynomialBase(i, c) {}
    /** @brief Construct variable term.
     * @param i Variable index.
     * @param c Scaling coefficient.
     */
    explicit Polynomial(const unsigned int i, const double c = 1.0) : PolynomialBase(i, c) {}
    /** @brief Construct a constant polynomial.
     * @param c Constant value.
     */
    Polynomial(const double c) : PolynomialBase(c) {}
    /** @brief Construct from a DACE polynomial object.
     * @param da Source DACE polynomial.
     */
    Polynomial(const DACE::DA &da) : PolynomialBase(da) {}

    /** @brief Return the first-order coefficient vector.
     * @return Linear coefficients.
     */
    DACE::AlgebraicVector<double> linear() const;
    /** @brief Return the gradient vector.
     * @return Gradient polynomial vector.
     */
    DACE::AlgebraicVector<Polynomial<Taylor>> gradient() const;

    /** @brief Move assignment operator.
     * @param p Source polynomial.
     * @return Reference to this object.
     */
    Polynomial<Taylor>& operator=(Polynomial<Taylor> &&p) noexcept;
    /** @brief Copy assignment operator.
     * @param p Source polynomial.
     * @return Reference to this object.
     */
    Polynomial<Taylor>& operator=(const Polynomial<Taylor> &p);
    /** @brief Assign from scalar constant.
     * @param c Constant value.
     * @return Reference to this object.
     */
    Polynomial<Taylor>& operator=(const double c);

    /** @brief In-place addition with another polynomial.
     * @param p Right-hand side polynomial.
     * @return Reference to this object.
     */
    Polynomial<Taylor>& operator+=(const Polynomial<Taylor> &p);
    /** @brief In-place addition with scalar.
     * @param c Right-hand side scalar.
     * @return Reference to this object.
     */
    Polynomial<Taylor>& operator+=(const double c);
    /** @brief In-place subtraction with another polynomial.
     * @param p Right-hand side polynomial.
     * @return Reference to this object.
     */
    Polynomial<Taylor>& operator-=(const Polynomial<Taylor> &p);
    /** @brief In-place subtraction with scalar.
     * @param c Right-hand side scalar.
     * @return Reference to this object.
     */
    Polynomial<Taylor>& operator-=(const double c);
    /** @brief In-place multiplication with another polynomial.
     * @param p Right-hand side polynomial.
     * @return Reference to this object.
     */
    Polynomial<Taylor>& operator*=(const Polynomial<Taylor> &p);
    /** @brief In-place multiplication with scalar.
     * @param c Right-hand side scalar.
     * @return Reference to this object.
     */
    Polynomial<Taylor>& operator*=(const double c);
    /** @brief In-place division by another polynomial.
     * @param p Right-hand side polynomial.
     * @return Reference to this object.
     */
    Polynomial<Taylor>& operator/=(const Polynomial<Taylor> &p);
    /** @brief In-place division by scalar.
     * @param c Right-hand side scalar.
     * @return Reference to this object.
     */
    Polynomial<Taylor>& operator/=(const double c);

    /** @brief Unary negation.
     * @return Negated polynomial.
     */
    Polynomial<Taylor> operator-() const;

    /** @brief Add two Taylor polynomials.
     * @param p1 Left polynomial.
     * @param p2 Right polynomial.
     * @return Sum polynomial.
     */
    friend Polynomial<Taylor> operator+(const Polynomial<Taylor> &p1, const Polynomial<Taylor> &p2);
    /** @brief Add scalar to polynomial.
     * @param p Polynomial term.
     * @param c Scalar term.
     * @return Sum polynomial.
     */
    friend Polynomial<Taylor> operator+(const Polynomial<Taylor> &p, const double c);
    /** @brief Add polynomial to scalar.
     * @param c Scalar term.
     * @param p Polynomial term.
     * @return Sum polynomial.
     */
    friend Polynomial<Taylor> operator+(const double c, const Polynomial<Taylor> &p);

    /** @brief Subtract two Taylor polynomials.
     * @param p1 Left polynomial.
     * @param p2 Right polynomial.
     * @return Difference polynomial.
     */
    friend Polynomial<Taylor> operator-(const Polynomial<Taylor> &p1, const Polynomial<Taylor> &p2);
    /** @brief Subtract scalar from polynomial.
     * @param p Polynomial term.
     * @param c Scalar term.
     * @return Difference polynomial.
     */
    friend Polynomial<Taylor> operator-(const Polynomial<Taylor> &p, const double c);
    /** @brief Subtract polynomial from scalar.
     * @param c Scalar term.
     * @param p Polynomial term.
     * @return Difference polynomial.
     */
    friend Polynomial<Taylor> operator-(const double c, const Polynomial<Taylor> &p);

    /** @brief Multiply two Taylor polynomials.
     * @param p1 Left polynomial.
     * @param p2 Right polynomial.
     * @return Product polynomial.
     */
    friend Polynomial<Taylor> operator*(const Polynomial<Taylor> &p1, const Polynomial<Taylor> &p2);
    /** @brief Multiply polynomial by scalar.
     * @param p Polynomial term.
     * @param c Scalar term.
     * @return Product polynomial.
     */
    friend Polynomial<Taylor> operator*(const Polynomial<Taylor> &p, const double c);
    /** @brief Multiply scalar by polynomial.
     * @param c Scalar term.
     * @param p Polynomial term.
     * @return Product polynomial.
     */
    friend Polynomial<Taylor> operator*(const double c, const Polynomial<Taylor> &p);

    /** @brief Divide one Taylor polynomial by another.
     * @param p1 Numerator polynomial.
     * @param p2 Denominator polynomial.
     * @return Quotient polynomial.
     */
    friend Polynomial<Taylor> operator/(const Polynomial<Taylor> &p1, const Polynomial<Taylor> &p2);
    /** @brief Divide polynomial by scalar.
     * @param p Numerator polynomial.
     * @param c Denominator scalar.
     * @return Quotient polynomial.
     */
    friend Polynomial<Taylor> operator/(const Polynomial<Taylor> &p, const double c);
    /** @brief Divide scalar by polynomial.
     * @param c Numerator scalar.
     * @param p Denominator polynomial.
     * @return Quotient polynomial.
     */
    friend Polynomial<Taylor> operator/(const double c, const Polynomial<Taylor> &p);

    /** @brief Multiply matching monomials coefficient-wise.
     * @param p Right-hand side polynomial.
     * @return Product polynomial.
     */
    Polynomial<Taylor> multiplyMonomials(const Polynomial<Taylor> &p) const;
    /** @brief Divide by a variable power when possible.
     * @param var Variable index.
     * @param p Power to divide out.
     * @return Reduced polynomial.
     */
    Polynomial<Taylor> divide(const unsigned int var, const unsigned int p = 1) const;
    /** @brief Partial derivative.
     * @param i 1-based variable index.
     * @return Derived polynomial.
     */
    Polynomial<Taylor> deriv(const unsigned int i) const;
    /** @brief Sequential mixed derivative.
     * @param ind Derivative counts by variable.
     * @return Derived polynomial.
     */
    Polynomial<Taylor> deriv(const std::vector<unsigned int> ind) const;
    /** @brief Antiderivative.
     * @param i 1-based variable index.
     * @return Integrated polynomial.
     */
    Polynomial<Taylor> integ(const unsigned int i) const;
    /** @brief Sequential mixed antiderivative.
     * @param ind Integration counts by variable.
     * @return Integrated polynomial.
     */
    Polynomial<Taylor> integ(const std::vector<unsigned int> ind) const;
    /** @brief Keep terms in an order window.
     * @param min Minimum order to keep.
     * @param max Maximum order to keep (0 means max DACE order).
     * @return Trimmed polynomial.
     */
    Polynomial<Taylor> trim(const unsigned int min, const unsigned int max = 0) const;
    /** @brief Truncate to current DACE truncation order.
     * @return Truncated polynomial.
     */
    Polynomial<Taylor> trunc() const;
    /** @brief Round coefficients to nearest integers.
     * @return Rounded polynomial.
     */
    Polynomial<Taylor> round() const;
    /** @brief Coefficient-wise modulo operation.
     * @param p Modulus value.
     * @return Modulo polynomial.
     */
    Polynomial<Taylor> mod(const double p) const;
    /** @brief Integer power.
     * @param p Integer exponent.
     * @return Powered polynomial.
     */
    Polynomial<Taylor> pow(const int p) const;
    /** @brief Real power.
     * @param p Real exponent.
     * @return Powered polynomial.
     */
    Polynomial<Taylor> pow(const double p) const;
    /** @brief Principal integer root.
     * @param p Root order.
     * @return Root polynomial.
     */
    Polynomial<Taylor> root(const int p = 2) const;
    /** @brief Multiplicative inverse.
     * @return Inverse polynomial.
     */
    Polynomial<Taylor> minv() const;
    /** @brief Square of polynomial.
     * @return Squared polynomial.
     */
    Polynomial<Taylor> sqr() const;
    /** @brief Square root of polynomial.
     * @return Square-root polynomial.
     */
    Polynomial<Taylor> sqrt() const;
    /** @brief Inverse square root of polynomial.
     * @return Inverse-square-root polynomial.
     */
    Polynomial<Taylor> isrt() const;
    /** @brief Cubic root of polynomial.
     * @return Cubic-root polynomial.
     */
    Polynomial<Taylor> cbrt() const;
    /** @brief Inverse cubic root of polynomial.
     * @return Inverse-cubic-root polynomial.
     */
    Polynomial<Taylor> icrt() const;
    /** @brief Hypotenuse-like combination.
     * @param p Other polynomial.
     * @return Result polynomial.
     */
    Polynomial<Taylor> hypot(const Polynomial<Taylor> &p) const;
    /** @brief Exponential of polynomial.
     * @return Exponential polynomial.
     */
    Polynomial<Taylor> exp() const;
    /** @brief Natural logarithm of polynomial.
     * @return Logarithm polynomial.
     */
    Polynomial<Taylor> log() const;
    /** @brief Logarithm with explicit base.
     * @param b Logarithm base.
     * @return Logarithm polynomial.
     */
    Polynomial<Taylor> logb(const double b = 10.0) const;
    /** @brief Base-10 logarithm of polynomial.
     * @return Base-10 logarithm polynomial.
     */
    Polynomial<Taylor> log10() const;
    /** @brief Base-2 logarithm of polynomial.
     * @return Base-2 logarithm polynomial.
     */
    Polynomial<Taylor> log2() const;
    /** @brief Sine of polynomial.
     * @return Sine polynomial.
     */
    Polynomial<Taylor> sin() const;
    /** @brief Cosine of polynomial.
     * @return Cosine polynomial.
     */
    Polynomial<Taylor> cos() const;
    /** @brief Tangent of polynomial.
     * @return Tangent polynomial.
     */
    Polynomial<Taylor> tan() const;
    /** @brief Arcsine of polynomial.
     * @return Arcsine polynomial.
     */
    Polynomial<Taylor> asin() const;
    /** @brief Arccosine of polynomial.
     * @return Arccosine polynomial.
     */
    Polynomial<Taylor> acos() const;
    /** @brief Arctangent of polynomial.
     * @return Arctangent polynomial.
     */
    Polynomial<Taylor> atan() const;
    /** @brief Two-argument arctangent with denominator polynomial.
     * @param p Denominator polynomial.
     * @return Angle polynomial.
     */
    Polynomial<Taylor> atan2(const Polynomial<Taylor> &p) const;
    /** @brief Hyperbolic sine of polynomial.
     * @return Hyperbolic-sine polynomial.
     */
    Polynomial<Taylor> sinh() const;
    /** @brief Hyperbolic cosine of polynomial.
     * @return Hyperbolic-cosine polynomial.
     */
    Polynomial<Taylor> cosh() const;
    /** @brief Hyperbolic tangent of polynomial.
     * @return Hyperbolic-tangent polynomial.
     */
    Polynomial<Taylor> tanh() const;
    /** @brief Inverse hyperbolic sine of polynomial.
     * @return Inverse-hyperbolic-sine polynomial.
     */
    Polynomial<Taylor> asinh() const;
    /** @brief Inverse hyperbolic cosine of polynomial.
     * @return Inverse-hyperbolic-cosine polynomial.
     */
    Polynomial<Taylor> acosh() const;
    /** @brief Inverse hyperbolic tangent of polynomial.
     * @return Inverse-hyperbolic-tangent polynomial.
     */
    Polynomial<Taylor> atanh() const;
    /** @brief Error function of polynomial.
     * @return Error-function polynomial.
     */
    Polynomial<Taylor> erf() const;
    /** @brief Complementary error function of polynomial.
     * @return Complementary-error-function polynomial.
     */
    Polynomial<Taylor> erfc() const;
    /** @brief Bessel J function.
     * @param n Order.
     * @return Bessel-J polynomial.
     */
    Polynomial<Taylor> BesselJFunction(const int n) const;
    /** @brief Bessel Y function.
     * @param n Order.
     * @return Bessel-Y polynomial.
     */
    Polynomial<Taylor> BesselYFunction(const int n) const;
    /** @brief Modified Bessel I function.
     * @param n Order.
     * @param scaled Use scaled variant when true.
     * @return Modified-Bessel-I polynomial.
     */
    Polynomial<Taylor> BesselIFunction(const int n, const bool scaled = false) const;
    /** @brief Modified Bessel K function.
     * @param n Order.
     * @param scaled Use scaled variant when true.
     * @return Modified-Bessel-K polynomial.
     */
    Polynomial<Taylor> BesselKFunction(const int n, const bool scaled = false) const;
    /** @brief Gamma function of polynomial.
     * @return Gamma polynomial.
     */
    Polynomial<Taylor> GammaFunction() const;
    /** @brief Natural logarithm of gamma function.
     * @return Log-gamma polynomial.
     */
    Polynomial<Taylor> LogGammaFunction() const;
    /** @brief Polygamma/Psi function.
     * @param n Order.
     * @return Psi polynomial.
     */
    Polynomial<Taylor> PsiFunction(const unsigned int n) const;
    
    /** @brief Convert to formatted textual representation.
     * @return Human-readable polynomial string.
     */
    std::string toString() const;

    /** @brief Compute a DACE interval bound.
     * @return Interval enclosure.
     */
    DACE::Interval bound() const;
    /** @brief Compute a tighter bound via additional evaluation.
     * @return Interval enclosure.
     */
    DACE::Interval tightBound() const;
    /** @brief Substitute a value into one variable.
     * @param var Variable index.
     * @param val Replacement value.
     * @return Substituted polynomial.
     */
    Polynomial<Taylor> plug(const unsigned int var, const double val = 0.0) const;
    /** @brief Evaluate monomials against polynomial-valued inputs.
     * @param values Polynomial argument values.
     * @return Evaluated scalar result.
     */
    double evalMonomials(const Polynomial<Taylor> &values) const;
    /** @brief Remap one variable index to another.
     * @param from Source variable index.
     * @param to Destination variable index.
     * @param val Scaling value.
     * @return Remapped polynomial.
     */
    Polynomial<Taylor> replaceVariable(const unsigned int from = 0, const unsigned int to = 0, const double val = 1.0) const;
    /** @brief Scale one variable.
     * @param var Variable index.
     * @param val Scaling factor.
     * @return Scaled polynomial.
     */
    Polynomial<Taylor> scaleVariable(const unsigned int var = 0, const double val = 1.0) const;
    /** @brief Affine translate one variable using x -> a*x + c.
     * @param var Variable index.
     * @param a Multiplicative factor.
     * @param c Additive constant.
     * @return Translated polynomial.
     */
    Polynomial<Taylor> translateVariable(const unsigned int var = 0, const double a = 1.0, const double c = 0.0) const;

    /** @brief Stream insertion operator.
     * @param out Output stream.
     * @param p Input polynomial.
     * @return Output stream reference.
     */
    friend std::ostream& operator<< (std::ostream &out, const Polynomial<Taylor> &p);
    /** @brief Stream extraction operator.
     * @param in Input stream.
     * @param p Output polynomial.
     * @return Input stream reference.
     */
    friend std::istream& operator>> (std::istream &in, Polynomial<Taylor> &p);

    /** @brief Create a random polynomial.
     * @param cm Coefficient magnitude scale.
     * @return Random polynomial.
     */
    static Polynomial<Taylor> random(const double cm);
    /** @brief Create identity polynomial for one variable.
     * @param var Variable index.
     * @return Identity polynomial.
     */
    static Polynomial<Taylor> identity(const unsigned int var);
    /** @brief Parse a polynomial from a single textual representation.
     * @param str Input text.
     * @return Parsed polynomial.
     */
    static Polynomial<Taylor> fromString(const std::string &str);
    /** @brief Parse a polynomial from multiple textual lines.
     * @param str Input text lines.
     * @return Parsed polynomial.
     */
    static Polynomial<Taylor> fromString(const std::vector<std::string> &str);
    /** @brief Read a polynomial from a stream.
     * @param is Input stream.
     * @return Parsed polynomial.
     */
    static Polynomial<Taylor> read(std::istream &is);
};

/** @brief Free-function exponential wrapper.
 * @param p Input polynomial.
 * @return Exponential polynomial.
 */
inline Polynomial<Taylor> exp(const Polynomial<Taylor> &p) { return p.exp(); }
/** @brief Free-function natural logarithm wrapper.
 * @param p Input polynomial.
 * @return Logarithm polynomial.
 */
inline Polynomial<Taylor> log(const Polynomial<Taylor> &p) { return p.log(); }
/** @brief Free-function base-10 logarithm wrapper.
 * @param p Input polynomial.
 * @return Base-10 logarithm polynomial.
 */
inline Polynomial<Taylor> log10(const Polynomial<Taylor> &p) { return p.log10(); }
/** @brief Free-function base-2 logarithm wrapper.
 * @param p Input polynomial.
 * @return Base-2 logarithm polynomial.
 */
inline Polynomial<Taylor> log2(const Polynomial<Taylor> &p) { return p.log2(); }
/** @brief Free-function sine wrapper.
 * @param p Input polynomial.
 * @return Sine polynomial.
 */
inline Polynomial<Taylor> sin(const Polynomial<Taylor> &p) { return p.sin(); }
/** @brief Free-function cosine wrapper.
 * @param p Input polynomial.
 * @return Cosine polynomial.
 */
inline Polynomial<Taylor> cos(const Polynomial<Taylor> &p) { return p.cos(); }
/** @brief Free-function tangent wrapper.
 * @param p Input polynomial.
 * @return Tangent polynomial.
 */
inline Polynomial<Taylor> tan(const Polynomial<Taylor> &p) { return p.tan(); }
/** @brief Free-function arcsine wrapper.
 * @param p Input polynomial.
 * @return Arcsine polynomial.
 */
inline Polynomial<Taylor> asin(const Polynomial<Taylor> &p) { return p.asin(); }
/** @brief Free-function arccosine wrapper.
 * @param p Input polynomial.
 * @return Arccosine polynomial.
 */
inline Polynomial<Taylor> acos(const Polynomial<Taylor> &p) { return p.acos(); }
/** @brief Free-function arctangent wrapper.
 * @param p Input polynomial.
 * @return Arctangent polynomial.
 */
inline Polynomial<Taylor> atan(const Polynomial<Taylor> &p) { return p.atan(); }
/** @brief Free-function hyperbolic sine wrapper.
 * @param p Input polynomial.
 * @return Hyperbolic-sine polynomial.
 */
inline Polynomial<Taylor> sinh(const Polynomial<Taylor> &p) { return p.sinh(); }
/** @brief Free-function hyperbolic cosine wrapper.
 * @param p Input polynomial.
 * @return Hyperbolic-cosine polynomial.
 */
inline Polynomial<Taylor> cosh(const Polynomial<Taylor> &p) { return p.cosh(); }
/** @brief Free-function hyperbolic tangent wrapper.
 * @param p Input polynomial.
 * @return Hyperbolic-tangent polynomial.
 */
inline Polynomial<Taylor> tanh(const Polynomial<Taylor> &p) { return p.tanh(); }
/** @brief Free-function inverse hyperbolic sine wrapper.
 * @param p Input polynomial.
 * @return Inverse-hyperbolic-sine polynomial.
 */
inline Polynomial<Taylor> asinh(const Polynomial<Taylor> &p) { return p.asinh(); }
/** @brief Free-function inverse hyperbolic cosine wrapper.
 * @param p Input polynomial.
 * @return Inverse-hyperbolic-cosine polynomial.
 */
inline Polynomial<Taylor> acosh(const Polynomial<Taylor> &p) { return p.acosh(); }
/** @brief Free-function inverse hyperbolic tangent wrapper.
 * @param p Input polynomial.
 * @return Inverse-hyperbolic-tangent polynomial.
 */
inline Polynomial<Taylor> atanh(const Polynomial<Taylor> &p) { return p.atanh(); }
/** @brief Free-function error function wrapper.
 * @param p Input polynomial.
 * @return Error-function polynomial.
 */
inline Polynomial<Taylor> erf(const Polynomial<Taylor> &p) { return p.erf(); }
/** @brief Free-function complementary error-function wrapper.
 * @param p Input polynomial.
 * @return Complementary-error-function polynomial.
 */
inline Polynomial<Taylor> erfc(const Polynomial<Taylor> &p) { return p.erfc(); }
/** @brief Free-function square-root wrapper.
 * @param p Input polynomial.
 * @return Square-root polynomial.
 */
inline Polynomial<Taylor> sqrt(const Polynomial<Taylor> &p) { return p.sqrt(); }
/** @brief Free-function cubic-root wrapper.
 * @param p Input polynomial.
 * @return Cubic-root polynomial.
 */
inline Polynomial<Taylor> cbrt(const Polynomial<Taylor> &p) { return p.cbrt(); }

namespace DACE {
    /** @brief Return the constant term of a Taylor polynomial.
     * @param p Input polynomial.
     * @return Constant coefficient.
     */
    inline double cons(const Polynomial<Taylor> &p) {
        return p.cons();
    }
    
    /** @brief Specialization of vector inversion for Taylor polynomial vectors. */
    template<>
    AlgebraicVector<Polynomial<Taylor>> AlgebraicVector<Polynomial<Taylor>>::invert() const;
}

#endif /* TAYLOR_H_ */

