// Copyright (c) 2026 Adam Evans (@a-ev)

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include "dace/dace.h"
#include <iostream>
#include <vector>
#include <string>
#include <initializer_list>
#include <cmath>
#include <functional>

/**
 * @brief Common base class that wraps a DACE polynomial object.
 *
 * This class exposes setup, coefficient access, evaluation, and utility
 * methods that are shared by specific polynomial families.
 */
class PolynomialBase {
protected:
    DACE::DA m_da;

public:
    /**
     * @brief Initialize the DACE engine with expansion order and variable count.
     * @param ord Maximum polynomial order.
     * @param nvar Number of independent variables.
     */
    static void init(const unsigned int ord, const unsigned int nvar) {
        DACE::DA::init(ord, nvar);
    }
    
    /**
     * @brief Check whether the DACE engine is initialized.
     * @return True when DACE is initialized, otherwise false.
     */
    static bool isInitialized() {
        return DACE::DA::isInitialized();
    }
    
    /**
     * @brief Get DACE version numbers.
     * @param maj Output major version.
     * @param min Output minor version.
     * @param patch Output patch version.
     */
    static void version(int &maj, int &min, int &patch) {
        DACE::DA::version(maj, min, patch);
    }
    
    /** @brief Validate runtime and compile-time DACE version compatibility. */
    static void checkVersion() {
        DACE::DA::checkVersion();
    }
    
    /**
     * @brief Set DACE truncation epsilon.
     * @param eps New epsilon value.
     * @return Previous epsilon value.
     */
    static double setEps(const double eps) {
        return DACE::DA::setEps(eps);
    }
    
    /**
     * @brief Get the current DACE truncation epsilon.
     * @return Current epsilon value.
     */
    static double getEps() {
        return DACE::DA::getEps();
    }
    
    /**
     * @brief Get the machine epsilon used by DACE.
     * @return Machine epsilon value.
     */
    static double getEpsMac() {
        return DACE::DA::getEpsMac();
    }
    
    /**
     * @brief Get the maximum order configured in DACE.
     * @return Maximum polynomial order.
     */
    static unsigned int getMaxOrder() {
        return DACE::DA::getMaxOrder();
    }
    
    /**
     * @brief Get the maximum number of variables supported by DACE.
     * @return Maximum variable count.
     */
    static unsigned int getMaxVariables() {
        return DACE::DA::getMaxVariables();
    }
    
    /**
     * @brief Get the maximum number of monomials supported by DACE.
     * @return Maximum monomial count.
     */
    static unsigned int getMaxMonomials() {
        return DACE::DA::getMaxMonomials();
    }
    
    /**
     * @brief Get the current truncation order.
     * @return Current truncation order.
     */
    static unsigned int getTO() {
        return DACE::DA::getTO();
    }
    
    /**
     * @brief Set the truncation order.
     * @param ot New truncation order; when zero, uses the maximum order.
     * @return The previous truncation order.
     */
    static unsigned int setTO(const unsigned int ot = 0) {
        if (ot == 0) {
            return DACE::DA::setTO(getMaxOrder());
        }
        return DACE::DA::setTO(ot);
    }
    
    /**
     * @brief Push a truncation order on the DACE TO stack.
     * @param ot Truncation order to push; when zero, uses the maximum order.
     */
    static void pushTO(const unsigned int ot = 0) {
        if (ot == 0) {
            DACE::DA::pushTO(getMaxOrder());
        } else {
            DACE::DA::pushTO(ot);
        }
    }
    
    /** @brief Restore the previous truncation order from the DACE TO stack. */
    static void popTO() {
        DACE::DA::popTO();
    }
    
    /** @brief Print DACE memory usage diagnostics. */
    static void memdump() {
        DACE::DA::memdump();
    }

    /** @brief Construct a zero polynomial. */
    PolynomialBase() : m_da() {}
    
    /**
     * @brief Construct variable term.
     * @param i Variable index.
     * @param c Variable scaling coefficient.
     */
    explicit PolynomialBase(const int i, const double c = 1.0) : m_da(i, c) {}
    
    /**
     * @brief Construct variable term.
     * @param i Variable index.
     * @param c Variable scaling coefficient.
     */
    explicit PolynomialBase(const unsigned int i, const double c = 1.0) : m_da(i, c) {}
    
    /**
     * @brief Construct a constant polynomial.
     * @param c Constant value.
     */
    PolynomialBase(const double c) : m_da(c) {}
    
    /**
     * @brief Construct from an existing DACE polynomial object.
     * @param da Source DACE object.
     */
    PolynomialBase(const DACE::DA &da) : m_da(da) {}
    
    /** @brief Default destructor. */
    ~PolynomialBase() noexcept = default;

    /**
     * @brief Check whether any coefficient is NaN.
     * @return Non-zero when a NaN is present.
     */
    int isnan() const {
        return m_da.isnan();
    }
    
    /**
     * @brief Check whether any coefficient is infinite.
     * @return Non-zero when an infinite value is present.
     */
    int isinf() const {
        return m_da.isinf();
    }
    
    /**
     * @brief Get the constant term.
     * @return Constant coefficient.
     */
    double cons() const {
        return m_da.cons();
    }
    
    /**
     * @brief Get the coefficient for a monomial index.
     * @param jj Exponent/index vector.
     * @return Coefficient value at the requested index.
     */
    double getCoefficient(const std::vector<unsigned int> &jj) const {
        return m_da.getCoefficient(jj);
    }
    
    /**
     * @brief Set the coefficient for a monomial index.
     * @param jj Exponent/index vector.
     * @param coeff New coefficient value.
     */
    void setCoefficient(const std::vector<unsigned int> &jj, const double coeff) {
        m_da.setCoefficient(jj, coeff);
    }
    
    /**
     * @brief Get a monomial by linear storage position.
     * @param npos Monomial position.
     * @return Monomial object at that position.
     */
    DACE::Monomial getMonomial(const unsigned int npos) const {
        return m_da.getMonomial(npos);
    }
    
    /**
     * @brief Write a monomial by linear storage position.
     * @param npos Monomial position.
     * @param m Output monomial storage.
     */
    void getMonomial(const unsigned int npos, DACE::Monomial &m) const {
        m_da.getMonomial(npos, m);
    }
    
    /**
     * @brief Get all monomials present in the polynomial.
     * @return Vector of monomial entries.
     */
    std::vector<DACE::Monomial> getMonomials() const {
        return m_da.getMonomials();
    }

    /**
     * @brief Get the number of non-zero monomials.
     * @return Number of stored monomials.
     */
    unsigned int size() const {
        return m_da.size();
    }
    
    /**
     * @brief Get the maximum absolute coefficient.
     * @return Maximum absolute coefficient value.
     */
    double abs() const {
        return m_da.abs();
    }
    
    /**
     * @brief Compute a norm for the polynomial.
     * @param type Norm type selector.
     * @return Requested norm value.
     */
    double norm(const unsigned int type = 0) const {
        return m_da.norm(type);
    }
    
    /**
     * @brief Compute per-order norm values.
     * @param var Variable selector (0 for total order).
     * @param type Norm type selector.
     * @return Norm values grouped by order.
     */
    std::vector<double> orderNorm(const unsigned int var = 0, const unsigned int type = 0) const {
        return m_da.orderNorm(var, type);
    }
    
    /**
     * @brief Estimate convergence norms.
     * @param var Variable selector (0 for total order).
     * @param type Norm type selector.
     * @param nc Maximum order for the estimate.
     * @return Estimated norms.
     */
    std::vector<double> estimNorm(const unsigned int var = 0, const unsigned int type = 0, const unsigned int nc = 0) const {
        return m_da.estimNorm(var, type, nc);
    }
    
    /**
     * @brief Estimate convergence norms and fit errors.
     * @param err Output fitting-error values.
     * @param var Variable selector (0 for total order).
     * @param type Norm type selector.
     * @param nc Maximum order for the estimate.
     * @return Estimated norms.
     */
    std::vector<double> estimNorm(std::vector<double> &err, const unsigned int var = 0, const unsigned int type = 0, const unsigned int nc = 0) const {
        return m_da.estimNorm(err, var, type, nc);
    }
    
    /**
     * @brief Compute an interval bound over the default DACE domain.
     * @return Interval enclosure of polynomial values.
     */
    DACE::Interval bound() const {
        return m_da.bound();
    }
    
    /**
     * @brief Estimate convergence radius.
     * @param eps Tolerance target.
     * @param type Radius-estimation mode.
     * @return Estimated convergence radius.
     */
    double convRadius(const double eps, const unsigned int type = 1) const {
        return m_da.convRadius(eps, type);
    }

    /**
     * @brief Evaluate the polynomial at vector arguments.
     * @tparam T Evaluation scalar type.
     * @param args Input argument vector.
     * @return Evaluated value in type @p T.
     */
    template<class T> T eval(const std::vector<T> &args) const {
        return m_da.eval(args);
    }
    
    /**
     * @brief Evaluate the polynomial at array arguments.
     * @tparam T Evaluation scalar type.
     * @param args Input array pointer.
     * @param length Number of array entries.
     * @return Evaluated value in type @p T.
     */
    template<class T> T eval(const T args[], const unsigned int length) const {
        return m_da.eval(args, length);
    }
    
    /**
     * @brief Evaluate a single-variable polynomial.
     * @tparam T Evaluation scalar type.
     * @param arg Input scalar argument.
     * @return Evaluated value in type @p T.
     */
    template<class T> T evalScalar(const T &arg) const {
        return m_da.evalScalar(arg);
    }
    
    /**
     * @brief Compile the polynomial for faster repeated evaluation.
     * @return Compiled DACE representation.
     */
    DACE::compiledDA compile() const {
        return m_da.compile();
    }

    /**
     * @brief Serialize the polynomial to an output stream.
     * @param os Output stream.
     */
    void write(std::ostream &os) const {
        m_da.write(os);
    }

    /**
     * @brief Access the internal DACE object pointer.
     * @return Pointer to the wrapped DACE object.
     */
    DACE::DA* getDA() const {
        return const_cast<DACE::DA*>(&m_da);
    }
};

/** @brief Primary polynomial template, specialized by family tags. */
template<typename Family> class Polynomial;
/** @brief Alias template declaration for algebraic vectors. */
template<typename T> class AlgebraicVector;

#endif /* POLYNOMIAL_H_ */

