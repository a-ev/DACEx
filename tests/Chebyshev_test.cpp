#include "dacex/DACEx.h"
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include <cmath> 
#include <map>
#include <functional>

using namespace DACE;

static int g_failures = 0;

// =============================================================
//  Visual Styling
// =============================================================
const std::string CYAN  = "\033[36m";
const std::string GREEN = "\033[32m";
const std::string RED   = "\033[31m";
const std::string YELLOW = "\033[33m";
const std::string RESET = "\033[0m";

void printHeader(const std::string& title) {
    std::cout << "\n" << CYAN << "==========================================================" << RESET << std::endl;
    std::cout << "  " << title << std::endl;
    std::cout << CYAN << "==========================================================" << RESET << std::endl;
}

void printSubHeader(const std::string& title) {
    std::cout << "\n" << YELLOW << "--- " << title << " ---" << RESET << std::endl;
}

// =============================================================
//  Test Result Verification
// =============================================================

void verifyValue(const std::string& testName, double actual, double expected, double tolerance = 1e-10) {
    double error = std::abs(actual - expected);
    std::cout << "  " << testName << ": ";
    if (error <= tolerance) {
        std::cout << GREEN << "[PASS]" << RESET << " Actual: " << actual 
                  << " (Expected: " << expected << ")" << std::endl;
    } else {
        ++g_failures;
        std::cout << RED << "[FAIL]" << RESET << " Actual: " << actual 
                  << " Expected: " << expected << " Error: " << error << std::endl;
    }
}

void verifyCoefficient(const Polynomial<Chebyshev>& poly, 
                       const std::vector<unsigned int>& indices,
                       double expected, 
                       double tolerance = 1e-10) {
    double actual = poly.getCoefficient(indices);
    double error = std::abs(actual - expected);
    
    std::cout << "  T(";
    for(size_t i = 0; i < indices.size(); ++i) {
        std::cout << indices[i];
        if(i < indices.size() - 1) std::cout << ",";
    }
    std::cout << "): ";
    
    if (error <= tolerance) {
        std::cout << GREEN << "[PASS]" << RESET << " Actual: " << std::scientific 
                  << std::setprecision(6) << actual << std::fixed << std::endl;
    } else {
        ++g_failures;
        std::cout << RED << "[FAIL]" << RESET << " Actual: " << std::scientific 
                  << std::setprecision(6) << actual << " Expected: " << expected 
                  << " Error: " << error << std::fixed << std::endl;
    }
}

// =============================================================
//  Test Categories
// =============================================================

void testBasicArithmetic() {
    printHeader("Basic Arithmetic Operations");
    
    Polynomial<Chebyshev> x = DA(1);
    Polynomial<Chebyshev> y = DA(2);
    
    printSubHeader("Addition and Subtraction");
    {
        Polynomial<Chebyshev> p1 = 3.0 + 2.0*x;
        Polynomial<Chebyshev> p2 = 1.0 + x;
        Polynomial<Chebyshev> sum = p1 + p2;
        Polynomial<Chebyshev> diff = p1 - p2;
        
        std::cout << "\n(3 + 2x) + (1 + x) = 4 + 3x" << std::endl;
        verifyCoefficient(sum, {0,0}, 4.0);
        verifyCoefficient(sum, {1,0}, 3.0);
        
        std::cout << "\n(3 + 2x) - (1 + x) = 2 + x" << std::endl;
        verifyCoefficient(diff, {0,0}, 2.0);
        verifyCoefficient(diff, {1,0}, 1.0);
    }
    
    printSubHeader("Multiplication");
    {
        // T_m * T_n = 0.5*(T_{m+n} + T_{|m-n|})
        Polynomial<Chebyshev> T0 = 1.0;
        Polynomial<Chebyshev> T1 = x;
        Polynomial<Chebyshev> T2 = 2.0*x*T1 - T0;
        
        std::cout << "\nT1 * T1 = 0.5*(T2 + T0)" << std::endl;
        Polynomial<Chebyshev> T1_sq = T1 * T1;
        verifyCoefficient(T1_sq, {0,0}, 0.5);
        verifyCoefficient(T1_sq, {2,0}, 0.5);
        
        std::cout << "\nT1 * T2 = 0.5*(T3 + T1)" << std::endl;
        Polynomial<Chebyshev> T1_T2 = T1 * T2;
        verifyCoefficient(T1_T2, {1,0}, 0.5);
        verifyCoefficient(T1_T2, {3,0}, 0.5);
        
        std::cout << "\n(x + y)^2" << std::endl;
        Polynomial<Chebyshev> sq = (x + y) * (x + y);
        verifyCoefficient(sq, {0,0}, 1.0);
        verifyCoefficient(sq, {2,0}, 0.5);
        verifyCoefficient(sq, {0,2}, 0.5);
        verifyCoefficient(sq, {1,1}, 2.0);
    }
    
    printSubHeader("Division and Inversion");
    {
        Polynomial<Chebyshev> p = 1.0 - 0.3*x;
        Polynomial<Chebyshev> invP = Polynomial<Chebyshev>(1.0) / p;
        Polynomial<Chebyshev> identity = p * invP;
        
        std::cout << "\nP * (1/P) ≈ 1" << std::endl;
        verifyCoefficient(identity, {0,0}, 1.0, 1e-9);
        
        // Check that higher order terms are near zero
        double higherOrderSum = 0.0;
        for(unsigned int i = 1; i <= 4; ++i) {
            higherOrderSum += std::abs(identity.getCoefficient({i,0}));
        }
        verifyValue("Higher order terms sum", higherOrderSum, 0.0, 1e-9);
    }
    
    printSubHeader("Scalar Operations");
    {
        Polynomial<Chebyshev> p = 2.0 + 3.0*x;
        Polynomial<Chebyshev> scaled = p * 2.0;
        Polynomial<Chebyshev> divided = p / 2.0;
        
        std::cout << "\n2*(2 + 3x) = 4 + 6x" << std::endl;
        verifyCoefficient(scaled, {0,0}, 4.0);
        verifyCoefficient(scaled, {1,0}, 6.0);
        
        std::cout << "\n(2 + 3x)/2 = 1 + 1.5x" << std::endl;
        verifyCoefficient(divided, {0,0}, 1.0);
        verifyCoefficient(divided, {1,0}, 1.5);
    }
}

void testCalculus() {
    printHeader("Calculus Operations");
    
    Polynomial<Chebyshev> x = DA(1);
    Polynomial<Chebyshev> y = DA(2);
    Polynomial<Chebyshev> T0 = 1.0;
    Polynomial<Chebyshev> T1 = x;
    Polynomial<Chebyshev> T2 = 2.0*x*T1 - T0;
    Polynomial<Chebyshev> T3 = 2.0*x*T2 - T1;
    
    printSubHeader("Differentiation");
    {
        // d/dx T2 = 4*T1
        std::cout << "\nd/dx T2(x) = 4*T1(x)" << std::endl;
        Polynomial<Chebyshev> dT2 = T2.deriv(1);
        verifyCoefficient(dT2, {1,0}, 4.0);
        verifyCoefficient(dT2, {0,0}, 0.0);
        
        // d/dx T3 = 6*(T2 + 0.5*T0)
        std::cout << "\nd/dx T3(x) = 6*T2(x) + 3*T0" << std::endl;
        Polynomial<Chebyshev> dT3 = T3.deriv(1);
        verifyCoefficient(dT3, {0,0}, 3.0);
        verifyCoefficient(dT3, {2,0}, 6.0);
        
        // d/dx (T1*T2) = d/dx[0.5*(T3 + T1)] = 3*T2 + 2*T0
        std::cout << "\nd/dx (T1*T2) = 3*T2 + 2*T0" << std::endl;
        Polynomial<Chebyshev> T1_T2 = T1 * T2;
        Polynomial<Chebyshev> d_T1_T2 = T1_T2.deriv(1);
        verifyCoefficient(d_T1_T2, {0,0}, 2.0);
        verifyCoefficient(d_T1_T2, {2,0}, 3.0);
        
        // Mixed partial: d²/dxdy (xy) = 1
        std::cout << "\nd²/dxdy (T1(x)*T1(y)) = 1" << std::endl;
        Polynomial<Chebyshev> xy = x * y;
        Polynomial<Chebyshev> dxy_dx = xy.deriv(1);
        Polynomial<Chebyshev> dxy_dxdy = dxy_dx.deriv(2);
        verifyCoefficient(dxy_dxdy, {0,0}, 1.0);
    }
    
    printSubHeader("Integration");
    {
        // ∫T2 dx = T3/6 - T1/2
        std::cout << "\n∫T2 dx = T3/6 - T1/2" << std::endl;
        Polynomial<Chebyshev> intT2 = T2.integ(1);
        verifyCoefficient(intT2, {1,0}, -0.5);
        verifyCoefficient(intT2, {3,0}, 1.0/6.0);
        
        // ∫T1 dx = T2/4 + T0/4 (with constant of integration)
        std::cout << "\n∫T1 dx = T2/4 + T0/4" << std::endl;
        Polynomial<Chebyshev> intT1 = T1.integ(1);
        verifyCoefficient(intT1, {2,0}, 0.25);
        verifyCoefficient(intT1, {0,0}, 0.25);  // Constant of integration
    }
}

void testTrigonometric() {
    printHeader("Trigonometric Functions");
    
    Polynomial<Chebyshev> x = DA(1);
    Polynomial<Chebyshev> y = DA(2);
    
    printSubHeader("sin, cos, tan");
    {
        // Test 1: sin on polynomial with T2 component
        std::cout << "\nTest 1: sin(0.3 + 0.2*T2(x)) at x=0.5" << std::endl;
        Polynomial<Chebyshev> T2 = 2.0*x*x - 1.0;
        Polynomial<Chebyshev> p1 = 0.3 + 0.2*T2;  // Higher order polynomial
        Polynomial<Chebyshev> sinP1 = sin(p1);
        std::vector<double> args1 = {0.5, 0.0};
        double pVal1 = p1.eval(args1);
        double result1 = sinP1.eval(args1);
        double expected1 = std::sin(pVal1);
        verifyValue("sin(P) with T2", result1, expected1, 0.005);
        
        // Test 2: sin on bivariate polynomial
        std::cout << "\nTest 2: sin(0.2*T1(x) + 0.1*T1(y)) at x=0.6, y=0.4" << std::endl;
        Polynomial<Chebyshev> p2 = 0.2*x + 0.1*y;
        Polynomial<Chebyshev> sinP2 = sin(p2);
        std::vector<double> args2 = {0.6, 0.4};
        double pVal2 = p2.eval(args2);
        double result2 = sinP2.eval(args2);
        double expected2 = std::sin(pVal2);
        verifyValue("sin(P) bivariate", result2, expected2, 0.005);
        
        // Test 3: cos on polynomial with T3 component
        std::cout << "\nTest 3: cos(0.4 + 0.15*T3(x)) at x=0.7" << std::endl;
        Polynomial<Chebyshev> T3 = 4.0*x*x*x - 3.0*x;
        Polynomial<Chebyshev> p3 = 0.4 + 0.15*T3;
        Polynomial<Chebyshev> cosP3 = cos(p3);
        std::vector<double> args3 = {0.7, 0.0};
        double pVal3 = p3.eval(args3);
        double result3 = cosP3.eval(args3);
        double expected3 = std::cos(pVal3);
        verifyValue("cos(P) with T3", result3, expected3, 0.005);
        
        // Test 4: tan with mixed terms
        std::cout << "\nTest 4: tan(0.25*T2(x) + 0.1*T1(y)) at x=0.5, y=0.3" << std::endl;
        Polynomial<Chebyshev> p4 = 0.25*T2 + 0.1*y;
        Polynomial<Chebyshev> tanP4 = tan(p4);
        std::vector<double> args4 = {0.5, 0.3};
        double pVal4 = p4.eval(args4);
        double result4 = tanP4.eval(args4);
        double expected4 = std::tan(pVal4);
        verifyValue("tan(P) mixed terms", result4, expected4, 0.005);
        
        // Identity test: sin²+ cos² = 1
        std::cout << "\nTest 5: sin²(P) + cos²(P) = 1 for P with T2" << std::endl;
        Polynomial<Chebyshev> p5 = 0.3 + 0.15*T2;
        Polynomial<Chebyshev> sinP5 = sin(p5);
        Polynomial<Chebyshev> cosP5 = cos(p5);
        Polynomial<Chebyshev> sin2 = sinP5 * sinP5;
        Polynomial<Chebyshev> cos2 = cosP5 * cosP5;
        Polynomial<Chebyshev> sum = sin2 + cos2;
        std::vector<double> args5 = {0.6, 0.0};
        double result5 = sum.eval(args5);
        verifyValue("sin²+cos²", result5, 1.0, 0.005);
    }
    
    printSubHeader("asin, acos, atan");
    {
        Polynomial<Chebyshev> smallX = 0.3 * x;
        
        std::cout << "\nTesting asin(0.3*x) at x=0.5 (argument~0.15)" << std::endl;
        Polynomial<Chebyshev> asinP = asin(smallX);
        std::vector<double> args = {0.5, 0.0};
        double result = asinP.eval(args);
        double expected = std::asin(0.15);
        verifyValue("asin(0.15)", result, expected, 0.005);
        
        std::cout << "\nTesting acos(0.3*x) at x=0.5" << std::endl;
        Polynomial<Chebyshev> acosP = acos(smallX);
        result = acosP.eval(args);
        expected = std::acos(0.15);
        verifyValue("acos(0.15)", result, expected, 0.005);
        
        std::cout << "\nTesting atan(0.5*x) at x=0.6" << std::endl;
        Polynomial<Chebyshev> smallX2 = 0.5 * x;
        Polynomial<Chebyshev> atanP = atan(smallX2);
        args = {0.6, 0.0};
        result = atanP.eval(args);
        expected = std::atan(0.3);
        verifyValue("atan(0.3)", result, expected, 0.005);
    }
    
    printSubHeader("atan2");
    {
        Polynomial<Chebyshev> y_pos = 1.0 + 0.1*x;
        Polynomial<Chebyshev> x_pos = 1.0 + 0.05*x;
        Polynomial<Chebyshev> angle = atan2(y_pos, x_pos);
        
        std::cout << "\nTesting atan2(y, x) for quadrant I" << std::endl;
        std::vector<double> args = {0.0, 0.0};
        double result = angle.eval(args);
        double expected = std::atan2(1.0, 1.0);  // π/4
        verifyValue("atan2(1, 1)", result, expected, 1e-5);
        
        // Quadrant II
        Polynomial<Chebyshev> x_neg = -1.0 + 0.05*x;
        angle = atan2(y_pos, x_neg);
        result = angle.eval(args);
        expected = std::atan2(1.0, -1.0);  // 3π/4
        verifyValue("atan2(1, -1)", result, expected, 0.005);
    }
}

void testExponentialLogarithmic() {
    printHeader("Exponential and Logarithmic Functions");
    
    Polynomial<Chebyshev> x = DA(1);
    Polynomial<Chebyshev> y = DA(2);
    
    printSubHeader("exp and log");
    {
        // Test 1: exp on polynomial with T2 component
        std::cout << "\nTest 1: exp(0.3 + 0.1*T2(x)) at x=0.5" << std::endl;
        Polynomial<Chebyshev> T2 = 2.0*x*x - 1.0;
        Polynomial<Chebyshev> p1 = 0.3 + 0.1*T2;
        Polynomial<Chebyshev> expP1 = exp(p1);
        std::vector<double> args1 = {0.5, 0.0};
        double pVal1 = p1.eval(args1);
        double result1 = expP1.eval(args1);
        double expected1 = std::exp(pVal1);
        verifyValue("exp(P) with T2", result1, expected1, 0.005);
        
        // Test 2: exp on bivariate polynomial
        std::cout << "\nTest 2: exp(0.2*T1(x) + 0.15*T1(y)) at x=0.6, y=0.4" << std::endl;
        Polynomial<Chebyshev> p2 = 0.2*x + 0.15*y;
        Polynomial<Chebyshev> expP2 = exp(p2);
        std::vector<double> args2 = {0.6, 0.4};
        double pVal2 = p2.eval(args2);
        double result2 = expP2.eval(args2);
        double expected2 = std::exp(pVal2);
        verifyValue("exp(P) bivariate", result2, expected2, 0.005);
        
        // Test 3: log on polynomial with T3 component
        std::cout << "\nTest 3: log(1.5 + 0.2*T3(x)) at x=0.6" << std::endl;
        Polynomial<Chebyshev> T3 = 4.0*x*x*x - 3.0*x;
        Polynomial<Chebyshev> p3 = 1.5 + 0.2*T3;
        Polynomial<Chebyshev> logP3 = log(p3);
        std::vector<double> args3 = {0.6, 0.0};
        double pVal3 = p3.eval(args3);
        double result3 = logP3.eval(args3);
        double expected3 = std::log(pVal3);
        verifyValue("log(P) with T3", result3, expected3, 0.005);
        
        // Test 4: exp(log(P)) = P identity with higher-order terms
        std::cout << "\nTest 4: exp(log(P)) = P for P=1.2+0.15*T2(x)+0.1*T1(y)" << std::endl;
        Polynomial<Chebyshev> p4 = 1.2 + 0.15*T2 + 0.1*y;
        Polynomial<Chebyshev> logP4 = log(p4);
        Polynomial<Chebyshev> expLogP4 = exp(logP4);
        std::vector<double> args4 = {0.5, 0.3};
        double result4 = expLogP4.eval(args4);
        double expected4 = p4.eval(args4);
        verifyValue("exp(log(P))", result4, expected4, 0.005);
    }
    
    printSubHeader("log variants");
    {
        // Test 1: log10 with T2 component
        std::cout << "\nTest 1: log10(1.5 + 0.2*T2(x)) at x=0.5" << std::endl;
        Polynomial<Chebyshev> T2 = 2.0*x*x - 1.0;
        Polynomial<Chebyshev> p1 = 1.5 + 0.2*T2;
        Polynomial<Chebyshev> log10P1 = log10(p1);
        std::vector<double> args1 = {0.5, 0.0};
        double pVal1 = p1.eval(args1);
        double result1 = log10P1.eval(args1);
        double expected1 = std::log10(pVal1);
        verifyValue("log10(P) with T2", result1, expected1, 0.005);
        
        // Test 2: log2 with mixed terms
        std::cout << "\nTest 2: log2(1.3 + 0.1*T2(x) + 0.15*T1(y)) at x=0.6, y=0.4" << std::endl;
        Polynomial<Chebyshev> p2 = 1.3 + 0.1*T2 + 0.15*y;
        Polynomial<Chebyshev> log2P2 = log2(p2);
        std::vector<double> args2 = {0.6, 0.4};
        double pVal2 = p2.eval(args2);
        double result2 = log2P2.eval(args2);
        double expected2 = std::log2(pVal2);
        verifyValue("log2(P) mixed", result2, expected2, 0.005);
        
        // Test 3: logb with T3 component
        std::cout << "\nTest 3: logb(1.8 + 0.15*T3(x), 5) at x=0.7" << std::endl;
        Polynomial<Chebyshev> T3 = 4.0*x*x*x - 3.0*x;
        Polynomial<Chebyshev> p3 = 1.8 + 0.15*T3;
        Polynomial<Chebyshev> logbP3 = logb(p3, 5.0);
        std::vector<double> args3 = {0.7, 0.0};
        double pVal3 = p3.eval(args3);
        double result3 = logbP3.eval(args3);
        double expected3 = std::log(pVal3) / std::log(5.0);
        verifyValue("logb(P, 5) with T3", result3, expected3, 0.005);
    }
}

void testPowerRoot() {
    printHeader("Power and Root Functions");
    
    Polynomial<Chebyshev> x = DA(1);
    Polynomial<Chebyshev> y = DA(2);
    
    printSubHeader("sqrt and pow");
    {
        // Test 1: sqrt on polynomial with T2 component
        std::cout << "\nTest 1: sqrt(1.5 + 0.2*T2(x)) at x=0.5" << std::endl;
        Polynomial<Chebyshev> T2 = 2.0*x*x - 1.0;
        Polynomial<Chebyshev> p1 = 1.5 + 0.2*T2;
        Polynomial<Chebyshev> sqrtP1 = sqrt(p1);
        std::vector<double> args1 = {0.5, 0.0};
        double pVal1 = p1.eval(args1);
        double result1 = sqrtP1.eval(args1);
        double expected1 = std::sqrt(pVal1);
        verifyValue("sqrt(P) with T2", result1, expected1, 0.005);
        
        // Test 2: sqrt on bivariate polynomial
        std::cout << "\nTest 2: sqrt(1.3 + 0.15*T1(x) + 0.1*T1(y)) at x=0.6, y=0.4" << std::endl;
        Polynomial<Chebyshev> p2 = 1.3 + 0.15*x + 0.1*y;
        Polynomial<Chebyshev> sqrtP2 = sqrt(p2);
        std::vector<double> args2 = {0.6, 0.4};
        double pVal2 = p2.eval(args2);
        double result2 = sqrtP2.eval(args2);
        double expected2 = std::sqrt(pVal2);
        verifyValue("sqrt(P) bivariate", result2, expected2, 0.005);
        
        // Test 3: pow with T3 component
        std::cout << "\nTest 3: pow(1.4 + 0.15*T3(x), 2.5) at x=0.6" << std::endl;
        Polynomial<Chebyshev> T3 = 4.0*x*x*x - 3.0*x;
        Polynomial<Chebyshev> p3 = 1.4 + 0.15*T3;
        Polynomial<Chebyshev> powP3 = pow(p3, 2.5);
        std::vector<double> args3 = {0.6, 0.0};
        double pVal3 = p3.eval(args3);
        double result3 = powP3.eval(args3);
        double expected3 = std::pow(pVal3, 2.5);
        verifyValue("pow(P, 2.5) with T3", result3, expected3, 0.005);
        
        // Test 4: sqr with mixed terms
        std::cout << "\nTest 4: sqr(1.2 + 0.1*T2(x) + 0.15*T1(y))" << std::endl;
        Polynomial<Chebyshev> p4 = 1.2 + 0.1*T2 + 0.15*y;
        Polynomial<Chebyshev> sqrP4 = sqr(p4);
        std::vector<double> args4 = {0.5, 0.3};
        double pVal4 = p4.eval(args4);
        double result4 = sqrP4.eval(args4);
        double expected4 = pVal4 * pVal4;
        verifyValue("sqr(P) mixed", result4, expected4, 0.005);
        
        // Test 5: sqrt(sqr(P)) = P identity
        std::cout << "\nTest 5: sqrt(sqr(P)) = P for P=1.15+0.1*T2(x)" << std::endl;
        Polynomial<Chebyshev> p5 = 1.15 + 0.1*T2;
        Polynomial<Chebyshev> sqrP5 = sqr(p5);
        Polynomial<Chebyshev> roundTrip = sqrt(sqrP5);
        std::vector<double> args5 = {0.6, 0.0};
        double result5 = roundTrip.eval(args5);
        double expected5 = p5.eval(args5);
        verifyValue("sqrt(sqr(P))", result5, expected5, 0.005);
    }
    
    printSubHeader("Root functions");
    {
        Polynomial<Chebyshev> T2 = 2.0*x*x - 1.0;
        Polynomial<Chebyshev> T3 = 4.0*x*x*x - 3.0*x;
        
        // Test 1: cbrt with T2 component
        std::cout << "\nTest 1: cbrt(1.5 + 0.2*T2(x)) at x=0.5" << std::endl;
        Polynomial<Chebyshev> p1 = 1.5 + 0.2*T2;
        Polynomial<Chebyshev> cbrtP1 = cbrt(p1);
        std::vector<double> args1 = {0.5, 0.0};
        double pVal1 = p1.eval(args1);
        double result1 = cbrtP1.eval(args1);
        double expected1 = std::cbrt(pVal1);
        verifyValue("cbrt(P) with T2", result1, expected1, 0.005);
        
        // Test 2: root (5th root) with T3 component
        std::cout << "\nTest 2: root(1.6 + 0.15*T3(x), 5) at x=0.6" << std::endl;
        Polynomial<Chebyshev> p2 = 1.6 + 0.15*T3;
        Polynomial<Chebyshev> root5P2 = root(p2, 5.0);
        std::vector<double> args2 = {0.6, 0.0};
        double pVal2 = p2.eval(args2);
        double result2 = root5P2.eval(args2);
        double expected2 = std::pow(pVal2, 1.0/5.0);
        verifyValue("root5(P) with T3", result2, expected2, 0.005);
        
        // Test 3: isrt with bivariate polynomial
        std::cout << "\nTest 3: isrt(1.4 + 0.1*T1(x) + 0.15*T1(y)) at x=0.6, y=0.4" << std::endl;
        Polynomial<Chebyshev> p3 = 1.4 + 0.1*x + 0.15*y;
        Polynomial<Chebyshev> isrtP3 = isrt(p3);
        std::vector<double> args3 = {0.6, 0.4};
        double pVal3 = p3.eval(args3);
        double result3 = isrtP3.eval(args3);
        double expected3 = 1.0 / std::sqrt(pVal3);
        verifyValue("isrt(P) bivariate", result3, expected3, 0.005);
        
        // Test 4: icrt with mixed terms
        std::cout << "\nTest 4: icrt(1.7 + 0.12*T2(x) + 0.1*T1(y)) at x=0.5, y=0.3" << std::endl;
        Polynomial<Chebyshev> p4 = 1.7 + 0.12*T2 + 0.1*y;
        Polynomial<Chebyshev> icrtP4 = icrt(p4);
        std::vector<double> args4 = {0.5, 0.3};
        double pVal4 = p4.eval(args4);
        double result4 = icrtP4.eval(args4);
        double expected4 = 1.0 / std::cbrt(pVal4);
        verifyValue("icrt(P) mixed", result4, expected4, 0.005);
    }
}

void testHyperbolic() {
    printHeader("Hyperbolic Functions");
    
    Polynomial<Chebyshev> x = DA(1);
    Polynomial<Chebyshev> y = DA(2);
    
    printSubHeader("sinh, cosh, tanh");
    {
        Polynomial<Chebyshev> T2 = 2.0*x*x - 1.0;
        Polynomial<Chebyshev> T3 = 4.0*x*x*x - 3.0*x;
        
        // Test 1: sinh with T2 component
        std::cout << "\nTest 1: sinh(0.3 + 0.1*T2(x)) at x=0.5" << std::endl;
        Polynomial<Chebyshev> p1 = 0.3 + 0.1*T2;
        Polynomial<Chebyshev> sinhP1 = sinh(p1);
        std::vector<double> args1 = {0.5, 0.0};
        double pVal1 = p1.eval(args1);
        double result1 = sinhP1.eval(args1);
        double expected1 = std::sinh(pVal1);
        verifyValue("sinh(P) with T2", result1, expected1, 0.005);
        
        // Test 2: cosh with T3 component
        std::cout << "\nTest 2: cosh(0.25 + 0.15*T3(x)) at x=0.6" << std::endl;
        Polynomial<Chebyshev> p2 = 0.25 + 0.15*T3;
        Polynomial<Chebyshev> coshP2 = cosh(p2);
        std::vector<double> args2 = {0.6, 0.0};
        double pVal2 = p2.eval(args2);
        double result2 = coshP2.eval(args2);
        double expected2 = std::cosh(pVal2);
        verifyValue("cosh(P) with T3", result2, expected2, 0.005);
        
        // Test 3: tanh with bivariate polynomial
        std::cout << "\nTest 3: tanh(0.2*T1(x) + 0.15*T1(y)) at x=0.5, y=0.4" << std::endl;
        Polynomial<Chebyshev> p3 = 0.2*x + 0.15*y;
        Polynomial<Chebyshev> tanhP3 = tanh(p3);
        std::vector<double> args3 = {0.5, 0.4};
        double pVal3 = p3.eval(args3);
        double result3 = tanhP3.eval(args3);
        double expected3 = std::tanh(pVal3);
        verifyValue("tanh(P) bivariate", result3, expected3, 0.005);
        
        // Test 4: cosh² - sinh² = 1 identity
        std::cout << "\nTest 4: cosh²(P) - sinh²(P) = 1 for P with T2" << std::endl;
        Polynomial<Chebyshev> p4 = 0.3 + 0.12*T2;
        Polynomial<Chebyshev> sinhP4 = sinh(p4);
        Polynomial<Chebyshev> coshP4 = cosh(p4);
        Polynomial<Chebyshev> cosh2 = coshP4 * coshP4;
        Polynomial<Chebyshev> sinh2 = sinhP4 * sinhP4;
        Polynomial<Chebyshev> diff = cosh2 - sinh2;
        std::vector<double> args4 = {0.5, 0.0};
        double result4 = diff.eval(args4);
        verifyValue("cosh²-sinh²", result4, 1.0, 0.005);
    }
    
    printSubHeader("asinh, acosh, atanh");
    {
        Polynomial<Chebyshev> T2 = 2.0*x*x - 1.0;
        Polynomial<Chebyshev> T3 = 4.0*x*x*x - 3.0*x;
        
        // Test 1: asinh with T2 component
        std::cout << "\nTest 1: asinh(0.2 + 0.08*T2(x)) at x=0.5" << std::endl;
        Polynomial<Chebyshev> p1 = 0.2 + 0.08*T2;
        Polynomial<Chebyshev> asinhP1 = asinh(p1);
        std::vector<double> args1 = {0.5, 0.0};
        double pVal1 = p1.eval(args1);
        double result1 = asinhP1.eval(args1);
        double expected1 = std::asinh(pVal1);
        verifyValue("asinh(P) with T2", result1, expected1, 0.005);
        
        // Test 2: asinh with bivariate polynomial
        std::cout << "\nTest 2: asinh(0.15*T1(x) + 0.1*T1(y)) at x=0.6, y=0.4" << std::endl;
        Polynomial<Chebyshev> p2 = 0.15*x + 0.1*y;
        Polynomial<Chebyshev> asinhP2 = asinh(p2);
        std::vector<double> args2 = {0.6, 0.4};
        double pVal2 = p2.eval(args2);
        double result2 = asinhP2.eval(args2);
        double expected2 = std::asinh(pVal2);
        verifyValue("asinh(P) bivariate", result2, expected2, 0.005);
        
        // Test 3: acosh with T2 component
        std::cout << "\nTest 3: acosh(1.2 + 0.15*T2(x)) at x=0.5" << std::endl;
        Polynomial<Chebyshev> p3 = 1.2 + 0.15*T2;
        Polynomial<Chebyshev> acoshP3 = acosh(p3);
        std::vector<double> args3 = {0.5, 0.0};
        double pVal3 = p3.eval(args3);
        double result3 = acoshP3.eval(args3);
        double expected3 = std::acosh(pVal3);
        verifyValue("acosh(P) with T2", result3, expected3, 0.005);
        
        // Test 4: atanh with T3 component
        std::cout << "\nTest 4: atanh(0.15 + 0.08*T3(x)) at x=0.6" << std::endl;
        Polynomial<Chebyshev> p4 = 0.15 + 0.08*T3;
        Polynomial<Chebyshev> atanhP4 = atanh(p4);
        std::vector<double> args4 = {0.6, 0.0};
        double pVal4 = p4.eval(args4);
        double result4 = atanhP4.eval(args4);
        double expected4 = std::atanh(pVal4);
        verifyValue("atanh(P) with T3", result4, expected4, 0.005);
    }
}

void testUtilityFunctions() {
    printHeader("Utility Functions");
    
    Polynomial<Chebyshev> x = DA(1);
    Polynomial<Chebyshev> y = DA(2);
    
    // Future utility function tests can be added here
}

void testNormAnalysis() {
    printHeader("Norm and Analysis Functions");
    
    Polynomial<Chebyshev> x = DA(1);
    Polynomial<Chebyshev> T0 = 1.0;
    Polynomial<Chebyshev> T1 = x;
    Polynomial<Chebyshev> T2 = 2.0*x*T1 - T0;
    
    printSubHeader("Norm functions");
    {
        Polynomial<Chebyshev> p = 1.0 + 2.0*T1 + 3.0*T2;
        
        std::cout << "\nTest polynomial: 1 + 2*T1 + 3*T2" << std::endl;
        
        // abs() = max coefficient
        double absVal = abs(p);
        std::cout << "  abs(p): " << absVal << std::endl;
        verifyValue("abs(p)", absVal, 3.0);
        
        // norm(0) = max norm
        double norm0 = norm(p, 0);
        verifyValue("norm(p, 0)", norm0, 3.0);
        
        // norm(1) = sum of absolute values
        double norm1 = norm(p, 1);
        verifyValue("norm(p, 1)", norm1, 6.0);
        
        // norm(2) = Euclidean norm
        double norm2 = norm(p, 2);
        double expected = std::sqrt(1.0 + 4.0 + 9.0);
        verifyValue("norm(p, 2)", norm2, expected, 1e-10);
        
        // Member function style should give same results
        std::cout << "\nVerifying member function equivalence:" << std::endl;
        verifyValue("p.abs() vs abs(p)", p.abs(), absVal);
        verifyValue("p.norm(0) vs norm(p,0)", p.norm(0), norm0);
        verifyValue("p.norm(1) vs norm(p,1)", p.norm(1), norm1);
    }
    
    printSubHeader("orderNorm");
    {
        Polynomial<Chebyshev> p = 1.0 + 2.0*T1 + 3.0*T2;
        std::vector<double> ordNorm = orderNorm(p, 0, 0);
        
        std::cout << "\norderNorm(p, 0, 0) by total order:" << std::endl;
        std::cout << "  Order 0: " << ordNorm[0] << " (expected: 1.0)" << std::endl;
        verifyValue("orderNorm[0]", ordNorm[0], 1.0);
        std::cout << "  Order 1: " << ordNorm[1] << " (expected: 2.0)" << std::endl;
        verifyValue("orderNorm[1]", ordNorm[1], 2.0);
        std::cout << "  Order 2: " << ordNorm[2] << " (expected: 3.0)" << std::endl;
        verifyValue("orderNorm[2]", ordNorm[2], 3.0);
    }
    
    printSubHeader("bound");
    {
        Polynomial<Chebyshev> p = 1.0 + 0.5*T1;
        DACE::Interval bounds = bound(p);
        
        std::cout << "\nbound(1 + 0.5*T1) over [-1,1]" << std::endl;
        std::cout << "  Lower bound: " << bounds.m_lb << " (expected: 0.5)" << std::endl;
        std::cout << "  Upper bound: " << bounds.m_ub << " (expected: 1.5)" << std::endl;
        // T1 ranges from -1 to 1, so 1+0.5*T1 ranges from 0.5 to 1.5
        verifyValue("Lower bound", bounds.m_lb, 0.5, 0.1);
        verifyValue("Upper bound", bounds.m_ub, 1.5, 0.1);
    }
    
}

void testEvaluation() {
    printHeader("Evaluation Functions");
    
    Polynomial<Chebyshev> x = DA(1);
    Polynomial<Chebyshev> y = DA(2);
    
    printSubHeader("Numerical evaluation");
    {
        Polynomial<Chebyshev> p = 1.0 + 2.0*x + 3.0*y;
        
        std::cout << "\nEvaluating 1 + 2x + 3y at (0.5, 0.3)" << std::endl;
        std::vector<double> args = {0.5, 0.3};
        double result = p.eval(args);
        double expected = 1.0 + 2.0*0.5 + 3.0*0.3;
        verifyValue("eval", result, expected);
    }
    
    printSubHeader("Symbolic composition");
    {
        Polynomial<Chebyshev> T1_x = x;
        Polynomial<Chebyshev> p = T1_x * T1_x;  // x²
        
        Polynomial<Chebyshev> Q = 0.5 * x;  // Substitute 0.5x for x
        DACE::AlgebraicVector<Polynomial<Chebyshev>> args(2);
        args[0] = Q;
        args[1] = y;
        
        Polynomial<Chebyshev> composition = p.eval(args);
        
        std::cout << "\nP(Q) where P=x² and Q=0.5x" << std::endl;
        std::cout << "Expected: (0.5x)² = 0.25*x²" << std::endl;
        
        // x² in Chebyshev: T1² = 0.5*(T2 + T0)
        // (0.5x)² = 0.25*x² = 0.25*0.5*(T2 + T0) = 0.125*(T2 + T0)
        verifyCoefficient(composition, {0,0}, 0.125);
        verifyCoefficient(composition, {2,0}, 0.125);
    }
}

void testCompoundAssignment() {
    printHeader("Compound Assignment Operators");
    
    Polynomial<Chebyshev> x = DA(1);
    
    printSubHeader("+=, -=, *=, /=");
    {
        Polynomial<Chebyshev> p = 2.0 + x;
        
        std::cout << "\nTesting +=" << std::endl;
        p += 1.0;
        verifyCoefficient(p, {0,0}, 3.0);
        verifyCoefficient(p, {1,0}, 1.0);
        
        std::cout << "\nTesting -=" << std::endl;
        p -= x;
        verifyCoefficient(p, {0,0}, 3.0);
        verifyCoefficient(p, {1,0}, 0.0, 1e-15);
        
        std::cout << "\nTesting *=" << std::endl;
        p = 2.0 + x;
        p *= 2.0;
        verifyCoefficient(p, {0,0}, 4.0);
        verifyCoefficient(p, {1,0}, 2.0);
        
        std::cout << "\nTesting /=" << std::endl;
        p /= 2.0;
        verifyCoefficient(p, {0,0}, 2.0);
        verifyCoefficient(p, {1,0}, 1.0);
    }
}

void testMultivariateOperations() {
    printHeader("Multivariate Operations");
    
    Polynomial<Chebyshev> x = DA(1);
    Polynomial<Chebyshev> y = DA(2);
    
    printSubHeader("Mixed partial derivatives");
    {
        Polynomial<Chebyshev> p = x*y;  // T1(x)*T1(y)
        
        std::cout << "\nTesting ∂²/∂x∂y [T1(x)*T1(y)]" << std::endl;
        Polynomial<Chebyshev> dp_dx = p.deriv(1);
        Polynomial<Chebyshev> dp_dxdy = dp_dx.deriv(2);
        
        // d/dx[T1(x)*T1(y)] = T0(x)*T1(y)
        // d/dy[T0(x)*T1(y)] = T0(x)*T0(y) = 1
        verifyCoefficient(dp_dxdy, {0,0}, 1.0);
    }
    
    printSubHeader("Cross-term multiplication");
    {
        Polynomial<Chebyshev> px = 1.0 + x;
        Polynomial<Chebyshev> py = 1.0 + y;
        Polynomial<Chebyshev> prod = px * py;
        
        std::cout << "\nTesting (1+x)*(1+y)" << std::endl;
        verifyCoefficient(prod, {0,0}, 1.0);
        verifyCoefficient(prod, {1,0}, 1.0);
        verifyCoefficient(prod, {0,1}, 1.0);
        verifyCoefficient(prod, {1,1}, 1.0);
    }
}

// =============================================================
//  Main Execution
// =============================================================

int main(int argc, char** argv) {
    std::cout << "\n" << CYAN << "=========================================================" << RESET << std::endl;
    std::cout << CYAN << "  COMPREHENSIVE CHEBYSHEV POLYNOMIAL ALGEBRA TEST SUITE " << RESET << std::endl;
    std::cout << CYAN << "=========================================================" << RESET << std::endl;
    
    // Initialize DACE with order 10, 2 variables (higher order for accurate transcendentals)
    DA::init(10, 2);
    
    const std::map<std::string, std::function<void()>> tests = {
        {"basic-arithmetic", testBasicArithmetic},
        {"calculus", testCalculus},
        {"trigonometric", testTrigonometric},
        {"exponential-logarithmic", testExponentialLogarithmic},
        {"power-root", testPowerRoot},
        {"hyperbolic", testHyperbolic},
        {"utility", testUtilityFunctions},
        {"norm-analysis", testNormAnalysis},
        {"evaluation", testEvaluation},
        {"compound-assignment", testCompoundAssignment},
        {"multivariate", testMultivariateOperations}
    };

    try {
        if (argc == 2) {
            std::string selected = argv[1];
            auto it = tests.find(selected);
            if (it == tests.end()) {
                std::cerr << RED << "Unknown test category: " << selected << RESET << "\n";
                std::cerr << "Available categories:\n";
                for (const auto& kv : tests) {
                    std::cerr << "  - " << kv.first << "\n";
                }
                return 2;
            }
            it->second();
        } else {
            // Run all test categories
            for (const auto& kv : tests) {
                kv.second();
            }
        }
        
        // Final summary
        std::cout << "\n" << CYAN << "=========================================================" << RESET << std::endl;
        if (g_failures == 0) {
            std::cout << GREEN << "           ALL TEST CATEGORIES COMPLETED                 " << RESET << std::endl;
        } else {
            std::cout << RED << "           TEST FAILURES: " << g_failures << "                           " << RESET << std::endl;
        }
        std::cout << CYAN << "=========================================================" << RESET << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "\n" << RED << "ERROR: " << e.what() << RESET << std::endl;
        return 1;
    }
    
    return (g_failures == 0) ? 0 : 1;
}

