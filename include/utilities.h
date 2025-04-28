//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "macros.h"
#include "mesh.h"

class Empty {};

template <typename T>
concept HasDmdt = requires(T t) {
    { t.dm_dt(0.0, 0.0, 0.0) } -> std::same_as<Real>;
};

template <typename T>
concept HasDedt = requires(T t) {
    { t.deps_k_dt(0.0, 0.0, 0.0) } -> std::same_as<Real>;
};

template <typename T>
concept HasSigma = requires(T t) {
    { t.sigma0(0.0, 0.0) } -> std::same_as<Real>;
};

template <typename T>
concept HasU = requires(T t) {
    { t.u } -> std::same_as<Real>;
};

/*template <typename T>
concept HasMass = requires(T t) {
    { t.mass(0.0, 0.0, 0.0) } -> std::same_as<Real>;
};*/

#define MAKE_THIS_ODEINT_STATE(data, array_size)                             \
    using array_type = std::array<Real, array_size>;                         \
    using value_type = typename array_type::value_type;                      \
    using iterator = typename array_type::iterator;                          \
    using const_iterator = typename array_type::const_iterator;              \
                                                                             \
    constexpr size_t size() const noexcept { return array_size; }            \
    constexpr iterator begin() noexcept { return data.begin(); }             \
    constexpr iterator end() noexcept { return data.end(); }                 \
    constexpr const_iterator begin() const noexcept { return data.begin(); } \
    constexpr const_iterator end() const noexcept { return data.end(); }     \
    constexpr value_type& operator[](size_t i) noexcept { return data[i]; }  \
    constexpr const value_type& operator[](size_t i) const noexcept { return data[i]; }

void print_array(Array const& arr);

/********************************************************************************************************************
 * TEMPLATE FUNCTION: print (Variadic)
 * DESCRIPTION: Prints multiple arguments to standard output, separated by a space.
 *              This template function recursively prints all provided arguments.
 ********************************************************************************************************************/
inline void print() {  // Base case: terminates the recursion.
    std::cout << std::endl;
}

template <typename First, typename... Rest>
void print(First first, Rest... rest) {
    std::cout << first << " ";
    print(rest...);  // Recursive call to print the remaining arguments.
}

/********************************************************************************************************************
 * FUNCTION TYPE DEFINITIONS
 * DESCRIPTION: Defines convenient aliases for unary, binary, and ternary functions operating on Reals.
 *              These function types are used throughout the codebase for various mathematical operations
 *              and physical calculations.
 ********************************************************************************************************************/
using UnaryFunc = std::function<Real(Real)>;                // Function taking one Real argument
using BinaryFunc = std::function<Real(Real, Real)>;         // Function taking two Real arguments
using TernaryFunc = std::function<Real(Real, Real, Real)>;  // Function taking three Real arguments

/********************************************************************************************************************
 * NAMESPACE: func
 * DESCRIPTION: Contains inline constexpr lambda functions that return constant values.
 ********************************************************************************************************************/
namespace func {
    // Always returns 0 regardless of the input.
    inline constexpr auto zero_3d = [](Real /*phi*/, Real /*theta*/, Real /*t*/) constexpr noexcept { return 0.; };
    inline constexpr auto zero_2d = [](Real /*phi*/, Real /*theta*/) constexpr noexcept { return 0.; };
    // Always returns 1 regardless of the input.
    inline constexpr auto one_3d = [](Real /*phi*/, Real /*theta*/, Real /*t*/) constexpr noexcept { return 1.; };
    inline constexpr auto one_2d = [](Real /*phi*/, Real /*theta*/) constexpr noexcept { return 1.; };
}  // namespace func

/********************************************************************************************************************
 * FUNCTION: Basic Math Functions                                                                                   *
 * DESCRIPTION: Inline functions for specific power calculations, a step function, and unit conversion.             *
 ********************************************************************************************************************/
inline Real pow52(Real a) { return std::sqrt(a * a * a * a * a); }
inline Real pow43(Real a) { return std::cbrt(a * a * a * a); }
inline Real pow23(Real a) { return std::cbrt(a * a); }
inline Real stepFunc(Real x) { return x > 0 ? 1 : 0; }
inline constexpr Real eVtoHz(Real eV) { return eV / con::h; }

/********************************************************************************************************************
 * FUNCTION: Fast Math & Interpolation Prototypes                                                                   *
 * DESCRIPTION: Prototypes for fast power, logarithm, and various interpolation functions.                          *
 ********************************************************************************************************************/
Real fast_pow(Real a, Real b);
Real fast_log(Real x);
Real interp(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
Real eq_space_interp(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
Real loglog_interp(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);
Real eq_space_loglog_interp(Real x0, Array const& x, Array const& y, bool lo_extrap = false, bool hi_extrap = false);

/********************************************************************************************************************
 * FUNCTION: Root Finding (Bisection Method)                                                                        *
 * DESCRIPTION: Template function to find the root of a function using the bisection method.                        *
 ********************************************************************************************************************/
template <typename Fun>
auto root_bisect(Fun f, decltype(f(0)) low, decltype(f(0)) high, decltype(f(0)) eps = 1e-6) -> decltype(f(0)) {
    using Scalar = decltype(f(0));
    for (; (high - low) > std::fabs((high + low) * 0.5) * eps;) {
        Scalar mid = 0.5 * (high + low);
        if (f(mid) * f(high) > 0)
            high = mid;
        else
            low = mid;
    }
    return 0.5 * (high + low);
}

/********************************************************************************************************************
 * FUNCTION: Utility Templates                                                                                      *
 * DESCRIPTION: Template functions for computing the minimum and maximum of provided values.                        *
 ********************************************************************************************************************/
template <typename T>
T min(T value) {
    return value;
}

template <typename T, typename... Args>
T min(T first, Args... args) {
    return std::min(first, std::min(args...));
}

template <typename T>
T max(T value) {
    return value;
}

template <typename T, typename... Args>
T max(T first, Args... args) {
    return std::max(first, std::max(args...));
}

/********************************************************************************************************************
 * FUNCTION: Fast Exponential and Logarithm Functions                                                               *
 * DESCRIPTION: Inline functions that provide fast approximations of exponential and logarithm functions using      *
 *              alternative methods when EXTREME_SPEED is defined.                                                  *
 ********************************************************************************************************************/
inline Real fast_exp(Real x) {
#ifdef EXTREME_SPEED
    // if (std::isnan(x)) return std::numeric_limits<Real>::quiet_NaN();
    // if (x == std::numeric_limits<Real>::infinity()) return std::numeric_limits<Real>::infinity();
    // if (x == -std::numeric_limits<Real>::infinity()) return 0.0;

    constexpr Real ln2 = 0.6931471805599453;
    constexpr Real inv_ln2 = 1.4426950408889634;

    Real y = x * inv_ln2;
    int64_t k = static_cast<int64_t>(y + (y >= 0 ? 0.5 : -0.5));
    Real r = x - k * ln2;

    // Real p = 1.0 + r * (1.0 + r * (0.5 + r * (0.166666666666666 + r * 0.041666666666666664)));

    Real p = 1.0 + r * (1.0 + r * (0.5 + r * (0.166666666666666)));

    return std::ldexp(p, k);
#else
    return std::exp(x);
#endif
}

inline double fast_log(double x) {
#ifdef EXTREME_SPEED
    if (x <= 0.) return -std::numeric_limits<double>::infinity();
    if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
    if (x == std::numeric_limits<double>::infinity()) return std::numeric_limits<double>::infinity();

    uint64_t bits;
    std::memcpy(&bits, &x, sizeof(x));
    int64_t exponent = ((bits >> 52) & 0x7FF) - 1023;
    bits = (bits & 0x000FFFFFFFFFFFFFULL) | 0x3FF0000000000000ULL;
    double f;
    std::memcpy(&f, &bits, sizeof(f));
    double p = -1.49278 + (2.11263 + (-0.729104 + 0.10969 * f) * f) * f;
    return p + 0.6931471806 * exponent;
#else
    return std::log(x);
#endif
}

inline double fast_log2(double val) {
#ifdef EXTREME_SPEED
    int64_t* const exp_ptr = reinterpret_cast<int64_t*>(&val);
    int64_t x = *exp_ptr;
    int log2 = ((x >> 52) & 0x7FF) - 1023;  // extract exponent bits

    // Step 2: Normalize mantissa to [1.0, 2.0)
    x &= ~(0x7FFLL << 52);  // clear exponent bits
    x |= (1023LL << 52);    // set exponent to 0
    *exp_ptr = x;           // val is now normalized to [1.0, 2.0)
    double mantissa = val;

    // Step 3: Polynomial approximation of log2(mantissa) in [1, 2)
    double y = mantissa - 1.0;  // small value in [0, 1)
    double log2_mantissa =
        y * (1.44269504088896340736 + y * (-0.721347520444482371076 + y * (0.479381953382630073738)));

    /*double log2_mantissa =
        y * (1.44269504088896340736 +  // 1/ln(2)
        y * (-0.721347520444482371076 +
        y * (0.479381953382630073738 +
        y * (-0.360673760222241834679 +
        y * (0.288539008177138356808 +
        y * (-0.139304958445395653244))))));*/

    return log2_mantissa + log2;
#else
    return std::log2(val);
#endif
}

inline double fast_exp2(double x) {
#ifdef EXTREME_SPEED
    int int_part = (int)x;
    double frac_part = x - int_part;

    // Polynomial approximation for 2^frac_part where 0 <= frac_part < 1
    // 4th order polynomial gives good balance of speed and accuracy
    /*double poly = 1.0 + frac_part * (0.693147180559945 +
                                     frac_part * (0.240226506959101 +
                                                  frac_part * (0.0555041086648216 + frac_part *
       0.00961812910762848)));*/
    double poly =
        1.0 + frac_part * (0.693147180559945 + frac_part * (0.240226506959101 + frac_part * (0.0555041086648216)));

    // Combine with integer power of 2 using bit manipulation
    int64_t bits = ((int64_t)(int_part + 1023)) << 52;
    double factor = *reinterpret_cast<double*>(&bits);

    return (poly * factor);
#else
    return std::exp2(x);
#endif
}

inline Real fast_pow(Real a, Real b) { return fast_exp2(b * fast_log2(a)); }
