#pragma once
#include <array>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace GFlinalg {

template <class T, T modPol>
class BasicBinPolynomial;

template <class T, T modPol>
class PowBinPolynomial;

template <class T, T modPol>
class TableBinPolynomial;

template <class T>
class BasicGFElem;

template <class T>
class PowGFElem;

template <class T>
class TableGFElem;

namespace op {

/**
 * @return Position of the leading 1 in the polynomial.
 */
template <class T>
uint8_t leadElemPos(const T& pol, uint8_t startPos = 1, size_t order = sizeof(T) << 3) {
    uint8_t pos = 0;

    for (uint8_t i = startPos; i < order + 1; ++i) {
        pos = order - i;

        if ((pol >> pos) & 1)
            break;
    }

    return order - pos;
}

/**
 * @return The degree of modulus polynomial.
 */
template <class T>
constexpr uint8_t modPolDegree(const T& modPol) {
    uint8_t pos           = 0;
    uint8_t containerSize = sizeof(T) << 3;

    for (uint8_t i = 1; i < containerSize + 1; ++i) {
        pos = containerSize - i;

        if ((modPol >> pos) & 1)
            break;
    }

    return pos;
}

/**
 * Internal polynomial addition (equivalent to \c XOR).
 */
template <class Polynomial>
Polynomial polSum(const Polynomial& a, const Polynomial& b) {
    return Polynomial(a.val() ^ b.val());
}

/// Internal multiplication version 1
template <class Polynomial>
Polynomial polMulOld(const Polynomial& a, const Polynomial& b) {
    Polynomial res(0);

    for (size_t i = 0; i < a.gfOrder(); ++i)
        if ((b.value >> i) & 1)
            res.val() ^= a.val() << i;

    res.reduce();

    return res;
}

/// Internal multiplication version 2
template <class Polynomial>
Polynomial polMul(const Polynomial& a, const Polynomial& b) {
    Polynomial res(a);

    res.val() = 0;

    auto av   = a.val();
    auto bv   = b.val();

    while (bv > 0) {
        if (bv & 1)
            res.val() ^= av;

        bv >>= 1;
        av <<= 1;

        if (av & a.gfOrder())
            av ^= a.getMod();
    }

    return res;
}

/**
 * Binary power function. Operation is performed in the polynomial form.
 */
template <class Polynomial>
Polynomial pow(Polynomial val, size_t power) {
    Polynomial res(val);
    res.val() = 1;

    while (power) {
        if (power & 1) {
            res = res * val;
            --power;
        } else {
            val *= val;
            power >>= 1;
        }
    }

    return res;
}

template <class Polynomial>
Polynomial polDiv(const Polynomial& a, const Polynomial& b) {
    if (b.val() == 0)
        throw std::out_of_range("Division by zero");

    return polMul<Polynomial>(a, /*b ^(-1)*/ op::pow<Polynomial>(b, b.gfOrder() - 2));
}

/**
 * GF element will be written to the \c ostream in polynomial form with higher powers of \c x written first.
 *
 * @example <tt>GF{1,0,1,1} -> x^3 + x + 1</tt>
 */
template <class Polynomial>
std::ostream& operator<<(std::ostream& out, Polynomial pol) {
    size_t deg = pol.degree();

    bool empty = true;

    while (pol.val()) {
        if (deg > 0) {
            out.operator<<("x");

            if (deg > 1)
                out.operator<<('^') << deg;
        } else
            out.operator<<('1');

        empty = false;

        pol.val() ^= 1 << pol.degree();

        deg = pol.degree(deg - 1);

        if (pol.val())
            out.operator<<('+');
    }

    if (empty)
        out.operator<<('0');

    return out;
}

template <class T, T modPol>
struct LUTArrPair {
    static constexpr uint64_t order = (1U << modPolDegree<T>(modPol));

    /**
     * Lookup table, converts powers of the primitive element to polynomials.
     */
    std::array<T, (order - 1) << 1> indToPol;

    /**
     * Lookup table, converts polynomials to powers of the primitive element.
     */
    std::array<size_t, order> polToInd;

    LUTArrPair() : polToInd(), indToPol() {
        BasicBinPolynomial<T, modPol> counter{1};
        BasicBinPolynomial<T, modPol> modifier{2};

        for (size_t i = 0; i < order - 1; ++i) {
            indToPol[i]           = counter.val();
            polToInd[indToPol[i]] = i;

            counter *= modifier;
        }

        // This is to avoid % operations in math operators
        for (size_t i = order - 1; i < indToPol.size(); ++i)
            indToPol[i] = indToPol[i - order + 1];
    }

    LUTArrPair(const std::array<T, (order - 1) * 2>& alph, const std::array<size_t, order>& ind)
        : indToPol(alph), polToInd(ind) {}

    LUTArrPair(const std::pair<std::array<T, (order - 1) * 2>, std::array<size_t, order>>& val)
        : indToPol(val.first), polToInd(val.second) {}
};

template <class T>
struct LUTVectPair {
    std::vector<T> indToPol;
    std::vector<T> polToInd;

    size_t order;

    explicit LUTVectPair(const T& modPol) : order(1U << modPolDegree<T>(modPol)) {
        polToInd.resize(order);
        indToPol.resize((order - 1) << 1);

        BasicGFElem<T> counter{1, modPol};
        BasicGFElem<T> modifier{2, modPol};

        for (size_t i = 0; i < order - 1; ++i) {
            indToPol[i]           = counter.val();
            polToInd[indToPol[i]] = i;
            counter *= modifier;
        }

        // This is to avoid % operations in math operators
        for (size_t i = order - 1; i < indToPol.size(); ++i)
            indToPol[i] = indToPol[i - order + 1];
    }
    
    explicit LUTVectPair(const LUTVectPair&& lut)
        : polToInd(std::move(lut.polToInd)), indToPol(std::move(lut.indToPol)) {}
};
} // namespace op

using op::leadElemPos;
using op::pow;
using op::operator<<;
using op::LUTArrPair;
using op::LUTVectPair;
}