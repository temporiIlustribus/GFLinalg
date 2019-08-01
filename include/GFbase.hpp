#pragma once
#include <array>
#include <stdexcept>
#include <iostream>
#include <memory>

namespace GFlinalg {

#ifdef  USE_PRIM_POL_TABLE
    
    static uint64_t prim_poly[33] =
    {   /*  0 */    0U,
        /*  1 */    3U,
        /*  2 */    7U,
        /*  3 */    11U,
        /*  4 */    19U,
        /*  5 */    37U,
        /*  6 */    67U,
        /*  7 */    137U,
        /*  8 */    285U,
        /*  9 */    529U,
        /* 10 */    1033U,
        /* 11 */    2053U,
        /* 12 */    4179U,
        /* 13 */    8219U,
        /* 14 */    17475U,
        /* 15 */    32771‬U,
        /* 16 */    69643U,
        /* 17 */    131081‬U,
        /* 18 */    262273‬U,
        /* 19 */    524327‬U,
        /* 20 */    1048585‬U,
        /* 21 */    2097157‬U,
        /* 22 */    4194307‬U,
        /* 23 */    8388641‬U,
        /* 24 */    16777351‬U,
        /* 25 */    33554441‬U,
        /* 26 */    67108935U,
        /* 27 */    134217767U,
        /* 28 */    268435465U‬,
        /* 29 */    536870917U‬,
        /* 30 */    1082130439U,
        /* 31 */    2147483657U‬,
        /* 32 */    4299161607‬U
    };
#endif

    template <class T, T modPol>
    class  BasicBinPolynomial;
    template <class T, T modPol>
    class  PowBinPolynomial;
    template <class T, T modPol>
    class  TableBinPolynomial;
    template <class T>
    class BasicGFElem;
    template <class T>
    class PowGFElem;
    template <class T>
    class TableGFElem;
    namespace op {
        //! Get the position of the leading 1 in the polynomial (from the left)
        template<class T>
        uint8_t leadElemPos(const T& pol, uint8_t startPos = 1, size_t order = sizeof(T) << 3) {
            uint8_t pos = 0;
            for (uint8_t i = startPos; i < order + 1; ++i) {
                pos = order - i;
                if ((pol >> pos) & 1) break;
            }
            return order - pos;
        }
        //! Get the degree of the modulus polynomial
        template<class T>
        constexpr uint8_t modPolDegree(const T& modPol) {
            uint8_t pos = 0;
            uint8_t containerSize = sizeof(T) << 3;
            for (uint8_t i = 1; i < containerSize + 1; ++i) {
                pos = containerSize - i;
                if ((modPol >> pos) & 1) break;
            }
            return pos;
        }
        //! Internal polynomial addition (is equivalent to XOR)
        template<class Polynomial>
        Polynomial polSum(const Polynomial& a, const Polynomial& b) {
            Polynomial res(a.val() ^ b.val());
            return res;
        }
        //! Internal multiplication version 1
        template<class Polynomial>
        Polynomial polMulOld(const Polynomial& a, const Polynomial& b) {
            Polynomial res(0);
            for (size_t i = 0; i < a.gfOrder(); ++i) {
                if ((b.value >> i) & 1) {
                    res.val() ^= a.val() << i;
                }
            }
            res.reduce();
            return res;
        }
        //! Internal multiplication version 2
        template<class Polynomial>
        Polynomial polMul(const Polynomial &a, const Polynomial &b) {
            Polynomial res{a};
            res.val() = 0;
            auto av = a.val();
            auto bv = b.val();
            while (bv > 0) {
                if (bv & 1)
                    res.val() ^= av;
                bv >>= 1;
                av <<= 1;
                if (av & a.gfOrder())
                    av ^= a.modpol;
            }
            return res;
        }
        //! Binary power function in Galois field, performing operations in polynomial form
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
        //! Internal division
        template<class Polynomial>
        Polynomial polDiv(const Polynomial& a, const Polynomial& b) {
            // Get b^-1: a / b = a * b^-1
            if (b.val() == 0)
                throw std::out_of_range("Division by zero");
            Polynomial invB = op::pow<Polynomial>(b, b.gfOrder() - 2);
            return polMul<Polynomial>(a, invB);
        }
        //! Writes Galois Field element to std::ostream in polynomial form
        /*!
         * * \param out ostream that the elemnt will be written to
         * * \param pol GF element to write to the ostream
         *
         *  GF element will be written to the ostream in polynomial form with higher powers of x
         *  written first.
         *
         *  Example:
         *
         *  Let a be a GF element with coefficients {1,0,1,1};
         *
         *  out << a; - Results in x^3 + x + 1 being written to "out";
         */
        template <class Polynomial>
        std::ostream& operator << (std::ostream& out, Polynomial pol) {
            size_t deg = pol.degree();
            bool flag = false;
            while (pol.val()) {
                if (deg > 0) {
                    out.operator<<('x');
                    if (deg > 1) {
                        out.operator<<('^');
                        out.operator<<(deg);
                    }
                } else {
                    out.operator<<('1');
                }
                flag = true;
                pol.val() ^= (1 << pol.degree());
                deg = pol.degree(deg - 1);
                if (pol.val())
                    out.operator<<('+');
            }
            if (!flag)
                out.operator<<('0');
            return out;
        }

        //! Pair of arrays, used as LUTs
        template <class T, T modPol>
        struct  LUTArrPair {
            static constexpr uint64_t order = (1U << modPolDegree<T>(modPol));
            std::array<T, (order - 1) << 1> indToPol; /*!<Look-up table, converts powers of primitive element to polynomials */
            std::array<size_t, order> polToInd; /*!<Look-up table, converts polynomials to powers of primitive element*/
            /*!
             * Creates a pair of look-up arrays (LUTPair):
             * * indToPol: power of primitive element -> polynomial
             * * polToInd: polynomial -> power of primitive element
             */
            LUTArrPair() : polToInd(), indToPol() {
                BasicBinPolynomial<T, modPol> counter{1};
                BasicBinPolynomial<T, modPol> modifier{2};
                for (size_t i = 0; i < order - 1; ++i) {
                    indToPol[i] = counter.val();
                    polToInd[indToPol[i]] = i;
                    counter *= modifier;
                }
                // This is to avoid % operations in math operators
                for (size_t i = order - 1; i < indToPol.size(); ++i)
                    indToPol[i] = indToPol[i - order + 1];
            }
            LUTArrPair(const std::array<T, (order - 1) * 2>& alph, const std::array<size_t, order>& ind) :
                indToPol(alph), polToInd(ind) {
            }
            LUTArrPair(const std::pair<std::array<T, (order - 1) * 2>, std::array<size_t, order>>& val) :
                indToPol(val.first), polToInd(val.second) {
            }
        };
        //! Pair of vectors, used as LUTs
        template <class T>
        struct LUTVectPair {
            std::vector<T> indToPol;
            std::vector<T> polToInd;
            size_t order;
            /*!
             * Creates a pair of look-up vectors (LUTPair):
             * * indToPol: power of primitive element -> polynomial
             * * polToInd: polynomial -> power of primitive element
             */
            LUTVectPair(const T& modPol) : order(1U << modPolDegree<T>(modPol)) {
                polToInd.resize(order);
                indToPol.resize((order - 1) << 1);
                BasicGFElem<T> counter{1, modPol};
                BasicGFElem<T> modifier{2, modPol};
                for (size_t i = 0; i < order - 1; ++i) {
                    indToPol[i] = counter.val();
                    polToInd[indToPol[i]] = i;
                    counter *= modifier;
                }
                // This is to avoid % operations in math operators
                for (size_t i = order - 1; i < indToPol.size(); ++i) {
                    indToPol[i] = indToPol[i - order + 1];
                }
            }
            LUTVectPair(const LUTVectPair&& lut) : polToInd(std::move(lut.polToInd)),
                indToPol(std::move(lut.indToPol)) {
            }
        };
    }
    using op::pow;
    using op::leadElemPos;
    using op::operator<<;
    using op::LUTArrPair;
    using op::LUTVectPair;
}