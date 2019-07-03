#pragma once
#include <iostream>
#include <memory>
#include <array>
#include <vector>

namespace GFlinalg {
    template <class T, size_t SZ, T modPol>
    class  BaseBinPolynomial {
    public:
        constexpr static size_t order = 1 << SZ;
    protected:
        using vectPair = std::pair<std::vector<T>, std::vector<size_t>>;
        using GFtable = std::array<std::array <T, order>, order>;
        T value;
        uint8_t sz;

        // Internal multiplication version 1
        BaseBinPolynomial polMul(const BaseBinPolynomial& a, const BaseBinPolynomial& b) const {
            BaseBinPolynomial<T, SZ> res(0);
            uint8_t leadPos = sz - leadElemPos(a.value);
            uint8_t modLeadPos = sz - leadElemPos(modPol);
            for (size_t i = 0; i < b.sz; ++i) {
                if ((b.value >> i) & 1) {
                    res.value ^= a.value << i;
                }
            }
            res.reduce();
            return res;
        }

        // Internal addition
        BaseBinPolynomial polSum(const BaseBinPolynomial& a, const BaseBinPolynomial& b) const {
            BaseBinPolynomial res(a.value ^ b.value);
            return res;
        }

        // Internal division
        BaseBinPolynomial polDiv(const BaseBinPolynomial& a, const BaseBinPolynomial& b) const {
            auto invB = galoisPow(b, SZ-2);
            return polMul(a, invB);
        }

        // Get the position of the leading 1 in the polynomial
        uint8_t leadElemPos(const T& pol) const {
            uint8_t pos = 0;
            for (uint8_t i = 1; i < sz; ++i) {
                if (pol >> (sz - i) & 1) {
                    pos = i;
                    break;
                }
            }
            return pos;
        }

        void reduce(const T& modulus) {
            auto pos = leadElemPos(modulus);
            // Reduce by modulus
            uint8_t i = 1;
            while (value >= modulus) {
                if ((value >> (sz - i)) & 1)
                    value ^= modulus << (pos - i);
                ++i;
            }
            if (value >= (1 << (sz - pos)))
                value ^= modulus;
        }

        BaseBinPolynomial galoisPow(const BaseBinPolynomial val, size_t power) {
            BaseBinPolynomial res = 1;
            while (power) {
                if (power & 1) {
                    res *= val;
                } else {
                    val *= val;
                    power >>= 1;
                }
            }
            return res;
        }

    public:

        BaseBinPolynomial() : value(0), sz(sizeof(T) << 3) {}
        BaseBinPolynomial(const T& val) : value(val), sz(sizeof(val) << 3) {
            reduce();
        }
        BaseBinPolynomial(int val) : value(static_cast<T>(val)), sz(sizeof(T) << 3) {}

        T getVal() { return value; }
        T& val() { return value; }
        size_t size() { return sz; }
        size_t GFsize() { return SZ; }

        void reduce() {
            reduce(modPol);
        }
        template <class T1>
        operator T1() { return static_cast<T1>(value); }



    };

    template <class T, size_t SZ, T modPol>
    class BasicBinPolynomial : public BaseBinPolynomial<T, SZ, modPol> {
    public:
        using BaseBinPolynomial<T, SZ, modPol>::BaseBinPolynomial;
        BasicBinPolynomial(const BaseBinPolynomial<T, SZ, modPol>& pol) : BaseBinPolynomial<T, SZ, modPol>(pol) {}
        BasicBinPolynomial operator * (const BasicBinPolynomial& other) const {
            return this->polMul(*this, other);
        }

        BasicBinPolynomial& operator *= (const BasicBinPolynomial& other) {
            *this = this->polSum(*this, other);
            return *this;
        }

        BasicBinPolynomial operator + (const BasicBinPolynomial& other) const {
            return this->polSum(*this, other);
        }
        
        BasicBinPolynomial& operator += (const BasicBinPolynomial& other) {
            *this = this->polSum(*this, other);
            return *this;
        }
    };

    template <class T, size_t SZ, T modPol>
    class PowBinPolynomial : public BaseBinPolynomial<T, SZ, modPol> {
    public:
        using BaseBinPolynomial<T, SZ, modPol>::value;
        using BaseBinPolynomial<T, SZ, modPol>::order;
        using BaseBinPolynomial<T, SZ, modPol>::polSum;
        using BaseBinPolynomial<T, SZ, modPol>::vectPair;
    private:
        static vectPair alphaToIndex;
    public:

        PowBinPolynomial(const BasicBinPolynomial<T, SZ, modPol>& pol) : BaseBinPolynomial<T, SZ, modPol>(pol.getVal()) {}
        PowBinPolynomial(const BaseBinPolynomial<T, SZ, modPol>& pol) : BaseBinPolynomial<T, SZ, modPol>(pol) {}
        /*
        Creates a pair of vectors (vectPair):
            first: power of primitive element -> polynomial
            second: polynomial -> power of primitive element
        */
        static constexpr auto makeAlphaToIndex() {
            std::vector<T> alpha(order - 1);
            std::vector<size_t> index(order);
            T counter = 1;
            for (size_t i = 0; i < order - 1; ++i) {
                alpha[i] = BasicBinPolynomial<T, SZ, modPol>(counter).getVal();
                index[alpha[i]] = i;
                counter <<= 1;
            }
            return std::make_pair(alpha, index);
        }

        // Debug function
        static void printAlpha() {
            for (size_t i = 0; i < order - 1; ++i) {
                std::cout << static_cast<size_t>(alphaToIndex.first[i]);
                std::cout << " : " << static_cast<size_t>(alphaToIndex.second[alphaToIndex.first[i]]) << "\n";
            }
        }

        /*
        Returns a pair of vectors (vectPair):
            first: power of primitive element -> polynomial
            second: polynomial -> power of primitive element
        */
        static vectPair getAlphaToIndex() {
            return alphaToIndex;
        }

        PowBinPolynomial operator * (const PowBinPolynomial& other) const {
            if (this->value == 0 || other.value == 0)
                return PowBinPolynomial(0);
            return PowBinPolynomial(alphaToIndex.first[(alphaToIndex.second[this->value] +
                                                        alphaToIndex.second[other.value]) % (order - 1)]);
        }

        PowBinPolynomial& operator *= (const PowBinPolynomial& other) {
            *this->val = alphaToIndex.first[(alphaToIndex.second[*this->value] +
                                             alphaToIndex.second[other.value]) % (order - 1)];
            return *this;
        }

        PowBinPolynomial operator / (const PowBinPolynomial& other) const {
            return PowBinPolynomial(alphaToIndex.first[alphaToIndex.second[this->value] -
                                    alphaToIndex.second[other.value] % (order - 1)]);
        }

        PowBinPolynomial& operator /= (const PowBinPolynomial& other) {
            *this->val = alphaToIndex.first[alphaToIndex.second[this->value] -
                alphaToIndex.second[other.value] % (order - 1)];
            return *this;
        }

        PowBinPolynomial operator + (const PowBinPolynomial& other) const {
            return PowBinPolynomial(polSum(this->value, other.value));
        }

        PowBinPolynomial operator += (const PowBinPolynomial& other) {
            *this = *this + other;
            return *this;
        }
    };


    template <class T, size_t SZ, T modPol>
    class TableBinPolynomial : public BasicBinPolynomial<T, SZ, modPol> {
    private:
        static T mulTable[];
        static T divTable[];
    public:
        using BaseBinPolynomial<T, SZ, modPol>::order;
        using BaseBinPolynomial<T, SZ, modPol>::GFtable;

        TableBinPolynomial(const BasicBinPolynomial<T, SZ, modPol>& pol) : BasicBinPolynomial<T, SZ, modPol>(pol) {}
        TableBinPolynomial(const PowBinPolynomial<T, SZ, modPol>& pol) : BasicBinPolynomial<T, SZ, modPol>(pol.getVal()) {}

        /*
        Creates multiplication table (GFtable):
            Table[pol1][pol2] = pol1 * pol2;
        */
        static constexpr GFtable makeMulTable() {
            GFtable temp;
            // Note: requires checking whether we have traversed through all elements
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = 0; j < order; ++j) {
                    BasicBinPolynomial<T, SZ, modPol> a(i);
                    BasicBinPolynomial<T, SZ, modPol> b(j);
                    temp[i][j] = (a * b).getValue();
                }
            }
            return temp;
        }

        /*
        Creates division table (GFtable) :
            Table[pol1][pol2] = pol1 / pol2;

        Note: multiplication table is required to create this table
        */
        static constexpr GFtable makeDivTable() {
            GFtable temp;
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = 0; j < order; ++j) {
                    BasicBinPolynomial<T, SZ, modPol> a(i);
                    BasicBinPolynomial<T, SZ, modPol> b(j);
                    temp[mulTable[a.value][b.value]][a.value] = b.value;
                    temp[mulTable[a.value][b.value]][b.value] = a.value;
                }
            }
            return temp;
        }

        using BasicBinPolynomial<T, SZ, modPol>::BasicBinPolynomial;

        TableBinPolynomial operator * (const TableBinPolynomial& other) const {
            return mulTable[*this->getVal(), other.getVal()];
        }

        TableBinPolynomial& operator *= (const TableBinPolynomial& other) const {
            *this->val = *this * other;
            return *this;
        }

        TableBinPolynomial operator / (const TableBinPolynomial& other) const {
            return divTable[*this->getVal(), other.getVal()];
        }

        TableBinPolynomial& operator /= (const TableBinPolynomial& other) const {
            *this->val = *this / other;
            return *this;
        }

        TableBinPolynomial operator + (const TableBinPolynomial& other) const {
            return TableBinPolynomial(*this + other);
        }

        TableBinPolynomial operator += (const TableBinPolynomial& other) const {
            *this = *this + other;
            return *this;
        }
    };
}