#pragma once
#include <iostream>
#include <memory>
#include <vector>

namespace GFlinalg {
    template <class T, size_t SZ>
    class  baseBinPolynomial {
    public:
        constexpr static size_t order = 1 << SZ;
    protected:
        using GFtable = std::array<std::array <T, order>, order>;
        T value;
        uint8_t sz;

        // Internal multiplication version 1
        baseBinPolynomial polMul(const baseBinPolynomial& a, const baseBinPolynomial& b) const {
            baseBinPolynomial<T,SZ> res(0);
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
        // TODO:
        // Internal multiplication version 2
        //baseBinPolynomial polMulAlt(const baseBinPolynomial& a, const baseBinPolynomial& b) const {
        //    baseBinPolynomial res(0);
        //    auto x = a.value;
        //    auto y = b.value;
        //    uint8_t modLeadPos = sz - leadElemPos(modPol);
        //    do {
        //        res ^= (y & 1);
        //
        //    } while (y << = 1)
        //            return res;
        //}

        // Internal addition
        baseBinPolynomial polSum(const baseBinPolynomial& a, const baseBinPolynomial& b) const {
            baseBinPolynomial res(a.value ^ b.value);
            return res;
        }

        // Get the position of the leading 1 in the polynomial
        uint8_t leadElemPos(const T& pol) const {
            uint8_t pos = 0;
            for (uint8_t i = 1; i < (sizeof(T) << 3); ++i) {
                if (pol >> ((sizeof(T) << 3) - i) & 1) {
                    pos = i;
                    break;
                }
            }
            return pos;
        }

        void reduce(T& modulus) {
            auto pos = leadElemPos(modulus);
            // Reduce by modulus
            uint8_t i = 1;
            while (value >= modulus) {
                if ((value >> (sz - i)) & 1)
                    value ^= modulus << (pos - i);
                ++i;
            }
        }

    public:
        static T modPol;
        baseBinPolynomial() : value(0), sz(sizeof(T) << 3) {}
        baseBinPolynomial(const T& val) : value(val), sz(sizeof(val) << 3) {
            reduce();
        }
        baseBinPolynomial(int val) : value(static_cast<T>(val)), sz(sizeof(T) << 3) {}

        T getVal() { return value; }
        T& val() { return value; }
        size_t size() { return sz; }
        size_t GFsize() { return SZ; }

        void reduce() {
            reduce(modPol);
        }

        operator int() { return static_cast<int>(value); }
        operator uint64_t() { return static_cast<uint64_t>(value); }
        operator T() { return value; }


    };


    //////////////////
    //              //
    //  TODO  LIST  //
    //              //
    //////////////////

    // Note: binary operators for the binPolynomial compute everything runtime
    // TableBased version will use tables to do all of the operations in O(1) time 
    // (where binPolynomial takes O(n^2) for GF(2^n))
    // The TableBased_binPolynomial class will contain 3 tables

    template <class T, size_t SZ>
    class BasicBinPolynomial : public baseBinPolynomial<T, SZ> {
    public:
        using baseBinPolynomial<T,SZ>::baseBinPolynomial;
        BasicBinPolynomial(const baseBinPolynomial<T, SZ>& pol) : baseBinPolynomial<T,SZ>(pol) {}
        BasicBinPolynomial operator * (const BasicBinPolynomial& other) const {
            return this->polMul(*this, other);
        }
        BasicBinPolynomial operator + (const BasicBinPolynomial& other) const {
            return this->polSum(*this, other);
        }
        BasicBinPolynomial& operator *= (const BasicBinPolynomial& other) {
            *this = this->polSum(*this, other);
            return *this;
        }
        BasicBinPolynomial& operator += (const BasicBinPolynomial& other) {
            *this = this->polSum(*this, other);
            return *this;
        }
    };

    template <class T, size_t SZ>
    class TableBinPolynomial : public BasicBinPolynomial<T, SZ> {
    private:
        static T mulTable[];

    public:
        using baseBinPolynomial<T, SZ>::order;
        using baseBinPolynomial<T, SZ>::GFtable;
        TableBinPolynomial(const BasicBinPolynomial<T, SZ>& pol) : BasicBinPolynomial<T, SZ>(pol) {}
        static constexpr GFtable makeTable() {
            GFtable temp;
            // Note: requires checking whethere we have traversed through all elements
            for (size_t i = 0; i < order; ++i) {
                for (size_t j = 0; j < order; ++j) {
                    BasicBinPolynomial<T, SZ> a(i);
                    BasicBinPolynomial<T, SZ> b(j);
                    temp[i][j] = (a * b).getValue();
                }
            }
            mulTable = temp;
            return mulTable;
        }
        using BasicBinPolynomial<T,SZ>::BasicBinPolynomial;
        
        TableBinPolynomial operator * (const TableBinPolynomial& other) const {
            return mulTable[*this->getVal(), other.getVal()];
        }

        TableBinPolynomial& operator *= (const TableBinPolynomial& other) const {
            *this->val = *this * other;
            return *this;
        }

        TableBinPolynomial operator + (const TableBinPolynomial& other) const {
            return TableBinPolynomial(other + *this);
        }

        TableBinPolynomial operator += (const TableBinPolynomial& other) const {
            *this = *this + other;
            return *this;
        }
    };

}