#pragma once
#include <iostream>
#include <memory>
#include <vector>

namespace GFlinalg {

    template <class T>
    class  binPolyniomial {
    private:
        T value;
        uint8_t sz;

        // Internal multiplication
        binPolyniomial<T> polMul(const binPolyniomial<T>& a, const binPolyniomial<T>& b) const {
            binPolyniomial<T> res(0);
            {
                T mask = 1;
                for (size_t i = 0; i < b.sz; ++i) {
                    res.value ^= a.value * (b.value & (mask));
                    mask <<= 1;
                }
            }
            res.reduce();
            return res;
        }
        // Internal addition
        binPolyniomial<T> polSum(const binPolyniomial<T>& a, const binPolyniomial<T>& b) const {
            binPolyniomial<T> res(a.value ^ b.value);
            res.reduce();
            return res;
        }

        uint8_t leadElemPos(T& pol) {
            uint8_t pos = 0;
            // Adjust modulus polynomial - zero leftmost bit
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
        binPolyniomial<T>() : value(0), sz(1) {}
        binPolyniomial<T>(const T& val) : value(val), sz(sizeof(val) << 3) {}
        //binPolyniomial<T>(const binPolynomial<T>& pol) : value(pol.value), sz(pol.sz) {}
        binPolyniomial<T>(int val) : value(static_cast<T>(val)), sz(sizeof(T) << 3) {}

        T getVal() { return value; }
        T& val() { return value; }
        size_t size() { return sz; }
        
        void reduce() {
            reduce(modPol);
        }

        operator int() { return static_cast<int>(value); }
        operator uint64_t() { return static_cast<uint64_t>(value); }

        binPolyniomial<T> operator * (binPolyniomial<T>& pol) const {
            return polMul(*this, pol);
        }

        binPolyniomial<T> operator + (binPolyniomial<T>& pol) const {
            return polSum(*this, pol);
        }
        binPolyniomial<T> operator *= (binPolyniomial<T>& pol) {
            *this = *this * pol;
            return *this;
        }

        binPolyniomial<T> operator += (binPolyniomial<T>& pol) {
            *this = *this + pol;
            return *this;
        }
    };


    //////////////////
    //              //
    //  TODO  LIST  //
    //              //
    //////////////////

    // Note: binary operators for the polynomials currently compute everything at runtime.
    // Final version will use tables to do all of the operations in O(1) time (currently O(n^2) for GF(2^n))
    // The binPolynomial class will contain all 3 tables as static fields

    // Inversion Table
    template<size_t N, class Data, class T>
    class inverseTable {
    public:
        static constexpr T& data = inverseTable<N - 1, binPolynomialInverse(N - 1), Data>::data;
    };

    template<class Data, class T>
    class inverseTable<0, Data, T> {
    public:
        static constexpr T data[] = {Data};
    };


    // Addition table
    template<size_t N, class Data, class T>
    class sumTable {
    public:
        static constexpr T data[];
    };

    // Multiplication table
    template<size_t N, class Data, class T>
    class mulTable {
    public:
        static constexpr T data[];
    };


}