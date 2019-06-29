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

        // Internal multiplication procedure
        binPolyniomial<T> polMul(binPolyniomial<T>& a, binPolyniomial<T>& b,
                                 binPolyniomial<T>& modulus) {
            binPolyniomial<T> res(0);
            {
                binPolyniomial<T> mask(1);
                for (size_t i = 0; i < b.size(); ++i) {
                    res.value ^= a.value * (b.value & (mask.value));
                    mask.value <<= 1;
                }
            }

            res.reduce();

            return res;
        }

        binPolyniomial<T> polSum(binPolyniomial<T>& a, binPolyniomial<T>& b,
                                 binPolyniomial<T>& modulus) {


        }

        uint8_t leadElemPos(binPolyniomial<T> pol) {
            uint8_t pos = 0;
            // Adjust modulus polynomial - zero leftmost bit
            for (uint8_t i = 1; i < pol.sz; ++i) {
                if (pol.value >> (pol.sz - i) & 1) {
                    pos = i;
                    break;
                }
            }
            return pos;
        }

    public:
        static T modPol;

        binPolyniomial(T& val) : value(val), sz(sizeof(val) << 3) {}
        binPolyniomial() : value(0), sz(1) {}
        binPolyniomial(int val) : value(static_cast<T>(val)), sz(sizeof(T) << 3) {}

        T getVal() { return value; }
        T& val() { return value; }
        size_t size() { return sz; }
        static void setModulus(T& val) {
            modPol = val;
        }
        static void setModulus(int val) {
            modPol = static_cast<T>(val);
        }
        void reduce(binPolyniomial<T>& modulus) {
            auto pos = leadElemPos(modulus);
            // Reduce by modulus
            uint8_t i = 1;
            while (value >= modulus.value) {
                if ((value >> (sz-i)) & 1)
                    value ^= modulus.value << (pos - i);
                ++i;
            }
            //this->value ^= temp.value;
        }

        void reduce() {
            binPolyniomial<T> temp(modPol);
            reduce(temp);
        }

        operator int() { return static_cast<int>(value); }
        operator uint64_t() { return static_cast<uint64_t>(value); }
        binPolyniomial<T> operator * (binPolyniomial<T>& pol) {
            binPolyniomial<T> temp(modPol);
            return polMul(this, pol, temp);
        }

        binPolyniomial<T> operator + (binPolyniomial<T>& pol) {
            binPolyniomial<T> temp(modPol);
            return polSum(this, pol, temp);
        }
    
    
    };


    //////////////////
    //              //
    //  TODO  LIST  //
    //              //
    //////////////////


    // Таблица инверсий
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


    // Таблица сложения
    template<size_t N, class Data, class T>
    class sumTable {
    public:
        static constexpr T& data;
    };

    template<class Data, class T>
    class sumTable<0, Data, T> {
    public:
        static constexpr T data[] = {Data};
    };


    // Таблица умножения
    template<size_t N, class Data, class T>
    class mulTable {
    public:
        static constexpr T& data;
    };

    template<class Data, class T>
    class mulTable<0, Data, T> {
    public:
        static constexpr T data[] = {Data};
    };

}