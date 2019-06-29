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
    public:
        operator int() { return static_cast<int>(value); }
        operator uint64_t() { return static_cast<uint64_t>(value); }
        operator uint32_t() { return static_cast<uint32_t>(value); }
        T getVal() { return value; }
        T& val() { return value; }
        T size() { return sz; }
        binPolyniomial(T& val): value(val), sz(sizeof(val) << 3) {}
        binPolyniomial(int val) : value(static_cast<T>(val)), sz(sizeof(val) << 3) {}

    }

    template<class T>
    T binPolynomialMul(T& a, T& b, T& modulus) {
        T mask = 1;
        T res = 0;
        for (size_t i = 0; i < b.size; ++i) {
            res ^= a << b & (mask)
            mask <<= 1;
        }
        std::cout << res;
        return res;
    }

    template<class T>
    T binPolynomialSum(T& a, T& b, T& modulus) {

    }

    template<class T>
    T binPolynomialInverse(T& a, T& b, T& modulus) {

    }

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
        static constexpr T& data = sumTable<N - 1, binPolynomialSum(N - 1), Data>::data;
    };

    template<class Data, class T>
    class sumTable<0, Data, T> {
    public:
        static constexpr T data[] = {Data};
    };


    //Таблица умножения
    template<size_t N, class Data, class T>
    class sumTable {
    public:
        static constexpr T& data = sumTable<N - 1, binPolynomialMul(N - 1), Data>::data;
    };

    template<class Data, class T>
    class sumTable<0, Data, T> {
    public:
        static constexpr T data[] = {Data};
    };

}