#include <benchmark/benchmark.h>
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include "GFlinalg.hpp"
typedef GFlinalg::BasicBinPolynomial<uint8_t, 11> basicPol8;
typedef GFlinalg::PowBinPolynomial<uint8_t, 11> powPol8;
typedef GFlinalg::TableBinPolynomial<uint8_t, 11> tablePol8;

typedef GFlinalg::BasicBinPolynomial<uint32_t, 37> basicPol32;
typedef GFlinalg::PowBinPolynomial<uint32_t, 37> powPol32;
typedef GFlinalg::TableBinPolynomial<uint32_t, 37> tablePol32;
// Comment this if you are having problems building project
template<>
powPol8::ArrayPair powPol8::alphaToIndex = powPol8::makeAlphaToIndex();
template<>
tablePol8::GFtable tablePol8::mulTable = tablePol8::makeMulTable();
template<>
tablePol8::GFtable tablePol8::divTable = tablePol8::makeInvMulTable();

template <class Pol>
static void BM_Reduction(benchmark::State& state) {
    Pol testVal(0);
    for (auto _ : state) {
        state.PauseTiming();
        size_t value = rand();
        testVal.val() = value;
        state.ResumeTiming();
        testVal.reduce();
    }
}
BENCHMARK_TEMPLATE(BM_Reduction, basicPol8);
BENCHMARK_TEMPLATE(BM_Reduction, powPol8);
BENCHMARK_TEMPLATE(BM_Reduction, tablePol8);
BENCHMARK_TEMPLATE(BM_Reduction, basicPol32);
BENCHMARK_TEMPLATE(BM_Reduction, powPol32);
BENCHMARK_TEMPLATE(BM_Reduction, tablePol32);

template <class Pol>
static void BM_Addition(benchmark::State& state) {
    Pol temp(0);
    for (auto _ : state) {
        state.PauseTiming();
        Pol a(rand());
        Pol b(rand());
        state.ResumeTiming();
        temp = a + b;
    }
}
BENCHMARK_TEMPLATE(BM_Addition, basicPol8);
BENCHMARK_TEMPLATE(BM_Addition, powPol8);
BENCHMARK_TEMPLATE(BM_Addition, tablePol8);
BENCHMARK_TEMPLATE(BM_Addition, basicPol32);
BENCHMARK_TEMPLATE(BM_Addition, powPol32);
BENCHMARK_TEMPLATE(BM_Addition, tablePol32);

template <class Pol>
static void BM_Mul(benchmark::State& state) {
    Pol temp(0);
    for (auto _ : state) {
        state.PauseTiming();
        Pol a(rand());
        Pol b(rand());
        state.ResumeTiming();
        temp = a * b;
    }
}
BENCHMARK_TEMPLATE(BM_Mul, basicPol8);
BENCHMARK_TEMPLATE(BM_Mul, powPol8);
BENCHMARK_TEMPLATE(BM_Mul, tablePol8);
BENCHMARK_TEMPLATE(BM_Mul, basicPol32);
BENCHMARK_TEMPLATE(BM_Mul, powPol32);
BENCHMARK_TEMPLATE(BM_Mul, tablePol32);

template <class Pol>
static void BM_Div(benchmark::State& state) {
    Pol temp(0);
    for (auto _ : state) {
        state.PauseTiming();
        Pol a(rand());
        Pol b(rand());
        state.ResumeTiming();
        temp = a / b;
    }
}

BENCHMARK_TEMPLATE(BM_Div, basicPol8);
BENCHMARK_TEMPLATE(BM_Div, powPol8);
BENCHMARK_TEMPLATE(BM_Div, tablePol8);
BENCHMARK_TEMPLATE(BM_Div, basicPol32);
BENCHMARK_TEMPLATE(BM_Div, powPol32);
BENCHMARK_TEMPLATE(BM_Div, tablePol32);

BENCHMARK_MAIN();