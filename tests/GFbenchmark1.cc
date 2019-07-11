#include <benchmark/benchmark.h>
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include "GFlinalg.hpp"
typedef GFlinalg::BasicBinPolynomial<uint8_t, 11> basicPol8;
typedef GFlinalg::PowBinPolynomial<uint8_t, 11> powPol8;
typedef GFlinalg::TableBinPolynomial<uint8_t, 11> tablePol8;

typedef GFlinalg::BasicBinPolynomial<uint16_t, 19> basicPol16;
typedef GFlinalg::PowBinPolynomial<uint16_t, 19> powPol16;
typedef GFlinalg::TableBinPolynomial<uint16_t, 19> tablePol16;


typedef GFlinalg::BasicBinPolynomial<uint32_t, 37> basicPol32;
typedef GFlinalg::PowBinPolynomial<uint32_t, 37> powPol32;
typedef GFlinalg::TableBinPolynomial<uint32_t, 37> tablePol32;
// Comment this if you are having problems building project
template<>
const powPol8::LUTPair powPol8::alphaToIndex = powPol8::makeAlphaToIndex();
template<>
const tablePol8::GFtable tablePol8::mulTable = tablePol8::makeMulTable();
template<>
const tablePol8::GFtable tablePol8::divTable = tablePol8::makeInvMulTable();

template<>
const powPol16::LUTPair powPol16::alphaToIndex = powPol16::makeAlphaToIndex();
template<>                                     
const tablePol16::GFtable tablePol16::mulTable = tablePol16::makeMulTable();
template<>                                     
const tablePol16::GFtable tablePol16::divTable = tablePol16::makeInvMulTable();

template<>
const powPol32::LUTPair powPol32::alphaToIndex = powPol32::makeAlphaToIndex();
template<>                                     
const tablePol32::GFtable tablePol32::mulTable = tablePol32::makeMulTable();
template<>                                     
const tablePol32::GFtable tablePol32::divTable = tablePol32::makeInvMulTable();

template <class Pol>
static void BM_Reduction(benchmark::State& state) {
    Pol testVal(0);
    std::uniform_int_distribution<uint32_t> uid(0, (1 << Pol::gfDegree()) - 1);
    std::default_random_engine rd;
    for (auto _ : state) {
        testVal.val() = uid(rd);
        benchmark::DoNotOptimize(
            testVal.reduce()
        );
    }
}
BENCHMARK_TEMPLATE(BM_Reduction, basicPol8);
BENCHMARK_TEMPLATE(BM_Reduction, powPol8);
BENCHMARK_TEMPLATE(BM_Reduction, tablePol8);
BENCHMARK_TEMPLATE(BM_Reduction, basicPol16);
BENCHMARK_TEMPLATE(BM_Reduction, powPol16);
BENCHMARK_TEMPLATE(BM_Reduction, tablePol16);
BENCHMARK_TEMPLATE(BM_Reduction, basicPol32);
BENCHMARK_TEMPLATE(BM_Reduction, powPol32);
BENCHMARK_TEMPLATE(BM_Reduction, tablePol32);

template <class Pol>
static void BM_Addition(benchmark::State& state) {
    Pol temp(0);
    std::uniform_int_distribution<uint32_t> uid(0, (1 << Pol::gfDegree()) - 1);
    std::default_random_engine rd;
    for (auto _ : state) {
        Pol a(uid(rd));
        Pol b(uid(rd));
        benchmark::DoNotOptimize(
            temp = a + b
        );
    }
}
BENCHMARK_TEMPLATE(BM_Addition, basicPol8);
BENCHMARK_TEMPLATE(BM_Addition, powPol8);
BENCHMARK_TEMPLATE(BM_Addition, tablePol8);
BENCHMARK_TEMPLATE(BM_Addition, basicPol16);
BENCHMARK_TEMPLATE(BM_Addition, powPol16);
BENCHMARK_TEMPLATE(BM_Addition, tablePol16);
BENCHMARK_TEMPLATE(BM_Addition, basicPol32);
BENCHMARK_TEMPLATE(BM_Addition, powPol32);
BENCHMARK_TEMPLATE(BM_Addition, tablePol32);

template <class Pol>
static void BM_Mul(benchmark::State& state) {
    Pol temp(0);
    std::uniform_int_distribution<uint32_t> uid(0, (1 << Pol::gfDegree()) - 1);
    std::default_random_engine rd;
    for (auto _ : state) {
        Pol a(uid(rd));
        Pol b(uid(rd));
        benchmark::DoNotOptimize(
            temp = a * b
        );
    }
}

template <class Pol>
static void BM_MulAlt(benchmark::State& state) {
    Pol temp(0);
    std::uniform_int_distribution<uint32_t> uid(1, Pol::gfOrder() - 1);
    std::default_random_engine rd;
    for (auto _ : state) {
        Pol a(uid(rd));
        Pol b(uid(rd));
        auto temp(a * b);
        while (temp != a) {
            temp *= b;
        }
        benchmark::DoNotOptimize(temp);
    }
}
BENCHMARK_TEMPLATE(BM_Mul, basicPol8);
BENCHMARK_TEMPLATE(BM_Mul, powPol8);
BENCHMARK_TEMPLATE(BM_Mul, tablePol8);
BENCHMARK_TEMPLATE(BM_Mul, basicPol16);
BENCHMARK_TEMPLATE(BM_Mul, powPol16);
BENCHMARK_TEMPLATE(BM_Mul, tablePol16);
BENCHMARK_TEMPLATE(BM_Mul, basicPol32);
BENCHMARK_TEMPLATE(BM_Mul, powPol32);
BENCHMARK_TEMPLATE(BM_Mul, tablePol32);

BENCHMARK_TEMPLATE(BM_MulAlt, basicPol32);
BENCHMARK_TEMPLATE(BM_MulAlt, powPol32);
BENCHMARK_TEMPLATE(BM_MulAlt, tablePol32);

template <class Pol>
static void BM_Div(benchmark::State& state) {
    Pol temp(0);
    std::uniform_int_distribution<uint32_t> uid1(1, (1 << Pol::gfDegree()) - 1);
    std::uniform_int_distribution<uint32_t> uid2(0, (1 << Pol::gfDegree()) - 1);
    std::default_random_engine rd;
    for (auto _ : state) {
        Pol a(uid2(rd));
        Pol b(uid1(rd));
        benchmark::DoNotOptimize(
            temp = a / b
        );
    }
}

BENCHMARK_TEMPLATE(BM_Div, basicPol8);
BENCHMARK_TEMPLATE(BM_Div, powPol8);
BENCHMARK_TEMPLATE(BM_Div, tablePol8);
BENCHMARK_TEMPLATE(BM_Div, basicPol16);
BENCHMARK_TEMPLATE(BM_Div, powPol16);
BENCHMARK_TEMPLATE(BM_Div, tablePol16);
BENCHMARK_TEMPLATE(BM_Div, basicPol32);
BENCHMARK_TEMPLATE(BM_Div, powPol32);
BENCHMARK_TEMPLATE(BM_Div, tablePol32);


static void BM_RandomTime(benchmark::State& state) {
    std::uniform_int_distribution<uint32_t> uid(0, 255);
    std::default_random_engine rd;
    for (auto _ : state) {
        benchmark::DoNotOptimize(uid(rd));
    }
}
BENCHMARK(BM_RandomTime);

BENCHMARK_MAIN();
