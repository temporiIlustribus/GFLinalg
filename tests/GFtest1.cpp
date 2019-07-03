#include <iostream>
#include <random>
#include <vector>
#include <string>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "GFlinalg.hpp"

typedef GFlinalg::PowBinPolynomial<uint8_t, 3, 11> powPol;
typedef GFlinalg::BaseBinPolynomial<uint8_t, 3, 11> basePol;
typedef GFlinalg::BasicBinPolynomial<uint8_t, 3, 11> basicPol;
template<>
powPol::arrayPair powPol::alphaToIndex = powPol::makeAlphaToIndex();

TEST_CASE("Basic reduction and data access", "[baseBinPolynomial]") {
    SECTION("Reduction") {
        REQUIRE(basePol(10).getVal() == 1);
        REQUIRE(basePol(11).getVal() == 0);
        REQUIRE(basePol(1).getVal() == 1);
        REQUIRE(basePol(42).getVal() == 6);
        REQUIRE(basePol(9).getVal() == 2);
    }
    SECTION("Data access") {
        basePol a(10);
        REQUIRE(a.getVal() == 1);
        REQUIRE(a.val() == 1);
        a.val() += 2;
        REQUIRE(a.getVal() == 3);
        REQUIRE(a.size() == 8);
        REQUIRE(a.GFsize() == 3);
    }
}

TEST_CASE("Basic arithmetic", "[BasicBinPolynomial]") {
    SECTION("Reduction") {
        REQUIRE(basicPol(10).getVal() == 1);
        REQUIRE(basicPol(11).getVal() == 0);
        REQUIRE(basicPol(1).getVal() == 1);
        REQUIRE(basicPol(42).getVal() == 6);
        REQUIRE(basicPol(9).getVal() == 2);
    }
    SECTION("Data access") {
        basicPol a(10);
        REQUIRE(a.getVal() == 1);
        REQUIRE(a.val() == 1);
        a.val() += 2;
        REQUIRE(a.getVal() == 3);
        REQUIRE(a.size() == 8);
        REQUIRE(a.GFsize() == 3);
    }
    SECTION("Addition") {
        basicPol a(10);
        basicPol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a+b).getVal() == 0);
        a.val() = 42;
        b.val() = 0;
        a.reduce();
        REQUIRE((basicPol{0} + basicPol{42}).getVal() == basicPol{42}.getVal());
        b.val() = 3;
        REQUIRE((a + b).getVal() == 5);
        b.val() = 10;
        REQUIRE((a+b).getVal() == 4);
        REQUIRE((a += b).getVal() == 4);
        REQUIRE(a.getVal() == 4);
        a.val() = 8;
        b.val() = 3;
        a.reduce();
        REQUIRE((a += b).getVal() == 0);
        REQUIRE(a.getVal() == 0);
    }
    SECTION("Multiplication") {
        basicPol a(10);
        basicPol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a * b).getVal() == 1);
        a.val() = 42;
        b.val() = 42;
        a.reduce();
        b.reduce();
        REQUIRE((a * b).getVal() == 2);
        a.val() = 0;
        REQUIRE((a * b).getVal() == 0);
        a.val() = 3;
        b.val() = 3;
        REQUIRE((a * b).getVal() == 5);
        a.val() = 7;
        b.val() = 4;
        REQUIRE((a * b).getVal() == 1);
        a.val() = 5;
        b.val() = 3;
        REQUIRE((a * b).getVal() == 4);
        REQUIRE((a *= b).getVal() == 4);
        REQUIRE(a.getVal() == 4);
    }
    SECTION("Division") {
        basicPol a(10);
        basicPol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a / b).getVal() == 1);
        REQUIRE((basicPol(2) / basicPol(6)).getVal() == 6);
        REQUIRE((basicPol(6) / basicPol(6)).getVal() == 1);
        REQUIRE((basicPol(10) / basicPol(7)).getVal() == 4);
        REQUIRE((basicPol(10) / basicPol(4)).getVal() == 7);
        REQUIRE((basicPol(4) / basicPol(5)).getVal() == 3);
        REQUIRE((basicPol(4) / basicPol(8)).getVal() == 5);
    }
}

TEST_CASE("Pow arithmetic", "[PowBinPolynomial]") {
    SECTION("Reduction") {
        REQUIRE(powPol(10).getVal() == 1);
        REQUIRE(powPol(11).getVal() == 0);
        REQUIRE(powPol(1).getVal() == 1);
        REQUIRE(powPol(42).getVal() == 6);
        REQUIRE(powPol(9).getVal() == 2);
    }
    SECTION("Data access") {
        powPol a(10);
        REQUIRE(a.getVal() == 1);
        REQUIRE(a.val() == 1);
        a.val() += 2;
        REQUIRE(a.getVal() == 3);
        REQUIRE(a.size() == 8);
        REQUIRE(a.GFsize() == 3);
    }
    SECTION("Addition") {
        powPol a(10);
        powPol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a + b).getVal() == 0);
        a.val() = 42;
        b.val() = 0;
        a.reduce();
        REQUIRE((a + b).getVal() == a.getVal());
        b.val() = 3;
        REQUIRE((a + b).getVal() == 3);
        b.val() = 10;
        REQUIRE((a + b).getVal() == 4);
        REQUIRE((a += b).getVal() == 4);
        REQUIRE(a.getVal() == 4);
        a.val() = 8;
        b.val() = 3;
        a.reduce();
        REQUIRE((a += b).getVal() == 0);
        REQUIRE(a.getVal() == 0);
    }
    SECTION("Multiplication") {
        powPol a(10);
        powPol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a * b).getVal() == 1);
        a.val() = 42;
        b.val() = 42;
        a.reduce();
        b.reduce();
        REQUIRE((a * b).getVal() == 2);
        a.val() = 0;
        REQUIRE((a * b).getVal() == 0);
        a.val() = 3;
        b.val() = 3;
        REQUIRE((a * b).getVal() == 5);
        a.val() = 7;
        b.val() = 4;
        REQUIRE((a * b).getVal() == 1);
        a.val() = 5;
        b.val() = 3;
        REQUIRE((a * b).getVal() == 4);
        REQUIRE((a *= b).getVal() == 4);
        REQUIRE(a.getVal() == 4);
    }
    SECTION("Division") {
        powPol a(10);
        powPol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a / b).getVal() == 1);
        REQUIRE((powPol(2) / powPol(6)).getVal() == 6);
        REQUIRE((powPol(6) / powPol(6)).getVal() == 1);
        REQUIRE((powPol(10) / powPol(7)).getVal() == 4);
        REQUIRE((powPol(10) / powPol(4)).getVal() == 7);
        REQUIRE((powPol(4) / powPol(5)).getVal() == 3);
        REQUIRE((powPol(4) / powPol(8)).getVal() == 5);
    }
}

