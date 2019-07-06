#include <iostream>
#include <random>
#include <vector>
#include <string>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "GFlinalg.hpp"

typedef GFlinalg::PowBinPolynomial<uint8_t, 11> powPol;
typedef GFlinalg::BasicBinPolynomial<uint8_t, 11> basicPol;
typedef GFlinalg::TableBinPolynomial<uint8_t, 11> tablePol;
// Comment this if you are having problems building project
template<>
powPol::ArrayPair powPol::alphaToIndex = powPol::makeAlphaToIndex();
template<>
tablePol::GFtable tablePol::mulTable = tablePol::makeMulTable();
template<>
tablePol::GFtable tablePol::divTable = tablePol::makeInvMulTable();

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
        REQUIRE(basicPol(a).getVal() == 1);
        REQUIRE(basicPol(a.getVal() + 2).getVal() == 3);
        REQUIRE(a.gfSize() == 3);
    }
    SECTION("Addition") {
        basicPol a(10);
        basicPol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a+b).getVal() == 0);
        REQUIRE(basicPol(10) + basicPol(1) == basicPol(0));
        REQUIRE(basicPol(42) + basicPol(5) == basicPol(3));
        REQUIRE(basicPol(42) + basicPol(0) == basicPol(94));
        REQUIRE(basicPol(8) + basicPol(3) == basicPol(0));
        REQUIRE((a += basicPol(6)) == basicPol(7));
        REQUIRE(a == basicPol(7));
        REQUIRE(a + basicPol(17) == basicPol(0));
    }
    SECTION("Multiplication") {
        basicPol a(10);
        basicPol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a * b).getVal() == 1);
        REQUIRE(basicPol(42) * basicPol(42) == basicPol(2));
        REQUIRE(basicPol(42) * basicPol(0) == basicPol(0));
        REQUIRE(basicPol(3) * basicPol(3) == basicPol(5));
        REQUIRE(basicPol(7) * basicPol(4) == basicPol(1));
        REQUIRE(basicPol(5) * basicPol(3) == basicPol(4));
        REQUIRE((a *= basicPol(40)) == basicPol(4));
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
    SECTION("Galois Power") {
        REQUIRE(GFlinalg::pow(basicPol(10), 2) == basicPol(1));
        REQUIRE(GFlinalg::pow(basicPol(15), 3) == basicPol(5));
        REQUIRE(GFlinalg::pow(basicPol(3), 3) == basicPol(4));
        REQUIRE(GFlinalg::pow(basicPol(42), 7) == basicPol(1));
        REQUIRE(GFlinalg::pow(basicPol(42), 8) == basicPol(42));
        REQUIRE(basicPol(42) * GFlinalg::pow(basicPol(42), 6) == basicPol(1));
    }
    SECTION("Compare operators") {
        REQUIRE(basicPol(42) > basicPol(5));
        REQUIRE(basicPol(5) > basicPol(4));
        REQUIRE(basicPol(10) < basicPol(2));
        REQUIRE(basicPol(42) < basicPol(17));
        REQUIRE(basicPol(42) > basicPol(128));
        REQUIRE(basicPol(176) >= basicPol(0));
        REQUIRE(basicPol(176) <= basicPol(0));
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
        REQUIRE(a.gfSize() == 3);
    }
    SECTION("Addition") {
        powPol a(10);
        powPol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a + b).getVal() == 0);
        REQUIRE(powPol(10) + powPol(1) == powPol(0));
        REQUIRE(powPol(42) + powPol(5) == powPol(3));
        REQUIRE(powPol(42) + powPol(0) == powPol(94));
        REQUIRE(powPol(8) + powPol(3) == powPol(0));
        REQUIRE((a += powPol(6)) == powPol(7));
        REQUIRE(a == powPol(7));
        REQUIRE(a + powPol(17) == powPol(0));
    }
    SECTION("Multiplication") {
        powPol a(10);
        powPol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a * b).getVal() == 1);
        REQUIRE(powPol(42) * powPol(42) == powPol(2));
        REQUIRE(powPol(42) * powPol(0) == powPol(0));
        REQUIRE(powPol(3) * powPol(3) == powPol(5));
        REQUIRE(powPol(7) * powPol(4) == powPol(1));
        REQUIRE(powPol(5) * powPol(3) == powPol(4));
        REQUIRE((a *= powPol(40)) == powPol(4));
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
    SECTION("Galois Power") {
        REQUIRE(GFlinalg::pow(powPol(10), 2) == powPol(1));
        REQUIRE(GFlinalg::pow(powPol(15), 3) == powPol(5));
        REQUIRE(GFlinalg::pow(powPol(3), 3) == powPol(4));
        REQUIRE(GFlinalg::pow(powPol(42), 7) == powPol(1));
        REQUIRE(GFlinalg::pow(powPol(42), 8) == powPol(42));
        REQUIRE(powPol(42) * GFlinalg::pow(powPol(42), 6) == powPol(1));
    }
    SECTION("Compare operators") {
        REQUIRE(powPol(42) > powPol(5));
        REQUIRE(powPol(5) > powPol(4));
        REQUIRE(powPol(10) < powPol(2));
        REQUIRE(powPol(42) < powPol(17));
        REQUIRE(powPol(42) > powPol(128));
        REQUIRE(powPol(176) >= powPol(0));
        REQUIRE(powPol(176) <= powPol(0));
    }
}

TEST_CASE("Table arithmetic", "[TableBinPolynomial]") {
    SECTION("Reduction") {
        REQUIRE(tablePol(10).getVal() == 1);
        REQUIRE(tablePol(11).getVal() == 0);
        REQUIRE(tablePol(1).getVal() == 1);
        REQUIRE(tablePol(42).getVal() == 6);
        REQUIRE(tablePol(9).getVal() == 2);
    }
    SECTION("Data access") {
        tablePol a(10);
        REQUIRE(a.getVal() == 1);
        REQUIRE(a.val() == 1);
        a.val() += 2;
        REQUIRE(a.getVal() == 3);
        REQUIRE(a.gfSize() == 3);
    }
    SECTION("Addition") {
        tablePol a(10);
        tablePol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a + b) == tablePol(0));
        REQUIRE(tablePol(10) + tablePol(1) == tablePol(0));
        REQUIRE(tablePol(42) + tablePol(5) == tablePol(3));
        REQUIRE(tablePol(42) + tablePol(0) == tablePol(94));
        REQUIRE(tablePol(8) + tablePol(3) == tablePol(0));
        REQUIRE((a += tablePol(6)) == tablePol(7));
        REQUIRE(a == tablePol(7));
        REQUIRE(a + tablePol(17) == tablePol(0));
    }
    SECTION("Multiplication") {
        tablePol::makeMulTable();
        tablePol a(10);
        tablePol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a * b).getVal() == 1);
        REQUIRE(tablePol(42) * tablePol(42) == tablePol(2));
        REQUIRE(tablePol(42) * tablePol(0) == tablePol(0));
        REQUIRE(tablePol(3) * tablePol(3) == tablePol(5));
        REQUIRE(tablePol(7) * tablePol(4) == tablePol(1));
        REQUIRE(tablePol(5) * tablePol(3) == tablePol(4));
        REQUIRE((a *= tablePol(40)) == tablePol(4));
        REQUIRE(a.getVal() == 4);
    }
    SECTION("Division") {
        tablePol::makeInvMulTable();
        tablePol a(10);
        tablePol b(1);
        REQUIRE(a.getVal() == 1);
        REQUIRE(b.getVal() == 1);
        REQUIRE((a / b).getVal() == 1);
        REQUIRE((tablePol(2) / tablePol(6)).getVal() == 6);
        REQUIRE((tablePol(6) / tablePol(6)).getVal() == 1);
        REQUIRE((tablePol(10) / tablePol(7)).getVal() == 4);
        REQUIRE((tablePol(10) / tablePol(4)).getVal() == 7);
        REQUIRE((tablePol(4) / tablePol(5)).getVal() == 3);
        REQUIRE((tablePol(4) / tablePol(8)).getVal() == 5);
    }
    SECTION("Galois Power") {
        REQUIRE(GFlinalg::pow(tablePol(10), 2) == tablePol(1));
        REQUIRE(GFlinalg::pow(tablePol(15), 3) == tablePol(5));
        REQUIRE(GFlinalg::pow(tablePol(3), 3) == tablePol(4));
        REQUIRE(GFlinalg::pow(tablePol(42), 7) == tablePol(1));
        REQUIRE(GFlinalg::pow(tablePol(42), 8) == tablePol(42));
        REQUIRE(tablePol(42) * GFlinalg::pow(tablePol(42), 6) == tablePol(1));
    }
    SECTION("Compare operators") {
        REQUIRE(tablePol(42) > tablePol(5));
        REQUIRE(tablePol(5) > tablePol(4));
        REQUIRE(tablePol(10) < tablePol(2));
        REQUIRE(tablePol(42) < tablePol(17));
        REQUIRE(tablePol(42) > tablePol(128));
        REQUIRE(tablePol(176) >= tablePol(0));
        REQUIRE(tablePol(176) <= tablePol(0));
    }    
}
