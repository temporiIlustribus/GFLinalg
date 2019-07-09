#include <iostream>
#include <random>
#include <vector>
#include <string>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "GFlinalg.hpp"

typedef GFlinalg::BasicBinPolynomial<uint8_t, 11> basicPol;
typedef GFlinalg::PowBinPolynomial<uint8_t, 11> powPol;
typedef GFlinalg::TableBinPolynomial<uint8_t, 11> tablePol;
//
typedef GFlinalg::BasicGFElem<uint8_t> basicElem;
typedef GFlinalg::PowGFElem<uint8_t> powElem;
// Comment this if you are having problems building project
template<>
const powPol::LUTPair powPol::alphaToIndex = powPol::makeAlphaToIndex();
template<>
const tablePol::GFtable tablePol::mulTable = tablePol::makeMulTable();
template<>
const tablePol::GFtable tablePol::divTable = tablePol::makeInvMulTable();
//
std::shared_ptr<GFlinalg::PowGFElem<uint8_t>::LUTPair> LUT(GFlinalg::PowGFElem<uint8_t>::LUTPair::makeLUT(8, 11));
#define lut1 LUT.get()

TEST_CASE("Basic arithmetic", "[BasicBinPolynomial]") {
    SECTION("Reduction") {
        REQUIRE(basicPol(10).val() == 1);
        REQUIRE(basicPol(11).val() == 0);
        REQUIRE(basicPol(1).val() == 1);
        REQUIRE(basicPol(42).val() == 6);
        REQUIRE(basicPol(9).val() == 2);
    }
    SECTION("Data access") {
        basicPol a(10);
        REQUIRE(a.val() == 1);
        REQUIRE(a.val() == 1);
        REQUIRE(basicPol(a).val() == 1);
        REQUIRE(basicPol(a.val() + 2).val() == 3);
        REQUIRE(a.gfDegree() == 3);
    }
    SECTION("Addition") {
        basicPol a(10);
        basicPol b(1);
        REQUIRE((tablePol{1}+a).val() == 0);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a + b).val() == 0);
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
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((tablePol{1} *a).val() == 1);
        REQUIRE((a * b).val() == 1);
        REQUIRE(basicPol(42) * basicPol(42) == basicPol(2));
        REQUIRE(basicPol(42) * basicPol(0) == basicPol(0));
        REQUIRE(basicPol(3) * basicPol(3) == basicPol(5));
        REQUIRE(basicPol(7) * basicPol(4) == basicPol(1));
        REQUIRE(basicPol(5) * basicPol(3) == basicPol(4));
        REQUIRE((a *= basicPol(40)) == basicPol(4));
        REQUIRE(a.val() == 4);
    }
    SECTION("Division") {
        basicPol a(10);
        basicPol b(1);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a / b).val() == 1);
        REQUIRE((basicPol(2) / basicPol(6)).val() == 6);
        REQUIRE((basicPol(6) / basicPol(6)).val() == 1);
        REQUIRE((basicPol(10) / basicPol(7)).val() == 4);
        REQUIRE((basicPol(10) / basicPol(4)).val() == 7);
        REQUIRE((basicPol(4) / basicPol(5)).val() == 3);
        REQUIRE((basicPol(4) / basicPol(8)).val() == 5);
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
        REQUIRE(powPol(10).val() == 1);
        REQUIRE(powPol(11).val() == 0);
        REQUIRE(powPol(1).val() == 1);
        REQUIRE(powPol(42).val() == 6);
        REQUIRE(powPol(9).val() == 2);
    }
    SECTION("Data access") {
        powPol a(10);
        REQUIRE(a.val() == 1);
        REQUIRE(a.val() == 1);
        a.val() += 2;
        REQUIRE(a.val() == 3);
        REQUIRE(a.gfDegree() == 3);
    }
    SECTION("Addition") {
        powPol a(10);
        powPol b(1);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a + b).val() == 0);
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
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a * b).val() == 1);
        REQUIRE(powPol(42) * powPol(42) == powPol(2));
        REQUIRE(powPol(42) * powPol(0) == powPol(0));
        REQUIRE(powPol(3) * powPol(3) == powPol(5));
        REQUIRE(powPol(7) * powPol(4) == powPol(1));
        REQUIRE(powPol(5) * powPol(3) == powPol(4));
        REQUIRE((a *= powPol(40)) == powPol(4));
        REQUIRE(a.val() == 4);
    }
    SECTION("Division") {
        powPol a(10);
        powPol b(1);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a / b).val() == 1);
        REQUIRE((powPol(2) / powPol(6)).val() == 6);
        REQUIRE((powPol(6) / powPol(6)).val() == 1);
        REQUIRE((powPol(10) / powPol(7)).val() == 4);
        REQUIRE((powPol(10) / powPol(4)).val() == 7);
        REQUIRE((powPol(4) / powPol(5)).val() == 3);
        REQUIRE((powPol(4) / powPol(8)).val() == 5);
    }
    SECTION("Galois Power") {
        REQUIRE(GFlinalg::pow(powPol(10), 2) == powPol(1));
        REQUIRE(GFlinalg::pow(powPol(15), 3) == powPol(5));
        REQUIRE(GFlinalg::pow(powPol(3), 3) == powPol(4));
        REQUIRE(GFlinalg::pow(powPol(42), 7) == powPol(1));
        REQUIRE(GFlinalg::pow(powPol(42), 8) == powPol(42));
        REQUIRE(GFlinalg::pow(powPol(42), 6) * powPol(42) == powPol(1));
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
        REQUIRE(tablePol(10).val() == 1);
        REQUIRE(tablePol(11).val() == 0);
        REQUIRE(tablePol(1).val() == 1);
        REQUIRE(tablePol(42).val() == 6);
        REQUIRE(tablePol(9).val() == 2);
    }
    SECTION("Data access") {
        tablePol a(10);
        REQUIRE(a.val() == 1);
        REQUIRE(a.val() == 1);
        a.val() += 2;
        REQUIRE(a.val() == 3);
        REQUIRE(a.gfDegree() == 3);
    }
    SECTION("Addition") {
        tablePol a(10);
        tablePol b(1);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
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
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a * b).val() == 1);
        REQUIRE(tablePol(42) * tablePol(42) == tablePol(2));
        REQUIRE(tablePol(42) * tablePol(0) == tablePol(0));
        REQUIRE(tablePol(3) * tablePol(3) == tablePol(5));
        REQUIRE(tablePol(7) * tablePol(4) == tablePol(1));
        REQUIRE(tablePol(5) * tablePol(3) == tablePol(4));
        REQUIRE((a *= tablePol(40)) == tablePol(4));
        REQUIRE(a.val() == 4);
    }
    SECTION("Division") {
        tablePol::makeInvMulTable();
        tablePol a(10);
        tablePol b(1);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a / b).val() == 1);
        REQUIRE((tablePol(2) / tablePol(6)).val() == 6);
        REQUIRE((tablePol(6) / tablePol(6)).val() == 1);
        REQUIRE((tablePol(10) / tablePol(7)).val() == 4);
        REQUIRE((tablePol(10) / tablePol(4)).val() == 7);
        REQUIRE((tablePol(4) / tablePol(5)).val() == 3);
        REQUIRE((tablePol(4) / tablePol(8)).val() == 5);
    }
    SECTION("Galois Power") {
        REQUIRE(GFlinalg::pow(tablePol(10), 2) == tablePol(1));
        REQUIRE(GFlinalg::pow(tablePol(15), 3) == tablePol(5));
        REQUIRE(GFlinalg::pow(tablePol(3), 3) == tablePol(4));
        REQUIRE(GFlinalg::pow(tablePol(42), 7) == tablePol(1));
        REQUIRE(GFlinalg::pow(tablePol(42), 8) == tablePol(42));
        REQUIRE(GFlinalg::pow(tablePol(42), 6) * tablePol(42) == tablePol(1));
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

TEST_CASE("Basic single template param arithmetic", "[BasicGFElem]") {
    SECTION("Reduction") {
        REQUIRE(static_cast<size_t>(basicElem(10, 11).val()) == 1);
        REQUIRE(basicElem(11, 11).val() == 0);
        REQUIRE(basicElem(1, 11).val() == 1);
        REQUIRE(basicElem(42, 11).val() == 6);
        REQUIRE(basicElem(9, 11).val() == 2);
    }
    SECTION("Data access") {
        basicElem a(10, 11);
        REQUIRE(a.val() == 1);
        REQUIRE(basicElem(a).val() == 1);
        REQUIRE(basicElem(a.val() + 2, 11).val() == 3);
        REQUIRE(a.gfDegree() == 3);
    }
    SECTION("Addition") {
        basicElem a(10, 11);
        basicElem b(1, 11);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a + b).val() == 0);
        REQUIRE(basicElem(10, 11) + basicElem(1, 11) == basicElem(0, 11));
        REQUIRE(basicElem(42, 11) + basicElem(5, 11) == basicElem(3, 11));
        REQUIRE(basicElem(42, 11) + basicElem(0, 11) == basicElem(94, 11));
        REQUIRE(basicElem(8, 11) + basicElem(3, 11) == basicElem(0, 11));
        REQUIRE((a += basicElem(6, 11)) == basicElem(7, 11));
        REQUIRE(a == basicElem(7, 11));
        REQUIRE(a + basicElem(17, 11) == basicElem(0, 11));
    }
    SECTION("Multiplication") {
        basicElem a(10, 11);
        basicElem b(1, 11);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a * b).val() == 1);
        REQUIRE(basicElem(42, 11) * basicElem(42, 11) == basicElem(2, 11));
        REQUIRE(basicElem(42, 11) * basicElem(0, 11) == basicElem(0, 11));
        REQUIRE(basicElem(3, 11) * basicElem(3, 11) == basicElem(5, 11));
        REQUIRE(basicElem(7, 11) * basicElem(4, 11) == basicElem(1, 11));
        REQUIRE(basicElem(5, 11) * basicElem(3, 11) == basicElem(4, 11));
        REQUIRE((a *= basicElem(40, 11)) == basicElem(4, 11));
        REQUIRE(a.val() == 4);
    }
    SECTION("Division") {
        basicElem a(10, 11);
        basicElem b(1, 11);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a / b).val() == 1);
        REQUIRE((basicElem(2, 11) / basicElem(6, 11)).val() == 6);
        REQUIRE((basicElem(6, 11) / basicElem(6, 11)).val() == 1);
        REQUIRE((basicElem(10, 11) / basicElem(7, 11)).val() == 4);
        REQUIRE((basicElem(10, 11) / basicElem(4, 11)).val() == 7);
        REQUIRE((basicElem(4, 11) / basicElem(5, 11)).val() == 3);
        REQUIRE((basicElem(4, 11) / basicElem(8, 11)).val() == 5);
    }
    SECTION("Galois Power") {
        REQUIRE(GFlinalg::pow(basicElem(10, 11), 2) == basicElem(1, 11));
        REQUIRE(GFlinalg::pow(basicElem(15, 11), 3) == basicElem(5, 11));
        REQUIRE(GFlinalg::pow(basicElem(3, 11), 3) == basicElem(4, 11));
        REQUIRE(GFlinalg::pow(basicElem(42, 11), 7) == basicElem(1, 11));
        REQUIRE(GFlinalg::pow(basicElem(42, 11), 8) == basicElem(42, 11));
        REQUIRE(basicElem(42, 11) * GFlinalg::pow(basicElem(42, 11), 6) == basicElem(1, 11));
    }
    SECTION("Compare operators") {
        REQUIRE(basicElem(42, 11) > basicElem(5, 11));
        REQUIRE(basicElem(5, 11) > basicElem(4, 11));
        REQUIRE(basicElem(10, 11) < basicElem(2, 11));
        REQUIRE(basicElem(42, 11) < basicElem(17, 11));
        REQUIRE(basicElem(42, 11) > basicElem(128, 11));
        REQUIRE(basicElem(176, 11) >= basicElem(0, 11));
        REQUIRE(basicElem(176, 11) <= basicElem(0, 11));

    }
}

TEST_CASE("Pow single template param arithmetic", "[PowGFElem]") {
    SECTION("Reduction") {
        REQUIRE(static_cast<size_t>(powElem(10, 11, lut1).val()) == 1);
        REQUIRE(powElem(11, 11, lut1).val() == 0);
        REQUIRE(powElem(1, 11, lut1).val() == 1);
        REQUIRE(powElem(42, 11, lut1).val() == 6);
        REQUIRE(powElem(9, 11, lut1).val() == 2);
    }
    SECTION("Data access") {
        powElem a(10, 11, lut1);
        REQUIRE(a.val() == 1);
        REQUIRE(powElem(a, lut1).val() == 1);
        REQUIRE(powElem(a.val() + 2, 11, lut1).val() == 3);
        REQUIRE(a.gfDegree() == 3);
    }
    SECTION("Addition") {
        powElem a(10, 11);
        powElem b(1, 11);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a + b).val() == 0);
        REQUIRE(powElem(10, 11) + powElem(1, 11) == powElem(0, 11));
        REQUIRE(powElem(42, 11) + powElem(5, 11) == powElem(3, 11));
        REQUIRE(powElem(42, 11) + powElem(0, 11) == powElem(94, 11));
        REQUIRE(powElem(8, 11) + powElem(3, 11) == powElem(0, 11));
        REQUIRE((a += powElem(6, 11)) == powElem(7, 11));
        REQUIRE(a == powElem(7, 11));
        REQUIRE(a + powElem(17, 11) == powElem(0, 11));
    }
    SECTION("Multiplication") {
        powElem a(10, 11, lut1);
        powElem b(1, 11, lut1);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a * b).val() == 1);
        REQUIRE(powElem(42, 11, lut1) * powElem(42, 11, lut1) == powElem(2, 11, lut1));
        REQUIRE(powElem(42, 11, lut1) * powElem(0, 11, lut1) == powElem(0, 11, lut1));
        REQUIRE(powElem(3, 11, lut1) * powElem(3, 11, lut1) == powElem(5, 11, lut1));
        REQUIRE(powElem(7, 11, lut1) * powElem(4, 11, lut1) == powElem(1, 11, lut1));
        REQUIRE(powElem(5, 11, lut1) * powElem(3, 11) == powElem(4, 11, lut1));
        REQUIRE((a *= powElem(40, 11, lut1)) == powElem(4, 11, lut1));
        REQUIRE(a.val() == 4);
    }
    SECTION("Division") {
        powElem a(10, 11, lut1);
        powElem b(1, 11, lut1);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a / b).val() == 1);
        REQUIRE((powElem(2, 11, lut1) / powElem(6, 11, lut1)).val() == 6);
        REQUIRE((powElem(6, 11, lut1) / powElem(6, 11, lut1)).val() == 1);
        REQUIRE((powElem(10, 11, lut1) / powElem(7, 11, lut1)).val() == 4);
        REQUIRE((powElem(10, 11, lut1) / powElem(4, 11, lut1)).val() == 7);
        REQUIRE((powElem(4, 11, lut1) / powElem(5, 11, lut1)).val() == 3);
        REQUIRE((powElem(4, 11, lut1) / powElem(8, 11, lut1)).val() == 5);
    }
    SECTION("Galois Power") {
        REQUIRE(GFlinalg::pow(powElem(10, 11, lut1), 2) == powElem(1, 11, lut1));
        REQUIRE(GFlinalg::pow(powElem(15, 11, lut1), 3) == powElem(5, 11, lut1));
        REQUIRE(GFlinalg::pow(powElem(3, 11, lut1), 3) == powElem(4, 11, lut1));
        REQUIRE(GFlinalg::pow(powElem(42, 11, lut1), 7) == powElem(1, 11, lut1));
        REQUIRE(GFlinalg::pow(powElem(42, 11, lut1), 8) == powElem(42, 11, lut1));
        REQUIRE(powElem(42, 11, lut1) * GFlinalg::pow(powElem(42, 11, lut1), 6) == powElem(1, 11, lut1));
    }
    SECTION("Compare operators") {
        REQUIRE(powElem(42, 11, lut1) > powElem(5, 11, lut1));
        REQUIRE(powElem(5, 11, lut1) > powElem(4, 11, lut1));
        REQUIRE(powElem(10, 11, lut1) < powElem(2, 11, lut1));
        REQUIRE(powElem(42, 11, lut1) < powElem(17, 11, lut1));
        REQUIRE(powElem(42, 11, lut1) > powElem(128, 11, lut1));
        REQUIRE(powElem(176, 11, lut1) >= powElem(0, 11, lut1));
        REQUIRE(powElem(176, 11, lut1) <= powElem(0, 11, lut1));
    }
}




