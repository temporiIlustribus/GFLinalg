#include <iostream>
#include <random>
#include <vector>
#include <string>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "GFTPlinalg.hpp"
#include "GFSPlinalg.hpp"

using accessor = GFlinalg::Accessor<uint8_t>;

typedef GFlinalg::BasicBinPolynomial<uint8_t, 11> basicPol;
typedef GFlinalg::PowBinPolynomial<uint8_t, 11> powPol;
typedef GFlinalg::TableBinPolynomial<uint8_t, 11> tablePol;
//
typedef GFlinalg::BasicGFElem<uint8_t> basicElem;
typedef GFlinalg::PowGFElem<uint8_t> powElem;
typedef GFlinalg::TableGFElem<uint8_t> tableElem;
// Comment this if you are having problems building project
template<>
const GFlinalg::LUTArrPair<uint8_t,11> powPol::alphaToIndex{};
template<>
const tablePol::GFtable tablePol::mulTable = tablePol::makeMulTable();
template<>
const tablePol::GFtable tablePol::divTable = tablePol::makeInvMulTable();
//
GFlinalg::LUTVectPair<uint8_t> LUT{11};
GFlinalg::LUTVectPair<uint8_t> const* lut1(&LUT);
std::vector<uint8_t> mul_table(tableElem::makeMulTable(11));
std::vector<uint8_t> const* mul = &mul_table;
std::vector<uint8_t> div_table(tableElem::makeInvMulTable(mul, 11));

std::vector<uint8_t> const* div1 = &div_table;

TEMPLATE_TEST_CASE("Basic arithmetic", "[template]", basicPol, powPol, tablePol) {
    SECTION("Reduction") {
        REQUIRE(TestType(10).val() == 1);
        REQUIRE(TestType(11).val() == 0);
        REQUIRE(TestType(1).val() == 1);
        REQUIRE(TestType(42).val() == 6);
        REQUIRE(TestType(9).val() == 2);
    }
    SECTION("Data access") {
        TestType a(10);
        REQUIRE(a.val() == 1);
        REQUIRE(a.val() == 1);
        REQUIRE(TestType(a).val() == 1);
        REQUIRE(TestType(a.val() + 2).val() == 3);
        REQUIRE(a.gfDegree() == 3);
    }
    SECTION("Addition") {
        TestType a(10);
        TestType b(1);
        REQUIRE((tablePol{1}+a).val() == 0);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a + b) == TestType(0));
        REQUIRE(TestType(10) + TestType(1) == TestType(0));
        REQUIRE(TestType(42) + TestType(5) == TestType(3));
        REQUIRE(TestType(42) + TestType(0) == TestType(94));
        REQUIRE(TestType(8) + TestType(3) == TestType(0));
        REQUIRE((a += TestType(6)) == TestType(7));
        REQUIRE(a == TestType(7));
        REQUIRE(a + TestType(17) == TestType(0));
    }
    SECTION("Multiplication") {
        TestType a(10);
        TestType b(1);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((tablePol{1} *a).val() == 1);
        REQUIRE((a * b) == TestType(1));
        REQUIRE(TestType(42) * TestType(42) == TestType(2));
        REQUIRE((TestType(42) * TestType(0)) == TestType(0));
        REQUIRE(TestType(3) * TestType(3) == TestType(5));
        REQUIRE(TestType(7) * TestType(4) == TestType(1));
        REQUIRE(TestType(5) * TestType(3) == TestType(4));
        REQUIRE((a *= TestType(40)) == TestType(4));
        REQUIRE(a.val() == 4);
    }
    SECTION("Multiplication chaining") {
        TestType a(3);
        TestType b(5);
        auto res{a * b};
        while (res != a) {
            res *= b;
            REQUIRE(res.degree() < TestType::gfOrder());
        }
    }
    SECTION("Division") {
        TestType a(10);
        TestType b(1);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a / b) == TestType(1));
        REQUIRE((TestType(2) / TestType(6)) == TestType(6));
        REQUIRE((TestType(6) / TestType(6)) == TestType(1));
        REQUIRE((TestType(10) / TestType(7)) == TestType(4));
        REQUIRE((TestType(10) / TestType(4)) == TestType(7));
        REQUIRE((TestType(4) / TestType(5)) == TestType(3));
        REQUIRE((TestType(4) / TestType(8)) == TestType(5));
    }
    SECTION("Galois Power") {
        REQUIRE(GFlinalg::pow(TestType(10), 2) == TestType(1));
        REQUIRE(GFlinalg::pow(TestType(15), 3) == TestType(5));
        REQUIRE(GFlinalg::pow(TestType(3), 3) == TestType(4));
        REQUIRE(GFlinalg::pow(TestType(42), 7) == TestType(1));
        REQUIRE(GFlinalg::pow(TestType(42), 8) == TestType(42));
        REQUIRE(TestType(42) * GFlinalg::pow(TestType(42), 6) == TestType(1));
    }
    SECTION("Compare operators") {
        REQUIRE(TestType(42) > TestType(5));
        REQUIRE(TestType(5) > TestType(4));
        REQUIRE(TestType(10) < TestType(2));
        REQUIRE(TestType(42) < TestType(17));
        REQUIRE(TestType(42) > TestType(128));
        REQUIRE(TestType(176) >= TestType(0));
        REQUIRE(TestType(176) <= TestType(0));
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

TEST_CASE("Table single template param arithmetic", "[TableGFElem]") {
    SECTION("Reduction") {
        REQUIRE(static_cast<size_t>(tableElem(10, 11, mul, div1).val()) == 1);
        REQUIRE(tableElem(11, 11, mul, div1).val() == 0);
        REQUIRE(tableElem(1, 11, mul, div1).val() == 1);
        REQUIRE(tableElem(42, 11, mul, div1).val() == 6);
        REQUIRE(tableElem(9, 11, mul, div1).val() == 2);
    }
    SECTION("Data access") {
        tableElem a(10, 11, mul, div1);
        REQUIRE(a.val() == 1);
        REQUIRE(tableElem(a, mul, div1).val() == 1);
        REQUIRE(tableElem(a.val() + 2, 11, mul, div1).val() == 3);
        REQUIRE(a.gfDegree() == 3);
    }
    SECTION("Addition") {
        tableElem a(10, 11);
        tableElem b(1, 11);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a + b).val() == 0);
        REQUIRE(tableElem(10, 11) + tableElem(1, 11) == tableElem(0, 11));
        REQUIRE(tableElem(42, 11) + tableElem(5, 11) == tableElem(3, 11));
        REQUIRE(tableElem(42, 11) + tableElem(0, 11) == tableElem(94, 11));
        REQUIRE(tableElem(8, 11) + tableElem(3, 11) == tableElem(0, 11));
        REQUIRE((a += tableElem(6, 11)) == tableElem(7, 11));
        REQUIRE(a == tableElem(7, 11));
        REQUIRE(a + tableElem(17, 11) == tableElem(0, 11));
    }
    SECTION("Multiplication") {
        tableElem a(10, 11, mul, div1);
        tableElem b(1, 11, mul, div1);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a * b).val() == 1);
        REQUIRE(tableElem(42, 11, mul, div1) * tableElem(42, 11, mul, div1) == tableElem(2, 11, mul, div1));
        REQUIRE(tableElem(42, 11, mul, div1) * tableElem(0, 11, mul, div1) == tableElem(0, 11, mul, div1));
        REQUIRE(tableElem(3, 11, mul, div1) * tableElem(3, 11, mul, div1) == tableElem(5, 11, mul, div1));
        REQUIRE(tableElem(7, 11, mul, div1) * tableElem(4, 11, mul, div1) == tableElem(1, 11, mul, div1));
        REQUIRE(tableElem(5, 11, mul, div1) * tableElem(3, 11) == tableElem(4, 11, mul, div1));
        REQUIRE((a *= tableElem(40, 11, mul, div1)) == tableElem(4, 11, mul, div1));
        REQUIRE(a.val() == 4);
    }
    SECTION("Division") {
        tableElem a(10, 11, mul, div1);
        tableElem b(1, 11, mul, div1);
        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a / b).val() == 1);
        REQUIRE((tableElem(2, 11, mul, div1) / tableElem(6, 11, mul, div1)).val() == 6);
        REQUIRE((tableElem(6, 11, mul, div1) / tableElem(6, 11, mul, div1)).val() == 1);
        REQUIRE((tableElem(10, 11, mul, div1) / tableElem(7, 11, mul, div1)).val() == 4);
        REQUIRE((tableElem(10, 11, mul, div1) / tableElem(4, 11, mul, div1)).val() == 7);
        REQUIRE((tableElem(4, 11, mul, div1) / tableElem(5, 11, mul, div1)).val() == 3);
        REQUIRE((tableElem(4, 11, mul, div1) / tableElem(8, 11, mul, div1)).val() == 5);
    }
    SECTION("Galois power") {
        REQUIRE(GFlinalg::pow(tableElem(10, 11, mul, div1), 2) == tableElem(1, 11, mul, div1));
        REQUIRE(GFlinalg::pow(tableElem(15, 11, mul, div1), 3) == tableElem(5, 11, mul, div1));
        REQUIRE(GFlinalg::pow(tableElem(3, 11, mul, div1), 3) == tableElem(4, 11, mul, div1));
        REQUIRE(GFlinalg::pow(tableElem(42, 11, mul, div1), 7) == tableElem(1, 11, mul, div1));
        REQUIRE(GFlinalg::pow(tableElem(42, 11, mul, div1), 8) == tableElem(42, 11, mul, div1));
        REQUIRE(tableElem(42, 11, mul, div1) * GFlinalg::pow(tableElem(42, 11, mul, div1), 6) == tableElem(1, 11, mul, div1));
    }
    SECTION("Compare operators") {
        REQUIRE(tableElem(42, 11, mul, div1) > tableElem(5, 11, mul, div1));
        REQUIRE(tableElem(5, 11, mul, div1) > tableElem(4, 11, mul, div1));
        REQUIRE(tableElem(10, 11, mul, div1) < tableElem(2, 11, mul, div1));
        REQUIRE(tableElem(42, 11, mul, div1) < tableElem(17, 11, mul, div1));
        REQUIRE(tableElem(42, 11, mul, div1) > tableElem(128, 11, mul, div1));
        REQUIRE(tableElem(176, 11, mul, div1) >= tableElem(0, 11, mul, div1));
        REQUIRE(tableElem(176, 11, mul, div1) <= tableElem(0, 11, mul, div1));
    }
}

TEST_CASE("GF elements storage", "[Accessor]") {
    SECTION("Empty") {
        REQUIRE(accessor().empty());
    }

    SECTION("Insertion") {
        accessor a;

        REQUIRE(accessor().empty());

        REQUIRE(!a.tryInsert(5));
        REQUIRE(a.tryInsert(basicElem(5, 11)));
        REQUIRE(a.size() == 1);

        a.clear();

        REQUIRE (a.empty());
    }

    SECTION("Iteration") {
        accessor a;

        REQUIRE(a.tryInsert(basicElem(5, 11)));
        REQUIRE(a.tryInsert(basicElem(6, 11)));
        REQUIRE(a.tryInsert(basicElem(7, 11)));
        REQUIRE(a.tryInsert(basicElem(8, 11)));
        REQUIRE(a.tryInsert(basicElem(9, 11)));

        int i = 5;

        using Elem = typename accessor::Value;

        for (Elem& elem : a) {
            REQUIRE(elem == i);
            ++i;
        }
    }
}