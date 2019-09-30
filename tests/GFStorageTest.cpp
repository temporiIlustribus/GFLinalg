#include "catch.hpp"
#include "GFStorage.h"

using T = uint8_t;
using Elem = GFlinalg::BasicGFElem<T>;
using Ref = GFlinalg::GFElemRef<Elem>;
using Ptr = GFlinalg::GFElemPtr<Elem>;

template<size_t R, size_t C>
using Matrix = GFlinalg::MatrixEngine<Elem, R, C>;

TEST_CASE("Smart reference", "[GFElemRef]") {
    SECTION("Data access") {
        Elem a(10, 11);

        Ref ref(a.val(), a.getState());

        REQUIRE(ref.val() == 1);
        REQUIRE(Ref(a).val() == 1);
        REQUIRE(a.gfDegree() == 3);
    }

    SECTION("Addition") {
        Elem eA(10, 11);
        Elem eB(1, 11);
        Ref a(eA);
        Ref b(eB);

        REQUIRE(a.val() == 1);
        REQUIRE(b.val() == 1);
        REQUIRE((a + b).val() == 0);

        REQUIRE(a + Elem(7, 11) == Elem(6, 11));
    }
    SECTION("Multiplication") {
        Elem eA(10, 11);
        Ref a(eA);

        REQUIRE((a *= Elem(40, 11)) == Elem(4, 11));
        REQUIRE(a.val() == 4);
    }
    SECTION("Division") {
        Elem eA(10, 11);
        Elem eB(1, 11);
        Ref a(eA);
        Ref b(eB);

        REQUIRE((a / b).val() == 1);

        Elem e2(2, 11);
        Elem e6(6, 11);
        Ref e2r(e2);
        Ref e6r(e6);

        REQUIRE((e2r / e6r).val() == 6);
    }
}

TEST_CASE("Smart pointer", "[GFElemPtr]") {
    SECTION("Conversion") {
        Elem a(10, 11);

        Ref ref(a.val(), a.getState());
        Ptr ptr = ref;

        REQUIRE(ptr->val() == ref.val());
        REQUIRE(*ptr == ref);
    }
}

TEST_CASE("Matrix engine", "[MatrixEngine]") {
    SECTION("Data access") {
        Elem a(10, 11);

        Matrix<2, 4> m(a.getState());

        REQUIRE(m.rows() == 2);
        REQUIRE(m.columns() == 4);
        REQUIRE(m.size() == 8);

        REQUIRE(m(0, 0) != a);

        m(0, 0) = a;

        REQUIRE(m(0, 0) == a);
    }
}