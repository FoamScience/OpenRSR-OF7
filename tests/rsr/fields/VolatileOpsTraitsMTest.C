#include <type_traits>
#include <vector>
#include "catch.H"
#include "error.H"
#include "fvCFD.H"
#include "VolatileDimensionedFieldOpsTraitsM.H"
#include "scalar.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define CHECK_FUNCTIONS(T, opName, op, opFunc, func)                           \
CHECK(decltype(opName##Helper::operator##opName##Works(res,a,b))::value == op);\
CHECK(decltype(opName##Helper::opFunc##Works(res, a ,b))::value == func);      \
REQUIRE_NOTHROW(opName##Helper::opFunc##Best(res, a, b))

#define CHECH_FUNCTIONS_RETURN(T, opName, opFunc, syntax)                      \
    using rT1 = std::remove_reference<decltype(syntax)>::type;                 \
    using rT2 = std::remove_reference                                          \
        <decltype(opName##Helper::opFunc##Best(res, a, b))>::type;             \
    REQUIRE( std::is_same<rT1, rT2>::value )

#define ADD_OPFUNC(T)                                                          \
void add(T& res, const T& a, const T& b)                                       \
{                                                                              \
    add(static_cast<Field<scalar>&>(res), a, b);                               \
}

using namespace Foam;

struct A { double a = 0; };
A operator+(const A& aa, const A& bb)
{
    A res;
    res.a = aa.a + bb.a;
    return res;
}

struct B { double b = 1; };
void add(B& res, const B& aa, const B& bb)
{
    res.b = aa.b + bb.b;
}

// Inherits from B and convertible to A
struct C : public B
{
    double c = 2;
    C() = default;
    C(const A& aa) : c(aa.a) {}

    operator A() { A a; a.a = c; return a; }
};

BEST_BINARY_PRODUCT(+, Plus, add, );

SCENARIO("Functional traits for volatile fields")
{
    GIVEN("A set of sample classes")
    {
        WHEN("Only the syntax (res = x op y) works")
        {
            A res, a, b;
            THEN("(x op y) is called when calling opFunc##Best()")
            {
                CHECK_FUNCTIONS(A, Plus, 1, add, 0);
                CHECH_FUNCTIONS_RETURN(C, Plus, add, (res = a+b));
            }
        }
        WHEN("Only the syntax opFunc(res, x, y) works")
        {
            B res, a, b;
            THEN("opFunc(res, x, y) is called when calling opFunc##Best()")
            {
                CHECK_FUNCTIONS(B, Plus, 0, add, 1);
                CHECH_FUNCTIONS_RETURN(C, Plus, add, add(res,a,b));
            }
        }
        WHEN("A Class is convertible to classes that work with both functions")
        {
            C res, a, b;
            res = a + b;
            THEN("Using the operator syntax is preferred")
            {
                CHECK_FUNCTIONS(C, Plus, 1, add, 1);
                CHECH_FUNCTIONS_RETURN(C, Plus, add, (res = a+b));
            }
        }
    }
}

// ************************************************************************* //
