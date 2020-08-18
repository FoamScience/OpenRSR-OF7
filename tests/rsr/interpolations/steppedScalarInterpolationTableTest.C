#include "catch.H"
#include "basicInterpolationTable.H"
#include "scalarList.H"
#include "error.H"

using namespace Foam;

SCENARIO("Stepped interpolation of a list of scalar values")
{
    // Setup required/sane entries for a canonical dictionary
    dictionary dict;
    dict.add<word>("interpolationType", "stepped"); // Default
    dict.add<bool>("periodic", false);
    dict.add<fileName>("file", "testData/testData.dat");

    FatalError.dontThrowExceptions();
	GIVEN("A valid linear interpolation table object for scalars")
	{
		auto lit = basicInterpolationTable<scalar>::New(dict);
		WHEN("interpolate() gets invoked")
		{
			std::vector<scalar> vec{1.0, 1.5, 1.9, 3.0, 3.1};
			std::vector<scalarList> calced(vec.size());
			for(unsigned int ci = 0; ci<vec.size(); ci++)
			{
				calced[ci] = lit->interpolate(vec[ci]);
			}

			THEN("it returns consistent interpolated values")
			{
				std::vector<scalarList> expected(vec.size());

				expected[0] = lit->values()[0].second();
				expected[1] = lit->values()[0].second();
				expected[2] = lit->values()[0].second();
				expected[3] = lit->values()[1].second();
				expected[4] = lit->values()[lit->values().size()-1].second();

                REQUIRE(expected == calced);
			}
		}
	}
}
