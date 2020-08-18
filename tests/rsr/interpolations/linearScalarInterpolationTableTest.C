#include "catch.H"
#include "basicInterpolationTable.H"
#include "scalarList.H"
#include "error.H"

using namespace Foam;

SCENARIO("Linear interpolation of a list of scalar values")
{
    // Setup required/sane entries for a canonical dictionary
    dictionary dict;
    dict.add<word>("interpolationType", "linear"); // Default
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
				expected[1] = scalarList{4.5, 1, 3.5, 5, 3, 3};
				expected[2] = scalarList{8.1, 1, 4.7, 6.6, 2.2, 4.6};
				expected[3] = scalarList
				{6.27272727, 6.45454545, 7.72727273, 7.90909091, 1.09090909, 2.27272727};
				expected[4] = lit->values()[lit->values().size()-1].second();

                // Compare vector of scalarLists with custom predicate
                REQUIRE_THAT(
                    calced,
                    Catch::Matchers::Predicate<std::vector<scalarList>>
                    (
                        [&expected]
                        (const std::vector<scalarList>& actual) -> bool {
                        for (unsigned long i = 0; i < actual.size(); ++i) {
                            forAll(actual[i], ci)
                            {
                                if
                                (
                                    actual[i][ci] != Approx(expected[i][ci])
                                )
                                {
                                    return false;
                                }
                            }
                        }
                        return true; }
                    )
                );
			}
		}
	}
}
