
#include "catch.H"
#include "scalarList.H"
#include "dictionary.H"
#include "sampleInterpolationTable.H"

using namespace Foam;

SCENARIO("Generic Interpolation Construction from list of values")
{
	scalarList sl0 {0, 1, 2, 3, 4};
	scalarList sl1 {9, 1, 5, 7, 2, 5};
	scalarList sl2 {5, 2, 4};
	scalarList sl3 {6, 7, 8, 8, 1, 2};

    GIVEN("A faulty list of Tuple2s where time is not monotically increasing")
    {
        WHEN("Instatiating time series with the decreasing-time")
		{
			List<Tuple2<scalar, scalarList>> lt
			{
				Tuple2<scalar, scalarList>{1.0, sl0},
				Tuple2<scalar, scalarList>{10.0, sl1},
				Tuple2<scalar, scalarList>{2.3, sl2},
				Tuple2<scalar, scalarList>{2.1, sl3},
			};

			THEN("Instantiation must fail")
			{
				REQUIRE_THROWS(sampleInterpolationTable<scalar>(lt, false));
			}
		}
		WHEN("Instatiating time series from values (duplicated time entry)")
		{
			List<Tuple2<scalar, scalarList>> lt
			{
				Tuple2<scalar, scalarList>{1.0, sl0},
				Tuple2<scalar, scalarList>{10.0, sl1},
				Tuple2<scalar, scalarList>{10.0, sl2},
				Tuple2<scalar, scalarList>{2.1, sl3},
			};

			THEN("Instantiation fails")
			{
			    REQUIRE_THROWS(sampleInterpolationTable<scalar>(lt, false));
			}
		}
	}

	GIVEN("A valid list of Tuple2 values")
	{
		List<Tuple2<scalar, scalarList>> lt
		{
			Tuple2<scalar, scalarList>{1.0, sl0},
			Tuple2<scalar, scalarList>{2.0, sl1},
			Tuple2<scalar, scalarList>{2.3, sl2},
			Tuple2<scalar, scalarList>{3.0, sl3},
		};

		WHEN("Instatiating time series from list of values")
		{
			sampleInterpolationTable<scalar> interpTable(lt);
			THEN("lookup() and projectTime must return cosistent results")
			{
				scalar t = 2.2;
				scalarList intValue = interpTable.interpolate(t);
			
                // Projected time
				REQUIRE(intValue[0] == t);
                // Index of next element
				REQUIRE(intValue[1] == 2.0);
			}
			THEN("lookup() doesn't accept less than startTime values")
			{
				scalar t = lt[0].first() - 1.0;
				REQUIRE_THROWS(interpTable.interpolate(t));
			}
			THEN("By default, time after endTime is not accepted")
			{
				scalar t = lt[lt.size()-1].first() + 1.0;
				REQUIRE_THROWS(interpTable.interpolate(t));
			}
		}
		WHEN("Instatiating time series from values (periodicity enabled)")
		{
			sampleInterpolationTable<scalar> interpTable(lt, true);
			THEN("After-endTime projection must return consistent time")
			{
				std::vector<scalar> vec
				{
					lt[lt.size()-1].first()+1.1,
					lt[lt.size()-1].first()+2.0,
					lt[lt.size()-1].first()+2.8,
					lt[lt.size()-1].first()+3.0,
					lt[lt.size()-1].first()+10.0,
					lt[lt.size()-1].first()+25.4
				};

				std::vector<scalar> res(vec.size());
				for(unsigned int ci = 0; ci<res.size(); ci++)
				{
					res[ci] = interpTable.interpolate(vec[ci])[0];
				}
				std::vector<scalar> expected{1.1, 2.0, 2.8, 3.0, 1.0, 1.4};
				REQUIRE_THAT(res, Catch::Matchers::Approx(expected).margin(1e-6));
			}
		}
	}
}

SCENARIO("Generic Interpolation Construction from a dictionary")
{
    GIVEN("A dictionary for interpolation table")
    {
        // Setup required/sane entries for a canonical dictionary
        dictionary dict;
        dict.add<word>("interpolationType", "sampleInterpolationTable");
        dict.add<bool>("isPeriodic", false);
        dict.add<word>("file", "testData/testData.csv");
        dict.add<bool>("hasHeaderLine", true);
        dict.add<label>("timeColumn", 0);
        dict.add<label>("valueColumns", 1);

        WHEN("Dictionary requests non-existent interpolation type")
        {
            dictionary newDict(dict);
            newDict.set("interpolationType", word("banana"));
            THEN("Interpolation table construction must fail")
            {
				REQUIRE_THROWS(basicInterpolationTable<scalar>::New(dict));
            }
        }
    }
}
