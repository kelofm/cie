// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/testing/inc/essentials.hpp"
#include "packages/commandline/inc/ArgParse.hpp"
#include "packages/stl_extension/inc/zip.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"


namespace cie::utils {


CIE_TEST_CASE("ArgParse", "[commandline]")
{
    CIE_TEST_CASE_INIT("ArgParse")

    {
        CIE_TEST_REQUIRE_NOTHROW(ArgParse());
        ArgParse parser;

        CIE_TEST_CHECK_NOTHROW(parser.addPositional("pos"));
        CIE_TEST_CHECK_NOTHROW(parser.addKeyword({"-s", "--start"}));
        CIE_TEST_CHECK_NOTHROW(parser.addFlag({"-k"}));
        CIE_TEST_CHECK_NOTHROW(parser.addFlag({"-r"}));
        CIE_TEST_CHECK_NOTHROW(parser.addKeyword({"-v"}));

        const int argc = 6;
        char const* argv[argc] = {
            "a_positional_argument",
            "--start",
            "a_keyword_argument",
            "-k",
            "-v",
            "1"
        };

        ArgParse::Results results;
        CIE_TEST_CHECK_NOTHROW(results = parser.parseArguments(argc,argv));

        const std::map<std::string,DynamicArray<std::string>> referenceResults {
            {"pos", {"a_positional_argument"}},
            {"s", {"a_keyword_argument"}},
            {"k", {"true"}},
            {"r", {"false"}},
            {"v", {"1"}}
        };

        for (const auto& r_pair : referenceResults)
        {
            const auto pairRange = results[r_pair.first];

            CIE_TEST_CHECK(!pairRange.empty());

            auto it = pairRange.begin();
            auto it_ref = r_pair.second.begin();
            const auto it_end = pairRange.end();
            for (; it!=it_end; ++it, ++it_ref)
            {
                CIE_TEST_REQUIRE(it_ref != r_pair.second.end());
                CIE_TEST_CHECK(it->second == *it_ref);
            }
            CIE_TEST_CHECK(it_ref == r_pair.second.end());
        }

        CIE_TEST_CHECK_NOTHROW(results.get<int>("v"));
        CIE_TEST_CHECK(results.get<int>("v") == 1);
        CIE_TEST_CHECK_NOTHROW(results.get<Size>("v"));
        CIE_TEST_CHECK(results.get<Size>("v") == 1);
        CIE_TEST_CHECK_NOTHROW(results.get<double>("v"));
        CIE_TEST_CHECK(results.get<double>("v") == 1.0);
        CIE_TEST_CHECK_NOTHROW(results.get<std::string>("v"));
        CIE_TEST_CHECK(results.get<std::string>("v") == "1");
    }

    {
        CIE_TEST_REQUIRE_NOTHROW(ArgParse(""));
        ArgParse parser("ParserTest");

        CIE_TEST_CHECK_NOTHROW(parser.addPositional("pos", "this is a positional argument"));
        CIE_TEST_CHECK_NOTHROW(parser.addKeyword({"-s", "--start"}, 4));
        CIE_TEST_CHECK_NOTHROW(parser.addFlag({"-k"}));
        CIE_TEST_CHECK_NOTHROW(parser.addFlag({"-r"}));
        CIE_TEST_CHECK_NOTHROW(parser.addKeyword(
            {"-v"},
            [](const ArgParse::ValueView& r_value) -> bool {return r_value.size() == 1 && *r_value.begin() == '1';}
        ));

        const int argc = 9;
        char const* argv[argc] = {
            "a_positional_argument",
            "--start",
            "a_keyword_argument_belonging_to_start",
            "another_one",
            "and_another_one",
            "and a last one",
            "-k",
            "-v",
            "1"
        };

        ArgParse::Results results;
        CIE_TEST_CHECK_NOTHROW(results = parser.parseArguments(argc,argv));
        CIE_TEST_CHECK_NOTHROW(parser.help(OutputStream::cout));
        CIE_TEST_CHECK_NOTHROW(OutputStream::cout << results << std::endl);
    }
}


} // namespace cie::utils
