// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp"


namespace cie::ct {


CIE_TEST_CASE("Match", "[compile_time]")
{

    CIE_TEST_CASE_INIT("Match")

    {
        CIE_TEST_CASE_INIT("Any")

        CIE_TEST_CHECK(Match<int>::Any<char,void,double,int>);
        CIE_TEST_CHECK(Match<int>::Any<int>);
        CIE_TEST_CHECK(!Match<int>::Any<void,float,double>);
        CIE_TEST_CHECK(!Match<int>::Any<unsigned int>);
        CIE_TEST_CHECK(!Match<int>::Any<>);

        CIE_TEST_CHECK(Match<int,char>::Any<void,double,int>);
        CIE_TEST_CHECK(Match<int,char>::Any<int>);
        CIE_TEST_CHECK(!Match<int,char>::Any<void,float,double>);
        CIE_TEST_CHECK(!Match<int,char>::Any<unsigned int>);
        CIE_TEST_CHECK(!Match<int,char>::Any<>);

        CIE_TEST_CHECK(!Match<>::Any<int>);
        CIE_TEST_CHECK(!Match<int>::Any<>);
        CIE_TEST_CHECK(!Match<>::Any<>);
    }

    {
        CIE_TEST_CASE_INIT("None")

        CIE_TEST_CHECK(!Match<int>::None<char,void,double,int>);
        CIE_TEST_CHECK(!Match<int>::None<int>);
        CIE_TEST_CHECK(Match<int>::None<void,float,double>);
        CIE_TEST_CHECK(Match<int>::None<unsigned int>);
        CIE_TEST_CHECK(Match<int>::None<>);

        CIE_TEST_CHECK(!Match<int,char>::None<void,double,int>);
        CIE_TEST_CHECK(!Match<int,char>::None<int>);
        CIE_TEST_CHECK(Match<int,char>::None<void,float,double>);
        CIE_TEST_CHECK(Match<int,char>::None<unsigned int>);
        CIE_TEST_CHECK(Match<int,char>::None<>);

        CIE_TEST_CHECK(Match<>::None<>);
    }

    {
        CIE_TEST_CASE_INIT("All")

        CIE_TEST_CHECK(Match<int>::All<int>);
        CIE_TEST_CHECK(Match<int,float>::All<int,float>);
        CIE_TEST_CHECK(Match<int,float>::All<float,int>);
        CIE_TEST_CHECK(Match<int>::All<int,int>);
        CIE_TEST_CHECK(Match<int,int>::All<int>);
        CIE_TEST_CHECK(Match<int,int>::All<int,float>);
        CIE_TEST_CHECK(Match<int>::All<int,float>);

        CIE_TEST_CHECK(!Match<int>::All<float>);
        CIE_TEST_CHECK(!Match<int,void>::All<int,float>);
        CIE_TEST_CHECK(!Match<int,float>::All<int>);

        CIE_TEST_CHECK(!Match<>::All<void>);
        CIE_TEST_CHECK(!Match<void>::All<>);
        CIE_TEST_CHECK(Match<>::All<>);
    }

    {
        CIE_TEST_CASE_INIT("Identical")

        CIE_TEST_CHECK(Match<int>::Identical<int>);
        CIE_TEST_CHECK(Match<void>::Identical<void>);
        CIE_TEST_CHECK(Match<int,void,float>::Identical<int,void,float>);

        CIE_TEST_CHECK(!Match<int>::Identical<unsigned int>);
        CIE_TEST_CHECK(!Match<int>::Identical<int,void>);
        CIE_TEST_CHECK(!Match<int,void>::Identical<int>);
        CIE_TEST_CHECK(!Match<int,void>::Identical<void,int>);

        CIE_TEST_CHECK(Match<>::Identical<>);
        CIE_TEST_CHECK(!Match<int>::Identical<>);
        CIE_TEST_CHECK(!Match<>::Identical<int>);
    }
}


CIE_TEST_CASE("MatchTuple", "[compile_time]")
{

    CIE_TEST_CASE_INIT("MatchTuple")

    {
        CIE_TEST_CASE_INIT("Any")

        CIE_TEST_CHECK(MatchTuple<std::tuple<int>>::Any<std::tuple<char,void,double,int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int>>::Any<std::tuple<int>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int>>::Any<std::tuple<void,float,double>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int>>::Any<std::tuple<unsigned int>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int>>::Any<std::tuple<>>);

        CIE_TEST_CHECK(MatchTuple<std::tuple<int,char>>::Any<std::tuple<void,double,int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int,char>>::Any<std::tuple<int>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int,char>>::Any<std::tuple<void,float,double>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int,char>>::Any<std::tuple<unsigned int>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int,char>>::Any<std::tuple<>>);

        CIE_TEST_CHECK(!MatchTuple<std::tuple<>>::Any<std::tuple<>>);
    }

    {
        CIE_TEST_CASE_INIT("None")

        CIE_TEST_CHECK(!MatchTuple<std::tuple<int>>::None<std::tuple<char,void,double,int>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int>>::None<std::tuple<int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int>>::None<std::tuple<void,float,double>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int>>::None<std::tuple<unsigned int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int>>::None<std::tuple<>>);

        CIE_TEST_CHECK(!MatchTuple<std::tuple<int,char>>::None<std::tuple<void,double,int>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int,char>>::None<std::tuple<int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int,char>>::None<std::tuple<void,float,double>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int,char>>::None<std::tuple<unsigned int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int,char>>::None<std::tuple<>>);

        CIE_TEST_CHECK(MatchTuple<std::tuple<>>::None<std::tuple<>>);
    }

    {
        CIE_TEST_CASE_INIT("All")

        CIE_TEST_CHECK(MatchTuple<std::tuple<int>>::All<std::tuple<int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int,float>>::All<std::tuple<int,float>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int,float>>::All<std::tuple<float,int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int>>::All<std::tuple<int,int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int,int>>::All<std::tuple<int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int,int>>::All<std::tuple<int,float>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int>>::All<std::tuple<int,float>>);

        CIE_TEST_CHECK(!MatchTuple<std::tuple<int>>::All<std::tuple<float>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int,long int>>::All<std::tuple<int,float>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int,float>>::All<std::tuple<int>>);

        CIE_TEST_CHECK(!MatchTuple<std::tuple<>>::All<std::tuple<long int>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<long int>>::All<std::tuple<>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<>>::All<std::tuple<>>);
    }

    {
        CIE_TEST_CASE_INIT("Identical")

        CIE_TEST_CHECK(MatchTuple<std::tuple<int>>::Identical<std::tuple<int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<long int>>::Identical<std::tuple<long int>>);
        CIE_TEST_CHECK(MatchTuple<std::tuple<int,long int,float>>::Identical<std::tuple<int,long int,float>>);

        CIE_TEST_CHECK(!MatchTuple<std::tuple<int>>::Identical<std::tuple<unsigned int>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int>>::Identical<std::tuple<int,long int>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int,long int>>::Identical<std::tuple<int>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int,long int>>::Identical<std::tuple<long int,int>>);

        CIE_TEST_CHECK(MatchTuple<std::tuple<>>::Identical<std::tuple<>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<int>>::Identical<std::tuple<>>);
        CIE_TEST_CHECK(!MatchTuple<std::tuple<>>::Identical<std::tuple<int>>);
    }
}


} // namespace cie::ct