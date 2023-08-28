// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- Internal Includes ---
#include "packages/io/inc/json.hpp"
#include "cmake_variables.hpp"

// --- STL Includes ---
#include <sstream>
#include <fstream>


namespace cie::io {


template <Size Index>
void checkJSONObject( JSONObject& r_json )
{
    // Existence
    CIE_TEST_CHECK( r_json.hasKey("bool") );
    CIE_TEST_CHECK( r_json.hasKey("string") );
    CIE_TEST_CHECK( r_json.hasKey("int") );
    CIE_TEST_CHECK( r_json.hasKey("float") );
    CIE_TEST_CHECK( r_json.hasKey("intArray") );
    CIE_TEST_CHECK( r_json.hasKey("stringArray") );
    CIE_TEST_CHECK( r_json.hasKey("mixedArray") );
    CIE_TEST_CHECK( r_json.hasKey("object") );

    // Type
    CIE_TEST_CHECK( r_json["bool"].is<bool>() );
    CIE_TEST_CHECK( r_json["string"].is<std::string>() );
    CIE_TEST_CHECK( r_json["int"].is<int>() );
    CIE_TEST_CHECK( r_json["float"].is<float>() );
    CIE_TEST_CHECK( r_json["intArray"].is<StaticArray<int,2>>() );
    CIE_TEST_CHECK( r_json["stringArray"].is<std::vector<std::string>>() );
    CIE_TEST_CHECK( r_json["mixedArray"].isArray() );
    CIE_TEST_CHECK( r_json["object"].isObject() );

    // False type
    CIE_TEST_CHECK( !r_json["bool"].is<int>() );
    CIE_TEST_CHECK( !r_json["string"].is<std::vector<float>>() );
    CIE_TEST_CHECK( !r_json["int"].is<std::string>() );
    CIE_TEST_CHECK( !r_json["intArray"].is<int>() );
    CIE_TEST_CHECK( !r_json["intArray"].is<std::vector<float>>() );
    CIE_TEST_CHECK( !r_json["stringArray"].is<std::string>() );
    CIE_TEST_CHECK( !r_json["stringArray"].is<std::vector<int>>() );
    CIE_TEST_CHECK( !r_json["mixedArray"].is<std::vector<float>>() );
    CIE_TEST_CHECK( !r_json["mixedArray"].is<std::vector<std::string>>() );
    CIE_TEST_CHECK( !r_json["object"].is<std::string>() );
    CIE_TEST_CHECK( !r_json["object"].isArray() );

    // Single value
    CIE_TEST_CHECK( r_json["bool"].as<bool>() );
    CIE_TEST_CHECK( r_json["string"].as<std::string>() == "std::string" );
    CIE_TEST_CHECK( r_json["int"].as<int>() == 12 );
    CIE_TEST_CHECK( r_json["float"].as<float>() == Approx(0.25) );
    CIE_TEST_CHECK( r_json["float"].as<double>() == Approx(0.25) );

    CIE_TEST_CHECK( r_json["mixedArray"][0].as<double>() == Approx(0.1) );
    CIE_TEST_CHECK( r_json["mixedArray"][1].as<std::string>() == "tenth" );

    // Array value
    std::vector<int>          intVector;
    StaticArray<int,2>         intArray;
    std::vector<float>        floatVector;
    StaticArray<float,2>       floatArray;
    std::vector<std::string>  stringVector;
    StaticArray<std::string,2> stringArray;

    CIE_TEST_CHECK_NOTHROW( intVector = r_json["intArray"].as<std::vector<int>>() );
    CIE_TEST_REQUIRE( intVector.size() == 2 );
    CIE_TEST_CHECK( intVector[0] == 1 );
    CIE_TEST_CHECK( intVector[1] == 0 );

    CIE_TEST_CHECK_NOTHROW( intArray = r_json["intArray"].as<StaticArray<int,2>>() );
    CIE_TEST_CHECK( intArray[0] == 1 );
    CIE_TEST_CHECK( intArray[1] == 0 );

    CIE_TEST_CHECK_NOTHROW( floatVector = r_json["intArray"].as<std::vector<float>>() );
    CIE_TEST_REQUIRE( floatVector.size() == 2 );
    CIE_TEST_CHECK( floatVector[0] == 1 );
    CIE_TEST_CHECK( floatVector[1] == 0 );

    CIE_TEST_CHECK_NOTHROW( floatArray = r_json["intArray"].as<StaticArray<float,2>>() );
    CIE_TEST_CHECK( floatArray[0] == 1 );
    CIE_TEST_CHECK( floatArray[1] == 0 );

    CIE_TEST_CHECK_NOTHROW( stringVector = r_json["stringArray"].as<std::vector<std::string>>() );
    CIE_TEST_REQUIRE( stringVector.size() == 2 );
    CIE_TEST_CHECK( stringVector[0] == "one" );
    CIE_TEST_CHECK( stringVector[1] == "zero" );

    CIE_TEST_CHECK_NOTHROW( stringArray = r_json["stringArray"].as<StaticArray<std::string,2>>() );
    CIE_TEST_CHECK( stringArray[0] == "one" );
    CIE_TEST_CHECK( stringArray[1] == "zero" );

    // New assignment
    JSONObject copy = r_json;

    CIE_TEST_CHECK_NOTHROW( copy.add("test", "testString") );
    CIE_TEST_CHECK( copy.hasKey("test") );
    CIE_TEST_CHECK( copy["test"].is<std::string>() );
    CIE_TEST_CHECK( copy["test"].as<std::string>() == "testString" );

    CIE_TEST_CHECK_NOTHROW( copy.add("not quite pi", std::vector<int> {3, 1, 4}) );
    CIE_TEST_CHECK( copy.hasKey("not quite pi") );
    CIE_TEST_CHECK( copy["not quite pi"].is<StaticArray<int,3>>() );
    StaticArray<int,3> intArray3;
    CIE_TEST_CHECK_NOTHROW( intArray3 = copy["not quite pi"].as<StaticArray<int,3>>() );
    CIE_TEST_CHECK( intArray3[0] == 3 );
    CIE_TEST_CHECK( intArray3[1] == 1 );
    CIE_TEST_CHECK( intArray3[2] == 4 );

    // Overwrite
    CIE_TEST_CHECK_NOTHROW( copy["test"].set( "another test" ) );
    CIE_TEST_CHECK( copy["test"].as<std::string>() == "another test" );

    // Attempted overwrite
    CIE_TEST_CHECK_THROWS( copy.add( "test", "invalid" ) );

    // Type mismatch
    CIE_TEST_CHECK_THROWS( copy["test"].set( std::vector<int> {0} ) );

    // Objects
    JSONObject object;
    CIE_TEST_CHECK_NOTHROW( object = copy["object"] );
    CIE_TEST_CHECK( object.hasKey("name") );
    CIE_TEST_CHECK( object.hasKey("power") );
    CIE_TEST_CHECK( object["name"].as<std::string>() == "mya-nee" );
    CIE_TEST_CHECK( object["power"].as<Size>() == 9001 );

    CIE_TEST_CHECK_NOTHROW( object = copy["object"].as<JSONObject>() );
    CIE_TEST_CHECK( object.hasKey("name") );
    CIE_TEST_CHECK( object.hasKey("power") );
    CIE_TEST_CHECK( object["name"].as<std::string>() == "mya-nee" );
    CIE_TEST_CHECK( object["power"].as<Size>() == 9001 );

    const std::string newContent = R"({"key" : "value"})";
    CIE_TEST_CHECK_NOTHROW( copy.add( "new object", JSONObject( newContent )) );
    CIE_TEST_CHECK( copy.hasKey( "new object" ) );
    CIE_TEST_CHECK( copy["new object"].is<JSONObject>() );
    CIE_TEST_CHECK( copy["new object"].hasKey( "key" ) );
    CIE_TEST_CHECK( copy["new object"]["key"].is<std::string>() );
    CIE_TEST_CHECK( copy["new object"]["key"].as<std::string>() == "value" );
}


CIE_TEST_CASE( "JSONObject", "[io]" )
{
    CIE_TEST_CASE_INIT( "JSONObject" )

    const std::string contents = R"({
        "bool"        : true,
        "string"      : "std::string",
        "int"         : 12,
        "float"       : 0.25,
        "intArray"    : [1, 0],
        "stringArray" : ["one", "zero"],
        "mixedArray"  : [0.1, "tenth"],
        "object" : {
            "name"  : "mya-nee",
            "power" : 9001
        }
    })";

    {
        JSONObject json;
        CIE_TEST_REQUIRE_NOTHROW( json = JSONObject(contents) );

        {
            CIE_TEST_CASE_INIT( "string constructor" )
            checkJSONObject<0>( json );
        } // string constructor

        {
            CIE_TEST_CASE_INIT( "string move constructor" )
            JSONObject local;
            CIE_TEST_REQUIRE_NOTHROW( local = JSONObject(std::string(contents)) );
            checkJSONObject<1>( local );
        } // string move constructor

        {
            CIE_TEST_CASE_INIT( "copy constructor" )
            JSONObject local(json);
            checkJSONObject<2>( local );
        } // copy constructor

        {
            CIE_TEST_CASE_INIT( "copy assignment operator" )
            JSONObject local;
            CIE_TEST_REQUIRE_NOTHROW( local = json );
            checkJSONObject<3>( local );
        } // copy assignment operator

        {
            CIE_TEST_CASE_INIT( "move constructor" )
            JSONObject local;
            CIE_TEST_REQUIRE_NOTHROW( local = JSONObject( JSONObject(contents) ) );
            checkJSONObject<4>( local );
        } // move constructor

        {
            CIE_TEST_CASE_INIT( "stream constructor" )

            std::stringstream stream;
            stream << contents;

            JSONObject local;
            CIE_TEST_REQUIRE_NOTHROW( local = JSONObject(stream) );
            checkJSONObject<5>( json );
        } // stream constructor

        {
            CIE_TEST_CASE_INIT( "construction from file" )

            const std::filesystem::path filePath = getOutputPath()/ "JSONObject_test.json";

            std::ofstream outputFile( filePath );
            outputFile << contents;
            outputFile.close();

            JSONObject local;
            CIE_TEST_REQUIRE_NOTHROW( local = JSONObject(filePath) );
            checkJSONObject<6>( local );
        } // construction from file
    }



}


} // namespace cie::io