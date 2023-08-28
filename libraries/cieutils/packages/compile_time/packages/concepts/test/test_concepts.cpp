// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <vector>
#include <array>
#include <deque>
#include <map>
#include <unordered_map>
#include <utility>

namespace cie {


template <concepts::Numeric T>
void numericTest( T value )
{}

template <concepts::RawPointer T>
void pointerTest( T ptr )
{}

template <class ValueType, class PtrType>
void classPtrTest( PtrType ptr )
requires concepts::RawPointer<PtrType,ValueType>
{}

template <concepts::NonRawPointer T>
void nonPtrTest( T value )
{}

template <concepts::Iterator IteratorType>
void iteratorTest( IteratorType it )
{}

template <class ValueType, class IteratorType>
void classIteratorTest( IteratorType it )
requires concepts::Iterator<IteratorType, ValueType>
{}

template <concepts::Container ContainerType>
void containerTest( ContainerType container )
{}

template <class ValueType, class ContainerType>
void classContainerTest( ContainerType container )
requires concepts::Container<ContainerType,ValueType>
{}

template <concepts::PointerContainer ContainerType>
void pointerContainerTest( ContainerType container )
{}

template <concepts::NonPointerContainer ContainerType>
void nonPtrContainerTest( ContainerType container )
{}

template <concepts::IteratorContainer ContainerType>
void iteratorContainerTest( ContainerType container )
{}




CIE_TEST_CASE( "Numeric Concepts", "[concepts]" )
{
    int     a = 0;
    Size    b = 0;
    float   c = 0;
    double  d = 0;

    numericTest(a);
    numericTest(b);
    numericTest(c);
    numericTest(d);
}


CIE_TEST_CASE( "Ptr Concepts", "[concepts]" )
{
    int     a = 0;
    void*   b = nullptr;

    pointerTest( &a );
    pointerTest( b );
    classPtrTest<int>( &a );
    nonPtrTest( a );
}


CIE_TEST_CASE( "Iterator Concepts", "[concepts]" )
{
    std::array<int,2>           a;
    std::vector<int>            b;
    std::deque<int>             c;
    std::map<int,int>           d;
    std::unordered_map<int,int> e;

    iteratorTest(a.begin());
    iteratorTest(b.begin());
    iteratorTest(c.begin());
    iteratorTest(d.begin());
    iteratorTest(e.begin());

    classIteratorTest<int>(a.begin());
    classIteratorTest<int>(b.begin());
    classIteratorTest<int>(c.begin());
    classIteratorTest<std::pair<const int,int>>(d.begin());
    classIteratorTest<std::pair<const int,int>>(e.begin());
}


CIE_TEST_CASE( "Container Concepts", "[concepts]" )
{
    std::array<int,2>                       a{};
    std::vector<int>                        b{};
    std::deque<int>                         c{};
    std::map<int,int>                       d{};
    std::unordered_map<int,int>             e{};
    std::deque<int*>                        f{};
    std::deque<std::deque<int>::iterator>   g{};

    containerTest( a );
    containerTest( b );
    containerTest( c );
    containerTest( d );
    containerTest( e );
    containerTest( f );
    containerTest( g );

    classContainerTest<int>( a );
    classContainerTest<int>( b );
    classContainerTest<int>( c );
    classContainerTest<std::pair<const int,int>>( d );
    classContainerTest<std::pair<const int,int>>( e );
    classContainerTest<int*>( f );
    classContainerTest<std::deque<int>::iterator>( g );

    pointerContainerTest( f );
    iteratorContainerTest( g );

    nonPtrContainerTest( a );
    nonPtrContainerTest( b );
    nonPtrContainerTest( c );
    nonPtrContainerTest( d );
    nonPtrContainerTest( e );
}


} // namespace cie