// --- Internal Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/trees/inc/abstree.hpp"

// --- STL Includes ---
#include <vector>
#include <memory>
#include <iostream>

namespace cie::utils {

struct TestTree : public AbsTree<TestTree,std::vector,std::shared_ptr<TestTree>>
{
    using Base = AbsTree<TestTree,std::vector,std::shared_ptr<TestTree>>;
    TestTree(Size level) : Base(level) {}
};


CIE_TEST_CASE( "AbsTree", "[trees]" )
{
    TestTree root(0);

    // Create function that splits nodes in two new ones
    // if a specified depth has not yet been reached
    Size targetDepth = 5;
    Size counter = 0;
    auto split = [targetDepth,&counter]( TestTree* node ) -> bool
    {
        if (node->level() < targetDepth)
        {
            ++counter;
            node->children().push_back( std::make_shared<TestTree>(node->level()+1) );
            node->children().push_back( std::make_shared<TestTree>(node->level()+1) );
            return true;
        }
        return false;
    };

    CIE_TEST_REQUIRE_NOTHROW( root.visit(split) );
    CIE_TEST_CHECK( counter == targetDepth );
}


}