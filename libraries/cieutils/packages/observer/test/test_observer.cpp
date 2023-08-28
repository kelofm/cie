// --- Internal Imports ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/observer/inc/Observer.hpp"
#include "packages/observer/inc/Subject.hpp"

// --- STD Imports ---
#include <iostream>
#include <typeinfo>

namespace cie::utils::observer {


// Observer wrapper
struct TestObserver : public Observer
{
    TestObserver(int* counter) : Observer(), _counter(counter) {++(*_counter);}
    ~TestObserver() {--(*_counter);}
    void onSubjectChange( ) override {}

    int* _counter;
};
using TestObserverPtr = std::shared_ptr<TestObserver>;


// Subject wrapper
struct TestSubject : public Subject
{
    TestSubject(int* counter) : Subject(), _counter(counter) {++(*_counter);}
    ~TestSubject() {--(*_counter);}

    int* _counter;
};
using TestSubjectPtr = std::shared_ptr<TestSubject>;



CIE_TEST_CASE( "Observer", "[observer]" )
{
    CIE_TEST_CASE_INIT( "Observer" )

    // Create counters
    int numberOfSubjects    = 0;
    int numberOfObservers   = 0;

    {
        // Create p_subject
        auto p_subject = SubjectPtr(
            new TestSubject(&numberOfSubjects)
        );

        // Create Observers
        auto p_observer1 = ObserverPtr(
            new TestObserver(&numberOfObservers)
        );

        auto p_observer2 = ObserverPtr(
            new TestObserver(&numberOfObservers)
        );

        // Check initial state of wrappers
        CIE_TEST_REQUIRE( numberOfSubjects   == 1 );
        CIE_TEST_REQUIRE( numberOfObservers  == 2 );

        CIE_TEST_CHECK( p_observer1->subject().get() == nullptr );
        CIE_TEST_CHECK_NOTHROW( p_observer1->onSubjectChange( ) );
        
        CIE_TEST_CHECK_NOTHROW( p_observer1->setSubject(p_subject) );
        CIE_TEST_CHECK( p_observer1->subject().get() == p_subject.get() );
        CIE_TEST_CHECK_NOTHROW( p_observer1->onSubjectChange( ) );

        // Check attaching a new observer
        int observerID = -1;
        CIE_TEST_CHECK_NOTHROW( observerID = p_subject->attachObserver( p_observer1 ) );
        CIE_TEST_CHECK( observerID == 0 );
        CIE_TEST_CHECK_NOTHROW( observerID = p_subject->attachObserver( p_observer2 ) );
        CIE_TEST_CHECK( observerID == 1 );

        // Check attaching an existing observer
        CIE_TEST_CHECK_NOTHROW( observerID = p_subject->attachObserver( p_observer1 ) );
        CIE_TEST_CHECK_NOTHROW( observerID = p_subject->observerID( p_observer1 ) );
        CIE_TEST_CHECK( observerID == 0 );

        CIE_TEST_CHECK_NOTHROW( observerID = p_subject->attachObserver( p_observer2 ) );
        CIE_TEST_CHECK_NOTHROW( observerID = p_subject->observerID( p_observer2 ) );
        CIE_TEST_CHECK( observerID == 1 );

        // Check detaching existing observer
        CIE_TEST_CHECK_NOTHROW( p_observer1->detach() );
        CIE_TEST_CHECK( p_subject->observerID( p_observer1 ) == -1 );
        CIE_TEST_CHECK( p_subject->observerID( p_observer2 ) == 0 );

        // Check detaching non-existent observer
        CIE_TEST_CHECK_NOTHROW( p_observer1->detach() );
        CIE_TEST_CHECK( p_subject->observerID( p_observer1 ) == -1 );
        CIE_TEST_CHECK( p_subject->observerID( p_observer2 ) == 0 );

        // Check detaching observer from p_subject
        CIE_TEST_CHECK_NOTHROW( p_subject->detachObserver( p_observer2 ) );
        CIE_TEST_CHECK( p_subject->observerID( p_observer2 ) == -1 );
        CIE_TEST_CHECK( p_subject->observerID( p_observer1 ) == -1 );

        // Check casting
        auto p_castSubject  = p_subject->as<TestSubject>();
        CIE_TEST_CHECK( typeid(p_subject) == typeid(SubjectPtr) );
        CIE_TEST_CHECK( typeid(p_castSubject) == typeid(TestSubjectPtr) );

        auto p_castObserver = p_observer1->as<TestObserver>();
        CIE_TEST_CHECK( typeid(p_observer1) == typeid(ObserverPtr) );
        CIE_TEST_CHECK( typeid(p_castObserver) == typeid(TestObserverPtr) );
    }

    // Check final state of wrappers
    CIE_TEST_REQUIRE( numberOfSubjects   == 0 );
    CIE_TEST_REQUIRE( numberOfObservers  == 0 );
}


} // namespace cie::utils