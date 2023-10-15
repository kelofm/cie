#ifndef CIEUTILS_Observer_HPP
#define CIEUTILS_Observer_HPP

// --- STL Imports ---
#include <memory>


namespace cie::utils::observer {

class Observer;
class Subject;
using ObserverPtr    = std::shared_ptr<Observer>;
using SubjectPtr     = std::shared_ptr<Subject>;


/// @addtogroup cieutils
class Observer : public std::enable_shared_from_this<Observer>
{
public:
    Observer();
    Observer( SubjectPtr subject );
    Observer( const ObserverPtr& copy );
    ObserverPtr operator=( const ObserverPtr& copy );
    virtual ~Observer();

    virtual void onSubjectChange() = 0;
    void detach();

    template <class DerivedObserver>
    std::shared_ptr<DerivedObserver> as();

    void setSubject( SubjectPtr subject );
    SubjectPtr subject();
    const SubjectPtr subject() const;

protected:
    std::weak_ptr<Subject>   _subject;
};


} // namespace cie::utils::observer

#include "packages/observer/impl/Observer_impl.hpp"

#endif