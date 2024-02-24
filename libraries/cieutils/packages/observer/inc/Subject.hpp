#ifndef CIEUTILS_Subject_HPP
#define CIEUTILS_Subject_HPP

// --- Internal Includes ---
#include "packages/observer/inc/Observer.hpp"

// --- STL Includes ---
#include <memory>
#include <vector>


namespace cie::utils::observer {


/// @addtogroup cieutils
class Subject : public std::enable_shared_from_this<Subject>
{
public:
    Subject() = default;
    Subject( const SubjectPtr& copy );
    SubjectPtr operator=( const SubjectPtr& copy );
    virtual ~Subject();

    void notifyObservers();
    int attachObserver( ObserverPtr observer );
    void detachObserver( ObserverPtr observer );
    int observerID( ObserverPtr observer );

    const std::vector<ObserverPtr>& observers() const;

    template <class DerivedSubject>
    std::shared_ptr<DerivedSubject> as();

protected:
    std::vector<ObserverPtr>     _observers;
};


} // namespace cie::utils::observer

#include "packages/observer/impl/Subject_impl.hpp"

#endif
