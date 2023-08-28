#ifndef CIE_UTILS_MACROS_TYPEDEFS_HPP
#define CIE_UTILS_MACROS_TYPEDEFS_HPP

// --- STL Includes ---
#include <memory>


/// Define smart pointers to the class as member typedefs.
#define CIE_DEFINE_CLASS_POINTERS(CLASS_NAME)                       \
    using SharedPointer = std::shared_ptr<CLASS_NAME>;              \
    using UniquePointer = std::unique_ptr<CLASS_NAME>;              \
    using WeakPointer   = std::weak_ptr<CLASS_NAME>;                \
    using SharedConstPointer = std::shared_ptr<const CLASS_NAME>;   \
    using UniqueConstPointer = std::unique_ptr<const CLASS_NAME>;   \
    using WeakConstPointer = std::weak_ptr<const CLASS_NAME>;


/// Explicitly declare the class' default constructor as default.
#define CIE_DEFINE_CLASS_DEFAULT_CONSTRUCTOR(CLASS_NAME) \
    /** Default default constructor.*/                   \
    CLASS_NAME() = default;


/** Explicitly declare the class' move constructor
 *  and move assignment operator as default.
 */
#define CIE_DEFINE_CLASS_DEFAULT_MOVES(CLASS_NAME)                  \
    /** Default move constructor.*/                                 \
    CLASS_NAME(CLASS_NAME&& r_rhs) noexcept = default;              \
    /** Default move assignment operator.*/                         \
    CLASS_NAME& operator=(CLASS_NAME&& r_rhs) noexcept = default;


/** Explicitly declare the class' copy constructor
 *  and copy assignment operator as default.
 */
#define CIE_DEFINE_CLASS_DEFAULT_COPIES(CLASS_NAME)           \
    /** Default copy constructor.*/                           \
    CLASS_NAME(const CLASS_NAME& r_rhs) = default;            \
    /** Default copy assignment operator.*/                   \
    CLASS_NAME& operator=(const CLASS_NAME& r_rhs) = default;


/** Explicitly declare the class' default constructor,
 *  move constructor, move assignment operator,
 *  copy constructor, and copy assignment operator
 *  as default.
 */
#define CIE_DEFINE_CLASS_DEFAULTS(CLASS_NAME)        \
    CIE_DEFINE_CLASS_DEFAULT_CONSTRUCTOR(CLASS_NAME) \
    CIE_DEFINE_CLASS_DEFAULT_MOVES(CLASS_NAME)       \
    CIE_DEFINE_CLASS_DEFAULT_COPIES(CLASS_NAME)



#endif
