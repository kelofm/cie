#ifndef CIE_GL_COLOR_HPP
#define CIE_GL_COLOR_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/logging/inc/OutputStream.hpp"

// --- STL Includes ---
#include <cstdint>
#include <initializer_list>
#include <ostream>


namespace cie {


/// @ingroup cieutils
class Color
{
public:
    using ValueType = std::uint8_t;

    using DataType = std::uint32_t;

public:
    Color();

    Color(ValueType v0, ValueType v1, ValueType v2, ValueType alpha);

    Color(ValueType v0, ValueType v1, ValueType v2);

    Color(const std::initializer_list<ValueType>& r_values);

    Color(Color&& r_rhs) = default;

    Color(const Color& r_rhs) = default;

    Color& operator=(Color&& r_rhs) = default;

    Color& operator=(const Color& r_rhs) = default;

    virtual ~Color() = default;

    std::string string() const;

    virtual std::string ANSI() const = 0;

    virtual std::string ANSIBackground() const = 0;

    virtual DataType hex() const = 0;

protected:
    template <Size ValueIndex>
    requires (0<=ValueIndex && ValueIndex<4)
    ValueType get() const;

    template <Size ValueIndex>
    requires (0<=ValueIndex && ValueIndex<4)
    void set(ValueType value);

private:
    void setValues(ValueType v0, ValueType v1, ValueType v2, ValueType v3 = 0xff);

private:
    DataType _data;
}; // class Color


std::ostream& operator<<(std::ostream& r_stream, const Color& r_color);

OutputStream& operator<<(OutputStream& r_stream, const Color& r_color);


class RGBAColor : public Color
{
public:
    using Color::ValueType;

    using Color::DataType;

public:
    using Color::Color;

    using Color::operator=;

    ValueType getRed() const;

    ValueType getGreen() const;

    ValueType getBlue() const;

    ValueType getAlpha() const;

    void setRed(ValueType value);

    void setGreen(ValueType value);

    void setBlue(ValueType value);

    void setAlpha(ValueType value);

    std::string ANSI() const override;

    std::string ANSIBackground() const override;

    DataType hex() const override;

public:
    /* --- Basic Colors --- */

    static const RGBAColor White;

    static const RGBAColor Black;

    /* --- TUM Palette --- */

    static const RGBAColor TUMBlack;

    static const RGBAColor TUMBlue;

    static const RGBAColor TUMWebBlue;

    static const RGBAColor TUMLightBlue;

    static const RGBAColor TUMLighterBlue;

    static const RGBAColor TUMLightestBlue;

    static const RGBAColor TUMOrange;

    static const RGBAColor TUMGreen;

    static const RGBAColor TUMGray;

    static const RGBAColor TUMLightGray;

    /* --- System Palette --- */

    static const RGBAColor Error;

    static const RGBAColor Warning;
}; // class RGBAColor


} // namespace cie

#endif