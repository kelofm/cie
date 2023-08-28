// ---- Utility Includes ---
#include "packages/exceptions/inc/exception.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- Internal Includes ---
#include "packages/types/inc/Color.hpp"


namespace cie {


template <Size ByteIndex>
constexpr void copyByte(const Color::DataType source, Color::DataType& r_destination)
{
    const Color::DataType shift = 8 * ByteIndex;
    const Color::DataType byteMask = 0xff << shift;
    const Color::DataType inverseByteMask = ~byteMask;
    r_destination &= inverseByteMask; // clear the target byte
    r_destination |= source & byteMask; // copy the target byte from source
}


Color::Color() : Color(0, 0, 0, 0xff)
{
}


Color::Color(ValueType v0, ValueType v1, ValueType v2, ValueType alpha)
{
    this->setValues(v0, v1, v2, alpha);
}


Color::Color(ValueType v0, ValueType v1, ValueType v2)
    : Color(v0, v1, v2, 0xff)
{
}


Color::Color(const std::initializer_list<ValueType>& r_values)
{
    const auto begin = r_values.begin();
    switch (r_values.size())
    {
        case 3:
        {
            this->setValues(*begin, *(begin + 1), *(begin + 2));
            break;
        }
        case 4:
        {
            this->setValues(*begin, *(begin + 1), *(begin + 2), *(begin + 3));
            break;
        }
        default:
        {
            CIE_THROW(
                OutOfRangeException,
                "Color takes 3 or 4 values, but " + std::to_string(r_values.size()) + " were given"
            )
            break;
        }
    }
}


void Color::setValues(ValueType v0, ValueType v1, ValueType v2, ValueType v3)
{
    DataType tmp;
    tmp = DataType(v0); copyByte<0>(tmp, _data);
    tmp = DataType(v1) << 8; copyByte<1>(tmp, _data);
    tmp = DataType(v2) << 16; copyByte<2>(tmp, _data);
    tmp = DataType(v3) << 24; copyByte<3>(tmp, _data);
}


template <Size ValueIndex>
requires (0<=ValueIndex && ValueIndex<4)
Color::ValueType Color::get() const
{
    constexpr Color::ValueType shift = 8 * ValueIndex;
    Color::DataType output = 0;
    copyByte<ValueIndex>(_data, output);
    return Color::ValueType(output >> shift);
}


template <Size ValueIndex>
requires (0<=ValueIndex && ValueIndex<4)
void Color::set(Color::ValueType value)
{
    constexpr Color::ValueType shift = 8 * ValueIndex;
    Color::DataType offsetValue = Color::DataType(value) << shift;
    copyByte<ValueIndex>(offsetValue, _data);
}


template Color::ValueType Color::get<0>() const;


template Color::ValueType Color::get<1>() const;


template Color::ValueType Color::get<2>() const;


template Color::ValueType Color::get<3>() const;


template void Color::set<0>(Color::ValueType value);


template void Color::set<1>(Color::ValueType value);


template void Color::set<2>(Color::ValueType value);


template void Color::set<3>(Color::ValueType value);


std::string Color::string() const
{
    return std::to_string(static_cast<unsigned int>(this->get<0>())) + ";"
           + std::to_string(static_cast<unsigned int>(this->get<1>())) + ";"
           + std::to_string(static_cast<unsigned int>(this->get<2>())) + ";"
           + std::to_string(static_cast<unsigned int>(this->get<3>()));
}


std::ostream& operator<<(std::ostream& r_stream, const Color& r_color)
{
    return r_stream << r_color.ANSI();
}


OutputStream& operator<<(OutputStream& r_stream, const Color& r_color)
{
    return r_stream << r_color.ANSI();
}


RGBAColor::ValueType RGBAColor::getRed() const
{
    return this->get<0>();
}


RGBAColor::ValueType RGBAColor::getGreen() const
{
    return this->get<1>();
}


RGBAColor::ValueType RGBAColor::getBlue() const
{
    return this->get<2>();
}


RGBAColor::ValueType RGBAColor::getAlpha() const
{
    return this->get<3>();
}


void RGBAColor::setRed(ValueType value)
{
    this->set<0>(value);
}


void RGBAColor::setGreen(ValueType value)
{
    this->set<1>(value);
}


void RGBAColor::setBlue(ValueType value)
{
    this->set<2>(value);
}


void RGBAColor::setAlpha(ValueType value)
{
    this->set<3>(value);
}


std::string RGBAColor::ANSI() const
{
    std::string ansi;
    ansi.reserve(19);
    ansi += "\x1b[38;2;";
    ansi += std::to_string(static_cast<unsigned int>(this->getRed()));
    ansi += ';';
    ansi += std::to_string(static_cast<unsigned int>(this->getGreen()));
    ansi += ';';
    ansi += std::to_string(static_cast<unsigned int>(this->getBlue()));
    ansi += 'm';
    return ansi;
}


std::string RGBAColor::ANSIBackground() const
{
    std::string ansi;
    ansi.reserve(19);
    ansi += "\x1b[48;2;";
    ansi += std::to_string(static_cast<unsigned int>(this->getRed()));
    ansi += ';';
    ansi += std::to_string(static_cast<unsigned int>(this->getGreen()));
    ansi += ';';
    ansi += std::to_string(static_cast<unsigned int>(this->getBlue()));
    ansi += 'm';
    return ansi;
}


RGBAColor::DataType RGBAColor::hex() const
{
    DataType output = 0x0;
    copyByte<2>(DataType(this->getRed()) << 16, output); // red
    copyByte<1>(DataType(this->getGreen()) << 8, output); // green
    copyByte<0>(DataType(this->getBlue()), output); // blue
    return output;
}


const RGBAColor RGBAColor::White(0xff, 0xff, 0xff);


const RGBAColor RGBAColor::Black(0x0, 0x0, 0x0);


const RGBAColor RGBAColor::TUMBlack(51, 51, 51);


const RGBAColor RGBAColor::TUMBlue(0, 82, 147); // #005293


const RGBAColor RGBAColor::TUMWebBlue(48, 112, 179);


const RGBAColor RGBAColor::TUMLightBlue(0, 101, 189);


const RGBAColor RGBAColor::TUMLighterBlue(100, 160, 200);


const RGBAColor RGBAColor::TUMLightestBlue(152, 198, 234);


const RGBAColor RGBAColor::TUMOrange(227, 114, 34); // #e07222


const RGBAColor RGBAColor::TUMGreen(162, 173, 0); // #a1ad00


const RGBAColor RGBAColor::TUMLightGray(218, 215, 203);


const RGBAColor RGBAColor::Error(0xd8, 0x0, 0x0c);


const RGBAColor RGBAColor::Warning(0xff, 0xcc, 0xcc);


} // namespace cie