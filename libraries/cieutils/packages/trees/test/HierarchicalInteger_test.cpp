// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/trees/inc/HierarchicalInteger.hpp"


namespace cie {


CIE_TEST_CASE("HierarchicalInteger", "[trees][!mayfail]")
{
    CIE_TEST_CASE_INIT("HierarchicalInteger")

    {
        CIE_TEST_CASE_INIT("2D")
        const unsigned short Dimension = 2;
        using HI = HierarchicalInteger<Dimension, unsigned>;

        CIE_TEST_CHECK_NOTHROW(HI(0));
        CIE_TEST_CHECK_NOTHROW(HI(1));
        CIE_TEST_CHECK_NOTHROW(HI(std::numeric_limits<HI::Value>::max()));

        {
            CIE_TEST_CASE_INIT("getCell")
            CIE_TEST_CHECK(HI::getCell(0) == 0);
            CIE_TEST_CHECK(HI::getCell(1) == 0);
            CIE_TEST_CHECK(HI::getCell(2) == 0);
            CIE_TEST_CHECK(HI::getCell(3) == 0);
            CIE_TEST_CHECK(HI::getCell(4) == 1);
            CIE_TEST_CHECK(HI::getCell(5) == 1);
            CIE_TEST_CHECK(HI::getCell(6) == 1);
            CIE_TEST_CHECK(HI::getCell(7) == 1);
            CIE_TEST_CHECK(HI::getCell(8) == 2);
            CIE_TEST_CHECK(HI::getCell(9) == 2);
            CIE_TEST_CHECK(HI::getCell(10) == 2);
            CIE_TEST_CHECK(HI::getCell(11) == 2);
            CIE_TEST_CHECK(HI::getCell(12) == 3);
            CIE_TEST_CHECK(HI::getCell(13) == 3);
            CIE_TEST_CHECK(HI::getCell(14) == 3);
            CIE_TEST_CHECK(HI::getCell(15) == 3);
            CIE_TEST_CHECK(HI::getCell(16) == 4);
            CIE_TEST_CHECK(HI::getCell(17) == 4);
            CIE_TEST_CHECK(HI::getCell(18) == 4);
            CIE_TEST_CHECK(HI::getCell(19) == 4);
            CIE_TEST_CHECK(HI::getCell(20) == 5);
        } // getCell

        {
            CIE_TEST_CASE_INIT("getParentCell")
            CIE_TEST_CHECK(HI::getParentCell(0) == 0);
            CIE_TEST_CHECK(HI::getParentCell(1) == 0);
            CIE_TEST_CHECK(HI::getParentCell(2) == 0);
            CIE_TEST_CHECK(HI::getParentCell(3) == 0);
            CIE_TEST_CHECK(HI::getParentCell(4) == 0);
            CIE_TEST_CHECK(HI::getParentCell(5) == 1);
            CIE_TEST_CHECK(HI::getParentCell(6) == 1);
            CIE_TEST_CHECK(HI::getParentCell(7) == 1);
            CIE_TEST_CHECK(HI::getParentCell(8) == 1);
            CIE_TEST_CHECK(HI::getParentCell(9) == 2);
            CIE_TEST_CHECK(HI::getParentCell(10) == 2);
            CIE_TEST_CHECK(HI::getParentCell(11) == 2);
            CIE_TEST_CHECK(HI::getParentCell(12) == 2);
            CIE_TEST_CHECK(HI::getParentCell(13) == 3);
            CIE_TEST_CHECK(HI::getParentCell(14) == 3);
            CIE_TEST_CHECK(HI::getParentCell(15) == 3);
            CIE_TEST_CHECK(HI::getParentCell(16) == 3);
            CIE_TEST_CHECK(HI::getParentCell(17) == 4);
            CIE_TEST_CHECK(HI::getParentCell(18) == 4);
            CIE_TEST_CHECK(HI::getParentCell(19) == 4);
            CIE_TEST_CHECK(HI::getParentCell(20) == 4);
            CIE_TEST_CHECK(HI::getParentCell(21) == 5);
        } // getParentCell

        {
            CIE_TEST_CASE_INIT("getLevel")
            CIE_TEST_CHECK(HI::getLevel(0) == 0);
            CIE_TEST_CHECK(HI::getLevel(1) == 0);
            CIE_TEST_CHECK(HI::getLevel(2) == 0);
            CIE_TEST_CHECK(HI::getLevel(3) == 0);
            CIE_TEST_CHECK(HI::getLevel(4) == 1);
            CIE_TEST_CHECK(HI::getLevel(5) == 1);
            CIE_TEST_CHECK(HI::getLevel(6) == 1);
            CIE_TEST_CHECK(HI::getLevel(7) == 1);
            CIE_TEST_CHECK(HI::getLevel(8) == 1);
            CIE_TEST_CHECK(HI::getLevel(9) == 1);
            CIE_TEST_CHECK(HI::getLevel(10) == 1);
            CIE_TEST_CHECK(HI::getLevel(11) == 1);
            CIE_TEST_CHECK(HI::getLevel(12) == 1);
            CIE_TEST_CHECK(HI::getLevel(13) == 1);
            CIE_TEST_CHECK(HI::getLevel(14) == 1);
            CIE_TEST_CHECK(HI::getLevel(15) == 1);
            CIE_TEST_CHECK(HI::getLevel(16) == 1);
            CIE_TEST_CHECK(HI::getLevel(17) == 1);
            CIE_TEST_CHECK(HI::getLevel(18) == 1);
            CIE_TEST_CHECK(HI::getLevel(19) == 1);
            CIE_TEST_CHECK(HI::getLevel(20) == 2);
        } // getLevel

        {
            CIE_TEST_CASE_INIT("getLocalNode")
            CIE_TEST_CHECK(HI::getLocalNode(0) == 0);
            CIE_TEST_CHECK(HI::getLocalNode(1) == 1);
            CIE_TEST_CHECK(HI::getLocalNode(2) == 2);
            CIE_TEST_CHECK(HI::getLocalNode(3) == 3);
            CIE_TEST_CHECK(HI::getLocalNode(4) == 0);
            CIE_TEST_CHECK(HI::getLocalNode(5) == 1);
            CIE_TEST_CHECK(HI::getLocalNode(6) == 2);
            CIE_TEST_CHECK(HI::getLocalNode(7) == 3);
            CIE_TEST_CHECK(HI::getLocalNode(8) == 0);
            CIE_TEST_CHECK(HI::getLocalNode(9) == 1);
            CIE_TEST_CHECK(HI::getLocalNode(10) == 2);
            CIE_TEST_CHECK(HI::getLocalNode(11) == 3);
            CIE_TEST_CHECK(HI::getLocalNode(12) == 0);
            CIE_TEST_CHECK(HI::getLocalNode(13) == 1);
            CIE_TEST_CHECK(HI::getLocalNode(14) == 2);
            CIE_TEST_CHECK(HI::getLocalNode(15) == 3);
            CIE_TEST_CHECK(HI::getLocalNode(16) == 0);
            CIE_TEST_CHECK(HI::getLocalNode(17) == 1);
            CIE_TEST_CHECK(HI::getLocalNode(18) == 2);
            CIE_TEST_CHECK(HI::getLocalNode(19) == 3);
            CIE_TEST_CHECK(HI::getLocalNode(20) == 0);

            CIE_TEST_CHECK(HI::getLocalNode(72) == 0);
            CIE_TEST_CHECK(HI::getLocalNode(73) == 1);
            CIE_TEST_CHECK(HI::getLocalNode(74) == 2);
            CIE_TEST_CHECK(HI::getLocalNode(75) == 3);
        } // getLocalNode

        {
            CIE_TEST_CASE_INIT("getLocalCell")
            CIE_TEST_CHECK(HI::getLocalCell(0) == 0);
            CIE_TEST_CHECK(HI::getLocalCell(1) == 0);
            CIE_TEST_CHECK(HI::getLocalCell(2) == 1);
            CIE_TEST_CHECK(HI::getLocalCell(3) == 2);
            CIE_TEST_CHECK(HI::getLocalCell(4) == 3);
            CIE_TEST_CHECK(HI::getLocalCell(5) == 0);
            CIE_TEST_CHECK(HI::getLocalCell(6) == 1);
            CIE_TEST_CHECK(HI::getLocalCell(7) == 2);
            CIE_TEST_CHECK(HI::getLocalCell(8) == 3);
            CIE_TEST_CHECK(HI::getLocalCell(9) == 0);
            CIE_TEST_CHECK(HI::getLocalCell(10) == 1);
            CIE_TEST_CHECK(HI::getLocalCell(11) == 2);
            CIE_TEST_CHECK(HI::getLocalCell(12) == 3);
            CIE_TEST_CHECK(HI::getLocalCell(13) == 0);
            CIE_TEST_CHECK(HI::getLocalCell(14) == 1);
            CIE_TEST_CHECK(HI::getLocalCell(15) == 2);
            CIE_TEST_CHECK(HI::getLocalCell(16) == 3);
            CIE_TEST_CHECK(HI::getLocalCell(17) == 0);
            CIE_TEST_CHECK(HI::getLocalCell(18) == 1);
            CIE_TEST_CHECK(HI::getLocalCell(19) == 2);
            CIE_TEST_CHECK(HI::getLocalCell(20) == 3);
            CIE_TEST_CHECK(HI::getLocalCell(21) == 0);
        } // getLocalCell

        {
            CIE_TEST_CASE_INIT("getLocal")

            CIE_TEST_CHECK(HI::getLocal(0u, 0) == 0);
            CIE_TEST_CHECK(HI::getLocal(2u, 1) == 0b10);
            CIE_TEST_CHECK(HI::getLocal(16u, 2) == 0b1100);
            CIE_TEST_CHECK(HI::getLocal(65u, 3) == 0b101101);
        }

        {
            CIE_TEST_CASE_INIT("decay")

            // Cell 0
            CIE_TEST_CHECK(HI(0).value() == 0);
            CIE_TEST_CHECK(HI(1).value() == 1);
            CIE_TEST_CHECK(HI(2).value() == 2);
            CIE_TEST_CHECK(HI(3).value() == 3);

            // Cell 1
            CIE_TEST_CHECK(HI(4).value() == 0);
            CIE_TEST_CHECK(HI(5).value() == 5);
            CIE_TEST_CHECK(HI(6).value() == 6);
            CIE_TEST_CHECK(HI(7).value() == 7);

            // Cell 2
            CIE_TEST_CHECK(HI(8).value() == 5);
            CIE_TEST_CHECK(HI(9).value() == 1);
            CIE_TEST_CHECK(HI(10).value() == 10);
            CIE_TEST_CHECK(HI(11).value() == 11);

            // Cell 3
            CIE_TEST_CHECK(HI(12).value() == 6);
            CIE_TEST_CHECK(HI(13).value() == 7);
            CIE_TEST_CHECK(HI(14).value() == 2);
            CIE_TEST_CHECK(HI(15).value() == 15);

            // Cell 4
            CIE_TEST_CHECK(HI(16).value() == 7);
            CIE_TEST_CHECK(HI(17).value() == 11);
            CIE_TEST_CHECK(HI(18).value() == 15);
            CIE_TEST_CHECK(HI(19).value() == 3);

            // Cell 14
            CIE_TEST_CHECK(HI(56).value() == 31);
            CIE_TEST_CHECK(HI(57).value() == 7);
            CIE_TEST_CHECK(HI(58).value() == 55);
            CIE_TEST_CHECK(HI(59).value() == 59);

            // Cell 17
            HI(68);
            CIE_TEST_CHECK(HI(68).value() == 7);
            CIE_TEST_CHECK(HI(69).value() == 47);
            CIE_TEST_CHECK(HI(70).value() == 59);
            CIE_TEST_CHECK(HI(71).value() == 71);

            // Cell 19
            CIE_TEST_CHECK(HI(80).value() == 71);
            CIE_TEST_CHECK(HI(81).value() == 75);
            CIE_TEST_CHECK(HI(82).value() == 79);
            CIE_TEST_CHECK(HI(83).value() == 3);
        } // decay
    }
}


} // namespace cie
