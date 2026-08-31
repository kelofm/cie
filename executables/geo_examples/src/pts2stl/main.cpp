// --- GEO Includes ---
#include "packages/triangulation/inc/meshgenerator.hpp"
#include "packages/io/inc/STLIO.hpp"

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <format>


using Scalar = double;


int main(int argc, const char** argv) {
    cie::utils::ArgParse parser("poly2stl");
    parser
        .addPositional("input-path")
        .addPositional("output-path")
        .addFlag({"-h", "--help"});
    cie::utils::ArgParse::Results arguments;

    try {
        arguments = parser.parseArguments(argc - 1, argv + 1);
    } catch (cie::Exception& rException) {
        parser.help(std::cerr);
        std::cerr << "\n" << rException.what() << std::endl;
        return 1;
    }

    try {
        const std::filesystem::path inputPath  = arguments.get<std::filesystem::path>("input-path");
        const std::filesystem::path outputPath = arguments.get<std::filesystem::path>("output-path");

        using Point = std::array<Scalar,2>;
        std::vector<Point> vertices;

        {
            std::ifstream file(inputPath);
            std::string line, component;

            for (std::size_t iLine=0ul; std::getline(file, line); ++iLine) {
                const auto itDelimiter = std::find(
                    line.begin(),
                    line.end(),
                    ',');
                CIE_CHECK(
                    itDelimiter != line.end(),
                    std::format(
                        "no delimiter found in line {} ({})",
                        iLine, line))

                const std::string xString(line.begin(), itDelimiter);
                const std::string yString(itDelimiter + 1, line.end());
                vertices.emplace_back();
                vertices.back().front() = std::stod(xString);
                vertices.back().back()  = std::stod(yString);
            } // while getline
        }

        const cie::geo::Triangulation triangulation = cie::geo::triangulate(
            vertices,
            cie::geo::TriangulationParameters());

        std::vector<float> output;
        output.reserve(triangulation.second.size() * 3 * 2);

        auto itOut = std::back_inserter(output);
        for (const auto triplet : triangulation.second) {
            for (const auto iVertex : triplet) {
                *itOut++ = triangulation.first[iVertex].front();
                *itOut++ = triangulation.first[iVertex].back();
            }
        }

        if (!vertices.empty()) {
            std::ofstream file(outputPath, std::ios::binary);
            cie::io::STLIO::Output<float,2> io(file);
            io.execute(output, {});
        } else {
            std::cerr << "empty input file\n";
            return 1;
        }

    } catch (cie::Exception& rException) {
        std::cerr << rException.what() << std::endl;
        return 1;
    }

    return 0;
}