// --- GL Includes ---
#include "packages/utility/inc/paths.hpp"
#include "packages/scene/inc/Scene.hpp"
#include "packages/context/inc/GLFWContext.hpp"
#include "packages/camera/inc/Camera.hpp"
#include "packages/camera/inc/OrthographicProjection.hpp"

// --- CSG Includes ---
#include "packages/primitives/inc/Cube.hpp"
#include "packages/trees/inc/Cell.hpp"
#include "packages/trees/inc/SpaceTreeNode.hpp"
#include "packages/trees/inc/CartesianGridSampler.hpp"
#include "packages/trees/inc/MidPointSplitPolicy.hpp"
#include "packages/trees/inc/WeightedSplitPolicy.hpp"
#include "packages/trees/inc/ContiguousSpaceTree.hpp"

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"
#include "packages/concurrency/inc/ResourceQueue.hpp"
#include "packages/macros/inc/logging.hpp"

// --- STL Includes ---
#include <cstdlib>
#include <chrono>
#include <filesystem>
#include <set>
#include <variant>
#include <iostream>


namespace cie {


// -------------------------------------------------------------------------
// TYPE ALIASES
// -------------------------------------------------------------------------

using CoordinateType = float;
using PointType      = geo::Traits<2,CoordinateType>::Point;

// Generic tree
using Primitive      = geo::Box<2,CoordinateType>;
using Cell           = geo::Cell<Primitive>;
using Node           = geo::SpaceTreeNode<Cell,int>;
using NodePtr        = Node::SharedPointer;

using Sampler        = geo::CartesianGridSampler<Primitive>;
using Splitter       = geo::WeightedSplitPolicy<Node::SamplePointIterator,Node::value_iterator>;

// Contiguous tree
using ContiguousTree = geo::ContiguousSpaceTree<geo::Cube<2,CoordinateType>,Size,tags::Serial>;
using ContiguousNode = ContiguousTree::Node;


// -------------------------------------------------------------------------
// CLASS DECLARATION
// -------------------------------------------------------------------------


template <class TNode> requires false
void makeRoot(const Size samplingOrder) {}


template <class TNode> requires std::is_same_v<TNode,Node>
NodePtr makeRoot(const Size samplingOrder)
{
    Sampler::Resolution resolution;
    std::fill(resolution.begin(), resolution.end(), samplingOrder);

    if constexpr (concepts::Cube<typename Node::Cell::Primitive>) {
        return NodePtr(new Node(
            Node::sampler_ptr(new Sampler(resolution)),
            Node::split_policy_ptr(new Splitter),
            0,
            PointType {0.0, 0.0},
            1.0
        ));
    } else if constexpr (concepts::Box<typename Node::Cell::Primitive>) {
        return NodePtr(new Node(
            Node::sampler_ptr(new Sampler(resolution)),
            Node::split_policy_ptr(new Splitter),
            0,
            PointType {0.0, 0.0},
            PointType {1.0, 1.0}
        ));
    } else {
        CIE_THROW(Exception, "Unsupported primitive!")
    }
}


template <class TTree> requires std::is_same_v<TTree,ContiguousTree>
std::unique_ptr<ContiguousTree> makeRoot(const Size samplingOrder)
{
    if constexpr (concepts::Cube<typename TTree::Geometry>) {
        return std::make_unique<ContiguousTree>(PointType {0.0, 0.0}, 1.0);
    } else if constexpr (concepts::Box<typename TTree::Geometry>) {
        return std::make_unique<ContiguousTree>(PointType {0.0, 0.0}, ContiguousTree {1.0, 1.0});
    } else {
        CIE_THROW(Exception, "Unsupported primitive!")
    }
}


template <class TPrimitive> requires concepts::Cube<TPrimitive>
void flattenGeometry(Ref<const TPrimitive> r_primitive, Ptr<typename TPrimitive::Coordinate> p_output)
{
    std::copy(r_primitive.base().begin(),
              r_primitive.base().end(),
              p_output);
    std::fill(p_output + TPrimitive::Dimension,
              p_output + 2 * TPrimitive::Dimension,
              r_primitive.length());
}


template <class TPrimitive> requires concepts::Box<TPrimitive>
void flattenGeometry(Ref<const TPrimitive> r_primitive, Ptr<typename TPrimitive::Coordinate> p_output)
{
    std::copy(r_primitive.base().begin(),
              r_primitive.base().end(),
              p_output);
    std::copy(r_primitive.lengths().begin(),
              r_primitive.lengths().end(),
              p_output + TPrimitive::Dimension);
}


class BadAppleScene : public gl::Scene
{
private:
    struct Settings
    {
        Size maxDepth;
        Size samplingOrder;
        Size windowWidth;
        Size windowHeight;
        Size bufferSize;
    }; // struct Settings

public:
    BadAppleScene(const Size depth,
                  const Size samplingOrder,
                  const Size width,
                  const Size height,
                  const Size bufferSize,
                  const bool useGenericTree,
                  utils::Logger& r_logger)
        : gl::Scene("BadAppleScene",
                    gl::makeVertexShader("coloredRectangleFrame::vertexShader"),
                    gl::makeGeometryShader("coloredRectangleFrame::geometryShader"),
                    gl::makeFragmentShader("coloredRectangleFrame::fragmentShader"),
                    r_logger),
          _settings{.maxDepth       = depth,
                    .samplingOrder  = samplingOrder,
                    .windowWidth    = width,
                    .windowHeight   = height,
                    .bufferSize     = bufferSize},
          _threadPool(mp::ThreadPoolSingleton::get())
    {
        auto p_camera = this->makeCamera<gl::Camera<gl::OrthographicProjection>>(r_logger);
        this->bindUniform("transformation", p_camera->transformationMatrix());
        this->_vertexData.resize(4 * intPow(4, _settings.maxDepth));
        if (useGenericTree) {
            this->_p_root = makeRoot<Node>(samplingOrder);
        } else {
            this->_p_root = makeRoot<ContiguousTree>(samplingOrder);
        }
    }


    gl::CameraPtr getCamera()
    {return *this->_cameras.begin();}


    Size updateVertexBuffer(Ref<const gl::Image> r_image, Ref<Node> r_root)
    {
        // Define target function
        const auto target = [&r_image, this](const PointType& r_point) -> int
        {return static_cast<int>(this->sampleImage(r_point[0], r_point[1], r_image)) - 128;};

        // Fill vertex data
        Size valueCount = 0;
        const std::function<void(const Node&)> collectVertexData = [&valueCount, &target, this](const Node& r_node) mutable -> void {
            const auto depth = r_node.level();
            if (r_node.isBoundary()) {
                std::scoped_lock lock(_mutex);
                if (_vertexData.size() < valueCount + 8 + 1) {
                    _vertexData.resize(valueCount + 8 + 1);
                }
                flattenGeometry<Node>(r_node, _vertexData.data() + valueCount);

                // Highlight boundary by dimming low depth split cells
                const float relativeDepth = (depth + 1) / float(_settings.maxDepth + 1);
                std::fill(_vertexData.data() + valueCount + 4,
                          _vertexData.data() + valueCount + 7,
                          relativeDepth * relativeDepth);
                _vertexData[valueCount+7] = 1.0f;

                valueCount += 8;
            }
        };

        r_root.scan(target,
                    collectVertexData,
                    _settings.maxDepth,
                    _threadPool);

        return valueCount;
    }

    Size updateVertexBuffer(Ref<const gl::Image> r_image, Ref<ContiguousTree> r_root)
    {
        // Define target function
        const auto target = [&r_image, &r_root, this] (Ref<const ContiguousNode> r_node) -> bool {
            PointType base, lengths;
            r_root.getNodeGeometry(r_node, base.begin(), lengths.begin());

            if constexpr (concepts::Cube<ContiguousTree::Geometry>) {
                lengths[1] = lengths[0];
            }

            PointType lengthIncrements {lengths[0] / (_settings.samplingOrder - 1),
                                        lengths[1] / (_settings.samplingOrder - 1)};
            bool first = true;
            bool test = false;
            for (int i=0; i<_settings.samplingOrder; ++i) {
                for (int j=0; j<_settings.samplingOrder; ++j) {
                    if (first) {
                        first = false;
                        test = this->sampleImage(base[0] + j * lengthIncrements[0],
                                                 base[1] + i * lengthIncrements[1],
                                                 r_image) < 128;
                    }
                    else if (test != (this->sampleImage(base[0] + j * lengthIncrements[0],
                                                        base[1] + i * lengthIncrements[1],
                                                        r_image) < 128)) {
                        return true;
                    }
                } // for j in samplingOrder
            } // for i in samplingOrder
            return false;
        };

        r_root.scan(target, _settings.maxDepth);

        // Collect leaf nodes to render
        Size valueCount = 0;
        const auto collectNodes = [&r_root, &valueCount, this](Ref<const ContiguousNode> r_node, Size depth) -> bool {
            if (r_node.isLeaf()) {
                //std::scoped_lock<std::mutex> lock(_mutex);
                if (_vertexData.size() < valueCount + 8 + 1) {
                    _vertexData.resize(valueCount + 8 + 1);
                }
                r_root.getNodeGeometry(r_node,
                                       _vertexData.begin() + valueCount,
                                       _vertexData.begin() + valueCount + 2);
                if constexpr (concepts::Cube<ContiguousTree::Geometry>) {
                    _vertexData[valueCount + 3] = _vertexData[valueCount + 2];
                }

                // Highlight boundary by dimming low depth split cells
                const float relativeDepth = (depth + 1) / float(_settings.maxDepth + 1);
                std::fill(_vertexData.data() + valueCount + 4,
                          _vertexData.data() + valueCount + 7,
                          relativeDepth * relativeDepth);
                _vertexData[valueCount+7] = 1.0f;

                valueCount += 8;
            }
            return depth <= _settings.maxDepth;
        };

        r_root.visit(collectNodes);
        return valueCount;
    }

private:
    bool update_impl()
    {
        CIE_PROFILE_SCOPE

        // Check whether frames are available or there might be more frames coming in
        if (std::cin.eof() && _frameQueue.empty()) {
            return false;
        }

        for (Size i_buffer=_frameQueue.size(); i_buffer<_settings.bufferSize; ++i_buffer) {
            CIE_PROFILE_SCOPE
            // Check whether there's data to be read from stdin
            const bool newFrame = std::cin.peek() != std::cin.eof();

            // Queue loading a new frame if one is available
            if (newFrame) {
                std::string framePath;
                std::getline(std::cin, framePath);
                if (!framePath.empty()) {
                    _frameQueue.emplace(framePath);
                }
            }
        }

        if (!_frameQueue.empty()) {
            CIE_PROFILE_SCOPE
            const auto valueCount = std::visit([this](auto& rp_root){
                return this->updateVertexBuffer(_frameQueue.pop(), *rp_root);
            }, _p_root);

            this->_p_bufferManager->template writeToBoundBuffer<gl::VertexBuffer>(this->_vertexData.begin(),
                                                                                  valueCount);

            glDrawArrays(GL_POINTS, 0, valueCount);
        }

        return true;
    }

    template <concepts::Integer CT>
    static gl::Image::value_type sampleImage(CT x, CT y, const gl::Image& r_image)
    {
        CIE_OUT_OF_RANGE_CHECK(x < r_image.width())
        CIE_OUT_OF_RANGE_CHECK(y < r_image.height())
        Size index = y * r_image.width() + x;
        index *= r_image.numberOfChannels();
        return std::accumulate(r_image.data() + index,
                               r_image.data() + index + r_image.numberOfChannels(),
                               static_cast<gl::Image::value_type>(0)) / r_image.numberOfChannels();
    }

    template <class CT>
    requires (!concepts::Integer<CT>)
    static gl::Image::value_type sampleImage(CT x, CT y, const gl::Image& r_image)
    {
        return BadAppleScene::sampleImage(
            Size(x * (r_image.width()-1)),
            Size(y * (r_image.height()-1)),
            r_image
        );
    }

public:
    Settings _settings;

    mp::ThreadPool<> _threadPool;

    mp::ResourceQueue<gl::Image> _frameQueue;

    std::variant<NodePtr,std::unique_ptr<ContiguousTree>> _p_root;

    gl::VertexBuffer::data_container_type _vertexData;

    std::mutex _mutex;
};


utils::ArgParse::Results badAppleArguments(int argc, char const* argv[])
{
    utils::ArgParse parser;
    parser.addKeyword(
        {"--depth", "-d"},
        utils::ArgParse::validatorFactory<Size>(),
        utils::ArgParse::DefaultValue {"9"},
        "Max depth of the quadtree."
    );
    parser.addKeyword(
        {"--sampling-order", "-s"},
        utils::ArgParse::validatorFactory<Size>(),
        utils::ArgParse::DefaultValue {"12"},
        "Number of sample points per cell per direction."
    );
    parser.addKeyword(
        {"--width", "-w"},
        utils::ArgParse::validatorFactory<Size>(),
        utils::ArgParse::DefaultValue {"960"},
        "Output width in pixels."
    );
    parser.addKeyword(
        {"--height", "-h"},
        utils::ArgParse::validatorFactory<Size>(),
        utils::ArgParse::DefaultValue {"720"},
        "Output height in pixels."
    );
    parser.addKeyword(
        {"--buffer-size", "-b"},
        utils::ArgParse::validatorFactory<Size>(),
        utils::ArgParse::DefaultValue {"10"},
        "Max buffer size for reading frames."
    );
    parser.addFlag({"--generic-tree"}, "Use a generic tree for subdivision instead of a contiguous one.");
    parser.addFlag({"--save"}, "Output screenshots of each frame.");
    parser.addFlag({"--unlock-framerate"}, "Don't limit the framerate.");
    parser.addFlag({"--help"}, "Print this help and exit.");

    try {
        auto args = parser.parseArguments(argc - 1, argv + 1);
        if (args.get<Bool>("help")) {
            parser.help(OutputStream::cout);
        }
        return args;
    }
    catch (Exception& r_exception)
    {
        OutputStream::cerr << r_exception.what() << std::endl;
        parser.help(OutputStream::cout);
        exit(1);
    }
}


// -------------------------------------------------------------------------
// MAIN
// -------------------------------------------------------------------------

int main(int argc, char const* argv[])
{
    // Parse command line arguments
    auto args = badAppleArguments(argc, argv);
    if (args.get<Bool>("help")) {
        return 0;
    }

    const Size bufferSize    = args.get<Size>("buffer-size");
    const Size depth         = args.get<Size>("depth");
    const Size samplingOrder = args.get<Size>("sampling-order");
    const Size width         = args.get<Size>("width");
    const Size height        = args.get<Size>("height");

    // Graphics setup
    auto p_log = std::make_shared<utils::Logger>(getOutputPath() / "bad_apple.log");
    auto p_context = gl::GLFWContextSingleton::get(p_log);
    auto p_window = p_context.lock()->newWindow(width, height);
    auto p_scene = p_window->makeScene<BadAppleScene>(
        depth,
        samplingOrder,
        width,
        height,
        bufferSize,
        args.get<Bool>("generic-tree"),
        *p_log
    );

    auto p_camera = p_scene->getCamera();
    p_camera->setPose({0.5, 0.5, 0.5},
                      {0.0, 0.0, -1.0},
                      {0.0, 1.0, 0.0});
    p_camera->setClippingPlanes(0.3, 1.0);
    p_camera->setAspectRatio(double(p_window->getSize().first) / p_window->getSize().second);

    using Clock = std::chrono::high_resolution_clock;
    auto t0 = Clock::now();
    constexpr const auto period = std::chrono::microseconds(Size(1e6 / 30.0));

    const Bool writeScreenshots = args.get<Bool>("save");
    const Bool limitFrameRate = !args.get<Bool>("unlock-framerate");

    // Main loop
    Size frameIndex = 0;
    while (p_window->update()) {
        if (limitFrameRate) {
            // Delay the next frame to maintain the desired frame rate
            const auto now = Clock::now();
            const auto elapsed = now - t0;
            if (elapsed < period) {
                std::this_thread::sleep_for(period - elapsed);
                t0 = Clock::now();
            } else {
                t0 = now;
            }
        }

        // Write screenshots if requested
        if (writeScreenshots) {
            p_window->screenshot("bad_apple_output_" + std::to_string(frameIndex) + ".png");
        }

        ++frameIndex;
    }

    return 0;
}


} // namespace cie



int main(int argc, char const* argv[])
{
    return cie::main(argc, argv);
}
