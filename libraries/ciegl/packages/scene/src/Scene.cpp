// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- Internal Includes ---
#include "packages/scene/inc/Scene.hpp"
#include "packages/buffer/inc/BufferManager.hpp"
#include "packages/buffer/inc/Buffer.hpp"
#include "packages/utility/inc/GLError.hpp"

// --- STL Includes ---
#include <algorithm>
#include <limits>


namespace cie::gl {


Size Scene::_activeSceneID = std::numeric_limits<Size>().max();


Scene::Scene(const std::string& r_name,
             Shader::SharedPointer p_vertexShader,
             Shader::SharedPointer p_geometryShader,
             Shader::SharedPointer p_fragmentShader,
             utils::Logger& r_logger,
             VertexBuffer::SharedPointer p_vertexBuffer,
             ElementBuffer::SharedPointer p_elementBuffer) :
    utils::Loggee(r_logger, r_name),
    utils::IDObject<Size>(std::numeric_limits<Size>().max()),
    _cameras(),
    _p_vertexShader(p_vertexShader),
    _p_geometryShader(p_geometryShader),
    _p_fragmentShader(p_fragmentShader),
    _p_bufferManager(nullptr),
    _uniforms(),
    _vaoID(std::numeric_limits<GLuint>().max())
{
    CIE_BEGIN_EXCEPTION_TRACING

    //glEnable(GL_DEPTH_TEST);

    // Set shaders
    this->setVertexShader(p_vertexShader);
    this->setGeometryShader(p_geometryShader);
    this->setFragmentShader(p_fragmentShader);

    // New program
    this->setID(glCreateProgram());
    this->logID("Build new program", this->getID());

    // Shaders
    if (this->_p_vertexShader)
        glAttachShader(this->getID(), _p_vertexShader->getID());
    else
        this->logID("Unset vertex shader!",
                        this->getID(),
                        LOG_TYPE_ERROR);

    if (this->_p_geometryShader)
        glAttachShader(this->getID(), _p_geometryShader->getID());

    if (this->_p_fragmentShader)
        glAttachShader(this->getID(), _p_fragmentShader->getID());
    else
        this->logID("Unset fragment shader!",
                    this->getID(),
                    LOG_TYPE_ERROR);

    // Link and check for errors
    glLinkProgram(this->getID());

    GLint linkStatus;
    glGetProgramiv(this->getID(), GL_LINK_STATUS, &linkStatus);

    if (linkStatus != GL_TRUE)
    {
        char buffer[512];
        glGetProgramInfoLog(this->getID(),
                                512,
                                NULL,
                                buffer);
        this->logID("Program linking failed!\n" + std::string(buffer),
                        this->getID(),
                        LOG_TYPE_ERROR);
    }

    // Vertex Array Object
    // TODO: wtf? also, is it okay here?
    glGenVertexArrays(1, &this->_vaoID);
    this->logID("Create Vertex Array Object", this->_vaoID);

    // Set active program and bind vertex array object
    this->activate();

    // Check buffers
    //  - no buffers should be bound before binding the vertex array object
    this->setBufferManager(BufferManager::SharedPointer(new BufferManager(this->logger())));

    // Initialize and bind buffers
    if (!p_vertexBuffer)
        p_vertexBuffer = this->_p_bufferManager->makeBuffer<VertexBuffer>();
    else
        this->_p_bufferManager->addBuffer(p_vertexBuffer);

    if (!p_elementBuffer)
        p_elementBuffer = this->_p_bufferManager->makeBuffer<ElementBuffer>();
    else
        this->_p_bufferManager->addBuffer(p_elementBuffer);

    this->_p_bufferManager->bindBuffer(p_vertexBuffer, true);
    this->_p_bufferManager->bindBuffer(p_elementBuffer, true);

    // Map shader attributes
    for (const auto& r_attribute : this->_p_vertexShader->attributes())
    {
        GLint attributeID = glGetAttribLocation(this->getID(),
                                                r_attribute.name().c_str());

        checkGLErrors(*this, "Could not find attribute: " + r_attribute.name());

        glVertexAttribPointer(attributeID,
                              r_attribute.size(),
                              GL_FLOAT,
                              GL_FALSE,
                              r_attribute.stride() * sizeof(float),
                              (void*)(r_attribute.offset() * sizeof(float)));

        glEnableVertexAttribArray(attributeID);

        checkGLErrors(*this, "Failed to enable attribute: " + r_attribute.name());
    }

    // Bind fragment shader outputs
    for (const auto& r_output : this->_p_fragmentShader->outputs())
        glBindFragDataLocation(this->getID(),
                                r_output.getID(),
                                r_output.name().c_str());

    // Create placeholders for uniforms
    for (const auto& r_uniform : this->_p_vertexShader->uniforms())
        this->_uniforms.emplace_back(
            new GLUniformPlaceholder(r_uniform, this->getID())
        );

    if (this->_p_geometryShader)
        for (const auto& r_uniform : this->_p_geometryShader->uniforms())
            this->_uniforms.emplace_back(
                new GLUniformPlaceholder(r_uniform, this->getID())
            );

    for (const auto& r_uniform : this->_p_fragmentShader->uniforms())
        this->_uniforms.emplace_back(
            new GLUniformPlaceholder(r_uniform, this->getID())
        );

    // Create empty textures
    for (const auto& r_texture : this->_p_vertexShader->textures())
        this->_textures.emplace_back(
            new Texture(r_texture)
        );

    for (const auto& r_texture : this->_p_fragmentShader->textures())
        this->_textures.emplace_back(
            new Texture(r_texture)
        );

    // Check gl errors
    if (glGetError() != 0)
        this->logID("Failed to make GL program",
                     this->getID(),
                     LOG_TYPE_ERROR);

    // Validate program
    glValidateProgram(this->getID());

    GLint programStatus;
    glGetProgramiv(this->getID(), GL_VALIDATE_STATUS, &programStatus);

    if (programStatus == GL_TRUE)
        this->logID("Program successfully validated", this->getID());
    else
    {
        GLsizei messageLength;
        GLchar buffer[512];
        glGetProgramInfoLog(this->getID(),
                             sizeof(buffer) / sizeof(GLchar),
                             &messageLength,
                             buffer);
        this->log("Program validation failed | ID_" + std::to_string(this->getID()),
                   LOG_TYPE_ERROR);
    }

    CIE_END_EXCEPTION_TRACING
}


Scene::~Scene()
{
    this->logID("Delete program", this->getID());
    glDeleteProgram(this->getID());

    this->logID("Delete Vertex Array Object", this->_vaoID);
    glDeleteVertexArrays(1, &_vaoID);
}


bool Scene::update()
{
    CIE_BEGIN_EXCEPTION_TRACING

    this->activate(false);

    // Update all cameras
    for (auto& rp_camera : this->_cameras)
        rp_camera->update();

    // Check whether the number of vertex buffers
    // and element buffers is identical
    Size numberOfBuffers = this->_p_bufferManager->elementBuffers().size();

    if (this->_p_bufferManager->vertexBuffers().size() != numberOfBuffers) [[unlikely]]
        CIE_THROW(OutOfRangeException, "Buffer size mismatch")

    // Save current bound buffers for resetting after the update
    auto p_boundVertexBuffer = this->_p_bufferManager->template getBoundBuffer<VertexBuffer>().lock();
    auto p_boundElementBuffer = this->_p_bufferManager->template getBoundBuffer<ElementBuffer>().lock();

    // Update uniforms
    for (const auto& rp_uniform : this->_uniforms)
        rp_uniform->write();

    checkGLErrors(*this, "Error updating uniforms");

    // Loop through buffers and draw their contents
    auto it_vertexBuffer    = this->_p_bufferManager->vertexBuffers().begin();
    auto it_elementBuffer   = this->_p_bufferManager->elementBuffers().begin();
    auto it_vertexBufferEnd = this->_p_bufferManager->vertexBuffers().end();

    bool output = true;
    for (; it_vertexBuffer!=it_vertexBufferEnd; ++it_vertexBuffer, ++it_elementBuffer)
    {
        this->_p_bufferManager->bindBuffer(
            *it_vertexBuffer,
            false
        );
        this->_p_bufferManager->bindBuffer(
            *it_elementBuffer,
            false
        );

        // Draw
        output = output && this->update_impl();

        checkGLErrors(*this,
                       "Error drawing from buffer");
    }

    // Reset bound buffers
    this->_p_bufferManager->bindBuffer(p_boundVertexBuffer, false);
    this->_p_bufferManager->bindBuffer(p_boundElementBuffer, false);

    return output;

    CIE_END_EXCEPTION_TRACING
}


void Scene::activate(bool log)
{
    if (!this->isActive())
    {
        glUseProgram(this->getID());
        if (log)
            this->logID("Set active program", this->getID());

        glBindVertexArray(this->_vaoID);
        if (log)
            this->logID("Bind Vertex Array Object", this->_vaoID);

        Scene::_activeSceneID = this->getID();

        checkGLErrors(*this, "Scene activation failed");
    }
}


bool Scene::isActive() const
{
    return this->getID() == Scene::_activeSceneID;
}


void Scene::addCamera(CameraPtr p_camera)
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto it_camera = std::find(_cameras.begin(),
                               _cameras.end(),
                               p_camera);

    CIE_CHECK(it_camera == _cameras.end(),
              "Camera '" << p_camera->name() << "' is already part of scene '" << this->name() << "'")

    _cameras.push_back(p_camera);

    CIE_END_EXCEPTION_TRACING
}


void Scene::removeCamera(CameraPtr p_camera)
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto it_camera = std::find(
        _cameras.begin(),
        _cameras.end(),
        p_camera
    );

    if (it_camera == _cameras.end())
        CIE_THROW(OutOfRangeException, "Attempt to remove non-existent camera")

    _cameras.erase(it_camera);

    CIE_END_EXCEPTION_TRACING
}


const typename Scene::camera_container& Scene::cameras() const
{
    return _cameras;
}


CameraPtr Scene::camera(Size index)
{
    CIE_OUT_OF_RANGE_CHECK(index < this->_cameras.size())

    auto itp_camera = this->_cameras.begin();
    std::advance(itp_camera, index);

    return *itp_camera;
}


Scene& Scene::setVertexShader(Shader::SharedPointer p_vertexShader)
{
    _p_vertexShader = p_vertexShader;
    return *this;
}


Scene& Scene::setGeometryShader(Shader::SharedPointer p_geometryShader)
{
    _p_geometryShader = p_geometryShader;
    return *this;
}


Scene& Scene::setFragmentShader(Shader::SharedPointer p_fragmentShader)
{
    _p_fragmentShader = p_fragmentShader;
    return *this;
}


const Shader::SharedPointer& Scene::vertexShader() const
{
    return _p_vertexShader;
}


const Shader::SharedPointer& Scene::geometryShader() const
{
    return _p_geometryShader;
}


const Shader::SharedPointer& Scene::fragmentShader() const
{
    return _p_fragmentShader;
}


Scene& Scene::setBufferManager(BufferManager::SharedPointer p_bufferManager)
{
    _p_bufferManager = p_bufferManager;
    return *this;
}


const BufferManager::SharedPointer& Scene::bufferManager() const
{
    return _p_bufferManager;
}


const Scene::uniform_container& Scene::uniforms() const
{
    return this->_uniforms;
}


const Scene::texture_container& Scene::textures() const
{
    return this->_textures;
}


void Scene::bindUniform(const std::string& r_name,
                         const glm::mat4& r_uniform)
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Find uniform in the internal list of uniforms
    GLUniformPtr& rp_uniform = this->findUniform(r_name);

    // Check uniform type
    CIE_CHECK(
        rp_uniform->properties().type == GL_FLOAT_MAT4,
        "Attempt to bind mismatching uniform types: " + r_name
        + " | type: " + std::to_string(rp_uniform->properties().type)
    )

    // Replace the uniform
    rp_uniform.reset(
        new FloatMat4GLUniform(*rp_uniform,
                                this->getID(),
                                r_uniform)
    );

    CIE_END_EXCEPTION_TRACING
}


void Scene::bindUniform(const std::string& r_name,
                         const glm::dvec3& r_uniform)
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Find uniform in the internal list of uniforms
    GLUniformPtr& rp_uniform = this->findUniform(r_name);

    // Check uniform type
    CIE_CHECK(
        rp_uniform->properties().type == GL_FLOAT_VEC3,
        "Attempt to bind mismatching uniform types: " + r_name
        + " | type: " + std::to_string(rp_uniform->properties().type)
    )

    // Replace the uniform
    rp_uniform.reset(
        new FloatVec3GLUniform<glm::dvec3>(*rp_uniform,
                                            this->getID(),
                                            r_uniform)
    );

    CIE_END_EXCEPTION_TRACING
}


void Scene::bindUniform(const std::string& r_name,
                         const glm::vec3& r_uniform)
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Find uniform in the internal list of uniforms
    GLUniformPtr& rp_uniform = this->findUniform(r_name);

    // Check uniform type
    CIE_CHECK(
        rp_uniform->properties().type == GL_FLOAT_VEC3,
        "Attempt to bind mismatching uniform types: " + r_name
        + " | type: " + std::to_string(rp_uniform->properties().type)
    )

    // Replace the uniform
    rp_uniform.reset(
        new FloatVec3GLUniform<glm::vec3>(*rp_uniform,
                                           this->getID(),
                                           r_uniform)
    );

    CIE_END_EXCEPTION_TRACING
}


void Scene::bindUniform(const std::string& r_name,
                         const GLfloat& r_uniform)
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Find uniform in the internal list of uniforms
    GLUniformPtr& rp_uniform = this->findUniform(r_name);

    // Check uniform type
    CIE_CHECK(
        rp_uniform->properties().type == GL_FLOAT,
        "Attempt to bind mismatching uniform types: " + r_name
        + " | type: " + std::to_string(rp_uniform->properties().type)
    )

    // Replace the uniform
    rp_uniform.reset(
        new FloatGLUniform(*rp_uniform,
                            this->getID(),
                            r_uniform)
    );

    CIE_END_EXCEPTION_TRACING
}


GLUniformPtr& Scene::findUniform(const std::string& r_name)
{
    auto itp_uniform = std::find_if(
        this->_uniforms.begin(),
        this->_uniforms.end(),
        [&r_name](const GLUniformPtr& rp_uniform) -> bool
        {
            if (rp_uniform->name() == r_name)     [[unlikely]]
                return true;
            else                                    [[likely]]
                return false;
        }
    );

    CIE_CHECK(
        itp_uniform != this->_uniforms.end(),
        "Uniform not found: " + r_name
    )

    return *itp_uniform;
}


TexturePtr& Scene::findTexture(const std::string& r_name)
{
    auto it = std::find_if(
        this->_textures.begin(),
        this->_textures.end(),
        [&r_name](const TexturePtr& rp_uniform) -> bool
        {
            if (rp_uniform->name() == r_name)     [[unlikely]]
                return true;
            else                                    [[likely]]
                return false;
        }
    );

    CIE_CHECK(
        it != this->_textures.end(),
        "Uniform not found: " + r_name
    )

    return *it;
}


} // namespace cie::gl