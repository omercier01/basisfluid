
#ifndef SHADERPIPELINE_H
#define SHADERPIPELINE_H

#define GLM_FORCE_RADIANS
#include "GL/glew.h"
#include "glm/gtc/type_ptr.hpp"

#include <unordered_map>


class ShaderPipeline
{
public:

    enum class PrimitiveStructure {ARRAY, ELEMENT, TRANSFORM_FEEDBACK, WORK_GROUP};

    class ShaderProgram
    {
    public:
        GLenum _shaderType;
        GLuint _glidShaderProgram;
    public:
        ShaderProgram(GLenum shaderType,
                      std::initializer_list<std::string> src);
        ~ShaderProgram() {}
    };

public:
    ShaderProgram  *_vertexShader,
                   *_tessControlShader,
                   *_tessEvalShader,
                   *_geometryShader,
                   *_fragmentShader,
                   *_computeShader;
    GLuint _glidProgramPipeline;
    GLuint _glidVao;

    int _nbVerticesPerPrimitive;
    GLenum _primitiveType;
    int _nbPrimitives;

public:
    ShaderPipeline();
    void UseShaders(std::initializer_list<ShaderProgram*> shaders);
    void RemoveShader(GLenum shaderStage);
    static GLbitfield ShaderTypeEnumToBitField(GLenum shaderType);
    static unsigned int NumberOfComponentsInType(GLenum type);
    GLenum GetAttribDataType(GLuint program, std::string attribName);
    ShaderProgram* GetShader(GLenum shaderType);
};


#endif // SHADERPIPELINE_H
