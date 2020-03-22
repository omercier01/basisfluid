
#include "ShaderPipeline.h"

#include "Application.h"

#include <glm/glm.hpp>

class ObstacleShaderPipeline : public ShaderPipeline {
public:
    ObstacleShaderPipeline() {

        const char* vertexSource = R"(
            #version 430 core

            in layout(location=0) vec2 pos;
            
            uniform mat4 vpMat;
            
            out gl_PerVertex {
                vec4 gl_Position;
            };
            
            void main() {
                gl_Position = vpMat * vec4(pos, 0, 1);
            }
        )";

        const char* fragmentSource = R"(
            #version 430 core
            
            uniform vec4 color;
            
            out vec4 fragColor;
            
            void main() {
                fragColor = color;
            }
        )";

        UseShaders({
            new ShaderPipeline::ShaderProgram(GL_VERTEX_SHADER  , {vertexSource  }),
            new ShaderPipeline::ShaderProgram(GL_FRAGMENT_SHADER, {fragmentSource})
            });

        bufferLineSegmentsVertices_loc = glGetAttribLocation(_vertexShader->_glidShaderProgram, "pos");
        bufferLineSegmentsVertices_offset = 0;
        bufferLineSegmentsVertices_nbComponents =
            NumberOfComponentsInType(GetAttribDataType(_vertexShader->_glidShaderProgram, "pos"));
        glEnableVertexArrayAttrib(_glidVao, bufferLineSegmentsVertices_loc);


        lineColor_shaderType = GL_FRAGMENT_SHADER;
        lineColor_loc = glGetUniformLocation(GetShader(lineColor_shaderType)->_glidShaderProgram, "color");

        vpMat_shaderType = GL_VERTEX_SHADER;
        vpMat_loc = glGetUniformLocation(GetShader(vpMat_shaderType)->_glidShaderProgram, "vpMat");

        _nbVerticesPerPrimitive = 2;
        _primitiveType = GL_LINES;
    }

    void Execute() {
        glBindProgramPipeline(_glidProgramPipeline);
        glBindVertexArray(_glidVao);
        glVertexArrayVertexBuffer(
            _glidVao,
            bufferLineSegmentsVertices_loc,
            app->_obstacleLines->_metadataBuffer.bufferId,
            0,
            app->_obstacleLines->_metadataBuffer.nbElementsPerComponent *
            SizeOfEnumType(app->_obstacleLines->_metadataBuffer.dataType));
        glVertexArrayAttribFormat(
            _glidVao,
            bufferLineSegmentsVertices_loc,
            bufferLineSegmentsVertices_nbComponents,
            app->_obstacleLines->_metadataBuffer.dataType,
            GL_FALSE,
            bufferLineSegmentsVertices_offset);

        glm::vec4 lineColor = { 0,0,1,1 };
        glProgramUniform4fv(GetShader(lineColor_shaderType)->_glidShaderProgram, lineColor_loc, 1, &lineColor[0]);

        glm::mat4 temp = app->_viewProjMat;
        glProgramUniformMatrix4fv(GetShader(vpMat_shaderType)->_glidShaderProgram, vpMat_loc, 1, GL_FALSE, glm::value_ptr(temp));

        glEnable(GL_DEPTH_TEST);

        _nbPrimitives = app->_obstacleLines->_metadataBuffer.nbElements / 2;

        glDrawArrays(_primitiveType, 0, _nbPrimitives * _nbVerticesPerPrimitive);
    }

    GLuint bufferLineSegmentsVertices_loc;
    GLuint bufferLineSegmentsVertices_offset;
    GLuint bufferLineSegmentsVertices_nbComponents;

    GLint lineColor_loc;
    GLenum lineColor_shaderType;

    GLint vpMat_loc;
    GLenum vpMat_shaderType;
};

