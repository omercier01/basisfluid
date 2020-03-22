
#include "ShaderPipeline.h"

#include "Application.h"

#include <glm/glm.hpp>

class ParticleShaderPipeline : public ShaderPipeline {
public:
    ParticleShaderPipeline() {

        const char* vertexSource = R"(
            #version 430 core
            in layout(location=0) vec2 pos;
            //in layout(location=1) vec3 color;
            
            uniform float size;
            uniform mat4 vpMat;
            
            out gl_PerVertex {
                vec4 gl_Position;
                float gl_PointSize;
            };
            out vec3 vColor;
            void main() {
                gl_Position = vpMat * vec4(pos, 0, 1);
                gl_PointSize = size;
                //vColor = color;
                vColor = vec3(0,0,0);
            }
        )";

        const char* fragmentSource = R"(
            #version 430 core
            
            in vec3 vColor;
                in flat uint bShow;
            
            out vec4 outColor;
            void main() {
                if(distance(gl_PointCoord, vec2(0.5,0.5)) < 0.5) {
                    outColor = vec4(vColor,0.1);
                } else {
                    discard;
                }
            }
        )";

        UseShaders({
            new ShaderPipeline::ShaderProgram(GL_VERTEX_SHADER  , {vertexSource  }),
            new ShaderPipeline::ShaderProgram(GL_FRAGMENT_SHADER, {fragmentSource})
            });

        bufferPositions_loc = glGetAttribLocation(_vertexShader->_glidShaderProgram, "pos");
        bufferPositions_offset = 0;
        bufferPositions_nbComponents =
            NumberOfComponentsInType(GetAttribDataType(_vertexShader->_glidShaderProgram, "pos"));
        glEnableVertexArrayAttrib(_glidVao, bufferPositions_loc);

        pointSize_shaderType = GL_VERTEX_SHADER;
        pointSize_loc = glGetUniformLocation(GetShader(pointSize_shaderType)->_glidShaderProgram, "size");

        vpMat_shaderType = GL_VERTEX_SHADER;
        vpMat_loc = glGetUniformLocation(GetShader(vpMat_shaderType)->_glidShaderProgram, "vpMat");

        _nbVerticesPerPrimitive = 1;
        _primitiveType = GL_POINTS;
    }

    void Execute() {
        glBindProgramPipeline(_glidProgramPipeline);
        glBindVertexArray(_glidVao);
        glVertexArrayVertexBuffer(
            _glidVao,
            bufferPositions_loc,
            app->_partPos->_metadataBuffer.bufferId,
            0,
            app->_partPos->_metadataBuffer.nbElementsPerComponent *
                SizeOfEnumType(app->_partPos->_metadataBuffer.dataType));
        glVertexArrayAttribFormat(
            _glidVao,
            bufferPositions_loc,
            bufferPositions_nbComponents,
            app->_partPos->_metadataBuffer.dataType,
            GL_FALSE,
            bufferPositions_offset);

        glProgramUniform1f(GetShader(pointSize_shaderType)->_glidShaderProgram, pointSize_loc, 4.f);

        glm::mat4 temp = app->_viewProjMat;
        glProgramUniformMatrix4fv(GetShader(vpMat_shaderType)->_glidShaderProgram, vpMat_loc, 1, GL_FALSE, glm::value_ptr(temp));


        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
        glDisable(GL_DEPTH_TEST);

        _nbPrimitives = app->_partPos->_metadataBuffer.nbElements;

        glDrawArrays(_primitiveType, 0, _nbPrimitives * _nbVerticesPerPrimitive);
    }

    GLuint bufferPositions_loc;
    GLuint bufferPositions_offset;
    GLuint bufferPositions_nbComponents;

    GLuint bufferColors_loc;
    GLuint bufferColors_offset;
    GLuint bufferColors_nbComponents;

    GLint pointSize_loc;
    GLenum pointSize_shaderType;

    GLint vpMat_loc;
    GLenum vpMat_shaderType;
};

