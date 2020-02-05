
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
            new ShaderPipeline::ShaderProgram(GL_VERTEX_SHADER  , {vertexSource  }, true),
            new ShaderPipeline::ShaderProgram(GL_FRAGMENT_SHADER, {fragmentSource}, true)
            });

        //pipelineParticles->connect_bufferPositions_to_attribute("pos");
        bufferPositions_loc = glGetAttribLocation(_vertexShader->_glidShaderProgram, "pos");
        bufferPositions_offset = 0;
        bufferPositions_nbComponents =
            NumberOfComponentsInType(GetAttribDataType(_vertexShader->_glidShaderProgram, "pos"));
        glEnableVertexArrayAttrib(_glidVao, bufferPositions_loc);

        ////pipelineParticles->connect_bufferColors_to_attribute("color");
        //bufferColors_loc = glGetAttribLocation(_vertexShader->_glidShaderProgram, "color");
        //bufferColors_offset = 0;
        //bufferColors_nbComponents = NumberOfComponentsInType(GetAttribDataType(_vertexShader->_glidShaderProgram, "color"));
        //glEnableVertexArrayAttrib(_glidVao, bufferColors_loc);

        //pipelineParticles->connect_pointSize_to_uniform(GL_VERTEX_SHADER, "size");
        pointSize_shaderType = GL_VERTEX_SHADER;
        pointSize_loc = glGetUniformLocation(GetShader(pointSize_shaderType)->_glidShaderProgram, "size");

        //pipelineParticles->connect_vpMat_to_uniform(GL_VERTEX_SHADER, "vpMat");
        vpMat_shaderType = GL_VERTEX_SHADER;
        vpMat_loc = glGetUniformLocation(GetShader(vpMat_shaderType)->_glidShaderProgram, "vpMat");

        //pipelineParticles->in_pointSize.receive(4);

        _nbVerticesPerPrimitive = 1;
        _primitiveType = GL_POINTS;
        //_nbPrimitives.set(particles->metadataBuffer.nbElements);
    }

    void Execute() {
        glBindProgramPipeline(_glidProgramPipeline);
        glBindVertexArray(_glidVao);
        //bufferPositions.refresh();
        //bufferPositions.set(in_bufferPositions.pullUpdate());
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

        //glVertexArrayVertexBuffer(_glidVao, bufferColors_loc, bufferColors.get().bufferId, 0, bufferColors.get().nbElementsPerComponent * sizeOfEnumType(bufferColors.get().dataType));
        //glVertexArrayAttribFormat(_glidVao, bufferColors_loc, bufferColors_nbComponents, bufferColors.get().dataType, GL_FALSE, bufferColors_offset);

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



//bool Application::Init_ObstacleShader() {
//
//
//
//    return true;
//}
