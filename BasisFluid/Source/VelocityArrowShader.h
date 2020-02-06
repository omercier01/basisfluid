
#include "ShaderPipeline.h"

#include "Application.h"

#include <glm/glm.hpp>

class VelocityArrowShaderPipeline : public ShaderPipeline {
public:
    VelocityArrowShaderPipeline() {

        const char* vertexSource = R"(
            #version 430 core
            
            in layout(location=0) vec2 pos;
            in layout(location=1) vec2 vec;
            
            uniform float pointSize;
            
            out gl_PerVertex {
                vec4 gl_Position;
                float gl_PointSize;
            };
            out layout(location = 1) vec2 vecOut;
            
            void main() {
                gl_Position = vec4(pos, 0, 1);
                gl_PointSize = pointSize;
                vecOut = vec;
            }
        )";

        const char* geometrySource = R"(
            #version 430 core
            
            layout(points) in;
            layout(triangle_strip, max_vertices=18) out;
    
            in gl_PerVertex {
                vec4 gl_Position;
            } gl_in[];
            in layout(location=1) vec2 vec[];

            uniform float lengthFactor;
            uniform mat4 vpMat;

            out gl_PerVertex
            {
                vec4 gl_Position;
            };
    
            void main() {
                vec2 begin = gl_in[0].gl_Position.xy;
                vec2 end = begin + lengthFactor*vec[0];
                vec2 right = lengthFactor*vec2(vec[0].y,-vec[0].x);
                vec2 branchRight = -1.0/sqrt(2)/4.0*lengthFactor*vec2(vec[0].x-vec[0].y, vec[0].x+vec[0].y);
                vec2 branchLeft = vec2(branchRight.y, -branchRight.x);
                        
                float th = 0.1;//length(end-begin);
                float factStem = 0.5;
                        
                gl_Position = vpMat * vec4(begin-factStem*th*right, 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(begin+factStem*th*right, 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(end-factStem*th*right, 0, 1);
                EmitVertex();
                EndPrimitive();
                        
                gl_Position = vpMat * vec4(end-factStem*th*right, 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(begin+factStem*th*right, 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(end+factStem*th*right, 0, 1);
                EmitVertex();
                EndPrimitive();
                        
                gl_Position = vpMat * vec4(end+branchRight+th*branchLeft, 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(end+branchRight-2*th*branchLeft, 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(end, 0, 1);
                EmitVertex();
                EndPrimitive();
                        
                gl_Position = vpMat * vec4(end+th*(begin-end), 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(end+branchRight-2*th*branchLeft, 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(end-2*th*(begin-end), 0, 1);
                EmitVertex();
                EndPrimitive();
                        
                gl_Position = vpMat * vec4(end+branchLeft-2*th*branchRight, 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(end+branchLeft+th*branchRight, 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(end-2*th*(begin-end), 0, 1);
                EmitVertex();
                EndPrimitive();
                        
                gl_Position = vpMat * vec4(end-2*th*(begin-end), 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(end+branchLeft+th*branchRight, 0, 1);
                EmitVertex();
                gl_Position = vpMat * vec4(end, 0, 1);
                EmitVertex();
                EndPrimitive();
            }
        )";


        const char* fragmentSource = R"(
            #version 430 core

            out vec4 color;
        
            uniform vec3 arrowColor;
            
            void main() {
                color = vec4(arrowColor, 1);
            }
        )";


        UseShaders({
            new ShaderPipeline::ShaderProgram(GL_VERTEX_SHADER  , {vertexSource  }, true),
            new ShaderPipeline::ShaderProgram(GL_GEOMETRY_SHADER, {geometrySource}, true),
            new ShaderPipeline::ShaderProgram(GL_FRAGMENT_SHADER, {fragmentSource}, true)
            });

        //pipelineArrows->connect_bufferStart_to_attribute("pos");
        bufferStart_loc = glGetAttribLocation(_vertexShader->_glidShaderProgram, "pos");
        bufferStart_offset = 0;
        bufferStart_nbComponents =
            NumberOfComponentsInType(GetAttribDataType(_vertexShader->_glidShaderProgram, "pos"));
        glEnableVertexArrayAttrib(_glidVao, bufferStart_loc);

        //pipelineArrows->connect_bufferVectors_to_attribute("vec");
        bufferVectors_loc = glGetAttribLocation(_vertexShader->_glidShaderProgram, "vec");
        bufferVectors_offset = 0;
        bufferVectors_nbComponents =
            NumberOfComponentsInType(GetAttribDataType(_vertexShader->_glidShaderProgram, "vec"));
        glEnableVertexArrayAttrib(_glidVao, bufferVectors_loc);

        //pipelineArrows->connect_arrowColor_to_uniform(GL_FRAGMENT_SHADER, "arrowColor");
        arrowColor_shaderType = GL_FRAGMENT_SHADER;
        arrowColor_loc = glGetUniformLocation(GetShader(arrowColor_shaderType)->_glidShaderProgram, "arrowColor");

        //pipelineArrows->connect_arrowLengthFactor_to_uniform(GL_GEOMETRY_SHADER, "lengthFactor");
        arrowLengthFactor_shaderType = GL_GEOMETRY_SHADER;
        arrowLengthFactor_loc = glGetUniformLocation(GetShader(arrowLengthFactor_shaderType)->_glidShaderProgram, "lengthFactor");

        //pipelineArrows->connect_vpMat_to_uniform(GL_GEOMETRY_SHADER, "vpMat");
        vpMat_shaderType = GL_GEOMETRY_SHADER;
        vpMat_loc = glGetUniformLocation(GetShader(vpMat_shaderType)->_glidShaderProgram, "vpMat");

        _nbVerticesPerPrimitive = 1;
        _primitiveType = GL_POINTS;
        //_nbPrimitives.set(bufferGridPoints->metadataBuffer.nbElements);

        /*createPullConnection(&velocityField->out_gridNodeLocationsCpu,
            &bufferGridPoints->in_metadataCpu)->activate();

        createPullConnection(&velocityField->out_vectorValuesCpu,
            &bufferArrows->in_metadataCpu)->activate();


        createPullConnection(&bufferGridPoints->out_metadataBuffer,
            &pipelineArrows->in_bufferStart)->activate();


        createPullConnection(&bufferArrows->out_metadataBuffer,
            &pipelineArrows->in_bufferVectors)->activate();

        createPullConnection(&cam->out_viewProjMatrix,
            &pipelineArrows->in_vpMat)->activate();*/

    }

    void Execute() {


        glBindProgramPipeline(_glidProgramPipeline);
        glBindVertexArray(_glidVao);

        glVertexArrayVertexBuffer(
            _glidVao,
            bufferStart_loc,
            app->_bufferGridPoints->_metadataBuffer.bufferId,
            0,
            app->_bufferGridPoints->_metadataBuffer.nbElementsPerComponent *
            SizeOfEnumType(app->_bufferGridPoints->_metadataBuffer.dataType)
        );
        glVertexArrayAttribFormat(
            _glidVao,
            bufferStart_loc,
            bufferStart_nbComponents,
            app->_bufferGridPoints->_metadataBuffer.dataType,
            GL_FALSE,
            bufferStart_offset
        );

        glVertexArrayVertexBuffer(
            _glidVao,
            bufferVectors_loc,
            app->_bufferArrows->_metadataBuffer.bufferId,
            0,
            app->_bufferArrows->_metadataBuffer.nbElementsPerComponent *
            SizeOfEnumType(app->_bufferArrows->_metadataBuffer.dataType)
        );
        glVertexArrayAttribFormat(
            _glidVao,
            bufferVectors_loc,
            bufferVectors_nbComponents,
            app->_bufferArrows->_metadataBuffer.dataType,
            GL_FALSE,
            bufferVectors_offset
        );

        //pipelineArrows->in_arrowColor.receive(vec3(1, 0, 0));
        {
            vec3 temp = { 1,0,0 };
            glProgramUniform3f(GetShader(arrowColor_shaderType)->_glidShaderProgram, arrowColor_loc, temp.x, temp.y, temp.z);
        }

        //pipelineArrows->in_arrowLengthFactor.receive(0.1f);
        glProgramUniform1f(
            GetShader(arrowLengthFactor_shaderType)->_glidShaderProgram,
            arrowLengthFactor_loc,
            app->_velocityArrowFactor
        );

        {
            glm::mat4 temp = app->_viewProjMat;
            glProgramUniformMatrix4fv(GetShader(vpMat_shaderType)->_glidShaderProgram, vpMat_loc, 1, GL_FALSE, glm::value_ptr(temp));
        }

        glDisable(GL_DEPTH_TEST);

        //_nbPrimitives.set(bufferGridPoints->metadataBuffer.nbElements);
        _nbPrimitives = app->_bufferArrows->_metadataBuffer.nbElements;

        glDrawArrays(_primitiveType, 0, _nbPrimitives * _nbVerticesPerPrimitive);

    }

    GLuint bufferStart_loc;
    GLuint bufferStart_offset;
    GLuint bufferStart_nbComponents;

    GLuint bufferVectors_loc;
    GLuint bufferVectors_offset;
    GLuint bufferVectors_nbComponents;

    GLint arrowColor_loc;
    GLenum arrowColor_shaderType;

    GLint arrowLengthFactor_loc;
    GLenum arrowLengthFactor_shaderType;

    GLint vpMat_loc;
    GLenum vpMat_shaderType;

};


