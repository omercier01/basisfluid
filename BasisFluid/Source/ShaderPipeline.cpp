
#include "ShaderPipeline.h"

#include <fstream>
#include <string>
#include <iostream>

using namespace std;

ShaderPipeline::ShaderProgram::ShaderProgram(
    GLenum shaderType,
    std::initializer_list<string> srcFilenames,
    bool readSrcFilenamesAsSourceCode)
{
    _shaderType = shaderType;

    size_t nbSources = srcFilenames.size();
    const char ** sourceCodes = new const char *[nbSources];
    int count = 0;
    string* sourceCodeStrings;

    if (readSrcFilenamesAsSourceCode) {
        for (initializer_list<string>::iterator iSrc = srcFilenames.begin();
            iSrc != srcFilenames.end(); iSrc++)
        {
            sourceCodes[count++] = iSrc->c_str();
        }

    }
    else {
        // Read the Vertex Shader code from the files
        // based on http://www.opengl-tutorial.org/beginners-tutorials/tutorial-2-the-first-triangle/

        sourceCodeStrings = new string[nbSources];

        for (initializer_list<string>::iterator iSrc = srcFilenames.begin();
            iSrc != srcFilenames.end(); iSrc++)
        {
            ifstream sourceCodeStream(*iSrc, ios::in);
            if (sourceCodeStream.is_open())
            {
                string line = "";
                while (getline(sourceCodeStream, line)) {
                    sourceCodeStrings[count] += line + "\n";
                }
                sourceCodeStream.close();
            }
            else {
                cout << "Error loading shader : " << *iSrc << endl;
            }
            sourceCodes[count] = sourceCodeStrings[count].c_str();
            count++;
        }

    }

    GLint result;
    int infoLength;
    char* info;

    const GLuint shader = glCreateShader(shaderType);
    glShaderSource(shader, GLsizei(nbSources), sourceCodes, NULL);
    glCompileShader(shader);
    glGetShaderiv(shader, GL_COMPILE_STATUS, &result);
    if (!result) {
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLength);
        info = new char[infoLength];
        glGetShaderInfoLog(shader, infoLength, NULL, info);
        cerr << "Compiler error in shader :" << endl;
        cerr << info << endl;
    }
    else {
        _glidShaderProgram = glCreateProgram();
        glProgramParameteri(_glidShaderProgram, GL_PROGRAM_SEPARABLE, GL_TRUE);
        glAttachShader(_glidShaderProgram, shader);
        glLinkProgram(_glidShaderProgram);
        glDetachShader(_glidShaderProgram, shader);
        glDeleteShader(shader);

        glGetProgramiv(_glidShaderProgram, GL_LINK_STATUS, &result);
        if (!result) {
            glGetProgramiv(_glidShaderProgram, GL_INFO_LOG_LENGTH, &infoLength);
            info = new char[infoLength];
            glGetProgramInfoLog(_glidShaderProgram, infoLength, NULL, info);
            cerr << "Linker error in shader :" << endl;
            cerr << info << endl;
        }
    }


}


GLenum ShaderPipeline::GetAttribDataType(GLuint program, string attribName)
{
    // need to loop over all attributes becuase glGetActiveAttrib apparently does not use the order of the vertex attrib (i.e. the order of glGetAttribLocation).
    GLint nbActiveAttribs;
    glGetProgramiv(program, GL_ACTIVE_ATTRIBUTES, &nbActiveAttribs);
    GLsizei tempLength; GLint tempSize; GLenum tempType; GLchar* tempName;
    tempName = new GLchar[attribName.size() + 1]; // +1 so if the returned name has the same beginning but is longer, we'll discard it.
    for (int i = 0; i < nbActiveAttribs; i++) {
        glGetActiveAttrib(program, i, GLsizei(attribName.size() + 1), &tempLength, &tempSize, &tempType, tempName);
        if (!strcmp(attribName.c_str(), tempName)) { // if same name
            return tempType;
        }
    }

    // if not found
    cout << "ShaderPipeline::GetAttribDataType : Unknown attribute name." << endl;
    return 0;
}



// adapted from https://chromium.googlesource.com/angle/angle/+/chromium/2027/src/common/utilities.cpp
unsigned int ShaderPipeline::NumberOfComponentsInType(GLenum type)
{
    switch (type)
    {
    case GL_BOOL:
    case GL_FLOAT:
    case GL_INT:
    case GL_SAMPLER_2D:
    case GL_SAMPLER_3D:
    case GL_SAMPLER_CUBE:
    case GL_SAMPLER_2D_ARRAY:
    case GL_INT_SAMPLER_2D:
    case GL_INT_SAMPLER_3D:
    case GL_INT_SAMPLER_CUBE:
    case GL_INT_SAMPLER_2D_ARRAY:
    case GL_UNSIGNED_INT_SAMPLER_2D:
    case GL_UNSIGNED_INT_SAMPLER_3D:
    case GL_UNSIGNED_INT_SAMPLER_CUBE:
    case GL_UNSIGNED_INT_SAMPLER_2D_ARRAY:
    case GL_SAMPLER_2D_SHADOW:
    case GL_SAMPLER_CUBE_SHADOW:
    case GL_SAMPLER_2D_ARRAY_SHADOW:
    case GL_UNSIGNED_INT:
        return 1;
    case GL_BOOL_VEC2:
    case GL_FLOAT_VEC2:
    case GL_INT_VEC2:
    case GL_UNSIGNED_INT_VEC2:
        return 2;
    case GL_INT_VEC3:
    case GL_FLOAT_VEC3:
    case GL_BOOL_VEC3:
    case GL_UNSIGNED_INT_VEC3:
        return 3;
    case GL_BOOL_VEC4:
    case GL_FLOAT_VEC4:
    case GL_INT_VEC4:
    case GL_UNSIGNED_INT_VEC4:
    case GL_FLOAT_MAT2:
        return 4;
    case GL_FLOAT_MAT2x3:
    case GL_FLOAT_MAT3x2:
        return 6;
    case GL_FLOAT_MAT2x4:
    case GL_FLOAT_MAT4x2:
        return 8;
    case GL_FLOAT_MAT3:
        return 9;
    case GL_FLOAT_MAT3x4:
    case GL_FLOAT_MAT4x3:
        return 12;
    case GL_FLOAT_MAT4:
        return 16;
    default:
        cout << "ShaderPipeline::numberOfComponentsInType(...) : Unknown data type." << endl;
        return NULL;
        break;
    }
    return 0;
}


GLbitfield ShaderPipeline::ShaderTypeEnumToBitField(
    GLenum shaderType)
{
    GLbitfield shaderBit;
    switch (shaderType) {
    case GL_VERTEX_SHADER:
        shaderBit = GL_VERTEX_SHADER_BIT;
        break;
    case GL_TESS_CONTROL_SHADER:
        shaderBit = GL_TESS_CONTROL_SHADER_BIT;
        break;
    case GL_TESS_EVALUATION_SHADER:
        shaderBit = GL_TESS_EVALUATION_SHADER_BIT;
        break;
    case GL_GEOMETRY_SHADER:
        shaderBit = GL_GEOMETRY_SHADER_BIT;
        break;
    case GL_FRAGMENT_SHADER:
        shaderBit = GL_FRAGMENT_SHADER_BIT;
        break;
    case GL_COMPUTE_SHADER:
        shaderBit = GL_COMPUTE_SHADER_BIT;
        break;
    }
    return shaderBit;
}

ShaderPipeline::ShaderProgram *ShaderPipeline::GetShader(GLenum shaderType)
{
    switch (shaderType) {
    case GL_VERTEX_SHADER:
        return _vertexShader;
        break;
    case GL_TESS_CONTROL_SHADER:
        return _tessControlShader;
        break;
    case GL_TESS_EVALUATION_SHADER:
        return _tessEvalShader;
        break;
    case GL_GEOMETRY_SHADER:
        return _geometryShader;
        break;
    case GL_FRAGMENT_SHADER:
        return _fragmentShader;
        break;
    case GL_COMPUTE_SHADER:
        return _computeShader;
        break;
    default:
        cout << "ProgramPipeline::ShaderProgram::getShader(...) : Unknown shader type." << endl;
        return NULL;
        break;
    }
}


ShaderPipeline::ShaderPipeline()
{
    _vertexShader = NULL;
    _tessControlShader = NULL;
    _tessEvalShader = NULL;
    _geometryShader = NULL;
    _fragmentShader = NULL;
    _computeShader = NULL;
    glGenProgramPipelines(1, &_glidProgramPipeline);
    glBindProgramPipeline(_glidProgramPipeline);

    glGenVertexArrays(1, &_glidVao);
    glBindVertexArray(_glidVao);
}

void ShaderPipeline::UseShaders(
    std::initializer_list<ShaderPipeline::ShaderProgram *> shaders)
{
    for (initializer_list<ShaderProgram*>::iterator iShader = shaders.begin();
        iShader != shaders.end(); iShader++)
    {
        GLbitfield shaderBit = ShaderTypeEnumToBitField((*iShader)->_shaderType);

        // make sure we don't mix compute shaders with other shader types.
        if ((*iShader)->_shaderType == GL_COMPUTE_SHADER) {
            glUseProgramStages(_glidProgramPipeline, GL_ALL_SHADER_BITS, 0);
        }
        else {
            glUseProgramStages(_glidProgramPipeline, GL_COMPUTE_SHADER_BIT, 0);
        }

        // store shader ID
        switch ((*iShader)->_shaderType) {
        case GL_VERTEX_SHADER:
            _vertexShader = *iShader;
            break;
        case GL_TESS_CONTROL_SHADER:
            _tessControlShader = *iShader;
            break;
        case GL_TESS_EVALUATION_SHADER:
            _tessEvalShader = *iShader;
            break;
        case GL_GEOMETRY_SHADER:
            _geometryShader = *iShader;
            break;
        case GL_FRAGMENT_SHADER:
            _fragmentShader = *iShader;
            break;
        case GL_COMPUTE_SHADER:
            _computeShader = *iShader;
            break;
        }

        // attach shader shage
        glUseProgramStages(_glidProgramPipeline, shaderBit,
            (*iShader)->_glidShaderProgram);
    }

    // check for validation errors.
    GLint result;
    int infoLength;
    char* info;
    glValidateProgramPipeline(_glidProgramPipeline);
    glGetProgramPipelineiv(_glidProgramPipeline, GL_VALIDATE_STATUS, &result);
    if (!result) {
        glGetProgramPipelineiv(_glidProgramPipeline, GL_INFO_LOG_LENGTH, &infoLength);
        info = new char[infoLength];
        glGetProgramPipelineInfoLog(_glidProgramPipeline, infoLength, NULL, info);
        cerr << "Validation error in program pipeline :" << endl;
        cerr << info << endl;
    }

}

void ShaderPipeline::RemoveShader(GLenum shaderStage)
{
    glUseProgramStages(_glidProgramPipeline,
        ShaderTypeEnumToBitField(shaderStage), 0);

    // store shader ID
    switch (shaderStage) {
    case GL_VERTEX_SHADER:
        _vertexShader = 0;
        break;
    case GL_TESS_CONTROL_SHADER:
        _tessControlShader = 0;
        break;
    case GL_TESS_EVALUATION_SHADER:
        _tessEvalShader = 0;
        break;
    case GL_GEOMETRY_SHADER:
        _geometryShader = 0;
        break;
    case GL_FRAGMENT_SHADER:
        _fragmentShader = 0;
        break;
    case GL_COMPUTE_SHADER:
        _computeShader = 0;
        break;
    }
}












