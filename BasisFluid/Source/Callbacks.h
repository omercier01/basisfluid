

#ifndef CALLBACKS_H
#define CALLBACKS_H

#include <GL/glew.h>
#include <string>

void CallbackWindowClosed();
void CallbackKey(int key, int /*scancode*/, int action, int /*mods*/);
void CallbackMouseScroll(double /*xOffset*/, double yOffset);
void CallbackDebug(GLenum source, GLenum /*type*/, GLuint /*id*/, GLenum severity, std::string message);

#endif /*CALLBACKS_H*/