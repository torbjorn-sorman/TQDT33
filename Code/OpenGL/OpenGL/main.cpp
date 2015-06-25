
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <Windows.h>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/GLU.h>


int main(int argc, char **argv)
{
    
    GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
    GLchar *str = "";
    glShaderSource(shader, 1, &str, NULL);
 

    printf("Press any key...");
    getchar();
    return 0;
}