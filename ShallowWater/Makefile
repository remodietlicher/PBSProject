THIRDPARTY = thirdparty_osx
GLEW_INCLUDE = /usr/local/Cellar/glew/1.10.0/include
GLEW_LIB = /usr/local/Cellar/glew/1.10.0/lib/
GLFW_INCLUDE = $(THIRDPARTY)/glfw/include
GLFW_LIB = $(THIRDPARTY)/glfw/
GLM_INCLUDE = $(THIRDPARTY)/glm
STBIMG_INCLUDE = $(THIRDPARTY)/stb_image

LIBDIR = -L$(GLEW_LIB) -L$(GLFW_LIB)
INCLUDEDIR = -I$(GLEW_INCLUDE) -I$(GLFW_INCLUDE) -I$(GLM_INCLUDE) -I$(STBIMG_INCLUDE)

LIBS = -lGLEW -lglfw -lGLEW

myOpenGL: source/*.cpp source/tdogl/*.cpp

	clang++ $(LIBDIR) -o $@ $^ $(INCLUDEDIR) $(LIBS) -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo

# .c.o:
# 	clang++ -c -o $@ $< -I$(GLEW_INCLUDE) -I$(GLFW_INCLUDE)

clean:
	rm -f *.o; rm myOpenGL
