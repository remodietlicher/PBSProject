/*
 main
 
 Copyright 2012 Thomas Dalling - http://tomdalling.com/

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

// third-party libraries
#include <GL/glew.h>
#include <GL/glfw.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// standard C++ libraries
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <cmath>

// tdogl classes
#include "tdogl/Program.h"
#include "tdogl/Texture.h"
#include "tdogl/Camera.h"

// Shallow Water
#include "sw_grid.h"
#include "sws.h"

// constants
const glm::vec2 SCREEN_SIZE(800, 600);

// globals for GPU
tdogl::Program* gProgram = NULL;
tdogl::Camera gCamera;
GLuint gVAO      = 0;
GLuint gVertexXY = 0; // Vertex positions
GLuint gVertexZ  = 0; // The water
GLuint gVertexG  = 0; // The ground
GLuint gNormal   = 0; // the normal of the water
GLuint gIndex    = 0; // the triangle mesh indices
GLuint gBoxVertXY= 0;
GLuint gBoxVertZ = 0;
GLuint gBoxIndex = 0;

// globals for CPU
GLfloat *vertexXY; // vertex positons in CPU RAM

// loads the vertex shader and fragment shader, and links them to make the global gProgram
static void LoadShaders() {
    std::vector<tdogl::Shader> shaders;
    shaders.push_back(tdogl::Shader::shaderFromFile("C:/Users/remo/PBSProject/ShallowWater/resources/vertex-shader.txt", GL_VERTEX_SHADER));
    shaders.push_back(tdogl::Shader::shaderFromFile("C:/Users/remo/PBSProject/ShallowWater/resources/fragment-shader.txt", GL_FRAGMENT_SHADER));
    gProgram = new tdogl::Program(shaders);
}

void computeNormals(GLfloat *vertexXY, GLfloat *vertexZ, GLfloat *normalData, Sw_grid *grid){

    /*
        computes the normals for vertexXY & vertexZ
        saves into normalData
    */


    int nCellsX = grid->xRes-1;
    int nCellsY = grid->yRes-1;
    int nVertsX = nCellsX + 1;
    int nVertsY = nCellsY + 1;
    
    for (int i = 0; i<nVertsX; i++)
        for (int j = 0; j<nVertsY; j++){
            int xm1 = (i==0)?0:i-1;
            int xp1 = (i==nVertsX-1)?nVertsX-1:i+1;
            int ym1 = (j==0)?0:j-1;
            int yp1 = (j==nVertsY-1)?nVertsY-1:j+1;
            
            int xm1ym1 = (xm1 + ym1*nVertsX);
            int xp1ym1 = (xp1 + ym1*nVertsX);
            int xm1yp1 = (xm1 + yp1*nVertsX);
            int xp1yp1 = (xp1 + yp1*nVertsX);
            glm::vec3 v1 = glm::vec3(
                                     vertexXY[xp1yp1*2 + 0] - vertexXY[xm1ym1*2 + 0],
                                     vertexXY[xp1yp1*2 + 1] - vertexXY[xm1ym1*2 + 1],
                                     vertexZ[xp1yp1] - vertexZ[xm1ym1]
            );
            glm::vec3 v2 = glm::vec3(
                                     vertexXY[xm1yp1*2 + 0] - vertexXY[xp1ym1*2 + 0],
                                     vertexXY[xm1yp1*2 + 1] - vertexXY[xp1ym1*2 + 1],
                                     vertexZ[xm1yp1] - vertexZ[xp1ym1]
            );
            glm::vec3 normal = glm::cross(v1, v2);
            int index = (i + j*nVertsX)*3;
            for (int k=0; k<3; k++)
                normalData[index + k] = normal[k];
        }
}

static void MakeTriangleMesh(Sw_grid *grid){

    /*
        builds the triangle mesh
    */

    int nCellsX = grid->xRes-1;
    int nCellsY = grid->yRes-1;
    int nCellsTot = nCellsX*nCellsY;
    int nVertsX = nCellsX + 1;
    int nVertsY = nCellsY + 1;
    int nVertsTot = nVertsX*nVertsY;
    int nTriangles = nCellsTot * 2; // 2 triangles per cell
    int nIndices = nTriangles*3;

    vertexXY  = new GLfloat[nVertsTot*2];

    // make and bind the VAO
    glGenVertexArrays(1, &gVAO);
    glBindVertexArray(gVAO);
    
    for (int j = 0; j<nVertsY; j++)
    for (int i = 0; i<nVertsX; i++){
        int index = (i + j*nVertsX)*2;
        vertexXY[index + 0] = (GLfloat)i * grid->dx;
        vertexXY[index + 1] = (GLfloat)j * grid->dx;
    }
    
    GLuint *indexData = new GLuint[nIndices];
    for (int j = 0; j<nCellsY ; j++)
    for (int i = 0; i<nCellsX ; i++){
        // lower left triangle
        int index = (i+j*nCellsX)*6;
        int mm = i   + j*nVertsX;
        int pm = i+1 + j*nVertsX;
        int mp = i   + (j+1)*nVertsX;
        int pp = i+1 + (j+1)*nVertsX;
        indexData[index + 0] = mm;
        indexData[index + 1] = pm;
        indexData[index + 2] = mp;

        // upper right triangle
        indexData[index + 3] = pm;
        indexData[index + 4] = pp;
        indexData[index + 5] = mp;
    }
    
    glGenBuffers(1, &gVertexXY);
    glBindBuffer(GL_ARRAY_BUFFER, gVertexXY);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*nVertsTot*2, vertexXY, GL_STATIC_DRAW);
    
    glGenBuffers(1, &gIndex);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndex);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*nIndices, indexData, GL_STATIC_DRAW);

    // unbind the VAO
    glBindVertexArray(0); //CANNOT DO: Open Gl reports 'invalid operation'

    delete indexData;
}

// loads a cube into the VAO and VBO globals: gVAO and gVBO
static void LoadHeigthGPU(Sw_grid *grid) {

    int nCellsX = grid->xRes-1;
    int nCellsY = grid->yRes-1;
    int nVertsX = nCellsX + 1;
    int nVertsY = nCellsY + 1;
    int nVertsTot = nVertsX*nVertsY;

    // make and bind the VAO
    glGenVertexArrays(1, &gVAO);
    glBindVertexArray(gVAO);
    
    GLfloat *vertexZ = grid->oldFields[HEIGHT];

    GLfloat *vertexG = grid->oldFields[GROUND];

    // Computing Normals
    GLfloat *normalData = new GLfloat[nVertsTot*3];
    computeNormals(vertexXY, vertexZ, normalData, grid);
    
    glGenBuffers(1, &gVertexZ);
    glBindBuffer(GL_ARRAY_BUFFER, gVertexZ);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*nVertsTot, vertexZ, GL_STATIC_DRAW);

    glGenBuffers(1, &gVertexG);
    glBindBuffer(GL_ARRAY_BUFFER, gVertexG);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*nVertsTot, vertexG, GL_STATIC_DRAW);
    
    glGenBuffers(1, &gNormal);
    glBindBuffer(GL_ARRAY_BUFFER, gNormal);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*nVertsTot*3, normalData, GL_STATIC_DRAW);
    
    // unbind the VAO
    glBindVertexArray(0); //CANNOT DO: Open Gl reports 'invalid operation'

}

static void LoadBoxGPU(Box *box) {

//    int nTriangles = 12;
    int nBoxIndices = 12*3;

    // make and bind the VAO
    glGenVertexArrays(1, &gVAO);
    glBindVertexArray(gVAO);

    GLfloat *boxvertsXY = new GLfloat[8*2];
    GLfloat *boxvertsZ  = new GLfloat[8];
    for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
            for(int l=0; l<2; l++){
                Vector3f pos = box->x + (box->x0*pow(-1.f, i)+box->y0*pow(-1.f, j)+box->z0*pow(-1.f, l))*0.5f;
                boxvertsXY[(i + j*2 + l*4)*2 + 0] = pos.x();
                boxvertsXY[(i + j*2 + l*4)*2 + 1] = pos.y();
                boxvertsZ[  i + j*2 + l*4    + 0] = pos.z();
            }

    GLuint indexData0[] = {
        0,1,2,
        1,3,2,
        5,1,3,
        5,3,7,
        4,5,7,
        4,7,6,
        6,7,3,
        6,3,2,
        4,0,2,
        4,2,6,
        4,5,1,
        4,1,0
    };

    GLuint *indexData = new GLuint[nBoxIndices];
    for (int i = 0; i < nBoxIndices; ++i)
    {
        indexData[i] = indexData0[i];
    }

    glGenBuffers(1, &gBoxVertXY);
    glBindBuffer(GL_ARRAY_BUFFER, gBoxVertXY);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*8*2, boxvertsXY, GL_STATIC_DRAW);

    glGenBuffers(1, &gBoxVertZ);
    glBindBuffer(GL_ARRAY_BUFFER, gBoxVertZ);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*8, boxvertsZ, GL_STATIC_DRAW);

    glGenBuffers(1, &gBoxIndex);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gBoxIndex);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*nBoxIndices, indexData, GL_STATIC_DRAW);

    // unbind the VAO
    glBindVertexArray(0); //CANNOT DO: Open Gl reports 'invalid operation'

    // delete boxvertsXY, boxvertsZ;

}

// draws a single frame
static void Render(Sw_grid *grid) {

    int nCellsX = grid->xRes-1;
    int nCellsY = grid->yRes-1;
    int nCellsTot = nCellsX*nCellsY;
    int nTriangles = nCellsTot * 2; // 2 triangles per cell
    int nIndices = nTriangles*3;
    int nBoxIndices = 12*3;

    // clear everything
    glClearColor(0, 0, 0, 1); // black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // bind the program (the shaders)
    gProgram->use();
    glBindVertexArray(gVAO);

    // set the "camera" uniform
    gProgram->setUniform("camera", gCamera.matrix());
    
    // VERTEX XY
    glEnableVertexAttribArray(gProgram->attrib("vertXY"));
    glBindBuffer(GL_ARRAY_BUFFER, gVertexXY);
    glVertexAttribPointer(
                          gProgram->attrib("vertXY"),                  // attribute
                          2,                  // size = size of vector to pass to vertex shader
                          GL_FLOAT,           // type
                          GL_FALSE,           // normalized?
                          0,                  // stride
                          (void*)0            // array buffer offset
                          );
    
    // VERTEX Z = HEIGHT
    glEnableVertexAttribArray(gProgram->attrib("vertZ"));
    glBindBuffer(GL_ARRAY_BUFFER, gVertexZ);
    glVertexAttribPointer(
                          gProgram->attrib("vertZ"),                  // attribute
                          1,                  // size = size of vector to pass to vertex shader
                          GL_FLOAT,           // type
                          GL_FALSE,           // normalized?
                          0,                  // stride
                          (void*)0            // array buffer offset
                          );

    // INDEX
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndex);
    
    // NORMALS
    glEnableVertexAttribArray(gProgram->attrib("normal"));
    glBindBuffer(GL_ARRAY_BUFFER, gNormal);
    glVertexAttribPointer(
                          gProgram->attrib("normal"),                  // attribute
                          3,                  // size = size of vector to pass to vertex shader
                          GL_FLOAT,           // type
                          GL_FALSE,           // normalized?
                          0,                  // stride
                          (void*)0            // array buffer offset
                          );
    
    glDrawElements(
                   GL_TRIANGLES,      // mode
                   nIndices,          // count = number of indices
                   GL_UNSIGNED_INT,   // type
                   0                  // element array buffer offset
                   );

    glDisableVertexAttribArray(gProgram->attrib("vertZ"));
    glDisableVertexAttribArray(gProgram->attrib("normal"));

    // VERTEX G = HEIGHT
    glEnableVertexAttribArray(gProgram->attrib("vertZ"));
    glBindBuffer(GL_ARRAY_BUFFER, gVertexG);
    glVertexAttribPointer(
                          gProgram->attrib("vertZ"),                  // attribute
                          1,                  // size = size of vector to pass to vertex shader
                          GL_FLOAT,           // type
                          GL_FALSE,           // normalized?
                          0,                  // stride
                          (void*)0            // array buffer offset
                          );

    // INDEX
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndex);
    
    glDrawElements(
                   GL_TRIANGLES,      // mode
                   nIndices,          // count = number of indices
                   GL_UNSIGNED_INT,   // type
                   0                  // element array buffer offset
                   );
    
    
    glDisableVertexAttribArray(gProgram->attrib("vertZ"));
    glDisableVertexAttribArray(gProgram->attrib("vertXY"));

    /////////////////// BOX //////////////
    glEnableVertexAttribArray(gProgram->attrib("vertXY"));
    glBindBuffer(GL_ARRAY_BUFFER, gBoxVertXY);
    glVertexAttribPointer(
                          gProgram->attrib("vertXY"),                  // attribute
                          2,                  // size = size of vector to pass to vertex shader
                          GL_FLOAT,           // type
                          GL_FALSE,           // normalized?
                          0,                  // stride
                          (void*)0            // array buffer offset
                          );

    glEnableVertexAttribArray(gProgram->attrib("vertZ"));
    glBindBuffer(GL_ARRAY_BUFFER, gBoxVertZ);
    glVertexAttribPointer(
                          gProgram->attrib("vertZ"),                  // attribute
                          1,                  // size = size of vector to pass to vertex shader
                          GL_FLOAT,           // type
                          GL_FALSE,           // normalized?
                          0,                  // stride
                          (void*)0            // array buffer offset
                          );

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gBoxIndex);
    
    glDrawElements(
                   GL_TRIANGLES,      // mode
                   nBoxIndices,          // count = number of indices
                   GL_UNSIGNED_INT,   // type
                   0                  // element array buffer offset
                   );

    glDisableVertexAttribArray(gProgram->attrib("vertZ"));
    glDisableVertexAttribArray(gProgram->attrib("vertXY"));

    // unbind the VAO, the program
    glBindVertexArray(0);
    gProgram->stopUsing();
    
    // swap the display buffers (displays what was just drawn)
    glfwSwapBuffers();
}


// update the scene based on the time elapsed since last update
void Update(float secondsElapsed, bool &bSolving, bool &bRain, bool &bTsunami, bool &bBlob, bool &bBeachRise, bool &bMoveBoxDown, bool &bMoveBoxUp) {

    //move position of camera based on WASD keys, and XZ keys for up and down
    const float moveSpeed = 2.0; //units per second
    if(glfwGetKey('S')){
        gCamera.offsetPosition(secondsElapsed * moveSpeed * -gCamera.forward());
    } else if(glfwGetKey('W')){
        gCamera.offsetPosition(secondsElapsed * moveSpeed * gCamera.forward());
    }
    if(glfwGetKey('A')){
        gCamera.offsetPosition(secondsElapsed * moveSpeed * -gCamera.right());
    } else if(glfwGetKey('D')){
        gCamera.offsetPosition(secondsElapsed * moveSpeed * gCamera.right());
    }
    if(glfwGetKey('Z')){
        gCamera.offsetPosition(secondsElapsed * moveSpeed * -glm::vec3(0,1,0));
    } else if(glfwGetKey('X')){
        gCamera.offsetPosition(secondsElapsed * moveSpeed * glm::vec3(0,1,0));
    }

    if(glfwGetKey('P')){
        bSolving = true;
    }
    if(glfwGetKey('L')){
        bSolving = false;
    }
    if(glfwGetKey('R')){
        bRain = true;
    }
    if(glfwGetKey('T')){
        bTsunami = true;
    }
    if(glfwGetKey('B')){
        bBeachRise = true;
    }

	if(glfwGetKey('H')){
        bMoveBoxDown = true;
    }	
	
	if(glfwGetKey('U')){
        bMoveBoxUp = true;
    }

    //rotate camera based on mouse movement
    const float mouseSensitivity = 0.1;
    int mouseX, mouseY;
    glfwGetMousePos(&mouseX, &mouseY);
    gCamera.offsetOrientation(mouseSensitivity * mouseY, mouseSensitivity * mouseX);
    glfwSetMousePos(0, 0); //reset the mouse, so it doesn't go out of the window

    //increase or decrease field of view based on mouse wheel
    const float zoomSensitivity = -0.2;
    float fieldOfView = gCamera.fieldOfView() + zoomSensitivity * (float)glfwGetMouseWheel();
    if(fieldOfView < 5.0f) fieldOfView = 5.0f;
    if(fieldOfView > 130.0f) fieldOfView = 130.0f;
    gCamera.setFieldOfView(fieldOfView);
    glfwSetMouseWheel(0);
}

class FPSMeasure {
    int nbFrames;
    float lastTime;
    float printFrequency;
    
public:
    FPSMeasure(float _printFrequency = 1.0):nbFrames(0), printFrequency(_printFrequency){
        lastTime = glfwGetTime();
    }
    
    void measure_fps(){
        double currentTime = glfwGetTime();
        nbFrames++;
        if ( currentTime - lastTime >= printFrequency ){ // If last prinf() was more than 1 sec ago
            // printf and reset timer
            float fps = double(nbFrames)/printFrequency;
            printf("%f fps\n", fps);
            nbFrames = 0;
            lastTime += printFrequency;
        }
    };
};

void initialConditions(Sw_grid *grid){
    const int xRes = grid->xRes;
    const int yRes = grid->yRes;
    
    float *eta    = grid->oldFields[ETA];
    float *height = grid->oldFields[HEIGHT];
    float *ground = grid->oldFields[GROUND];
    float *velx   = grid->oldFields[VELX];
    float *vely   = grid->oldFields[VELY];
    float *neweta    = grid->newFields[ETA];
    float *newheight = grid->newFields[HEIGHT];
    float *newground = grid->newFields[GROUND];
    float *newvelx   = grid->newFields[VELX];
    float *newvely   = grid->newFields[VELY];

    for (int i = 0; i < xRes; ++i)
    for (int j = 0; j < yRes; ++j){
        int index = i + j*xRes;
        //if (i > xRes/4 && i < xRes*3/4 && j > yRes/4 && j < yRes*3/4)
        //    eta[index] = 0.2;
        //else
        eta[index] = 1;
        ground[index] = 0.0;
        height[index] = ground[index] + eta[index];
        velx[index] = 0.0;
        vely[index] = 0.0;
        
        neweta   [index] = 0.0;
        newheight[index] = 0.0;
        newground[index] = 0.0;
        newvelx  [index] = 0.0;
        newvely  [index] = 0.0;
    }
}

void initialConditionBeach(Sw_grid *grid){
    const int xRes = grid->xRes;
    const int yRes = grid->yRes;
    
    float *eta    = grid->oldFields[ETA];
    float *height = grid->oldFields[HEIGHT];
    float *ground = grid->oldFields[GROUND];
    float *velx   = grid->oldFields[VELX];
    float *vely   = grid->oldFields[VELY];

    for (int i = 0; i < xRes; ++i)
    for (int j = 0; j < yRes; ++j){
        int index = i + j*xRes;
        
        if (j < yRes*1/4)
            eta[index] = 0.1;
        else
            eta[index] = 0.0;

        if (j > yRes*3/4)
            ground[index] = 0.2 * (float)(j-yRes*3/4)*grid->dx;
        else
            ground[index] = 0.0;
        height[index] = ground[index] + eta[index];
        velx[index] = 0.0;
        vely[index] = 0.0;
    }
}

/* SPECIAL EFFECTS */

void rain(Sw_grid *grid){

    int ii = rand() % grid->xRes;
    int jj = rand() % grid->yRes;
    float *eta = grid->oldFields[ETA];
//    float *height = grid->oldFields[HEIGHT];
    float drop = 0.05;
    for (int i = -1; i < 2; ++i)
    for (int j = -2; j < 2; ++j) {
            eta[i+ii + (j+jj)*grid->xRes] += drop;
    }
    
}

void tsunami(Sw_grid *grid){
    float *eta = grid->oldFields[ETA];
    float *height = grid->oldFields[HEIGHT];
//    float *vely = grid->oldFields[VELY];
    for (int i = 0; i < grid->xRes; ++i)
    for (int j = 0; j < grid->yRes/10; ++j){
        // vely[   i + j*grid->xRes] += 0.01;
        height[i + j*grid->xRes] += (1.0 - (float)((float)j/(float)grid->yRes*10.0))*0.01;
        eta[i + j*grid->xRes] += (1.0 - (float)((float)j/(float)grid->yRes*10.0))*0.01;
    }
    printf("TSUNAMI!!!\n");
}

void beachRise(Sw_grid *grid){
    float *ground = grid->oldFields[GROUND];
    for (int i = 0; i < grid->xRes; ++i)
    for (int j = 0; j < grid->yRes/4; ++j){
        ground[i + j*grid->xRes] += (1.0 - (float)((float)j/(float)grid->yRes*4.0))*0.01;
    }
}

void beachRise2(Sw_grid *grid){
    float *ground = grid->oldFields[GROUND];
    for (int i = 0; i < grid->xRes; ++i)
    for (int j = grid->yRes*3/4; j < grid->yRes; ++j){
        ground[i + j*grid->xRes] += (float)(j-grid->yRes*3/4)/(float)(grid->yRes/4)*0.01;
    }
}

void blob(Sw_grid *grid){
    float *eta = grid->oldFields[ETA];
//    float *height = grid->oldFields[GROUND];
    int xRes = grid->xRes;
    int yRes = grid->yRes;
    for (int i = 0; i < xRes; ++i)
    for (int j = 0; j < yRes; ++j){
        // if (i > xRes*4/10 && i < xRes*6/10 && j > yRes*4/10 && j < yRes*6/10)
            eta[   i + j*xRes] -= 0.01;
    }
}

void initGLFW(){
    // initialise GLFW
    if(!glfwInit())
        throw std::runtime_error("glfwInit failed");

    // open a window with GLFW
    glfwOpenWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 3);
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 2);
    glfwOpenWindowHint(GLFW_WINDOW_NO_RESIZE, GL_TRUE);
    if(!glfwOpenWindow((int)SCREEN_SIZE.x, (int)SCREEN_SIZE.y, 8, 8, 8, 8, 16, 0, GLFW_WINDOW))
        throw std::runtime_error("glfwOpenWindow failed. Can your hardware handle OpenGL 3.2?");

    // GLFW settings
    glfwDisable(GLFW_MOUSE_CURSOR);
    glfwSetMousePos(0, 0);
    glfwSetMouseWheel(0);

    // initialise GLEW
    glewExperimental = GL_TRUE; //stops glew crashing on OSX :-/
    if(glewInit() != GLEW_OK)
        throw std::runtime_error("glewInit failed");

    // GLEW throws some errors, so discard all the errors so far
    while(glGetError() != GL_NO_ERROR) {}

    // print out some info about the graphics drivers
    std::cout << "OpenGL version: " << glGetString(GL_VERSION) << std::endl;
    std::cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
    std::cout << "Vendor: " << glGetString(GL_VENDOR) << std::endl;
    std::cout << "Renderer: " << glGetString(GL_RENDERER) << std::endl;

    // make sure OpenGL version 3.2 API is available
    if(!GLEW_VERSION_3_2)
        throw std::runtime_error("OpenGL 3.2 API is not available.");

    // OpenGL settings
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

}

// the program starts here
void AppMain() {
    float size = 20;
    Sw_grid *grid = new Sw_grid(size, size, 1.0/size);

    // initialConditionBeach(grid);
    initialConditions(grid);
    Box *box = new Box(Vector3f(0.5f, 0.5f, 1.0f), 0.1, Vector3f(0.2f, 0.0f, 0.0f),  Vector3f(0.0f, 0.2f, 0.0f),  Vector3f(0.0f, 0.0f, 0.2f));

    // initiate solver and bind grid
    SWSolver sw_solver;
    sw_solver.setGrid(grid);
    SWRBSolver swrb_solver(grid, box);

    // initialize GLFW & GLEW
    initGLFW();

    // load vertex and fragment shaders into opengl
    LoadShaders();

    // setup gCamera
    gCamera.setPosition(glm::vec3(2,2,2));
    gCamera.setViewportAspectRatio(SCREEN_SIZE.x / SCREEN_SIZE.y);
    gCamera.lookAt(glm::vec3(0,0,0));
    
    // create buffer and fill it with the points of the triangle
    MakeTriangleMesh(grid);
    
    FPSMeasure fps;
    float sim_time = 0.0;
    float phys_time = 0.0;
    bool bSolving   = false;
    bool bRain      = false;
    bool bTsunami   = false;
    bool bBlob      = false;
    bool bBeachRise = false;
	bool bMoveBoxDown = false;
	bool bMoveBoxUp = false;

    // run while the window is open
    double lastTime = glfwGetTime();
    while(glfwGetWindowParam(GLFW_OPENED)){

        // update the scene based on the time elapsed since last update
        double thisTime = glfwGetTime();
        float dt = thisTime - lastTime;
        Update(dt, bSolving, bRain, bTsunami, bBlob, bBeachRise, bMoveBoxDown, bMoveBoxUp);
        
        // if((int)(thisTime*10)%2 == 0){
        //     float potEnergy = solver.computePotentialEnergy();
        //     float kinEnergy = solver.computeKineticEnergy();
        //     printf("Kinetic Energy = %f, potEnergy = %f, total Energy = %f\n", kinEnergy, potEnergy, kinEnergy+potEnergy);
        // }

        if (bRain) {
            rain(grid);
            bRain = false;
        }
        if (bTsunami) {
            tsunami(grid);
            bTsunami = false;
        }
        if (bBeachRise) {
            beachRise(grid);
            // beachRise2(grid);
            bBeachRise = false;
        }
        if (bBlob) {
            blob(grid);
            bBlob = false;
        }

		if (bMoveBoxDown) {
			box->x -= Vector3f(0.0f, 0.0f, 0.01f);
			bMoveBoxDown = false;
		}

		if (bMoveBoxUp) {
			box->x += Vector3f(0.0f, 0.0f, 0.01f);
			bMoveBoxUp = false;
		}

        if(bSolving){
            sim_time += dt;

            float *eta = grid->oldFields[ETA];
            int xRes = grid->xRes;
            int yRes = grid->yRes;
            float H = 0.;
            for (int i = 0; i < xRes - 1; i++)
            for (int j = 0; j < yRes - 1; j++)
                H += eta[i+j*xRes];
            H /= (float)(xRes*yRes);
            float maxdt = 0.001;//0.1*grid->dx / sqrt(sw_solver.g*H);
            
            //while(phys_time < sim_time){
                
                float dtPhys = std::min(maxdt, sim_time - phys_time);
//                sw_solver.advanceTimestep(dtPhys);
                swrb_solver.advanceTimestep(dtPhys);
                phys_time += dtPhys;
            //}
			//bSolving = false;
        }
        lastTime = thisTime;

        // Copy height field to GPU and compute normals
        LoadHeigthGPU(grid);
        LoadBoxGPU(box);

        // draw one frame
        Render(grid);

        // check for errors
        GLenum error = glGetError();
        if(error != GL_NO_ERROR)
            std::cerr << "OpenGL Error " << error << ": " << (const char*)gluErrorString(error) << std::endl;

        // measure frames per second
        fps.measure_fps();
        
        //exit program if escape key is pressed
        if(glfwGetKey(GLFW_KEY_ESC))
            glfwCloseWindow();
    }

    // clean up and exit
    glfwTerminate();
}


int main(int argc, char *argv[]) {
    try {
        AppMain();
    } catch (const std::exception& e){
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
