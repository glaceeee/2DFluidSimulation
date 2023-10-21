#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <windows.h>
#include <string>

#define IX(i,j) ((i) + (N+2)*(j)) //index
#define SWAP(x0,x) {double* tmp = x0; x0 = x; x = tmp;} //swap 2 pointers
#define LERP(val1,val2,x) (val1 + x*(val2-val1)) //linear interpolate macro used in advect

/* Creates and initializes the passed GLFWwindow object, also calls glewInit() and makes sure it is GLEW_OK */
void initializeWindow(GLFWwindow*& window, int width, int height, const char* name) {
    if (!glfwInit()) {
        std::cout << "glfwInit() didn't work" << std::endl;
        return;
    }
    
    window = glfwCreateWindow(width, height, name, NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        std::cout << "passed GLFWwindow* reference doesn't exist" << std::endl;
        return;
    }

    glfwMakeContextCurrent(window);

    /* initialize glew, so we can use all the other openGL functions that were introduced past version 1.1 */
    if (glewInit() != GLEW_OK) {
        std::cout << "glewInit != GLEW_OK" << std::endl;
        return;
    }  
}

/* Initializes all the necessary arrays used in the simulation */
void initializeArrays(double*& u, double*& u_prev, double*& v, double*& v_prev, double*& dens, double*& dens_prev, const int& size) {
    u = new double[size];
    u_prev = new double[size];
    v = new double[size];
    v_prev = new double[size];
    dens = new double[size];
    dens_prev = new double[size];

    for (int i = 0; i < size; i++) {
        u[i] = 0;
        u_prev[i] = 0;
        v[i] = 0;
        v_prev[i] = 0;
        dens[i] = 0;
        dens_prev[i] = 0;
    }
}

/* Goes through inputArray of certain values, i.e. velocities, densities, and puts those in the given valueArray, boundary cells included*/
void addSources(double* x, double* sourceArray, int size, double dt) {
    for (int i = 0; i < size; i++) {
        x[i] += dt*sourceArray[i];
    }
}

/* Sets the boundaries. Free-slip condition */
void set_bnd(int N, int b, double* x) {
    for (int i = 1; i <= N; i++) {
        x[IX(0, i)] = (b == 1) ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = (b == 1) ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i, 0)] = (b == 2) ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N + 1)] = (b == 2) ? -x[IX(i, N)] : x[IX(i, N)];
    }
    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(N + 1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N+1, 1)]);
    x[IX(0, N + 1)] = 0.5 * (x[IX(0, N)] + x[IX(1, N+1)]);
    x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N+1, N)]);
}

/* Handles diffusion. In one time interval the passed quantity should gradually become the same as its 4 adjacent cells. */
void diffuse(double* x, double* x0, double diff, int N, int b, double dt) {
    double a = dt*diff*N*N; //Controls rate of diffusion. The bigger, the quicker things diffuse. Multiplication by dt*N*N to convert units
    for (int k = 0; k < 20; k++) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                /* This equation makes sure of multiple things:
                    1) The relations between x0(i,j) and x(i,j) are constant such that if x0(i,j) > x(i,j), x(i,j) < x0(i,j) and vice versa.
                       This works for all values of a, since this equation finds a value for x(i,j) that it needs to remain accurate to its'
                       surroundings, based on the current values of the surrounding cells and a and its' relation to x0(i,j).
                    2) The cells' updated values influence each other and they adjust each other based on the others' needs.
                       In other words: cell(i,j) updates (now cell(i,j)'s value is correct from its point of view) 
                       => cell(i,j+1) updates next based on the value of cell(i,j). Now cell(i,j+1) updated its value so its right from its point of view.
                       now next time we get to cell(i,j), cell(i,j) will see the update in value in cell(i,j+1), that it, itself caused
                       among other cells and update itself so its value fits from its point of view again. Now cell(i,j+1) needs to update its value again
                       so it fits from its point of view again. It does this until at some point the values start to converge, i.e. cell(i,j)
                       is more and more satisfied with cell(i,j+1)'s value and doesnt need to change itself much, and vice versa.
                    3) That the process converges. If you were to determine the matrix for all the N*N unknowns in this system. One unknown for each cell's value.
                       Then what you would see is that in every calculation for every cell the coefficient of x[IX(i,j)], that we are solving for here,
                       is (1+4*a), while the 4 adjacent cells' coefficients are a. The remaining unknowns' coefficients are 0.
                       Thus the matrix is diagonally dominant, since (1+4*a) > a+a+a+a+0+0+0+...+0. Therefore, since we're using gauss-seidel to determine
                       the solution to this system of linear equations, this system will always converge, and we will always get the values we're looking for. */
                x[IX(i,j)] = (x0[IX(i,j)] + a * (x[IX(i+1,j)] + x[IX(i-1,j)] + x[IX(i,j+1)] + x[IX(i,j-1)])) / (1+4*a);
            }
        }
        /* When we're diffusing densities, this line is irrelevant, however when diffusing velocities it is crucial.
           This line makes sure that whenever we update a cell, that is next to a boundary cell,'s velocity,
           we immediately reflect that cell's velocity on the boundary cell, so that the boundary cell can now,
           in the next gauss-seidel iteration, steer the velocity of the cell away from the boundary, so the cell points more in the
           boundaries tangential direction. */
       set_bnd(N, b, x);
    }
}

/* Moves the fluid along its vectors (roughly), using semi-lagranian advection */
void advect(double* x, double* x0, double* u, double* v, int N, int b, double dt) {
   double particleX, particleY, particleX_INT, particleY_INT, particleX_FRACT, particleY_FRACT, lerpX1, lerpX2, interpolatedValue;

    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {

            /* Determine "particle's" position */
            particleX = i - u[IX(i,j)] * dt * N;
            particleY = j - v[IX(i,j)] * dt * N;
            /* Handles overflow by capping the max value of particleX and particleY at the respective row/columns boundary cell.
               This way no matter how big the velocity at a certain point is, particleX will always be in the same row and particleY
               will never be out of bounds. */
            if(particleX < 0.5) particleX = 0.5; if(particleX > N+0.5) particleX = N+0.5; //handle x-overflow
            if(particleY < 0.5) particleY = 0.5; if(particleY > N+0.5) particleY = N+0.5; //handle y-overflow
            particleX_FRACT = modf(particleX, &particleX_INT);
            particleY_FRACT = modf(particleY, &particleY_INT);
            particleY_INT++;

            // Interpolate the values of 4 cells to approximate the particle's value
            lerpX1 = LERP(x0[IX(int(particleX_INT), int(particleY_INT))], x0[IX(int(particleX_INT+1), int(particleY_INT))], particleX_FRACT);
            lerpX2 = LERP(x0[IX(int(particleX_INT), int(particleY_INT-1))], x0[IX(int(particleX_INT+1), int(particleY_INT)-1)], particleX_FRACT);
            x[IX(i,j)] = LERP(lerpX2, lerpX1, particleY_FRACT);

            particleX = 0, particleY = 0, particleX_INT = 0, particleY_INT = 0, particleX_FRACT = 0, particleY_FRACT = 0, lerpX1 = 0, lerpX2 = 0, interpolatedValue = 0;
        }
    }
    /* vector field changed thus update boundary values */
    set_bnd(N, b, x);
}

/* Clears any divergence and thus forces incompressibility */
void project(int N, double* u, double* v, double* p, double* div) {
    int i, j, k;
    double h = 1.0/N; //partly here so that things converge better/faster/more when computing the scalar field p. also here to keep units consistent, but i dont get that part tbh

    /* First calculate divergence with discrete formula */
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
            div[IX(i, j)] = 0.5*(u[IX(i+1, j)] - u[IX(i-1, j)] + v[IX(i, j+1)] - v[IX(i, j-1)]);
            p[IX(i, j)] = 0;
        }
    }

    /* Make divergence consistent along boundaries, i.e. just copy over the values of the already present values to the boundary cells 
       and also initialize all boundary cells in p to 0 */
    set_bnd(N, 0, div); set_bnd(N, 0, p);

    /* Now using Gauss-Seidel again calculate a scalar-field p, such that del*del*p = del*currentVelocityField
       i.e. such that div(del*p) = div(currentVelocityField). The gradient of p has no curl and exactly the
       divergence div(currentVelocityField). Thus since del*p is irrotational and has the same exact divergence as
       currentVelocityField at every point (x,y) in a continous context and (i,j) in a discrete context and is therefore
       the exact irrotational vector field we need to subtract from currentVelocityField to obtain a divergence free vector field
       and keep our fluid incompressible and mass conservant. This will also converge since, again, the matrix of this system is
       diagonally dominant. Definition of diagonal dominance: abs(aii) = sum_of_all(abs(aij)). */
    for (k = 0; k < 20; k++) {
        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                p[IX(i, j)] = (p[IX(i-1, j)] + p[IX(i+1, j)] + p[IX(i, j-1)] + p[IX(i, j+1)] - div[IX(i, j)]) / 4;
            }
        }
        /* as approximation let the boundary cells be the same value as the cells that surround them, better than just having them all be 0
           more consistency and also more expected behaviour */
        set_bnd(N, 0, p);
    }

    /* Now calculate the gradient (i.e. del*p) of our newly determined scalar field p */
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
            u[IX(i, j)] -= 0.5 * (p[IX(i+1, j)] - p[IX(i-1, j)]);
            v[IX(i, j)] -= 0.5 * (p[IX(i, j+1)] - p[IX(i, j-1)]);
        }
    }

    /* vector field changed thus update boundary values accordingly */
    set_bnd(N, 1, u); set_bnd(N, 2, v);
}

/* Advanced the velocites by one timestep */
void velocity_step(int N, double dt, double diff, double* u, double* u_prev, double* v, double* v_prev) {
    addSources(u, u_prev, (N+2) * (N+2), dt); addSources(v, v_prev, (N+2) * (N+2), dt);
    SWAP(u_prev, u); diffuse(u, u_prev, diff, N, 1, dt); 
    SWAP(v_prev, v); diffuse(v, v_prev, diff, N, 2, dt);
    project(N, u, v, u_prev, v_prev);
    SWAP(u_prev, u); SWAP(v_prev, v); 
    advect(u, u_prev, u_prev, v_prev, N, 1, dt); advect(v, v_prev, u_prev, v_prev, N, 2, dt);
    project(N, u, v, u_prev, v_prev);
}

/* Advances the densities by one timestep */
void density_step(int N, double dt, double diff, double* dens, double* dens_prev, double* u, double* v) {
    addSources(dens, dens_prev, (N+2) * (N+2), dt);
    SWAP(dens_prev, dens); 
    diffuse(dens, dens_prev, diff, N, 0, dt);
    SWAP(dens_prev, dens); 
    advect(dens, dens_prev, u, v, N, 0, dt);
}

/* Draws the grid based on the values in the passed 'densityArray'. The densityArray includes boundary cells. It is up to the caller to adjust the array to their need */
void drawGrid(double* densityArray, int N, double gridSize) {
    int arrayIndex;
    double gridSizeTimes2 = gridSize * 2;
    double x = -1.0;
    double y = -1.0;

    glBegin(GL_QUADS);
    for (int j = 1; j < N+1; j++) {
        for (int i = 1; i < N+1; i++) {
            arrayIndex = IX(i,j);
            glColor3d(densityArray[arrayIndex], densityArray[arrayIndex], densityArray[arrayIndex]);
            glVertex2d(x, y);
            glVertex2d(x + gridSizeTimes2, y);
            glVertex2d(x + gridSizeTimes2, y + gridSizeTimes2);
            glVertex2d(x, y + gridSizeTimes2);
            x += gridSizeTimes2;
        }
        x = -1.0;
        y += gridSizeTimes2;
    }
    glEnd();
}

/* Handles user interface. Mouse input: Left & drag to add sources, Right & drag to move fluid (i.e. only add velocities) */
// TO-DO: mouse movement in negative part of coordinate system, no reflection but same as with positive direcion 
void getSourceValuesFromUI(int N, double dt, int width, double dyeIntensity, int brushSize, double* dens_prev, double* u_prev, double* v_prev, double x_prev, double y_prev, double x, double y, GLFWwindow* window) {
    double h = 1.0 / (double)N;
    bool leftMousePressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    bool rightMousePressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
    int screenX_prev = (int)(x_prev * ((double)N/(double)width));
    int screenY_prev = (int)((width - y_prev) * ((double)N / (double)width));
    int screenX = (int)(x * ((double)N/(double)width));
    int screenY = (int)((width - y) * ((double)N / (double)width));
    int signX = (0 < (screenX - screenX_prev)) - (0 > (screenX - screenX_prev));
    int signY = (0 < (screenY - screenY_prev)) - (0 > (screenY - screenY_prev));
    int stepsInX;
    stepsInX = (signX > 0 && screenX - screenX_prev < 0) ? abs(screenX - screenX_prev) + N + 1 : screenX - screenX_prev + 1;
   // stepsInX = (signX < 0 && screenX - screenX_prev > 0) ? screenX - screenX_prev + 1 : screenX - screenX_prev + 1;
    int x_step = screenX_prev; //geometry dash reference?!?!
    int y_step = screenY_prev;
    int x_stepSave = x_step;
    int y_stepSave = y_step;
    //double slope = ((double)(abs(screenX - screenX_prev)) != 0) ? ((double)(screenY - screenY_prev) / (double)(abs(screenX - screenX_prev))) : 3;
    double slope = ((double)(screenY - screenY_prev) / (double)(abs(screenX - screenX_prev)));

    /* Make fluid loop from one end of a boundary too another. This does not affect the free-slip boundary condition of our actual boundary cells
       and is more here to make moving your cursor along the screen easier and more satisfying. For example once your cursor hits the left most
       boundary we'll move the user's cursor to x=200*/
    //screenX = (screenX > N) ? (screenX % N) : screenX;
    //screenX = (screenX < 1) ? (screenX % N)+N : screenX;
    //screenY = (screenY > N) ? (screenY % N) : screenY;
    //screenY = (screenY < 1) ? (screenY % N)+N : screenY;

    /* Clear the arrays so we can add new source values to them */
    for (int i = 0; i < (N + 2) * (N + 2); i++) {
        u_prev[i] = 0;
        v_prev[i] = 0;
        dens_prev[i] = 0;
    }

    ///* If left mouse button pressed, add density AND velocity (along linear path of cursor) */
    //if (leftMousePressed) {
    //    /* Handle vertical interpolation (i.e. slopes of inf), cursor movement is considered purely vertical at angle of roughly 0 with a given tolerance */
    //    if (isinf(slope)) {
    //        std::cout << "slope is inf" << std::endl;
    //        while (y_step != screenY+signY) {
    //            v_prev[IX(screenX, y_step)] = ((double)(screenY - screenY_prev)*h) / dt;
    //            dens_prev[IX(screenX, y_step)] = dyeIntensity;
    //            y_step += signY;
    //        }
    //    }
    //    /* Handle non-vertical interpolation */
    //    else {
    //        /* Only check for x_step, to ensure that angles of 0 (i.e. horizontal cursor movements) are also drawn/recorded */
    //        for(int i = 0; i < stepsInX; i++) {
    //           // for (int i = 0; i < brushSize; i++) {
    //           //     if (x_step + i <= N + 1 || y_step + i <= N+1)
    //                u_prev[IX(x_step, y_step)] = ((double)(signX*(stepsInX-1))*h) / dt;
    //                v_prev[IX(x_step, y_step)] = ((double)(screenY - screenY_prev)*h) / dt;
    //                dens_prev[IX(x_step, y_step)] = dyeIntensity;
    //                //std::cout << "u:" << std::endl;
    //                //std::cout << u_prev[IX(x_step, y_step)] << std::endl;
    //                //std::cout << "v:" << std::endl;
    //                //std::cout << v_prev[IX(x_step, y_step)] << std::endl;
    //                x_step += signX;
    //                x_step = (x_step > N) ? (x_step % N) : x_step;
    //                x_step = (x_step < 1) ? (x_step % N) + N : x_step;
    //           //}

    //            if (x_step != stepsInX) {
    //                /* Always advance one step in direction set by signX and then follow with the necessary amount of steps in the y-direction. More below */
    //                    y_stepSave = y_step;
    //                    //x_stepSave = x_step;

    //                /* Handle cases where x_step is in between screenX and screenX_prev but the two variables are on opposite horizontal sides of the screen */
    //                    switch(signX){
    //                    /* 1) Cursor went from right to left? */
    //                    case(1):
    //                        y_step = (x_step < screenX_prev && x_step < screenX)
    //                            ? std::floor(slope * (abs(x_step + N - screenX_prev))) + screenY_prev
    //                            : std::floor(slope * (abs(x_step - screenX_prev))) + screenY_prev;
    //                        break;
    //                    /* 2) Cursor went from left to right */
    //                   /* case(-1):
    //                        y_step = (x_step < screenX_prev && x_step < screenX)
    //                            ? std::floor(slope * (abs(x_step - N - screenX_prev))) + screenY_prev
    //                            : std::floor(slope * (abs(x_step - screenX_prev))) + screenY_prev;
    //                        break;*/
    //                    default:
    //                        std::cout << "Invalid signX" << std::endl;
    //                        break;
    //                    }

    //                    y_step = std::floor(slope * (abs(x_step - screenX_prev))) + screenY_prev;

    //                    //if (x_step < screenX_prev && x_step < screenX) { y_step = std::floor(slope * (abs(x_step + N - screenX_prev))) + screenY_prev; }
    //                    ///* 2) Cursor went from left to right? */
    //                    //else if (x_step > screenX_prev && x_step > screenX) { y_step = std::floor(slope * (abs(x_step - N - screenX_prev))) + screenY_prev; }
    //                    ///* 3) Cursor is didn't switch walls */
    //                    //else { y_step = std::floor(slope * (abs(x_step - screenX_prev))) + screenY_prev; }

    //                /* Handle slopes of pi/4 (45°) < gradient < inf */
    //                    if (abs(y_step - y_stepSave) > 1) {
    //                        for (int y_step_step = 0; y_step_step < abs(y_step - y_stepSave)-1; y_step_step++) {
    //                            if (signY > 0) { 
    //                                u_prev[IX(x_step, y_stepSave + y_step_step)] = ((double)(screenX - screenX_prev) * h) / dt;
    //                                v_prev[IX(x_step, y_stepSave + y_step_step)] = ((double)(screenY - screenY_prev) * h) / dt;
    //                                dens_prev[IX(x_step, y_stepSave + y_step_step)] = dyeIntensity; 
    //                            }
    //                            if (signY < 0) { 
    //                                u_prev[IX(x_step, y_stepSave - y_step_step)] = ((double)(screenX - screenX_prev) * h) / dt;
    //                                v_prev[IX(x_step, y_stepSave - y_step_step)] = ((double)(screenY - screenY_prev) * h) / dt;
    //                                dens_prev[IX(x_step, y_stepSave - y_step_step)] = dyeIntensity; 
    //                            }
    //                        }
    //                    }
    //            }
    //        }
    //    }
    //}

    /* If right mouse button pressed, add density AND velocity (along linear path of cursor) */
    if (leftMousePressed) {
        /* Handle vertical interpolation (i.e. slopes of inf), cursor movement is considered purely vertical at angle of roughly 0 with a given tolerance */
        if (isinf(slope)) {
            while (y_step != screenY + signY) {
                v_prev[IX(screenX, y_step)] = ((double)(screenY - screenY_prev) * h) / dt;
                dens_prev[IX(screenX, y_step)] = dyeIntensity;
                y_step += signY;
            }
        }
        /* Handle non-vertical interpolation */
        else {
            /* Only check for x_step, to ensure that angles of 0 (i.e. horizontal cursor movements) are also drawn/recorded */
            while (x_step != screenX + signX) {
                u_prev[IX(x_step, y_step)] = ((double)(screenX - screenX_prev) * h) / dt;
                v_prev[IX(x_step, y_step)] = ((double)(screenY - screenY_prev) * h) / dt;
                dens_prev[IX(x_step, y_step)] = dyeIntensity;
                x_step += signX;

                if (x_step != screenX + signX) {
                    /* Always advance one step in direction set by signX and then follow with the necessary amount of steps in the y-direction. More below */
                    y_stepSave = y_step;
                    y_step = std::floor(slope * (abs(x_step - screenX_prev))) + screenY_prev;

                    /* Handle slopes of pi/4 (45°) < gradient < inf */
                    if (abs(y_step - y_stepSave) > 1) {
                        for (int y_step_step = 0; y_step_step < abs(y_step - y_stepSave) - 1; y_step_step++) {
                            if (signY > 0) {
                                u_prev[IX(x_step, y_stepSave + y_step_step)] = ((double)(screenX - screenX_prev) * h) / dt;
                                v_prev[IX(x_step, y_stepSave + y_step_step)] = ((double)(screenY - screenY_prev) * h) / dt;
                                dens_prev[IX(x_step, y_stepSave + y_step_step)] = dyeIntensity;
                            }
                            if (signY < 0) {
                                u_prev[IX(x_step, y_stepSave - y_step_step)] = ((double)(screenX - screenX_prev) * h) / dt;
                                v_prev[IX(x_step, y_stepSave - y_step_step)] = ((double)(screenY - screenY_prev) * h) / dt;
                                dens_prev[IX(x_step, y_stepSave + y_step_step)] = dyeIntensity;
                            }
                        }
                    }
                }
            }
        }
    }

    /* If right mouse button pressed, add density AND velocity (along linear path of cursor) */
    if (rightMousePressed) {
        /* Handle vertical interpolation (i.e. slopes of inf), cursor movement is considered purely vertical at angle of roughly 0 with a given tolerance */
        if (isinf(slope)) {
            while (y_step != screenY + signY) {
                v_prev[IX(screenX, y_step)] = ((double)(screenY - screenY_prev) * h) / dt;
                y_step += signY;
            }
        }
        /* Handle non-vertical interpolation */
        else {
            /* Only check for x_step, to ensure that angles of 0 (i.e. horizontal cursor movements) are also drawn/recorded */
            while (x_step != screenX + signX) {
                u_prev[IX(x_step, y_step)] = ((double)(screenX - screenX_prev) * h) / dt;
                v_prev[IX(x_step, y_step)] = ((double)(screenY - screenY_prev) * h) / dt;
                x_step += signX;

                if (x_step != screenX + signX) {
                    /* Always advance one step in direction set by signX and then follow with the necessary amount of steps in the y-direction. More below */
                    y_stepSave = y_step;
                    y_step = std::floor(slope * (abs(x_step - screenX_prev))) + screenY_prev;

                    /* Handle slopes of pi/4 (45°) < gradient < inf */
                    if (abs(y_step - y_stepSave) > 1) {
                        for (int y_step_step = 0; y_step_step < abs(y_step - y_stepSave) - 1; y_step_step++) {
                            if (signY > 0) {
                                u_prev[IX(x_step, y_stepSave + y_step_step)] = ((double)(screenX - screenX_prev) * h) / dt;
                                v_prev[IX(x_step, y_stepSave + y_step_step)] = ((double)(screenY - screenY_prev) * h) / dt;
                            }
                            if (signY < 0) {
                                u_prev[IX(x_step, y_stepSave - y_step_step)] = ((double)(screenX - screenX_prev) * h) / dt;
                                v_prev[IX(x_step, y_stepSave - y_step_step)] = ((double)(screenY - screenY_prev) * h) / dt;
                            }
                        }
                    }
                }
            }
        }
    }
}

/* The main simulation loop. Call this to execute the 2D eulerian fluid simulation. Function returns when passed GLFWwindow* object is closed. */
void simulate(int N, int width, double gridSize, double dt, double* dens, double* dens_prev, double* u, double* u_prev, double* v, double* v_prev, double diff, double dyeIntensity, int brushSize, GLFWwindow* window) {
    /* Used for tracking user input */
    double cursorX_prev = 0, cursorY_prev = 0, cursorX = 0, cursorY = 0;

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window)) {
        cursorX_prev = cursorX;
        cursorY_prev = cursorY;
        glfwGetCursorPos(window, &cursorX, &cursorY);
        glClear(GL_COLOR_BUFFER_BIT);
        getSourceValuesFromUI(N, dt, width, dyeIntensity, brushSize, dens_prev, u_prev, v_prev, abs(cursorX_prev), abs(cursorY_prev), abs(cursorX), abs(cursorY), window);

            int i = 4;
            int j = 100;
            u_prev[IX(i, j)] = 1;
            dens_prev[IX(i, j)] = 5;
            v_prev[IX(i, j)] = 0;

             i = 196;
             j = 100;
            u_prev[IX(i, j)] = -0.015;
            dens_prev[IX(i, j)] = 2;
            v_prev[IX(i, j)] = 0;

        velocity_step(N, dt, diff, u, u_prev, v, v_prev);
        density_step(N, dt, diff, dens, dens_prev, u, v);
        drawGrid(dens, N, gridSize);
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}
/* hello future me
   if you ever decide to come back to this cause youre feeling especially bored or need an ego boost or whatever you left off in the getSourceValuesFromUI subroutine.
   the commented out part in that subroutine handles left mouse inputs appropriatly. the part that isnt commented out was there to implement the ability to swipe your
   cursor across the screen in any direction and have the fluid just sort of wrap around the window (however the boundary conditions are still no-slip and NOT periodic).
   that part is not finished yet.
   here are a few other things you might want to add since youre coming back to this now:
    - vorticity confinement (look it up incase you dont remember but i think you will cause im pretty sure youll be nostalgic for this time of your life by then) 
    - colors
    - non-square window support
    - make the 'brushSize' parameter in the sourceValuesFromUI routine actually do what its supposed to (increase the radius of where youre adding densities and velocities) 
    - obstacles: handle them and add the ability to place them wherever you please
    - give the user the ability to add constant source values wherever they please, and direction and velocity, could even come from the boundary cells
    - make this whole program run on the gpu, i.e. make it do all the calculations on the gpu
   ok byeee (also i hope you didnt kill yourself yet that would be a bummer) */
int main()
{
    GLFWwindow* window;
    int width = 800;
    int height = 800;
  //  double screenRatio = (double)width / (double)height;
    const int N = 200;
    const double dt = 0.5;
    double gridSize = 1.0/(double)N;
    const int totalAmountOfCells = (N+2) * (N+2);
    double *u = nullptr, *u_prev = nullptr, *v = nullptr, *v_prev = nullptr, *dens = nullptr, *dens_prev = nullptr, *p = nullptr, *div = nullptr;

    initializeArrays(u, u_prev, v, v_prev, dens, dens_prev, totalAmountOfCells);
    initializeWindow(window, width, height, "ouroborous 41-100%");
    simulate(N, width, gridSize, dt, dens, dens_prev, u, u_prev, v, v_prev, 0.000000002299, 6, 3, window);
    glfwTerminate();

    delete[] u;
    delete[] v;
    delete[] u_prev;
    delete[] v_prev;
    delete[] dens;
    delete[] dens_prev;
    return 0;
}