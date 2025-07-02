#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <windows.h>
#include <string>
#include <ctime>
#include <iomanip>

#define ITERATIONS 50
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
    for (int k = 0; k < ITERATIONS; k++) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                x[IX(i,j)] = (x0[IX(i,j)] + a * (x[IX(i+1,j)] + x[IX(i-1,j)] + x[IX(i,j+1)] + x[IX(i,j-1)])) / (1+4*a);
            }
        }
        /* When we're diffusing densities, this line is irrelevant, however when diffusing velocities it is crucial.
           This line makes sure that whenever we update a cell, that is next to a boundary cell's velocity,
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
    for (k = 0; k < ITERATIONS; k++) {
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

/* The main simulation loop. Call this to execute the 2D eulerian fluid simulation. Function returns when passed GLFWwindow* object is closed. */
void simulate(int N, int width, double gridSize, double dt, double* dens, double* dens_prev, double* u, double* u_prev, double* v, double* v_prev, double diff, double dyeIntensity, int brushSize, GLFWwindow* window) {
    /* Used for tracking user input */
    double cursorX_prev = 0, cursorY_prev = 0, cursorX = 0, cursorY = 0;

    int frames = 0;
    double avg[30] = {0};
    int it = 0;
    time_t start = time(0);
    time_t prev = time(0);
    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window) && time(0) - start < 30) {
        if (time(0) - prev >= 1) {
            printf("time(0)-start: %f\n", difftime(time(0), start));
            printf("time(0)-prev: %f\n", difftime(time(0), prev));
            printf("Cells/s: %d\n", frames * (N * N));
            avg[it++] = (frames * (N * N)) / (difftime(time(0), prev));
            printf("Avg[%d]: %f\n", it - 1, avg[it - 1]);
            frames = 0;
            prev = time(0);
        }
        if (time(0) - start >= 30) break;
        glClear(GL_COLOR_BUFFER_BIT);

        int i = 200;
        int j = 255;
        u_prev[IX(i, j)] = 0.5;
        dens_prev[IX(i, j)] = 5;
        v_prev[IX(i, j)] = 0;

        i = 300;
        j = 255;
        u_prev[IX(i, j)] = -0.5;
        dens_prev[IX(i, j)] = 5;
        v_prev[IX(i, j)] = 0;

        velocity_step(N, dt, diff, u, u_prev, v, v_prev);
        density_step(N, dt, diff, dens, dens_prev, u, v);
        drawGrid(dens, N, gridSize);

        for (int i = 0; i < (N + 2) * (N + 2); i++) {
            u_prev[i] = 0;
            v_prev[i] = 0;
            dens_prev[i] = 0;
        }
        
        glfwSwapBuffers(window);
        glfwPollEvents();
        frames++;
    }

    time_t end = time(0);
    double total_avg = 0;
    for (int i = 0; i < it; i++) {
        total_avg += avg[i];
    }
    printf("it: %d\n", it);
    printf("Avg Cells/s: %f\nTotal_Avg: %f\nSeconds: %f\n", total_avg / (double)it, total_avg, difftime(end, start));
    fflush(stdout);
}

int main()
{
    GLFWwindow* window;
    int width = 1440;
    int height = 1440;
  //  double screenRatio = (double)width / (double)height;
    const int N = 544;
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
