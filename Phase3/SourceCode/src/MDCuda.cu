#include "MDCuda.h"

// -------------  CUDA ---------------

#define NUM_THREADS_PER_BLOCK 512

// -------------  CUDA ---------------

// Number of particles
int N = 5000;

// -------------  CUDA ---------------

__device__ int d_N;
__device__ double d_PE;

// -------------  CUDA ---------------

double PE;
double NA = 6.022140857e23;
double kBSI = 1.38064852e-23;  // m^2*kg/(s^2*K)

//  Size of box, which will be specified in natural units
double L;

//  Initial Temperature in Natural Units
double Tinit;  //2;

//  Vectors!
const int MAXPART=5001;
//  Position
double* r = (double *) malloc(MAXPART*3*sizeof(double));
//  Velocity
double* v= (double *) malloc(MAXPART*3*sizeof(double));
//  Acceleration
double* a= (double *) malloc(MAXPART*3*sizeof(double));

// -------------  CUDA ---------------

double *d_r, *d_a; 
double aux = N * 3 * sizeof(double);

// -------------  CUDA ---------------

char *atype = (char *)malloc(3 * sizeof(char));

//  Function prototypes

//  initialize positions on simple cubic lattice, also calls function to initialize velocities
void initialize();  

//  update positions and velocities using Velocity Verlet algorithm 
//  print particle coordinates to file for rendering via VMD or other animation software
//  return 'instantaneous pressure'
double VelocityVerlet(double dt, FILE *fp);  

//  Compute Force using F = -dV/dr
//  solve F = ma for use in Velocity Verlet
//  Compute total potential energy from particle coordinates
__global__ void computeAccelerationsPotentialGPU(double *d_a, double *d_r, double *d_Pot);

//  Numerical Recipes function for generation gaussian distribution
double gaussdist();

//  Initialize velocities according to user-supplied initial Temperature (Tinit)
void initializeVelocities();

//  Compute mean squared velocity from particle velocities and total kinetic energy from particle mass and velocities
double MeanSquaredVelocityKinetic();

__device__ double atomicAddDouble(double* address, double val); 

void computeAccelerationsPotential();

int main(){

    int i, NumTime;
    double dt, Vol, Temp, Press, Pavg = 0, Tavg = 0, rho, VolFac, TempFac, PressFac, timefac, KE, mvs, gc, Z;
    char prefix[1000], tfn[1000], ofn[1000], afn[1000];
    FILE *tfp, *ofp, *afp;
    
    scanf("%s",prefix);
    strcpy(tfn,prefix);
    strcat(tfn,"_traj.xyz");
    strcpy(ofn,prefix);
    strcat(ofn,"_output.txt");
    strcpy(afn,prefix);
    strcat(afn,"_average.txt");
    
    scanf("%s",atype);
    
    if (strcmp(atype,"He")==0) {
        
        VolFac = 1.8399744000000005e-29;
        PressFac = 8152287.336171632;
        TempFac = 10.864459551225972;
        timefac = 1.7572698825166272e-12;
        
    }
    else if (strcmp(atype,"Ne")==0) {
        
        VolFac = 2.0570823999999997e-29;
        PressFac = 27223022.27659913;
        TempFac = 40.560648991243625;
        timefac = 2.1192341945685407e-12;
        
    }
    else if (strcmp(atype,"Ar")==0) {
        
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
        
    }
    else if (strcmp(atype,"Kr")==0) {
        
        VolFac = 4.5882712000000004e-29;
        PressFac = 59935428.40275003;
        TempFac = 199.1817584391428;
        timefac = 8.051563913585078e-13;
        
    }
    else if (strcmp(atype,"Xe")==0) {
        
        VolFac = 5.4872e-29;
        PressFac = 70527773.72794868;
        TempFac = 280.30305642163006;
        timefac = 9.018957925790732e-13;
        
    }
    else {
        
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
        strcpy(atype,"Ar");
        
    }

    scanf("%lf",&Tinit);

    if (Tinit<0.) {
        printf("\n  !!!!! ABSOLUTE TEMPERATURE MUST BE A POSITIVE NUMBER!  PLEASE TRY AGAIN WITH A POSITIVE TEMPERATURE!!!\n");
        exit(0);
    }

    if (N>=MAXPART) {
        printf("\n\n\n  MAXIMUM NUMBER OF PARTICLES IS %i\n\n  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY \n\n", MAXPART);
        exit(0);
    }

    // Convert initial temperature from kelvin to natural units
    Tinit /= TempFac;
    
    scanf("%lf",&rho);

    // -------------  CUDA ---------------

    // Copy N to the device variable d_N
    cudaMemcpyToSymbol(d_N, &N, sizeof(int));

    // -------------  CUDA ---------------

    Vol = N/(rho*NA);
    
    Vol /= VolFac;

    if (Vol<N) {
        printf("\n\n\n  YOUR DENSITY IS VERY HIGH!\n\n");
        printf("  THE NUMBER OF PARTICLES IS %i AND THE AVAILABLE VOLUME IS %f NATURAL UNITS\n",N,Vol);
        printf("  SIMULATIONS WITH DENSITY GREATER THAN 1 PARTCICLE/(1 Natural Unit of Volume) MAY DIVERGE\n");
        printf("  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY AND RETRY\n\n");
        exit(0);
    }

    // Length of the box in natural units:
    L = cbrt(Vol);
    
    //  Files that we can write different quantities to
    tfp = fopen(tfn,"w");    //  The MD trajectory, coordinates of every particle at each timestep
    ofp = fopen(ofn,"w");    //  Output of other quantities (T, P, gc, etc) at every timestep
    afp = fopen(afn,"w");    //  Average T, P, gc, etc from the simulation
    
    NumTime = 200;
    dt = 0.5e-14/timefac;

    if (strcmp(atype,"He")==0) {
        dt = 0.2e-14/timefac;
        NumTime=50000;
    }

    initialize();
    
    // -------------  CUDA ---------------

    computeAccelerationsPotential();

    // -------------  CUDA ---------------

    fprintf(tfp,"%i\n",N);
    fprintf(ofp,"  time (s)              T(t) (K)              P(t) (Pa)           Kinetic En. (n.u.)     Potential En. (n.u.) Total En. (n.u.)\n");

    for (i=0; i<NumTime+1; i++) {
        
        Press = VelocityVerlet(dt, tfp);
        Press *= PressFac;

        mvs = MeanSquaredVelocityKinetic()/N;
        KE = MeanSquaredVelocityKinetic()*0.5;

        Temp = mvs/3 * TempFac;

        gc = NA*Press*(Vol*VolFac)/(N*Temp);
        Z  = Press*(Vol*VolFac)/(N*kBSI*Temp);
        
        Tavg += Temp;
        Pavg += Press;
        
        fprintf(ofp,"  %12.12e  %12.12f  %12.12f %12.12f  %12.12f  %12.12f \n",i*dt*timefac,Temp,Press,KE, PE, KE+PE);
    }
    
    // Because we have calculated the instantaneous temperature and pressure,
    // we can take the average over the whole simulation here
    Pavg /= NumTime;
    Tavg /= NumTime;
    Z = Pavg*(Vol*VolFac)/(N*kBSI*Tavg);
    gc = NA*Pavg*(Vol*VolFac)/(N*Tavg);

    fprintf(afp,"  Total Time (s)      T (K)               P (Pa)      PV/nT (J/(mol K))         Z           V (m^3)              N\n");
    fprintf(afp," --------------   -----------        ---------------   --------------   ---------------   ------------   -----------\n");
    fprintf(afp,"  %12.12e  %12.12f       %12.12f     %12.12f       %12.12f        %12.12e         %i\n",i*dt*timefac,Tavg,Pavg,gc,Z,Vol*VolFac,N);
    
    printf("\n  AVERAGE TEMPERATURE (K):                 %15.5f\n",Tavg);
    printf("\n  AVERAGE PRESSURE  (Pa):                  %15.5f\n",Pavg);
    printf("\n  PV/nT (J * mol^-1 K^-1):                 %15.5f\n",gc);
    printf("\n  PERCENT ERROR of pV/nT AND GAS CONSTANT: %15.5f\n",100*fabs(gc-8.3144598)/8.3144598);
    printf("\n  THE COMPRESSIBILITY (unitless):          %15.5f \n",Z);
    printf("\n  TOTAL VOLUME (m^3):                      %10.5e \n",Vol*VolFac);
    printf("\n  NUMBER OF PARTICLES (unitless):          %i \n", N);
    
    free(r);
    free(v);
    free(a);
    free(atype);
    fclose(tfp);
    fclose(ofp);
    fclose(afp);
    
    return 0;
}

void initialize() {
    int n, p, i, j, k;
    double pos;
    
    // Number of atoms in each direction
    n = int(ceil(cbrt(N)));
    
    //  spacing between atoms along a given direction
    pos = L / n;
    
    //  index for number of particles assigned positions
    p = 0;
    double xPos, yPos, halfPos;
    halfPos = 0.5 * pos;

    //  initialize positions
    for (i=0; i<n; i++) {
        xPos = i*pos + halfPos;
        for (j=0; j<n; j++) {
            yPos = j*pos + halfPos;
            for (k=0; k<n; k++) {
                if (p<N*3) {
                    r[p++] = xPos;
                    r[p++] = yPos;
                    r[p++] = k*pos + halfPos;
                }
            }
        }
    }
    
    // Call function to initialize velocities
    initializeVelocities();
}

//  Function to calculate the averaged velocity squared
double MeanSquaredVelocityKinetic() {
    double v2 = 0;

    for (int i = 0; i < N * 3; i += 5) {
        v2 += v[i] * v[i] + v[i + 1] * v[i + 1] + v[i + 2] * v[i + 2] + v[i + 3] * v[i + 3] + v[i + 4] * v[i + 4];
    }

    return v2;
}

// -------------  CUDA ---------------

// __device__ double atomicAddDouble(double* address, double val) {
//     unsigned long long int* address_as_ull = (unsigned long long int*)address;
//     unsigned long long int old = *address_as_ull, assumed;

//     do {
//         assumed = old;
//         old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +
//                        __longlong_as_double(assumed)));
//     } while (assumed != old);

//     return __longlong_as_double(old);
// }

//   Uses the derivative of the Lennard-Jones potential to calculate
//   the forces on each atom.  Then uses a = F/m to calculate the
//   accelleration of each atom. 
__global__ void computeAccelerationsPotentialGPU(double *d_a, double *d_r, double *d_Pot) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    __shared__ double shared_d_r[NUM_THREADS_PER_BLOCK * 3];

    // Each thread loads the values of d_r into shared memory
    for (int k = 0; k < 3; ++k) {
        shared_d_r[threadIdx.x * 3 + k] = d_r[i * 3 + k];
    }

    if (i < d_N) {
        double Pot = 0.0;
        double d_a_aux[3] = {0.0, 0.0, 0.0};

        for (int j = 0; j < d_N; j++) {

            if (i != j) {

                double vals[3];
                double rij[3];
                double rSqd = 0;

                rij[0] = shared_d_r[threadIdx.x * 3] - d_r[j * 3];
                rij[1] = shared_d_r[threadIdx.x * 3 + 1] - d_r[j * 3 + 1];
                rij[2] = shared_d_r[threadIdx.x * 3 + 2] - d_r[j * 3 + 2];

                rSqd = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

                double rSqd3 = rSqd*rSqd*rSqd;
                double rSqd6 = rSqd3*rSqd3;
                Pot += ((1-rSqd3)/(rSqd6));
            
                double f = ((48 - 24*rSqd3)/(rSqd6*rSqd));

                vals[0] = rij[0] * f; 
                vals[1] = rij[1] * f; 
                vals[2] = rij[2] * f;

                d_a_aux[0] += vals[0];
                d_a_aux[1] += vals[1];
                d_a_aux[2] += vals[2];
            }
        }
        d_Pot[i] = Pot;

        d_a[i * 3] = d_a_aux[0];
        d_a[i * 3 + 1] = d_a_aux[1];
        d_a[i * 3 + 2] = d_a_aux[2];
    }
}

void computeAccelerationsPotential() {

    double Pot = 0.0;
    double h_Pot[N];
    double* d_Pot;

    int size = N * 3; 

    for (int i = 0; i < size; i++) { 
        a[i] = 0;
    }

    for (int i = 0; i < N; i++)
        h_Pot[i] = 0;

    cudaMalloc((void**)&d_r, aux);
    cudaMalloc((void**)&d_a, aux);
    cudaMalloc((void**)&d_Pot, N * sizeof(double));
    checkCUDAError("Mem Allocation");

    // Copy all working Data to Device
    cudaMemset(d_a, 0, aux);
    cudaMemset(d_Pot, 0,  N * sizeof(double));
    cudaMemcpy(d_r, r, aux, cudaMemcpyHostToDevice);
    checkCUDAError("Memcpy Host -> Device");

    int blocos = (N + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK;  
    computeAccelerationsPotentialGPU<<<blocos, NUM_THREADS_PER_BLOCK>>>(d_a, d_r, d_Pot);
    cudaDeviceSynchronize();
    checkCUDAError("Error in computeAccelerationsPotentialGPU");

    // copy the output to the host
    cudaMemcpy(a, d_a, aux, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_Pot, d_Pot, N * sizeof(double), cudaMemcpyDeviceToHost);
    checkCUDAError("Memcpy Device -> Host");

    for (int i = 0; i < N; i++) {
        Pot += h_Pot[i];
    }

    // free the device memory
    cudaFree(d_r);
    cudaFree(d_a);
    cudaFree(d_Pot);
    checkCUDAError("Free Mem");

    PE = Pot * 4;
}

// -------------  CUDA ---------------

double VelocityVerlet(double dt, FILE *fp) {
    
    int i;
    double psum = 0., temp1, temp2, dt1 = 0.5 * dt;

    for (i=0; i<N*3; i += 2) {
        temp1 = a[i] * dt1;
        r[i] += (v[i] + temp1) * dt;
        v[i] += temp1;

        temp2 = a[i+1] * dt1;
        r[i+1] += (v[i+1] + temp2) * dt;
        v[i+1] += temp2;
    }

    // -------------  CUDA ---------------

    computeAccelerationsPotential();

    // -------------  CUDA ---------------

    for (i=0; i<N*3; i += 2) {
        v[i] += a[i] * dt1;
        v[i+1] += a[i+1] * dt1;
    }
    
    // Elastic walls
    for (i=0; i<N*3; i += 2) {
        if (r[i]<0. || r[i]>=L) {
            v[i] *=-1.;
            psum += fabs(v[i]);
        }

        if (r[i+1]<0. || r[i+1]>=L) {
            v[i+1] *=-1.;
            psum += fabs(v[i+1]);
        }
    }
    
    return psum/(3*L*L*dt);
}

void initializeVelocities() {
    
    double vCM[3] = {0, 0, 0};
    double vSqdSum=0.;
    double lambda;
    
    for (int i=0; i < N*3; i += 2) {
        //  Pull a number from a Gaussian Distribution
        v[i] = gaussdist();
        v[i+1] = gaussdist();
    }
    
    for (int i=0; i<N; i++) {
        vCM[0] += v[i*3+0];
        vCM[1] += v[i*3+1];
        vCM[2] += v[i*3+2];
    }

    vCM[0] /= N;
    vCM[1] /= N;
    vCM[2] /= N;
    
    for (int i=0; i<N; i++) {
        v[i*3+0] -= vCM[0];
        v[i*3+1] -= vCM[1];
        v[i*3+2] -= vCM[2];
    }

    for (int i = 0; i < N * 3; i += 5) {
        vSqdSum += v[i] * v[i] + v[i + 1] * v[i + 1] + v[i + 2] * v[i + 2] + v[i + 3] * v[i + 3] + v[i + 4] * v[i + 4];
    }
    
    lambda = sqrt( 3*(N-1)*Tinit/vSqdSum);
    
    for (int i=0; i<N*3; i +=2) {
        v[i] *= lambda;
        v[i+1] *= lambda;
    }
}

//  Numerical recipes Gaussian distribution number generator
double gaussdist() {
    static bool available = false;
    static double gset;
    double fac, rsq, v1, v2, returnValue;
    if (!available) {
        do {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;
        
        returnValue =  v2*fac;
    } else {
        
        available = false;
        returnValue = gset;
        
    }
    return returnValue;
}
