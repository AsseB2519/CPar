#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

// Number of particles
int N = 2160;

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
void computeAccelerationsPotential();

//  Numerical Recipes function for generation gaussian distribution
double gaussdist();

//  Initialize velocities according to user-supplied initial Temperature (Tinit)
void initializeVelocities();

//  Compute mean squared velocity from particle velocities and total kinetic energy from particle mass and velocities
double MeanSquaredVelocityKinetic();

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
    computeAccelerationsPotential();
    
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
    
    for (int i=0; i<N*3; i++) {
        v2 += v[i]*v[i];
    }
    
    return v2;
}

//   Uses the derivative of the Lennard-Jones potential to calculate
//   the forces on each atom.  Then uses a = F/m to calculate the
//   accelleration of each atom. 
void computeAccelerationsPotential() {
    int i, j, k, pos1, pos2;
    double Pot=0., f, rSqd, rSqd3, rSqd6, auxrij;
    double rij[3]; // position of i relative to j
    
    for (i = 0; i < N*3; i++) {  // set all accelerations to zero
        a[i] = 0;
    }

    for (i = 0; i < N-1; i++) {   // loop over all distinct pairs i,j
        pos1 = i*3;
        for (j = i+1; j < N; j++) {
            pos2 = j*3;
            rSqd = 0;
            rij[0]=r[pos1] - r[pos2];
            rij[1]=r[pos1+1] - r[pos2+1];
            rij[2]=r[pos1+2] - r[pos2+2];
            rSqd = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];

            rSqd3 = rSqd*rSqd*rSqd;
            rSqd6=rSqd3*rSqd3;
            Pot+=((1-rSqd3)/(rSqd6));
            
            //  From derivative of Lennard-Jones with sigma and epsilon set equal to 1 in natural units!
            f = ((48 - 24*rSqd3)/(rSqd6*rSqd));

            for (k = 0; k < 3; k++) {
                //  from F = ma, where m = 1 in natural units!
                auxrij= rij[k] * f;
                a[pos1+k] += auxrij;
                a[pos2+k] -= auxrij;
            }
        }
    }
    PE = Pot*8;
}

// returns sum of dv/dt*m/A (aka Pressure) from elastic collisions with walls
double VelocityVerlet(double dt, FILE *fp) {
    
    int i;
    double psum = 0., temp, dt1 = 0.5 * dt;

    //  Update positions and velocity with current velocity and acceleration
    for (i=0; i<N*3; i++) {
        temp = a[i] * dt1;
        r[i] += (v[i]+temp)*dt;
        v[i] += temp;
    }

    //  Update accellerations from updated positions
    computeAccelerationsPotential();

    //  Update velocity with updated acceleration
    for (i=0; i<N*3; i++) {
        v[i] += a[i]*dt1;
    }
    
    // Elastic walls
    for (i=0; i<N*3; i++) {
        if (r[i]<0. || r[i]>=L) {
            v[i] *=-1.;
            psum += fabs(v[i]);
        }
    }
    
    return psum/(3*L*L*dt);
}

void initializeVelocities() {
    
    int i, j;
    double vCM[3] = {0, 0, 0};
    double vSqdSum=0.;
    double lambda;
    
    for (i=0; i<N*3; i++) {
        //  Pull a number from a Gaussian Distribution
        v[i] = gaussdist();
    }
    
    // Vcm = sum_i^N  m*v_i/  sum_i^N  M
    // Compute center-of-mas velocity according to the formula above
    for (i=0; i<N; i++) {
        for (j=0; j<3; j++) {
            vCM[j] += v[i*3+j];
        }
    }
    
    for (i=0; i<3; i++) vCM[i] /= N;
    
    //  Subtract out the center-of-mass velocity from the
    //  velocity of each particle... effectively set the
    //  center of mass velocity to zero so that the system does
    //  not drift in space!
    for (i=0; i<N; i++) {
        for (j=0; j<3; j++) {
            v[i*3+j] -= vCM[j];
        }
    }
    
    //  Now we want to scale the average velocity of the system
    //  by a factor which is consistent with our initial temperature, Tinit
    for (i=0; i<N*3; i++) {
            vSqdSum += v[i]*v[i];
    }
    
    lambda = sqrt( 3*(N-1)*Tinit/vSqdSum);
    
    for (i=0; i<N*3; i++) {
        v[i] *= lambda;
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