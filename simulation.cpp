#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "gillespie.h"
#include "Pearson.h"

const float Nmax = 200000;

///////// initialisation

long MCeq = 1000;
long MCsteps = 20000;


void generateDataFromUD(){
    
    float bT = 0,  bC = 0,  beta = 0,  kTpl = 0,  kCpl = 0,  kTmin = 0,  kCmin = 0,  kappaT = 0,  kappaC = 0,  dT = 0,  dC = 0, sigmaT = 0, sigmaC = 0, delta = 0,  kpl = 0,  kmin = 0,  deltaP = 0,  alphaT = 0,  deltaT = 0, alphaC = 0, deltaC = 0;
    
    // reading parameters
    float params [21] = { bT,  bC,  beta,  kTpl,  kCpl,  kTmin,  kCmin,  kappaT,  kappaC,  dT,  dC, sigmaT, sigmaC, delta,  kpl,  kmin,  deltaP,  alphaT,  deltaT, alphaC, deltaC};
    char s[20];
    
    FILE * paramsFile = fopen("parameters","r"); // read mode
    int i= 0;
    
    if (paramsFile != NULL){
        while (fgetc(paramsFile) != EOF)
        {
            fscanf(paramsFile, "%s\t%f\n", s, &params[i]);
            i ++;
        }
    }
    
    fclose(paramsFile);
    
    
    // finish reading parameters
    
    FILE * pFile = fopen ("alll1", "w");
    
    if( Nmax  - params[15] - params[16] - params[14] > mT + mC + cT + cC + mu)
    {
        // equilibrating
        Tstop = MCeq;
        gillespie(params, pFile);
        
        // data recording at the steady state
        Tstop = MCsteps;
        gillespie(params, pFile);
        
    }
    
    fclose(pFile);
    
}


/////////////////////////-------------- main

int main(){
    
    srand (time(NULL));

    
    float pearsonCyclesP = 0;
    float pearsonCyclesM = 0;
    
    float pearsonVarM = 0;
    float pearsonVarP = 0;
    
    float pTCycle = 0;
    float variancePCycle = 0;
    
    float CPCycle = 0;
    float varianceCPCycle = 0;

    
    // averaging results over Ncycles simulations
    float Ncycles = 0;
    
    while(Ncycles < 10){
        
        generateDataFromUD();
        
        pearsonVarP = Pearson(-1);
        pearsonVarM = Pearson(1);
        
        pearsonCyclesP += pearsonVarP;
        pearsonCyclesM += pearsonVarM;
        
        pTCycle += meanpT;
        variancePCycle += variancepT;
        
        CPCycle += meanCP;
        varianceCPCycle += varianceCP;
        
        Ncycles ++;
    }
    

    pearsonCyclesP = pearsonCyclesP/Ncycles;
    pearsonCyclesM = pearsonCyclesM/Ncycles;
    
    pTCycle = pTCycle/Ncycles;
    variancePCycle = variancePCycle/Ncycles;
    
    CPCycle = CPCycle/Ncycles;
    varianceCPCycle = varianceCPCycle/Ncycles;
    
    printf("C_P\t\tVar(C_P)\tP_T\t\tVar(P_T)\tP(P_T, P_C)\tP(m_T, m_C)\n ");
    printf("%f\t%f\t%f\t%f\t%f\t%f\n ", CPCycle, varianceCPCycle, pTCycle, variancePCycle, pearsonCyclesP, pearsonCyclesM);
    
    
    return 0;
    
}
