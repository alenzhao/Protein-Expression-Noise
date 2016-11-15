//
//  gillespie.h
//
//
//  Created by Araks  Martirosyan on 14/11/16.
//
//

#ifndef ____gillespie__
#define ____gillespie__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "mt64.h"
#include "rando3.h"

#endif /* defined(____gillespie__) */

float Tstop = 30000;
float dt = 1;
long mT = 0, mC = 0, mu = 0, cT = 0, cC = 0, pC = 0, pT = 0, cP = 0;

float meanpT = 0;
float meanCP = 0;

float variancepT = 0;
float varianceCP = 0;

//Gillespie Algorithm
int gillespie ( float * params, FILE * pFile){
    
    float bT, bC, beta, kTpl, kCpl, kTmin, kCmin, kappaT, kappa2, dT, dC, sigmaT, sigmaC, delta, kpl, kmin, deltaP, alphaT, deltaT, alphaC, deltaC;
    
    // reading parameters
    bT =  params[0];
    bC = params[1];
    beta = params[2];
    kTpl = params[3];
    kCpl = params[4];
    kTmin = params[5];
    kCmin = params[6];
    kappaT = params[7];
    kappa2 = params[8];
    dT = params[9];
    dC = params[10];
    sigmaT = params[11];
    sigmaC = params[12];
    delta = params[13];
    kpl = params[14];
    kmin = params[15];
    deltaP = params[16];
    alphaT = params[17];
    deltaT = params[18];
    alphaC = params[19];
    deltaC = params[20];
    
   // for(int j = 1; j<21; j++){
   //     printf("%f\n", params[j]);
   // }
    
    init_genrand64(rand());
    
    unsigned long meanN = 0;
    
    float epsilon = DBL_EPSILON + 0.000000001;
    int stop = 0, i = 0;
    float  a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a19, a20, a21, a22, a23, a24, a25, a26, a27, omega;
    
    long mTt, mCt, mut, cTt, cCt;
    
    meanpT = 0;
    meanCP = 0;
    
    variancepT = 0;
    varianceCP = 0;
    
    
    float u = 0.0;
    float timeT = 0.01;
    float timeTflow = 0.0;
    float r1 = 0.0;
    
    // calculating the probabilities of events
    
    a1 = bT;
    
    a2 = dT * mT;
    a3 = kTpl * mT * mu;
    a4 = kTmin * cT;
    
    a5 = bC;
    a6 = dC * mC ;
    a7 = kCpl * mC * mu;
    a8 = kCmin * cC;
    
    a9 = beta;
    a10 = delta * mu;
    a11 = kappaT * cT;
    a12 = kappa2 * cC;
    
    a19 = sigmaT * cT;
    a20 = sigmaC * cC;
    
    a21 = mT * alphaT;
    a22 = mC * alphaC;
    a23 = deltaT * pT;
    a24 = deltaC * pC;
    
    a25 = kpl * pT * pC;
    a26 = deltaP * cP;
    a27 = kmin * cP;
    
    
    omega = a1 + a2 + a3 + a4 + a5 + a6 + a7+ a8 + a9 + a10 + a11 + a12 + a19 + a20 + a21 + a22 + a23 + a24 + a25 + a26 + a27;
    
    
    while ( timeTflow < Tstop && !( fabs(omega) < epsilon )) { // start while loop
        
        //random number in [1, 0]
        r1 =  genrand64_real1();
        
        //random number in [1, 0]
        u = genrand64_real1();
        
        // next reaction's time period
        timeT = (float) -log(u)/ omega;
        
        // choosing next reaction
        r1 = r1*omega;
        
        // ceRNAT
        if(0 <= r1  && r1  <= a1) {
            mT ++;
        }
        else if ( r1   <= a1 + a2){
            mT -- ;
        }
        //unbinding
        else if ( r1   <=  a1 + a2 + a3 ){
            mT -- ;
            mu -- ;
            cT ++ ;
        }
        else if (  r1   <= a1 + a2 + a3 + a4 ){
            
            mT ++ ;
            mu ++ ;
            cT -- ;
            
        }
        // ceRNAC
        else if (  r1   <= a1 + a2 + a3 + a4 + a5 ) {
            mC ++ ;
        }
        
        else if (  r1   <= a1 + a2 + a3 + a4 + a5 + a6){
            mC -- ;
            
        }
        // binding
        else if (   r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7){
            mC -- ;
            mu -- ;
            cC ++ ;
            
        }
        
        else if ( r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 ){
            mC ++;
            mu ++;
            cC --;
            
        }
        //miRNA
        else if( r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 ) {
            mu ++;
        }
        else if ( r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 ){
            mu -- ;
        }
        //recycling
        else if(  r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 ) {
            mu ++ ;
            cT --;
            
        }
        else if ( r1  <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 ){
            mu ++;
            cC --;
        }
        
        //c degradation
        else if( r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a19 ) {
            cT -- ;
        }
        
        else if ( r1   <=  a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a19 + a20 ){
            
            cC -- ;
            
        }
        // protein
        else if ( r1   <=  a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a19 + a20 + a21){
            
            pT ++;
            
        }
        else if ( r1   <=  a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a19 + a20 + a21 + a22){
            
            pC ++ ;
            
        }
        else if ( r1   <=  a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a19 + a20 + a21 + a22 + a23){
            
            pT -- ;
            
        }
        else if ( r1   <=  a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a19 + a20 + a21 + a22 + a23 + a24){
            
            pC -- ;
            
        }
        //protein complex
        else if ( r1   <=  a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a19 + a20 + a21 + a22 + a23 + a24 + a25){
            
            pC -- ;
            pT -- ;
            cP ++ ;
            
        }
        else if ( r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a19 + a20 + a21 + a22 + a23 + a24 + a25 + a26){
            cP -- ;
            
        }
        else if ( r1   <= a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a19 + a20 + a21 + a22 + a23 + a24 + a25 + a26 + a27){
            cP --;
            pT ++;
            pC ++;
            
        }
        
        timeTflow += timeT;
        
        // recording
        
        if(timeTflow > dt * i && mT>=0 && mC>=0   ) {
            
            fprintf( pFile, "tau:\t%d\t%li\t%li\t%li\t%li\t%li\t%li\t%f\t%li\t%li\n", i, mC, mT, cP, mu, pT, pC, 10.0, cT + mT, cC + mC);
            
            variancepT += pT*pT;
            meanpT += pT;
            
            meanCP += cP;
            varianceCP += cP * cP;
            
            meanN ++;
            i++;
        }
        
        a1 = bT ;
        a2 = dT * mT ;
        a3 = kTpl * mT * mu ;
        a4 = kTmin * cT ;
        
        a5 = bC;
        a6 = dC * mC ;
        a7 = kCpl * mC * mu ;
        a8 = kCmin * cC ;
        
        a9 = beta;
        a10 = delta * mu ;
        a11 = kappaT * cT ;
        a12 = kappa2 * cC ;
        
        a19 = sigmaT * cT ;
        a20 = sigmaC * cC ;
        
        a21 = mT * alphaT;
        a22 = mC * alphaC;
        a23 = deltaT * pT;
        a24 = deltaC * pC;
        
        a25 = kpl * pT * pC;
        a26 = deltaP * cP;
        a27 = kmin * cP;
        
        omega =  a1 + a2 + a3 + a4 + a5 + a6 + a7+ a8 + a9 + a10 + a11 + a12 + a19 + a20 + a21 + a22 + a23 + a24 + a25 + a26 + a27;
        // will define the timeT and next action probability
        
    } //end of  while loop
    
    
    // calculating average
    meanpT = meanpT/meanN;
    meanCP = meanCP/meanN;
    
    //calculating variance
    variancepT = variancepT/meanN - meanpT * meanpT;
    varianceCP = varianceCP/meanN - meanCP * meanCP;
    
    return 0;
    
}
