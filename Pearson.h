#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <float.h>
#include <vector>
#include <iostream>     // std::cout


using std::vector;

#define _USE_MATH_DEFINES

//calculating Pearson correlation coefficient
float Pearson(int switching) // -1 for proteins and 1 for mRNAs
{
    
    float sigmap1 = 0;
    float sigmap2 = 0;
    
    float avp1 = 0;
    float avp2 = 0;
    
    long Nentries = 0;
    
    float pearsonCoeff = 0;
    
    long p1p2 = 0;
    
    long p1 = 0, p2=0, length, m1, m2, mu, c1, c2, c1m1, c2m2;
    
    float timeMu, var;
    
    FILE * pFile;
    
    char s[50];
    
    /////pay attention on filename
    pFile = fopen ("alll1", "r");
    
    if (pFile != NULL){
        
        //  printf("kuku%i\t", var<varmax );
        while (fgetc(pFile) != EOF)
        {
            
            if(switching < 0){ // between p1 p2
                fscanf(pFile, "%s\t%li\t%li\t%li\t%li\t%li\t%li\t%li\t%f\t%li\t%li\n", s, &length, &m1, &m1, &m2, &mu, &c1, &c2, &var, &c1m1, &c2m2);}
            else if (switching == 1) { // between m1 m2
                fscanf(pFile, "%s\t%li\t%li\t%li\t%li\t%li\t%li\t%li\t%f\t%li\t%li\n", s, &length, &c1, &c2, &m2, &mu, &m1, &m2, &var, &c1m1, &c2m2);
            }
            
            Nentries ++;
            avp1 += c1;
            avp2 += c2;
            
            sigmap1 += c1 * c1;
            sigmap2 += c2 * c2;
            
            p1p2 += c1 * c2;
            
            //   printf("p1p2: %li, %li, %li\n", p1p2, p1, p2);
            
        }
    }
    fclose (pFile);
    
    avp1 = (float) avp1/Nentries;
    avp2 = (float) avp2/Nentries;
    
    
    sigmap1 = (float) sqrt(sigmap1/Nentries - avp1*avp1);
    sigmap2 = (float) sqrt(sigmap2/Nentries - avp2*avp2);
    
    
    pearsonCoeff =  p1p2/Nentries - avp1* avp2;
    pearsonCoeff /= sigmap1 * sigmap2;
    
    return pearsonCoeff;
    
}



