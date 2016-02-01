/* nsmcore.c - Core NSM solver. Generates one trajectory of the process.  */

/* A. Hellander 2012-06-15 (Revision) */
/* P. Bauer and S. Engblom 2012-04-10 (Revision) */
/* A. Hellander 2009-11-24 (Revision) */
/* J. Cullhed 2008-06-18 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "propensities.h"
#include "nsm.h"
#include "binheap.h"
#include "report.h"
#include  "outputwriter.h"
#include "hdf5.h"
#include "hdf5_hl.h"





void nsm_core(const size_t *irD,const size_t *jcD,double *prD,
              const int *u0,
              const size_t *irN,const size_t *jcN,const int *prN,
              const size_t *irG,const size_t *jcG,
              const double *tspan,const size_t tlen,
              const double *vol,const double *data,const int *sd,
              const size_t Ncells,
              const size_t Mspecies,const size_t Mreactions,
              const size_t dsize,int report_level,
              const size_t *irK,const size_t *jcK,const double *prK,
              urdme_output_writer *writer,
              double theta,
              double Csize)
/* Specification of the inputs:
 
 Ncells
 Number of subvolumes.
 
 Mspecies
 Number of species.
 
 Hence Ndofs = Ncells*Mspecies.
 
 Mreactions
 Total number of reactions.
 
 dsize
 Size of data vector sent to propensities.
 
 tlen
 Number of sampling points in time.
 
 report_level
 The desired degree of feedback during simulations. 0, 1, and 2 are
 currently supported options.
 
 Diffusion matrix D. Double sparse (Ndofs X Ndofs).
 Macroscopic diffusion matrix. D(i,j) is the diffusion rate from dof #j to
 dof #i. This matrix uses the CSR-format and not CSC because fast access to
 rows is needed.
 
 Initial state vector u0. Integer (Mspecies X Ncells).
 Gives the initial copy number of the species in each subvolume.
 
 Stochiometric matrix N. Integer sparse (Mspecies X Nreactions).
 N(:,j) describes how reaction j changes the number of species.
 
 Dependency graph G. Integer sparse (Mreactions X Mspecies+Mreactions).
 G(i,Mspecies+j) is non-zero if executing reaction j means that reaction i
 needs to be re-evaluated. The first Mspecies columns of G similarily cover
 diffusion events.
 
 tspan. Double vector.
 Output times. tspan[0] is the start time and tspan[length(tspan)-1] is the
 stop time.
 
 vol. Double vector (length Ncells).
 vol[i] gives the volume of cell #i.
 
 data. Double matrix (dsize X Ncells).
 Generalized data matrix, data(:,j) gives a data vector for cell #j.
 
 sd. Integer vector (length Ncells).
 Subdomain number. sd[i] is the subdomain of cell #i. The vector sd can also
 be used to separate boundaries, line segments and points.
 
 Format of sparse matrices:
 G, N and S are sparse matrices in compressed column format (CCS). D is sparse
 but in compressed row format (CRS), or equivalently, a transposed matrix in
 CCS format.
 jcD, irD, prD (double *)
 jcN, irN, prN (int *)
 jcG, irG (int *)
 
 Propensities:
 a vector of function pointers (length Mreactions) is input by
 linking with the prototypes in propensities.h and function
 definitions in a user-specified .c-file. The type of this vector is
 PropensityFun which defines the input to a property function. See
 propensities.h for more details.
 
 Ordering of the dofs:
 Dof #i is located in cell #(i/Mspecies), and the dofs located in
 cell #j is u0(:,j). Thus, u0 is understood as a matrix of size
 Mspecies X Ncells.
 
 The output is a matrix U (Ndofs X length(tspan)).
 U(:,j) contains the state of the system at tspan(j).
 */
{
    double told,tt = tspan[0];
    double rdelta,rrdelta;
    double rand,cum,old;
    double *srrate,*rrate;
    double *sdrate,*Ddiag;
    double *rtimes;
    double old_rrate = 0.0,old_drate = 0.0;
    
    double totrate;
    
    int *node,*heap,*xx;
    long int total_reactions = 0;
    long int total_diffusion = 0;
    int dof,col;
    
    int subvol,event,re,spec,errcode = 0;
    size_t i,j,it = 0;
    size_t to_node,to_vol = 0;
    const size_t Ndofs = Ncells*Mspecies;


    printf("\ntheta = %f\n",theta);
    printf("Csize = %f\n",Csize);
    
/****************** Crowding Agent (Random) *************/
    double geo[Ncells];
    double rn;
    //double theta = 0.90; // Crowding Density
    //double Csize = 0.3; // 0.7: small crowding agent; 0.3: large crowding agent.
    int geo_temp = 0;
//    int sum = 0;
//    do {
//        sum = 0;
        for(int geo_index = 0;geo_index<Ncells;geo_index++){
           rn = drand48();
           geo_temp = 1-theta+rn;
	   geo[geo_index] = floor(geo_temp);
//	   if(geo_temp == 0) geo[geo_index]+= (1.00-theta);
	   if(geo_temp == 0) geo[geo_index]+= Csize;
//            if(geo[geo_index] == 0) sum+=1;
//        }
    }
//    while (sum != Ncells*(theta)); // Rejection Sampling
/********************************************************/
/********************* CHROMATIN ************************/
#define ARRAYSIZE(x)  (sizeof(x)/sizeof(*(x)))
   FILE *file;
   float array[Ncells][1];
   size_t iChrm, jChrm, kChrm;
   char buffer[BUFSIZ], *ptr;
   float Chrm[Ncells], MAX_Pixel = 190.00;
   file = fopen("/usr/local/pyurdme/pyurdme/urdme/src/nsm/avg.txt","r");
   if(file != NULL){
       for ( iChrm = 0; fgets(buffer, sizeof buffer, file); ++iChrm )
           for ( jChrm = 0, ptr = buffer; jChrm < ARRAYSIZE(*array); ++jChrm, ++ptr )
               array[iChrm][jChrm] = (float)strtof(ptr, &ptr);
       fclose(file);  
       for (size_t idxx = 0; idxx<Ncells; idxx++)
    //	 Chrm[idxx] = MAX_Pixel - array[idxx][0];
          { Chrm[idxx] = (MAX_Pixel-array[idxx][0])/MAX_Pixel;
	}
	//printf("\n%f  %f  %f",idxx,array[idxx][0],Chrm[idxx]);}
    }else{
       for (size_t idxx = 0; idxx<Ncells; idxx++)
           Chrm[idxx] = 1;
    }

for (int idxxx = 0; idxxx<Ncells; idxxx++)
      // Chrm[idxx] = MAX_Pixel - array[idxx][0];
      printf("\n%d %f %f",idxxx,array[idxxx][0],Chrm[idxxx]);
    
printf("\n\nNCells = %u\n\n",Ncells);
      // Chrm[idxx] = 66/(5+array[idxx][0]);
//	 Chrm[idxx] = 1;
/********************************************************/
    PropensityFun *rfun;
    rfun = ALLOC_propensities();
    
    ReportFun report;
    if (report_level)
        report = &reportFun1;
    else
        report = NULL;
//    
    /* Set xx to the initial state. xx will always hold the current solution. */
    xx = (int *)malloc(Ndofs*sizeof(int));
    memcpy(xx,u0,Ndofs*sizeof(int));
    
    /* Create reaction rate matrix (Mreactions X Ncells) and total rate
     vector. In rrate we store all propensities for chemical rections,
     and in srrate the sum of propensities in every subvolume. */
    rrate = (double *)malloc(Mreactions*Ncells*sizeof(double));
    srrate = (double *)malloc(Ncells*sizeof(double));
    
    
    /* Calculate the propensity for every reaction and every
     subvolume. Store the sum of the reaction intensities in each
     subvolume in srrate. */
    for (i = 0; i < Ncells; i++) {
        srrate[i] = 0.0;
        for (j = 0; j < Mreactions; j++) {
            //rrate[i*Mreactions+j] =
            //(*rfun[j])(&xx[i*Mspecies],tt,vol[i],&data[i*dsize],sd[i],i,xx,irK,jcK,prK);
            //srrate[i] += rrate[i*Mreactions+j];
            rrate[i*Mreactions+j] =
            (*rfun[j])(&xx[i*Mspecies],tt,vol[i],&data[i*dsize],sd[i]);
            srrate[i] += rrate[i*Mreactions+j];
        }
    }
    
    /* Total diffusion rate vector (length Mcells). It will hold
     the total diffusion rates in each subvolume. */
    sdrate = (double *)malloc(Ncells*sizeof(double));
    
    /* The diagonal value of the D-matrix is used frequently. For
     efficiency, we store the negative of D's diagonal in Ddiag. */
    
    Ddiag = (double *)malloc(Ndofs*sizeof(double));
    for (i = 0; i < Ndofs; i++) {
        Ddiag[i] = 0.0;
        size_t tmp_ndx;
        for (j = jcD[i]; j < jcD[i+1]; j++)
            if (irD[j] != i){
                prD[j] = prD[j]*geo[i/Mspecies]*Chrm[i/Mspecies];
                Ddiag[i] += prD[j];
            }else{
                tmp_ndx = j;
            }
        prD[tmp_ndx] = -1*Ddiag[i];
            //if (irD[j] == i) Ddiag[i] = -prD[j]*geo[i/Mspecies]*Chrm[i/Mspecies];
            //for (i = jcD[col], cum = 0.0; i < jcD[col+1]; i++)
            //if (irD[i] != col && (cum += prD[i]*geo[irD[i]/Mspecies]*Chrm[irD[i]/Mspecies]) > rand)
        
    }
    
    /* Calculate the total diffusion rate for each subvolume. */
    for(i = 0; i < Ncells; i++) {
        sdrate[i] = 0.0;
        for(j = 0; j < Mspecies; j++)
        sdrate[i] += Ddiag[i*Mspecies+j]*xx[i*Mspecies+j];
    }
    
    /* Create binary (min)heap. */
    rtimes = (double *)malloc(Ncells*sizeof(double));
    node = (int *)malloc(Ncells*sizeof(int));
    heap = (int *)malloc(Ncells*sizeof(int));
    
    /* Calculate times to next event (reaction or diffusion)
     in each subvolume and initialize heap. */
    for (i = 0; i < Ncells; i++) {
        rtimes[i] = -log(1.0-drand48())/(srrate[i]+sdrate[i])+tspan[0];
        heap[i] = node[i] = i;
    }
    initialize_heap(rtimes,node,heap,Ncells);
    
    
    /* Main loop. */
    for (;;) {
        
        /* Get the subvolume in which the next event occurred.
         This subvolume is on top of the heap. */
        told = tt;
        tt   = rtimes[0];
        subvol = node[0];
        
        /* Store solution if the global time counter tt has passed the
         next time is tspan. */
        
        if (tt >= tspan[it] || isinf(tt)) {
            
            for (; it < tlen && (tt >= tspan[it] || isinf(tt)); it++) {
                
                if (report){
                    report(tspan[it],tspan[0],tspan[tlen-1],total_diffusion,total_reactions,0,report_level);
                }
                
                write_state(writer,xx);
                
            }
            
            /* If the simulation has reached the final time, flush the buffer and exit. */
            if (it >= tlen){
                flush_buffer(writer);
                break;
            }
            
        }
        
        /* First check if it is a reaction or a diffusion event. */
        totrate = srrate[subvol]+sdrate[subvol];
        rand = drand48();
        
        
        if (rand*totrate <= srrate[subvol]) {
            /* Reaction event. */
            event = 0;
            
            /* a) Determine the reaction re that did occur (direct SSA). */
            rand *= totrate;
            for (re = 0, cum = rrate[subvol*Mreactions]; re < Mreactions && rand > cum; re++, cum += rrate[subvol*Mreactions+re])
            ;
            
            if (re >= Mreactions){
                re--;
            }
            
            /* b) Update the state of the subvolume subvol and sdrate[subvol]. */
            for (i = jcN[re]; i < jcN[re+1]; i++) {
                xx[subvol*Mspecies+irN[i]] += prN[i];
                if (xx[subvol*Mspecies+irN[i]] < 0) errcode = 1;
                sdrate[subvol] += Ddiag[subvol*Mspecies+irN[i]]*prN[i];
            }
            
            /* c) Recalculate srrate[subvol] using dependency graph. */
            for (i = jcG[Mspecies+re], rdelta = 0.0; i < jcG[Mspecies+re+1]; i++) {
                old = rrate[subvol*Mreactions+irG[i]];
                j = irG[i];
                //rdelta +=
                //(rrate[subvol*Mreactions+j] =
                // (*rfun[j])(&xx[subvol*Mspecies],tt,vol[subvol],&data[subvol*dsize],sd[subvol],subvol,xx,irK,jcK,prK)
                // )-old;
                rdelta +=
                (rrate[subvol*Mreactions+j] =
                 (*rfun[j])(&xx[subvol*Mspecies],tt,vol[subvol],&data[subvol*dsize],sd[subvol])
                 )-old;
            }
            srrate[subvol] += rdelta;
            
            total_reactions++; /* counter */
        }
        else {
            /* Diffusion event. */
            event = 1;
            
            /* a) Determine which species... */
            rand *= totrate;
            rand -= srrate[subvol];
            
            for (spec = 0, dof = subvol*Mspecies, cum = Ddiag[dof]*xx[dof];
                 spec < Mspecies && rand > cum;
                 spec++, cum += Ddiag[dof+spec]*xx[dof+spec]);
            
            
            /* b) and then the direction of diffusion. */
            col = dof+spec;
            rand = drand48()*Ddiag[col];
            
            /* Search for diffusion direction. */
            for (i = jcD[col], cum = 0.0; i < jcD[col+1]; i++)
            //if (irD[i] != col && (cum += prD[i]*geo[irD[i]/Mspecies]*Chrm[irD[i]/Mspecies]) > rand)
            if (irD[i] != col && (cum += prD[i]) > rand)
            break;
            
            /* paranoia fix: */
            if (i >= jcD[col+1]){
               i--;
            }
            
            to_node = irD[i];
            to_vol = to_node/Mspecies;
            
            /* c) Execute the diffusion event (check for negative elements). */
            //if (geo[to_vol]>0){
	    xx[subvol*Mspecies+spec]--;
            if (xx[subvol*Mspecies+spec] < 0){
                    errcode = 1;
            }
            xx[to_node]++;
           // } 
            /* Save reaction and diffusion rates. */
            old_rrate = srrate[to_vol];
            old_drate = sdrate[to_vol];
            
            /* Recalculate the reaction rates using dependency graph G. */
            if (Mreactions > 0){
                for (i = jcG[spec], rdelta = 0.0, rrdelta = 0.0; i < jcG[spec+1]; i++) {
                    
                    j = irG[i];
                    old = rrate[subvol*Mreactions+j];
                    
                    //rdelta +=
                    //  (rrate[subvol*Mreactions+j] =
                    //    (*rfun[j])(&xx[subvol*Mspecies],tt,vol[subvol],&data[subvol*dsize],sd[subvol],subvol,xx,irK,jcK,prK)
                    //  )-old;
                    rdelta +=
                      (rrate[subvol*Mreactions+j] =
                        (*rfun[j])(&xx[subvol*Mspecies],tt,vol[subvol],&data[subvol*dsize],sd[subvol])
                      )-old;
                    old = rrate[to_vol*Mreactions+j];
                    
                    //rrdelta += 
					//  (rrate[to_vol*Mreactions+j] = 
                    //    (*rfun[j])(&xx[to_vol*Mspecies],tt,vol[to_vol],&data[to_vol*dsize],sd[to_vol],to_vol,xx,irK,jcK,prK)
                    //  )-old;
                    rrdelta += 
					  (rrate[to_vol*Mreactions+j] = 
                        (*rfun[j])(&xx[to_vol*Mspecies],tt,vol[to_vol],&data[to_vol*dsize],sd[to_vol])
                      )-old;
                }
                
                srrate[subvol] += rdelta;
                srrate[to_vol] += rrdelta;
            }
            
            /* Adjust diffusion rates. */
           //if (geo[to_vol]>0){
	    sdrate[subvol] -= Ddiag[subvol*Mspecies+spec];
            sdrate[to_vol] += Ddiag[to_vol*Mspecies+spec];
            
            total_diffusion++; /* counter */
           // }
        }
        
        /* Compute time to new event for this subvolume. */
        totrate = srrate[subvol]+sdrate[subvol];  
        if (totrate > 0.0)
        rtimes[0] = -log(1.0-drand48())/totrate+tt;
        else
        rtimes[0] = INFINITY;
        
        /* Update the heap. */
        update(0,rtimes,node,heap,Ncells);
        
        /* If it was a diffusion event, also update the other affected
         node. */
        if (event) {
            totrate = srrate[to_vol]+sdrate[to_vol];      
            if (totrate > 0.0) {
                if (!isinf(rtimes[heap[to_vol]]))
                rtimes[heap[to_vol]] = 
                (old_rrate+old_drate)/totrate*(rtimes[heap[to_vol]]-tt)+tt;
                else
                /* generate a new waiting time */
                rtimes[heap[to_vol]] = -log(1.0-drand48())/totrate+tt;
            } 
            else
            rtimes[heap[to_vol]] = INFINITY;
            
            update(heap[to_vol],rtimes,node,heap,Ncells);
        } 
        
        /* Check for error codes. */
        if (errcode) {
            /* Report the error that occurred. */
            if (report)
                report(tt,tspan[0],tspan[tlen-1],total_diffusion,total_reactions,errcode,report_level);
            /* Cannot continue. Clear this solution and exit. */
            break;
        }
    }
    
    
    FREE_propensities(rfun);
    free(heap);
    free(node);
    free(rtimes);
    free(Ddiag);
    free(sdrate);
    free(srrate);
    free(rrate);
    free(xx); 
    
}
