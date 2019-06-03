
/* $Header: /disks/hermosa/public/NEURON/CA1-99/bp/lib/TEST/newshiftsyn.c 2000/02/25 poirazi */
/* Adding a temporal offset for stimulation of trains */

/* $Header: /disks/newport/freeway/brannon/rs/p-and-s/SourceCode/C/newshiftsyn.c,v 1.7 1999/03/24 18:19:21 brannon Exp $ */

/* $Log: newshiftsyn.c,v $
 * Revision 1.7  1999/03/24 18:19:21  brannon
 * Ok, now NEURON should be happy.
 *
 * Revision 1.6  1999/03/24 18:09:25  brannon
 * The end of newshiftsyn files is not compatible with NEURON's vector
 * scanf() function. Working on this.
 *
 * Revision 1.5  1998/10/22 10:05:03  brannon
 * The compiler croaked on nested comments so I fixed my attempt to
 * simply comment out and copy the old code by removing the old code and
 * associated comments.
 *
 * Revision 1.4  1998/10/22 10:02:58  brannon
 * I changed the sprintf which created filenames to write the p_seed in
 * integer as opposed to float from.
 *
 * Revision 1.3  1998/05/07 19:33:44  brannon
 * rm shiftsyn.log -> rm -f shiftsyn.log
 * My runs were held up because it was waiting for a reply to removal of
 * a previous shiftsyn.log
 *
 * Revision 1.2  1998/04/08 18:33:29  brannon
 * Write a \n after each dt of synaptic input.
 *
 * Revision 1.1  1998/04/08 18:32:48  brannon
 * Initial revision
 *
 * Revision 1.12  1997/04/15 19:35:59  niebur
 * *** empty log message ***
 *
 * Revision 1.11  1995/02/01  01:55:28  ernst
 * Eliminated finite loop condition: if s=1, the restart in the main loop after
 * enforce_refrac failed (cond=1) was ineffective, because make_train_offset
 * and make-spike trains each time gave identical results. This happened because
 * the sigmas in them contained terms 1-s and thus became zero.
 *
 * This has been corrected by restarting the spiketrain now COMPLETELY
 * (ie, all synapses) when that occurs.
 *
 * PS this condition occurred for " shiftsyn testfile 100 500 0.1 100 1 0.2 6074"
 *
 * Revision 1.10  1995/01/28  22:00:46  ernst
 * Writing now a log file at the beginning of the run which is
 * deleted at its successful completion.
 *
 * Revision 1.9  1995/01/28  08:04:27  ernst
 * Now enforcing refractory period (in new routine enforce_refrac).
 * Also found out that the Numerical recipes generator has a problem
 * if its seed is too large (> 54000, apparently). The reason is
 * unknown, but I do not allow such large seeds anymore.
 *
 * Revision 1.8  1995/01/20  22:13:26  ernst
 * Added more debugging message in context with the consistency check in write_to_disk.
 *
 * Revision 1.7  1995/01/18  02:41:29  ernst
 * Added a global variable ("errfile") which is identical to outfilename and
 * which is printed out each time there is an error message. This way, it will
 * be possible to determine which run makes a problem in a series of many runs
 * in which the shell might garble the order of messages.
 *
 * Revision 1.6  1995/01/18  02:01:03  ernst
 * removed debugging messages
 *
 * Revision 1.5  1994/10/29  02:25:07  ernst
 * Fixed a problem in  make_spike_trains: the sigma of the Gaussian there
 * was dependend on the total time. There is no reason this should be
 * the case.
 *
 * Revision 1.4  1993/10/04  08:11:56  ernst
 * Introduced computation and writing to a file of
 * the total activity at a given time (variable totact).
 *
 * Revision 1.3  1993/09/23  01:17:08  ernst
 * Introduced larger scatter for shift between spike trains.
 * The scatter is now in units of n_tot, NOT in units of the
 * average spiking frequency.
 * */

/* $Id: newshiftsyn.c,v 1.7 1999/03/24 18:19:21 brannon Exp $ */

/***************************************************************************/
/*                                                                         */
/* Purpose of this program: make time-dependent input for use with NEURON  */
/*                                                                         */
/* Input: explained before 'main'                                          */
/*                                                                         */
/* Output: file:                                                           */
/*                                                                         */
/*         syndata<parameters> contains spike times for the given number   */
/*         of synapses (read by NEURON),                                   */
/*                                                                         */
/* Note: the random number generator is defined in /cit/ernst/num-rec.h    */
/*                                                                         */
/****************************************************************************/

#include "/usr/include/stdlib.h"
#include "/usr/include/string.h"
//#include "/usr/include/math.h"
#include "ken.h"
#include "num-rec.h"



#define RAND_MAX 2147483647
#define myrand ( (double) random() / RAND_MAX)

char errfile[1000]; /* global var for debugging purposes*/
/********************************************************************/
/*                                                                  */
/* Allocate array space for times (doubles)                         */
/*                                                                  */
/********************************************************************/
double *alloctimes(int nelements)
{
int i;
double *p;
p = calloc(nelements,sizeof(double));

if(!p)
  {
    printf("allocation failure in allocspikes, working on %s\n; aborting",
	   errfile);
    exit(1);
  }
for(i=0;i<nelements;i++)
  p[i]=0;

return(p);
}

/********************************************************************/
/*                                                                  */
/* Allocate array space for spikes (ints)                           */
/*                                                                  */
/********************************************************************/
int *allocspikes(int nelements)
{
int i;
int *p;
p = calloc(nelements,sizeof(int));

if(!p)
  {
    printf("allocation failure in allocspikes, working on %s\n; aborting",
	   errfile);
    exit(1);
  }
for(i=0;i<nelements;i++)
  p[i]=0;

return(p);
}

/********************************************************************/
/*                                                                  */
/* Read input from console                                          */
/*                                                                  */
/********************************************************************/
void read_input(int argc,char *argv[],int *p_n_syn,double *p_t_tot,double *p_dt,double *p_av_rate,double *p_p,double *p_s,int *p_seed, double *p_offset, char *outfilename,char *totactname)
{
float f_t_tot,f_dt,f_av_rate,f_p,f_s,f_offset;

  if(argc != 10)
    {
      printf
	("Usage: shiftsyn outfile n_syn t_tot dt av_rate s p seed offset\n (Error file is %s)\n",errfile);
      exit(1);
    }
  else
    {

    

          sscanf(argv[2],"%i",p_n_syn);
          sscanf(argv[3],"%g",&f_t_tot);
          sscanf(argv[4],"%g",&f_dt);
          sscanf(argv[5],"%g",&f_av_rate);
          sscanf(argv[6],"%g",&f_s);
          sscanf(argv[7],"%g",&f_p);
          sscanf(argv[8],"%i",p_seed);
          *p_offset = atof(argv[9]); 
     
    }
/* (translating the floats into doubles)*/

if(f_p<0 || f_p>1 || f_p<0 || f_s>1 )
   {
     printf("p and s must be from [0..1];, working on %s EXIT!\n\n",errfile);
     exit(1);
   }

*p_t_tot = f_t_tot;  
*p_dt = f_dt ;		
*p_p=f_p;
*p_s=f_s;
  
/*

Brannon: I dont have the faintest idea why p_seed would be cast as float here. 
But I do have a non-faint idea as to why I want to make it integer.


*/

/* making the complete outfile name (NB: dt is in ms):*/
sprintf(outfilename,"%s-%d-%.2f-%.2f-%.2f-%.2f-%.2f-%d-%.2f",\
	argv[1],*p_n_syn,f_t_tot,f_dt,f_av_rate,f_s,f_p,*p_seed,*p_offset);

/* also write the complete outfile name to global variable which is
   used in error messages as a signature:*/
sprintf(errfile,"%s-%d-%.2f-%.2f-%.2f-%.2f-%.2f-%d-%.2f",\
	argv[1],*p_n_syn,f_t_tot,f_dt,f_av_rate,f_s,f_p,*p_seed,*p_offset);

/* making the name for the total activity file:*/
sprintf(totactname,"%s-%d-%.2f-%.2f-%.2f-%.2f-%.2f-%d-%.2f",\
	"totact",*p_n_syn,f_t_tot,f_dt,f_av_rate,f_s,f_p,*p_seed,*p_offset);


/* t_tot and dt are in ms, therefore also the rate has to be in spikes */
/* per ms and not spikes per second: */  
*p_av_rate = f_av_rate/1000.;	

}

/********************************************************************/
/*                                                                  */
/* make_train_offset:                                               */
/* Compute the offset for the whole spiketrain.                     */
/* A Gaussian random number with sigma s/(1-s)                      */
/*                                                                  */
/* s=0  <==> all spiketrains are completely periodic,               */
/* s=1  <==> all spiketrains are completely aperiodic,              */
/* s=0.5 <==> spiketrains are distributed with sigma= */
/*                                                                  */
/********************************************************************/

double make_train_offset(double p, double s, double av_rate, \
			 double t_tot,int *p_seed)
{
  double offset,delta_t;
  double sigma_p;
  
  delta_t = 1./av_rate;

  if(s > APRECISION && s<=1.)		/* finite value for s */
    {
      sigma_p = p*(1-s)/s * delta_t; /* determine sigma of gaussian */
      offset= sigma_p*gasdev(p_seed);
    }
  else if(fabs(s) < APRECISION)	     /* s is zero ==> "infinite sigma" */
    {
      offset= p* t_tot * ran1(p_seed);    /* pick random offset  */
    }
  else				/* error! */
    nrerror("error in make_train_offset: s must be 0 <= s <= 1");

return(offset);
}
/********************************************************************/
/*                                                                  */
/* dsort:                                                           */
/* sorting routine, from numrec.                                    */
/* Identical to routine sort there, exc. for one difference (shown) */
/*                                                                  */
/********************************************************************/

void dsort(n,ra)
int n;
double ra[];	
       /* declaration of ra and rra is the only difference to numrec's sort */
{
        int l,j,ir,i;
        double rra;

        l=(n >> 1)+1;
        ir=n;
        for (;;) {
                if (l > 1)
                        rra=ra[--l];
                else {
                        rra=ra[ir];
                        ra[ir]=ra[1];
                        if (--ir == 1) {
                                ra[1]=rra;
                                return;
                        }
                }
                i=l;
                j=l << 1;
                while (j <= ir) {
                        if (j < ir && ra[j] < ra[j+1]) ++j;
                        if (rra < ra[j]) {
                                ra[i]=ra[j];
                                j += (i=j);
                        }
                        else j=ir+1;
                }
                ra[i]=rra;
        }
}

/********************************************************************/
/*                                                                  */
/* sort_train:                                                      */
/* bring spikes in range [0..t_tot] and                             */
/* sort spiketrain in ascending order.                              */
/*                                                                  */
/********************************************************************/

sort_train(double *spiketrain,int length_train,double t_tot)
{
  int i;

  for (i=0;i<length_train; i++)
    {
      /*bring in range [0 .. t_tot]:*/

      spiketrain[i]=fmod(spiketrain[i],t_tot);

      while(spiketrain[i]<0)
	  spiketrain[i] += t_tot;

      if(spiketrain[i]>t_tot || spiketrain[i] <0)
	{
	  printf("sort_train: This cannot happen\n, working on %s",errfile);
	  exit(1);
	}
    }

/* DEBUG:*/
/*
  printf("\nNext train, in range 0-%g, unordered (working on %s):\n",
  t_tot,errfile);
  for (i=0;i<length_train; i++)
    {
      printf("%i %g\n",i,(float)spiketrain[i]);
    }
*/

/* use numerical recipes heapsort routine. ATTENTION: subtract 1 from
the pointer to the array spiketrain because the routine expects array
indices 1 to n, not 0 to n-1. See NRC p. 15 (unit offsets vs. zero
offsets) */

  dsort(length_train,spiketrain-1); 
}



/********************************************************************
 *                                                                  *
 * enforce_refrac:                                                  *
 * enforce refractory period.                                       *
 *                                                                  *
 * Method: If a spike is closer than t_ref to the next one,         *
 * the next one is pushed by t_refrac. This MAY lead the last       *
 * spike(s) to be pushed beyond t_tot. If that happens, we push     *
 * back the last spikes by t_ref/2 (not t_ref; this guarantees      *
 * that the order remains conserved).                               * 
 *                                                                  *
 * returns 1 if successful, -1 if not                               *
 *                                                                  *
 ********************************************************************/

int enforce_refrac(double *spiketrain, int syn, int length_train, 
		   double t_ref, double dt,double t_tot)
{
  int i,j,done;
  if(t_ref <= dt) {
    printf("ERROR1 in enforce_refrac, working on %s\n; aborting\n",
	   errfile);
    exit(1);
  }
  
  done=0;
  while( !done) {   /* loop goes to l_t-2; see below for very last interval*/
    for (i=0;i<length_train-2; i++)  
      if(spiketrain[i+1]-spiketrain[i] <= t_ref) {   /* if too close */
	spiketrain[i+1] += t_ref;                    /* push next spike*/
	j=i+1;             /*if spike order violated, keep pushing till ok:*/
	if(j+1 >= length_train) 
	  return(-1);          /* failure; make new spike train*/
	while((spiketrain[j]) > spiketrain[j+1]) {
	  spiketrain[j+1] += t_ref;
	  j++;
	  if(j+1 >= length_train) /*make sure j-loop stays w/in array bounds*/
	    return(-1);          /*else: failure; make new spike train*/
	}
      }
	
    /*Now spikes are separated by at least t_ref. Check if sth. fell off end:*/
    done=1;
    for (i=1;i<length_train; i++) 
      if(spiketrain[i] > t_tot) {   /*if yes, put the last one back a bit:*/
	spiketrain[i-1] -= t_ref/2.; /*subtract only t_r/2 to keep ordered*/
	if(spiketrain[i-1]<0)      /*failure; make new spike train*/
	  return(-1);
	done=0; /*have to do the while loop again (possibly more than once)*/
      }         
    }

/* The only interval that has not been checked is the very last one.
   In the unlikely case this is a problem, we would have to do a
   major re-ordering (spikes are already ordered!).
   Therefore, if this happens, declare failure, return with (-1) and make 
   an entirely new spike train:*/
  
  if(spiketrain[length_train-1]-spiketrain[length_train-2] <= t_ref)
    return(-1);
  else {
    /*Debug:*/
    /*
       printf("\nOrdered train %i, in range 0-%g (working on %s):\n",
       syn,t_tot,errfile);
       for (i=0;i<length_train; i++)
       {
       printf("%i %g\n",i,(float)spiketrain[i]);
       }
     */
    return(1);
  }
}

/********************************************************************/
/*                                                                  */
/* append_train:                                                    */
/* append spiketrain to alltrains.                                  */
/* In this array, spiketrains are arranged one after the other.     */
/* ie., first the complete train for synapse 1, then for syn. 2 etc */
/*                                                                  */
/********************************************************************/

append_train(double *sortedtrain,double *alltrains, int syn,\
	     int length_train)
{
int i;

for (i=0;i<length_train; i++)
  {
    /*store all synapses together in array:*/
    alltrains[syn*length_train+i]=sortedtrain[i];
  }

}

/********************************************************************/
/*                                                                  */
/* write_to_disk:                                                   */
/* write spike trains to disk and                                   */
/* write the total activity to a separate file.                     */
/* also do consistency checking.                                    */
/*                                                                  */
/********************************************************************/

void write_to_disk(double *alltrains,FILE *outfile,FILE *f_totact,\
		   int length_train,double t_tot,int n_syn,double dt,double offset)
{
  double t;
  int nspike,syn,outsyn;
  int *indnexttime; 
  double *nexttime;
  int totact;
  FILE *problem_file;
  char consist_filename[1000];

  /*To avoid that we have to read the array for each synapse
    completely at each time step, we store the next time at which a
    spike occurs in nexttime and its index in the array in
    indnexttime:*/
  
  indnexttime = allocspikes(n_syn); 
  nexttime = alloctimes(n_syn); 

  for(syn=0;syn<n_syn;syn++)	/* find first spike for each synapse */
    {
      indnexttime[syn]=0;     
      nexttime[syn]=alltrains[syn*length_train]+offset; /* time of first spike */
    }

  for(t=0;t<=t_tot; t+=dt)
    {
      totact=0;                   /*total activity at this time*/
      for (syn=0;syn<n_syn;syn++)
	{
	  if(fabs(t-nexttime[syn])<dt)	/*this synapse has a spike at time t */
	    {
	      fprintf(outfile,"1"); /* write spike to disk */
	      indnexttime[syn] ++;   /* next spike of this synapse: */
	      nexttime[syn] = offset + alltrains[syn*length_train+indnexttime[syn]];        
              totact ++;            /*one spike more in total activity*/
	      if(nexttime[syn] > t_tot)
		{
		  nexttime[syn] = t_tot+dt;    
                }
            }
	  else			/* synapse has no spike at time t */
	    {
	      fprintf(outfile,"0");
	    }
	  if (syn<(n_syn-1)) {
	    fprintf(outfile," ");
	  }
	}
      if (t < (t_tot-dt)) {
	fprintf(outfile,"\n");
      }

      /*write total activity for this time:*/
      //      fprintf(f_totact,"%g %d\n",t,totact);
    }
  
  /*make sure all spikes have been found (consistency check):*/
  //  for(syn=0;syn<n_syn;syn++)
  // {
  //    if(indnexttime[syn] != length_train)
  //	{
  //	  printf("Working on %s: Someth. wrong w. synapse %d:\n",errfile,syn);
  //	  printf("Last spike: %i != length of spiketrain: %i\n",
  //		 indnexttime[syn],length_train);

  //	  sprintf(consist_filename,"ERROR_%s",errfile);
  //	  problem_file = fopen(consist_filename,"a");
  //
  //	  fprintf(problem_file,"\n\n*******\n\nn_syn: %i\n",n_syn);
  //	  fprintf(problem_file,"length_train: %i\n",length_train);
  //	  fprintf(problem_file,"t_tot: %i\n",t_tot);
  //	  fprintf(problem_file,"syn: %i\n\n",syn);
 
  //	  for (outsyn=0;outsyn<n_syn;outsyn++)
  //	    for(nspike=0;nspike<length_train;nspike++)
  //	      fprintf(problem_file,"%i %i: %g\n",
  //		      outsyn,nspike,alltrains[outsyn*length_train+nspike]);

  //	  fclose (problem_file);
  //	}
  //  }   
}

/********************************************************************/
/*                                                                  */
/* make_spike_train:                                                */
/* add shift to individual spikes                                   */
/* given by gaussian with sigma_s                                   */
/*                                                                  */
/********************************************************************/

void make_spike_trains(double *spiketrain,double *gentrain,double p,\
		       double s,double av_rate,double t_tot,double toff,\
		       int *p_seed)
{
  double t,delta_t;
  double sigma_s;
  int i;
  double alpha;
  
  alpha=2.;    
  delta_t = 1./av_rate; 
  i=0;
  
  /*sigma_s = t_tot;   this is wrong, no reason for t_tot dependence!*/ 
  sigma_s = (1-p)*(1-s)*alpha*delta_t; /* determine sigma of gaussian */

  for(t=0.;t<t_tot;t+=delta_t) /* make spike train with gaussian */
    {
      spiketrain[i] = gentrain[i] + toff + sigma_s*gasdev(p_seed);
      i++;
    }
}

/********************************************************************/
/*                                                                  */
/* make_generator_train:                                            */
/* make master spike train                                          */
/*                                                                  */
/********************************************************************/

void make_generator_train(double *gentrain,double t_tot,\
			  double av_rate,double p,int *p_seed)
{
  double t,delta_t;
  double sigma_g;
  int i;
  
  delta_t = 1./av_rate;
  i=0;

  if(p >APRECISION && p<=1.)		/* finite value for p */
    {
      sigma_g = (1-p)/p * delta_t; /* determine sigma of gaussian */

      for(t=0.;t<t_tot;t+=delta_t) /* make spike train with gaussian */
	{
	  gentrain[i] = t+sigma_g*gasdev(p_seed);
	  i++;
	}
    }
  else if(fabs(p) < APRECISION)	     /* p is zero ==> "infinite sigma" */
    {
      for(t=0.;t<t_tot;t+=delta_t) 
	{
	  gentrain[i] = t_tot * ran1(p_seed);    /* pick random spike times  */
	  i++;
	}
    }
  else				/* error! */
    nrerror("error in make_generator_train: p must be 0 <= p <= 1");
}

      
/********************************************************************/
/*                                                                  */
/* Main.                                                            */
/* Parameters:                                                      */
/* n_syn    - Number of synapses                                    */
/* t_tot    - total simulation time (milliseconds)                  */
/* dt       - simulation time step  (milliseconds)                  */
/* av_rate  - average spike rate  (spikes per second)               */
/* p        - periodicity parameter, 0<p<1                          */
/* s        - synchronicity parameter, 0<s<1                        */
/* offset   - spike temporal offset (milliseconds)                  */
/*                                                                  */
/********************************************************************/

main(int argc, char *argv[])
{
  int n_syn,syn;
  int length_train;
  int start_offset;
  double t_tot,dt,av_rate, offset;
  FILE *outfile,*f_totact,*logfile;
  char outfilename[1000],totactname[1000];
  double *spiketrain, *gentrain, *alltrains;
  double p,s;
  double toff;
  int seed,cond;
  int i;
  double t_ref=1.; /*refractory period, in ms*/

 

  strcpy(errfile,"shiftsyn.err");

  read_input(argc, argv,&n_syn,&t_tot,&dt,&av_rate,&p,&s,&seed,&offset,outfilename,totactname);
  
  // printf( "%g\n", offset);

  OPENFILE(logfile,"shiftsyn.log","w",exit(1));
  fprintf(logfile, "Right now, shiftsyn is working on %s\n",errfile);
  fclose(logfile);
  
  if(seed>50000) {
    printf("seed too large; problem for Num. Recipes random generator\n");
    exit(1);
  }

  if(t_tot*av_rate<=1)
    nrerror("run too short (less than one spike on average), EXIT");
  else
  length_train = floor(t_tot*av_rate);
  spiketrain = alloctimes(length_train);
  gentrain = alloctimes(length_train);
  alltrains =  alloctimes(n_syn*length_train);
  outfile = fopen(outfilename, "w");
 //  f_totact = fopen(totactname, "w");
  
   
 make_generator_train(gentrain,t_tot,av_rate,p,&seed);
  
  for(syn=0;syn<n_syn;syn++)       /* loop over all synapses */
    {
      /* compute global offset for whole spiketrain of this synapse: */
      toff = make_train_offset(p,s,av_rate,t_tot,&seed);
      
      /* add shifts to individual spike trains */
      make_spike_trains(spiketrain,gentrain,p,s,av_rate,t_tot,toff,&seed);
      
      /* sort spiketrain */
      sort_train(spiketrain,length_train,t_tot);
      
      /* enforce refractory period */
      cond=enforce_refrac(spiketrain,syn,length_train,t_ref,dt,t_tot);
 
     
      if (cond == 1)  {             /*normal case*/
	/* append to alltrains:*/
	append_train(spiketrain,alltrains,syn,length_train);
      }
      else  if (cond == -1) {       /*make a new spiketrain for same syn*/
	if(s<0.9) {            
	  syn --;   /*if s!=1 it is enough to restart this synapse */
	}
	else { /*if s==1 need also new gen. train and all new synapses */  
	  make_generator_train(gentrain,t_tot,av_rate,p,&seed);
	  syn=-1;
	}
      }
      else {
	printf("ERROR in main: this cannot happen\n");
	exit(1);
      }
    }
  
  /* All spiketimes are now computed. Write them to disk:*/
  
  write_to_disk(alltrains,outfile,f_totact,length_train,t_tot,n_syn,dt,offset);
  
  /* Debugging output: */
  /*
     for(i=0;i<n_syn*length_train; i++)
     {
     printf("Working on %s: %d %f\n",errfile,i,alltrains[i]); 
     printf("%f \n",alltrains[i]);
     }
     */

/* successfully completed, rm logfile: */
   system("rm -f shiftsyn.log");
}





