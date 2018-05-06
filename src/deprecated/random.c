/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*                              random.c
*
* This is a file containing a series of functions written by John
* Skilling for generating various kinds of random numbers. Charlie 
* lifted these out to plug them into jelly, but if a progrma were to use
* both bayesys and jelly there may be problems with multiple
* definitions. Here, then, is my solution - rename the functions in this
* file sufficiently hard so that it can exist as a standalone utility to
* be called from fortran from anywhere. 
*
* So from fortran:
*
*      integer seed,state(4)
* 	 integer idummy
*      real    rdummy
*      real*8    gdummy
*      real*8  ddummy
*
*      seed = -1
*      call RandomInit(state,seed)
*      idummy = RandomInteger(state)
*      rdummy = 1.0*RandomDouble(state)
*      gdummy = RandomGaussian(state)
*      ddummy = RandomDouble(state)
*
*  Haxxed by pjm in summer 2002
*
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Filename:  bayesubs.c
 * 
 * Purpose:   Random and Peano utility subroutines for BayeSys3.
 * 
 * History:   Random.c  17 Nov 1994 - 28 Jan 2002
 *            Peano.c   10 Apr 2001
 * 
 *            Copyright (c) 1994-2002, Maximum Entropy Data Consultants Ltd.
 *-----------------------------------------------------------------------------
 */
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

typedef unsigned *Rand_t;

//#define DEBUG 1

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function:  RanInit
 *
 * Purpose:   Initialise generator for other Random procedures
 *
 * Method:    Hashing based on ran4 in "Numerical Recipes"
 *
 *            Set 64-bit counter (Rand[1] high, Rand[0] low)  =  (seed, 0)
 *
 * History:   John Skilling   15 Jan 2002
 *-----------------------------------------------------------------------------
 */
int   craninit(         //   O  seed, either from input or time
Rand_t   Rand,         //   O  Random generator state
int      seed)         // I    Seed: +ve = value, -ve = time seed
{
	
    static const unsigned   Z = (unsigned)(-1) >> 1;
    unsigned    j, k;

    k = 1;
    for( j = 0; k; ++j )
        k += k;
    if( seed < 0 )
        seed = (int)(time(NULL) & Z);  // still OK after A.D.2030
    Rand[0] = 0;
    Rand[1] = (unsigned)seed;
    Rand[2] = Z;                       // Rangauss storage switch
#ifdef DEBUG
    fprintf(stdout, "Raninit %d %d\n", seed, Rand[2]);
#endif
    Rand[3] = Rand[1];
    return seed;
}


int   randominit_(         //   O  seed, either from input or time
Rand_t   Rand,         //   O  Random generator state
int      *seed)    
{
	return craninit(Rand, *seed);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function:  Ranint
 *
 * Purpose:   Random integer sample, in [-2^31,2^31).
 *
 * Method:              (64-bit hashing) XOR (simple congruential)
 *            Hashing from ran4 in "Numerical Recipes", with re-seeding.
 *            Simple congruential should make any non-randomness undetectable.
 *
 * Note:      Period = 2^64
 *
 * History:   John Skilling   28 Jan 2002
 *-----------------------------------------------------------------------------
 */
int cranint(          //   O  Value
Rand_t   Rand)         // I O  Random generator state
{
	unsigned i, j, k, m, n;

    n = ++Rand[0];
    m = n ? Rand[1] : ++Rand[1];   // Increment seed after 2^32 calls
    k = n ^ 0xbaa96887;
    j = k >> 16;
    k &= 0xffff;
    i = (j - k) * (j + k);
    m ^= (((i >> 16) | (i << 16)) ^ 0xb4f0c4a7) + k * j;
    k = m ^ 0x1e17d32c;
    j = k >> 16;
    k &= 0xffff;
    i = (j - k) * (j + k);
    n ^= (((i >> 16) | (i << 16)) ^ 0x178b0f3c) + k * j;
    k = n ^ 0x03bcdc3c;
    j = k >> 16;
    k &= 0xffff;
    i = (j - k) * (j + k);
    m ^= (((i >> 16) | (i << 16)) ^ 0x96aa3a59) + k * j;
    k = m ^ 0x0f33d1b2;
    j = k >> 16;
    k &= 0xffff;
    i = (j - k) * (j + k);
    n ^= (((i >> 16) | (i << 16)) ^ 0xaa5835b9) + k * j;
    Rand[3] = Rand[3] * 0x0019660d + 0x3c6ef35f;

//#ifdef DEBUG
//    fprintf(stdout, "RanInt %d\n", (int)(n ^ Rand[3]));
//#endif    
    return (int)(n ^ Rand[3]);
}

int randominteger_(Rand_t   Rand)  
{ 
    return cranint(Rand); 
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function:  Randouble
 *
 * Purpose:   Random double-precision floating point sample, inside (0,1).
 *
 * Method:    Development of Ranint.
 *
 * History:   John Skilling   6 May 1995, 3 Dec 1995, 24 Aug 1996
 *-----------------------------------------------------------------------------
 */
double crandouble(      //   O  Value
Rand_t   Rand)         // I O  Random generator state
{
    static const double STEP = 1.0 / ((unsigned)(-1) + 1.0);
    double   r;
    do  r = ((double)(unsigned)cranint(Rand) + 0.5) * STEP;
    while( r >= 1.0 );

#ifdef DEBUG
   fprintf(stdout, "randouble_ %e\n",  r);
#endif

    return r;
}

double randomdouble_(Rand_t   Rand)  
{ 
    return crandouble(Rand); 
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Function:  Rangauss
 * 
 * Purpose:   Sample from Gaussian N(0,1)
 * 
 * Notes: (1) Requires previous call to RanInit to initialise the generator
 *            stored in Rand[0..1] AND the Rangauss switch in Rand[2]
 *        (2) Maximum outlier extent cannot quite reach 7 standard deviations
 *        (3) Generates 2 deviates every other call
 * 
 * History:   JS         22 Jan 1994   Box-Muller
 *                       19 Oct 1995
 *-----------------------------------------------------------------------------
 */
double crangauss(       //   O  Value
Rand_t   Rand)         // I O  Random generator state
{
    static const unsigned   Z = (unsigned)(-1) >> 1;
    static const double     F = (double)(((unsigned)(-1) >> 4) + 1);
    double a, r, g1, g2;

	double ret = 0.0;
	
	//fprintf(stdout, "Rand2 %d\n", Rand[2]);

    if( Rand[2] == Z )
    {
// Generate 2 normal deviates (g1*a) and (g2*a)
        do
        {
            g1 = 2.0 * crandouble(Rand) - 1.0;                      // (-1,1)
            g2 = 2.0 * crandouble(Rand) - 1.0;                      // (-1,1)
            r = g1 * g1 + g2 * g2;
        } while( r >= 1.0 );
        a = sqrt(-2.0 * log(r) / r);                       // 0 < a < 3.0e10
// Store (int) representation of 2nd normal deviate
//                           (|g2*a|<7) * 2^28  is well within  (-2^31,2^31)
        Rand[2] = (unsigned)((int)(g2 * a * F));
	//fprintf(stdout, "rangauss_ %p %e\n", Rand, g1*a);
       ret = (g1 * a);
    }
    else
    {
// Recover 2nd normal deviate from (int) representation
        g2 = (double)((int)Rand[2]) / F;
        Rand[2] = Z;              // Outside range of (int) stored deviate
	//
        ret = g2;
    }

#ifdef DEBUG
  	fprintf(stdout, "rangauss_ %e \n", ret);
#endif

//    if (!finite(ret))
  //  	return rangauss_(Rand);
	//else 
		return ret;
}

double randomgaussian_(Rand_t   Rand)  
{ 
	return crangauss(Rand); 
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
