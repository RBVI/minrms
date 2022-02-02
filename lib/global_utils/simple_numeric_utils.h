#ifndef _SIMPLE_NUMERIC_UTILS_H
#define _SIMPLE_NUMERIC_UTILS_H

#include <cstdlib>
#include <cmath>
#include <cassert>
using namespace std;

#include "simple_numeric_types.h"


template<class T>
inline T SQR(T x)  { return x*x; }

template<class T>
inline T ABS(T x)  { return ( ((x) > 0) ? (x) : -(x) ); }

template<class T>
inline T MAX(T a, T b) { return ( ((a) > (b)) ? (a) : (b) ); }

template<class T>
inline T MIN(T a, T b) { return ( ((a) < (b)) ? (a) : (b) ); }



// Pick the correct version of floor() or ceil to use:
// (decision should happen at compile time)

#define FLOOR(x)        \
        ( sizeof (x) == sizeof(float )        ?        floorf(x) \
        : sizeof (x) == sizeof(double)        ?        floor(x)  \
        : floorl(x) )

#define CEIL(x)        \
        ( sizeof (x) == sizeof(float )        ?        ceilf(x)  \
        : sizeof (x) == sizeof(double)        ?        ceil(x)   \
        : ceill(x))

// All compilers agree that floor(1.5) is 1.0.
// On some (but not all!!) compilers floor(-1.5) returns -1.0.
// UGH. This bugs the crap out of me.  In my humble
// opinion, I consider this broken.  I think it should return -2.0;
#define FLOOR_IS_BROKEN = (floor(-1.5) > -1.5)

// Similarly, on some compilers ceil(-1.5) returns -2.0.
#define CEIL_IS_BROKEN = (ceil(-1.5) < -1.5)

// The "fixed" versions of floor() and ceil()
// are type independent, and 

// Don't use this: It ends up evaluating X twice.
//#define ROUND_DOWN(X) {(X<0)?-static_cast<Integer>(-X):static_cast<Integer>(X)}


template<class T>
inline T ROUND_DOWN(T x)
{
  T x_rounded = FLOOR(x);
  return ((x_rounded > x) ? (x_rounded-1.0) : x_rounded);
}


template<class T>
inline T ROUND_UP(T x)
{
  T x_rounded = CEIL(x);
  return ((x_rounded < x) ? (x_rounded+1.0) : x_rounded);
}



template<class T>
inline T ROUND_NEAREST(T x)
{ return ROUND_DOWN(x + 0.5); }


// special cases (examples)
// ROUND_NEAREST_CHEAP(4.5) = 5, but
// ROUND_NEAREST_CHEAP(-4.5) = -5  (not -4)

template<class T>
inline T ROUND_NEAREST_CHEAP(T x)
{ 
  if (x > 0.0)
    return FLOOR(x+0.5);
  else
    return FLOOR(x-0.5);
}  




//Truncates an integer to a multiple of another integer (always rounding down)
inline long long int ROUND_DOWN_DIV(long long n, long long d)
{ return
    ((n < 0) && (n%d != 0))
    ? n/d-1
    : n/d;
}



//Returns what's left over after truncation
inline long long ROUND_DOWN_REMAINDER(long long n, long long d)
{ long long answer = n%d;
 if ((n < 0) && (answer != 0))
   answer += d;
 assert((0 <= answer) && (answer < d));
 return answer;
}



//Truncates an integer to a multiple of another integer (always rounding up)
inline long long int ROUND_UP_DIV(long long n, long long d)
{
  return
    ((n%d == 0)
     ?
     n/d
     :
     ((n>0) ? (n/d)+1 : n/d));
}





template<class T>
inline bool EqualToTolerance(T a, T b, T tolerance)
{ return (((a+b) == 0.0) ?
          ((a==b) ? true : false)
          :
          ((2.0*ABS((a - b) / (a + b))) < tolerance)); }



template<class T>
inline bool EqualTo1PartIn10_4(T a, T b)
{ return EqualToTolerance(a, b, (T)1e-4); }



template<class T>
inline bool EqualTo1PartIn10_13(T a, T b)
{ return EqualToTolerance(a, b, (T)1e-13); }



// The next function is a hack I wrote.
// I'm not sure if it is useful.
template<class T>
inline bool EqualToToleranceOrZero(T a, T b, T tolerance)
{ return (((a == 0) || (b == 0))
          ?
          true
          :
          (((a+b) == 0.0)
           ?
           false
           :
           ((2.0*ABS((a - b) / (a + b))) < tolerance)
           )
          );
}



template<class T>
inline bool EqualTo1PartIn10_4orZero(T a, T b)
{ return EqualToToleranceOrZero(a, b, (T)1e-4); }



template<class T>
inline bool EqualTo1PartIn10_13orZero(T a, T b)
{ return EqualToToleranceOrZero(a, b, (T)1e-13); }


// Choose the seed for the random-number generator
inline void INIT_RANDOM_SEED(long n) { srand48( n ); }


//RANDOM_INT(n) returns uniformly distributed integers ranging from 0 to n-1.
//The maximum allowed value of "n" is RAND_MAX, which (I think) is usually 2^31

inline long RANDOM_INT(unsigned long n)
{
  assert(n < RAND_MAX);
  return lrand48() % n;
}


//RANDOM_REAL_0_1() returns uniformly distributed "Real" numbers over [0, 1)

inline double RANDOM_REAL_0_1()
{
  return drand48();
}




// RANDOM_EXPONENTIAL() returns a number in the range [0,infinity)
// with an exponential probability distribution  P(x) = exp(-x)
//    (Scale the final answer to get different decay rates.
//     P(x) = (1/K) exp(-K x))


inline double RANDOM_EXPONENTIAL()
{
  double u = RANDOM_REAL_0_1(); // generates a random number in the range [0,1)
  u = 1.0 - u;                  // After this command, u is in the range (0,1]
                                // We want u>0 because log(u) not defined at 0
  return -log(u);
}



//RANDOM_GAUSSIAN() returns a random number in the range
//                            (-infinity, +infinity)
//with distrubition of P(x) proportional to exp( - x^2 / 2 )
//   (the "standard normal distribution".)
//The number should be Gaussian distributed, and have a variance, <x^2> = 1.

inline double RANDOM_GAUSSIAN()
{
  double u = RANDOM_REAL_0_1(); // generates a random number in the range [0,1)
  u = 1.0 - u;                  // After this command, u is in the range (0,1]
                                // We want u>0 because log(u) not defined at 0

  // The following line should generate r with a probability density of
  // P(r) proportional to  r * exp(- r^2/2)

  double r =  sqrt( - log (u) * 2.0 );
  //cerr << "u=" << u << ", log(u)=" << log(u) << ", sqrt(-log(u))="
  //     << sqrt(-log(u)) << endl;


  // Now generate another random number between 0 and 2*pi
  double theta = (2.0*M_PI) * RANDOM_REAL_0_1();

  // Convert from polar (r,theta) to cartesian coordinates (x,y)
  // The "x" and "y" values should both be Gaussian distributed.
  // (Probability density P(x) should be proportional to exp(-lambda x^2))

  //double return_val = r*cos(theta);
  //cerr << "    theta=" << theta << ", cos(theta)=" << cos(theta)
  //     << ", r=" << r << ", r*cos(theta)=" << return_val << endl;

  return r * cos( theta );

} //RANDOM_GAUSSIAN()





//RANDOM_GAUSSIAN_VERSION_COEFF_NUMERATOR()
//returns a random number in the range (-infinity, +infinity)
//with distrubition of P(x) proportional to exp( - A x^2 )
//The numbers should be Gaussian distributed.

inline double RANDOM_GAUSSIAN_COEFF_NUMERATOR(double A)
{
  double u = RANDOM_REAL_0_1(); // generates a random number in the range [0,1)
  u = 1.0 - u;                  // After this command, u is in the range (0,1]
                                // We want u>0 because log(u) not defined at 0

  // The following line should generate r with a probability density of
  // P(r) proportional to  r * exp(-A r^2)

  double r =  sqrt( - log (u) / A );
  //cerr << "A=" << A << endl;
  //cerr << "u=" << u << ", log(u)=" << log(u) << ", sqrt(-log(u))="
  //     << sqrt(-log(u)) << endl;


  // Now generate another random number between 0 and 2*pi
  double theta = (2.0*M_PI) * RANDOM_REAL_0_1();

  // Convert from polar (r,theta) to cartesian coordinates (x,y)
  // The "x" and "y" values should both be Gaussian distributed.
  // (Probability density P(x) should be proportional to exp(-lambda x^2))


  //double return_val = r*cos(theta);
  //cerr << "    theta=" << theta << ", cos(theta)=" << cos(theta) << ", r=" << r << ", r*cos(theta)=" << return_val << endl;

  return r * cos( theta );

} //RANDOM_GUASSIAN_VERSION_COEFF_NUMERATOR()



//RANDOM_GAUSSIAN_COEFF_DENOMINATOR() returns a random number in the range
//(-infinity, +infinity)
//with distrubition of P(x) proportional to exp( - x^2 / B )
//The numbers should be Gaussian distributed.

inline double RANDOM_GAUSSIAN_COEFF_DENOMINATOR(double B)
{
  double u = RANDOM_REAL_0_1(); // generates a random number in the range [0,1)
  u = 1.0 - u;                  // After this command, u is in the range (0,1]
                                // We want u>0 because log(u) not defined at 0

  // The following line should generate r with a probability density of
  // P(r) proportional to  r * exp(- r^2/B)

  double r =  sqrt( - log (u) * B );
  //cerr << "B=" << B << endl;
  //cerr << "u=" << u << ", log(u)=" << log(u) << ", sqrt(-log(u))="
  //     << sqrt(-log(u)) << endl;


  // Now generate another random number between 0 and 2*pi
  double theta = (2.0*M_PI) * RANDOM_REAL_0_1();

  // Convert from polar (r,theta) to cartesian coordinates (x,y)
  // The "x" and "y" values should both be Gaussian distributed.
  // (Probability density P(x) should be proportional to exp(-lambda x^2))

  //double return_val = r*cos(theta);
  //cerr << "    theta=" << theta << ", cos(theta)=" << cos(theta)
  //     << ", r=" << r << ", r*cos(theta)=" << return_val << endl;

  return r * cos( theta );

} //RANDOM_GUASSIAN_COEFF_DENOMINATOR()







/* FindMedian()
 * finds the median value in a set of integers (rounded up).
 * If the set is the empty set, it returns n_min-1.
 * Because of the funny way the "set" is passed to the function,
 * I doubt anyone else will have a use for this function.
 * The set itself should entirely fit within the interval [n_min, n_max].
 * The "aNinSet" array indicates whether a given integer
 * in the range [n_min, n_max] belongs to the set.
 * An integer "n" belongs to the set iff (aNinSet[n-n_min] == true).
 * The aNinSet array should be large enough to store a flag
 * for each "n" in the whole interval from n_min to n_max
 *  (that is, it should contain (1 + n_max - n_min) entries). */

inline int FindMedian(bool *aNinSet, int n_min, int n_max)
{
  int n;
  int num_in_set = 0;
  for(n= n_min;
      (n <= n_max);
      ++n)
  {
    if (aNinSet[n-n_min])
      ++num_in_set;
  }
  if (num_in_set == 0)
    return (n_min - 1);

  int half_num_in_set = (num_in_set+1)/2;
                                         
  n = n_min - 1;
  int count = 0;
  do {
    ++n;
    if (aNinSet[n-n_min])
      ++count;
  } while (count < half_num_in_set);

  assert(aNinSet[n-n_min] == true);
  return n;
} /* FindMedian() */


/* AccumDist()
 *  When AccumDist() is called, the interval from "x_min" to "x_max"
 * is partitioned into "num_bins" equal-sized subintervals.
 * The array of integers "aBins" stores counters, one for each sub-interval.
 * If "x" lies in the range x_min <= x < x_max, then this function finds the
 * sub-interval enclosing x, and increments the counter in that
 * sub-interval's bin.
 *     Returns:
 * The function returns true if x lies inside the interval [x_min,x_max]
 * (if one of the counters was incremented), and false otherwise. */

inline bool
AccumDist(double x,
          long *aBins,
          long  num_bins,
          double x_min,
          double x_max)
{
  assert( x_min != x_max );
  assert(num_bins >= 1);
  int which_bin =
    static_cast<int>(floor(((x-x_min)/(x_max-x_min))
                           *
                           static_cast<double>(num_bins)));
  if ((which_bin >= 0) && (which_bin < num_bins)) {
    ++(aBins[which_bin]);
    return true;
  }
  else return false;
}



// The "random_shuffle" STL algorithm does not suffice,
// because they do not supply a method to specify the seed
// of the default random number generator.
template<class RandomAccessIterator>
inline void
RANDOM_SHUFFLE(RandomAccessIterator pA,
               RandomAccessIterator pB)
{
  for (RandomAccessIterator pI = pA; pI != pB ; ++pI)
  {
    RandomAccessIterator pJ = RANDOM_INT((pI-pA)+1) + pA;
    swap(*pI, *pJ);
  }// for (i = 0 to n-1)
}// inline Shuffle()




/*
//(Not needed anymore.)
// RANDOM_SHUFFLE() takes an array of objects of type T, and a size parameter
// and randomly reorders its contents.
template<class T>
inline void
RANDOM_SHUFFLE(T *aT,             //the array
               unsigned long n)   //size of the array.
{
  assert(n >= 2); //must have at least 2 elements
  for (int i = 0; i < n ; ++i)
  {
    int temp;
    int j = RANDOM_INT(i+1);
    // aT[i] <--swap--> aT[j]
    temp = aT[i];
    aT[i] = aT[j];
    aT[j] = temp;
  }// for (i = 0 to n-1)
}// inline RANDOM_SHUFFLE()
*/



#endif //#ifndef _SIMPLE_NUMERIC_UTILS_H

