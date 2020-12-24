/*  Header for xoroshiro128+ 

Written in 2016-2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */
// #ifndef XOROSHIRO128PLUS_H
// #define XOROSHIRO128PLUS_H

#include <stdint.h>

/* This is xoroshiro128+ 1.0, our best and fastest small-state generator
   for floating-point numbers. We suggest to use its upper bits for
   floating-point generation, as it is slightly faster than
   xoroshiro128**. It passes all tests we are aware of except for the four
   lower bits, which might fail linearity tests (and just those), so if
   low linear complexity is not considered an issue (as it is usually the
   case) it can be used to generate 64-bit outputs, too; moreover, this
   generator has a very mild Hamming-weight dependency making our test
   (http://prng.di.unimi.it/hwd.php) fail after 5 TB of output; we believe
   this slight bias cannot affect any application. If you are concerned,
   use xoroshiro128** or xoshiro256+.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. 

   NOTE: the parameters (a=24, b=16, b=37) of this version give slightly
   better results in our test than the 2016 version (a=55, b=14, c=36).
*/


// static uint64_t s[2];

uint64_t next(uint64_t *s);

uint64_t generate_random_index(uint64_t *state, double maxindexd);

/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

// void jump(void);


//  This is the long-jump function for the generator. It is equivalent to
//    2^96 calls to next(); it can be used to generate 2^32 starting points,
//    from each of which jump() will generate 2^32 non-overlapping
//    subsequences for parallel distributed computations. 

// void long_jump(void);
// // #endif;
