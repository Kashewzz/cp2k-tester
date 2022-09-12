

#include <algorithm>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <cmath>

#include "f77blas.h"

#define MAX_M 64
#define MAX_N 64
#define CACHELINE 64
#define rzero   0.0
#if (CACHELINE == 64)
  #define MulByCachelen(N_) ( (N_) << 6 )
  #define DivByCachelen(N_) ( (N_) >> 6 )
#endif
#define ALIGNMENT CACHELINE

/** >1: number of locks, =1: omp critical, =0: atomic */
#define CP2K_SYNCHRONIZATION 0
// ensures sufficient parallel slack
#define CP2K_MIN_NPARALLEL 240
// ensures amortized atomic overhead
#define CP2K_MIN_NLOCAL 160

#if !defined(ITYPE)
# define ITYPE double
#endif

#if !defined(MAX_SIZE)
# define MAX_SIZE ((MAX_M) * (MAX_N))
#endif

unsigned int isqrt_u64(unsigned long long x)
{
  unsigned long long b; unsigned int y = 0, s;
  for (s = 0x80000000/*2^31*/; 0 < s; s >>= 1) {
    b = y | s; y |= (b * b <= x ? s : 0);
  }
  return y;
}




#define UPDIV(N, MULT) (((N) + ((MULT) - 1)) / (MULT))
#define MAX(A, B) ((A) < (B) ? (B) : (A))
#define ABS(A) (0 <= (A) ? (A) : -(A))
#define UP2(N, NPOT) (((N) + ((NPOT) - 1)) & ~((NPOT) - 1))
#define ALIGN(POINTER, ALIGNMENT/*POT*/) ((POINTER) + (UP2((uintptr_t)(POINTER), ALIGNMENT) - ((uintptr_t)(POINTER))) / sizeof(*(POINTER)))
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define DELTA(T0, T1) ((T0) < (T1) ? ((T1) - (T0)) : ((T0) - (T1)))
#define FEQ(A, B) ((A) == (B))

#define AlignPtr(vp) \
   (void*) (CACHELINE + MulByCachelen(DivByCachelen((size_t) (vp))))

void dzero (const int N, double *X, const int incX)
/*
 * X <- 0
 */
{
   int i, n;
   if (incX == 1)
   {
	  n = N; //SHIFT;
	  i = n >> 5;
	  if (i)
	  {
		 n -= (i << 5);
		 do
		 {
			*X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = X[7] = X[8] = X[9] =
			X[10] = X[11] = X[12] = X[13] = X[14] = X[15] = X[16] = X[17] =
			X[18] = X[19] = X[20] = X[21] = X[22] = X[23] = X[24] = X[25] =
			X[26] = X[27] = X[28] = X[29] = X[30] = X[31] = rzero;
			X += 32;
		 }
		 while(--i);
	  }
	  if (n >> 4) /* >= 16 */
	  {
		 *X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = X[7] = X[8] = X[9] =
		 X[10] = X[11] = X[12] = X[13] = X[14] = X[15] = rzero;
		 X += 16;
		 n -= 16;
	  }
	  if (n >> 3) /* >= 8 */
	  {
		 *X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = X[7] = rzero;
		 X += 8;
		 n -= 8;
	  }
	  switch(n)
	  {
		 case 1:
			*X = rzero;
			break;
		 case 2:
			*X = X[1] = rzero;
			break;
		 case 3:
			*X = X[1] = X[2] = rzero;
			break;
		 case 4:
			*X = X[1] = X[2] = X[3] = rzero;
			break;
		 case 5:
			*X = X[1] = X[2] = X[3] = X[4] = rzero;
			break;
		 case 6:
			*X = X[1] = X[2] = X[3] = X[4] = X[5] = rzero;
			break;
		 case 7:
			*X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = rzero;
			break;
		 default:;
	  }
   }
   else
   {
     #ifdef TREAL
		 for (i=N; i; i--, X += incX) *X = rzero;
     #else
		 for (n=incX<<1, i=N; i; i--, X += n) *X = X[1] = rzero;
     #endif
   }
}
double flushcache(int size)
/*
 * flush cache by reading enough mem; note that if the compiler gets
 * really smart, may be necessary to make vp a global variable so it
 * can't figure out it's not being modified other than during setup;
 * the fact that ATL_dzero is external will confuse most compilers
 */
{
  static void *vp=NULL;
  static int N = 0;
  double *cache;
  double dret=0.0;
  int i;

  if (size < 0) /* flush cache */
  {
     assert(vp);
     cache = (double*) AlignPtr(vp);
     if (N > 0) for (i=0; i != N; i++) dret += cache[i];
  } else if (size > 0) {
     vp = malloc(CACHELINE + size);
     assert(vp);
     N = size / sizeof(double);
     cache = (double*) AlignPtr(vp);
     dzero(N, cache, 1);
  }
  return(dret);
}

template<typename T>
void add(T* dst, const T* src, int nrows, int ncols, int ld_src = 0)
{
  const int ld = (0 == ld_src ? ncols : ld_src);
  {
    for (int i = 0; i < nrows; ++i) {
      for (int j = 0; j < ncols; ++j) {
        const T value = src[i*ld+j];
        dst[i*ncols+j] += value;
      }
    }
  }
}


unsigned long long timer_tick_rtc(void)
{
  unsigned long long result;
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  result = 1000000000ULL * t.tv_sec + t.tv_nsec;
  return result;
}

unsigned long long timer_tick(void)
{
  unsigned long long result;

  result = timer_tick_rtc();

  return result;
}

double timer_duration_rtc(unsigned long long tick0, unsigned long long tick1)
{
  double result = (double)DELTA(tick0, tick1);
  result *= 1E-9;
  return result;
}

double timer_duration(unsigned long long tick0, unsigned long long tick1)
{
  double result;
  result = timer_duration_rtc(tick0, tick1);
  return result;
}


int main(int argc, char* argv[])
{
 

    // flush cache
    double ret = flushcache(1);

    typedef ITYPE T;
    int m = 1 < argc ? std::atoi(argv[1]) : 23;
    const int q = ((1ULL << 30) / (3 * m * m * sizeof(T)));
    const int r = 2 < argc ? (0 < std::atoi(argv[2]) ? std::atoi(argv[2]) : ('+' == *argv[2]
      ? (q << std::strlen(argv[2])) : ('-' == *argv[2]
      ? (q >> std::strlen(argv[2])) : 0))) : 0;
    const int t = 3 < argc ? (0 < std::atoi(argv[3]) ? std::atoi(argv[3]) : ('+' == *argv[3]
      ? ((CP2K_MIN_NLOCAL) << std::strlen(argv[3])) : ('-' == *argv[3]
      ? ((CP2K_MIN_NLOCAL) >> std::strlen(argv[3])) : -1))) : -1;
    int k = 5 < argc ? std::atoi(argv[5]) : m;
    int n = 4 < argc ? std::atoi(argv[4]) : k;
    char transa = 'N', transb = 'N';
    ITYPE alpha = 1, beta = 1;

    const int csize = m * n;
    if ((MAX_SIZE) < csize) {
      throw "The size M x N is exceeding MAX_SIZE!";
    }

    const int asize = m * k, bsize = k * n, aspace = ALIGNMENT / sizeof(T);
    const int s = 0 < r ? r : ((2ULL << 30) / ((asize + bsize) * sizeof(T))); // 2 GByte
    const int u = 0 < t ? t : static_cast<int>(isqrt_u64(s * CP2K_MIN_NLOCAL / CP2K_MIN_NPARALLEL));
    const size_t bwsize = static_cast<size_t>((s * (asize + bsize)/*load*/ + UPDIV(s, u) * csize * 2/*accumulate*/) * sizeof(T));
    const double gflops = 2.0 * s * m * n * k * 1E-9, scale = 1.0 / s;
    const char ops[] = "FLOPS";
    const char *const env_check = getenv("CHECK");
    const double check = ABS(NULL == env_check ? 0 : atof(env_check));

    struct raii { // avoid std::vector (first-touch init. causes NUMA issue)
      T *a, *b, *c;
      raii(int asize_, int bsize_, int csize_)
        : a(new T[static_cast<size_t>(asize_)]), b(new T[static_cast<size_t>(bsize_)])
        , c(new T[static_cast<size_t>(csize_)]) {}
      ~raii() { delete[] a; delete[] b; delete[] c; }
    } buffer(s * asize + aspace - 1, s * bsize + aspace - 1, csize);
    T *const a = ALIGN(buffer.a, ALIGNMENT);
    T *const b = ALIGN(buffer.b, ALIGNMENT);
    T * /*const*/ c = buffer.c; // no alignment, but thread-local array will be aligned


  // Initialize a and b matrix
    for (int i = 0; i < s; ++i) {
        int NCOLS = k;
        int LD = m;
        for (int ii = 0; ii < NCOLS; ++ii) { 
            for (int jj = 0; jj < LD; ++jj) { 
                const int kk = ii * LD + jj; 
                ((ITYPE*)(a + i * asize))[kk] = 0.5 - ((double)rand())/((double)RAND_MAX); 
            } 
        } 
    }

    for (int i = 0; i < s; ++i) {
        int NCOLS = n;
        int LD = k;
        for (int ii = 0; ii < NCOLS; ++ii) { 
            for (int jj = 0; jj < LD; ++jj) { 
                const int kk = ii * LD + jj; 
                ((ITYPE*)(b + i * bsize))[kk] = 0.5 - ((double)rand())/((double)RAND_MAX); 
            } 
        } 
    }

    {
      fprintf(stdout, "m=%lli n=%lli k=%lli size=%lli memory=%.1f MB (%s)\n\n",
        static_cast<long long>(m), static_cast<long long>(n), static_cast<long long>(k), static_cast<long long>(s),
        1.0 * (s * (asize + bsize) * sizeof(T)) / (1 << 20), 8 == sizeof(T) ? "DP" : "SP");

      struct raii_expect { // avoid std::vector (first-touch init. causes NUMA issue)
        T *expect;
        explicit raii_expect(int size): expect(0 < size ? new T[static_cast<size_t>(size)] : 0) {}
        ~raii_expect() { delete[] expect; }
      } expect_buffer(FEQ(0, check) ? 0 : csize);
      T *const expect = (0 == expect_buffer.expect ? c : expect_buffer.expect);
      const T zero = 0;

      { // LAPACK/BLAS3 (warmup BLAS Library)
        std::fill_n(expect, csize, zero);

        for (int i = 0; i < s; i += u) {
          T tmp[MAX_SIZE] = { 0 }; // make sure that stacksize is covering the problem size
          const T *ai = a + i * asize, *bi = b + i * bsize;
          for (int j = 0; j < MIN(u, s - i); ++j) {
            const T *const aij = ai + asize, *const bij = bi + bsize;
            // dgemm_(&transa, &transb, m, n, k,
            //   &alpha, ai, &m, bi, &k, &beta, tmp, &m);
            ai = aij;
            bi = bij;
          }
          add(expect, tmp, m, n); // atomic
        }
      }

      { // LAPACK/BLAS3 (reference)
        // fprintf(stdout, "LAPACK/BLAS...\n");
        std::fill_n(c, csize, zero);
        const unsigned long long start = timer_tick();
        for (int i = 0; i < s; i += u) {
          T tmp[MAX_SIZE] = { 0 }; // make sure that stacksize is covering the problem size
          T *ai = a + i * asize, *bi = b + i * bsize;
          for (int j = 0; j < MIN(u, s - i); ++j) {
            T *aij = ai + asize, *bij = bi + bsize;
            dgemm_(&transa, &transb, &m, &n, &k,
              &alpha, ai, &m, bi, &k, &beta, tmp, &m);
            ai = aij;
            bi = bij;
          }
          add(c, tmp, m, n); // atomic
        }
        const double duration = timer_duration(start, timer_tick());
        char* form;
        form = "%4d   %4d   %4d   %5.3f        %5.3f         %5.2f      %5.3f\n";
        (void) fprintf( stdout, "%s", "\n");
        (void) fprintf( stdout, "%s%s",
                        "   M      N      K   ",
                        " GFLOP/s   Bandwidth[GB/s]   calls/s[Hz]    time[ms]\n" );
        (void) fprintf(stdout, form, m, n, k, gflops / duration, bwsize / (duration * (1 << 30)), s / duration, 1000.0 * duration);
        (void) fprintf( stdout, "%s", "\n");
        // if (0 < duration) {
        //   fprintf(stdout, "\tperformance: %.1f G%s/s\n", gflops / duration, ops);
        //   fprintf(stdout, "\tbandwidth: %.1f GB/s\n", bwsize / (duration * (1 << 30)));
        //   fprintf(stdout, "\tcalls/s: %.0f Hz\n", s / duration);
        // }
        // fprintf(stdout, "\tduration: %.0f ms\n", 1000.0 * duration);
      }

    }
  
  return 0;
}
