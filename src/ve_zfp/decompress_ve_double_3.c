#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <inttypes.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _FTRACE
#include <ftrace.h>
#endif

#include "zfp.h"
#include "ve_common.h"

#define EBITS      11                       /* number of exponent bits */
#define NBMASK     0xaaaaaaaaaaaaaaaaUL     /* negabinary mask         */
#define DIMS       3
#define BLOCK_SIZE (1 << (2 * DIMS))        /* values per block        */
#define EBIAS      ((1 << (EBITS - 1)) - 1) /* exponent bias           */

#define UInt      uint64_t
#define Int       int64_t
#define Scalar    double
#define FABS      fabs
#define FLOOR     floor
#define LOG2      log2
#define POW       pow
#define INV_XFORM inv_xform_3
#define PERM      perm_3

#include "ve_bitstream.c"
#include "inv_xform_3.c"
#include "scatter_vl_3.c"
#include "decode_block_vl.c"

static void decompress_ve_double_3_vl(zfp_stream* stream, zfp_field* field, uint64_t bmin, uint64_t bmax)
{
  const int vl = (bmax - bmin);
  double fblock[vl * BLOCK_SIZE] __attribute__((aligned(8)));
  bitstream lstream; // Local bitstream for this thread.

  lstream.bits = 0;
  lstream.buffer = 0;
  lstream.ptr = (uint64_t *) ((uint8_t*) stream->stream->begin + bmin * (stream->maxbits / CHAR_BIT));;
  lstream.begin = lstream.ptr;

  decode_block_vl(stream, &lstream, fblock, bmin, bmax);

  scatter_vl_3(field, fblock, bmin, bmax);
}

void decompress_ve_double_3(zfp_stream* stream, zfp_field* field)
{
  uint64_t i, nx, ny, nz, bx, by, bz, blocks, bmin, bmax;
  int it, nt;

  nx = field->nx; bx = (nx + 3) / 4;
  ny = field->ny; by = (ny + 3) / 4;
  nz = field->nz; bz = (nz + 3) / 4;
  blocks = bx * by * bz;

  if (zfp_stream_compression_mode(stream) != zfp_mode_fixed_rate) {
    fprintf(stderr, "%s:%s: only fixed rate mode is supported.\n", __FILE__, __FUNCTION__);
    return;
  }
#ifdef _FTRACE
    ftrace_region_begin("DECOMP");
#endif
  #pragma omp parallel private(i,nt,it,bmin,bmax)
  {
    #ifdef _OPENMP
      nt = omp_get_num_threads();
      it = omp_get_thread_num();
    #else
      nt = 1;
      it = 0;
    #endif

    bmin = blocks * it / nt;
    bmax = blocks * (it + 1) / nt;

    for (i = bmin; i < bmax; i += VL) {
      uint64_t n = (bmax - i) < VL ? (bmax - i) : VL;
      decompress_ve_double_3_vl(stream, field, i, i + n);
    }
  }
#ifdef _FTRACE
    ftrace_region_end("DECOMP");
#endif

  stream->stream->ptr = (uint64_t *) ((uint8_t*) stream->stream->begin + blocks * (stream->maxbits / CHAR_BIT));
}
