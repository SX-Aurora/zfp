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

#define EBITS      8                        /* number of exponent bits */
#define NBMASK     0xaaaaaaaau              /* negabinary mask         */
#define DIMS       2
#define BLOCK_SIZE (1 << (2 * DIMS))        /* values per block        */
#define EBIAS      ((1 << (EBITS - 1)) - 1) /* exponent bias           */

#define UInt      uint32_t
#define Int       int32_t
#define Scalar    float
#define FABS      fabsf
#define FLOOR     floorf
#define LOG2      log2f
#define POW       powf
#define INV_XFORM inv_xform_2
#define PERM      perm_2

#include "ve_bitstream.c"
#include "inv_xform_2.c"
#include "scatter_vl_2.c"
#include "decode_block_vl.c"

static void decompress_ve_float_2_vl(zfp_stream* stream, zfp_field* field, uint64_t bmin, uint64_t bmax)
{
  const int vl = (bmax - bmin);
  float fblock[vl * BLOCK_SIZE] __attribute__((aligned(8)));
  bitstream lstream; // Local bitstream for this thread.

  lstream.bits = 0;
  lstream.buffer = 0;
  lstream.ptr = (uint64_t *) ((uint8_t*) stream->stream->begin + bmin * (stream->maxbits / CHAR_BIT));;
  lstream.begin = lstream.ptr;

  decode_block_vl(stream, &lstream, fblock, bmin, bmax);

  scatter_vl_2(field, fblock, bmin, bmax);
}

void decompress_ve_float_2(zfp_stream* stream, zfp_field* field)
{
  uint64_t i, nx, ny, bx, by, blocks, bmin, bmax;
  int it, nt;

  nx = field->nx; bx = (nx + 3) / 4;
  ny = field->ny; by = (ny + 3) / 4;
  blocks = bx * by;

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
      decompress_ve_float_2_vl(stream, field, i, i + n);
    }
  }
#ifdef _FTRACE
    ftrace_region_end("DECOMP");
#endif

  stream->stream->ptr = (uint64_t *) ((uint8_t*) stream->stream->begin + blocks * (stream->maxbits / CHAR_BIT));
}
