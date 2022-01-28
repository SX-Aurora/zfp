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
#define DIMS       4
#define BLOCK_SIZE (1 << (2 * DIMS))        /* values per block        */
#define EBIAS      ((1 << (EBITS - 1)) - 1) /* exponent bias           */

#define UInt      uint32_t
#define Int       int32_t
#define Scalar    float
#define FABS      fabsf
#define FLOOR     floorf
#define LOG2      log2f
#define POW       powf
#define FWD_XFORM fwd_xform_4
#define PERM      perm_4

#include "ve_bitstream.c"
#include "fwd_xform_4.c"
#include "gather_vl_4.c"
#include "encode_block_vl.c"

static void compress_ve_float_4_vl(zfp_stream* stream, const zfp_field* field, uint64_t bmin, uint64_t bmax)
{
  const int vl = (bmax - bmin);
  float fblock_in[vl * BLOCK_SIZE];
  bitstream lstream; // Local bitstream for this thread.

  lstream.bits = 0;
  lstream.buffer = 0;
  lstream.ptr = (uint64_t *) ((uchar*) stream->stream->begin + bmin * (stream->maxbits / CHAR_BIT));;
  lstream.begin = lstream.ptr;

  gather_vl_4(field, fblock_in, bmin, bmax);

  encode_block_vl(stream, &lstream, fblock_in, bmin, bmax);
}

void compress_ve_float_4(zfp_stream* stream, const zfp_field* field)
{
  uint64_t i, nx, ny, nz, nw, bx, by, bz, bw, blocks, bmin, bmax;
  int it, nt;

  nx = field->nx; bx = (nx + 3) / 4;
  ny = field->ny; by = (ny + 3) / 4;
  nz = field->nz; bz = (nz + 3) / 4;
  nw = field->nw; bw = (nw + 3) / 4;
  blocks = bx * by * bz * bw;

  if (zfp_stream_compression_mode(stream) != zfp_mode_fixed_rate) {
    fprintf(stderr, "%s:%s: only fixed rate mode is supported.\n", __FILE__, __FUNCTION__);
    return;
  }

#ifdef _FTRACE
    ftrace_region_begin("COMP");
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
      compress_ve_float_4_vl(stream, field, i, i + n);
    }
  }
#ifdef _FTRACE
    ftrace_region_end("COMP");
#endif

  stream->stream->ptr = (uint64_t *) ((uchar*) stream->stream->begin + blocks * (stream->maxbits / CHAR_BIT));
}
