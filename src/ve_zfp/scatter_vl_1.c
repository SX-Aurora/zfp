static void scatter_vl_1(const zfp_field* field, Scalar *fblock, uint64_t bmin, uint64_t bmax)
{
  const int vl = (bmax - bmin);

  Scalar *data = (Scalar *) field->data;

  const uint64_t nx = field->nx;
  const uint64_t bx = (nx + 3) / 4;
  const uint64_t sx = field->sx ? field->sx : 1;
  uint64_t blk_idx;
  int partial_scatter = 0;

  // TODO: can be make simpler as only last block might need scatter copy.
  if (nx % 4 != 0) {
    // Possible partial scatter.
    #pragma _NEC shortloop
    for (blk_idx = bmin; blk_idx < bmax; blk_idx++) {
      int64_t b = blk_idx;
      int64_t x;
      x = 4 * b;
      if (x + 4 > nx) partial_scatter++;
    }
  }

  if (!partial_scatter) {
#ifdef _FTRACE
    ftrace_region_begin("DECOMP - SCATTER");
#endif
    for (blk_idx = bmin; blk_idx < bmax; blk_idx++) {
      uint64_t b = blk_idx; 
      uint64_t x;
      x = 4 * b;
      Scalar *p = &data[x * sx]; 
  
      /* scatter block to strided array */
      int ix;
      #pragma _NEC unroll_completely
      for (ix = 0; ix < 4; ix++, p += sx) {
        *p = fblock[ix * vl + (blk_idx - bmin)];
      }
    }
#ifdef _FTRACE
    ftrace_region_end("DECOMP - SCATTER");
#endif
  } else {
#ifdef _FTRACE
    ftrace_region_begin("DECOMP - PARTIAL SCATTER");
#endif
    #pragma _NEC ivdep
    #pragma _NEC shortloop
    for (blk_idx = bmin; blk_idx < bmax; blk_idx++) {
      uint64_t b = blk_idx;
      uint64_t x;
      x = 4 * b;

      int64_t ix;
      #pragma _NEC unroll_completely
      for (ix = 0; ix < 4; ix++) {
        if (x+ix < nx) {
          data[(x+ix) * sx] = fblock[ix * vl + (blk_idx - bmin)];
        }
      }
    }
#ifdef _FTRACE
    ftrace_region_end("DECOMP - PARTIAL SCATTER");
#endif
  }
}
