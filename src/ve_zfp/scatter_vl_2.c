static void scatter_vl_2(const zfp_field* field, Scalar *fblock, uint64_t bmin, uint64_t bmax)
{
  const int vl = (bmax - bmin);

  Scalar *data = (Scalar *) field->data;

  const uint64_t nx = field->nx;
  const uint64_t ny = field->ny;
  const uint64_t bx = (nx + 3) / 4;
  const uint64_t by = (ny + 3) / 4;
  const uint64_t sx = field->sx ? field->sx : 1;
  const uint64_t sy = field->sy ? field->sy : nx;
  uint64_t blk_idx;
  int partial_scatter = 0;

  if (nx % 4 != 0 || ny % 4 != 0) {
    // Possible partial scatter.
    #pragma _NEC shortloop
    for (blk_idx = bmin; blk_idx < bmax; blk_idx++) {
      int64_t b = blk_idx;
      int64_t x, y;
      x = 4 * (b % bx);
      b /= bx;
      y = 4 * b;
      if (x + 4 > nx || y + 4 > ny) partial_scatter++;
    }
  }

  if (!partial_scatter) {
#ifdef _FTRACE
    ftrace_region_begin("DECOMP - SCATTER");
#endif
    for (blk_idx = bmin; blk_idx < bmax; blk_idx++) {
      uint64_t b = blk_idx; 
      uint64_t x, y;
      x = 4 * (b % bx);
      b /= bx;
      y = 4 * b;
      Scalar *p = &data[x * sx + y * sy]; 
  
      /* scatter block to strided array */
      int ix, iy;
      #pragma _NEC unroll_completely
      for (iy = 0; iy < 4; iy++, p += sy - 4 * sx) {
        #pragma _NEC unroll_completely
        for (ix = 0; ix < 4; ix++, p += sx) {
          *p = fblock[(4 * iy + ix) * vl + (blk_idx - bmin)];
        }
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
      uint64_t x, y;
      x = 4 * (b % bx);
      b /= bx;
      y = 4 * b;

      int64_t ix, iy;
      #pragma _NEC unroll_completely
      for (iy = 0; iy < 4; iy++) {
        #pragma _NEC unroll_completely
        for (ix = 0; ix < 4; ix++) {
          if (x+ix < nx && y+iy < ny) {
            data[(x+ix) * sx + (y+iy) * sy] = fblock[(4 * iy + ix) * vl + (blk_idx - bmin)];
          }
        }
      }
    }
#ifdef _FTRACE
    ftrace_region_end("DECOMP - PARTIAL SCATTER");
#endif
  }
}
