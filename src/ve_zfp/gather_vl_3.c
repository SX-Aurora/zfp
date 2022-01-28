// n = 0, 1, 2, 3, 4
// n==4, nothing to copy
// n==3,                       p3=p0
// n==2,              p2=p1    p3=p0
// n==1,      p1=p0   p2=p1=p0 p3=p0
// n==0, p0=0 p1=p0=0 p2=p1=p0 p3=p0=0
//
//     nx-4   nx-3   nx-2   nx-1   nx+0   nx+1   nx+2   nx+3   nx+4
//     --+------+------+------+------+------+------+------+-----+-----
//
//    ix 0      1      2      3 .............................. nothing to do
//    
//    ix        0      1      2      3
//   xx-nx     -3     -2     -1      0
//    ofs                           -3
//    n=x+3-nx => 0
//
//    ix               0      1      2      3
//   xx-nx            -2     -1      0     +1
//    ofs                           -1     -3
//    n=x+3-nx => 1
//
//    ix                      0      1      2      3
//   xx-nx                   -1      0     +1     +2
//    ofs                           -1     -2     -3
//    n=x+3-nx => 2
//
//    ix                             0      1      2      3 .. not possible
//
// tab_ofs[3  x+3-nx][4  ix] = { {0,0,0,-3}, {0,0,-1,-3}, {0,-1,-2,-3} }
//
static const int tab_ofs[3][4] = { {0,0,0,-3}, {0,0,-1,-3}, {0,-1,-2,-3} };

static void gather_vl_3(const zfp_field* field, Scalar *fblock_in, uint64_t bmin, uint64_t bmax)
{
  const int vl = (bmax - bmin);

  Scalar *data = (Scalar *) field->data;

  int ix, iy, iz;
  int64_t xx[vl * 4];
  int64_t yy[vl * 4];
  int64_t zz[vl * 4];

  uint64_t blk_idx;
  const int64_t nx = field->nx;
  const int64_t ny = field->ny;
  const int64_t nz = field->nz;
  const int64_t sx = field->sx ? field->sx : 1;
  const int64_t sy = field->sy ? field->sy : nx;
  const int64_t sz = field->sz ? field->sz : nx * ny;

  const int64_t bx = (nx + 3) / 4;
  const int64_t by = (ny + 3) / 4;
  const int64_t bz = (nz + 3) / 4;
  int partial_gather = 0;

  if (nx % 4 != 0 || ny % 4 != 0 || nz % 4 != 0) {
    // Possible partial gather.
    for (blk_idx = bmin; blk_idx < bmax; blk_idx++) {
      int64_t b = blk_idx;
      int64_t x, y, z;
      x = 4 * (b % bx);
      b /= bx;
      y = 4 * (b % by);
      b /= by;
      z = 4 * b;
      if (x + 4 > nx || y + 4 > ny || z + 4 > nz) partial_gather++;
    }
  }

  if (!partial_gather) {
#ifdef _FTRACE
    ftrace_region_begin("COMP - GATHER");
#endif
    #pragma _NEC ivdep
    #pragma _NEC shortloop
    for (blk_idx = bmin; blk_idx < bmax; blk_idx++) {
      int64_t b = blk_idx;
      int64_t x, y, z;
      x = 4 * (b % bx);
      b /= bx;
      y = 4 * (b % by);
      b /= by;
      z = 4 * b;
      const Scalar *p = &data[x * sx + y * sy + z * sz];

      #pragma _NEC unroll_completely
      for (iz = 0; iz < 4; iz++, p += sz - 4 * sy) {
        #pragma _NEC unroll_completely
        for (iy = 0; iy < 4; iy++, p += sy - 4 * sx) {
          #pragma _NEC unroll_completely
          for (ix = 0; ix < 4; ix++, p += sx) {
            fblock_in[(16 * iz + 4 * iy + ix) * vl + (blk_idx - bmin)] = *p;
          }
        }
      }
    }
#ifdef _FTRACE
    ftrace_region_end("COMP - GATHER");
#endif
  } else {
#ifdef _FTRACE
    ftrace_region_begin("COMP - PARTIAL GATHER1");
#endif
    #pragma _NEC ivdep
    #pragma _NEC shortloop
    for (blk_idx = bmin; blk_idx < bmax; blk_idx++) {
      int64_t b = blk_idx;
      int64_t x, y, z;
      int i;
      x = 4 * (b % bx);
      b /= bx;
      y = 4 * (b % by);
      b /= by;
      z = 4 * b;
      #pragma _NEC unroll_completely
      for (i = 0; i < 4; i++) {
        int o = blk_idx - bmin + i * vl;
        xx[o] = x + i;
        if (xx[o] >= nx) xx[o] += tab_ofs[x+3-nx][i];
        yy[o] = y + i;
        if (yy[o] >= ny) yy[o] += tab_ofs[y+3-ny][i];
        zz[o] = z + i;
        if (zz[o] >= nz) zz[o] += tab_ofs[z+3-nz][i];
      }
    }
#ifdef _FTRACE
    ftrace_region_end("COMP - PARTIAL GATHER1");
#endif

#ifdef _FTRACE
    ftrace_region_begin("COMP - PARTIAL GATHER2");
#endif
    #pragma _NEC ivdep
    #pragma _NEC shortloop
    for (blk_idx = bmin; blk_idx < bmax; blk_idx++) {
      int64_t x, y, z;
      #pragma _NEC unroll_completely
      for (iz = 0; iz < 4; iz++) {
        z = zz[blk_idx - bmin + iz * vl];
        #pragma _NEC unroll_completely
        for (iy = 0; iy < 4; iy++) {
          y = yy[blk_idx - bmin + iy * vl];
          #pragma _NEC unroll_completely
          for (ix = 0; ix < 4; ix++) {
            x = xx[blk_idx - bmin + ix * vl];
            fblock_in[(16 * iz + 4 * iy + ix) * vl + (blk_idx - bmin)] = data[x * sx + y * sy + z * sz];
          }
        }
      }
    }
#ifdef _FTRACE
    ftrace_region_end("COMP - PARTIAL GATHER2");
#endif
  }
}
