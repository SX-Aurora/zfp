/* --------------------------------------------------------------------------------------------------------------------------------------------- */
/* forward decorrelating 3D transform */
static void fwd_xform_3(Int* iblock, int blk_idx, const int vl)
{
  int x, y, z;
  Int xx, yy, zz, ww;

  /* transform along x */
  #pragma _NEC unroll_completely
  for (z = 0; z < 4; z++) {
    #pragma _NEC unroll_completely
    for (y = 0; y < 4; y++) {
      xx = iblock[blk_idx + (0 + 4 * y + 16 * z) * vl];
      yy = iblock[blk_idx + (1 + 4 * y + 16 * z) * vl];
      zz = iblock[blk_idx + (2 + 4 * y + 16 * z) * vl];
      ww = iblock[blk_idx + (3 + 4 * y + 16 * z) * vl];

      /*
      ** non-orthogonal transform
      **        ( 4  4  4  4) (xx)
      ** 1/16 * ( 5  1 -1 -5) (yy)
      **        (-4  4  4 -4) (zz)
      **        (-2  6 -6  2) (ww)
      */
      xx += ww;
      xx >>= 1;
      ww -= xx;
      zz += yy;
      zz >>= 1;
      yy -= zz;
      xx += zz;
      xx >>= 1;
      zz -= xx;
      ww += yy;
      ww >>= 1;
      yy -= ww;
      ww += yy >> 1;
      yy -= ww >> 1;

      iblock[blk_idx + (0 + 4 * y + 16 * z) * vl] = xx;
      iblock[blk_idx + (1 + 4 * y + 16 * z) * vl] = yy;
      iblock[blk_idx + (2 + 4 * y + 16 * z) * vl] = zz;
      iblock[blk_idx + (3 + 4 * y + 16 * z) * vl] = ww;
    }
  }

  /* transform along y */
  #pragma _NEC unroll_completely
  for (x = 0; x < 4; x++) {
    #pragma _NEC unroll_completely
    for (z = 0; z < 4; z++) {
      xx = iblock[blk_idx + (x + 4 * 0 + 16 * z) * vl];
      yy = iblock[blk_idx + (x + 4 * 1 + 16 * z) * vl];
      zz = iblock[blk_idx + (x + 4 * 2 + 16 * z) * vl];
      ww = iblock[blk_idx + (x + 4 * 3 + 16 * z) * vl];

      /*
      ** non-orthogonal transform
      **        ( 4  4  4  4) (xx)
      ** 1/16 * ( 5  1 -1 -5) (yy)
      **        (-4  4  4 -4) (zz)
      **        (-2  6 -6  2) (ww)
      */
      xx += ww;
      xx >>= 1;
      ww -= xx;
      zz += yy;
      zz >>= 1;
      yy -= zz;
      xx += zz;
      xx >>= 1;
      zz -= xx;
      ww += yy;
      ww >>= 1;
      yy -= ww;
      ww += yy >> 1;
      yy -= ww >> 1;

      iblock[blk_idx + (x + 4 * 0 + 16 * z) * vl] = xx;
      iblock[blk_idx + (x + 4 * 1 + 16 * z) * vl] = yy;
      iblock[blk_idx + (x + 4 * 2 + 16 * z) * vl] = zz;
      iblock[blk_idx + (x + 4 * 3 + 16 * z) * vl] = ww;
    }
  }

  /* transform along z */
  #pragma _NEC unroll_completely
  for (y = 0; y < 4; y++) {
    #pragma _NEC unroll_completely
    for (x = 0; x < 4; x++) {
      xx = iblock[blk_idx + (x + 4 * y + 16 * 0) * vl];
      yy = iblock[blk_idx + (x + 4 * y + 16 * 1) * vl];
      zz = iblock[blk_idx + (x + 4 * y + 16 * 2) * vl];
      ww = iblock[blk_idx + (x + 4 * y + 16 * 3) * vl];

      /*
      ** non-orthogonal transform
      **        ( 4  4  4  4) (x)
      ** 1/16 * ( 5  1 -1 -5) (y)
      **        (-4  4  4 -4) (z)
      **        (-2  6 -6  2) (w)
      */
      xx += ww;
      xx >>= 1;
      ww -= xx;
      zz += yy;
      zz >>= 1;
      yy -= zz;
      xx += zz;
      xx >>= 1;
      zz -= xx;
      ww += yy;
      ww >>= 1;
      yy -= ww;
      ww += yy >> 1;
      yy -= ww >> 1;

      iblock[blk_idx + (x + 4 * y + 16 * 0) * vl] = xx;
      iblock[blk_idx + (x + 4 * y + 16 * 1) * vl] = yy;
      iblock[blk_idx + (x + 4 * y + 16 * 2) * vl] = zz;
      iblock[blk_idx + (x + 4 * y + 16 * 3) * vl] = ww;
    }
  }
}
