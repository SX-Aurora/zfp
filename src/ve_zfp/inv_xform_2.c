/* ------------------------------------------------------------------------------------------------ */
/* inverse decorrelating 2D transform */
static void inv_xform_2(Int* iblock, int blk_idx, const int vl)
{
  int x, y;
  Int xx, yy, zz, ww;

  /* transform along y */
  #pragma _NEC unroll_completely
  for (x = 0; x < 4; x++) {
    xx = iblock[blk_idx + (x + 4 * 0) * vl];
    yy = iblock[blk_idx + (x + 4 * 1) * vl];
    zz = iblock[blk_idx + (x + 4 * 2) * vl];
    ww = iblock[blk_idx + (x + 4 * 3) * vl];
    
    /*
    ** non-orthogonal transform
    **       ( 4  6 -4 -1) (xx)
    ** 1/4 * ( 4  2  4  5) (yy)
    **       ( 4 -2  4 -5) (zz)
    **       ( 4 -6 -4  1) (ww)
    */
    yy += ww >> 1; ww -= yy >> 1;
    yy += ww; ww <<= 1; ww -= yy;
    zz += xx; xx <<= 1; xx -= zz;
    yy += zz; zz <<= 1; zz -= yy;
    ww += xx; xx <<= 1; xx -= ww;
    
    iblock[blk_idx + (x + 4 * 0) * vl] = xx;
    iblock[blk_idx + (x + 4 * 1) * vl] = yy;
    iblock[blk_idx + (x + 4 * 2) * vl] = zz;
    iblock[blk_idx + (x + 4 * 3) * vl] = ww;
  }

  /* transform along x */
  #pragma _NEC unroll_completely
  for (y = 0; y < 4; y++) {
    xx = iblock[blk_idx + (0 + 4 * y) * vl];
    yy = iblock[blk_idx + (1 + 4 * y) * vl];
    zz = iblock[blk_idx + (2 + 4 * y) * vl];
    ww = iblock[blk_idx + (3 + 4 * y) * vl];
    
    /*
    ** non-orthogonal transform
    **       ( 4  6 -4 -1) (xx)
    ** 1/4 * ( 4  2  4  5) (yy)
    **       ( 4 -2  4 -5) (zz)
    **       ( 4 -6 -4  1) (ww)
    */
    yy += ww >> 1; ww -= yy >> 1;
    yy += ww; ww <<= 1; ww -= yy;
    zz += xx; xx <<= 1; xx -= zz;
    yy += zz; zz <<= 1; zz -= yy;
    ww += xx; xx <<= 1; xx -= ww;
    
    iblock[blk_idx + (0 + 4 * y) * vl] = xx;
    iblock[blk_idx + (1 + 4 * y) * vl] = yy;
    iblock[blk_idx + (2 + 4 * y) * vl] = zz;
    iblock[blk_idx + (3 + 4 * y) * vl] = ww;
  }
}
