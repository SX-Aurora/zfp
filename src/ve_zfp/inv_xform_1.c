/* ------------------------------------------------------------------------------------------------ */
/* inverse decorrelating 1D transform */
static void inv_xform_1(Int* iblock, int blk_idx, const int vl)
{
  Int xx, yy, zz, ww;

  /* transform along x */
  xx = iblock[blk_idx + 0 * vl];
  yy = iblock[blk_idx + 1 * vl];
  zz = iblock[blk_idx + 2 * vl];
  ww = iblock[blk_idx + 3 * vl];
    
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

  iblock[blk_idx + 0 * vl] = xx;
  iblock[blk_idx + 1 * vl] = yy;
  iblock[blk_idx + 2 * vl] = zz;
  iblock[blk_idx + 3 * vl] = ww;
}
