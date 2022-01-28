/* --------------------------------------------------------------------------------------------------------------------------------------------- */
/* forward decorrelating 1D transform */
static void fwd_xform_1(Int* iblock, int blk_idx, const int vl)
{
  Int xx, yy, zz, ww;

  /* transform along x */
  xx = iblock[blk_idx + 0 * vl];
  yy = iblock[blk_idx + 1 * vl];
  zz = iblock[blk_idx + 2 * vl];
  ww = iblock[blk_idx + 3 * vl];

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

  iblock[blk_idx + 0 * vl] = xx;
  iblock[blk_idx + 1 * vl] = yy;
  iblock[blk_idx + 2 * vl] = zz;
  iblock[blk_idx + 3 * vl] = ww;

}
