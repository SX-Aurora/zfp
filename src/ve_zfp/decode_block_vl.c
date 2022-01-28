static void decode_block_vl(zfp_stream* zfp, bitstream *lstream, Scalar* fblock, uint64_t bmin, uint64_t bmax)
{
  const int vl = (bmax - bmin);

  Int iblock[BLOCK_SIZE * VL] __attribute__((aligned(8)));
  UInt ublock[BLOCK_SIZE * VL] __attribute__((aligned(8)));

  int blk_idx;
  int i;
  int flag_planes;
  int flag_rle;
  int flag_x;
  uint32_t intprec, k;

  uint32_t e[VL];
  #pragma _NEC vreg(e)
  int emax[VL];
  #pragma _NEC vreg(emax)
  int maxprec[VL];
  #pragma _NEC vreg(maxprec)
  int lbits[VL];
  #pragma _NEC vreg(lbits)
  uint32_t kmin[VL];
  #pragma _NEC vreg(kmin)
  Scalar tmp_ldexp[VL];
  #pragma _NEC vreg(tmp_ldexp)
  uint32_t m[VL];
  #pragma _NEC vreg(m)
  uint32_t n[VL];
  #pragma _NEC vreg(n)
  uint64_t x[VL];
  #pragma _NEC vreg(x)
  int state[VL];
  #pragma _NEC vreg(state)

  uint32_t bs_bits[VL];
  #pragma _NEC vreg(bs_bits)
  uint64_t bs_buffer[VL];
  #pragma _NEC vreg(bs_buffer)
  uint64_t *bs_ptr[VL] __attribute__((aligned(8)));

#ifdef _FTRACE
  ftrace_region_begin("DECOMP - DECODE EXP");
#endif
  /* Initialize bitstream for each block */
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    bs_bits[blk_idx] = 0;
    bs_buffer[blk_idx] = 0;
    bs_ptr[blk_idx] = (uint64_t *) ((uint8_t*) lstream->begin + (blk_idx) * (zfp->maxbits / CHAR_BIT));
  }

  /* test if block has nonzero values */
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    e[blk_idx] = lstream_read_bit(bs_bits, bs_buffer, bs_ptr, blk_idx);
  }

  /* decode common exponent */
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    if (e[blk_idx]) {
      emax[blk_idx] = (int)lstream_read_bits(EBITS, bs_bits, bs_buffer, bs_ptr, blk_idx) - EBIAS;
    }
  }

  /*
   * maximum number of bit planes to decode
   *
   * Apply computation for all blocks even empty one:
   * . maxprec not used for empty block.
   * . prevent to build vector mask
   */
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    maxprec[blk_idx] = MIN(zfp->maxprec, (uint32_t)MAX(0, emax[blk_idx] - zfp->minexp + 2 * (DIMS + 1)));
  }

  /* initialize data array to all zeros */
  for (i = 0; i < BLOCK_SIZE * vl; i++) {
    ublock[i] = 0.0;
  }
  
  /* decode integer coefficients */
  intprec = CHAR_BIT * (uint32_t)sizeof(UInt);
  k = intprec;
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    if (e[blk_idx]) {
      kmin[blk_idx] = intprec > maxprec[blk_idx] ? intprec - maxprec[blk_idx] : 0;
      lbits[blk_idx] = zfp->maxbits - (1 + EBITS);
      n[blk_idx] = 0;
      if (lbits[blk_idx] && k > kmin[blk_idx]) flag_planes++;
    }
  }
#ifdef _FTRACE
  ftrace_region_end("DECOMP - DECODE EXP");
#endif

#ifdef _FTRACE
  ftrace_region_begin("DECOMP - BITSTREAM");
#endif
#if DIMS!=4
  /* decode one bit plane at a time from MSB to LSB */
  while (flag_planes) {
    flag_planes = 0;
    k--;

    /* decode first n bits of bit plane #k */
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      m[blk_idx] = MIN(n[blk_idx], lbits[blk_idx]);
      lbits[blk_idx] -= m[blk_idx];
      if (e[blk_idx]) {
        x[blk_idx] = lstream_read_bits(m[blk_idx], bs_bits, bs_buffer, bs_ptr, blk_idx);
      }
    }

    /* unary run-length decode remainder of bit plane */
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      state[blk_idx] = !!e[blk_idx];
      flag_rle += state[blk_idx];
    }

    while (flag_rle) {
      flag_rle = 0;
      #pragma _NEC shortloop
      for (blk_idx = 0; blk_idx < vl; blk_idx++) {
        if (n[blk_idx] == BLOCK_SIZE || lbits[blk_idx] == 0) state[blk_idx] = 0;
        if (state[blk_idx] == 1) {
          // If bit == 0, state = 0, remaining bits are 0, stop decoding for this plane.
          // If bit == 1, state = 1, next bits are supposed to be 0, decode them until '1' is reached and output 1.
          state[blk_idx] = lstream_read_bit(bs_bits, bs_buffer, bs_ptr, blk_idx);
          lbits[blk_idx]--;
        }
        if (state[blk_idx] >= 1) {
          if (n[blk_idx] < BLOCK_SIZE - 1 && lbits[blk_idx]) {
            // If bit == 0, state = 1 + 1 = 2, continue to read '0' bit.
            // If bit == 1, state = 0 + 1 = 1, stop reading bit and output 1.
            state[blk_idx] = 1 + !lstream_read_bit(bs_bits, bs_buffer, bs_ptr, blk_idx);
            lbits[blk_idx]--;
          }
          n[blk_idx]++;
        }
        if (state[blk_idx] == 2 && lbits[blk_idx] == 0) {
          // We are reading '0' but reach end of available bits, we will have to output last bit '1' at the correct place.
          n[blk_idx]++;
        }
        if (state[blk_idx] == 1 || (state[blk_idx] == 2 && (n[blk_idx] == BLOCK_SIZE || lbits[blk_idx] == 0))) {
          // state == 1, we were reading '0' and read '1', so output 1.
          // state == 2, we are reading '0' and reach end of block or available bits, so output 1.
          x[blk_idx] += (uint64_t)1 << (n[blk_idx]-1);
        }
        flag_rle += state[blk_idx];
      } // for(blk_idx)
    } // while(flag_rle)

    /* deposit bit plane from x */
    flag_x = 1;
    i = 0;
    while (flag_x) {
      flag_x = 0;
      #pragma _NEC shortloop
      for (blk_idx = 0; blk_idx < vl; blk_idx++) {
        ublock[blk_idx + i * vl] += (UInt)(x[blk_idx] & 1u) << k;
      }
      i++;
      #pragma _NEC shortloop
      for (blk_idx = 0; blk_idx < vl; blk_idx++) {
        if (e[blk_idx]) {
          x[blk_idx] = x[blk_idx] >> 1;
          if (x[blk_idx]) flag_x++;
        }
      }
    }

    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      if (e[blk_idx]) {
        if (lbits[blk_idx] && k > kmin[blk_idx]) flag_planes++;
      }
    }
  } // while (flag_planes)
#else
  /* decode one bit plane at a time from MSB to LSB */
  while (flag_planes) {
    flag_planes = 0;
    k--;

    /* decode first n bits of bit plane #k */
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      if (e[blk_idx]) {
        m[blk_idx] = MIN(n[blk_idx], lbits[blk_idx]);
        lbits[blk_idx] -= m[blk_idx];
      }
    }

    uint32_t max_m = 0;
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      if (e[blk_idx] && m[blk_idx] > max_m) max_m = m[blk_idx];
    }

    for (i = 0; i < max_m; i++) {
      #pragma _NEC shortloop
      for (blk_idx = 0; blk_idx < vl; blk_idx++) {
        if (e[blk_idx] && i < m[blk_idx]) {
          if (lstream_read_bit(bs_bits, bs_buffer, bs_ptr, blk_idx))
            ublock[blk_idx + i * vl] += (UInt)1 << k;
        }
      }
    }

    /* unary run-length decode remainder of bit plane */
    flag_rle = 0;
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      state[blk_idx] = !!e[blk_idx];
      flag_rle += state[blk_idx];
    }

    while (flag_rle) {
      flag_rle = 0;
      #pragma _NEC ivdep
      #pragma _NEC shortloop
      for (blk_idx = 0; blk_idx < vl; blk_idx++) {
        if (n[blk_idx] == BLOCK_SIZE || lbits[blk_idx] == 0) state[blk_idx] = 0;
        if (state[blk_idx] == 1) {
          // If bit == 0, state = 0, remaining bits are 0, stop decoding for this plane.
          // If bit == 1, state = 1, next bits are supposed to be 0, decode them until '1' is reached and output 1.
          state[blk_idx] = lstream_read_bit(bs_bits, bs_buffer, bs_ptr, blk_idx);
          lbits[blk_idx]--;
        }
        if (state[blk_idx] >= 1) {
          if (n[blk_idx] < BLOCK_SIZE - 1 && lbits[blk_idx]) {
            // If bit == 0, state = 1 + 1 = 2, continue to read '0' bit.
            // If bit == 1, state = 0 + 1 = 1, stop reading bit and output 1.
            state[blk_idx] = 1 + !lstream_read_bit(bs_bits, bs_buffer, bs_ptr, blk_idx);
            lbits[blk_idx]--;
          }
          n[blk_idx]++;
        }
        if (state[blk_idx] == 2 && lbits[blk_idx] == 0) {
          // We are reading '0' but reach end of available bits, we will have to output last bit '1' at the correct place.
          n[blk_idx]++;
        }
        if (state[blk_idx] == 1 || (state[blk_idx] == 2 && (n[blk_idx] == BLOCK_SIZE || lbits[blk_idx] == 0))) {
          // state == 1, we were reading '0' and read '1', so output 1.
          // state == 2, we are reading '0' and reach end of block or available bits, so output 1.
          ublock[blk_idx + (n[blk_idx]-1) * vl] += (UInt)1 << k;
        }
        flag_rle += state[blk_idx];
      } // for(blk_idx)
    } // while (flag_rle)

    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      if (e[blk_idx]) {
        if (lbits[blk_idx] && k > kmin[blk_idx]) flag_planes++;
      }
    }

  } // while(flag_planes)
#endif
#ifdef _FTRACE
  ftrace_region_end("DECOMP - BITSTREAM");
#endif
  /*
   * reorder unsigned coefficients and convert to signed integer
   *
   * Apply computation for all blocks even empty one:
   * . ublock was set to zero for empty block
   * . iblock will be equal to zero for empty block
   * . prevent to build vector mask
   */
#ifdef _FTRACE
  ftrace_region_begin("DECOMP - TRANSFORM");
#endif
  for (i = 0; i < BLOCK_SIZE; i++) { 
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      iblock[blk_idx + PERM[i] * vl] = (Int)(((ublock[blk_idx + i * vl]) ^ NBMASK) - NBMASK);
    }
  }

  /*
   * perform decorrelating transform
   *
   * Apply computation for all blocks even empty one:
   * . iblock was set to zero for empty block
   * . prevent to build vector mask
   */
  #pragma _NEC ivdep
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    INV_XFORM(iblock, blk_idx, vl);
  }

  /* perform inverse block-floating-point transform */
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    tmp_ldexp[blk_idx] = 0.0;
  }  
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    if (e[blk_idx]) {
      /* map integer x relative to exponent e to floating-point number */
      // ldexpf is not supported by NEC compilers.
      // ldexpf(1.0f, emax[blk_idx] - (CHAR_BIT * (int)sizeof(float) - 2));
      tmp_ldexp[blk_idx] = POW(2.0, (Scalar)(emax[blk_idx] - (CHAR_BIT * (int)sizeof(Scalar) - 2)));
    }
  }

  /*
   * Apply computation for all blocks even empty one:
   * . iblock == 0 for empty block, this set fblock to 0.0f
   * . prevent to build vector mask
   * . prevent to deal with empty block.
   */ 
  for (i = 0; i < BLOCK_SIZE; i++) {
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      /* compute power-of-two scale factor s */
      /* compute p-bit float x = s*y where |y| <= 2^(p-2) - 1 */
      fblock[blk_idx + i * vl] = (Scalar)(tmp_ldexp[blk_idx] * iblock[blk_idx + i * vl]);
    }
  }
#ifdef _FTRACE
  ftrace_region_end("DECOMP - TRANSFORM");
#endif
}
