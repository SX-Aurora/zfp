static void encode_block_vl(zfp_stream* zfp, bitstream *lstream, const Scalar* fblock, int64_t bmin, int64_t bmax)
{
  const int vl = (bmax - bmin);

  Int  iblock[BLOCK_SIZE * VL] __attribute__((aligned(8)));
  UInt ublock[BLOCK_SIZE * VL] __attribute__((aligned(8)));

  Scalar max[VL];
  #pragma _NEC vreg(max)
  Scalar tmp_ldexp[VL];
  #pragma _NEC vreg(tmp_ldexp)
  int emax[VL];
  #pragma _NEC vreg(emax)
  int maxprec[VL];
  #pragma _NEC vreg(maxprec)
  uint32_t e[VL];
  #pragma _NEC vreg(e)
  uint32_t bs_bits[VL];
  #pragma _NEC vreg(bs_bits)
  uint64_t bs_buffer[VL];
  #pragma _NEC vreg(bs_buffer)
  uint64_t *bs_ptr[VL] __attribute__((aligned(8)));
  int lbits[VL];
  #pragma _NEC vreg(lbits)
  uint64_t x[VL];
  #pragma _NEC vreg(x)
  uint32_t m[VL];
  #pragma _NEC vreg(m)
  uint32_t n[VL];
  #pragma _NEC vreg(n)
  int state[VL];
  #pragma _NEC vreg(state)
  uint32_t kmin[VL];
  #pragma _NEC vreg(kmin)
  uint32_t tmp_bits[VL];
  #pragma _NEC vreg(tmp_bits)
  uint64_t tmp_buffer[VL];
  #pragma _NEC vreg(tmp_buffer)
#if DIMS==4
  int c[VL];
  #pragma _NEC vreg(c)
#endif

  int flag_planes;
  int flag_rle;
  int blk_idx;
  int i;
  uint32_t intprec;
  uint32_t k;

#ifdef _FTRACE
  ftrace_region_begin("COMP - TRANSFORM");
#endif
  /* Compute maximum exponent */
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    max[blk_idx] = 0.0;
  }

  for (i = 0; i < BLOCK_SIZE; i++) {
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      Scalar f = FABS(fblock[blk_idx + i * vl]);
      if (max[blk_idx] < f)
        max[blk_idx] = f;
    }
  }

  /* normalized maximum exponent for x >= 0 */
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    if (max[blk_idx] > 0) {
      int e;
      // frexpf is not supported by NEC compilers.
      // frexpf(max[blk_idx], &e);
      e = (int)(1 + FLOOR(LOG2(max[blk_idx])));
      emax[blk_idx] = MAX(e, 1 - EBIAS);
    } else {
      emax[blk_idx] = -EBIAS;
    }
  }

  /* maximum number of bit planes to encode */
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    maxprec[blk_idx] = MIN(zfp->maxprec, (uint32_t)MAX(0, emax[blk_idx] - zfp->minexp + 2 * (DIMS + 1)));
  }

  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    e[blk_idx] = maxprec[blk_idx] ? emax[blk_idx] + EBIAS : 0;
  }

  /* Initialize bitstream for each block */
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    bs_bits[blk_idx] = 0;
    bs_buffer[blk_idx] = 0;
    bs_ptr[blk_idx] = (uint64_t *) ((uint8_t*) lstream->begin + (blk_idx) * (zfp->maxbits / CHAR_BIT));
  }

#if DIMS!=4
  /*
   * Encode common exponent; LSB indicates that exponent is nonzero
   *
   * Note: . With fix bit rate we always write exponent at beginning of the
   *         local buffer (bits per block = 64 values x N bits = N x 64 bits).
   *       . Max exponent size is EBITS = 8 bits (don't update stream,
   *         update just local buffer).
   */
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    if (e[blk_idx]) {
      bs_buffer[blk_idx] = 2 * e[blk_idx] + 1;
      bs_bits[blk_idx] = 1 + EBITS;
      lbits[blk_idx] = zfp->maxbits - (1 + EBITS);
    } else {
      lbits[blk_idx] = zfp->maxbits;
    }
  }
#else
  /*
   * Encode common exponent; LSB indicates that exponent is nonzero
   *
   */
  #pragma _NEC ivdep
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    if (e[blk_idx]) {
      lbits[blk_idx] = zfp->maxbits - (1 + EBITS);
      lstream_write_bits(2 * e[blk_idx] + 1, 1 + EBITS, bs_bits, bs_buffer, bs_ptr, blk_idx);
    } else {
      lbits[blk_idx] = zfp->maxbits;
    }
  }
#endif

  /* encode integer block */
  /* map floating-point number x to integer relative to exponent e */
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    if (e[blk_idx]) {
      // ldexpf is not supported by NEC compilers.
      // ldexpf(1.0f, (CHAR_BIT * (int)sizeof(float) - 2) - emax[blk_idx]);
      tmp_ldexp[blk_idx] = POW(2.0, (Scalar)((CHAR_BIT * (int32_t)sizeof(Scalar) - 2) - emax[blk_idx]));
    }
  }

  /*
   * Perform forward block-floating-point transform
   *
   * Apply computation for all blocks even empty one:
   * . iblock not used for empty block
   * . prevent to build vector mask
   */
  for (i = 0; i < BLOCK_SIZE; i++) {
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      iblock[blk_idx + i * vl] = (Int)(tmp_ldexp[blk_idx] * fblock[blk_idx + i * vl]);
    }
  }

  /*
   * Perform decorrelating transform
   *
   * Apply computation for all blocks even empty one:
   * . iblock not used for empty block
   * . prevent to build vector mask
   */
  #pragma _NEC ivdep
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    FWD_XFORM(iblock, blk_idx, vl);
  }

  /*
   * Reorder signed coefficients and convert to unsigned integer
   *
   * Apply computation for all blocks even empty one:
   * . ublock not used for empty block
   * . prevent to build vector mask
   */
  for (i = 0; i < BLOCK_SIZE; i++) {
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      /* map two's complement signed integer to negabinary unsigned integer */
      ublock[blk_idx + i * vl] = ((UInt)iblock[blk_idx + PERM[i] * vl] + NBMASK) ^ NBMASK;
    }
  }

  /* encode integer coefficients */
  intprec = CHAR_BIT * (uint32_t)sizeof(UInt);
  k = intprec;
  flag_planes = 0;
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    if (e[blk_idx]) {
      kmin[blk_idx] = intprec > maxprec[blk_idx] ? intprec - maxprec[blk_idx] : 0;
      n[blk_idx] = 0;
      if (lbits[blk_idx] && k > kmin[blk_idx]) flag_planes++;
    }
  }
#ifdef _FTRACE
  ftrace_region_end("COMP - TRANSFORM");
#endif

#ifdef _FTRACE
  ftrace_region_begin("COMP - BITSTREAM");
#endif
#if DIMS!=4
  /* Bitstream encoding part for 1, 2 and 3D arrays */
  /* encode one bit plane at a time from MSB to LSB */
  while (flag_planes) {
    flag_planes = 0;
    k--;

    /* step 1: extract bit plane #k to x */
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      x[blk_idx] = 0;
      tmp_buffer[blk_idx] = 0;
      tmp_bits[blk_idx] = 0;
    }
    for (i = 0; i < BLOCK_SIZE; i++) {
      /* 
       * Apply computation for all blocks even empty one:
       * . ublock not used for empty block
       * . prevent to build vector mask
       */
      #pragma _NEC shortloop
      for (blk_idx = 0; blk_idx < vl; blk_idx++) {
        x[blk_idx] += (uint64_t)((ublock[blk_idx + i * vl] >> k) & 1u) << i;
      }
    }

    /* step 2: encode first n bits of bit plane */
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      m[blk_idx] = MIN(n[blk_idx], lbits[blk_idx]);
      lbits[blk_idx] -= m[blk_idx];
      if (e[blk_idx]) {
        x[blk_idx] = lstream_write_bits(x[blk_idx], m[blk_idx], bs_bits, bs_buffer, bs_ptr, blk_idx);
      }
    }

    /* step 3: unary run-length encode remainder of bit plane */
    flag_rle = 0;
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      state[blk_idx] = !!e[blk_idx];
      flag_rle += state[blk_idx];
    }

    int cnt = 0;
    while (flag_rle) {
      flag_rle = 0;
      #pragma _NEC ivdep
      #pragma _NEC shortloop
      for (blk_idx = 0; blk_idx < vl; blk_idx++) {
        if (n[blk_idx] == BLOCK_SIZE || lbits[blk_idx] == 0) state[blk_idx] = 0;
        if (state[blk_idx] == 1) {
          // If there are only 0 remaining, output 0 and stop, otherwise output 1 and continue.
          uint32_t bit = !!x[blk_idx];
          tmp_buffer[blk_idx] += (uint64_t)bit << tmp_bits[blk_idx];
          tmp_bits[blk_idx]++;
          state[blk_idx] = bit;
          lbits[blk_idx]--;
        }
        if (state[blk_idx] >= 1) {
          if (n[blk_idx] < BLOCK_SIZE - 1 && lbits[blk_idx]) {
            // Output current bit
            // If bit = 0, state = 1 + 1 = 2, continue without checking remaining bits.
            // If bit = 1, state = 1 + 0, continue and check remaining bits.
            uint32_t bit = x[blk_idx] & 1u;
            tmp_buffer[blk_idx] += (uint64_t)bit << tmp_bits[blk_idx];
            tmp_bits[blk_idx]++;
            state[blk_idx] = 1 + !bit;
            lbits[blk_idx]--;
          }
          x[blk_idx] = x[blk_idx] >> 1;
          n[blk_idx]++;
        }
        flag_rle += state[blk_idx];
      } // for (blk_idx)
      cnt += 2;
      if (cnt == 64) {
        // In the worst case, 64 bits are ready to be output
        #pragma _NEC ivdep
        #pragma _NEC shortloop
        for (blk_idx = 0; blk_idx < vl; blk_idx++) {
          if (tmp_bits[blk_idx]) {
            lstream_write_bits(tmp_buffer[blk_idx], tmp_bits[blk_idx], bs_bits, bs_buffer, bs_ptr, blk_idx);
            tmp_bits[blk_idx] = 0;
            tmp_buffer[blk_idx] = 0;
          }
        }
        cnt = 0;
      }
    } // while (flag_rle)

    // Output bits for this plane #k
    #pragma _NEC ivdep
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      if (tmp_bits[blk_idx]) {
        lstream_write_bits(tmp_buffer[blk_idx], tmp_bits[blk_idx], bs_bits, bs_buffer, bs_ptr, blk_idx);
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
  /* Bitstream encoding part for 4D arrays */
  /* encode one bit plane at a time from MSB to LSB */
  while (flag_planes) {
    flag_planes = 0;
    k--;

    /* step 1: encode first n bits of bit plane #k */
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      tmp_buffer[blk_idx] = 0;
      tmp_bits[blk_idx] = 0;
      if (e[blk_idx]) {
        m[blk_idx] = MIN(n[blk_idx], lbits[blk_idx]);
        lbits[blk_idx] -= m[blk_idx];
      }
    }
    int cnt = 0;
    int max_m = 0;
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      if (e[blk_idx] && m[blk_idx] > max_m) max_m = m[blk_idx];
    } 
    for (i = 0; i < max_m; i++) {
      #pragma _NEC ivdep
      #pragma _NEC shortloop
      for (blk_idx = 0; blk_idx < vl; blk_idx++) {
        if (e[blk_idx] && i < m[blk_idx]) {
          uint32_t bit = (ublock[blk_idx + i * vl] >> k) & 1u;
          tmp_buffer[blk_idx] += (uint64_t)bit << tmp_bits[blk_idx];
          tmp_bits[blk_idx]++;
        }
      }
      cnt++;
      if (cnt == 62) {
        // Flush if we reach 62 bits, so we keep room for rle encoding.
        #pragma _NEC ivdep
        #pragma _NEC shortloop
        for (blk_idx = 0; blk_idx < vl; blk_idx++) {
          if (tmp_bits[blk_idx]) {
            lstream_write_bits(tmp_buffer[blk_idx], tmp_bits[blk_idx], bs_bits, bs_buffer, bs_ptr, blk_idx);
            tmp_bits[blk_idx] = 0;
            tmp_buffer[blk_idx] = 0;
          }
        }
        cnt = 0;
      }
    }

    /* step 2: count remaining one-bits in bit plane */
    int min_m = BLOCK_SIZE;
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      c[blk_idx] = 0;
      if (e[blk_idx] && m[blk_idx] < min_m) min_m = m[blk_idx];
    }
    for (i = min_m; i < BLOCK_SIZE; i++) {
      #pragma _NEC shortloop
      for (blk_idx = 0; blk_idx < vl; blk_idx++) {
        if (i >= m[blk_idx]) {
          c[blk_idx] += (ublock[blk_idx + i * vl] >> k) & 1u;
        }
      }
    }

    /* step 3: unary run-length encode remainder of bit plane */
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
          lbits[blk_idx]--;
          // Output 0 if there are only '0' remaining, state = 0, stop
          // Output 1 if there are at least one '1' remaining, state = 1, continue.
          uint32_t bit = !!c[blk_idx];
          tmp_buffer[blk_idx] += (uint64_t)bit << tmp_bits[blk_idx];
          tmp_bits[blk_idx]++;
          state[blk_idx] = bit;
        }
        if (state[blk_idx] >= 1) {
          if (n[blk_idx] < BLOCK_SIZE - 1 && lbits[blk_idx]) {
            lbits[blk_idx]--;
            // Output current bit
            // If bit = 0, state = 1 + 1 = 2, continue to output '0' without checking remaining bits.
            // If bit = 1, state = 1 + 0, continue and check remaining bits.
            uint32_t bit = (ublock[blk_idx + n[blk_idx] * vl] >> k) & 1u;
            tmp_buffer[blk_idx] += (uint64_t)bit << tmp_bits[blk_idx];
            tmp_bits[blk_idx]++;
            state[blk_idx] = 1 + !bit;
          }
          if (state[blk_idx] != 2) c[blk_idx]--;
          n[blk_idx]++; // Either we have output '0', either we have output '1', so increment position in bit plane.
        }
        flag_rle += state[blk_idx];
      } // for(blk_idx)
      cnt += 2;
      if (cnt >= 63 && flag_rle) {
        // At least 63 bits are ready to be output, flush them to have room for next iteration.
        #pragma _NEC ivdep
        #pragma _NEC shortloop
        for (blk_idx = 0; blk_idx < vl; blk_idx++) {
          if (tmp_bits[blk_idx]) {
            lstream_write_bits(tmp_buffer[blk_idx], tmp_bits[blk_idx], bs_bits, bs_buffer, bs_ptr, blk_idx);
            tmp_bits[blk_idx] = 0;
            tmp_buffer[blk_idx] = 0;
          }
        }
        cnt = 0;
      }
    } // while(flag_rle)

    // Output bits for this plane #k
    #pragma _NEC ivdep
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      if (tmp_bits[blk_idx]) {
        lstream_write_bits(tmp_buffer[blk_idx], tmp_bits[blk_idx], bs_bits, bs_buffer, bs_ptr, blk_idx);
      }
    }

    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      if (e[blk_idx]) {
        if (lbits[blk_idx] && k > kmin[blk_idx]) flag_planes++;
      }
    }
  } // while(flag_planes)
#endif
#ifdef _FTRACE
  ftrace_region_end("COMP - BITSTREAM");
#endif
  
#ifdef _FTRACE
  ftrace_region_begin("COMP - FLUSH BITSTREAM");
#endif
  /*
   * Write at least minbits bits by padding with zeros
   *
   * We support only fixed bit rate, zfp->minbits == zfp->maxbits >= 1
   * We support only encoded block size multiple of 64 bits (zfp->maxbits % 64 == 0)
   */
  /*
   * Flush remaining bits in buffer if any.
   * For empty blocks, bs_bits[blk_idx] == 0
   */
  #pragma _NEC ivdep
  #pragma _NEC shortloop
  for (blk_idx = 0; blk_idx < vl; blk_idx++) {
    if (bs_bits[blk_idx]) {
      lbits[blk_idx] = lbits[blk_idx] - (WSIZE - bs_bits[blk_idx]);
      *bs_ptr[blk_idx]++ = bs_buffer[blk_idx];
      bs_buffer[blk_idx] = 0;
      bs_bits[blk_idx] = 0;
    }
  }
  /*
   * Pad with zeros.
   * For empty block, lbits[blk_idx] == zfp->maxbits
   */
  #pragma _NEC novector
  for (i = 0; i < zfp->maxbits; i += WSIZE) {
    #pragma _NEC ivdep
    #pragma _NEC shortloop
    for (blk_idx = 0; blk_idx < vl; blk_idx++) {
      if (lbits[blk_idx] > 0) {
        *bs_ptr[blk_idx]++ = 0;
        lbits[blk_idx] -= WSIZE;
      }
    }
  }
#ifdef _FTRACE
  ftrace_region_end("COMP - FLUSH BITSTREAM");
#endif
}
