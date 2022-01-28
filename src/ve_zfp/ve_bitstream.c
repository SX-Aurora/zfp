/* ------------------------------------------------------------------------------------------------ */
/* bit stream structure (opaque to caller) */
struct bitstream {
  uint32_t bits;   /* number of buffered bits (0 <= bits < WSIZE) */
  uint64_t buffer; /* buffer for incoming/outgoing bits (buffer < 2^bits) */
  uint64_t* ptr;   /* pointer to next uint64_t to be read/written */
  uint64_t* begin; /* beginning of stream */
  uint64_t* end;   /* end of stream (currently unused) */
};

/* ------------------------------------------------------------------------------------------------ */
/* write 0 <= n <= 64 bits */
static uint64_t lstream_write_bits(uint64_t value, uint32_t n, uint32_t bs_bits[VL], uint64_t bs_buffer[VL], uint64_t *bs_ptr[VL], int32_t blk_idx)
{
  /* append bit string to buffer */
  bs_buffer[blk_idx] += (uint64_t)(value << bs_bits[blk_idx]);
  bs_bits[blk_idx] += n;
  /* is buffer full? */
  if (bs_bits[blk_idx] >= WSIZE) {
    /* 1 <= n <= 64; decrement n to ensure valid right shifts below */
    value >>= 1;
    n--;
    /* assert: 0 <= n < 64; WSIZE <= s->bits <= WSIZE + n */
//  do {
      /* output WSIZE bits while buffer is full */
      bs_bits[blk_idx] -= WSIZE;
      /* assert: 0 <= s->bits <= n */
      *bs_ptr[blk_idx]++ = bs_buffer[blk_idx];
      /* assert: 0 <= n - s->bits < 64 */
      bs_buffer[blk_idx] = (uint64_t)(value >> (n - bs_bits[blk_idx]));
//  } while (sizeof(bs_buffer[blk_idx]) < sizeof(value) && bs_bits[blk_idx] >= WSIZE);
  }
  /* assert: 0 <= s->bits < WSIZE */
  bs_buffer[blk_idx] &= ((uint64_t)1 << bs_bits[blk_idx]) - 1;
  /* assert: 0 <= n < 64 */
  return value >> n;
}

/* ------------------------------------------------------------------------------------------------ */
/* read 0 <= n <= 64 bits */
static uint64_t lstream_read_bits(uint32_t n, uint32_t bs_bits[VL], uint64_t bs_buffer[VL], uint64_t *bs_ptr[VL], int32_t blk_idx)
{
  uint64_t value = bs_buffer[blk_idx];
  if (bs_bits[blk_idx] < n) {
    /* keep fetching wsize bits until enough bits are buffered */
    // n is at most 64 bits
//  do {
      /* assert: 0 <= s->bits < n <= 64 */
      bs_buffer[blk_idx] = *bs_ptr[blk_idx]++;
      value += (uint64_t)bs_buffer[blk_idx] << bs_bits[blk_idx];
      bs_bits[blk_idx] += WSIZE;
//  } while (sizeof(bs_buffer[blk_idx]) < sizeof(value) && bs_bits[blk_idx] < n);
    /* assert: 1 <= n <= s->bits < n + wsize */
    bs_bits[blk_idx] -= n;
    if (!bs_bits[blk_idx]) {
      /* value holds exactly n bits; no need for masking */
      bs_buffer[blk_idx] = 0;
    }
    else {
      /* assert: 1 <= s->bits < wsize */
      bs_buffer[blk_idx] >>= WSIZE - bs_bits[blk_idx];
      /* assert: 1 <= n <= 64 */
      value &= ((uint64_t)2 << (n - 1)) - 1;
    }
  }
  else {
    /* assert: 0 <= n <= s->bits < wsize <= 64 */
    bs_bits[blk_idx] -= n;
    bs_buffer[blk_idx] >>= n;
    value &= ((uint64_t)1 << n) - 1;
  }
  return value;
}

/* ------------------------------------------------------------------------------------------------ */
/* read single bit (0 or 1) */
static uint32_t lstream_read_bit(uint32_t bs_bits[VL], uint64_t bs_buffer[VL], uint64_t *bs_ptr[VL], int32_t blk_idx)
{
  uint32_t bit;
  if (!bs_bits[blk_idx]) {
    bs_buffer[blk_idx] = *bs_ptr[blk_idx]++;
    bs_bits[blk_idx] = WSIZE;
  }
  bs_bits[blk_idx]--;
  bit = (uint32_t)bs_buffer[blk_idx] & 1u;
  bs_buffer[blk_idx] >>= 1;
  return bit;
}
