#ifndef _VE_FUNCTIONS_H
#define _VE_FUNCTIONS_H

extern void compress_ve_float_1(zfp_stream* stream, const zfp_field* field);
extern void compress_ve_float_2(zfp_stream* stream, const zfp_field* field);
extern void compress_ve_float_3(zfp_stream* stream, const zfp_field* field);
extern void compress_ve_float_4(zfp_stream* stream, const zfp_field* field);

extern void decompress_ve_float_1(zfp_stream* stream, zfp_field* field);
extern void decompress_ve_float_2(zfp_stream* stream, zfp_field* field);
extern void decompress_ve_float_3(zfp_stream* stream, zfp_field* field);
extern void decompress_ve_float_4(zfp_stream* stream, zfp_field* field);

extern void compress_ve_double_1(zfp_stream* stream, const zfp_field* field);
extern void compress_ve_double_2(zfp_stream* stream, const zfp_field* field);
extern void compress_ve_double_3(zfp_stream* stream, const zfp_field* field);
extern void compress_ve_double_4(zfp_stream* stream, const zfp_field* field);

extern void decompress_ve_double_1(zfp_stream* stream, zfp_field* field);
extern void decompress_ve_double_2(zfp_stream* stream, zfp_field* field);
extern void decompress_ve_double_3(zfp_stream* stream, zfp_field* field);
extern void decompress_ve_double_4(zfp_stream* stream, zfp_field* field);

#endif // _VE_FUNCTIONS_H
