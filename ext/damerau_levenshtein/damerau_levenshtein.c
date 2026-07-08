#include "ruby.h"
#include <limits.h>
#include <stdint.h>

VALUE DamerauLevenshteinBinding = Qnil;

void Init_damerau_levenshtein();

VALUE method_internal_distance(VALUE self, VALUE _s, VALUE _t, VALUE _block_size, VALUE _max_distance);

void Init_damerau_levenshtein() {
	DamerauLevenshteinBinding = rb_define_module("DamerauLevenshteinBinding");
	rb_define_method(DamerauLevenshteinBinding, "internal_distance", method_internal_distance, 4);
}

VALUE method_internal_distance(VALUE self, VALUE _s, VALUE _t, VALUE _block_size, VALUE _max_distance){
  int i, i1, j, j1, k, half_tl, cost, *d, distance, del, ins, subs, transp, block;
  int half_sl;
  int stop_execution = 0;
  int min = 0;
  int current_distance = 0;
  int pure_levenshtein = 0;
  int s_len, t_len, width, height;
  long long *s, *t;
  size_t width_z, matrix_len;

  Check_Type(_s, T_ARRAY);
  Check_Type(_t, T_ARRAY);

  int block_size = NUM2INT(_block_size);
  int max_distance = NUM2INT(_max_distance);
  if (block_size < 0) rb_raise(rb_eArgError, "block_size must not be negative");
  if (max_distance < 0) rb_raise(rb_eArgError, "max_distance must not be negative");

  if (RARRAY_LEN(_s) >= INT_MAX || RARRAY_LEN(_t) >= INT_MAX)
    rb_raise(rb_eArgError, "input arrays are too large");
  s_len = (int) RARRAY_LEN(_s);
  t_len = (int) RARRAY_LEN(_t);

  if (block_size == 0) {
    pure_levenshtein = 1;
    block_size = 1;
  }


  if (s_len == 0) return INT2NUM(t_len);
  if (t_len == 0) return INT2NUM(s_len);
  //case of lengths 1 must present or it will break further in the code
  if (s_len == 1 && t_len == 1 && RARRAY_AREF(_s, 0) != RARRAY_AREF(_t, 0)) return INT2NUM(1);

  //the length difference alone already exceeds the threshold; skip the
  //matrix work and report over-max the same way the early-stop path does
  if ((s_len > t_len ? s_len - t_len : t_len - s_len) > max_distance)
    return INT2NUM(max_distance + 1);

  s = ALLOC_N(long long, s_len);
  t = ALLOC_N(long long, t_len);

  //RARRAY_AREF (not a cached RARRAY_PTR) on each iteration: the ALLOC_N
  //calls above and NUM2LL below may trigger GC
  for (i=0; i < s_len; i++) s[i] = NUM2LL(RARRAY_AREF(_s, i));
  for (i=0; i < t_len; i++) t[i] = NUM2LL(RARRAY_AREF(_t, i));

  width = s_len + 1;
  height = t_len + 1;
  width_z = (size_t)width;
  if ((size_t)width > SIZE_MAX / sizeof(int) / (size_t)height) {
    xfree(s);
    xfree(t);
    rb_raise(rb_eNoMemError, "distance matrix is too large");
  }
  matrix_len = width_z * (size_t)height;

  //one-dimentional representation of 2 dimentional array len(s)+1 * len(t)+1
  d = ALLOC_N(int, matrix_len);
  //populate 'vertical' row starting from the 2nd position (first one is filled already)
  for(i = 0; i < height; i++){
    d[(size_t)i*width_z] = i;
  }

  //fill up array with scores
  for(i = 1; i<width; i++){
    d[i] = i;
    if (stop_execution == 1) break;
    current_distance = 10000;
    for(j = 1; j<height; j++){

      cost = 1;
      if(s[i-1] == t[j-1]) cost = 0;

      half_sl = s_len/2;
      half_tl = t_len/2;

      block = block_size < half_sl || half_sl == 0 ? block_size : half_sl;
      block = block < half_tl || half_tl == 0 ? block : half_tl;

      while (block >= 1){
        int swap1 = 1;
        int swap2 = 1;
        i1 = i - (block * 2);
        j1 = j - (block * 2);
        //the loops below read s[i1 .. i-block-1] and t[i-block .. i-1]
        //(and vice versa), so skip them when those ranges fall outside
        //either array
        if (i1 < 0 || i > t_len) {
          swap1 = 0;
        } else {
          for (k = i1; k < i1 + block; k++) {
            if (s[k] != t[k + block]){
              swap1 = 0;
              break;
            }
          }
        }
        if (j1 < 0 || j > s_len) {
          swap2 = 0;
        } else {
          for (k = j1; k < j1 + block; k++) {
            if (t[k] != s[k + block]){
              swap2 = 0;
              break;
            }
          }
        }

        del = d[(size_t)j*width_z + i - 1] + 1;
        ins = d[(size_t)(j-1)*width_z + i] + 1;
        min = del;
        if (ins < min) min = ins;
        //if (i == 2 && j==2) return INT2NUM(swap2+5);
        if (pure_levenshtein == 0 && i >= 2*block && j >= 2*block && swap1 == 1 && swap2 == 1){
          transp = d[(size_t)(j - block * 2) * width_z + i - block * 2] + cost + block -1;
          if (transp < min) min = transp;
          block = 0;
        } else if (block == 1) {
          subs = d[(size_t)(j-1)*width_z + i - 1] + cost;
          if (subs < min) min = subs;
        }
        block--;
      }
      d[(size_t)j*width_z + i] = min;
      if (current_distance > d[(size_t)j*width_z + i]) current_distance = d[(size_t)j*width_z + i];
    }
    if (current_distance > max_distance) {
      stop_execution = 1;
    }
  }
  distance = d[matrix_len - 1];
  if (stop_execution == 1) distance = current_distance;

  xfree(d);
  xfree(s);
  xfree(t);
  return INT2NUM(distance);
}
