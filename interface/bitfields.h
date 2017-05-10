#ifndef BITFIELDS_H
#define BITFIELDS_H


struct ap_uint3_t {
  // three-bit unsigned field,
  // allowed values are 0...7
  unsigned int b : 3;
};

struct ap_uint5_t {
  // three-bit unsigned field,
  // allowed values are 0...7
  unsigned int b : 5;
};

struct ap_uint11_t {
  // three-bit unsigned field,
  // allowed values are 0...7
  unsigned int b : 11;
};

#endif
