/*********************************************************************
 * Filename:   md5.h
 * Author:     Brad Conte (brad AT bradconte.com)
 * Copyright:
 * Disclaimer: This code is presented "as is" without any guarantees.
 * Details:    Defines the API for the corresponding MD5 implementation.
 *********************************************************************/

#ifndef MD5_H
#define MD5_H

/*************************** HEADER FILES ***************************/
#include <cstddef>

/****************************** MACROS ******************************/
constexpr int MD5_BLOCK_SIZE = 16;  // MD5 outputs a 16 byte digest

/**************************** DATA TYPES ****************************/
using BYTE = unsigned char;  // 8-bit byte
using WORD = unsigned int;   // 32-bit word, change to "long" for 16-bit machines

using MD5_CTX = struct {
  BYTE data[64];
  WORD datalen;
  unsigned long long bitlen;
  WORD state[4];
};

/*********************** FUNCTION DECLARATIONS **********************/
void md5_init(MD5_CTX *ctx);
void md5_update(MD5_CTX *ctx, const BYTE data[], size_t len);
void md5_final(MD5_CTX *ctx, BYTE hash[]);
void md5_file(const char filename[], char hashout[33]);

#endif  // MD5_H
