#ifndef MD5_H
#define MD5_H

#include <cstddef>
#include <cstdint>
#include <string>

constexpr int MD5_BLOCK_SIZE = 16;  // MD5 outputs a 16 byte digest

using BYTE = unsigned char;  // 8-bit byte
using WORD = unsigned int;  // 32-bit word, change to "long" for 16-bit machines

using MD5_CTX = struct {
  BYTE data[64];
  WORD datalen;
  std::uint64_t bitlen;
  WORD state[4];
};

void md5_init(MD5_CTX *ctx);
void md5_update(MD5_CTX *ctx, const BYTE data[], size_t len);
void md5_final(MD5_CTX *ctx, BYTE hash[]);
void md5_file(const std::string &filename, char hashout[33]);

#endif  // MD5_H
