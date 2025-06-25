// Based on the public domain MD5 implementation by Brad Conte
// https://github.com/B-Con/crypto-algorithms/blob/master/md5.c

#include "md5.h"

#include <array>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <span>
#include <string>

#include "sn3d.h"

namespace {
constexpr int MD5_BLOCK_SIZE = 16;  // MD5 outputs a 16 byte digest

using BYTE = unsigned char;  // 8-bit byte
using WORD = std::uint32_t;  // 32-bit word, change to "long" for 16-bit machines

using MD5_CTX = struct {
  std::array<BYTE, 64> data;
  WORD datalen;
  std::uint64_t bitlen;
  WORD state[4];
};

constexpr auto F(const WORD x, const WORD y, const WORD z) -> WORD { return (x & y) | (~x & z); }
constexpr auto G(const WORD x, const WORD y, const WORD z) -> WORD { return (x & z) | (y & ~z); }
constexpr auto H(const WORD x, const WORD y, const WORD z) -> WORD { return x ^ y ^ z; }
constexpr auto I(const WORD x, const WORD y, const WORD z) -> WORD { return y ^ (x | ~z); }

constexpr auto ff_func(WORD a, const WORD b, const WORD c, const WORD d, const WORD m, const int s, const WORD t)
    -> WORD {
  a += F(b, c, d) + m + t;
  a = b + std::rotl(a, s);
  return a;
}

constexpr auto gg_func(WORD a, const WORD b, const WORD c, const WORD d, const WORD m, const int s, const WORD t)
    -> WORD {
  a += G(b, c, d) + m + t;
  a = b + std::rotl(a, s);
  return a;
}

constexpr auto hh_func(WORD a, const WORD b, const WORD c, const WORD d, const WORD m, const int s, const WORD t)
    -> WORD {
  a += H(b, c, d) + m + t;
  a = b + std::rotl(a, s);
  return a;
}

constexpr auto ii_func(WORD a, const WORD b, const WORD c, const WORD d, const WORD m, const int s, const WORD t)
    -> WORD {
  a += I(b, c, d) + m + t;
  a = b + std::rotl(a, s);
  return a;
}

/*********************** FUNCTION DEFINITIONS ***********************/
constexpr void md5_transform(MD5_CTX *ctx, std::span<const BYTE> data) {
  WORD a = 0;
  WORD b = 0;
  WORD c = 0;
  WORD d = 0;
  WORD m[16];
  WORD i = 0;
  WORD j = 0;

  // MD5 specifies big endian byte order, but this implementation assumes a little
  // endian byte order CPU. Reverse all the bytes upon input, and re-reverse them
  // on output (in md5_final()).
  for (i = 0, j = 0; i < 16; ++i, j += 4) {
    m[i] = (data[j]) + (data[j + 1] << 8) + (data[j + 2] << 16) + (data[j + 3] << 24);
  }

  a = ctx->state[0];
  b = ctx->state[1];
  c = ctx->state[2];
  d = ctx->state[3];

  a = ff_func(a, b, c, d, m[0], 7, 0xd76aa478);
  d = ff_func(d, a, b, c, m[1], 12, 0xe8c7b756);
  c = ff_func(c, d, a, b, m[2], 17, 0x242070db);
  b = ff_func(b, c, d, a, m[3], 22, 0xc1bdceee);
  a = ff_func(a, b, c, d, m[4], 7, 0xf57c0faf);
  d = ff_func(d, a, b, c, m[5], 12, 0x4787c62a);
  c = ff_func(c, d, a, b, m[6], 17, 0xa8304613);
  b = ff_func(b, c, d, a, m[7], 22, 0xfd469501);
  a = ff_func(a, b, c, d, m[8], 7, 0x698098d8);
  d = ff_func(d, a, b, c, m[9], 12, 0x8b44f7af);
  c = ff_func(c, d, a, b, m[10], 17, 0xffff5bb1);
  b = ff_func(b, c, d, a, m[11], 22, 0x895cd7be);
  a = ff_func(a, b, c, d, m[12], 7, 0x6b901122);
  d = ff_func(d, a, b, c, m[13], 12, 0xfd987193);
  c = ff_func(c, d, a, b, m[14], 17, 0xa679438e);
  b = ff_func(b, c, d, a, m[15], 22, 0x49b40821);

  a = gg_func(a, b, c, d, m[1], 5, 0xf61e2562);
  d = gg_func(d, a, b, c, m[6], 9, 0xc040b340);
  c = gg_func(c, d, a, b, m[11], 14, 0x265e5a51);
  b = gg_func(b, c, d, a, m[0], 20, 0xe9b6c7aa);
  a = gg_func(a, b, c, d, m[5], 5, 0xd62f105d);
  d = gg_func(d, a, b, c, m[10], 9, 0x02441453);
  c = gg_func(c, d, a, b, m[15], 14, 0xd8a1e681);
  b = gg_func(b, c, d, a, m[4], 20, 0xe7d3fbc8);
  a = gg_func(a, b, c, d, m[9], 5, 0x21e1cde6);
  d = gg_func(d, a, b, c, m[14], 9, 0xc33707d6);
  c = gg_func(c, d, a, b, m[3], 14, 0xf4d50d87);
  b = gg_func(b, c, d, a, m[8], 20, 0x455a14ed);
  a = gg_func(a, b, c, d, m[13], 5, 0xa9e3e905);
  d = gg_func(d, a, b, c, m[2], 9, 0xfcefa3f8);
  c = gg_func(c, d, a, b, m[7], 14, 0x676f02d9);
  b = gg_func(b, c, d, a, m[12], 20, 0x8d2a4c8a);

  a = hh_func(a, b, c, d, m[5], 4, 0xfffa3942);
  d = hh_func(d, a, b, c, m[8], 11, 0x8771f681);
  c = hh_func(c, d, a, b, m[11], 16, 0x6d9d6122);
  b = hh_func(b, c, d, a, m[14], 23, 0xfde5380c);
  a = hh_func(a, b, c, d, m[1], 4, 0xa4beea44);
  d = hh_func(d, a, b, c, m[4], 11, 0x4bdecfa9);
  c = hh_func(c, d, a, b, m[7], 16, 0xf6bb4b60);
  b = hh_func(b, c, d, a, m[10], 23, 0xbebfbc70);
  a = hh_func(a, b, c, d, m[13], 4, 0x289b7ec6);
  d = hh_func(d, a, b, c, m[0], 11, 0xeaa127fa);
  c = hh_func(c, d, a, b, m[3], 16, 0xd4ef3085);
  b = hh_func(b, c, d, a, m[6], 23, 0x04881d05);
  a = hh_func(a, b, c, d, m[9], 4, 0xd9d4d039);
  d = hh_func(d, a, b, c, m[12], 11, 0xe6db99e5);
  c = hh_func(c, d, a, b, m[15], 16, 0x1fa27cf8);
  b = hh_func(b, c, d, a, m[2], 23, 0xc4ac5665);

  a = ii_func(a, b, c, d, m[0], 6, 0xf4292244);
  d = ii_func(d, a, b, c, m[7], 10, 0x432aff97);
  c = ii_func(c, d, a, b, m[14], 15, 0xab9423a7);
  b = ii_func(b, c, d, a, m[5], 21, 0xfc93a039);
  a = ii_func(a, b, c, d, m[12], 6, 0x655b59c3);
  d = ii_func(d, a, b, c, m[3], 10, 0x8f0ccc92);
  c = ii_func(c, d, a, b, m[10], 15, 0xffeff47d);
  b = ii_func(b, c, d, a, m[1], 21, 0x85845dd1);
  a = ii_func(a, b, c, d, m[8], 6, 0x6fa87e4f);
  d = ii_func(d, a, b, c, m[15], 10, 0xfe2ce6e0);
  c = ii_func(c, d, a, b, m[6], 15, 0xa3014314);
  b = ii_func(b, c, d, a, m[13], 21, 0x4e0811a1);
  a = ii_func(a, b, c, d, m[4], 6, 0xf7537e82);
  d = ii_func(d, a, b, c, m[11], 10, 0xbd3af235);
  c = ii_func(c, d, a, b, m[2], 15, 0x2ad7d2bb);
  b = ii_func(b, c, d, a, m[9], 21, 0xeb86d391);

  ctx->state[0] += a;
  ctx->state[1] += b;
  ctx->state[2] += c;
  ctx->state[3] += d;
}

constexpr void md5_init(MD5_CTX *ctx) {
  ctx->datalen = 0;
  ctx->bitlen = 0;
  ctx->state[0] = 0x67452301;
  ctx->state[1] = 0xEFCDAB89;
  ctx->state[2] = 0x98BADCFE;
  ctx->state[3] = 0x10325476;
}

constexpr void md5_update(MD5_CTX *ctx, std::span<const BYTE> data) {
  size_t i = 0;

  for (i = 0; i < data.size(); ++i) {
    ctx->data[ctx->datalen] = data[i];
    ctx->datalen++;
    if (ctx->datalen == 64) {
      md5_transform(ctx, ctx->data);
      ctx->bitlen += 512;
      ctx->datalen = 0;
    }
  }
}

constexpr void md5_final(MD5_CTX *ctx, std::span<BYTE, MD5_BLOCK_SIZE> hash) {
  size_t i = 0;

  i = ctx->datalen;

  // Pad whatever data is left in the buffer.
  if (ctx->datalen < 56) {
    ctx->data[i++] = 0x80;
    while (i < 56) {
      ctx->data[i++] = 0x00;
    }
  } else if (ctx->datalen >= 56) {
    ctx->data[i++] = 0x80;
    while (i < 64) {
      ctx->data[i++] = 0x00;
    }
    md5_transform(ctx, ctx->data);
    memset(ctx->data.data(), 0, 56);
  }

  // Append to the padding the total message's length in bits and transform.
  ctx->bitlen += ctx->datalen * 8;
  ctx->data[56] = ctx->bitlen;
  ctx->data[57] = ctx->bitlen >> 8;
  ctx->data[58] = ctx->bitlen >> 16;
  ctx->data[59] = ctx->bitlen >> 24;
  ctx->data[60] = ctx->bitlen >> 32;
  ctx->data[61] = ctx->bitlen >> 40;
  ctx->data[62] = ctx->bitlen >> 48;
  ctx->data[63] = ctx->bitlen >> 56;
  md5_transform(ctx, ctx->data);

  // Since this implementation uses little endian byte ordering and MD uses big endian,
  // reverse all the bytes when copying the final state to the output hash.
  for (i = 0; i < 4; ++i) {
    hash[i] = (ctx->state[0] >> (i * 8)) & 0x000000ff;
    hash[i + 4] = (ctx->state[1] >> (i * 8)) & 0x000000ff;
    hash[i + 8] = (ctx->state[2] >> (i * 8)) & 0x000000ff;
    hash[i + 12] = (ctx->state[3] >> (i * 8)) & 0x000000ff;
  }
}

}  // anonymous namespace

// added by Luke Shingles
auto md5_file(const std::string &filename) -> std::string {
  MD5_CTX ctx;
  md5_init(&ctx);

  FILE *infile = fopen(filename.c_str(), "r");

  assert_always(infile != nullptr);

  std::array<BYTE, 1024> buffer{};

  size_t numbytes = 1;
  while (numbytes != 0 && feof(infile) == 0) {
    numbytes = fread(buffer.data(), sizeof(BYTE), buffer.size(), infile);
    assert_always(ferror(infile) == 0);
    md5_update(&ctx, std::span{buffer}.first(numbytes));
  }

  fclose(infile);

  std::array<BYTE, MD5_BLOCK_SIZE> hashbytes{};
  md5_final(&ctx, hashbytes);

  std::string hashout(2 * MD5_BLOCK_SIZE, '0');  // Initialize with zeros
  for (int j = 0; j < MD5_BLOCK_SIZE; j++) {
    snprintf(&hashout[2 * j], (2 * MD5_BLOCK_SIZE) + 1 - (2 * j), "%02x", hashbytes[j]);
  }

  return hashout;
}

void md5_test() {
  MD5_CTX ctx;
  md5_init(&ctx);

  constexpr BYTE buffer[] = "md5 test string\n";

  md5_update(&ctx, std::span{buffer}.first(sizeof(buffer) - 1));

  std::array<BYTE, MD5_BLOCK_SIZE> hashbytes{};
  md5_final(&ctx, hashbytes);

  std::string hashout(2 * MD5_BLOCK_SIZE, '0');  // Initialize with zeros
  for (int j = 0; j < MD5_BLOCK_SIZE; j++) {
    snprintf(&hashout[2 * j], (2 * MD5_BLOCK_SIZE) + 1 - (2 * j), "%02x", hashbytes[j]);
  }
  assert_always(hashout == "4b0cff9625b0501c7b9ccc6569113ddf");
}
