#ifndef MD5_H
#define MD5_H

#include <cstddef>
#include <string>

auto md5_file(const std::string &filename) -> std::string;
void md5_test();

#endif  // MD5_H
