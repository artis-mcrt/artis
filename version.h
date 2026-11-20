// The git metadata and the compile-time options of the build. The Makefile writes the
// definitions into the generated version.cc. The git metadata changes at almost every git
// operation. A definition in this header makes sn3d.cc and exspec.cc compile again each time.
// Keep the definitions in version.cc.

#ifndef VERSION_H
#define VERSION_H

#include <string_view>

extern const char* const GIT_VERSION;
extern const char* const GIT_BRANCH;
extern const char* const GIT_STATUS;
// the file that the symlink artisoptions.h points to, or "not a symlink"
extern const char* const ARTISOPTIONS_SYMLINK_TARGET;
// the first line of the text of "mpicxx --version"
extern const char* const COMPILER_VERSION_LINE;
extern const char* const COMPILE_FLAGS;
// the lines of artisoptions.h, separated by newlines
extern const std::string_view ARTISOPTIONS_LINES;

#endif  // VERSION_H
