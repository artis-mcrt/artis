// The git metadata of the build. The Makefile writes the definitions into the generated
// version.cc. The values change at almost every git operation. A definition in this header
// makes sn3d.cc and exspec.cc compile again each time, so keep the definitions in version.cc.

#ifndef VERSION_H
#define VERSION_H

extern const char* const GIT_VERSION;
extern const char* const GIT_BRANCH;
extern const char* const GIT_STATUS;

#endif  // VERSION_H
