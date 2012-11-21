#ifndef __INI_H__
#define __INI_H__

#include <stdio.h>

//Header file for the input.txt file reader.  The input file should be formated in the form   [section]s, name=value (whitespace stripped)
//Comments starting with ';' or '#'.
//The section is "" if name=value pair is parsed before any section heading.

//The input file reader works as follows. For each name=value pair, the call handler function is called with given conf pointer as well as section, name, and value. Handler should return (not as usually) nonzero(this is a line information)  on success, zero on error.

//Returns 0 on success, line number of first error on parse error, -1 on file open error.

int ini_parse(const char* filename,
              int (*handler)(void* config, const char* section,
                             const char* name, const char* value),
              void* conf);

//Same as ini_parse(), but takes a FILE* instead of filename. This doesn not  close the file when it's finished. The caller have to do that.
int ini_parse_file(FILE* file,
                   int (*handler)(void* config, const char* section,
                                  const char* name, const char* value),
                   void* conf);

// Maximum line length for any line in input file
#ifndef INI_MAX_LINE
#define INI_MAX_LINE 200
#endif

#endif /* __INI_H__ */
