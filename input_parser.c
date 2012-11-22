// A very simple  file parser for the config file
#ifndef INI_C
#define INI_C
#include <stdio.h>
#include <ctype.h>
#include <string.h>
  
#include <stdlib.h>
#endif
  
//This set a upper maximum for the config parser . MAX_SECTION defines the maximum number of sections. MAX_NAME defines the maximum number of names per section.
#define MAX_SECTION 50
#define MAX_NAME 50
  
// Strip whitespace chars off end of given string, in place. Return s.
static char *
rstrip (char *s) 
{
  char *p = s + strlen (s); //pointer to whitespace 
  while (p > s && isspace (*--p)) // go through the string via pointer arithmetic 
    *p = '\0';//de
  return s;
}

 
// Return a pointer to the first non-whitespace char in given string.
static char *
lskip (const char *s) 
{
  while (*s && isspace (*s)) //
    s++;
  return (char *) s;
}  

// This function return a pointer to the first char c or ';' comment in given string, or pointer to
// null at the end of the string if neither found. ';' must be prefixed by a whitespace
// character to register as a comment.
static char *
find_char_or_comment (const char *s, char c) 
{
  int was_whitespace = 0;
  while (*s && *s != c && !(was_whitespace && *s == ';'))
    {
      was_whitespace = isspace (*s);
      s++;
    }
  return (char *) s;
}  

//Version of strncpy that ensures dest (size bytes) is null-terminated.
static char *
strncpy0 (char *dest, const char *src, size_t size) 
{
  strncpy (dest, src, size);
  dest[size - 1] = '\0';
  return dest;
}

 
// See documentation in header file.
  int
ini_parse_file (FILE * file,
		int (*handler) (void *, const char *, const char *,
				 const char *), void *conf) 
{
  
//Uses a fair bit of stack (use heap instead if you need to) 
  char *line;
  char section[MAX_SECTION] = "";
  char prev_name[MAX_NAME] = "";
   char *start;
  char *end;
  char *name;
  char *value;
  int lineno = 0;
  int error = 0;
   
    line = (char *) malloc (INI_MAX_LINE);
  if (!line)
    {
      return -2;
    }
  
    
    // Scan through file line by line
    while (fgets (line, INI_MAX_LINE, file) != NULL)
    {
      lineno++;
       start = line;
      start = lskip (rstrip (start));
       if (*start == ';' || *start == '#')
	{
	  
	    // I add this to allow # as  comments at start of line 
	}
      
      else if (*start == '[')
	{
	  
	    //Find the section via [  A. The end is set by ]
	    end = find_char_or_comment (start + 1, ']');
	  if (*end == ']')
	    {
	      *end = '\0';
	      strncpy0 (section, start + 1, sizeof (section));
	      *prev_name = '\0';
	    }
	  
	  else if (!error)
	    {
	      
		// If there is no ] in the line with [ the parser will crash. 
		error = lineno;
	    }
	}
      
      else if (*start && *start != ';')
	{
	  
	    //Find the values by point to = if there is a text in front of it 
	    end = find_char_or_comment (start, '=');
	  if (*end != '=')
	    {
	      end = find_char_or_comment (start, ':');
	    }
	  if (*end == '=' || *end == ':')
	    {
	      *end = '\0';
	      name = rstrip (start);
	      value = lskip (end + 1);
	      end = find_char_or_comment (value, '\0');
	      if (*end == ';')
		*end = '\0';
	      rstrip (value);
	       
		// Valid name[=:]value pair found, call handler 
		strncpy0 (prev_name, name, sizeof (prev_name));
	      if (!handler (conf, section, name, value) && !error)
		error = lineno;
	    }
	  
	  else if (!error)
	    {
	      
		// No = or : found on name value line. This will lead to a error 
		error = lineno;
	    }
	}
    }
   
    free (line);
  
    return error;
}

 
  int
ini_parse (const char *filename,
	   int (*handler) (void *, const char *, const char *, const char *),
	   void *conf) 
{
  FILE * file;
  int error;
   file = fopen (filename, "r");
  if (!file)
    return -1;
  error = ini_parse_file (file, handler, conf);
  fclose (file);
  return error;
}


