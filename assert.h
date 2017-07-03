#undef assert

#ifdef NDEBUG

#define	assert(e)	((void)0)

#else

extern int printout(const char *restrict format, ...);
#define assert(e) if (!(e)) { printout("%s:%u: failed assertion `%s' in function %s\n", __FILE__, __LINE__, #e, __PRETTY_FUNCTION__); abort(); }

#endif