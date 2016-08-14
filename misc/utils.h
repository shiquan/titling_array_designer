#ifndef UTILS_COMMON_HEADER
#define UTILS_COMMON_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include <time.h>

#define check_double_free(p) do {\
        void **_pp = (void**)&(p);                                      \
        if (_pp==NULL || *_pp == NULL) {				\
	    fprintf(stderr,"[error] double free %s, %d", __FUNCTION__, __LINE__); \
	    exit(EXIT_FAILURE);						\
	}								\
    } while(0)

#define ignore_free(p) do \
	void **_pp = (void**)&(p);					\
	if (*_pp!=NULL && _pp != NULL) {				\
	    free(*_pp);							\
	    *_pp = NULL;						\
	}								\
    } while(0)

#define safe_free(p) do				\
    {						\
	void **_pp = (void**)&(p);					\
        if (_pp==NULL || *_pp == NULL) {				\
	    fprintf(stderr,"[error] double free %s, %d", __FUNCTION__, __LINE__); \
	    exit(EXIT_FAILURE);						\
	}								\
	free(*_pp);							\
        *_pp = NULL;                                                    \
    } while(0)

#define check_mem(p) do				\
    {						\
	void **_pp = (void**)&p;		\
	if (_pp == NULL || *_pp == NULL) {				\
	    fprintf(stderr, "[memory out] func: %s, line: %d\n", __FUNCTION__, __LINE__);\
	    exit(EXIT_FAILURE);						\
	}								\
    }while(0)

#define str_errno() (errno == 0 ? "None" : strerror(errno))

#define clear_errno() do \
    {\
	fprintf(stderr, "%s\n", str_errno());\
	errno = 0;\
    }while(0)

#define error(_line, ...) do						\
    {									\
	fprintf(stderr, "[error] [func: %s, line: %d] " _line "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
	errno = 0;							\
	exit(EXIT_FAILURE);						\
    }while(0)

#define error_print(_line, ...) do					\
    {									\
	fprintf(stderr, "[error] [func: %s, line: %d] " _line "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    }while(0)

#define warnings(_line, ...) do						\
    {									\
	if (errno == 0) {						\
	    fprintf(stderr, "[warnings] " _line "\n", ##__VA_ARGS__);	\
	} else {							\
	    fprintf(stderr, "[warnings] Errno: %s. " _line "\n", str_errno(), ##__VA_ARGS__); \
	}								\
    }while(0)

#define debug_print(line, ...) do {\
	fprintf(stderr, "[ ** DEBUG ** func: %s, line: %d ] " line "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    } while(0)

#define LOG_print(_line, ...) do {\
	time_t second;\
	time(&second);\
	char _buff[100];							\
	strftime (_buff, 100, "%Y-%m-%d %H:%M:%S", localtime (&second));	\
	fprintf(stderr, "[%s] " _line "\n", _buff, ##__VA_ARGS__); \
    } while(0)

#define BE_SMART_STRING "Please DONOT post this error message on the forum or copy it into the emails. Try to figure out this issue by youself by reading the log information carefully and checking you input arguments."

#endif
