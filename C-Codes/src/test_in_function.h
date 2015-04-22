#ifndef _TEST_INT_FUNCTION_H_INCLUDED_
#define _TEST_INT_FUNCTION_H_INCLUDED_

#include <stdio.h>

#define TEST_IN_FUNCTION_START { printf("***Start test in %s[%d]:\nfunction: %s\n",__FILE__, __LINE__,__FUNCTION__);
#define TEST_IN_FUNCTION_END printf("End test in %s[%d]***\n",__FILE__, __LINE__);}

#endif
