//#define NTEMP1(x) NTEMP(x)
//#define COUNT_END1(x) COUNT_END(x)
size_t COUNT_END(FUN_NAME) = __COUNTER__;
size_t NTEMP(FUN_NAME) {
	return COUNT_END(FUN_NAME) + (NPASS>0?NPASS-1:0) - PREV_END - 1;
	#undef NPASS
	#define NPASS 0
	#undef NEXT_TEMP_INDEX
	#define NEXT_TEMP_INDEX __COUNTER__ - COUNT_END(FUN_NAME1) - 1
	#undef PREV_END
	#define PREV_END COUNT_END(FUN_NAME1)
}
#undef FUN_NAME
#undef FUN_NAME1