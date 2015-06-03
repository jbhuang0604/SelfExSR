#include "QPBO.h"

#ifdef _MSC_VER
#pragma warning(disable: 4661)
typedef __int32 int32_t;
typedef __int64 int64_t;
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

// Instantiations
template class QPBO<int32_t>;
template class QPBO<int64_t>;
template class QPBO<float>;
template class QPBO<double>;

template <> 
	inline void QPBO<int32_t>::get_type_information(char*& type_name, char*& type_format)
{
	type_name = "int32";
	type_format = "d";
}

template <> 
	inline void QPBO<int64_t>::get_type_information(char*& type_name, char*& type_format)
{
	type_name = "int64";
	type_format = "Ld";
}

template <> 
	inline void QPBO<float>::get_type_information(char*& type_name, char*& type_format)
{
	type_name = "float";
	type_format = "f";
}

template <> 
	inline void QPBO<double>::get_type_information(char*& type_name, char*& type_format)
{
	type_name = "double";
	type_format = "Lf";
}


