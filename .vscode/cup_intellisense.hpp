#pragma once

// Ensure CUP headers are parsed with a backend selected.
#if !defined(CUP_BACKEND_QUASI_STATIC) && !defined(CUP_BACKEND_MIE)
#define CUP_BACKEND_QUASI_STATIC
#endif

#include "../modules/cup/cup.hpp"
