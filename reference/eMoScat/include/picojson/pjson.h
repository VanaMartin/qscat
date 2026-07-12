
// This is the wrapper over PicoJSON -- please use it instead of direct including picojson.h

#ifndef _PICOJSON_PJSON_H
#define _PICOJSON_PJSON_H

#define PICOJSON_USE_INT64
#ifdef __CUDACC__
  #pragma diag_suppress boolean_controlling_expr_is_constant
  //#pragma diag_suppress code_is_unreachable
#endif

#include "picojson/picojson.h"

#ifdef __CUDACC__
  #pragma diag_default boolean_controlling_expr_is_constant
  //#pragma diag_default code_is_unreachable
#endif

#endif // _PICOJSON_PJSON_H