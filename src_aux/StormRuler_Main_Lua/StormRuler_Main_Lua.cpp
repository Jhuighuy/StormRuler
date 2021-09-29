// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// Copyright (C) 2021 Oleg Butakov
// 
// Permission is hereby granted, free of charge, to any person 
// obtaining a copy of this software and associated documentation 
// files (the "Software"), to deal in the Software without 
// restriction, including without limitation the rights  to use, 
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the  
// Software is furnished to do so, subject to the following 
// conditions:
// 
// The above copyright notice and this permission notice shall be 
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#include <StormRuler_API.h>

#include <iostream>
#include <type_traits>

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
} // extern "C"

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

template <typename T>
auto type_name(T&&) noexcept {
  std::string name;
#ifdef __clang__
  name = __PRETTY_FUNCTION__;
#elif defined(__GNUC__)
  name = __PRETTY_FUNCTION__;
#elif defined(_MSC_VER)
  name = __FUNCSIG__;
#endif
  return name;
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/// @{
inline void PushReturn(lua_State* L, lua_Number ret) {
  lua_pushnumber(L, ret);
}
inline void PushReturn(lua_State* L, lua_Integer ret) {
  lua_pushinteger(L, ret);
}
inline void PushReturn(lua_State* L, void* ret) {
  lua_pushlightuserdata(L, ret);
}
/// @}

/// @{
template<typename tType>
inline tType UnwrapArgument(lua_State* L, int iArg)  {
  static_assert( std::is_pointer_v<tType> );
  return reinterpret_cast<tType>(lua_touserdata(L, iArg));
}
template<>
inline SR_REAL UnwrapArgument<SR_REAL>(lua_State* L, int iArg) {
  return SR_REAL( luaL_checknumber(L, iArg) );
}
template<>
inline SR_INTEGER UnwrapArgument<SR_INTEGER>(lua_State* L, int iArg) {
  return SR_INTEGER( luaL_checkinteger(L, iArg) );
}
/// @}

/// @{
template<typename tFunc>
struct tLuaFunc;
template<typename tRet, typename... tArgs>
struct tLuaFunc<tRet(tArgs...)> {
  template<tRet(func)(tArgs...)>
  static int Call(lua_State* L) {
    if (lua_gettop(L) != sizeof...(tArgs)) {
      return luaL_error(L, "Unexpected amount of arguments.");
    }
    int iArg = 0;
    if constexpr (!std::is_void_v<tRet>) {
      tRet ret = func((++iArg, UnwrapArgument<tArgs>(L, iArg))...);
      PushReturn(L, ret);
      return 1;
    } else {
      func((++iArg, UnwrapArgument<tArgs>(L, iArg))...);
      return 0;
    }
  }
}; // tLuaFunc
#define SR_WrapLua(func) ( &tLuaFunc<decltype(func)>::Call<func> )
/// @}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

static const luaL_Reg gFuncList[] = {
    // ---------------------------------------- //
    { "InitMesh",    SR_WrapLua(SR_InitMesh)     },
    // ---------------------------------------- //
    { "AllocR",      SR_WrapLua(SR_AllocR)       },
    { "AllocC",      SR_WrapLua(SR_AllocC)       },
    { "AllocS",      SR_WrapLua(SR_AllocR)       },
    { "Alloc_Mold",  SR_WrapLua(SR_Alloc_MoldR)  },
    { "Free",        SR_WrapLua(SR_FreeR)        },
    // ---------------------------------------- //
    { "Fill",        SR_WrapLua(SR_FillR)        },
    { "Fill_Random", SR_WrapLua(SR_Fill_RandomR) },
    { "Set",         SR_WrapLua(SR_SetR)         },
    { "Scale",       SR_WrapLua(SR_ScaleR)       },
    { "Add",         SR_WrapLua(SR_AddR)         },
    { "Sub",         SR_WrapLua(SR_SubR)         },
    // Mul, Div
    { "FuncProd",    SR_WrapLua(SR_FuncProdR)    },
    // ---------------------------------------- //
    { "ApplyBCs",    SR_WrapLua(SR_ApplyBCsR)    },
    { "Grad",        SR_WrapLua(SR_GradR)        },
    { "Div",         SR_WrapLua(SR_DivR)         },
    { "DivGrad",     SR_WrapLua(SR_DivGradR)     },
    { "DivKGrad",    SR_WrapLua(SR_DivKGradR)    },
    // ---------------------------------------- //
    { nullptr,       nullptr                     },
  };

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

int main(int argc, char** argv) {
  lua_State* L = luaL_newstate();
  luaL_openlibs(L);

  // Exposing the API.
  lua_newtable(L);
  luaL_setfuncs(L, gFuncList, 0);
  lua_setglobal(L, "SR");

  // Executing the script.
  int ret = luaL_dofile(L, argv[1]);
  if (ret != 0) {
    std::cerr << "Error occurs while executing luaL_dofile, RET=" << ret << std::endl;
    std::cerr << "Error: " << lua_tostring(L, -1) << std::endl;
    return 1;
  }

  // Gracefully exiting.
  std::cout << "Exiting StormRuler Lua" << std::endl;
  lua_close(L);
  return 0;
}
