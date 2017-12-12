/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2012 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * @file rdl_config.h
 */

#ifndef RDL_CONFIG_H
#define RDL_CONFIG_H

/* #undef RDL_ENABLE_LOGGING */
/* #undef RDL_BUILD_STATIC */

// The headers then have to be able to work in two different modes:
// - dllexport when one is building the library,
// - dllimport for clients using the library.
//
// On Linux, set the visibility accordingly. If C++ symbol visibility
// is handled by the compiler, see: http://gcc.gnu.org/wiki/Visibility

// On Linux, for GCC >= 4, tag symbols using GCC extension.
#if __GNUC__ >= 4
# define RDL_DLLIMPORT __attribute__((visibility("default")))
# define RDL_DLLEXPORT __attribute__((visibility("default")))
# define RDL_DLLLOCAL  __attribute__((visibility("hidden")))
#else // if __GNUC__ >= 4

// Otherwise (GCC < 4 or another compiler is used), export everything.
# define RDL_DLLIMPORT
# define RDL_DLLEXPORT
# define RDL_DLLLOCAL
#endif // __GNUC__ >= 4

// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#ifdef rdl_EXPORTS
# define RDL_DLLAPI RDL_DLLEXPORT
#else // ifdef rdl_EXPORTS
# define RDL_DLLAPI RDL_DLLIMPORT
#endif // RDL_EXPORTS
#define RDL_LOCAL RDL_DLLLOCAL

#endif // ifndef RDL_CONFIG_H
