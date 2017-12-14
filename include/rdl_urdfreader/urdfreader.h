/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_URDFREADER_H__
#define __RDL_URDFREADER_H__

#include <rdl_dynamics/rdl_config.h>

namespace RobotDynamics
{
    struct Model;
    namespace Urdf
    {
        /**
         * @brief Read urdf from file path
         * @param filename
         * @param model
         * @param floating_base
         * @param verbose
         * @return
         */
        RDL_DLLAPI bool urdfReadFromFile(const char *filename, Model *model, bool floating_base, bool verbose = false);

        /**
         * @brief Read urdf from string contents
         * @param model_xml_string
         * @param model
         * @param floating_base
         * @param verbose
         * @return
         */
        RDL_DLLAPI bool urdfReadFromString(const char *model_xml_string, Model *model, bool floating_base, bool verbose = false);

        /**
         * @brief This will build a map of joint name to body name.
         * @param model_xml_string Urdf file contents
         * @param jointBodyMap Modified
         * @return
         *
         * @warning This will NOT give any information about a floating body/joint. The floating body will be
         * ignored since it's not moved by a joint called out in the urdf. Only joints/bodies in 'joint'/'link'
         * tags will be used to populate the map.
         */
        bool parseJointBodyNameMapFromString(const char *model_xml_string, std::map<std::string,std::string> &jointBodyMap);

        /**
         * @brief This will build a map of joint name to body name.
         * @param filename Filepath
         * @param jointBodyMap Modified
         * @return
         *
         * @warning This will NOT give any information about a floating body/joint. The floating body will be
         * ignored since it's not moved by a joint called out in the urdf. Only joints/bodies in 'joint'/'link'
         * tags will be used to populate the map.
         */
        bool parseJointBodyNameMapFromFile(const char *filename, std::map<std::string,std::string> &jointBodyMap);
    }
}

/* ___RDL_URDFREADER_H__ */
#endif // ifndef __RDL_URDFREADER_H__
