/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_FRAME_EXCEPTIONS__
#define __RDL_FRAME_EXCEPTIONS__

#include <stdexcept>
#include <exception>

namespace RobotDynamics
{
    /**
     * @class ReferenceFrameException
     * @ingroup reference_frame
     * @brief A custom exception for frame operations
     */
    class ReferenceFrameException : public std::exception
    {
public:

        /**
         * @brief Constructor
         * @param err
         */
        ReferenceFrameException(const std::string& err) : msg(err)
        {}

        virtual const char* what() const throw()
        {
            return msg.c_str();
        }

        std::string msg /**< Custom exception message*/;
    };
}
#endif // ifndef __RDL_FRAME_EXCEPTIONS__
