/***********************************************************************
    Copyright (c) 2000-2021,
    QuantumPlasma team, Max Planck Institut fur Kernphysik
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:

    * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.  

    * Redistributions in binary form must reproduce the above
        copyright notice, this list of conditions and the following
        disclaimer in the documentation and/or other materials provided
        with the distribution.  

    * Neither the name of the copyright holder nor the names of its
        contributors may be used to endorse or promote products derived
        from this software without specific prior written permission.

    * The software using and/or implementing this library must include
        a citation to [] in the official documentation.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
    OF THE POSSIBILITY OF SUCH DAMAGE.
*************************************************************************/

#include "SFQEDtoolkit_BaseInterface_Functions.h"

//inclusion needed to perform the file existence check
#include <sys/types.h>
#include <sys/stat.h>

//define environment variable name
// use string literals to initialize a const pointer to char
const char *variable_name = "SFQED_TOOLKIT_USER";

//make sure to pass an empty message string that will be
// filled with what the toolkit is doing
std::string find_path(std::string &message){
    //retrieve the path to the coefficients' location

    const char * val = std::getenv(variable_name);

    message = "";

    // invalid to assign nullptr to std::string
    if ( val == nullptr ) {
            //overwrite the message with informations about the location of 
            message = message + "Environment variable " + std::string(variable_name)
                    + " not found! Loading SFQED coefficients from root directory.\n";
               
            val = ".";
    }

    //build the final path to the coefficients
    std::string str_path = std::string(val) + std::string("/coefficients/");

    message = message + "Loading SFQED coefficients from " + str_path + ".\n";

    return str_path;
}


bool file_checker(std::string common_path, const char *phtn_file_names[], const int &entry_num, std::string &message){

    struct stat info;
    
    bool returner = true;

    for(int i = 0; i < entry_num; i++){

        const char *path_to_check = (std::string(common_path) + std::string(phtn_file_names[i])).c_str();
        
        if( stat(path_to_check, &info) != 0 || (info.st_mode & S_IFDIR) ){
            message = message + "File " + std::string(phtn_file_names[i]) + " not found!\n";
            returner = false;
        }

    }

    return returner;
}
