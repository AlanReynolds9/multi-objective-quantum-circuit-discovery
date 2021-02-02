// MOQCD1 (Multi-objective quantum circuit discovery - algorithm 1.)
//
// Copyright (c) 2020  Ian Whittley, Alan Reynolds (alanreynolds9@gmail.com)
//
// This file is part of MOQCD1.
//
// MOQCD1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// MOQCD1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with MOQCD1.  If not, see
// <http://www.gnu.org/licenses/>.

//----------------------------------------------------------------------------------------------------------------------

#ifndef UTILITY_IO_CHECKED_STREAM_H
#define UTILITY_IO_CHECKED_STREAM_H

/*! \file checked_stream.h
    \brief Wrappers for ifstream and ofstream that call a user defined 
           error function if the file they are used to access cannot be
           opened.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

namespace utils
{
  namespace io
  {
    /*! \struct ThrowException
        \brief Error handling class for checked_stream.

        Throws an exception of type std::runtime_error if a file
        cannot be opened. The filepath is included as part of the
        exception string.
    */
    struct ThrowException
    {
      /*! Error handling function.
      */
      void handle_open_error(std::string const& filepath)
      {
        throw std::runtime_error(std::string("Error opening file: ") + filepath);
      }
    };

    /*! \struct CallExit
        \brief Error handling class for checked_stream.

        Sends an error report to std::cerr if a file cannot be opened
        and then calls exit(1). The filepath is included in the error
        message.
    */
    struct CallExit
    {
      /*! Error handling function.
      */
      void handle_open_error(std::string const& filepath)
      {
        std::cerr << "Error opening file: " << filepath << std::endl;
        exit(1);
      }
    };

    /*! \class checked
        \brief A simple file stream wrapper that is templated on the 
              stream to be wrapped and the error handler. The default
              error handler is CallExit - simply because we have typedefs
              for the most common usage that explicitly use ThrowException.

        A simple wrapper for file streams that calls a user defined error
        handler if the file it is called upon to access cannot be opened
        for any reason. Note that both the constructor and 'open' function
        take a std::string argument rather than the standard char*. 
    */
    template<class Stream, class ErrorHandler = CallExit>
    class checked : public Stream, public ErrorHandler
    {
    public:
      /// default constructor
      checked() : Stream(){}

      /// construct an open stream to the file at the path provided or calls the specified error handler.
      checked(std::string const& filepath) : Stream()
      {
        if(filepath.empty())
          ErrorHandler::handle_open_error("<Empty file path>");
          
        open(filepath);
        if(!Stream::is_open())
          ErrorHandler::handle_open_error(filepath);
      }
      
      // Construct a stream with explicit file mode.
      // (Repeating code, but no way to provide default value for the second argument)
      checked(std::string const& filepath, std::ios_base::openmode mode) : Stream()
      {
        if (filepath.empty())
          ErrorHandler::handle_open_error("<Empty file path>");
          
        open(filepath, mode);
        if(!Stream::is_open())
          ErrorHandler::handle_open_error(filepath);
      }

      /// opens a stream to the file at the path provided or calls the specified error handler.
      // should we have an additional flags parameter?
      void open(std::string const& filepath)
      {
        if(Stream::is_open())
        {
          Stream::close();
          Stream::clear();
        }

        if(filepath.empty())
          ErrorHandler::handle_open_error("<Empty file path>");

        Stream::open(filepath.c_str());
        if(!Stream::is_open())
          ErrorHandler::handle_open_error(filepath);
      }

      // Open a stream with explicit file mode.
      // (Repeating code, but no way to provie default value for the second argument)
      void open(std::string const& filepath, std::ios_base::openmode mode)
      {
        if(Stream::is_open())
        {
          Stream::close();
          Stream::clear();
        }

        if(filepath.empty())
          ErrorHandler::handle_open_error("<Empty file path>");

        Stream::open(filepath.c_str(), mode);
        if(!Stream::is_open())
          ErrorHandler::handle_open_error(filepath);
      }
    };       

    /*! \var typedef checked<std::ifstream> checked_ifstream
        \brief Canonical form for a wrapped input file stream.
      
        typedef for a checked std::ifstream using exceptions as the error reporting mechanism.
    */
    typedef checked<std::ifstream, ThrowException> checked_ifstream; 

    /*! \var typedef checked<std::ofstream> checked_ofstream
        \brief Canonical form for a wrapped output file stream.
      
        typedef for a checked std::ofstream using exceptions as the error reporting mechanism.
    */
    typedef checked<std::ofstream, ThrowException> checked_ofstream;
  }
}

#endif
