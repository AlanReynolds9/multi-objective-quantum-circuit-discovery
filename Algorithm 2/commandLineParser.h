// MOQCD2 (Multi-objective quantum circuit discovery - algorithm 2.)
//
// Copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
//
// This file is part of MOQCD2.
//
// MOQCD2 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// MOQCD2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with MOQCD2.  If not, see
// <http://www.gnu.org/licenses/>.

//----------------------------------------------------------------------------------------------------------------------

#ifndef COMMAND_LINE_PARSER_H
#define COMMAND_LINE_PARSER_H

#include <stdexcept>
#include <vector>
#include <string>

class CommandLineError : public std::runtime_error
{
public:
  explicit CommandLineError(const std::string& msg) : runtime_error(msg) {}
};

class CommandLineParser
{
  // Class that takes the command line and breaks it into flags and values, which are returned to the user
public:
  CommandLineParser(const std::vector<std::string>& simpleFlags, const std::vector<std::string>& flagsWithValues);
  void readArguments(int argc, char* argv[]);

  bool finished() const;
  void next();
  std::string flag() const;
  std::string value() const;

private:
  std::vector<std::string> simple_flags;
  std::vector<std::string> flags_with_values;
  std::vector<std::string> arguments_;

  int arg_index;
  bool finished_;
  std::string flag_;
  std::string value_;
};

#endif
