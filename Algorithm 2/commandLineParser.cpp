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

#include <vector>
#include <string>
#include <iostream>
#include "commandLineParser.h"

using namespace std;

CommandLineParser::CommandLineParser(const vector<string>& simpleFlags, const vector<string>& flagsWithValues) :
simple_flags(simpleFlags),
flags_with_values(flagsWithValues)
{
}

void CommandLineParser::readArguments(int argc, char* argv[])
{
  // Resest the array and index
  arguments_.clear();
  arg_index = 0;
  finished_ = false;

  // Read in the command line arguments
  for (int i = 1; i < argc; ++i)  // Skip the program name
  {
    arguments_.push_back(string(argv[i]));
  }

  // Ready the first flag (and value if there is one)
  next();
}

bool CommandLineParser::finished() const
{
  return finished_;
}

void CommandLineParser::next()
{
  // Function to ready the next flag and the associated value if there is one.
  // Also called by 'readArguments()' to get the first flag and value.

  // Find the first argument that begins with a dash
  string argument = " ";
  while (argument[0] != '-' && arg_index < static_cast<int>(arguments_.size()))
  {
    argument = arguments_[arg_index++];
  }

  if (argument[0] != '-' || argument == "-")
  {
    // Either no flag is found, or we have found a lone dash, which tells us to stop reading the command line.
    finished_ = true;
    return;
  }

  // Flag is 'argument' with the initial dash stripped away.
  flag_ = argument.substr(1);
  if (find(simple_flags.begin(), simple_flags.end(), flag_) != simple_flags.end())
  {
    // A simple flag
    value_ = "";
  }
  else if (find(flags_with_values.begin(), flags_with_values.end(), flag_) != flags_with_values.end())
  {
    // A flag with a value
    if (arg_index < static_cast<int>(arguments_.size()))
    {
      value_ = arguments_[arg_index++];
    }
    else
    {
      // Flag should have an associated value, but none can be found
      throw CommandLineError("Flag -" + flag_ + " should have a value, but none can be found.");
    }
  }
  else
  {
    // An unknown flag.
    throw CommandLineError("Unknown command line flag -" + flag_ + ".");
  }
}

string CommandLineParser::flag() const
{
  return flag_;
}

string CommandLineParser::value() const
{
  return value_;
}

