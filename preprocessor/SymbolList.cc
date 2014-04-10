/*
 * Copyright (C) 2003-2014 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "SymbolList.hh"

void
SymbolList::addSymbol(const string &symbol)
{
  symbols.push_back(symbol);
}

void
SymbolList::writeOutput(const string &varname, ostream &output) const
{
  output << varname << " = char(";
  for (vector<string>::const_iterator it = symbols.begin();
       it != symbols.end(); ++it)
    {
      if (it != symbols.begin())
        output << ",";
      output << "'" << *it << "'";
    }
  output << ");" << endl;
}

void
SymbolList::clear()
{
  symbols.clear();
}
