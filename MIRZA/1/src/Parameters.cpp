// ===========================================================================
// Name        : Parameters.cpp
// Author      : Mohsen Khorshid
// Copyright   : University of Basel, 2010
// Description : Alignment model miRNA to mRNA target
// ===========================================================================
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "Parameters.h"
#include <string>
#include <math.h>
#include <vector>
#include <iostream>

namespace std {


void Parameters::Initialize() {
	//This block will set and initialize the Hybridization positional energies
	//TODO make sure the initialize values coming from Turner 2004 Nearest Neighbors DB

	//Initialize the Wa[Base_in_mRNA] with the base content of 3'UTRs from expressed representative RefSeqs
	Wa[A] = 0.277;
	Wa[C] = 0.209;
	Wa[G] = 0.214;
	Wa[T] = 0.300;

	//Wab[Base_in_mRNA][Base_in_miRNA]
	Wab[A][A] = 0.0;
	Wab[A][C] = 0.0;
	Wab[A][G] = 0.0;
	Wab[A][T] = 1.0;

	Wab[C][A] = 0.0;
	Wab[C][C] = 0.0;
	Wab[C][G] = 1.0;
	Wab[C][T] = 0.0;

	Wab[G][A] = 0.0;
	Wab[G][C] = 1.0;
	Wab[G][G] = 0.0;
	Wab[G][T] = 1.0;

	Wab[T][A] = 1.0;
	Wab[T][C] = 0.0;
	Wab[T][G] = 1.0;
	Wab[T][T] = 0.0;

}

Parameters::Parameters() {
	// TODO Auto-generated constructor stub
	this->Initialize();
}

Parameters::~Parameters() {
	// TODO Auto-generated destructor stub
}

}
