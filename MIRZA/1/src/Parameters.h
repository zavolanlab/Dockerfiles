// ===========================================================================
// Name        : Parameters.h
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

#include <string>
#include <math.h>
#include <vector>
#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#define A 0
#define C 1
#define G 2
#define T 3
#define U 3
#define N 4

#define NUMBER_OF_NUCLEOTIDES 4


namespace std {

class Parameters {

	typedef std::vector<double> double_vector;
	typedef std::vector<double_vector> double_2d_vector;

private:
	string mRNA;
	string miRNA;


public:

	double_vector E_hybrid; //Hybridization positions, where an Energy E_hybrid[i] is added for hybridizing position i of the miRNA
	double eE_GC;
	double eE_AT;
	double eE_GT;

	double E_sym;     //  Symmetric bases loop energy contribution
	double E_open;    //Energy E_open is added for every internal loop that is "OPENED"
	double E_miRNA_assymetric_loop;
	double  E_mRNA_assymetric_loop;

	double Wab[4][4]; // Base i from mRNA           hybridized to j from miRNA
	double Wa[4];     // Base i from mRNA is/is NOT hybridized to   the  miRNA (majinal probability)

	void Initialize();

	Parameters();
	virtual ~Parameters();
};

}

#endif /* PARAMETERS_H_ */
