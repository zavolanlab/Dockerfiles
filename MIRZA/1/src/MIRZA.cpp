// ===========================================================================
// Name        : MIRZA.cpp
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

#include "Hybridize.h"
#include <stdlib.h>
#include <string>
#include <math.h>
#include <time.h>


#define miRNA_FIXED_LENGTH 21
#define NUMBER_OF_ENERGY_PARAMETER_LENGTH 6
#define NOISE 0.05
#define ENERGY_SCALE 20


using namespace std;

//The Hybridizer
Hybridize Hybridizer = Hybridize();

int main(int argc, char* argv[]) {

	int mRNAlength = atoi(argv[4]);
	int updatePriors = 0;

	string PriorUpdate_ARGUMENT = argv[5];
	if(!(PriorUpdate_ARGUMENT.compare("update"))){updatePriors = 1;}

	//Initialize the Hybridizer and allocate the necessary memory
	Hybridizer.mRNA_Length = mRNAlength;
	Hybridizer.updatePriors = updatePriors;
	Hybridizer.Initialize();


	//microRNA expression levels
	string miRNA_expression_file = argv[1];

	//set the file Names for Hybridizer
	string my_mRNA_file          = argv[2];
	string my_miRNA_file         = argv[3];


	Hybridizer.Read_mRNA_fasta_file(my_mRNA_file);

	//Estimate the base content fractions in the mRNA sequences
	Hybridizer.Params.Initialize();

	cerr << "__________________________________________________" << endl;
	Hybridizer.Read_miRNA_fasta_file(my_miRNA_file);
	cerr << "__________________________________________________" << endl;

	//Setting the initial miRNA expression levels in the sample
	Hybridizer.Read_miRNA_expressions(miRNA_expression_file);
	cerr << "__________________________________________________" << endl;

	Hybridizer.Gaussian_Penalty = 0;

	// These are the parameters of the hybridization model
	double parameters[miRNA_FIXED_LENGTH + NUMBER_OF_ENERGY_PARAMETER_LENGTH] = {0};

	cout << "Setting the Energy Parameters of the Hybridizer with the following values: " << endl;
	parameters[0]= 2.38714;
	parameters[1]= 2.39491;
	parameters[2]= -3.61813;
	parameters[3]= -0.285659;
	parameters[4]= 1.01614;
	parameters[5]= 1.03229;

	//These are the 21 positional parameters
	parameters[6]= -3.24074;
	parameters[7]= -0.249397;
	parameters[8]= 3.45158;
	parameters[9]= 3.59005;
	parameters[10]= 0.609733;
	parameters[11]= 2.46235;
	parameters[12]= -0.0431051;
	parameters[13]= -0.815238;
	parameters[14]= -4.25303;
	parameters[15]= -3.08829;
	parameters[16]= -1.84676;
	parameters[17]= -4.55569;
	parameters[18]= -1.75914;
	parameters[19]= -1.53761;
	parameters[20]= -1.78762;
	parameters[21]= -1.46246;
	parameters[22]= -4.20649;
	parameters[23]= -1.79764;
	parameters[24]= -2.24505;
	parameters[25]= -4.15307;
	parameters[26]= -3.82612;

	for (unsigned int i = 0; i < miRNA_FIXED_LENGTH	+ NUMBER_OF_ENERGY_PARAMETER_LENGTH; i++) {
		cout << "init[" << i << "] = " << parameters[i] << ";" << endl;

		Hybridizer.Gaussian_Penalty += (parameters[i] / ENERGY_SCALE) * (parameters[i] / ENERGY_SCALE);
	}


	//Initialize the looping energies
	Hybridizer.eE_OPEN=0;
	Hybridizer.eE_SYM=0;
	Hybridizer.eE_mRNA_Assymetric_loops=0;
	Hybridizer.eE_miRNA_Assymetric_loops=0;

	//Initializing the miRNA priors relative to their expression in the cell
	Hybridizer.Initialize_miRNA_priors();

	//Initialize mRNA versus miRNA likelihood ratios
	Hybridizer.Initialize_miRNA_mRNA_likelihood_ratios();

	//Initialaize the values mRNA likelihood ratios
	Hybridizer.Initialize_mRNA_log_likelihood_ratios();

	//Set the energy parameters [with length of ENERGY_PARAMETER_LENGTH] that are being optimized and Calculate the Exponentials and Coefficients
	Hybridizer.Initialize_Global_Parameters_and_Prepare(parameters);

	//Start the calculation
	Hybridizer.Calculate_data_log_likelihood_ratio();


	return 0;
}
