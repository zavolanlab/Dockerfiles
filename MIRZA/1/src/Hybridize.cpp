// ===========================================================================
// Name        : Hybridize.cpp
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

using namespace std;

#include "Hybridize.h"
#include "Parameters.h"
#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>

Hybridize::Hybridize() {

	eE_OPEN = ZERO;
	eE_SYM = ZERO;

	eE_mRNA_Assymetric_loops = ZERO;
	eE_miRNA_Assymetric_loops = ZERO;

	BEST_HYBRID_i= ZERO;
	BEST_HYBRID_j= ZERO;
	BEST_HYBRID_SCORE = STOP;

	Gaussian_Penalty = ZERO;

	F_alpha= NULL;
	H_alpha=NULL;
	Best_Hybrid=NULL;
	Trace_Back=NULL;

	Initialize();
}

void Hybridize::Initialize_miRNA_priors() {

	memset(miRNA_Priors,ZERO,sizeof(miRNA_Priors));

	total_miRNA_expressions = 0; //This is for normalizing the miRNA expression
	for (unsigned int miRNA_id = 0; miRNA_id < miRNAs.Get_Size(); miRNA_id++) {
		total_miRNA_expressions += miRNA_Expression_Hash[miRNAs.GetSequence(miRNA_id).header];
	}


	//This is for normalizing the miRNA expression
	for (unsigned int miRNA_id = 0; miRNA_id < miRNAs.Get_Size(); miRNA_id++) {
		miRNA_Priors[miRNA_id] = miRNA_Expression_Hash[miRNAs.GetSequence(miRNA_id).header] / total_miRNA_expressions;
	}

	return;
}

void Hybridize::Initialize_mRNA_log_likelihood_ratios(){

	memset(mRNA_likelihood_ratios,ZERO,sizeof(mRNA_likelihood_ratios));
	return;
}

void Hybridize::Initialize_miRNA_mRNA_likelihood_ratios(){
	memset(Ratio,ZERO,sizeof(Ratio));
	return;
}


void Hybridize::Initialize_Global_Parameters_and_Prepare(double myParams[]){
        cerr << "Initializing the coefficients and energy arrays..." << endl;
        //Set the energy parameters [with length of ENERGY_PARAMETER_LENGTH] that are being optimized
        Params.eE_AT = exp(myParams[E_AT]);
        Params.eE_GC = exp(myParams[E_GC]);
        Params.eE_GT = 1; //exp(0);


        double eE_alpha[NUMBER_OF_NUCLEOTIDES] = { 0 };

        eE_alpha[A] = Params.Wa[T] * Params.eE_AT;
        eE_alpha[C] = Params.Wa[G] * Params.eE_GC;
        eE_alpha[G] = Params.Wa[C] * Params.eE_GC + Params.Wa[T] * Params.eE_GT;
        eE_alpha[T] = Params.Wa[A] * Params.eE_AT + Params.Wa[G] * Params.eE_GT;

        Params.Wab[T][G] = Params.Wa[T] * Params.eE_GT / eE_alpha[G];
        Params.Wab[C][G] = Params.Wa[C] * Params.eE_GC / eE_alpha[G];

        //Params.Wab[G][C] = 1;
        //Params.Wab[T][A] = 1;

        Params.Wab[A][T] = Params.Wa[A] * Params.eE_AT / eE_alpha[T];
        Params.Wab[G][T] = Params.Wa[G] * Params.eE_GT / eE_alpha[T];

        
        //Calculate the Exponentials
        Params.E_open = myParams[E_OPEN];
        Params.E_mRNA_assymetric_loop = myParams[E_mRNA_ASSYMETRIC_LOOP];
        Params.E_miRNA_assymetric_loop = myParams[E_miRNA_ASSYMETRIC_LOOP];
        Params.E_sym = myParams[E_SYMMETRIC_LOOP];

        eE_SYM = exp(Params.E_sym);
        eE_OPEN = exp(Params.E_open);
        eE_mRNA_Assymetric_loops = exp(Params.E_mRNA_assymetric_loop);
        eE_miRNA_Assymetric_loops= exp(Params.E_miRNA_assymetric_loop);

        //this block calculates the exponentials of
        // eE_hybridization =e^(E_hybrid[i]+E_alpha[A,C,G,T])

        unsigned  Number_of_miRNAs = miRNAs.Get_Size();


        memset(eE_hybridization, NOTHING, sizeof(eE_hybridization));

        string miRNA_sequence="";
        unsigned miRNA_Sequence_Length =0;


        for (unsigned int miRNA_id = 0; miRNA_id < Number_of_miRNAs; miRNA_id++) {
                miRNA_sequence = miRNAs.GetSequence(miRNA_id).sequence;
                miRNA_Sequence_Length = miRNAs.GetSequence(miRNA_id).sequence.length();

                for (unsigned int i = 0; i      < miRNA_Sequence_Length; i++) {
                        eE_hybridization[miRNA_id][i] = exp(
                                        myParams[miRNA_HYBRIDIZATION_ENERGY_PARAMETER_OFFSET + i] +
                                        eE_alpha[convert(miRNA_sequence[i])]);
                }
        }
	

        memset(Coefficients, 0, sizeof(Coefficients));
        //Coefficients[A][A] = 0.0;
        //Coefficients[A][C] = 0.0;
        //Coefficients[A][G] = 0.0;
        Coefficients[A][T] = Params.Wab[T][A] / Params.Wa[T];

        //Coefficients[C][A] = 0.0;
        //Coefficients[C][C] = 0.0;
        Coefficients[C][G] = Params.Wab[G][C] / Params.Wa[G];
        //Coefficients[C][T] = 0.0;

        //Coefficients[G][A] = 0.0;
        Coefficients[G][C] = Params.Wab[C][G] / Params.Wa[C];
        //Coefficients[G][G] = 0.0;
        Coefficients[G][T] = Params.Wab[T][G] / Params.Wa[T];

        Coefficients[T][A] = Params.Wab[A][T] / Params.Wa[A];
        //Coefficients[T][C] = 0.0;
        Coefficients[T][G] = Params.Wab[G][T] / Params.Wa[G];
        //Coefficients[T][T] = 0.0;

	cerr << "Initializing the coefficients and energy arrays...Done" << endl;
	return;
}

double Hybridize::Calculate_data_log_likelihood_ratio() {

	unsigned  int  Number_of_mRNAs =  mRNAs.Get_Size();
	unsigned  int Number_of_miRNAs = miRNAs.Get_Size();
	//Calculate
	for (unsigned int mRNA_id = 0; mRNA_id < Number_of_mRNAs; mRNA_id++) {
		mRNA_likelihood_ratios[mRNA_id] =0;
		for (unsigned int miRNA_id = 0; miRNA_id< Number_of_miRNAs; miRNA_id++) {
			resetMemory();
			Ratio[miRNA_id][mRNA_id] = CalculateHibridizationRatio_miRNA_mRNA(miRNA_id,mRNA_id);
			mRNA_likelihood_ratios[mRNA_id] += Ratio[miRNA_id][mRNA_id] * miRNA_Priors[miRNA_id];
		}
	}

	//Perfrom Expectation Maximization for updating the Priors and return back the Maximum Likelihood Ratio of the Data
	//Update miRNA priors and mRNA likelihood ratios based on the updated priors
	
	if (updatePriors) {
		Update_miRNA_priors();
	}

	for (unsigned int mRNA_id = 0; mRNA_id < Number_of_mRNAs; mRNA_id++) {
		for (unsigned int miRNA_id = 0; miRNA_id < Number_of_miRNAs; miRNA_id++) {
			cout << "+++++++++++++++++++++++" <<endl;
			cout << "likelihood ratio of mRNA ("
				<< mRNAs.GetSequence(mRNA_id).header
				<< ") "
				<< mRNAs.GetSequence(mRNA_id).sequence
				<< " and miRNA ("
				<< miRNAs.GetSequence(miRNA_id).header
				<< ") "
				<< miRNAs.GetSequence(miRNA_id).sequence
				<< " = "
				<< Ratio[miRNA_id][mRNA_id]
				<< "\tRatio*miRNA_Prior= "
				<< Ratio[miRNA_id][mRNA_id]*miRNA_Priors[miRNA_id]
				<< "\tmiRNA_Prior= "
				<< miRNA_Priors[miRNA_id]
				<< endl;
			cout << "+++++++++++++++++++++++" << endl;
		}
	}

	double data_log_likelihood_ratio = 0;
	for (unsigned int mRNA_id = 0; mRNA_id < Number_of_mRNAs; mRNA_id++) {
		data_log_likelihood_ratio += log(mRNA_likelihood_ratios[mRNA_id]);
	}
	cout << "Gaussian Penalty of the parameters = " << Gaussian_Penalty << endl;
	cout << "Log Likelihood Ratio of the data = " << data_log_likelihood_ratio << endl;

	data_log_likelihood_ratio -= Gaussian_Penalty;

	double return_value = -1 * data_log_likelihood_ratio;
	cout << "Return Value =" << return_value << endl;

	return return_value;
}

void Hybridize::Initialize() {

	F_alpha = new double**[miRNA_LENGTH];
	H_alpha = new double**[miRNA_LENGTH];
	Best_Hybrid = new double**[miRNA_LENGTH];
	Trace_Back = new int***[miRNA_LENGTH];
	for (int x = 0; x < miRNA_LENGTH; x++) {
		F_alpha[x] = new double*[mRNA_Length];
		H_alpha[x] = new double*[mRNA_Length];
		Best_Hybrid[x] = new double*[mRNA_Length];
		Trace_Back[x] = new int**[mRNA_Length];
		for (int y = 0; y < mRNA_Length; y++) {
			F_alpha[x][y] = new double[NUMBER_OF_HYBRIDIZATION_MODES];
			H_alpha[x][y] = new double[NUMBER_OF_HYBRIDIZATION_MODES];
			Best_Hybrid[x][y] = new double[NUMBER_OF_HYBRIDIZATION_MODES];
			Trace_Back[x][y] = new int*[NUMBER_OF_HYBRIDIZATION_MODES];
			for (int z = 0; z < NUMBER_OF_HYBRIDIZATION_MODES; z++) {
				Trace_Back[x][y][z] = new int[NUMBER_OF_TRACE_BACK_ELEMENTS];
			}
		}
	}
}

void Hybridize::resetMemory() {

	for (int x = 0; x < miRNA_LENGTH; x++) {
		for (int y = 0; y < mRNA_Length; y++) {
			memset(F_alpha[x][y], ZERO, sizeof(double) * NUMBER_OF_HYBRIDIZATION_MODES);
			memset(H_alpha[x][y], ZERO, sizeof(double) * NUMBER_OF_HYBRIDIZATION_MODES);
			memset(Best_Hybrid[x][y], ZERO, sizeof(double) * NUMBER_OF_HYBRIDIZATION_MODES);
			for (int z = 0; z < NUMBER_OF_HYBRIDIZATION_MODES; z++) {
				memset(Trace_Back[x][y][z], ZERO, sizeof(int) * NUMBER_OF_TRACE_BACK_ELEMENTS);
			}
		}
	}

	BEST_HYBRID_i = ZERO;
	BEST_HYBRID_j = ZERO;

	BEST_HYBRID_SCORE = STOP;

	// Set Z0 and Z1 to ZERO
	Z0= ZERO;
	H= ZERO;
	return;
}

double Hybridize::CalculateHibridizationRatio_miRNA_mRNA(unsigned int miRNA_id, unsigned int mRNA_id) {

	//This block will calculate the F_alpha[i][j][k] and H_alpha[i][j][k] where k = HYBRIDIZED_PAIR or UNMATCHED_SYMMETRIC_LOOP or UNMATCHED_ASYMMETRIC_mRNA or UNMATCHED_ASYMMETRIC_miRNA
	string  mRNA_sequence =  mRNAs.GetSequence( mRNA_id).reversed_sequence;
        string miRNA_sequence = miRNAs.GetSequence(miRNA_id).sequence;
        unsigned int miRNA_Seqeunce_length = miRNA_sequence.length();
        unsigned int  mRNA_Seqeunce_length =  mRNA_Length;


	for (unsigned int i = 0; i < miRNA_Seqeunce_length; i++) {
		int miRNA_base = convert(miRNA_sequence[i]);

		for (unsigned int j = 0; j < mRNA_Seqeunce_length; j++) {
			int mRNA_base =  convert(mRNA_sequence[j]);

			fill_the_F_and_H_matrices(miRNA_id, mRNA_id, Coefficients[miRNA_base][mRNA_base], i, j);
			fill_the_Best_Hybrid_matrix(miRNA_id, mRNA_id, Coefficients[miRNA_base][mRNA_base], i, j);

			//Update the Best Hybrid structure
			if(BEST_HYBRID_SCORE <= Best_Hybrid[i][j][HYBRIDIZED_PAIR]){
				BEST_HYBRID_i= i;
				BEST_HYBRID_j= j;
				BEST_HYBRID_SCORE = Best_Hybrid[i][j][HYBRIDIZED_PAIR];
			}

			Z0 += F_alpha[i][j][HYBRIDIZED_PAIR];
			 H += H_alpha[i][j][HYBRIDIZED_PAIR];
		}
	}

	double ratio = H/Z0;

	//Trace back and print out the best hybrid structure
	Trace_Back_and_Print_Hybridization(miRNA_id, mRNA_id, ratio);
	return ratio;
}

void Hybridize::fill_the_F_and_H_matrices(unsigned int miRNA_id, unsigned int mRNA_id, double coefficient, int i, int j) {

	if( i!=0 && j!=0){ //Standard case hybridizing base i+1th from miRNA with base j+1th of the mRNA
		F_alpha[i][j][HYBRIDIZED_PAIR] =               eE_hybridization[miRNA_id][i] * ( ONE + F_alpha[i - 1][j - 1][HYBRIDIZED_PAIR] + F_alpha[i - 1][j - 1][UNMATCHED_SYMMETRIC_LOOP] + F_alpha[i - 1][j - 1][UNMATCHED_ASYMMETRIC_mRNA] + F_alpha[i - 1][j - 1][UNMATCHED_ASYMMETRIC_miRNA]);
		H_alpha[i][j][HYBRIDIZED_PAIR] = coefficient * eE_hybridization[miRNA_id][i] * ( ONE + H_alpha[i - 1][j - 1][HYBRIDIZED_PAIR] + H_alpha[i - 1][j - 1][UNMATCHED_SYMMETRIC_LOOP] + H_alpha[i - 1][j - 1][UNMATCHED_ASYMMETRIC_mRNA] + H_alpha[i - 1][j - 1][UNMATCHED_ASYMMETRIC_miRNA]);

		F_alpha[i][j][UNMATCHED_SYMMETRIC_LOOP] = eE_SYM * (eE_OPEN * F_alpha[i - 1][j - 1][HYBRIDIZED_PAIR] + F_alpha[i - 1][j - 1][UNMATCHED_SYMMETRIC_LOOP]);
		H_alpha[i][j][UNMATCHED_SYMMETRIC_LOOP] = eE_SYM * (eE_OPEN * H_alpha[i - 1][j - 1][HYBRIDIZED_PAIR] + H_alpha[i - 1][j - 1][UNMATCHED_SYMMETRIC_LOOP]);

		F_alpha[i][j][UNMATCHED_ASYMMETRIC_mRNA] =   eE_mRNA_Assymetric_loops * (eE_OPEN * F_alpha[i][j - 1][HYBRIDIZED_PAIR] + F_alpha[i][j - 1][UNMATCHED_SYMMETRIC_LOOP] + F_alpha[i][j - 1][UNMATCHED_ASYMMETRIC_mRNA]);
		H_alpha[i][j][UNMATCHED_ASYMMETRIC_mRNA] =   eE_mRNA_Assymetric_loops * (eE_OPEN * H_alpha[i][j - 1][HYBRIDIZED_PAIR] + H_alpha[i][j - 1][UNMATCHED_SYMMETRIC_LOOP] + H_alpha[i][j - 1][UNMATCHED_ASYMMETRIC_mRNA]);

		F_alpha[i][j][UNMATCHED_ASYMMETRIC_miRNA] = eE_miRNA_Assymetric_loops * (eE_OPEN * F_alpha[i - 1][j][HYBRIDIZED_PAIR] + F_alpha[i - 1][j][UNMATCHED_SYMMETRIC_LOOP] + F_alpha[i - 1][j][UNMATCHED_ASYMMETRIC_miRNA]);
		H_alpha[i][j][UNMATCHED_ASYMMETRIC_miRNA] = eE_miRNA_Assymetric_loops * (eE_OPEN * H_alpha[i - 1][j][HYBRIDIZED_PAIR] + H_alpha[i - 1][j][UNMATCHED_SYMMETRIC_LOOP] + H_alpha[i - 1][j][UNMATCHED_ASYMMETRIC_miRNA]);

		return;
	}

	else if (i == 0 && j != 0) { //The boundary case when we start hybridization the first base of mRNA to somewhere in the middle of miRNA (some first bases of miRNA gets skipped)
		F_alpha[i][j][HYBRIDIZED_PAIR] =                eE_hybridization[miRNA_id][i];
		H_alpha[i][j][HYBRIDIZED_PAIR] =  coefficient * eE_hybridization[miRNA_id][i];


		F_alpha[i][j][UNMATCHED_SYMMETRIC_LOOP] = 0;
		H_alpha[i][j][UNMATCHED_SYMMETRIC_LOOP] = 0;

		F_alpha[i][j][UNMATCHED_ASYMMETRIC_mRNA] = eE_mRNA_Assymetric_loops * (eE_OPEN * F_alpha[i][j - 1][HYBRIDIZED_PAIR] + F_alpha[i][j - 1][UNMATCHED_SYMMETRIC_LOOP] + F_alpha[i][j - 1][UNMATCHED_ASYMMETRIC_mRNA]);
		H_alpha[i][j][UNMATCHED_ASYMMETRIC_mRNA] = eE_mRNA_Assymetric_loops * (eE_OPEN * H_alpha[i][j - 1][HYBRIDIZED_PAIR] + H_alpha[i][j - 1][UNMATCHED_SYMMETRIC_LOOP] + H_alpha[i][j - 1][UNMATCHED_ASYMMETRIC_mRNA]);

		F_alpha[i][j][UNMATCHED_ASYMMETRIC_miRNA] = 0;
		H_alpha[i][j][UNMATCHED_ASYMMETRIC_miRNA] = 0;


		return;
	}

	else if (j == 0 && i != 0) { //The boundary case when we start hybridization the first base of miRNA to somewhere in the middle of mRNA (some first bases of mRNA gets skipped)
		F_alpha[i][j][HYBRIDIZED_PAIR] =               eE_hybridization[miRNA_id][i];
		H_alpha[i][j][HYBRIDIZED_PAIR] = coefficient * eE_hybridization[miRNA_id][i];

		F_alpha[i][j][UNMATCHED_SYMMETRIC_LOOP] = 0;
		H_alpha[i][j][UNMATCHED_SYMMETRIC_LOOP] = 0;

		F_alpha[i][j][UNMATCHED_ASYMMETRIC_mRNA] = 0;
		H_alpha[i][j][UNMATCHED_ASYMMETRIC_mRNA] = 0;

		F_alpha[i][j][UNMATCHED_ASYMMETRIC_miRNA] = eE_miRNA_Assymetric_loops * (eE_OPEN * F_alpha[i - 1][j][HYBRIDIZED_PAIR] + F_alpha[i - 1][j][UNMATCHED_SYMMETRIC_LOOP] + F_alpha[i - 1][j][UNMATCHED_ASYMMETRIC_miRNA]);
		H_alpha[i][j][UNMATCHED_ASYMMETRIC_miRNA] = eE_miRNA_Assymetric_loops * (eE_OPEN * H_alpha[i - 1][j][HYBRIDIZED_PAIR] + H_alpha[i - 1][j][UNMATCHED_SYMMETRIC_LOOP] + H_alpha[i - 1][j][UNMATCHED_ASYMMETRIC_miRNA]);

		return;
	}

	else{ //The boundry case when we start hybridizing the first base of miRNA with the first base of mRNA

		F_alpha[i][j][HYBRIDIZED_PAIR] =               eE_hybridization[miRNA_id][i];
		H_alpha[i][j][HYBRIDIZED_PAIR] = coefficient * eE_hybridization[miRNA_id][i];

		F_alpha[i][j][UNMATCHED_SYMMETRIC_LOOP] = 0;
		H_alpha[i][j][UNMATCHED_SYMMETRIC_LOOP] = 0;

		F_alpha[i][j][UNMATCHED_ASYMMETRIC_mRNA] = 0;
		H_alpha[i][j][UNMATCHED_ASYMMETRIC_mRNA] = 0;

		F_alpha[i][j][UNMATCHED_ASYMMETRIC_miRNA] = 0;
		H_alpha[i][j][UNMATCHED_ASYMMETRIC_miRNA] = 0;

		return;
	}
}

void Hybridize::fill_the_Best_Hybrid_matrix(unsigned int miRNA_id, unsigned int mRNA_id, double coefficient, int i, int j){

	double MAXIMUM_for_HYBRIDIZED_PAIR = 0;
	double MAXIMUM_for_UNMATCHED_SYMMETRIC_LOOP = 0;
	double MAXIMUM_for_UNMATCHED_ASYMMETRIC_mRNA = 0;
	double MAXIMUM_for_UNMATCHED_ASYMMETRIC_miRNA = 0;

	int state;
	if( i!=0 && j!=0){
		//Best Hybrid and trace back for HYBRIDIZED_PAIR state
		MAXIMUM_for_HYBRIDIZED_PAIR = Get_Max_Hybrid_state(
				ONE,
				Best_Hybrid[i - 1][j - 1][HYBRIDIZED_PAIR],
				Best_Hybrid[i - 1][j - 1][UNMATCHED_SYMMETRIC_LOOP],
				Best_Hybrid[i - 1][j - 1][UNMATCHED_ASYMMETRIC_mRNA],
				Best_Hybrid[i - 1][j - 1][UNMATCHED_ASYMMETRIC_miRNA],
				state);
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_i]= i - 1;
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_j]= j - 1;
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_STATE] = state;


		//Best Hybrid and trace back for UNMATCHED_SYMMETRIC_LOOP state
		MAXIMUM_for_UNMATCHED_SYMMETRIC_LOOP = Get_Max_symmetric_loop_state(eE_OPEN *
				Best_Hybrid[i - 1][j - 1][HYBRIDIZED_PAIR],
				Best_Hybrid[i - 1][j - 1][UNMATCHED_SYMMETRIC_LOOP],
				state);
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_i]= i - 1;
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_j]= j - 1;
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_STATE] = state;

		//Best Hybrid and trace back for UNMATCHED_ASYMMETRIC_mRNA state
		MAXIMUM_for_UNMATCHED_ASYMMETRIC_mRNA = Get_Max_assymetric_mRNA_loop_state(eE_OPEN *
				Best_Hybrid[i][j - 1][HYBRIDIZED_PAIR],
				Best_Hybrid[i][j - 1][UNMATCHED_SYMMETRIC_LOOP],
				Best_Hybrid[i][j - 1][UNMATCHED_ASYMMETRIC_mRNA],
				state);
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_i]= i;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_j]= j - 1;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_STATE]= state;

		//Best Hybrid and trace back for UNMATCHED_ASYMMETRIC_miRNA state
		MAXIMUM_for_UNMATCHED_ASYMMETRIC_miRNA = Get_Max_assymetric_miRNA_loop_state(eE_OPEN *
				Best_Hybrid[i - 1][j][HYBRIDIZED_PAIR],
				Best_Hybrid[i - 1][j][UNMATCHED_SYMMETRIC_LOOP],
				Best_Hybrid[i - 1][j][UNMATCHED_ASYMMETRIC_miRNA],
				state);
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_i]= i - 1;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_j]= j;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_STATE]= state;
	}

	else if (i == 0 && j != 0) {
		//Best Hybrid and trace back for HYBRIDIZED_PAIR state
		MAXIMUM_for_HYBRIDIZED_PAIR = ONE;
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_i]= i;
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_j]= j - 1;
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_STATE] = BEGIN_OF_THE_HYBRID;

		//Best Hybrid and trace back for UNMATCHED_SYMMETRIC_LOOP state
		MAXIMUM_for_UNMATCHED_SYMMETRIC_LOOP = 0;
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_i]= i;
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_j]= j - 1;
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_STATE] = BEGIN_OF_THE_HYBRID;

		//Best Hybrid and trace back for UNMATCHED_ASYMMETRIC_miRNA state
		MAXIMUM_for_UNMATCHED_ASYMMETRIC_mRNA = Get_Max_assymetric_mRNA_loop_state(eE_OPEN *
				Best_Hybrid[i][j - 1][HYBRIDIZED_PAIR],
				Best_Hybrid[i][j - 1][UNMATCHED_SYMMETRIC_LOOP],
				Best_Hybrid[i][j - 1][UNMATCHED_ASYMMETRIC_mRNA],
				state);
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_i]= i;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_j]= j - 1;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_STATE] = state;

		//Best Hybrid and trace back for UNMATCHED_ASYMMETRIC_miRNA state
		MAXIMUM_for_UNMATCHED_ASYMMETRIC_miRNA = 0;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_i]= i;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_j]= j - 1;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_STATE] = BEGIN_OF_THE_HYBRID;
	}

	else if (i != 0 && j == 0) {
		//Best Hybrid and trace back for HYBRIDIZED_PAIR state
		MAXIMUM_for_HYBRIDIZED_PAIR = ONE;
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_i]= i - 1;
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_j]= j;
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_STATE] = BEGIN_OF_THE_HYBRID;

		//Best Hybrid and trace back for UNMATCHED_SYMMETRIC_LOOP state
		MAXIMUM_for_UNMATCHED_SYMMETRIC_LOOP = 0;
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_i]= i - 1;
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_j]= j;
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_STATE] = BEGIN_OF_THE_HYBRID;

		//Best Hybrid and trace back for UNMATCHED_ASYMMETRIC_mRNA state
		MAXIMUM_for_UNMATCHED_ASYMMETRIC_mRNA = 0;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_i]= i - 1;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_j]= j;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_STATE] = BEGIN_OF_THE_HYBRID;

		//Best Hybrid and trace back for UNMATCHED_ASYMMETRIC_miRNA state
		MAXIMUM_for_UNMATCHED_ASYMMETRIC_miRNA = Get_Max_assymetric_miRNA_loop_state(eE_OPEN *
				Best_Hybrid[i - 1][j][HYBRIDIZED_PAIR],
				Best_Hybrid[i - 1][j][UNMATCHED_SYMMETRIC_LOOP],
				Best_Hybrid[i - 1][j][UNMATCHED_ASYMMETRIC_miRNA],
				state);
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_i]= i - 1;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_j]= j;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_STATE] = state;
	}

	else{
		//Best Hybrid and trace back for HYBRIDIZED_PAIR state
		MAXIMUM_for_HYBRIDIZED_PAIR = ONE;
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_i]= i;
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_j]= j;
		Trace_Back[i][j][HYBRIDIZED_PAIR][TRACE_BACK_STATE] = BEGIN_OF_THE_HYBRID;

		//Best Hybrid and trace back for UNMATCHED_SYMMETRIC_LOOP state
		MAXIMUM_for_UNMATCHED_SYMMETRIC_LOOP = 0;
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_i]= i;
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_j]= j;
		Trace_Back[i][j][UNMATCHED_SYMMETRIC_LOOP][TRACE_BACK_STATE] = BEGIN_OF_THE_HYBRID;

		//Best Hybrid and trace back for UNMATCHED_ASYMMETRIC_mRNA state
		MAXIMUM_for_UNMATCHED_ASYMMETRIC_mRNA = 0;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_i]= i;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_j]= j;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_mRNA][TRACE_BACK_STATE] = BEGIN_OF_THE_HYBRID;

		//Best Hybrid and trace back for UNMATCHED_ASYMMETRIC_miRNA state
		MAXIMUM_for_UNMATCHED_ASYMMETRIC_miRNA = 0;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_i]= i;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_j]= j;
		Trace_Back[i][j][UNMATCHED_ASYMMETRIC_miRNA][TRACE_BACK_STATE] = BEGIN_OF_THE_HYBRID;

	}


	Best_Hybrid[i][j][HYBRIDIZED_PAIR]           = coefficient * eE_hybridization[miRNA_id][i] * MAXIMUM_for_HYBRIDIZED_PAIR;
	Best_Hybrid[i][j][UNMATCHED_SYMMETRIC_LOOP]  = eE_SYM                                      * MAXIMUM_for_UNMATCHED_SYMMETRIC_LOOP;
	Best_Hybrid[i][j][UNMATCHED_ASYMMETRIC_mRNA] = eE_mRNA_Assymetric_loops                    * MAXIMUM_for_UNMATCHED_ASYMMETRIC_mRNA;
	Best_Hybrid[i][j][UNMATCHED_ASYMMETRIC_miRNA]= eE_miRNA_Assymetric_loops                   * MAXIMUM_for_UNMATCHED_ASYMMETRIC_miRNA;
	return;
}

void Hybridize::Update_miRNA_priors(void){

	cerr << "Updating the priors and likelihood ratios" << endl;
	std::vector<double> miRNA_updated_priors;
	std::vector<double> mRNA_updated_likelihood_ratios;

	unsigned int Number_of_mRNAs=mRNAs.Get_Size();
	unsigned int Number_of_miRNAs=miRNAs.Get_Size();


	//Allocate memory and set to zero
	miRNA_updated_priors.resize(Number_of_miRNAs);
	miRNA_updated_priors.assign(Number_of_miRNAs, ZERO);
	mRNA_updated_likelihood_ratios.resize(Number_of_miRNAs);
	mRNA_updated_likelihood_ratios.assign(Number_of_mRNAs, ZERO);


	for (unsigned int miRNA_id = 0; miRNA_id < Number_of_miRNAs; miRNA_id++) {
		cout << "miRNA prior BEFORE\t" << miRNAs.GetSequence(miRNA_id).header << "\t" << miRNA_Priors[miRNA_id] << endl;
	}


	double delta=BIG_INTEGER;
	while (delta > ACCURACY) {
		// while the Relative mean difference >0.0001 then tru to update the miRNA priors
		delta = 0;
		for (unsigned int miRNA_id = 0; miRNA_id < Number_of_miRNAs; miRNA_id++) {
			miRNA_updated_priors[miRNA_id] = 0;
			for (unsigned int mRNA_id = 0; mRNA_id < Number_of_mRNAs; mRNA_id++) {
					miRNA_updated_priors[miRNA_id] += miRNA_Priors[miRNA_id] * (Ratio[miRNA_id][mRNA_id] / mRNA_likelihood_ratios[mRNA_id]);
			}
		}

		double normalization_factor = 0;
		for (unsigned int miRNA_id = 0; miRNA_id < Number_of_miRNAs; miRNA_id++) {
			normalization_factor += miRNA_updated_priors[miRNA_id];
		}

		//normalizing the updated priors and Updating the miRNA_Priors array
		for (unsigned int miRNA_id = 0; miRNA_id < Number_of_miRNAs; miRNA_id++) {

			//Normalize the new priors
			double x = miRNA_updated_priors[miRNA_id] / normalization_factor;

			//update the relative mean difference
			delta += (x - miRNA_Priors[miRNA_id]) * (x - miRNA_Priors[miRNA_id]);

			//update the priors
			miRNA_Priors[miRNA_id] = x;
		}
		//finish up the Euclidean distance = sqrt( sum((xi - xj)^2))
		delta =sqrt(delta);

		//update the mRNA likelihood ratio values
		for (unsigned int mRNA_id = 0; mRNA_id < Number_of_mRNAs; mRNA_id++) {
			mRNA_updated_likelihood_ratios[mRNA_id]=0;
			for (unsigned int miRNA_id = 0; miRNA_id < Number_of_miRNAs; miRNA_id++) {
				mRNA_updated_likelihood_ratios[mRNA_id]	+= Ratio[miRNA_id][mRNA_id] * miRNA_Priors[miRNA_id];
			}
			mRNA_likelihood_ratios[mRNA_id]	= mRNA_updated_likelihood_ratios[mRNA_id];
		}
	}

	for (unsigned int miRNA_id = 0; miRNA_id < Number_of_miRNAs; miRNA_id++) {
		cout << "miRNA prior AFTER\t" << miRNAs.GetSequence(miRNA_id).header << "\t" << miRNA_Priors[miRNA_id] << endl;
	}

	cerr << "Updating the priors and likelihood ratios...Done!" << endl;
	return;
}

void Hybridize::Trace_Back_and_Print_Hybridization(unsigned int miRNA_id,unsigned int mRNA_id, double Ratio){

	string mRNA_sequence = mRNAs.GetSequence(mRNA_id).reversed_sequence;
	string miRNA_sequence = miRNAs.GetSequence(miRNA_id).sequence;
	string miRNA_alignment = "";
	string mRNA_alignment = "";
	string alignment = "";
	int miRNA_length = miRNA_sequence.length() - 1;
	int mRNA_length = mRNA_sequence.length() - 1;
	int i = miRNA_length;
	int j = mRNA_length;
	int current_MODE;


	//print the dangling ends till the begining of the best end_position od the best hybrid
	int max = 0;

	if (i > j) {
		max = i;
	} else {
		max = j;
	}

	int end_position_of_best_hybrid = 0;
	if (BEST_HYBRID_i < BEST_HYBRID_j) {
		end_position_of_best_hybrid = BEST_HYBRID_i;
	} else {
		end_position_of_best_hybrid = BEST_HYBRID_j;
	}


	for (int z = max; z > end_position_of_best_hybrid; z--) {
		if(i > BEST_HYBRID_i || j > BEST_HYBRID_j){
			alignment = alignment + "#";
			if (i > BEST_HYBRID_i) {
				miRNA_alignment = miRNA_alignment + miRNA_sequence[i];
				i--;
			} else {
				miRNA_alignment = "-" + miRNA_alignment;
			}
			if (j > BEST_HYBRID_j) {
				mRNA_alignment = mRNA_alignment + mRNA_sequence[j];
				j--;
			} else {
				mRNA_alignment = "-" + mRNA_alignment;
			}
		}
	}

	//At this time, the i should be identical to BEST_HYBRID_i and the j should be identical to BEST_HYBRID_j
	if (i != BEST_HYBRID_i && j != BEST_HYBRID_j) {
		cerr << "The end of the hybrid is not in HYBRIDIZED_PAIR state...something is wrong\n";
		exit(1);
	}
	cout << "End   of Hybrid " << "\tmiRNA_position=" << (BEST_HYBRID_i + 1) << "\tmRNA_position=" << (BEST_HYBRID_j + 1) << endl;

	int trace_i = BEST_HYBRID_i;
	int trace_j = BEST_HYBRID_j;

	current_MODE = HYBRIDIZED_PAIR;
	
	while(current_MODE != BEGIN_OF_THE_HYBRID){
		int miRNA_base = convert(miRNA_sequence[trace_i]);
                int mRNA_base =  convert(mRNA_sequence[trace_j]);
		
		make_alignment(miRNA_id, mRNA_id, Coefficients[miRNA_base][mRNA_base], trace_i, trace_j, current_MODE, miRNA_alignment, alignment, mRNA_alignment);
		i = trace_i;
		j = trace_j;
		trace_i = Trace_Back[i][j][current_MODE][TRACE_BACK_i];
		trace_j = Trace_Back[i][j][current_MODE][TRACE_BACK_j];
		current_MODE = Trace_Back[i][j][current_MODE][TRACE_BACK_STATE];
	}
	cout << "Begin of Hybrid " << "\tmiRNA_position=" << (i + 1) << "\tmRNA_position=" << (j + 1) << endl;


	i--;
	j--;

	//Extend the alignments
	max = 0;
	if (i > j) {
		max = i;
	} else {
		max = j;
	}
	
	for (int z = max; z >= 0; z--) {
		if(i >= 0 || j >= 0){
			alignment = alignment + "#";
			if (i >= 0) {
				miRNA_alignment = miRNA_alignment + miRNA_sequence[i];
				i--;
			} else {
				miRNA_alignment = miRNA_alignment + "-";
			}
			if (j >= 0) {
				mRNA_alignment = mRNA_alignment + mRNA_sequence[j];
				j--;
			} else {
				mRNA_alignment = mRNA_alignment + "-";
			}
		}
	}

	//Since we back-tracked we need to reverse the alignment
	miRNA_alignment=sequence_reverse(miRNA_alignment);
	mRNA_alignment=sequence_reverse(mRNA_alignment);
	alignment=sequence_reverse(alignment);


	//Print out the alignment
	cout << "######## Hybridyzation ######## " << endl;
	cout << mRNAs.GetSequence(mRNA_id).header
		 <<	"\t<<< VERSUS >>>\t"
		 << miRNAs.GetSequence(miRNA_id).header
		 << "\t----> Ratio(H/Z0)= " << Ratio
		 << endl;
	cout << "miRNA             5' "
 		 << "\t"
 		 << miRNA_alignment
 		 << " 3'\t"
 		 << miRNA_sequence << endl;
 	cout << "A L I G N M E N T" << "\t" << alignment << endl;
 	cout << "mRNA   (reversed) 3' "
 		 << "\t"
 		 << mRNA_alignment
 		 << " 5'\t"
 		 << mRNA_sequence << endl;
 	cout << endl;


 	return;
}

void Hybridize::make_alignment(unsigned int miRNA_id, unsigned int mRNA_id, double coefficient, int i, int j, int MODE, string &miRNA_alignment, string &alignment, string &mRNA_alignment) {
	string mRNA_sequence = mRNAs.GetSequence(mRNA_id).reversed_sequence;
	string miRNA_sequence = miRNAs.GetSequence(miRNA_id).sequence;

	if(MODE == BEGIN_OF_THE_HYBRID){
		MODE = HYBRIDIZED_PAIR;
	}
	switch (MODE) {
	case HYBRIDIZED_PAIR:
		miRNA_alignment = miRNA_alignment + miRNA_sequence[i];
		mRNA_alignment = mRNA_alignment + mRNA_sequence[j];
		if (coefficient) {
			alignment = alignment + "|";
		} else {
			alignment = alignment + " ";
			cout << "How can you hybridize position "
					<< (i + 1) << " of miRNA "
					<< miRNAs.GetSequence(miRNA_id).header
					<< "(" << miRNA_sequence
					<< ") which is ["
					<< miRNA_sequence[i]
					<< "] to position " << (j + 1) << " of mRNA "
					<<  mRNAs.GetSequence(mRNA_id).header
					<< "(" << mRNA_sequence
					<< ") which is ["
					<< mRNA_sequence[j]
					<< "]"
					<< endl;
		}
		break;
	case UNMATCHED_SYMMETRIC_LOOP:
		miRNA_alignment = miRNA_alignment + miRNA_sequence[i];
		mRNA_alignment = mRNA_alignment + mRNA_sequence[j];
		alignment = alignment + "O";
		break;
	case UNMATCHED_ASYMMETRIC_mRNA:
		miRNA_alignment = miRNA_alignment + "-";
		mRNA_alignment = mRNA_alignment + mRNA_sequence[j];
		alignment = alignment + "v";
		break;
	case UNMATCHED_ASYMMETRIC_miRNA:
		miRNA_alignment = miRNA_alignment + miRNA_sequence[i];
		mRNA_alignment = mRNA_alignment + "-";
		alignment = alignment + "^";
		break;
		;
	}
	return;
}

string Hybridize::sequence_reverse(string s) {

	int begin = 0;
	int end = s.length() - 1;
	while (begin < end) {
		swap(s[begin++], s[end--]);
	}
	return s;
}

int Hybridize::convert(char a) {
	if (a == 'a' || a == 'A')
		return (int) A;
	else if (a == 'c' || a == 'C')
		return (int) C;
	else if (a == 'g' || a == 'G')
		return (int) G;
	else if (a == 'u' || a == 'U')
		return (int) T;
	else if (a == 't' || a == 'T')
		return (int) T;
	else
		return (int) N;
}

void Hybridize::Set_mRNA(string s) {
	mRNA = s;
	return;
}

string Hybridize::Get_mRNA() {
	return mRNA;
}

void Hybridize::Set_miRNA(string s) {
	miRNA = s;
	return;
}

string Hybridize::Get_miRNA() {
	return miRNA;
}

double Hybridize::Get_Max_Hybrid_state(double one,double hybrid,double symmetric_loop,double mRNA_assymetric_loop, double miRNA_assymetric_loop, int &STATE){

	double maximum = one;
	STATE = BEGIN_OF_THE_HYBRID;
	if (hybrid > maximum) {
		maximum = hybrid;
		STATE = HYBRIDIZED_PAIR;
	}
	if(symmetric_loop > maximum){
		maximum = symmetric_loop;
		STATE = UNMATCHED_SYMMETRIC_LOOP;
	}
	if (mRNA_assymetric_loop > maximum) {
		maximum = mRNA_assymetric_loop;
		STATE = UNMATCHED_ASYMMETRIC_mRNA;
	}
	if (miRNA_assymetric_loop > maximum) {
		maximum = miRNA_assymetric_loop;
		STATE = UNMATCHED_ASYMMETRIC_miRNA;
	}

	return maximum;
}

double Hybridize::Get_Max_symmetric_loop_state(double hybrid,double symmetric_loop, int &STATE){

	double maximum = hybrid;
	STATE = HYBRIDIZED_PAIR;
	if (symmetric_loop > maximum) {
		maximum = symmetric_loop;
		STATE = UNMATCHED_SYMMETRIC_LOOP;
	}

	return maximum;
}

double Hybridize::Get_Max_assymetric_mRNA_loop_state(double hybrid,double symmetric_loop,double mRNA_assymetric_loop,int &STATE){

	double maximum = hybrid;
	STATE = HYBRIDIZED_PAIR;
	if (symmetric_loop > maximum) {
		maximum = symmetric_loop;
		STATE = UNMATCHED_SYMMETRIC_LOOP;
	}
	if (mRNA_assymetric_loop > maximum) {
		maximum = mRNA_assymetric_loop;
		STATE = UNMATCHED_ASYMMETRIC_mRNA;
	}

	return maximum;
}

double Hybridize::Get_Max_assymetric_miRNA_loop_state(double hybrid,double symmetric_loop,double miRNA_assymetric_loop,int &STATE){

	double maximum = hybrid;
	STATE = HYBRIDIZED_PAIR;
	if (symmetric_loop > maximum) {
		maximum = symmetric_loop;
		STATE = UNMATCHED_SYMMETRIC_LOOP;
	}
	if (miRNA_assymetric_loop > maximum) {
		maximum = miRNA_assymetric_loop;
		STATE = UNMATCHED_ASYMMETRIC_miRNA;
	}

	return maximum;
}

void Hybridize::Get_miRNA_Sequences(string miRNA_fileName) {

	ifstream infile;
	infile.open(miRNA_fileName.c_str(), ios::in);
	if (!infile) {
		cerr << "Can't open input file " << miRNA_fileName << endl;
		exit(1);
	} else {
		cerr << "Reading Fasta File: " << miRNA_fileName << endl;
		number_of_miRNA_sequences=0;
		while (!infile.eof()) {
			string my_header = "";
			string my_seq = "";
			infile >> my_header;
			infile >> my_seq;
			miRNA_Sequences[my_header]=my_seq;
			number_of_miRNA_sequences++;
		}
		infile.close();
	}
	cerr << "Number of read sequences = " << number_of_miRNA_sequences << endl;
}

void Hybridize::Get_mRNA_Sequences(string mRNA_fileName) {

	ifstream infile;
	infile.open(mRNA_fileName.c_str(), ios::in);
	if (!infile) {
		cerr << "Can't open input file " << mRNA_fileName << endl;
		exit(1);
	} else {
		cerr << "Reading Fasta File: " << mRNA_fileName << endl;
		number_of_mRNA_sequences=0;
		while (!infile.eof()) {
			string my_header = "";
			string my_seq = "";
			infile >> my_header;
			infile >> my_seq;
			number_of_mRNA_sequences++;
			mRNA_Sequences[my_header]=my_seq;
			string reversed_sequence= sequence_reverse(my_seq);
			mRNA_Reversed_Sequences[my_header] = reversed_sequence;
		}
		infile.close();
	}
	cerr << "Number of read sequences = " << number_of_miRNA_sequences << endl;
}

void Hybridize::Read_mRNA_fasta_file(string fileName){

	this->mRNAs.Set_fileName(fileName);
	this->mRNAs.ReadFasta();

	return;
}

void Hybridize::Read_miRNA_fasta_file(string fileName){
	this->miRNAs.Set_fileName(fileName);
	this->miRNAs.ReadFasta();

	return;
}

void Hybridize::Read_miRNA_expressions(string fileName){
	//Get the miRNA expressions as their priors

	//Read Expression Profile
	ifstream infile;
	infile.open(fileName.c_str(), ios::in);
	if (!infile) {
		cerr << "Can't open input file " << fileName << endl;
		exit(1);
	} else {
		cerr << "Reading miRNA expression File: " << fileName << endl;
	}
	while (!infile.eof()) {
		string my_mirna_header = "";
		double my_mirna_expression;
		infile >> my_mirna_header;
		infile >> my_mirna_expression;
		my_mirna_header = '>' + my_mirna_header;
		miRNA_Expression_Hash[my_mirna_header] = my_mirna_expression;
	}
	return;
}

Hybridize::~Hybridize() {
}

