// ===========================================================================
// Name        : Hybridize.h
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
#include <string>
#include <iostream>
#include "Parameters.h"
#include "Sequence.h"
#include <map>
#include <vector>

#ifndef HYBRIDIZE_H_
#define HYBRIDIZE_H_

#define NOTHING -1
#define STOP -10

//Parameter length configurations
#define NUMBER_OF_HYBRIDIZATION_MODES 4 // Number of states which is (h, sym, mRNA, miRNA)
#define miRNA_HYBRIDIZATION_ENERGY_PARAMETER_OFFSET 6

//TRACE BACK elements
#define NUMBER_OF_TRACE_BACK_ELEMENTS 3
#define TRACE_BACK_i 0
#define TRACE_BACK_j 1
#define TRACE_BACK_STATE 2

//Noise
#define NOISE 0.05

//ONE
#define ONE 1
#define ZERO 0

//Nucleotides
#define A 0
#define C 1
#define G 2
#define T 3
#define U 3
#define N 4

//Global Energy parameters
#define E_AT  0
#define E_GC  1
#define E_OPEN  2
#define E_mRNA_ASSYMETRIC_LOOP  3
#define E_miRNA_ASSYMETRIC_LOOP  4
#define E_SYMMETRIC_LOOP  5

//Hybridization modes
#define BEGIN_OF_THE_HYBRID 5
#define HYBRIDIZED_PAIR 0
#define UNMATCHED_SYMMETRIC_LOOP 1
#define UNMATCHED_ASYMMETRIC_mRNA 2
#define UNMATCHED_ASYMMETRIC_miRNA 3
#define HYBRIDIZATION_TYPE 4

#define ACCURACY 0.0001
#define BIG_INTEGER 10000

#define TOTAL_NUMBER_OF_mRNAs 5000
#define TOTAL_NUMBER_OF_miRNAs	500
#define mRNA_LENGTH 50
#define miRNA_LENGTH 21

#define NUMBER_OF_NUCLEOTIDES 4
class Hybridize {
	typedef std::vector<double> double_vector;
	typedef std::vector<double_vector> double_2d_vector;
	typedef std::vector<double_2d_vector> double_3d_vector;
	typedef std::vector<double_3d_vector> double_4d_vector;

	typedef std::vector<int> int_vector;
	typedef std::vector<int_vector> int_2d_vector;
	typedef std::vector<int_2d_vector> int_3d_vector;
	typedef std::vector<int_3d_vector> int_4d_vector;

private:

	string mRNA;
	string miRNA;

	std::map<string, double> miRNA_Expression_Hash;


	std::map<string, string> mRNA_Sequences;
	int number_of_mRNA_sequences;

	std::map<string, string> miRNA_Sequences;
	std::map<string, string> mRNA_Reversed_Sequences;
	int number_of_miRNA_sequences;

	double total_miRNA_expressions;

	void fill_the_F_and_H_matrices(unsigned int miRNA_id, unsigned int mRNA_id, double coefficient, int i, int j);
	void fill_the_Best_Hybrid_matrix(unsigned int miRNA_id, unsigned int mRNA_id, double coefficient, int i, int j);
	void make_alignment(unsigned int miRNA_id, unsigned int mRNA_id, double coefficient, int i, int j, int MODE, string &miRNA_alignment, string &alignment, string &mRNA_alignment);

	int BEST_HYBRID_i;
	int BEST_HYBRID_j;
	double BEST_HYBRID_SCORE;
	int END_STATE;
	


public:

	int mRNA_Length;
	int updatePriors;

	Sequence mRNAs;
	Sequence miRNAs;
	Parameters Params;
	double ***F_alpha;
	double ***H_alpha;
	double ***Best_Hybrid;
	int    ****Trace_Back;
	double Ratio[TOTAL_NUMBER_OF_miRNAs][TOTAL_NUMBER_OF_mRNAs];
	double miRNA_Priors[TOTAL_NUMBER_OF_miRNAs];
	double mRNA_likelihood_ratios[TOTAL_NUMBER_OF_mRNAs];

	double Gaussian_Penalty;
	double Z0;
	double H;

	//Exponents tables
	double eE_hybridization[TOTAL_NUMBER_OF_miRNAs][miRNA_LENGTH];
	double Coefficients[NUMBER_OF_NUCLEOTIDES][NUMBER_OF_NUCLEOTIDES];
	double eE_SYM;
	double eE_OPEN;

	double eE_mRNA_Assymetric_loops;
	double eE_miRNA_Assymetric_loops;



	//some sequence modification functions
	string sequence_reverse(string );
	int convert(char );

	void Read_mRNA_fasta_file(string);
	void Read_miRNA_fasta_file(string);
	void Read_miRNA_expressions(string fileName);


	void Get_mRNA_Sequences(string mRNA_fileName);
	void Get_miRNA_Sequences(string miRNA_fileName);

	void Set_mRNA(string);
	string Get_mRNA();
	void Set_miRNA(string);         // It REVERSEs the sequence of miRNA when it sets the miRNA.
	string Get_miRNA();

	void Initialize_Global_Parameters_and_Prepare(double []);

	void Initialize(void);
	void resetMemory(void);
	void Estimate_the_mRNA_base_content(void);
	void Initialize_miRNA_mRNA_likelihood_ratios(void);
	void Initialize_miRNA_priors(void);
	void Initialize_mRNA_log_likelihood_ratios(void);
	double CalculateHibridizationRatio_miRNA_mRNA(unsigned int miRNA_id , unsigned int mRNA_id); //This function will populate the partition sum F0, Fa, Ha, Z0,Z1 and returns the Likelihood of ONE mRNA and miRNA hybridization

	void Update_miRNA_priors(void);
	double Calculate_data_log_likelihood_ratio(void);

	void Trace_Back_and_Print_Hybridization(unsigned int miRNA_id , unsigned int mRNA_id, double Ratio);

	double Get_Max_Hybrid_state(double one,double hybrid,double symmetric_loop,double mRNA_assymetric_loop, double miRNA_assymetric_loop, int &STATE);
	double Get_Max_symmetric_loop_state(double hybrid,double symmetric_loop, int &STATE);
	double Get_Max_assymetric_mRNA_loop_state(double hybrid,double symmetric_loop,double mRNA_assymetric_loop, int &STATE);
	double Get_Max_assymetric_miRNA_loop_state(double hybrid,double symmetric_loop,double miRNA_assymetric_loop, int &STATE);

	Hybridize();
	virtual ~Hybridize();
};

#endif /* HYBRIDIZE_H_ */
