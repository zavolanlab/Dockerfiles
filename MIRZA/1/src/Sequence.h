// ===========================================================================
// Name        : Sequence.h
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

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#define A 0
#define C 1
#define G 2
#define T 3
#define U 3
#define N 4
#define NOTHING -1

class Sequence {

struct SEQUENCE{
	string header;
	string sequence;
	string reversed_sequence;
	std::vector<int> Converted_Sequence;
	std::vector<int> Reversed_Converted_Sequence;
	SEQUENCE(){
		cerr << "header and sequence are empty, nothing inserted" <<endl;
		this->header="";
		this->sequence="";
		this->reversed_sequence="";


	}
	SEQUENCE(string h, string s){
		this->header=h;
		this->sequence=s;

		//Calculate and save the reversed sequence
		this->reversed_sequence=s;
		int begin = 0;
		int end = s.length() - 1;
		while (begin < end) {
			swap(reversed_sequence[begin++], reversed_sequence[end--]);
		}

		//Calculate the Converted sequence
		this->Converted_Sequence.resize(s.length(),NOTHING);
		this->Reversed_Converted_Sequence.resize(s.length(),NOTHING);

		for(unsigned int i=0;i<s.length();i++){
			if (s[i] == 'a' || s[i] == 'A')
				Converted_Sequence[i]= (int) A;
			else if (s[i] == 'c' || s[i] == 'C')
				Converted_Sequence[i]= (int) C;
			else if (s[i] == 'g' || s[i] == 'G')
				Converted_Sequence[i]= (int) G;
			else if (s[i] == 'u' || s[i] == 'U')
				Converted_Sequence[i]= (int) T;
			else if (s[i] == 't' || s[i] == 'T')
				Converted_Sequence[i]= (int) T;
			else
				Converted_Sequence[i]= (int) N;
		}

		//Calculate and save the reversed of the converted sequence
		begin = 0;
		end = s.length() - 1;
		while (begin < end) {
			Reversed_Converted_Sequence[end]=Converted_Sequence[begin];
			Reversed_Converted_Sequence[begin]=Converted_Sequence[end];
			begin++;
			end--;
		}
	}
};

private:
	string fileName;
	unsigned int Size;
	std::vector<SEQUENCE> SEQUENCEVector; //Array of Sequences

public:

	void AddSequence(string header, string sequence);
	SEQUENCE GetSequence(int );
	void     SetSequence(int ,string );
	//string sequence_reverse(string s);

	Sequence();
	virtual ~Sequence();
	void Set_fileName(string);
	string Get_fileName();
	void Set_Size(unsigned int);
	unsigned int Get_Size();
	void ReadFasta();
};



#endif /* SEQUENCE_H_ */
