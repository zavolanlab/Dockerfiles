// ===========================================================================
// Name        : Sequence.cpp
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
#include "Sequence.h"

using namespace std;

Sequence::Sequence() {
	// TODO Auto-generated constructor stub
	this->Set_Size(-1);
	this->Set_fileName("");
}

void Sequence::Set_fileName(string s) {
	this->fileName = s;
}
string Sequence::Get_fileName() {
	return this->fileName;
}
void Sequence::Set_Size(unsigned int s) {
	this->Size = s;
}
unsigned int Sequence::Get_Size() {
	return this->Size;
}


void Sequence::AddSequence(string header, string seq) {
	try {
		this->SEQUENCEVector.push_back(SEQUENCE(header,seq));
		this->Size++;
	} catch (int ExNum) {
		cerr << "While adding sequences, An exception with number =" << ExNum
				<< "has been caught " << endl;
		exit(1);
	}
}

Sequence::SEQUENCE Sequence::GetSequence(int i) {
	return this->SEQUENCEVector[i];
}

void Sequence::SetSequence(int i,string s) {
	SEQUENCEVector[i].sequence=s;
}

void Sequence::ReadFasta() {

	ifstream infile;
	infile.open(fileName.c_str(), ios::in);
	if (!infile) {
		cerr << "Can't open input file " << fileName << endl;
		exit(1);
	} else {
		cerr << "Reading Fasta File: " << fileName << endl;
		this->Size = -1;
		while (!infile.eof()) {
			string my_header = "";
			string my_seq = "";
			infile >> my_header;
			infile >> my_seq;
			AddSequence(my_header, my_seq);
		}
		infile.close();
	}
	cerr << "Number of read sequences = " << this->Size << endl;
}

Sequence::~Sequence() {
	// TODO Auto-generated destructor stub
}

