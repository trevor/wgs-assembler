/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005, J. Craig Venter Institute. All rights reserved.
 * Author: Brian Walenz
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

#ifndef PROFILEREADERRORS_H
#define PROFILEREADERRORS_H

static const char* rcsid_PROFILEREADERRORS_H = "$Id: profileReadErrors.h,v 1.3 2011-08-31 09:10:43 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <map>
#include <string>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_IID.h"
#include "FileUtils.h"
#include "IllegalStateException.h"
#include "IOException.h"
#include "RuntimeException.h"
#include "StringUtils.h"

#define SNAPPER_FILE_LINE_BUFFER_SIZE 1536

#define INITIAL_ERROR_MATRIX_SIZE 1536
#define INITIAL_ERROR_MATRIX_BUCKET_SIZE 1024
#define ERROR_MATRIX_RESIZE_RESERVE_FACTOR 2

#define readSnapperFileLine(snapperFile, line, lineNum) \
	line = new char[SNAPPER_FILE_LINE_BUFFER_SIZE]; \
	line = FileUtils::readLine(snapperFile, line, SNAPPER_FILE_LINE_BUFFER_SIZE); \
	lineNum++;

#define isBaseCall(base) \
	(base == 'a') || (base == 'A') || \
	(base == 'c') || (base == 'C') || \
	(base == 'g') || (base == 'G') || \
	(base == 't') || (base == 'T')

#define isBaseError(base) \
	(base == 'A') || (base == 'C') || \
	(base == 'G') || (base == 'T')

#define isBaseGap(base) \
	base == '-'

#define isBaseUnknown(base) \
	base == 'N'

typedef enum AlignmentErrorType
{
	UNKNOWN, MISMATCH, INSERTION, DELETION
};

typedef struct AlignmentError
{
	AS_IID iid;
	AlignmentErrorType type;
	
	AlignmentError(AS_IID iid = 0, AlignmentErrorType type = UNKNOWN)
	{
		this->iid = iid;
		this->type = type;
	}
};

typedef struct BasePosition
{
	size_t position;
	size_t readsAtPosition;
	vector<AlignmentError> errors;
	double readLengthPercent;
	
	BasePosition(size_t position)
	{
		this->position = position;
		this->readsAtPosition = 0;
		this->errors.reserve(INITIAL_ERROR_MATRIX_BUCKET_SIZE);
		this->readLengthPercent = -1;
	}
	
	void addError(AlignmentError error, uint16 readLength)
	{
		this->errors.push_back(error);
		
		this->readLengthPercent = (this->readLengthPercent >= -1) ? 
			(this->readLengthPercent + ((double)this->position / readLength)) / 2 : 
			((double)this->position / readLength);
	}
};

void writeOutput(const char* outputFile, map<AS_IID, uint16>& readMap, vector<BasePosition*>& errorMatrix);

void processReadAlignment(AS_IID readIID, uint16 readLength, const char* readSequence, const char* genomeSequence, 
	vector<BasePosition*>& errorMatrix);
void processSnapperFile(const char* snapperFile, map<AS_IID, uint16>& readMap, vector<BasePosition*>& errorMatrix);

BasePosition* getBaseErrorBucket(vector<BasePosition*>& errorMatrix, size_t base);
AS_IID getReadIID(const char* readDefLine);
uint16 getReadLength(string readInfoLine);

void printUsage(const char* executableName);

int main(int argc, char** argv);

#endif
