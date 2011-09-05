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

#ifndef READALIGNMENTPROFILER_H
#define READALIGNMENTPROFILER_H

static const char* rcsid_READALIGNMENTPROFILER_H = "$Id: ReadAlignmentProfiler.h,v 1.4 2011-09-05 16:49:44 mkotelbajcvi Exp $";

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>

using namespace std;

#include "AlignmentDataReader.h"
#include "AlignmentError.h"
#include "AlignmentErrorType.h"
#include "AS_global.h"
#include "BaseAlignment.h"
#include "ErrorUtils.h"
#include "FileUtils.h"
#include "ReadAlignment.h"
#include "StringUtils.h"

using namespace Utility;

namespace ReadAnalysis
{
	static const uint16 PROFILE_DATA_PERCENT_INCREMENT = 10;
	
	static const char* PROFILE_DATA_OUTPUT_SUMMARY_PREFIX = "#";
	
	class ReadAlignmentProfiler
	{
	public:
		ReadAlignmentProfiler();
		~ReadAlignmentProfiler();
		
		inline static bool isBaseCall(const char base)
		{
			return (base == 'a') || (base == 'A') || 
				(base == 'c') || (base == 'C') || 
				(base == 'g') || (base == 'G') || 
				(base == 't') || (base == 'T');
		}
		
		inline static bool isBaseError(const char base)
		{
			return (base == 'A') || (base == 'C') || 
				(base == 'G') || (base == 'T');
		}
		
		inline static bool isBaseGap(const char base)
		{
			return base == '-';;
		}
		
		inline static bool isBaseUnknown(const char base)
		{
			return base == 'N';;
		}
		
		void writeProfile(string path);
		void writeProfile(FILE* stream);
		void profileData(vector<ReadAlignment*>& data, AlignmentDataStats& dataStats);
		
		bool& getVerbose()
		{
			return this->verbose;
		}
		
		void setVerbose(bool verbose)
		{
			this->verbose = verbose;
		}
		
	protected:
		bool verbose;
		vector<ReadAlignment*> data;
		AlignmentDataStats dataStats;
		vector<BaseAlignment*> bases;
		FILE* stream;
		
		void profileBase(size_t readIndex, size_t baseIndex);
		void initBases();
	};
}

#endif
