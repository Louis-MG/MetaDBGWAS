/*
## Copyright (C) <2017>  <bioMerieux, Universite Claude Bernard Lyon 1,
## Centre National de la Recherche Scientifique>

## 1. This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU Affero General Public License as published
## by the Free Software Foundation version 3 of the  License and under the
## terms of article 2 below.
## 2. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
## or FITNESS FOR A PARTICULAR PURPOSE. See below the GNU Affero General
## Public License for more details.
## You should have received a copy of the GNU Affero General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 3. Communication to the public by any means, in particular in the form of
## a scientific paper, a poster, a slideshow, an internet page, or a patent,
## of a result obtained directly or indirectly by running this program must
## cite the following paper :
##  Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix,
##  Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic
##  Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017,
##  Cold Spring Harbor Labs Journals, doi:10.1101/113563.
##  (url: http://www.biorxiv.org/content/early/2017/03/03/113563)
## -------------------------------------------------------------------------

## Authors (alphabetically): Jacob L., Jaillard M., Lima L.
## -------------------------------------------------------------------------

## DISCLAIMER : This file was modified by Louis-Mael Gueguen. THe modified version of the program is published under the dual licensing AGPLv3, and zlib license.
*/

#ifndef _TOOL_map_reads_HPP_
#define _TOOL_map_reads_HPP_
#include <cstdlib>
/********************************************************************************/
#include <gatb/gatb_core.hpp>
/********************************************************************************/

//this struct contains the name of a kmer, and pattern is a vector of its absence/presence
struct SKmer {
    std::vector<int> pattern;
    bool corrected;
};


// declarations
SKmer process_line(const std::string& line_buffer);
//SKmer binarise_counts(SKmer& data); // TODO: supp
SKmer minor_allele_description(SKmer& data);
void write_bugwas_gemma(const std::string& outputFolder, const std::vector<std::vector<int>>& vector_of_unique_patterns, std::string& rawname, std::vector<std::string>& filenames, std::map<std::vector<int>, std::vector<int>>& map_unique_to_all);


class input_files_gen : public Tool
{
public:

    // Constructor
    input_files_gen ();

    // Actual job done by the tool is here
    [[noreturn]] void execute ();

    //overriding this in order to exit the tool when finding a problem with the arguments
    IProperties* run (int argc, char* argv[])
    {
        IProperties* toReturn = Tool::run(argc, argv);
        if (!toReturn)
            std::exit(1);
        return toReturn;
    }
};

/********************************************************************************/

#endif

