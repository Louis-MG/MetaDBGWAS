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
*/

#include "global.h"
#include "map_reads.hpp"
#include "Utils.h"
#include "PhenoCounter.h"
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <vector>
#include <set>
#include <map>

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
// The Tool constructor allows to give a name to our tool.
// This name appears when one gets help in the command line or in the final output
map_reads::map_reads ()  : Tool ("map_reads") //give a name to our tool
{
    populateParser(this);
}



/*********************************************************************
** METHOD  :
** PURPOSE : executes map_reads
** INPUT   : output folder name, unitigs
** OUTPUT  : bugwas_input.* files (4) and weight_correction file
** RETURN  :
** REMARKS : deletes graph, NodeIdToUnitigs obects as well as the tmpFolder
*********************************************************************/
void map_reads::execute ()
{
    //populates strains file if given
    checkParametersMapping(this);

    //get the parameters
    string referenceOutputFolder = getInput()->getStr(STR_OUTPUT); //keep the original output folder in a variable for file manipulation
    string outputFolder = stripLastSlashIfExists(getInput()->getStr(STR_OUTPUT))+string("/step1");
    // measuring time
    auto t1 = std::chrono::high_resolution_clock::now();

    //generate the bugwas input
    std::cout << "[Starting conversion of abundance matrix to presence/absence matrix ...]" << std::endl;
    std::string matrix = referenceOutputFolder + "matrix/reindeer_output/TRUC";//TODO: replace with filename; //input abundance matrix
    std::string output = "bugwas_input.all_rows.binary" ;//TODO:same ; //output presence/absence matrix
    // streams
    std::ifstream stream (matrix, std::ifstream::binary);
    std::ofstream outstream (outputFolder + output, std::ofstream::binary);
    std::ofstream weight_corr_track ("weight_correction", std::ofstream::binary);
    // other files names prefix
    size_t lastindex = output.find_last_of(".");
    std::string rawname = output.substr(0, lastindex);
    // must use stream.fail() with a switch case to now if the thing went fine : no
    // file operation causes C++ to stop.
    // see https://www.eecs.umich.edu/courses/eecs380/HANDOUTS/cppBinaryFileIO-2.html
    if (stream.fail()) {
        perror("ERROR: could not open the reindeer matrix file");
        std::exit(1);
    } else if (outstream.fail()) {
        perror("ERROR: could not open bugwas_input.all_rows.binary");
        std::exit(1);
    } else if (weight_corr_track.fail()) {
        perror("ERROR: could not open weight_correction");
        std::exit(1);
    }

    // variables
    std::vector<SKmer> vector_of_kmers;
    std::vector<std::vector<int>> vector_of_unique_patterns ;
    std::vector<std::string> filenames ;
    std::vector<std::vector<int>> unique_to_all ;
    std::map<std::vector<int>, std::vector<int>> map_unique_to_all ; //each time a unique is accountered, insert the pattern as a key and the n in the vector of values
    // reads the input data
    outstream << "ps ";
    // read the data line by line:
    std::string line_buffer;
    std::set<std::vector<int>> vector_set;
    int n = 0; // line counter
    while(std::getline(stream, line_buffer).good()) {
        if (line_buffer.starts_with("query")) {
            // we obtain the file names and store them:
            std::istringstream input(line_buffer);
            for (std::string word; std::getline(input, word, '\t'); ) {
                filenames.push_back(word);
            }
            // removes "query"
            filenames.erase(filenames.begin());
            // write the header of both all_rows and all_rows_unique :
            for (const std::string &i : filenames) {
                outstream << i << " ";
            }
            outstream << "\n" ;
        } else {
            // we write the body:
            // 1: parse the line and build the SKmer
            SKmer raw_data = process_line(line_buffer);
            // 2: change, if needed, the allele description of the SKmer
            SKmer data = minor_allele_description(raw_data);
            // 3: keep track of the change in allele description
            weight_corr_track << data.corrected << "\n";
            // next
            vector_of_kmers.push_back(data);
            outstream << n << " " ;
            for (const auto &i : data.pattern) {
                outstream << i << " " ;
            }
            outstream << "\n" ;
            // looks for the vector in the set of unique vectors :
            auto search = vector_set.find(data.pattern) ;
            // if not found, adds it to the vector of unique vectors
            if (search == vector_set.end()) {
                vector_of_unique_patterns.push_back(data.pattern);
                vector_set.emplace(data.pattern);
                std::vector<int> unitigs{n};
                map_unique_to_all.emplace(data.pattern, unitigs);
            } else {
                map_unique_to_all[data.pattern].push_back(n);
                // add n to the item (vector) pointed out above
            }
            n++;
        }
    }
    stream.close();
    outstream.close();
    weight_corr_track.close();

    // writes uniques and unique_to_all, gemma unique patterns to nb unitigs outputs
    write_bugwas_gemma(vector_of_unique_patterns, rawname, filenames, map_unique_to_all);

    // finishing measuring time
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::minutes>(t2 - t1);
    std::cout << "Conversion took " << ms_int.count() << "min\n";

    //after the mapping, free some memory that will not be needed anymore
    delete graph;

    //clean-up - saving some disk space
    //remove GATB's graph file
    remove((referenceOutputFolder+string("/graph.h5")).c_str());

    cerr << endl << "[Convertion process finished!]" << endl;
    cerr.flush();
}
