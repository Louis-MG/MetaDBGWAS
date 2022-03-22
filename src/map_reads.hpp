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

#ifndef _TOOL_map_reads_HPP_
#define _TOOL_map_reads_HPP_
#include <cstdlib>
/********************************************************************************/
#include <gatb/gatb_core.hpp>
/********************************************************************************/

//this struct contains the name of a kmer, and pattern is a vector of its absence/presence
struct SKmer {
    std::string name;
    std::vector<int> pattern;
    int corrected;
};


// declarations
SKmer process_line(const std::string& line_buffer);
void write_bugwas_gemma(const std::vector<std::vector<int>>& vector_of_unique_patterns, std::string& rawname, std::vector<std::string>& filenames, std::map<std::vector<int>, std::vector<int>>& map_unique_to_all);
SKmer minor_allele_description(SKmer data);


class map_reads : public Tool
{
public:

    // Constructor
    map_reads ();

    // Actual job done by the tool is here
    void execute ();

    //overriding this in order to exit the tool when finding a problem with the arguments
    IProperties* run (int argc, char* argv[])
    {
        IProperties* toReturn = Tool::run(argc, argv);
        if (!toReturn)
            std::exit(1);
        return toReturn;
    }
};

SKmer process_line(const std::string& line_buffer) {
    /*
     * this function processes lines by putting them in a structure than contains  the Kmer name, a vector of its absence/presence pattern. Ignores the corrected attribute.
     */
    std::vector<int> output_pattern;
    std::string kmer_name;
    std::istringstream input(line_buffer);
    //loop over tab-separated words in the line :
    for (std::string word; std::getline(input, word, '\t'); ) {
        if (word.starts_with(">")) { // if id resets line
            kmer_name = word;
        } else if (word.ends_with("*")) { // if abundance is 0
            output_pattern.push_back(0);
        } else {
            output_pattern.push_back(1); // if abundance is anything else than 0
        }
    }
    // Kmer output_struct{kmer_name, output_pattern};
    // output_pattern.clear();
    return {kmer_name, output_pattern};
}

SKmer minor_allele_description(SKmer data) {
    /*
     * this function changes the pattern of presence/absence of a SKmer into the minor allele description if needed, and changed the 'corrected' accordingly (1: did not change; -1: changed).
     */
    float sum ;
    std::vector<int> corr_vector;
    for (auto& n : data.pattern) {
        sum += n;
    }
    if (sum/(float)data.pattern.size() > 0.5) {
        for (int i = 0; i < data.pattern.size(); i++) {
            std::cout << i << std::endl;
            switch (data.pattern.at(i)) {
                case 0:
                    corr_vector.push_back(1);
                    break;
                case 1:
                    corr_vector.push_back(0);
                    break;
            }
        }
        data.corrected = -1;
        data.pattern = corr_vector;
    } else {
        data.corrected = 1;
    }
    return data;
}


void write_bugwas_gemma(const std::vector<std::vector<int>>& vector_of_unique_patterns, std::string& rawname, std::vector<std::string>& filenames, std::map<std::vector<int>, std::vector<int>>& map_unique_to_all) {
    /*
     * this function builds output files : unique_patterns, unique_to_all, and gemma_pattern_to_nb_unitigs, gemma_unitig_to_patterns.
     */
    std::ofstream outstream_unique (rawname+".unique_rows.binary", std::ofstream::binary);
    std::ofstream outstream_unique_to_all (rawname+".unique_rows_to_all_rows.binary", std::ofstream::binary);
    std::ofstream outstream_gemma_pattern_to_nb_unitigs ("gemma_input.pattern_to_nb_of_unitigs.binary", std::ofstream::binary);
    std::ofstream outstream_gemma_unitig_to_patterns ("gemma_input.unitig_to_pattern.binary", std::ofstream::binary);


    //error check
    if (outstream_unique.fail()) {
        perror("ERROR: could not open bugwas_input.unique_rows.binary");
        std::exit(1);
    } else if (outstream_unique_to_all.fail()) {
        perror("ERROR: could not open bugwas_inputunique_rows_to_all_rows.binary");
        std::exit(1);
    } else if (outstream_gemma_pattern_to_nb_unitigs.fail()) {
        perror("ERROR: could not open gemma_input.pattern_to_nb_of_unitigs.binary");
        std::exit(1);
    } else if (outstream_gemma_unitig_to_patterns.fail()) {
        perror("ERROR: could not open gemma_input.unitig_to_pattern.binary");
        std::exit(1);
    }

    //TODO: check that my understanding of the output files is ok

    // header for unique_pattern
    outstream_unique << "ps ";
    for (const std::string &i : filenames) {
        outstream_unique << i << " ";
    }
    outstream_unique << "\n";

    int n = 0;
    for (const auto &i : vector_of_unique_patterns) {
        // writes the unique patterns in their file
        outstream_unique << n << " ";
        for (const auto &j : i) {
            outstream_unique << j << " " ;
        }
        outstream_unique << "\n" ;
        // writes connection between unique pattern and all the unitigs each represents in unique_to_all, and the reciprocal in unitig_to_pattern
        for (const auto &j : map_unique_to_all[i]) {
            outstream_unique_to_all << j << " ";
            outstream_gemma_unitig_to_patterns << j << " " << n << "\n";
        }
        outstream_unique_to_all << "\n" ;

        // writes the number of unitigs that each unique pattern represents
        outstream_gemma_pattern_to_nb_unitigs << n << " " << map_unique_to_all[i].size() << "\n";

        n++;
    }
    outstream_unique.close();
    outstream_unique_to_all.close();
    outstream_gemma_pattern_to_nb_unitigs.close();
    outstream_gemma_unitig_to_patterns.close();
}

/********************************************************************************/

#endif

