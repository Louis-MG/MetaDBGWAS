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

#include "build_dbg.hpp"
#include "global.h"
#include "GraphOutput.h"

using namespace std;

// The Tool constructor allows to give a name to our tool.
// This name appears when one gets help in the command line or in the final output
build_dbg::build_dbg ()  : Tool ("build_dbg") //give a name to our tool
{
    populateParser(this);
}

class EdgeConstructionVisitor : public boost::static_visitor<>    {
private:
    const string& linear_seqs_name;
/*
 * This small part of the code refers to GraphOutput.h which contains functions to build the .edges and .nodes files. Only the Unitigs names must be given.
 * It takes care of writing the "FF" (Forward to Forward), "FR" (Forward to Reverse), and other annotations.
 */
public:
    EdgeConstructionVisitor (const string &linear_seqs_name) : linear_seqs_name(linear_seqs_name) {}
    template<size_t span>
    void operator() (GraphOutput<span>& graphOutput) const
    {
        graphOutput.open();
        graphOutput.load_nodes_extremities(linear_seqs_name);
        graphOutput.construct_graph(linear_seqs_name);
        graphOutput.close();
    }
};


/*********************************************************************
** METHOD  :
** PURPOSE : produce a debruijn graph with the unitigs
** INPUT   : fasta file with unitigs from Bcalm
** OUTPUT  : outputs some stats
** RETURN  : .edges and .nodes files
** REMARKS :
*********************************************************************/
void build_dbg::execute ()
{
    cerr << "Creating .edges and .nodes files ..." << endl;
    //get the parameters
    int kmerSize = getInput()->getInt(STR_KSKMER_SIZE);
    //gets the number of cores to use
    int nbCores = getInput()->getInt(STR_NBCORES);
    //creates variable where the unitigs file is
    string linear_seqs_name = getInput()->getStr(STR_PATH_TO_FASTA_FILES);
    //create the step1 folder in the outputfolder
    string outputFolder = stripLastSlashIfExists(getInput()->getStr(STR_OUTPUT)); //TODO: should I keep the files in the common output folder or in a separate one ?

    //Builds the DBG with unitigs using GATB
    graph = new Graph ; //TODO: changer en GaphUnitig (GraphUnitigs.hpp) when the code will be usable (god knows when)
    //test line to load the .h5 from bcalm :
    // *graph = gatb::core::debruijn::impl::Graph::load(outputFolder + "/unitigs/fof.h5");
    // line that works :
    *graph = gatb::core::debruijn::impl::Graph::create("-in %s -kmer-size %d -abundance-min 0 -out %s/graph -nb-cores %d", //TODO: same
                                                       linear_seqs_name.c_str(), kmerSize, outputFolder.c_str(), nbCores);

    //builds and outputs .nodes and .edges.dbg files, see GraphOutput.h for the inner code
    typedef boost::variant <
        GraphOutput<KMER_SPAN(0)>,
        GraphOutput<KMER_SPAN(1)>,
        GraphOutput<KMER_SPAN(2)>,
        GraphOutput<KMER_SPAN(3)>
    >  GraphOutputVariant;

    GraphOutputVariant graphOutput;
    if (kmerSize < KMER_SPAN(0))  {  graphOutput = GraphOutput<KMER_SPAN(0)>(graph, outputFolder+string("/graph")); }
    else if (kmerSize < KMER_SPAN(1))  {  graphOutput = GraphOutput<KMER_SPAN(1)>(graph, outputFolder+string("/graph")); }
    else if (kmerSize < KMER_SPAN(2))  {  graphOutput = GraphOutput<KMER_SPAN(2)>(graph, outputFolder+string("/graph")); }
    else if (kmerSize < KMER_SPAN(3))  {  graphOutput = GraphOutput<KMER_SPAN(3)>(graph, outputFolder+string("/graph")); }
    else { throw gatb::core::system::Exception ("Graph failure because of unhandled kmer size %d", kmerSize); }
    boost::apply_visitor (EdgeConstructionVisitor(linear_seqs_name),  graphOutput);

    //save disk space
    //remove(linear_seqs_name.c_str());

    //print some stats
    cout << "################################################################################" << endl;
    cout << "Stats: " << endl;
    cout << "Number of kmers: " << graph->getInfo()["kmers_nb_solid"]->getInt() << endl;
    cout << "Number of unitigs: " << getNbLinesInFile(outputFolder+string("/graph.nodes")) << endl;
    cout << "################################################################################" << endl;
}
