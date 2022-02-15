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

#include <string>
#include "Utils.h"

using namespace std;


const char* STR_MAX_NEIGHBOURHOOD = "-nh";
const char* STR_KSKMER_SIZE = "-k";
const char* STR_OUTPUT = "-output";
const char* STR_NBCORES = "-nb-cores";
const char* STR_KEEP_NA = "-keepNA";

//TODO: several questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
//const char* STR_COUNT_MODE = "-count";
//TODO: several questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option

//variables controlling where the executable is
string DBGWAS_lib = ""; //set on first command in main.cpp


//global vars used by both programs
Graph* graph;
vector< UnitigIdStrandPos >* nodeIdToUnitigId;

void populateParser (Tool *tool) {
  // We add some custom arguments for command line interface
  tool->getParser()->push_front (new OptionNoParam (STR_KEEP_NA, "Keep strains with phenotype NA.",  false));
  tool->getParser()->push_front (new OptionOneParam (STR_MAX_NEIGHBOURHOOD, "Denotes the neighbourhood to be considered around the significant unitigs.",  false, "5"));

  tool->getParser()->push_front (new OptionOneParam (STR_OUTPUT, "Path to the folder where the final and temporary files will be stored.",  false, "output"));
  tool->getParser()->push_front (new OptionOneParam (STR_KSKMER_SIZE, "K-mer size.",  false, "31"));

  //TODO: several questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
  //tool->getParser()->push_front (new OptionOneParam (STR_COUNT_MODE, "The count mode. If \"01\", then the count mode is seen as presence/absence. If \"Freq\", then the count mode is seen as the frequency",  false, "01"));
  //TODO: several questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
}
