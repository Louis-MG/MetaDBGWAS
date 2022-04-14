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

#include "global.h"
#include <string>

using namespace std;

const char* STR_PATH_TO_FASTA_FILES = "--files";
const char* STR_KSKMER_SIZE = "--kmer";
const char* STR_OUTPUT = "--output";
const char* STR_NBCORES = "--threads";
const char* STR_STRAINS_FILE = "--strains";
const char* STR_KEEP_NA = "--keepNA";
const char* STR_PHENOTYPE_THRESHOLD = "-phenoThreshold";

//global vars used by both programs
Graph* graph;
vector< Strain >* strains = NULL;
vector< UnitigIdStrandPos >* nodeIdToUnitigId;
bool keepNA = false;
double phenotypeThreshold = 0.0;

void populateParser (Tool *tool) {
  // We add some custom arguments for command line interface
  tool->getParser()->push_front (new OptionOneParam(STR_NBCORES, "Number of cores to use", false, "4"));
  tool->getParser()->push_front (new OptionNoParam (STR_KEEP_NA, "Keep strains with phenotype NA.",  false));
  tool->getParser()->push_front (new OptionOneParam (STR_OUTPUT, "Path to the folder where the final and temporary files will be stored.",  false, "output"));
  tool->getParser()->push_front (new OptionOneParam (STR_KSKMER_SIZE, "K-mer size.",  false, "31"));
  tool->getParser()->push_front (new OptionOneParam (STR_PATH_TO_FASTA_FILES, "Path to the folder where the input fasta files are stored", true, "./"));
  tool->getParser()->push_front (new OptionOneParam (STR_STRAINS_FILE, "A text file describing the strains containing 3 columns: 1) ID of the strain; 2) Phenotype (a real number or NA); 3) Path to a multi-fasta file containing the sequences of the strain. This file needs a header. Check the sample_example folder or https://gitlab.com/leoisl/dbgwas/raw/master/sample_example/strains for an example.",  true));
  tool->getParser()->push_front (new OptionOneParam (STR_PHENOTYPE_THRESHOLD, "Phenotype threshold. Values <= than this are considered Phenotype 0, and > than this are considered Phenotype 1. Used only on the visualisation.",  false, "0.0"));
}
