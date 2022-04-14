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

#include "Utils.h"
#include "PhenoCounter.h"
#include "global.h"

using namespace std;

int getNbLinesInFile(const string &filename) {
  std::ifstream file;
  openFileForReading(filename, file);

  // Number of lines in the file
  int n = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');

  file.close();

  return n;
}

void fatalError (const string &message) {
  cerr << endl << endl << "[FATAL ERROR] " << message << endl << endl;
  cerr.flush();
  exit(1);
}

//strips all last "/" if exists in the parameter
string stripLastSlashIfExists (string path) {
  while(!path.empty() && path.back()=='/')
    path.pop_back();
  return path;
}


void openFileForReading(const string &filePath, ifstream &stream) {
  stream.open(filePath);
  if (!stream.is_open()) {
    stringstream ss;
    ss << "Error opening file " << filePath;
    fatalError(ss.str());
  }
}

/*
 * Next lines are for generating inputs for step 2
 */


void checkParametersMapping(Tool *tool) {
    //check the strains file
    keepNA = tool->getInput()->get(STR_KEEP_NA) != 0;
    string strainsFile = tool->getInput()->getStr(STR_STRAINS_FILE);
    checkStrainsFile(strainsFile);

}


//this function also populates strains if needed
void checkStrainsFile(const string &strainsFile) {

    vector<Strain> localStrains;
    bool header=true;
    ifstream input;
    openFileForReading(strainsFile, input);
    set<string> allIds;

    for(string line; getline( input, line ); )
    {
        //parse header
        if (header) {
            header=false;
            continue;
        }

        //ignore empty lines
        if (line.size()==0)
            continue;

        //create the strain
        stringstream ss;
        ss << line;
        string id, pheno, path;
        ss >> id >> pheno >> path;

        //check for duplicated IDs
        if (allIds.find(id)!=allIds.end()) {
            stringstream ss;
            ss << "Duplicated IDs in " << strainsFile << ": " << id << endl;
            fatalError(ss.str());
        }
        allIds.insert(id);

        //check for disallowed phenotypes
        bool phenoIsNumber;
        double phenoAsNumber;
        std::tie(phenoIsNumber, phenoAsNumber) = is_number(pheno);

        //allowed phenotypes are only "NA" or numbers between 0 and 1
        if (!phenoIsNumber && pheno!="NA") {
            stringstream ss;
            ss << "Phenotype not allowed: " << pheno << " . The only allowed values for phenotypes are real numbers or NA." << endl;
            fatalError(ss.str());
        }

        //check if the path is ok
        ifstream file;
        openFileForReading(path, file);
        if (!file.is_open()) {
            stringstream ss;
            ss << "Error opening file " << path << " in " << strainsFile << endl;
            fatalError(ss.str());
        }
        file.close();


        //add the strain if it is different from NA
        if (pheno=="NA" && keepNA==false) {
            cerr << "[WARNING] Skipping strain " << id << " because its phenotype is NA and " << STR_KEEP_NA << " is not set." << endl;
        }else {
            Strain strain(id, pheno, path);
            localStrains.push_back(strain);
        }
    }
    input.close();

    //in the end, check if strain is null. If it is, populate it
    if (strains==NULL)
        strains = new vector<Strain>(localStrains);
}


//tries to parse s, and returns a pair<bool, double>
//the first value indicates if s was successfully parsed into a double
//the second value indicates the double (it is only valid if the first is true)
tuple<bool, float> is_number(const std::string& s) {
    float number;
    try {
        number = std::stod(s);
    }
    catch(...) {
        return make_tuple(false, number);
    }
    return make_tuple(true, number);
}

void openFileForWriting(const string &filePath, ofstream &stream) {
    stream.open(filePath);
    if (!stream.is_open()) {
        stringstream ss;
        ss << "Error opening file " << filePath;
        fatalError(ss.str());
    }
}


void createFolder(const string &path) {
    boost::filesystem::path folder(path.c_str());

    if (boost::filesystem::exists(folder))
        return;

    if (!boost::filesystem::create_directories(folder)) {
        stringstream ss;
        ss << "Could not create dir " << path << " - unknown reasons...";
        fatalError(ss.str());
    }
}


//Read all strings in the readsFile file and return them as a vector of strings
vector<string> getVectorStringFromFile(const string &readsFile) {
    vector<string> allReadFilesNames;
    string tempStr;

    ifstream readsFileStream;
    openFileForReading(readsFile, readsFileStream);
    while (getline(readsFileStream, tempStr)) {
        if (tempStr.size() > 0)
            allReadFilesNames.push_back(tempStr);
    }
    readsFileStream.close();

    return allReadFilesNames;
}


void Strain::createPhenotypeCounter(const string &filePath, vector< Strain >* strains) {
    PhenoCounter phenoCounter;
    for (const auto &strain : (*strains))
        phenoCounter.add(strain.phenotype, 1);
    //serialize phenoCounter
    ofstream phenoCounterFile;
    {
        openFileForWriting(filePath, phenoCounterFile);
        boost::archive::text_oarchive boostOutputArchive(phenoCounterFile);
        //serialization itself
        boostOutputArchive & phenoCounter;
    } //boostOutputArchive is closed on destruction
    phenoCounterFile.close();
}

