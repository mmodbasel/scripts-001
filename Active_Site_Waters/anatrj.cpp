/* 
 * File:   anatrj.cpp
 * Author: Martin Smiesko
 *
 * Created on January 2, 2017, 4:00 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "definitions.h"
#include "structures.h"

#define REF_ATOM_TYPE   "CA"
#define MAX_AS_RES_NUM  100

char    atom_giveDonAcc(ProteinAtom * pa2) {
        char    pa[5];
        sscanf(pa2->getAtomType(), "%s", pa);
        pa[4] = 0;
        if (pa2->getProtonNum() != CARBON) {
                // backbone
                if(strcmp(pa, "N") == 0) return DONOR;
                if(strcmp(pa, "O") == 0) return ACCEPTOR;
                if(strcmp(pa, "OXT") == 0) return ACCEPTOR;
                // GLN - GLU
                if(strcmp(pa, "OE1") == 0) return ACCEPTOR;
                if(strcmp(pa, "OE2") == 0) return ACCEPTOR;
                if(strcmp(pa, "NE2") == 0) return DONOR;
                // ASN - ASP
                if(strcmp(pa, "OD1") == 0) return ACCEPTOR;
                if(strcmp(pa, "OD2") == 0) return ACCEPTOR;
                if(strcmp(pa, "ND2") == 0) return DONOR;
                // THR, SER, CYS, TYR
                if(strcmp(pa, "OG1") == 0) return DON_ACC;
                if(strcmp(pa, "OG") == 0) return DON_ACC;
                if(strcmp(pa, "SG") == 0) return DON_ACC;
                if(strcmp(pa, "OH") == 0) return DON_ACC;
                // MET
                if(strcmp(pa, "SD") == 0) return ACCEPTOR;
                // TRP, LYS, ARG
                if(strcmp(pa, "NE1") == 0) return DONOR;
                if(strcmp(pa, "NZ") == 0) return DONOR;
                if(strcmp(pa, "NE") == 0) return DONOR;
                if(strcmp(pa, "NH1") == 0) return DONOR;
                if(strcmp(pa, "NH2") == 0) return DONOR;
                // HIS - always unsure due to protonation
                if(strcmp(pa, "ND1") == 0) return DON_ACC_UNSURE;
                if(strcmp(pa, "NE2") == 0) return DON_ACC_UNSURE;
        }
        return NO;
}

#define SIG_LENGTH       16
class RawLine{
protected:
        RawLine * previous;
        int     atomnr;
        char    signature[SIG_LENGTH];
        char    letter;
        double  occ;
public:
        RawLine *       getPrevious() { return previous; }
        int             getAtomNr() { return atomnr; }
        double          getOccupancy() { return occ; }
        char            getLetter() { return letter; }
        char *          getSignature() { return signature; }
        
        RawLine(int nr, char * sig, double o, RawLine * rwln) {
                atomnr = nr;
                strncpy(signature, sig, SIG_LENGTH-1);
                // termination character
                signature[SIG_LENGTH-1] = 0;
                // save letter separately
                letter = signature[5];
                // and make it the same : for all
                signature[5] = ':';
                occ = o;
                previous = rwln;
        }
};

class   Ignore {
protected:        
        Ignore *        previous;
        int             atomnr;
public:
        int             getNr() { return atomnr; }
        Ignore *        getPrevious() { return previous; }
        Ignore(int n, Ignore * p) { previous = p; atomnr = n; }
};

class   DBref {
protected:
        DBref *         previous;
        char            chain;
        char            dbrefStr[MAX_DBREF_LENGTH];
public:
        char            getChain() { return chain; }
        char *          getDBrefStr() { return dbrefStr; }
        DBref *         getPrevious() { return previous; }
        
        void            addDBrefStr(char * v) {
                                // if not the first one, insert a spacer
                                if(strlen(dbrefStr) != 0) strcat(dbrefStr, "_");
                                strcat(dbrefStr, v);

                                if(strlen(dbrefStr) > MAX_DBREF_LENGTH) {
                                        printf("Max length of a DBREF record exceeded!\n");
                                        exit(1);
                                }
                        }
        DBref(char c, char * dbr, DBref * p) {
                chain = c;
                previous = p;
                strcpy(dbrefStr, dbr);
        }
};

int     shallBeIgnored(int n, Ignore * ig) {
        while(ig) {
                if(n == ig->getNr()) return YES;
                ig = ig->getPrevious();
        }
        return NO;
}

void    readProteinPDB(char * filename, Residue * &lastResidue, Molecule * &lastMolecule, Water * &lastWater, RawAtom * &lastRawAtom, char readH = NO) {
	int             atomCounter = 0;
	int     	residueCounter = 0;
        int             waterCounter = 0;
        int             moleculeCounter = 0;
	int             resNum;
	int             lastNum  = -1;
	int             atomNum;
        unsigned char   N;
	char		letter = '-';	// amino acid one letter code + extended ones
	char		lineBuffer[255];
	char		resName[5];
	char		lastName[5];
	char		atomType[5];
        char		ts[10];         // temporary string
        char            element[3];
        char            chain;
        char            lastChain;
        char            altLoc;
        double          occupancy;
        double          Bfactor;
	Point		pP;		// pP holding protein atom's coordinates
        int             charge;
        Residue *       curResidue = NULL;
        Molecule *      curMolecule = NULL;
        Water *         curWater = NULL;
        char            flag;
        int             position;

        ProteinAtom *   pa;
        LigandAtom *    la;
        RawAtom *       ra;
        RawAtom *       ra2;
        
        lastResidue = NULL; 
        lastMolecule = NULL;
        lastWater = NULL;
        
        RawLine *       rl = NULL;
        RawLine *       rlLast = NULL;
        RawLine *       rl2 = NULL;

        Ignore *        ig = NULL;
        
        FILE *	input = fopen(filename, "rt");
	
	// open file for reading
	if (input == NULL) { printf("Error! Problem opening file %s.\n", filename); exit(1); }

        // P R L I M I N A R Y   C H E C K   F O R   O C C U P A N C I E S
        // read signatures
        while (fgets(lineBuffer, sizeof(lineBuffer), input) != NULL  &&  strstr(lineBuffer, "ENDMDL") != lineBuffer)
	{
                     if(strstr(lineBuffer, "ATOM") == lineBuffer) flag = ATOM_LINE;
                else if(strstr(lineBuffer, "HETATM") == lineBuffer) flag = HETATM_LINE;
                else flag = NO;
                     
                if(flag != NO) {
                        // atom number
                        if(flag == ATOM_LINE) { strncpy(ts, lineBuffer+4, 7); ts[7] = 0; }
                        else if(flag == HETATM_LINE) { strncpy(ts, lineBuffer+6, 5); ts[5] = 0; }
                        if(sscanf(ts, "%d", &atomNum) != 1) { printf("Error1 (atom num) reading PDB line : %s\n", lineBuffer); exit(1); };

                        // occupancy
                        strncpy(ts, lineBuffer+56, 4);
                        ts[4] = 0;
                        if(sscanf(ts, "%lf", &occupancy) != 1) { printf("Error8 (occupancy) reading PDB line : %s\n", lineBuffer); exit(1); }

                        // signature starts at offset 11 and is 15 chars long
                        rl = new RawLine(atomNum, lineBuffer+11, occupancy, rl);
                }
        }
        // backup the pointer
        rlLast = rl;
        
        // loop over signatures and find duplicates
        printf("List of ignored atom numbers (because of occupancy):\n");
        while(rl) {
                // rl2 = rlLast;
                rl2 = rl->getPrevious();
                while(rl2) {
                        // if not the same atom but the signature is the same
                        if(rl != rl2  &&  strcmp(rl->getSignature(), rl2->getSignature()) == 0) {
                                // keep rl
                                if(rl->getOccupancy() > rl2->getOccupancy()  ||  (rl->getOccupancy() == rl2->getOccupancy()  &&  rl->getLetter() < rl2->getLetter())) {
                                        ig = new Ignore(rl2->getAtomNr(), ig);
                                        printf("%d, ", rl2->getAtomNr());
                                }
                                // keep rl2
                                else if(rl->getOccupancy() < rl2->getOccupancy()  ||  (rl->getOccupancy() == rl2->getOccupancy()  &&  rl->getLetter() > rl2->getLetter())) {
                                        ig = new Ignore(rl->getAtomNr(), ig);
                                        printf("%d, ", rl->getAtomNr());
                                }
                        }
                        rl2 = rl2->getPrevious();
                }
                rl = rl->getPrevious();
        }
        printf("\nEnd.\n");
        
        // free the memory
        rl = rlLast;
        while(rl) {
                rl2 = rl;
                rl = rl->getPrevious();
                delete rl2;
        }
        
        // R E A D   F I N A L   D A T A
        rewind(input);
	// clear resName
	strcpy(lastName, "");
        lastChain = '-';
	
        // START OF READING COORDINATES
        // READ ATOM & HETATM
        while (fgets(lineBuffer, sizeof(lineBuffer), input) != NULL  &&  strstr(lineBuffer, "ENDMDL") != lineBuffer)
	{
                // STRUCTURAL DATA LINES
                     if(strstr(lineBuffer, "ATOM") == lineBuffer) flag = ATOM_LINE;
                else if(strstr(lineBuffer, "HETATM") == lineBuffer) flag = HETATM_LINE;
                else flag = NO;
                
                if(flag != NO) {
                        // atom number
                        if(flag == ATOM_LINE) { strncpy(ts, lineBuffer+4, 7); ts[7] = 0; }
                        else if(flag == HETATM_LINE) { strncpy(ts, lineBuffer+6, 5); ts[5] = 0; }
                        
                        if(sscanf(ts, "%d", &atomNum) != 1) { printf("Error1 (atom num) reading PDB line : %s\n", lineBuffer); exit(1); };
                        if(atomNum < 0) { printf("Error1 atom number smaller than zero at line : %s\n", lineBuffer); exit(1); };
                        
                        // atom type
                        // strncpy(atomType, lineBuffer+12, 4);
                        // atomType[4] = 0;
                        strncpy(ts, lineBuffer+12, 4);
                        ts[4] = 0;
                        if(sscanf(ts, "%s", atomType) != 1) { printf("Error2 (atom type) reading PDB line : %s\n", lineBuffer); exit(1); };

                        // alternative orientation
                        altLoc = lineBuffer[16];
                        
                        // res name
                        strncpy(resName, lineBuffer+17, 3);
                        resName[3] = 0;
                        // f(sscanf(ts, "%s", resName) != 1) { printf("Error3 (residue name) reading PDB line : %s\n", lineBuffer); exit(1); };
                        
                        // res num
                        strncpy(ts, lineBuffer+22, 4);
                        ts[4] = 0;
                        // if(sscanf(ts, "%d", &resNum) != 1) { printf("Error4 (residue number) reading PDB line : %s\n", lineBuffer); exit(1); }
                        if(sscanf(ts, "%d", &resNum) != 1) {
                                // printf("Could not read res num as a number, trying decoding...");
                                
                                     if(ts[0] == 'A') resNum = 10000;
                                else if(ts[0] >  'A') resNum = 10000 + (ts[0]-'A')*36*36*36;
                                else { printf("PROBLEM 1\n"); exit(2); };     
                                
                                if(ts[1] >= 'A') resNum = resNum + (ts[1]-'@'+10)*36*36;
                                else             resNum = resNum + (ts[1]-'0')*36*36;
                                     
                                if(ts[2] >= 'A') resNum = resNum + (ts[2]-'@'+10)*36;
                                else             resNum = resNum + (ts[2]-'0')*36;

                                if(ts[3] >= 'A') resNum = resNum + (ts[3]-'@'+10);
                                else             resNum = resNum + (ts[3]-'0');

                                // printf("Warning: ResNum transformed from %s to %d\n", ts, resNum);
                                // "Error4 (residue number) reading PDB line : %s\n", lineBuffer);
                                // exit(1);
                        }
                        
                        // chain ID : 1 char
                        chain = lineBuffer[21];

                        // X
                        strncpy(ts, lineBuffer+30, 8);
                        ts[8] = 0;
                        if(sscanf(ts, "%lf", &pP.x) != 1) { printf("Error5 (x-coord) reading PDB line : %s\n", lineBuffer); exit(1); }
                        // Y
                        strncpy(ts, lineBuffer+38, 8);
                        ts[8] = 0;
                        if(sscanf(ts, "%lf", &pP.y) != 1) { printf("Error6 (y-coord) reading PDB line : %s\n", lineBuffer); exit(1); }
                        // Z
                        strncpy(ts, lineBuffer+46, 8);
                        ts[8] = 0;
                        if(sscanf(ts, "%lf", &pP.z) != 1) { printf("Error7 (z-coord) reading PDB line : %s\n", lineBuffer); exit(1); }
                        
                        // occupancy
                        strncpy(ts, lineBuffer+56, 4);
                        ts[4] = 0;
                        if(sscanf(ts, "%lf", &occupancy) != 1) { printf("Error8 (occupancy) reading PDB line : %s\n", lineBuffer); exit(1); }
                        // if(occupancy < 0.0  ||  occupancy > 1.0) { printf("Error8 (occupancy) reading PDB line : %s\n", lineBuffer); exit(1); }
                        if(occupancy < 0.0  ||  occupancy > 1.0) { printf("Warning8 (occupancy) reading PDB line : %s\n", lineBuffer); occupancy = 1.0; }
                        
                        // B-factor
                        strncpy(ts, lineBuffer+60, 5);
                        ts[5] = 0;
                        if(sscanf(ts, "%lf", &Bfactor) != 1) { printf("Error9 (B-factor) reading PDB line : %s\n", lineBuffer); exit(1); }
                        
                        // element
                        strncpy(element, lineBuffer + 76, 2);
                        element[2] = '\0';
                        toUpperCase(element);
                        N = giveProtonNumber(element);
                        
                        // charge   e.g. " N1+"
                        strncpy(ts, lineBuffer + 78, 2);
                        ts[1] = '\0';
                        // if(sscanf(ts, "%hd", &charge) != 1) charge = 0; // { printf("Error11 (charge) reading PDB line : %s\n", lineBuffer); exit(1); }
                        if(sscanf(ts, "%d", &charge) != 1) charge = 0; // { printf("Error11 (charge) reading PDB line : %s\n", lineBuffer); exit(1); }
                        else {
                                if(lineBuffer[79] == '-') charge = -charge;
                                // printf("Non-zero charge %d found at line %s", charge, lineBuffer);
                        }
                        // printf("%7d %-5s %3s %4hd%8.3f%8.3f%8.3f\n", atomNum, atomType, resName, resNum, pP.x, pP.y, pP.z);
			

                        if((N == HYDROGEN  &&  readH == NO)  ||  shallBeIgnored(atomNum, ig) == YES) {}
                        else {
                                // if the name or number is different, make a new prot.residue
                                if(strcmp(lastName, resName) != 0  ||  lastNum != resNum  ||  lastChain != chain ){
                                        // set new name and number
                                        strcpy(lastName, resName);
                                        lastNum = resNum;
                                        lastChain = chain;

                                        if(flag == ATOM_LINE) {
                                                // encode into one letter amino acid code
                                                letter = residueName_C2L(resName);
                                                // create a new residue
                                                lastResidue = new Residue(resNum, letter, resName, chain, giveRotBonds(letter), lastResidue);
                                                // increase counter as there is a new residue
                                                residueCounter++;
                                                // printf("%s ", resName);
                                        }
                                        else if(flag == HETATM_LINE) {
                                                if(strcmp(resName, "HOH") == 0  ||  strcmp(resName, "WAT") == 0  ||  strcmp(resName, "DOD") == 0  ||  strcmp(resName, "T3P") == 0) {
                                                        lastWater = new Water(resNum, '~', resName, chain, lastWater);
                                                        waterCounter++;
                                                        letter = '~';
                                                        // printf("%s ", resName);
                                                }
                                                else {
                                                        lastMolecule = new Molecule(resNum, 'X', resName, chain, lastMolecule);
                                                        moleculeCounter++;
                                                        letter = 'X';
                                                        printf("new HET residue: %c%s %c %d\n", altLoc, resName, chain, resNum);
                                                }
                                        }
                                }

                                // add atom
                                if(flag == ATOM_LINE) {
                                        // add atom to the residue
                                        lastResidue->addAtom(atomNum, atomType, altLoc, pP, occupancy, Bfactor, N, charge);
                                }
                                else if(flag == HETATM_LINE) {
                                        // add water
                                        if(letter == '~') {
                                                // add atom to the water
                                                lastWater->addAtom(atomNum, atomType, altLoc, pP, occupancy, Bfactor, N, charge);
                                        }
                                        // add molecule
                                        else if (letter == 'X') {
                                                // add atom to the molecule
                                                lastMolecule->addAtom(atomNum, atomType, altLoc, pP, occupancy, Bfactor, N, charge);
                                        }
                                        else {
                                                printf("Unknown structural element.\n");
                                                exit(1);
                                        }
                                }
                                // increase the counter
                                atomCounter++;
                        }
		}
        }
        // END OF READING COORDINATES
        
        if(atomCounter == 0) {
                printf("No useful atoms were read... Strange.\n");
                exit(1);
        }
        if(residueCounter == 0) {
                printf("No protein atoms were read... Exit.\n");
                exit(1);
        }

        
        ProteinAtom *   pa2;
        LigandAtom *    la2;
        
        // initialize last pointer
        lastRawAtom = lastResidue->getLastAtom();
        
        // chain all atoms together
        curResidue = lastResidue;
        ra = NULL;
        while(curResidue) {
                pa = curResidue->getLastAtom();
                while(pa) {
                        // remember the last one for linking to other elements
                        ra = (RawAtom *) pa;
                        // if pointer to the previous is null, link to the previous one
                        if(pa->getPrevious() == NULL  &&  curResidue->getPrevious() != NULL) {
                                // ra->setPrevious((RawAtom *) curResidue->getPrevious()->getLastAtom());
                                pa2 = curResidue->getPrevious()->getLastAtom();
                                ra2 = (RawAtom *) pa2;
                                ra->setPrevious(ra2);
                        }
                        // set species
                        pa->setSpecies(PROTEIN);
                        
                        ra->setDAC(atom_giveDonAcc(pa));
                        // printf("->%s<- is DAC=%d\n", pa->getAtomType(), ra->getDAC());
                        
                        pa = pa->getPrevious();
                }
                curResidue = curResidue->getPrevious();
        }

        curMolecule = lastMolecule;
        if(curMolecule != NULL) ra->setPrevious((RawAtom *) curMolecule->getLastAtom());
        while(curMolecule) {
                la = curMolecule->getLastAtom();
                while(la) {
                        // if pointer to the previous is null, link to the previous one
                        ra = (RawAtom *) la;
                        // if(la->getPrevious() == NULL  &&  curMolecule->getPrevious() != NULL) la->setPrevious((RawAtom *) curMolecule->getPrevious()->getLastAtom());
                        if(la->getPrevious() == NULL  &&  curMolecule->getPrevious() != NULL) {
                                // ra->setPrevious((RawAtom *) curMolecule->getPrevious()->getLastAtom());
                                la2 = curMolecule->getPrevious()->getLastAtom();
                                ra2 = (RawAtom *) la2;
                                ra->setPrevious(ra2);
                        }

                        // remember the last one for linking to other elements
                        la = la->getPrevious();
                }
                curMolecule = curMolecule->getPrevious();
        }
        
        curWater = lastWater;
        if(curWater != NULL) ra->setPrevious((RawAtom *) curWater->getLastAtom());
        while(curWater) {
                pa = curWater->getLastAtom();
                while(pa) {
                        // remember the last one for linking to other elements
                        ra = (RawAtom *) pa;
                        // if pointer to the previous is null, link to the previous one
                        if(pa->getPrevious() == NULL  &&  curWater->getPrevious() != NULL) {
                                // ra->setPrevious((RawAtom *) curResidue->getPrevious()->getLastAtom());
                                pa2 = curWater->getPrevious()->getLastAtom();
                                ra2 = (RawAtom *) pa2;
                                ra->setPrevious(ra2);
                        }
                        // set species
                        pa->setSpecies(WATER);
                        
                        pa = pa->getPrevious();
                }
                curWater = curWater->getPrevious();
        }        
        
        curMolecule = lastMolecule;
        
        // R E A D   L I N K
        printf("Reading LINK records...\n");
	char		resName2[5];
	char		atomType2[5];
        char            chain2;
        int             resNum2;
        char            altLoc2;
        RawAtom *       raHelp;
        StructuralUnit *su;
        
        rewind(input);
        while (fgets(lineBuffer, sizeof(lineBuffer), input) != NULL) {
                
                if(strstr(lineBuffer, "LINK  ") == lineBuffer) {
                        // A T O M   1
                        // atom type
                        strncpy(atomType, lineBuffer+12, 4);
                        atomType[4] = 0;
                        // alt. location
                        altLoc = lineBuffer[16];
                        // residue name
                        strncpy(resName, lineBuffer+17, 3);
                        resName[3] = 0;
                        // chain
                        chain = lineBuffer[21];
                        // residue number
                        strncpy(ts, lineBuffer+22, 4);
                        ts[4] = 0;
                        sscanf(ts, "%d", &resNum);

                        // A T O M   2
                        // atom type
                        strncpy(atomType2, lineBuffer+42, 4);
                        atomType2[4] = 0;
                        // alt. location
                        altLoc2 = lineBuffer[46];
                        // residue name
                        strncpy(resName2, lineBuffer+47, 3);
                        resName2[3] = 0;
                        // chain
                        chain2 = lineBuffer[51];
                        // residue number
                        strncpy(ts, lineBuffer+52, 4);
                        ts[4] = 0;
                        sscanf(ts, "%d", &resNum2);
                        
                        printf("Link : %s %c %s %c %d  --- %s %c %s %c %d\n", atomType, altLoc, resName, chain, resNum, atomType2, altLoc2, resName2, chain2, resNum2);
                        
                        
                        
                        // loop over all residues and molecules to find the corresponding atom
                        raHelp = (RawAtom *) lastResidue->getLastAtom();
                        ra = ra2 = NULL;
                        while (raHelp  &&  (ra == NULL  ||  ra2 == NULL)) {
                                su = raHelp->getResidue();
                                
                                // check resNum, chain, altLoc and strings
                                if(su->getUnitNumber() == resNum  &&
                                   su->getChain() == chain  &&
                                   // raHelp->getAltLocation() == altLoc  &&
                                   strcmp(su->getResidueName(), resName) == 0  &&
                                   strcmp(raHelp->getAtomType(), atomType) == 0) { 
                                        ra = raHelp;
                                        // printf("a1 found\n");
                                }
                                // check resNum, chain, altLoc and strings
                                if(su->getUnitNumber() == resNum2  &&
                                   su->getChain() == chain2  &&
                                   // raHelp->getAltLocation() == altLoc2  &&
                                   strcmp(su->getResidueName(), resName2) == 0  &&
                                   strcmp(raHelp->getAtomType(), atomType2) == 0) { 
                                        ra2 = raHelp;
                                        // printf("a2 found\n");
                                }
                                
                                raHelp = raHelp->getPrevious();
                        }

                        
                        
                        
                        if(ra != NULL  &&  ra2 != NULL) {
                                     if( ra->getSpecies() == PROTEIN) pa = (ProteinAtom *) ra;
                                else if(ra2->getSpecies() == PROTEIN) pa = (ProteinAtom *) ra2;
                                
                                if(ra->getSpecies() == PROTEIN  &&  ra2->getSpecies() == PROTEIN) {
                                        printf("Link 2 protein atoms %s %d & %s %d\n", ra->getAtomType(), ra->getNumber(), ra2->getAtomType(), ra2->getNumber());
                                        // connect the two atoms
                                        pa->setNeighbor((ProteinAtom *) ra2);
                                }
                                else if(ra->getSpecies() == LIGAND  &&  ra2->getSpecies() == LIGAND) {
                                        if(ra->getResidue() == ra2->getResidue()) {
                                                // printf("Link 2 atoms from the same ligand %d & %d (%s)\n", ra->getNumber(), ra2->getNumber(), ra->getResidue()->getResidueName());
                                        }
                                        else {
                                                printf("Link atoms %d & %d of two ligands (%s & %s)\n", ra->getNumber(), ra2->getNumber(), ra->getResidue()->getResidueName(), ra2->getResidue()->getResidueName());
                                                ((LigandAtom *) ra)->setLinked((LigandAtom *) ra2);
                                                ((LigandAtom *) ra2)->setLinked((LigandAtom *) ra);
                                        }
                                        // connect the two atoms
                                        
                                        // ((LigandAtom *) ra)->setNeighbor((LigandAtom *) ra2);
                                        // ((LigandAtom *) ra)->setBond(((LigandAtom *) ra)->getNeighborCount()-1, SINGLE);
                                        ((LigandAtom *) ra)->setNeighborAndBond((LigandAtom *) ra2, SINGLE);
                                }

                                else if((ra->getSpecies() == PROTEIN  &&  ra2->getSpecies() == LIGAND)  ||  (ra2->getSpecies() == PROTEIN  &&  ra->getSpecies() == LIGAND)) {
                                        printf("Ligand-protein link at atoms %s %d (%s) & %s %d (%s)\n",
                                                ra->getAtomType(), ra->getNumber(), ra->getResidue()->getResidueName(), ra2->getAtomType(), ra2->getNumber(), ra2->getResidue()->getResidueName());
                                        // connect the two atoms
                                        // pa->setNeighbor((ProteinAtom *) ra2);

                                        if( ra->getSpecies() == LIGAND) { 
                                                ((Molecule *)  ra->getResidue())->setProteinBound(YES);
                                                ((LigandAtom *) ra)->setLinked((LigandAtom *) ra2);
                                        }
                                        if(ra2->getSpecies() == LIGAND) { 
                                                ((Molecule *) ra2->getResidue())->setProteinBound(YES);
                                                ((LigandAtom *) ra2)->setLinked((LigandAtom *) ra);
                                        }
                                }

                                else if((ra->getSpecies() == METAL  &&  ra2->getSpecies() == LIGAND)  ||  (ra2->getSpecies() == METAL  &&  ra->getSpecies() == LIGAND)) {
                                        // within one residue
                                        if(ra->getResidue() == ra2->getResidue()) {
                                                printf("Metal bond within one ligand %s %d & %s %d (%s)\n",
                                                        ra->getAtomType(), ra->getNumber(), ra2->getAtomType(), ra2->getNumber(), ra2->getResidue()->getResidueName());
                                                // connect the two atoms
                                                ((LigandAtom *) pa)->setNeighbor((LigandAtom *) ra2);

                                                ((Molecule *) ra->getResidue())->setMetalloligand(YES);
                                        }
                                        // between two residues
                                        else {
                                                printf("Metal bond to a cofactor at atoms %s %d (%s) & %s %d (%s)\n",
                                                        ra->getAtomType(), ra->getNumber(), ra->getResidue()->getResidueName(), ra2->getAtomType(), ra2->getNumber(), ra2->getResidue()->getResidueName());
                                                // connect the two atoms
                                                if( ra->getSpecies() == METAL) ((LigandAtom *) ra2)->setMetalNeighbor((LigandAtom *)ra);
                                                if(ra2->getSpecies() == METAL) ((LigandAtom *) ra)->setMetalNeighbor((LigandAtom *)ra2);
                                        }
                                }
                                else if((ra->getSpecies() == METAL  &&  ra2->getSpecies() == PROTEIN)  ||  (ra2->getSpecies() == METAL  &&  ra->getSpecies() == PROTEIN)) {
                                        printf("Metal-protein bond at atoms %s %d & %s %d\n",
                                                ra->getAtomType(), ra->getNumber(), ra2->getAtomType(), ra2->getNumber());
                                        if( ra->getSpecies() == METAL) ((Molecule *) ra->getResidue())->setProteinBound(YES);
                                        if(ra2->getSpecies() == METAL) ((Molecule *) ra2->getResidue())->setProteinBound(YES);
                                }

                                else if((ra->getSpecies() == METAL  &&  ra2->getSpecies() == WATER)  ||  (ra2->getSpecies() == METAL  &&  ra->getSpecies() == WATER)) {
                                        printf("Water bound to metal - atoms: %s %d & %s %d\n",
                                                ra->getAtomType(), ra->getNumber(), ra2->getAtomType(), ra2->getNumber());
                                        // to be done
                                        // if( ra->getSpecies() == WATER) ((Water *) ra->getResidue())->addBondToMetal();  // CONECT is duplicate, it is enough to scan once
                                        // if(ra2->getSpecies() == WATER) ((Water *) ra2->getResidue())->addBondToMetal();
                                }
                                else if(ra->getSpecies() == METAL  &&  ra2->getSpecies() == METAL) {
                                        printf("metal - metal bond: %s %d & %s %d\n",
                                                ra->getAtomType(), ra->getNumber(), ra2->getAtomType(), ra2->getNumber());
                                        ((LigandAtom *) ra2)->setMetalNeighbor((LigandAtom *)ra);
                                        ((LigandAtom *) ra)->setMetalNeighbor((LigandAtom *)ra2);
                                }
                                else {
                                        printf("WARNING: IGNORING LINK at atoms %d %s %s & %d %s %s\n",
                                                ra->getNumber(), ra->getAtomType(), ra->getResidue()->getResidueName(), ra2->getNumber(), ra2->getAtomType(), ra2->getResidue()->getResidueName());

                                        printf("species1=%d  species2=%d\n", ra->getSpecies(), ra2->getSpecies());
                                        // exit(1);
                                }
                        }
                        else {
                                printf("Could not find corresponding atoms - SKIPPED\n");
                                // exit(1);
                        }
                        // printf("LINK done\n");
                }
        }

        
        // R E A D   C O N E C T
        printf("Reading CONECT records...\n");
        rewind(input);
        while (fgets(lineBuffer, sizeof(lineBuffer), input) != NULL) {
                
                // if(strstr(lineBuffer, "LINK") == lineBuffer) {}
                
                if(strstr(lineBuffer, "CONECT") == lineBuffer) {
                        strncpy(ts, lineBuffer+6, 5);
                        ts[5] = 0;
                        sscanf(ts, "%d", &atomNum);

                        // loop over all residues and molecules to find the corresponding atom
                        ra2 = (RawAtom *) lastResidue->getLastAtom();
                        ra = NULL;
                        while (ra2  &&  ra == NULL) {
                                // if(ra->getNumber() == atomNum) printf("Connections from atom %d\n", atomNum);
                                if(ra2->getNumber() == atomNum) { 
                                        // printf("Connections from atom %d\n", atomNum);
                                        ra = ra2;
                                }
                                ra2 = ra2->getPrevious();
                        }
                        
                        if(ra != NULL) {
                        
                                pa = (ProteinAtom *) ra;        
                                position = 11;
                                do {
                                        strncpy(ts, lineBuffer + position, 5);
                                        ts[5] = 0;
                                        flag = sscanf(ts, "%d", &lastNum);
                                        // printf("->%s<- flag=%d : ", num, flag);
                                        if(flag == 1) {
                                                ra2 = (RawAtom *) lastResidue->getLastAtom();
                                                while (ra2) {
                                                        if(ra2->getNumber() == lastNum) {

                                                                if(ra->getSpecies() == PROTEIN  &&  ra2->getSpecies() == PROTEIN) {
                                                                        printf("Link 2 protein atoms %s %d & %s %d\n", ra->getAtomType(), ra->getNumber(), ra2->getAtomType(), ra2->getNumber());
                                                                        // connect the two atoms
                                                                        pa->setNeighbor((ProteinAtom *) ra2);
                                                                }
                                                                else if(ra->getSpecies() == LIGAND  &&  ra2->getSpecies() == LIGAND) {
                                                                        if(ra->getResidue() == ra2->getResidue()) {
                                                                                // printf("Link 2 atoms from the same ligand %d & %d (%s)\n", ra->getNumber(), ra2->getNumber(), ra->getResidue()->getResidueName());
                                                                        }
                                                                        else {
                                                                                printf("Link atoms %d & %d of two ligands (%s & %s)\n", ra->getNumber(), ra2->getNumber(), ra->getResidue()->getResidueName(), ra2->getResidue()->getResidueName());
                                                                                ((LigandAtom *) ra)->setLinked((LigandAtom *) ra2);
                                                                                ((LigandAtom *) ra2)->setLinked((LigandAtom *) ra);
                                                                        }
                                                                        // connect the two atoms
                                                                        ((LigandAtom *) ra)->setNeighborAndBond((LigandAtom *) ra2, SINGLE);
                                                                        ((LigandAtom *) ra2)->setNeighborAndBond((LigandAtom *) ra, SINGLE);
                                                                }

                                                                else if((ra->getSpecies() == PROTEIN  &&  ra2->getSpecies() == LIGAND)  ||  (ra2->getSpecies() == PROTEIN  &&  ra->getSpecies() == LIGAND)) {
                                                                        printf("Ligand-protein link at atoms %s %d (%s) & %s %d (%s)\n",
                                                                                ra->getAtomType(), ra->getNumber(), ra->getResidue()->getResidueName(), ra2->getAtomType(), ra2->getNumber(), ra2->getResidue()->getResidueName());
                                                                        // connect the two atoms
                                                                        pa->setNeighbor((ProteinAtom *) ra2);

                                                                        if( ra->getSpecies() == LIGAND) { 
                                                                                ((Molecule *)  ra->getResidue())->setProteinBound(YES);
                                                                                // ((Molecule *) ((LigandAtom *) ra)->getResidue())->setProteinBound(YES);
                                                                                ((LigandAtom *) ra)->setLinked((LigandAtom *) ra2);
                                                                        }
                                                                        if(ra2->getSpecies() == LIGAND) { 
                                                                                ((Molecule *) ra2->getResidue())->setProteinBound(YES);
                                                                                ((LigandAtom *) ra2)->setLinked((LigandAtom *) ra);
                                                                        }
                                                                }

                                                                else if((ra->getSpecies() == METAL  &&  ra2->getSpecies() == LIGAND)  ||  (ra2->getSpecies() == METAL  &&  ra->getSpecies() == LIGAND)) {
                                                                        // within one residue
                                                                        if(ra->getResidue() == ra2->getResidue()) {
                                                                                printf("Metal bond within one ligand %s %d & %s %d (%s)\n",
                                                                                        ra->getAtomType(), ra->getNumber(), ra2->getAtomType(), ra2->getNumber(), ra2->getResidue()->getResidueName());
                                                                                // connect the two atoms
                                                                                ((LigandAtom *) pa)->setNeighbor((LigandAtom *) ra2);

                                                                                ((Molecule *) ra->getResidue())->setMetalloligand(YES);
                                                                        }
                                                                        // between two residues
                                                                        else {
                                                                                printf("Metal bond to a cofactor at atoms %s %d (%s) & %s %d (%s)\n",
                                                                                       ra->getAtomType(), ra->getNumber(), ra->getResidue()->getResidueName(), ra2->getAtomType(), ra2->getNumber(), ra2->getResidue()->getResidueName());
                                                                                // connect the two atoms
                                                                                if( ra->getSpecies() == METAL) ((LigandAtom *) ra2)->setMetalNeighbor((LigandAtom *)ra);
                                                                                if(ra2->getSpecies() == METAL) ((LigandAtom *) ra)->setMetalNeighbor((LigandAtom *)ra2);
                                                                        }
                                                                }
                                                                else if((ra->getSpecies() == METAL  &&  ra2->getSpecies() == PROTEIN)  ||  (ra2->getSpecies() == METAL  &&  ra->getSpecies() == PROTEIN)) {
                                                                        printf("Metal-protein bond at atoms %s %d & %s %d\n",
                                                                               ra->getAtomType(), ra->getNumber(), ra2->getAtomType(), ra2->getNumber());
                                                                        if( ra->getSpecies() == METAL) ((Molecule *) ra->getResidue())->setProteinBound(YES);
                                                                        if(ra2->getSpecies() == METAL) ((Molecule *) ra2->getResidue())->setProteinBound(YES);
                                                                }

                                                                else if((ra->getSpecies() == METAL  &&  ra2->getSpecies() == WATER)  ||  (ra2->getSpecies() == METAL  &&  ra->getSpecies() == WATER)) {
                                                                        printf("Water bound to metal - atoms: %s %d & %s %d\n",
                                                                               ra->getAtomType(), ra->getNumber(), ra2->getAtomType(), ra2->getNumber());
                                                                        // to be done
                                                                        // if( ra->getSpecies() == WATER) ((Water *) ra->getResidue())->addBondToMetal();  // CONECT is duplicate, it is enough to scan once
                                                                        if(ra2->getSpecies() == WATER) ((Water *) ra2->getResidue())->addBondToMetal();
                                                                }
                                                                else if(ra->getSpecies() == METAL  &&  ra2->getSpecies() == METAL) {
                                                                        printf("metal - metal bond: %s %d & %s %d\n",
                                                                                ra->getAtomType(), ra->getNumber(), ra2->getAtomType(), ra2->getNumber());
                                                                        ((LigandAtom *) ra2)->setMetalNeighbor((LigandAtom *)ra);
                                                                        ((LigandAtom *) ra)->setMetalNeighbor((LigandAtom *)ra2);
                                                                }
                                                                else {
                                                                        printf("WARNING: IGNORING CONNECT at atoms %d %s %s & %d %s %s\n",
                                                                                ra->getNumber(), ra->getAtomType(), ra->getResidue()->getResidueName(), ra2->getNumber(), ra2->getAtomType(), ra2->getResidue()->getResidueName());
                                                                        
                                                                        printf("species1=%d  species2=%d\n", ra->getSpecies(), ra2->getSpecies());
                                                                        // exit(1);
                                                                }
                                                        }
                                                        ra2 = ra2->getPrevious();
                                                }
                                        }

                                        position = position + 5;
                                } while (flag == 1);
                        }
                }
        }

        
        // R E A D   C Y S T E I N E   B R I D G E S
        printf("Reading SSBOND records...\n");
        rewind(input);
        while (fgets(lineBuffer, sizeof(lineBuffer), input) != NULL) {
                
                if(strstr(lineBuffer, "SSBOND") == lineBuffer) {
                        strncpy(resName, lineBuffer+11, 3);
                        resName[3] = 0;

                        strncpy(resName2, lineBuffer+25, 3);
                        resName2[3] = 0;
                        
                        if(strcmp(resName, "CYS") == 0  &&  strcmp(resName2, "CYS") == 0) {
                                chain = lineBuffer[15];
                                chain2 = lineBuffer[29];

                                sscanf(lineBuffer+17, "%d", &resNum);
                                sscanf(lineBuffer+31, "%d", &resNum2);

                                curResidue = lastResidue;
                                while(curResidue) {
                                        if(curResidue->getID() == 'C') {
                                                if(curResidue->getChain() == chain  &&  curResidue->getUnitNumber() == resNum) {
                                                        ((Residue *) curResidue)->setID('c');
                                                        printf("Cysteine %d of chain %c marked as -S-S- bridged\n", resNum, chain);
                                                }
                                                if(curResidue->getChain() == chain2  &&  curResidue->getUnitNumber() == resNum2) {
                                                        ((Residue *) curResidue)->setID('c');
                                                        printf("Cysteine %d of chain %c marked as -S-S- bridged\n", resNum2, chain2);
                                                }
                                        }
                                        curResidue = curResidue->getPrevious();
                                }
                        }
                        else {
                                printf("\nWARNING: -S-S- bond for non-cysteine residues detected, but ignored!\n\n");
                        }
                }
        }
        
        
	printf("\n");

	
	// SUMMARY
	printf("Read PDB file %s : %d atoms, %d residues, %d molecules, %d waters\n", filename, atomCounter, residueCounter, moleculeCounter, waterCounter);
        
        
        // M I S S I N G   P R O T E I N   A T O M S
        rewind(input);
        flag = 0;
        while (fgets(lineBuffer, sizeof(lineBuffer), input) != NULL) {
                // find the header
                if(strstr(lineBuffer, "REMARK 470 MISSING ATOM") == lineBuffer) {
                        flag = 1;
                        printf("Missing protein atoms record found...\n");
                }
                // find second header
                else if(flag == 1  &&  strstr(lineBuffer, "REMARK 470   M RES CSSEQI  ATOMS") == lineBuffer) flag = 2;
                // data
                else if(flag == 2  &&  strstr(lineBuffer, "REMARK 470") == lineBuffer) {
                        strncpy(resName, lineBuffer + 15, 3);
                        chain = lineBuffer[19];
                        strncpy(ts, lineBuffer + 20, 4);
                        sscanf(ts, "%d", &resNum);
                        
                        printf("%s %c %d\n", resName, chain, resNum);
                        
                        ra = lastResidue->getLastAtom();
                        while(ra) {
                                // if the same and still useful, deusefulize it :)
                                if(ra->getResidue()->getChain() == chain  &&  strcmp(ra->getResidue()->getResidueName(), resName) == 0  &&  ra->getResidue()->getUnitNumber() == resNum  &&  ((Residue *) ra->getResidue())->getMissAtom() == 0) {
                                        ((Residue *) ra->getResidue())->setMissAtom();
                                        printf("Residue %s %c %d has missing atoms -> setting missing atom flag\n", ra->getResidue()->getResidueName(), ra->getResidue()->getChain(), ra->getResidue()->getUnitNumber());
                                }
                                ra = ra->getPrevious();
                        }
                }
        }
        
        
        
        // M I S S I N G   H E T E R O A T O M S
        rewind(input);
        flag = 0;
        while (fgets(lineBuffer, sizeof(lineBuffer), input) != NULL) {
                // find the header
                if(strstr(lineBuffer, "REMARK 610 MISSING HETEROATOM") == lineBuffer) {
                        flag = 1;
                        printf("Missing heteroatoms record found...\n");
                }
                // find second header
                else if(flag == 1  &&  strstr(lineBuffer, "REMARK 610   M RES C SSEQI") == lineBuffer) flag = 2;
                // data
                else if(flag == 2  &&  strstr(lineBuffer, "REMARK 610") == lineBuffer) {
                        strncpy(resName, lineBuffer + 15, 3);
                        chain = lineBuffer[19];
                        strncpy(ts, lineBuffer + 20, 5);
                        sscanf(ts, "%d", &resNum);
                        
                        printf("%s %c %d\n", resName, chain, resNum);
                        
                        ra = lastMolecule->getLastAtom();
                        while(ra) {
                                // if the same and still useful, deusefulize it :)
                                if(ra->getResidue()->getChain() == chain  &&  strcmp(ra->getResidue()->getResidueName(), resName) == 0  &&  ra->getResidue()->getUnitNumber() == resNum  &&  ((Molecule *) ra->getResidue())->getUseful() == YES) {
                                        ((Molecule *) ra->getResidue())->setUseful(NO);
                                        ((Molecule *) ra->getResidue())->setMissAtom(YES);
                                        printf("Molecule %s %c %d has missing atoms -> clearing useful flag\n", ra->getResidue()->getResidueName(), ra->getResidue()->getChain(), ra->getResidue()->getUnitNumber());
                                        
                                }
                                ra = ra->getPrevious();
                        }
                }
        }

        
        DBref *         db;
        DBref *         dbLast = NULL;
        int             i;
        
        // D B R E F
        rewind(input);
        flag = 0;
        while (fgets(lineBuffer, sizeof(lineBuffer), input) != NULL) {
                // if DBREF record
                if(strstr(lineBuffer, "DBREF") == lineBuffer) {
                        chain = lineBuffer[12];

                        strncpy(ts, lineBuffer + 33, 8);
                        for(i = 0; i < 8; i++) if(ts[i] == ' ') ts[i] = 0;
                        
                        flag = NO;
                        db = dbLast;
                        while(db) {
                                if(db->getChain() == chain) {
                                        flag = YES;
                                        db->addDBrefStr(ts);
                                }
                                db = db->getPrevious();
                        }
                        if(flag == NO) {
                                dbLast = new DBref(chain, ts, dbLast);
                        }
                }
        }
        // couple the dbrefs with mols
        curMolecule = lastMolecule;
        while(curMolecule) {
                db = dbLast;
                while(db) {
                        if(db->getChain() == curMolecule->getChain()) curMolecule->setDBrefStr(db->getDBrefStr());
                        db = db->getPrevious();
                }
                curMolecule = curMolecule->getPrevious();
        }

        
        
        // close the template input file
	if (fclose(input) != 0) { printf("Error! Closing file %s failed.\n", filename); exit(1); }
}

class WISR {
protected:
        int     rn;
        WISR *  previous;
public:
        WISR *  getPrevious() { return previous; }
        int     getRN() { return rn; }
        
        WISR(int RN, WISR * p) { rn = RN; previous = p; }
};		


int main(int argc, char *argv[]) {
	printf("anatrj ver. %4.2f by Martin Smiesko\n", VERSION);

        Molecule * 	ligand = NULL;
        Molecule * 	lastRawLigand = NULL;
        Water *         water = NULL;
        Water *         lastWater = NULL;
        Residue *       residue = NULL;
        Residue *       lastResidue = NULL;
        RawAtom *       lastRawAtom = NULL;
        LigandAtom *    la = NULL;
        ProteinAtom *   pa;
        int             i, j;
        // char            cmdln[255];
        int             ligNum = -1;
        char            ligNam[4];
        double          wisr = 10.0;     // water inclusion sphere radius
        char *          strng;
        WISR *          last_wis_res = NULL;
        WISR *          wis_res = NULL;
        Point           as_centroid = { 0, 0, 0 };
        int             flag;
        int             ligandInReach = NO;
        int             intruderInReach = NO;
        int             la_in_site = 0;
       
        // arg[0] - executable
        // arg[1] - PDB file
        // arg[2] - ResName
        // arg[3] - ResNumber
        // arg[4] - active site residue numbers
        // arg[5] - water inclusion sphere
        
        if(argc == 6  &&  strcmp(argv[1], "") != 0  &&  strcmp(argv[2], "") != 0  &&  strcmp(argv[3], "") != 0  &&  strcmp(argv[4], "") != 0  &&  strcmp(argv[5], "") != 0) {
                //  S T A R T   O F   L O A D I N G
                readProteinPDB(argv[1], lastResidue, lastRawLigand, lastWater, lastRawAtom);
                printf("\nProtein read.\n\n");
                // REF LIGAND NAME
                strncpy(ligNam, argv[2], 3);
                // REF LIGAND NUMBER
                sscanf(argv[3], "%d", &ligNum);
                // residue numbers defining the active site centroid                 
                // sscanf(argv[4], "%lf", &wisr);
                strng = strtok(argv[4], ",");
                i = j = 0;
                while(strng != NULL) {
                    if(sscanf(strng, "%d", &i) != 1) { printf("Problem reading WISRs\n"); exit(1); };
                    // printf("%d\n", i);
                    last_wis_res = new WISR(i, last_wis_res);
                    strng = strtok(NULL, ",");
                    j++;
                }

                // RADIUS
                sscanf(argv[5], "%lf", &wisr);
                
                printf("Analyzed ligand ID: %s %3d\n", ligNam, ligNum);
                // IsThereABigEnoughChain(residue);
                
                
                // getting centroid and check if all needed atoms present
                j = 0;
                wis_res = last_wis_res;
                while(wis_res) {
                    printf("Looking for atom %s of WISR Nr. %d\n", REF_ATOM_TYPE, wis_res->getRN());
                    
                    flag = NO;
                    // loop over protein residues
                    residue = lastResidue;
                    while(residue  &&  flag == NO) {
                        if(wis_res->getRN() == residue->getUnitNumber()) {
                            printf("Matching residue found\n");
                            
                            pa = residue->getLastAtom();
                            while(pa) {
                                // printf("Matching %s vs. %s\n", pa->getAtomType(), REF_ATOM_TYPE);
                                if(strcmp(pa->getAtomType(), REF_ATOM_TYPE) == 0) {
                                    printf("Matching atom found\n");
                                    as_centroid = as_centroid + pa->getCoord();
                                    flag = YES;
                                    j++;
                                }
                                pa = pa->getPrevious();
                            }
                        }
                        // next residue
                        residue = residue->getPrevious();
                    }
                    if(flag == NO) { printf("Atom needed for centroid not present in the PDB structure\n"); exit(1); };
                    
                    // next atom from the command line
                    wis_res = wis_res->getPrevious();
                }
                as_centroid.x = as_centroid.x / (double) j;
                as_centroid.y = as_centroid.y / (double) j;
                as_centroid.z = as_centroid.z / (double) j;
                printf("Centroid coords as follows (x, y, z) : %8.3f %8.3f %8.3f\n", as_centroid.x, as_centroid.y, as_centroid.z);
                


                // loop over all ligands
                ligand = lastRawLigand;
                while(ligand) {
                    // clear flag
                    flag = NO;
                    
                    printf("ligand %s %d\n", ligand->getResidueName(), ligand->getUnitNumber());
                    if(strcmp(ligand->getResidueName(), ligNam) == 0  &&  ligand->getUnitNumber() == ligNum) {
                        printf("Monitored ligand found\n");
                        flag = YES;
                    }
                    
                    // ignore heme
                    if(strcmp(ligand->getResidueName(), "HEM") != 0) 
                    {
                        la = ligand->getLastAtom();
                        while(la) {
                            if(giveDistance(la->getCoord(), as_centroid) <= wisr) {
                                if(flag == YES) {
                                    ligandInReach = YES;
                                    printf("Ligand in reach\n");
                                    if(la->getProtonNum() > HYDROGEN) la_in_site++;
                                }
                                else {
                                    intruderInReach = YES;
                                    printf("Intruder in reach\n");
                                }
                            }

                            la = la->getPrevious();
                        }
                    }
                    
                    ligand = ligand->getPrevious();
                }
                printf("\nLigand analysis done.\n\n");
                
                
                water = lastWater;
                i = 0;
                while(water) {
                    if(giveDistance(water->getO()->getCoord(), as_centroid) <= wisr) {
                        i++;
                    }
                    water = water->getPrevious();
                }
                printf("%s lir= %d laiac= %d iir= %d : %d waters within %8.3f Angstrom from the active site centroid.\n", argv[1], ligandInReach, la_in_site, intruderInReach, i, wisr);


                
		printf("anatrj.x finished gracefully.\n");
	}
	else { 
		printf("\nProgram to perform some analyses on the PDB structures.\n");
		printf("\nUsage: ./anatrj.x  PDBfile LigandName LigandNumber ActiveSiteResNumbers WaterInclusionSphereRadius\n");
                printf("       ./anatrj.x  frame_0001.pdb  PRC  918  110,112,120...  10.0\n");
                
		exit(1); 
	}
        
	return 0;
}
