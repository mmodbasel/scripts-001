/* 
 * File:   structures.h
 * Author: Martin Smiesko
 *
 * Created on June 19, 2014, 4:57 PM
 */

#ifndef STRUCTURES_H
#define	STRUCTURES_H

#include "definitions.h"
#include "points_vectors.h"

class   StructuralUnit;

class RawAtom {
protected:
        RawAtom *               previous;
	Point			coord;		// coordinates
	unsigned char		protonNum;	// proton number of the atom
	char			atomType[5];	// atom label = PDBX atom type
	double			radius;		// VdW radius in Angstroms
	double			charge;		// electrostatic charge
	int                     number;		// atom number as appears in the source PDB file
	char                    hybridization;	// hybridization state
        char                    dac;            // donor = 1, acceptor = 2, nothing = 0
	unsigned char           neighborCount;	// number of neighbors of the atom
        unsigned char           metalNeighborCount;	// number of neighbors of the atom
        StructuralUnit *        residue;        // residue pointer to which it belongs
        
        double                  occupancy;
        double                  Bfactor;
        // char                    chain;
        char                    altLocation;
        char                    species;        // protein, water, ligand, metal
public:
	void			setCoord(Point p) { coord = p; }
	void			setNumber(int n) { number = n; }
	void			setCharge(double q) { charge = q; }
        void                    setHybridization(char value) { hybridization = value; }
        void                    setDAC(char value) { dac = value; }
        void                    setResidue(StructuralUnit * r) { residue = r; }
        void                    setRadius(double r) { radius = r; }
        void                    setOccupancy(double value) { occupancy = value; }
        void                    setBfactor(double value) { Bfactor = value; }
        // void                    setChain(char value) { chain = value; }
        void                    setAltLocation(char value) { altLocation = value; }
        void                    setPrevious(RawAtom * value) { previous = value;  }
        void                    setSpecies(char value) { species = value; }

	int             	getNumber() { return number; }
	Point			getCoord() { return coord; }
	unsigned char		getProtonNum() { return protonNum; }
	double			getRadius() { return radius; }
	double			getCharge() { return charge; }
	char *			getAtomType() { return atomType; }
	char                    getHybridization() { return hybridization; }
        RawAtom *               getPrevious() { return previous; }
        char                    getDAC() { return dac; }
	unsigned char           getNeighborCount() { return neighborCount; }
        unsigned char           getMetalNeighborCount() { return metalNeighborCount; }
        StructuralUnit *        getResidue() { return residue; }
        double                  getOccupancy() { return occupancy; }
        double                  getBfactor() { return Bfactor; }
        // char                    getChain() { return chain; }
        char                    getAltLocation() { return altLocation; }
        char                    getSpecies() { return species; }
        
        void                    incNeighborCount() { neighborCount++; };

        // number  atom_type  alternate_location  x  y  z  occupancy  B-factor  proton_num  charge  index
        RawAtom(int n, char * at, char aL, Point p, double occ, double Bf, unsigned char pN, double Q) {
 		number = n;
		strcpy(atomType, at);
		coord = p;
		charge = Q;
                previous = NULL;

                // protonNum = giveProtonNumber(at);
                protonNum = pN;
                // radius = 0.0;
                radius = giveVdWRadius(pN);

                hybridization = UNCLEAR;
                dac = NO;
                neighborCount = 0;
                metalNeighborCount = 0;
                residue = NULL;
                
                altLocation = aL;
                occupancy = occ;
                Bfactor = Bf;
        };
};

class	Ring;

class	LigandAtom : public RawAtom {
	LigandAtom *	previous;	// pointer to the previous atom in the list
	LigandAtom *	next;		// pointer to the next atom in the list
	unsigned char	inRing;		// flag on, if atom already in a ring
	unsigned char	written;	// flag if the atom has been written in a ring searching algorithm
	LigandAtom *	neighbor[ATOM_MAX_NEIGHBORS];	// pointers to neighboring atoms
	char		bond[ATOM_MAX_NEIGHBORS];	// bond order, bonds in the same order like neighbors
	char		formalQ;	// formal charge
	Ring *		ring[ATOM_MAX_NEIGHBORS];	// pointers to rings, where this atom belong to
        unsigned char   nonPolarH;
        unsigned char   rotatableH;
        RawAtom *       metalNeighbor[METAL_MAX_NEIGHBORS];
        int             molNum;
        double          surfaceArea;
        double          polarArea;
        LigandAtom *    original;
        LigandAtom *    linked;         // atom linked to another residue
        char            neighT;         // count of neighbors as appeared in the template file
        
	public:
	// void		setPrevious(LigandAtom * a) { previous = a; }
	void		setNext(LigandAtom * a) { next = a; }
	void		setRingFlag(Ring * r) { ring[inRing] = r; inRing++; };
	void		setWritten() { written = YES; }
	void		clearWritten() { written = NO; }
	void		setNeighbor(LigandAtom * ptr) {
                                // security routine
                                int     i;
                                for(i = 0; i < neighborCount; i++)
                                        if(neighbor[i] == ptr) { 
                                                // printf("LIG:SET_NEIGH : already there!");
                                                return;
                                        }
                
				neighbor[neighborCount] = ptr;
                                bond[neighborCount] = UNCLEAR;
				neighborCount++;
				if(neighborCount > ATOM_MAX_NEIGHBORS) { printf("Error! Too many atoms connected to atom %3d\n",  neighborCount); exit(1); }
			}
	void		setNeighborAndBond(LigandAtom * ptr, char bondOrder) {
                                // security routine
                                int     i;
                                for(i = 0; i < neighborCount; i++) 
                                        if(neighbor[i] == ptr) { 
                                                // printf("LIG:SET_NEIGH : already there!");
                                                return;
                                        }
                
				neighbor[neighborCount] = ptr;
                                bond[neighborCount] = bondOrder;
				neighborCount++;
				if(neighborCount > ATOM_MAX_NEIGHBORS) { printf("Error! Too many atoms connected to atom %3d\n",  neighborCount); exit(1); }
			}
	void		setMetalNeighbor(RawAtom * ptr) {
                                // security routine
                                int     i;
                                for(i = 0; i < neighborCount; i++) 
                                        if(neighbor[i] == ptr) { 
                                                // printf("LIG(M):SET_NEIGH : already there!");
                                                return;
                                        }
                
				metalNeighbor[metalNeighborCount] = ptr;
				metalNeighborCount++;
				if(metalNeighborCount > METAL_MAX_NEIGHBORS) { printf("Error! Too many atoms connected to metal %3d\n",  metalNeighborCount); exit(1); }
			}        
	void		setBond(unsigned char n, char value) { bond[n] = value; }
	void		incFormalQ() { formalQ++; }
	void		decFormalQ() { formalQ--; }
        void            setVdW(double value) { radius = value; }
        void            setNonPolarH(unsigned char v) { nonPolarH = v; }
        void            setRotatableH(unsigned char v) { rotatableH = v; }
        void            setMolNum(int value) { molNum = value; }
        // void            setSurfaceArea(double value) { surfaceArea = value; }
        // void            setPolarArea(double value) { polarArea = value; }
        void            setOriginal(LigandAtom * value) { original = value; } 
        void            setLinked(LigandAtom * v) { linked = v; }
        // void            setAtomType(char * v) { strncpy(atomType, v, 4); atomType[4] = 0; } 
        void            setAtomType(char * v) { strcpy(atomType, v); }
        void            setNeighT(char v) { neighT = v; }
        
	unsigned char	getRingFlag() { return inRing; };	
	unsigned char	getWritten() { return written; }
	LigandAtom *	getPrevious() { return previous; }
	LigandAtom *	getNext() { return next; }
	LigandAtom *	getNeighbor(unsigned char i) { return neighbor[i]; }
        RawAtom *	getMetalNeighbor(unsigned char i) { return metalNeighbor[i]; }
	char		getBond(unsigned char n) { return bond[n]; }
	char		getFormalQ() { return formalQ; }
	Ring *		getNthRing(unsigned char n) { return ring[n]; }
        int             getMolNum() { return molNum; }
        double          getPolarArea() { return polarArea; }
        double          getSurfaceArea() { return surfaceArea; }
        LigandAtom *	getOriginal() { return original; }
        LigandAtom *    getLinked() { return linked; }
        char            getNeighT() { return neighT; }
        unsigned char	getRotatableH() { return rotatableH; };

        unsigned char	isNonPolarH() { return nonPolarH; };
        void            addSurfaceArea(double value) { surfaceArea += value; }
        void            addPolarArea(double value) { polarArea += value; }
        

        // LigandAtom(int n, char * at, Point p, double Q, LigandAtom * prev) : RawAtom(n, at, p, Q) {
        LigandAtom(int n, char * at, char aL, Point p, double occ, double Bf, unsigned char pN, double Q, LigandAtom * prev) : RawAtom(n, at, aL, p, occ, Bf, pN, Q) {
                for(neighborCount = 0; neighborCount < ATOM_MAX_NEIGHBORS; neighborCount++) {
                        neighbor[neighborCount] = NULL;
                        bond[neighborCount] = UNCLEAR;
                        ring[neighborCount] = NULL;
                }
                for(neighborCount = 0; neighborCount < METAL_MAX_NEIGHBORS; neighborCount++) {
                        metalNeighbor[neighborCount] = NULL;
                }
                previous = prev;
		if(prev != NULL) prev->setNext(this);
		next = NULL;
        	neighborCount = 0;
		inRing = NO;
		written = NO;
		formalQ = 0;
                nonPolarH = YES;
                this->setPrevious((RawAtom *) prev);
                
                if(isMetal(pN)) {
                        species = METAL;
                        // printf("METAL");
                }
                else species = LIGAND;
                dac = NO;
                molNum = -1;
                surfaceArea = polarArea = 0;
                original = NULL;
                linked = NULL;
                neighT = -1;
                rotatableH = NO;
	}
};

class	ProteinAtom : public RawAtom {
	ProteinAtom *	previous;	// pointer to the previous atom in the list
	ProteinAtom *	next;		// pointer to the next atom in the list
	// char		index;		// -3: H, O; -2: N, C; -1: Calpha; 0: Cbeta; 1: Gamma atoms; ...
        unsigned char   Hcount;
	ProteinAtom *	neighbor[ATOM_MAX_NEIGHBORS];	// pointers to neighboring atoms
        int             flag;           // multipurpose flag variable - clear before use :)   
	
	public:
	void		setNext(ProteinAtom * a) { next = a; }
	// void		setIndex(char i) { index = i; }
        void            setHcount(unsigned char value) { Hcount = value; }
        void            setVdW(double value) { radius = value; }
        void            setFlag(int value) { flag = value; }

	// char		getIndex() { return index; }
        unsigned char   getHcount() { return Hcount; }
	ProteinAtom *	getPrevious() { return previous; }
	ProteinAtom *	getNext() { return next; }
        ProteinAtom *	getNeighbor(unsigned char i) { return neighbor[i]; }
        int             getFlag() { return flag; }

	void		setNeighbor(ProteinAtom * ptr) {
                                // security routine
                                int     i;
                                for(i = 0; i < neighborCount; i++)
                                        if(neighbor[i] == ptr) {
                                                // printf("PROT:SET_NEIGH : already there!");
                                                return;
                                        }
                
				neighbor[neighborCount] = ptr;
				neighborCount++;
				if(neighborCount > ATOM_MAX_NEIGHBORS) { printf("Error! Too many atoms (%3d) connected to protein atom %s %d\n",  neighborCount, this->atomType, this->getNumber()); exit(1); }
			}
        // number  atom_type  alternate_location  x  y  z  occupancy  B-factor  proton_num  charge  index
        ProteinAtom(int n, char * at, char aL, Point p, double occ, double Bf, unsigned char pN, double Q, ProteinAtom * prev) : RawAtom(n, at, aL, p, occ, Bf, pN, Q) {
		previous = prev;
		if(prev != NULL) prev->setNext(this);
		next = NULL;
		// index = i;
                neighborCount = 0;
                Hcount = -1;
                this->setPrevious((RawAtom *) prev);
                // dac = atom_giveDonAcc(this);
                flag = 0;
	}
};

class	StructuralUnit {
protected:
	char			ID;			// one letter identifier of an amino acid or a residue, t for template, l for ligand
	int             	unitNumber;		// the unit number as read from the file
	unsigned char		activesite;		// active site flag
	ProteinAtom *           lastAtom;		// pointer to the last atom of the molecule
	ProteinAtom *		firstAtom;              // pointer to the first atom of the molecule
        char                    residueName[4];         // residue name long
        char                    chain;
public:
	void			setActiveSite() { activesite = YES; }
        void			clearActiveSite() { activesite = NO; }
	ProteinAtom *           getLastAtom() { return lastAtom; }
	ProteinAtom *           getFirstAtom() { return firstAtom; }
	int                     getUnitNumber() { return unitNumber; }
	char			getID() { return ID; }
        char *			getResidueName() { return residueName; }
	unsigned char		getActiveSite() { return activesite; }
	char                    getChain() { return chain; }
};

class	Residue : public StructuralUnit {
	Residue *		previous;	// pointer to the previous residue
	Residue *		next;		// pointer to the next residue
	Residue *		lastConf;	// pointer to the last conformation of this residue
	unsigned char		rotBonds;	// number of rotatable bonds
	int                     cloneN;		// number of clones of this residue
	int			bumpsToLigand;	// number of bumps a residue has with the ligand
	int			bumpsToProtein;	// number of bumps a residue has with all other residues including all atoms; -1 if not yet calculated
        Residue *               master;         // pointer to the master residue, NULL for the master itself
        int                     confID;         // conformation number
        unsigned char           bridged;        // flag indicating that the residue is H-bonded to a residue in the active site
        unsigned char           missAtom;       // missing atoms flag
	
	public:
	void			setPrevious(Residue * r) { previous = r; }
	void			setNext(Residue * r) { next = r; }
	void			setBumpsToLigand(int b) { bumpsToLigand = b; }
	void			setRotBonds(unsigned char n) { rotBonds = n; }
        void                    setMaster(Residue * r) { master = r; }
        void                    setConfID(int n) { confID = n; }
        void                    setBridged() { bridged = YES; }
        void                    setLastconf(Residue * r) { lastConf = r; }
        void                    setMissAtom() { missAtom++; }
        void                    setID(char c) { ID = c; }
        
	Residue *		getPrevious() { return previous; }
	Residue *		getNext() { return next; }
	Residue * 		getLastConf() { return lastConf; }
	unsigned char		getRotBonds() { return rotBonds; }
	// double			getTorsion(unsigned char n) { return torsion[n]; }
	int	 		getBumpsToLigand() { return bumpsToLigand; }
	unsigned int		getCloneN() { return cloneN; }
        Residue *               getMaster() { return master; }
        int                     getConfID() { return confID; }
        unsigned char           getBridged() { return bridged; }
        unsigned char           getMissAtom() { return missAtom; }

        void			decreaseCloneN(){ cloneN--; }
        void			addAtom(int num, char * at, char aL, Point p, double occ, double Bf, unsigned char protNum, double Q){
					// lastAtom = new ProteinAtom(num, at, aL, p, occ, Bf, protNum, Q, ind, lastAtom);
                                        lastAtom = new ProteinAtom(num, at, aL, p, occ, Bf, protNum, Q, lastAtom);
                                        lastAtom->setResidue(this);
                                        lastAtom->setSpecies(PROTEIN);
					// if it is the first atom
					if(lastAtom->getPrevious() == NULL) firstAtom = lastAtom;
                                }
     
        // constructor		
	Residue(int n, char l, char * rname, char ch, unsigned char rotBs, Residue * p) {
		unitNumber = n;
		ID = l;
		previous = p;
		if(p != NULL) p->setNext(this);
		next = NULL;
		rotBonds = rotBs;
		activesite = NO;
		lastConf = NULL;
		// atomCount = 0;
		cloneN = 0;
		lastAtom = NULL;
		// torsion[0] = 1000; torsion[1] = 1000; torsion[2] = 1000; torsion[3] = 1000;
		bumpsToProtein = -1;
		bumpsToLigand = 0;
                master = NULL;
                confID = 0;
                bridged = NO;
                strcpy(residueName, rname);
                residueName[3] = 0;
                chain = ch;
                missAtom = 0;
	}
};

class   ResidueList {
protected:
        ResidueList *           previous;
        Residue *		residue;
public:
	ResidueList *		getPrevious() { return previous; }
	Residue *		getResidue() { return residue; }

        ResidueList(ResidueList * p, Residue * r) { previous = p; residue = r; }
        ResidueList(Residue * r) { residue = r; }
};

class Water : public StructuralUnit {
	Water *		previous;
	Water *		next;
	unsigned char	structural;
        unsigned char	bridging;
        unsigned char   solvate;        // water which has 2 and more H-bonds to the same ligand
        unsigned char   status;         // YES = on, NO = off
        unsigned char   bond2W;         // bound to other water
        unsigned char   bond2L;         // bound to ligand
        unsigned char   bond2P;         // bound to protein
        unsigned char   bond2M;         // bound to metal
        unsigned char   flag;

	public:
                
        // number  atom_type  alternate_location  x  y  z  occupancy  B-factor  proton_num  charge  index
        void		addAtom(int num, char * at, char aL, Point p, double occ, double Bf, unsigned char protNum, double Q){                
	                        // lastAtom = new ProteinAtom(num, at, p, Q, i, lastAtom);
                                lastAtom = new ProteinAtom(num, at, aL, p, occ, Bf, protNum, Q, lastAtom);
                                lastAtom->setResidue(this);
                                lastAtom->setSpecies(WATER);
                                lastAtom->setDAC(DON_ACC);
                        }

	void		setStructural() { structural = YES; }
        void		setSolvate(unsigned char v) { solvate = v; }
        void		setBridging() { bridging = YES; }
	void		setNext(Water * w) { next = w; }
        void            setStatus(unsigned char s) { status = s; }
        void            setFlag(unsigned char v) { flag = v; }
        
        void            addBondToWater() { bond2W++; }
        void            addBondToLigand() { bond2L++; }
        void            addBondToProtein() { bond2P++; }
        void            addBondToMetal() { bond2M++; }
        
	unsigned char	getStructural() { return structural; }
        unsigned char	getSolvate() { return solvate; }
        unsigned char	getBridging() { return bridging; }
	Water * 	getPrevious() { return previous; }
	Water * 	getNext() { return next; }
        unsigned char   getStatus() { return status; }
        unsigned char   getBonds2Water() { return bond2W; }
        unsigned char   getBonds2Ligand() { return bond2L; }
        unsigned char   getBonds2Protein() { return bond2P; }
        unsigned char   getBonds2Metal() { return bond2M; }
        unsigned char   getFlag() { return flag; }
        ProteinAtom *   getO() {
                                     if(this->getLastAtom()->getProtonNum() == OXYGEN) return this->getLastAtom();
                                else if(this->getLastAtom()->getPrevious()->getProtonNum() == OXYGEN) return this->getLastAtom()->getPrevious();
                                else if(this->getLastAtom()->getPrevious()->getPrevious()->getProtonNum() == OXYGEN) return this->getLastAtom()->getPrevious()->getPrevious();
                                else { printf("Error getting water oxygen atom coords\n"); exit(1); }
                                return NULL;
                        }
	
        Water(int n, char l, char * rname, char ch, Water * p) {
		unitNumber = n;
                structural = NO;
		previous = p;
		lastAtom = NULL;
		next = NULL;
		if(p != NULL) p->setNext(this);
		ID = l;
                status = YES;
                strcpy(residueName, rname);
                residueName[3] = 0;
                chain = ch;
                bond2W = bond2L = bond2P = bond2M = 0;
                activesite = 0;
                bridging = NO;
                flag = 0;
                solvate = NO;
	}
};

class	Ring {
	unsigned char	size;
	unsigned char	lipophilic;	// set to 1 if all atoms are carbons or sulphur, and no substitutents
	unsigned char	linkerRing;	// set if the ring has more than 1 chain connected
	Point		centroid;
	Point		norm;
	LigandAtom *	member[MAX_RING_SIZE];
	Ring *		previous;
	Ring *		fusedto[MAX_RING_SIZE];
	unsigned char	aromatic;
	unsigned char	bondsdone;
	unsigned char	fusedCount;
	unsigned char	number;		// ring number as found by the ring finding routine
        unsigned char   dd;             // delocalized and double count
        unsigned char   ddouts;         // delocalized and double count
        
	public:
	void		setSize(unsigned char value) { size = value; }
	void		setCentroid(Point value) { centroid = value; }
	void		setNormal(Point v1) { norm = v1; }
	void		setMember(LigandAtom * value) { 
				if(size < MAX_RING_SIZE) {
                                        member[size] = value;
                                        size++;
                                }
                                else {
					printf("Error: Maximum ring size exceeded (n = %d).\n", MAX_RING_SIZE);
					exit(1);
				}
			}
	void		setFused(Ring * value) { 
				fusedto[fusedCount] = value;
				fusedCount++;
				if(fusedCount > MAX_RING_SIZE) {
					printf("Error: Maximum number of fused rings exceeded (n = %d).\n", MAX_RING_SIZE);
					exit(1);
				}
			}
	void		setLipophilic(unsigned char value) { lipophilic = value; };
	void		setLinkerRing(unsigned char value) { linkerRing = value; };
	void		setAromatic() { aromatic = YES; }
	void		setBondsDone() { bondsdone = YES; }
        void            setDelocAndDouble(unsigned char value) { dd = value; }
        void            setDelocAndDoubleOuts(unsigned char value) { ddouts = value; }
		
	unsigned char	getSize() { return size; }
	Point		getCentroid() { return centroid; }
	LigandAtom *	getMember(unsigned char value) { return member[value]; }
	unsigned char	getFusedCount() { return fusedCount; }
	Ring *		getFused(unsigned char value) { return fusedto[value]; }
	Ring *		getPrevious() { return previous; }
	unsigned char	getLipophilic() { return lipophilic; }
	unsigned char	getLinkerRing() { return linkerRing; }
	Point		getNorm() { return norm; }
	unsigned char	getAromatic() { return aromatic; }
	unsigned char	getBondsDone() { return bondsdone; }
	unsigned char	getNumber() { return number; }
        unsigned char	getDelocAndDouble() { return dd; }
        unsigned char	getDelocAndDoubleOuts() { return ddouts; }

	Ring(Ring * p) { 
		previous = p;
		size = 0;
		lipophilic = YES;
		linkerRing = NO;
		aromatic = NO;
		bondsdone = NO;
		fusedCount = 0;
                dd = ddouts = 0;
		// if this is the first ring, make 0, otherwise based on the previous one
		if(p == NULL) number = 0;
		else number = p->getNumber() + 1;
	}
};

template <class AtomType>
unsigned char	giveNeighborSum(AtomType * a, unsigned char protonNum) {
	unsigned char	i;
	unsigned char	sum = 0;
	
	for(i = 0; i < a->getNeighborCount(); i++) if(a->getNeighbor(i)->getProtonNum() == protonNum) sum++;
	return sum;
}

template <class AtomType>
unsigned char	isMeth(AtomType * a){
	// if the number of hydrogens is equal to number of attached atoms
	if(a->getProtonNum() == CARBON  &&  a->getNeighborCount() - 1 == giveNeighborSum(a, HYDROGEN)) return YES;
	else return NO;
}

template <class AtomType>
Point	giveCentroid(Ring * ring) {
	Point		center = { 0, 0, 0 };
	Point		p;
        AtomType *      a;
	unsigned char	size = ring->getSize();
	unsigned char	i, j, k;
	unsigned char	connections = 0;
	unsigned char	isInThisRing;
	
	// count all connections
	// for each atom of the ring
	for(i = 0; i < size; i++) {
		// analyze all its neighbor atoms
		a = ring->getMember(i);
		// if the belong this current ring
		for(j = 0; j < a->getNeighborCount(); j++){
			isInThisRing = NO;
                        // compare with all memebers of this ring
			for(k = 0; k < size; k++) if(a->getNeighbor(j) == ring->getMember(k)) isInThisRing = YES;

			// if it is not in the ring 
			if(isInThisRing == NO){
				k = a->getNeighbor(j)->getProtonNum();
                                // and it is not a halogend or hydrogen or methyl, then it must be a connection to some atom
				if(!(k == HYDROGEN  ||  k == FLUORINE  ||  k == CHLORINE  ||  k == BROMINE  ||  k == IODINE  ||  isMeth(a->getNeighbor(j)) == YES)) connections++;
			}
		}
	}
	
	for(i = 0; i < size; i++) {
		p = ring->getMember(i)->getCoord();
		center = center + p;
		// if not carbon or sulphur, clear the lipophilic flag
		if(!(ring->getMember(i)->getProtonNum() == CARBON || ring->getMember(i)->getProtonNum() == SULPHUR)) ring->setLipophilic(NO);
	}
	
	if(connections > 1) { ring->setLipophilic(NO); ring->setLinkerRing(YES); }
	
	center.x = center.x / size;
	center.y = center.y / size;
	center.z = center.z / size;
	return center;
}

Point	giveRingNormal(Ring * ring) {
	Point		center = { 0, 0, 0 };
	Point		p;
	unsigned char	size = ring->getSize();
	unsigned char	i, j;
	int             n = 0;
	Vector		v1, v2;
	Vector		V;
	
	p.x = 0; p.y = 0; p.z = 0;
	for(i = 0; i < size; i++) {
		p = ring->getMember(i)->getCoord();
		center = center + p;
	}
	center.x = center.x / size;
	center.y = center.y / size;
	center.z = center.z / size;
	
	
	p.x = 0; p.y = 0; p.z = 0;
	for(i = 0; i < size; i++) {
		for(j = i; j < size; j++) {
			if(i != j) {
				v1 = giveVector(center, ring->getMember(i)->getCoord());
				v2 = giveVector(center, ring->getMember(j)->getCoord());
				V = giveVectorProduct(v1, v2);
				p.x = p.x + V.x;
				p.y = p.y + V.y;
				p.z = p.z + V.z;
				n++;
			}
		}
	}
	p.x = p.x / (double)n;
	p.y = p.y / (double)n;
	p.z = p.z / (double)n;


	V.x = p.x; V.y = p.y; V.z = p.z;
	V = giveUnitVector(V);
	p.x = V.x; p.y = V.y; p.z = V.z;
	p = p + center;
	
	return p;
}

class	Centroid {
	private:
	Point		cPoint;
	Centroid *	previous;
	Point		ringNormal;
	Ring *		ring;
	
	public:
	void	setCPoint(Point pt) { cPoint = pt; }
	void	setNormal(Point pt) { ringNormal = pt; }
	void	setPrevious(Centroid * value) { previous = value; }
	
	Point		getCPoint() { return cPoint; }
	Centroid *	getPrevious() { return previous; }
	Point		getNormal() { return ringNormal; }
	Ring *		getRing() { return ring; }
	
	Centroid(Point p, Point rn, Ring * r, Centroid * c) { cPoint = p; ringNormal = rn; ring = r; previous = c; }
};

class	Molecule : public StructuralUnit  {
	private:
	LigandAtom *	lastAtom;
        // LigandAtom *	firstAtom;
	LigandAtom *	help;
	// char		name[5];
	Centroid *	centroid;
	Molecule *	previous;
	Ring *		ring;
	char		totalQ;
        unsigned char   flexible;
        int             rotBonds;
        double          selfIndex;
        double          MW;
        int             n_carbohyd_atoms;
        int             n_atoms;
        int             n_atoms_in_template;
        
        unsigned char   proteinBound;
        unsigned char   metalloligand;
        unsigned char   property;       // heme, modified aa, etc.
        unsigned char   carbohydrate;
        unsigned char   useful;
        unsigned char   wBridge;
        unsigned char   wClose;
        unsigned char   wBound;
        unsigned char   missAtom;
        
        char		dbrefStr[MAX_DBREF_LENGTH];
        char            fullResname[MAX_RES_NAME_LENGTH];
        char            molRootName[128];
        
	public:
	// void	addAtom(int num, char * name, Point p, double Q) {
        void            addAtom(int num, char * at, char aL, Point p, double occ, double Bf, unsigned char protNum, double Q){
                                // lastAtom = new ProteinAtom(num, at, aL, p, occ, Bf, protNum, Q, ind, lastAtom);
                                help = lastAtom;
                                lastAtom = new LigandAtom(num, at, aL, p, occ, Bf, protNum, Q, lastAtom);
                                lastAtom->setResidue(this);

                                if(help == NULL) firstAtom = (ProteinAtom *) lastAtom;
                                else help->setNext(lastAtom);
                                
                                if(protNum == OXYGEN  ||  protNum == OXYGEN  ||  protNum == SULPHUR) lastAtom->setDAC(DON_ACC_UNSURE);
                        }
	void	removeAllAtoms() {
			while(lastAtom){
				help = lastAtom;
				lastAtom = lastAtom->getPrevious();
				delete help;
			}
		}
	void	setCentroid(Centroid * value) { centroid = value; }
	void	setRing(Ring * value) { ring = value; }
	void	setPrevious(Molecule * value) { previous = value; }
	void	setTotalQ(char value) { totalQ = value; }
	void    setSelfIndex(double value) { selfIndex = value; }
        void    setRotBonds(int value) { rotBonds = value; }
        void    setProperty(unsigned char value) { property = value; }
        void    setProteinBound(unsigned char value) { proteinBound = value; }
        void    setMetalloligand(unsigned char value) { metalloligand = value; }
        void    setCarbohydrate() { carbohydrate = YES; }
        void    setUseful(unsigned char value) { useful = value; }
        void    setCloseWatersCount(unsigned char value) { wClose = value; }
        void    setBridgeWatersCount(unsigned char value) { wBridge = value; }
        void    setBoundWatersCount(unsigned char value) { wBound = value; }
        void    setAtomNumber(int value) { n_atoms = value; }
        void    setAtomNumberInTemplate(int value) { n_atoms_in_template = value; }
        void    setCBHAtomNumber(int value) { n_carbohyd_atoms = value; }
        void    setMissAtom(unsigned char value) { missAtom = value; }
        void    setDBrefStr(char * v) { strcpy(dbrefStr, v); }
        void    setFullResname(char * v) { strcpy(fullResname, v); }
        void    setMolRootName(char * v) { strcpy(molRootName, v); }
        
        void    addResname(char * v) { 
                        strcat(fullResname, v);
                        if(strlen(fullResname) >= sizeof(fullResname)) {
                                printf("Maximum string length (%d) for the full residue name exceeded\n", MAX_RES_NAME_LENGTH);
                                exit(1);
                        }
                }
        
        void	incTotalQ() { totalQ++; }
        void	decTotalQ() { totalQ--; }
        void    addFlexible() { flexible++; }
        void    addWeight(double w) { MW = MW + w; }

	LigandAtom *	getLastAtom() { return lastAtom; }
	Centroid *	getCentroid() { return centroid; }
	Ring *		getRing() { return ring; }
	Molecule *	getPrevious() { return previous; }
	char		getTotalQ() { return totalQ; }
        LigandAtom *	getFirstAtom() { return (LigandAtom *) firstAtom; }
        unsigned char   getFlexible() { return flexible; }
        double          getSelfIndex() { return selfIndex; }
        int             getRotBonds() { return rotBonds; }
        unsigned char   getProperty() { return property; }
        unsigned char   getProteinBound() { return proteinBound; }
        unsigned char   getMetalloligand() { return metalloligand; }
        unsigned char   getCarbohydrate() { return carbohydrate; }
        unsigned char   getUseful() { return useful; }
        unsigned char   getCloseWatersCount() { return wClose; }
        unsigned char   getBridgeWatersCount() { return wBridge; }
        unsigned char   getBoundWatersCount() { return wBound; }
        int             getAtomNumber() { return n_atoms; }
        int             getAtomNumberInTemplate() { return n_atoms_in_template; }
        int             getCBHAtomNumber() { return n_carbohyd_atoms; }
        double          getMolecularWeight() { return MW; }
        char *		getDBrefStr() { return dbrefStr; }
        unsigned char   getMissAtom() { return missAtom; }
        char *		getFullResname() { return fullResname; }
        char *		getMolRootName() { return molRootName; }
        
        Molecule(int n, char l, char * rname, char ch, Molecule * p) {
		unitNumber = n;
		lastAtom = NULL;
		strcpy(residueName, rname);
                residueName[3] = 0;
		previous = p;
                ID = l;
		totalQ = 0;
                flexible = 0;
                selfIndex = rotBonds = 0;
                property = 0;
                proteinBound = NO;
                chain = ch;
                metalloligand = NO;
                carbohydrate = NO;
                useful = YES;
                wBridge = wClose = wBound = 0;
                MW = n_carbohyd_atoms = n_atoms = 0;
                dbrefStr[0] = 0;
                missAtom = 0;
                n_atoms_in_template = 0;
                fullResname[0] = 0;
                strcpy(molRootName, "noname.dat");
	}
};

#endif	/* STRUCTURES_H */

