/* 
 * File:   myMMGBSA.cpp
 * Author: Martin Smiesko
 *
 * Created on January 22, 2019, 10:27 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define YES             1
#define NO              0

#define EHPHOB              -0.3
#define CUTOFF_NonBonded    12.0
#define MIN_EB_FOR_REPORTING    -0.1

#define	MAXATOMS	65536
#define	MAXNEIGHBORS	6
#define RESNAMELENGTH   5
#define RESNUMOFFSET    88
#define CHAINOFFSET     94
#define CHARGEOFFSET    101
#define RESNAMEOFFSET   119

#define SINGLE_BOND     1
#define DOUBLE_BOND     2
#define TRIPLE_BOND     3

#define	PI 			3.14159265358979
#define RAD_180                 3.14159265358979

#define RAD_90                  PI / 2.0
#define RAD_60                  PI / 3.0
#define RAD_45                  PI / 4.0
#define RAD_30                  PI / 6.0
#define RAD_20          	PI / 9.0
#define RAD_15                  PI / 12.0
#define RAD_10                  PI / 18.0
#define RAD_9                   PI / 20.0
#define RAD_6                   PI / 30.0
#define RAD_5                   PI / 36.0
#define RAD_3                   PI / 60.0
#define RAD_2                   PI / 90.0
#define RAD_1                   PI / 180.0

#define RAD_360                 2.0 * PI
#define RAD_150                 5.0 * RAD_30
#define RAD_125                 25.0 * RAD_5
#define RAD_120                 2.0 * RAD_60
#define RAD_105                 7.0 * RAD_15
#define RAD_100                 10.0 * RAD_10
#define RAD_80          	4.0 * RAD_20
#define RAD_75                  5.0 * RAD_15
#define RAD_50                  5.0 * RAD_10
#define RAD_40                  2.0 * RAD_20
#define RAD_35                  7.0 * RAD_5
#define RAD_25                  5.0 * RAD_5
#define RAD_8                   4.0 * RAD_2
#define RAD_4                   2.0 * RAD_2

#define HB_THRESHOLD_ANGLE  RAD_105


#define 	ACTINIUM	89
#define 	ALUMINIUM	13
#define 	AMERICIUM	95
#define 	ANTIMONY	51
#define 	ARGON           18
#define 	ARSENIC         33
#define 	ASTATINE	85
#define 	BARIUM          56
#define 	BERKELIUM	97
#define 	BERYLLIUM	4
#define 	BISMUTH         83
#define 	BOHRIUM         107
#define 	BORON           5
#define 	BROMINE         35
#define 	CADMIUM         48
#define 	CAESIUM         55
#define 	CALCIUM         20
#define 	CALIFORNIUM	98
#define 	CARBON          6
#define 	CERIUM          58
#define 	CHLORINE	17
#define 	CHROMIUM	24
#define 	COBALT          27
#define 	COPPER          29
#define 	CURIUM          96
#define 	DARMSTADTIUM	110
#define 	DUBNIUM         105
#define 	DYSPROSIUM	66
#define 	EINSTEINIUM	99
#define 	ERBIUM          68
#define 	EUROPIUM	63
#define 	FERMIUM         100
#define 	FLUORINE	9
#define 	FRANCIUM	87
#define 	GADOLINIUM	64
#define 	GALLIUM         31
#define 	GERMANIUM	32
#define 	GOLD            79
#define 	HAFNIUM         72
#define 	HASSIUM         108
#define 	HELIUM          2
#define 	HOLMIUM         67
#define 	HYDROGEN	1
#define 	INDIUM          49
#define 	IODINE          53
#define 	IRIDIUM         77
#define 	IRON            26
#define 	KRYPTON         36
#define 	LANTHANUM	57
#define 	LAWRENCIUM	103
#define 	LEAD            82
#define 	LITHIUM         3
#define 	LUTETIUM	71
#define 	MAGNESIUM	12
#define 	MANGANESE	25
#define 	MEITNERIUM	109
#define 	MENDELEVIUM	101
#define 	MERCURY         80
#define 	MOLYBDENUM	42
#define 	NEODYMIUM	60
#define 	NEON            10
#define 	NEPTUNIUM	93
#define 	NICKEL          28
#define 	NIOBIUM         41
#define 	NITROGEN	7
#define 	NOBELIUM	102
#define 	OSMIUM          76
#define 	OXYGEN          8
#define 	PALLADIUM	46
#define 	PHOSPHORUS	15
#define 	PLATINUM	78
#define 	PLUTONIUM	94
#define 	POLONIUM	84
#define 	POTASSIUM	19
#define 	PRASEODYMIUM	59
#define 	PROMETHIUM	61
#define 	PROTACTINIUM	91
#define 	RADIUM          88
#define 	RADON           86
#define 	RHENIUM         75
#define 	RHODIUM         45
#define 	RUBIDIUM	37
#define 	RUTHENIUM	44
#define 	RUTHERFORDIUM	104
#define 	SAMARIUM	62
#define 	SCANDIUM	21
#define 	SEABORGIUM	106
#define 	SELENIUM	34
#define 	SILICON         14
#define 	SILVER          47
#define 	SODIUM          11
#define 	STRONTIUM	38
#define 	SULPHUR         16
#define 	TANTALUM	73
#define 	TECHNETIUM	43
#define 	TELLURIUM	52
#define 	TERBIUM         65
#define 	THALLIUM	81
#define 	THORIUM         90
#define 	THULIUM         69
#define 	TIN             50
#define 	TITANIUM	22
#define 	TUNGSTEN	74
#define 	URANIUM         92
#define 	VANADIUM	23
#define 	XENON           54
#define 	YTTERBIUM	70
#define 	YTTRIUM         39
#define 	ZINC            30
#define 	ZIRCONIUM	40
#define         DEUTERIUM       1


double	giveVdWRadius(unsigned char N) {
        switch(N) {
            case	ACTINIUM	: return 2.00;	break;
            case	ALUMINIUM	: return 2.00;	break;
            case	AMERICIUM	: return 2.00;	break;
            case	ANTIMONY	: return 2.00;	break;
            case	ARGON		: return 1.88;	break;
            case	ARSENIC		: return 1.85;	break;
            case	ASTATINE	: return 2.00;	break;
            case	BARIUM		: return 2.00;	break;
            case	BERKELIUM	: return 2.00;	break;
            case	BERYLLIUM	: return 2.00;	break;
            case	BISMUTH		: return 2.00;	break;
            case	BOHRIUM		: return 2.00;	break;
            case	BORON		: return 2.00;	break;
            case	BROMINE		: return 1.85;	break;
            case	CADMIUM		: return 1.58;	break;
            case	CAESIUM		: return 2.00;	break;
            case	CALCIUM		: return 2.00;	break;
            case	CALIFORNIUM	: return 2.00;	break;
            case	CARBON		: return 1.70;	break;
            case	CERIUM		: return 2.00;	break;
            case	CHLORINE	: return 1.75;	break;
            case	CHROMIUM	: return 2.00;	break;
            case	COBALT		: return 2.00;	break;
            case	COPPER		: return 1.40;	break;
            case	CURIUM		: return 2.00;	break;
            case	DARMSTADTIUM	: return 2.00;	break;
            case	DUBNIUM		: return 2.00;	break;
            case	DYSPROSIUM	: return 2.00;	break;
            case	EINSTEINIUM	: return 2.00;	break;
            case	ERBIUM		: return 2.00;	break;
            case	EUROPIUM	: return 2.00;	break;
            case	FERMIUM		: return 2.00;	break;
            case	FLUORINE	: return 1.47;	break;
            case	FRANCIUM	: return 2.00;	break;
            case	GADOLINIUM	: return 2.00;	break;
            case	GALLIUM		: return 1.87;	break;
            case	GERMANIUM	: return 2.00;	break;
            case	GOLD		: return 1.66;	break;
            case	HAFNIUM		: return 2.00;	break;
            case	HASSIUM		: return 2.00;	break;
            case	HELIUM		: return 1.40;	break;
            case	HOLMIUM		: return 2.00;	break;
            case	HYDROGEN	: return 1.09;	break;
            case	INDIUM		: return 1.93;	break;
            case	IODINE		: return 1.98;	break;
            case	IRIDIUM		: return 2.00;	break;
            case	IRON		: return 2.00;	break;
            case	KRYPTON		: return 2.02;	break;
            case	LANTHANUM	: return 2.00;	break;
            case	LAWRENCIUM	: return 2.00;	break;
            case	LEAD		: return 2.02;	break;
            case	LITHIUM		: return 1.82;	break;
            case	LUTETIUM	: return 2.00;	break;
            case	MAGNESIUM	: return 1.73;	break;
            case	MANGANESE	: return 2.00;	break;
            case	MEITNERIUM	: return 2.00;	break;
            case	MENDELEVIUM	: return 2.00;	break;
            case	MERCURY		: return 1.55;	break;
            case	MOLYBDENUM	: return 2.00;	break;
            case	NEODYMIUM	: return 2.00;	break;
            case	NEON		: return 1.54;	break;
            case	NEPTUNIUM	: return 2.00;	break;
            case	NICKEL		: return 1.63;	break;
            case	NIOBIUM		: return 2.00;	break;
            case	NITROGEN	: return 1.55;	break;
            case	NOBELIUM	: return 2.00;	break;
            case	OSMIUM		: return 2.00;	break;
            case	OXYGEN		: return 1.52;	break;
            case	PALLADIUM	: return 1.63;	break;
            case	PHOSPHORUS	: return 1.80;	break;
            case	PLATINUM	: return 1.72;	break;
            case	PLUTONIUM	: return 2.00;	break;
            case	POLONIUM	: return 2.00;	break;
            case	POTASSIUM	: return 2.75;	break;
            case	PRASEODYMIUM	: return 2.00;	break;
            case	PROMETHIUM	: return 2.00;	break;
            case	PROTACTINIUM	: return 2.00;	break;
            case	RADIUM		: return 2.00;	break;
            case	RADON		: return 2.00;	break;
            case	RHENIUM		: return 2.00;	break;
            case	RHODIUM		: return 2.00;	break;
            case	RUBIDIUM	: return 2.00;	break;
            case	RUTHENIUM	: return 2.00;	break;
            case	RUTHERFORDIUM	: return 2.00;	break;
            case	SAMARIUM	: return 2.00;	break;
            case	SCANDIUM	: return 2.00;	break;
            case	SEABORGIUM	: return 2.00;	break;
            case	SELENIUM	: return 1.90;	break;
            case	SILICON		: return 2.10;	break;
            case	SILVER		: return 1.72;	break;
            case	SODIUM		: return 2.27;	break;
            case	STRONTIUM	: return 2.00;	break;
            case	SULPHUR		: return 1.80;	break;
            case	TANTALUM	: return 2.00;	break;
            case	TECHNETIUM	: return 2.00;	break;
            case	TELLURIUM	: return 2.06;	break;
            case	TERBIUM		: return 2.00;	break;
            case	THALLIUM	: return 1.96;	break;
            case	THORIUM		: return 2.00;	break;
            case	THULIUM		: return 2.00;	break;
            case	TIN		: return 2.17;	break;
            case	TITANIUM	: return 2.00;	break;
            case	TUNGSTEN	: return 2.00;	break;
            case	URANIUM		: return 1.86;	break;
            case	VANADIUM	: return 2.00;	break;
            case	XENON		: return 2.16;	break;
            case	YTTERBIUM	: return 2.00;	break;
            case	YTTRIUM		: return 2.00;	break;
            case	ZINC		: return 1.39;	break;
            case	ZIRCONIUM	: return 2.00;	break;
            default: break;
        }
        
        printf("VdW: Unsupported element proton number %d\n", N);
        exit(1);
        
	return 0;
}

class Point {
	public:
	double	x, y, z;
	
	Point	operator+(Point p) { 
		Point p0;
		p0.x = x + p.x;
		p0.y = y + p.y;
		p0.z = z + p.z;
		return p0;
	}
	
	Point	operator-(Point p) { 
		Point p0;
		p0.x = x - p.x;
		p0.y = y - p.y;
		p0.z = z - p.z;
		return p0;
	}
};

typedef struct {
    int     neighbor[MAXNEIGHBORS];
    int     bondorder[MAXNEIGHBORS];
    
} NeighRec;

class Atom {
    Point   coord;
    int     atype;
    double  charge;
    int     number;
    int     neig[MAXNEIGHBORS];
    int     bond[MAXNEIGHBORS];
    char    resname[RESNAMELENGTH];
    int     resnum;
    char    chain;
    Atom *  previous;

    int     neighborcount;
    char    moltype;
    int     element;
    double  vdw;

public:
    Point   getCoords()  { return coord; }
    int     getAtomType(){ return atype; }
    double  getCharge()  { return charge; }
    int     getNumber()  { return number; }
    char *  getResName() { return resname; }
    int     getResNum()  { return resnum; }
    char    getChain()   { return chain; }
    Atom *  getPrevious(){ return previous; }
    int     getNeig(int n) { return neig[n]; }
    int     getBond(int n) { return bond[n]; }
    char    getMolType() { return moltype; }
    int     getNeighCount() { return neighborcount; }
    int     getElement() { return element; }
    double  getRadius() { return vdw; }
    
    void    setMolType(char mt) { moltype = mt; }    
    void    setNeighbors(NeighRec nr) {
                int i, j = 0;
                for(i = 0; i < MAXNEIGHBORS; i++) {
                    neig[i] = nr.neighbor[i];
                    bond[i] = nr.bondorder[i];
                    if(neig[i] != 0) j++;
                }
                neighborcount = j;
            }

    Atom(Point Acoord, int Atype, double Acharge, int Anumber, Atom * Aprevious) {
        coord = Acoord;
        atype = Atype;
        charge = Acharge;
        number = Anumber;
        previous = Aprevious;
    }

    Atom(Point Acoord, int Atype, double Acharge, int Anumber, char * Resname, int Resnum, char Chain, Atom * Aprevious) {
        coord = Acoord;
        atype = Atype;
        charge = Acharge;
        number = Anumber;
        strncpy(resname, Resname, 4); resname[4] = 0;
        resnum = Resnum;
        chain = Chain;
        previous = Aprevious;
        
             if(atype <= 14) element = CARBON;
        else if(atype <= 23) element = OXYGEN;
        else if(atype <= 40) element = NITROGEN;
        else if(atype <= 48) element = HYDROGEN;
        else if(atype <= 52) element = SULPHUR;
        else if(atype == 53) element = PHOSPHORUS;
        else if(atype == 56) element = FLUORINE;
        else if(atype == 57) element = CHLORINE;
        else if(atype == 58) element = BROMINE;
        else if(atype == 59) element = IODINE;
        else if(atype == 60) element = SILICON;
        else if(atype == 66) element = SODIUM;
        else if(atype == 67) element = POTASSIUM;
        else if(atype == 70) element = CALCIUM;
        else if(atype == 72) element = MAGNESIUM;
        else if(atype == 157) {
            printf("WARNING: Atomtype %d assumed to be quarternary N, please check!\n", atype);
            element = NITROGEN;
        }
        else if(atype == 79) element = IRON;
        else if(atype == 80) element = IRON;
        else if(atype == 85) element = COPPER;
        else if(atype == 86) element = COPPER;
        else if(atype == 87) element = ZINC;
        // ions
        else if(atype == 151) element = CHLORINE;
        
        // else if(atype == ) element = ;
        else { printf("Unsupported atom type/element %d\n", atype); exit(1); };
        
        vdw = giveVdWRadius(element); 
    }
};

double  giveDistance(Point p1, Point p2) {
    p1.x = p1.x - p2.x;
    p1.y = p1.y - p2.y;
    p1.z = p1.z - p2.z;
    return (sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z)); 
}

int isLipoAtom(Atom * a) {
    if(a->getElement() == CARBON  ||  a->getElement() == CHLORINE  ||  a->getElement() == BROMINE  ||  a->getElement() == FLUORINE  ||  a->getElement() == IODINE) return YES;
    if(a->getElement() == SULPHUR  &&  a->getNeighCount() <= 2) return YES;
    return NO;
}

class Vector {
	public:
	double	x, y, z;
	
	Vector	operator+ (Vector v) { 
		Vector v0;
		v0.x = x + v.x;
		v0.y = y + v.y;
		v0.z = z + v.z;
		return v0;
	}
	
	Vector	operator- (Vector v) { 
		Vector v0;
		v0.x = x - v.x;
		v0.y = y - v.y;
		v0.z = z - v.z;
		return v0;
	}
};

Vector	point2vector (Point p) {
	Vector	v;
	v.x = p.x;
	v.y = p.y;
	v.z = p.z;
	return v;
}

double	giveVectorLengthSquared(Vector v1) { return (v1.x*v1.x + v1.y*v1.y + v1.z*v1.z); }

double giveAngle(Point A, Point B, Point C) {
	double	asq, bsq, csq, arg;
	asq = giveVectorLengthSquared(point2vector(B - C));
        bsq = giveVectorLengthSquared(point2vector(A - C));
	csq = giveVectorLengthSquared(point2vector(B - A));
        arg = (bsq - asq - csq) / (-2.0 * sqrt(asq * csq));
        if(arg < -1.0) arg = -1.0;
        else if(arg > 1.0) arg = 1.0;
        return acos(arg);
}

double  give_e_hphob(Atom * a1, Atom * a2, double d) {
    double  scale, e = 0;

    scale = 2.0 * (d - a1->getRadius() - a2->getRadius() - 2.0) / 3.0;
    
          if(scale <= -1.0) e = EHPHOB;
    else if (scale <   1.0) e = EHPHOB * (0.25 * scale * scale * scale  -  0.75 * scale  +  0.5);
    else                    e = 0;
    
    return e;
}

double  give_e_elstat(Atom * a1, Atom * a2, double d) {
    // normal
    // return (332.4 * a1->getCharge() * a2->getCharge() / d);
    
    if(d > CUTOFF_NonBonded) return 0;
    
    // attenuated   ( 1 - (R/CTOFNB)**2 )**2
    double f;
    f = d / CUTOFF_NonBonded;
    f = f * f;
    f = 1 - f;
    f = f * f;
    return (f * 332.4 * a1->getCharge() * a2->getCharge() / d);
}

Atom *  index2atom(Atom * lastatom, int index) {
    Atom * a = lastatom;
    while(a) {
        if(a->getNumber() == index) return a;
        a = a->getPrevious();
    }
    return NULL;
}

// according to Yeti force-field - linearity (X-H...Y) corrected H-bond energy function
double	giveHbondEnergy(Point X, Point H, Point Y, unsigned char Nx, unsigned char Ny, int lowCutOff = NO){
	double	C, D;
	double	cosine;
	double	radius;
	double	wellDepth;
	double	rOpt;
	double	angleXHY;
        Point   p;
	
	p = H - Y;
        radius = p.x * p.x  +  p.y * p.y  +  p.z * p.z;
        
        
	// Arbitrary 3 A limit for H-bonds to O & N, 4 A to S
	if((Ny < SULPHUR  &&  radius < 9.0)  ||  radius < 16.0){
		angleXHY = giveAngle(X, H, Y);
		if(angleXHY > HB_THRESHOLD_ANGLE) {
                        //make sqrt out of dist
                        radius = sqrt(radius);
                        
			     if(Nx == 8  &&  Ny == 8)  { wellDepth = -4.946; rOpt = 1.746; }
			else if(Nx == 8  &&  Ny == 7)  { wellDepth = -4.655; rOpt = 1.878; }
			else if(Nx == 8  &&  Ny == 16) { wellDepth = -1.746; rOpt = 2.535; }
			else if(Nx == 7  &&  Ny == 8)  { wellDepth = -4.073; rOpt = 1.877; }
			else if(Nx == 7  &&  Ny == 7)  { wellDepth = -3.491; rOpt = 2.003; }
			else if(Nx == 7  &&  Ny == 16) { wellDepth = -1.455; rOpt = 2.667; }
			else if(Nx == 16 &&  Ny == 8)  { wellDepth = -2.328; rOpt = 2.099; }
			else if(Nx == 16 &&  Ny == 7)  { wellDepth = -2.037; rOpt = 2.088; }
			else if(Nx == 16 &&  Ny == 16) { wellDepth = -1.164; rOpt = 3.009; }
			else { printf("Error: Wrong atom type(s) (Nx = %d, Ny = %d)for hydrogen bonding passed.\n", Nx, Ny); exit(1); }

                        // correction for too short distances sometimes seen in MD frames
                        if(lowCutOff == YES  &&  radius < rOpt) radius = rOpt;

			C = -5 * wellDepth * pow(rOpt, 12);
			D = -6 * wellDepth * pow(rOpt, 10);

			cosine = cos(angleXHY);
			cosine = cosine * cosine;
			
			return ( cosine * (C/pow(radius, 12) - D/pow(radius, 10)) );
		}
	}
	return 0;
}

int isHbondDonor(Atom * a) {}

int isHbondAcceptor(Atom * a) {
    int i;
    int flag;
         if(a->getElement() == OXYGEN) return YES;
    else if(a->getElement() == NITROGEN) {
        flag = 0;
        for(i = 0; i < a->getNeighCount(); i++) if(a->getBond(i) == DOUBLE_BOND) flag++;
            
             if(a->getNeighCount() == 1) return YES;
        else if(a->getNeighCount() == 2  &&  flag == 1) return YES;
        else if(a->getNeighCount() == 3  &&  flag == 0  &&  a->getAtomType() != 25) return YES;  // macromodel atom type is 25 is N-SP2
    }
    else if (a->getElement() == SULPHUR) {
        if(a->getNeighCount() <= 2) return YES;
    }
    return NO;
}

int main(int argc, char *argv[])
{
	FILE        *inDAT_file = NULL;
	char        lineBuffer[256];
	Atom *      a1 = NULL;
        Atom *      a2 = NULL;
        Atom *      aLast = NULL;
        
	Point       c;
	NeighRec    nrec;
        double      q;
        int         at;
        char        rnam[RESNAMELENGTH];
        int         rnum;
        char        chain;
        char        string[16];
        
	int         nOfAtoms;
	int         i, j;
        int         ligofi; // ligand of interest
        
        double      tE_hphob = 0;
        double      e_hphob = 0;
        double      tE_hphob_P = 0;
        double      tE_hphob_L = 0;
        double      tE_hphob_M = 0;
        double      tE_hphob_H = 0;
        double      tE_hphob_W = 0;

        double      tE_elstat = 0;
        double      e_elstat = 0;
        double      tE_elstat_P = 0;
        double      tE_elstat_L = 0;
        double      tE_elstat_M = 0;
        double      tE_elstat_H = 0;
        double      tE_elstat_W = 0;

        double      tE_Hbond = 0;
        double      e_Hbond = 0;
        double      tE_Hbond_P = 0;
        double      tE_Hbond_L = 0;
        double      tE_Hbond_M = 0;
        double      tE_Hbond_H = 0;
        double      tE_Hbond_W = 0;

        double      distance;
        Atom *      Xatom;
        Atom *      Yatom;
        
        int         CorrectShortHB = 0;
        
	if (argc == 5)
	{
		//  O P E N   I N P U T   F I L E   -  MACROMODEL .dat file
		printf("Opening source .mae file %s for reading... ", argv[1]);
		inDAT_file = fopen(argv[1], "rt");
		if (inDAT_file == NULL)  { printf("Error! Problem opening file.\n"); exit(1); }
		else printf("OK.\n");
                
                if(sscanf(argv[3], "%d", &ligofi) != 1) { printf("Error reading ligand number.\n"); exit(1); }
                if(sscanf(argv[4], "%d", &CorrectShortHB) != 1) { printf("Error reading CorrectShortHB. Use 0 to deactivate or 1 to activate.\n"); exit(1); }

		// read source dat file in order to decode matrix
		if(fgets(lineBuffer, sizeof(lineBuffer), inDAT_file) != NULL) {
			if(sscanf(lineBuffer, "%d", &nOfAtoms) == 1) {
				if(nOfAtoms > MAXATOMS) { printf("Error! Molecule is too big (%d atoms, max=%d)\n", nOfAtoms, MAXATOMS); exit(1); }
			}
			else { printf("Problem reading .DAT file header.\n"); exit(1); }
		}

/*
DAT FORMAT LINES:
  1    2   3     4 5     6 7     8 9    10 11   12 13     14          15          16       
  44   119 1     0 0     0 0     0 0     0 0     0 0   33.757072    9.479884   18.040653    28RA  21  0.46000  0.46000 ARG  HH21     
  16   131 1   142 1     0 0     0 0     0 0     0 0   37.662975   -1.068309   21.734814   906X   70 -0.58500 -0.58500 ACE           
*/              
                
		i = 0;
                while(fgets(lineBuffer, sizeof(lineBuffer), inDAT_file) != NULL  &&  i < nOfAtoms){
                        i++;

                        // get coords        type 1  2  3  4  5  6  7  8  9  10 11 12  x   y   z  
			if(sscanf(lineBuffer, "%d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf",
                            &at,
                            &nrec.neighbor[0], &nrec.bondorder[0],
                            &nrec.neighbor[1], &nrec.bondorder[1],
                            &nrec.neighbor[2], &nrec.bondorder[2],
                            &nrec.neighbor[3], &nrec.bondorder[3],
                            &nrec.neighbor[4], &nrec.bondorder[4],
                            &nrec.neighbor[5], &nrec.bondorder[5],
                            &c.x, &c.y, &c.z
                           ) != 16) 
			{
                            printf("Problem reading .DAT file data at line %d.\n", i); exit(1);
			}
                        
                        // get resname
                        strncpy(rnam, lineBuffer+RESNAMEOFFSET, RESNAMELENGTH - 1);
                        for(j = 0; j < RESNAMELENGTH; j++) if(rnam[j] == ' ') rnam[j] = 0;
                        rnam[RESNAMELENGTH - 1] = 0;

                        // charge
                        if(sscanf(lineBuffer+CHARGEOFFSET, "%lf", &q) != 1) {
                            printf("Problem reading .DAT file data at line %d (charge).\n", i); exit(1);
			}
                        
                        // resnum
                        strncpy(string, lineBuffer+RESNUMOFFSET, 6);
                        string[6] = 0;
                        if(sscanf(string, "%d", &rnum) != 1) {
                            printf("Problem reading .DAT file data at line %d (resnum).\n", i); exit(1);
			}
                        
                        // chain
                        chain = lineBuffer[CHAINOFFSET];
                        
                        aLast = new Atom(c, at, q, i, rnam, rnum, chain, aLast);
                        aLast->setNeighbors(nrec);
                        
                        // water
                        if(strcmp(rnam, "T3P") == 0  ||  strcmp(rnam, "WAT") == 0  ||  strcmp(rnam, "HOH") == 0) aLast->setMolType('w');
                        // membrane
                        else if(strcmp(rnam, "POPC") == 0) aLast->setMolType('m');
                        // heme
                        else if(strcmp(rnam, "HEM") == 0) aLast->setMolType('h');
                        // residues
                        else if(strcmp(rnam, "GLY") == 0  ||
                                strcmp(rnam, "ALA") == 0  ||
                                strcmp(rnam, "VAL") == 0  ||
                                strcmp(rnam, "LEU") == 0  ||
                                strcmp(rnam, "ILE") == 0  ||
                                strcmp(rnam, "PRO") == 0  ||
                                strcmp(rnam, "MET") == 0  ||
                                strcmp(rnam, "CYS") == 0  ||
                                strcmp(rnam, "PHE") == 0  ||
                                strcmp(rnam, "TYR") == 0  ||
                                strcmp(rnam, "TRP") == 0  ||
                                strcmp(rnam, "HIS") == 0  ||
                                strcmp(rnam, "HIP") == 0  ||
                                strcmp(rnam, "ASP") == 0  ||
                                strcmp(rnam, "ASH") == 0  ||
                                strcmp(rnam, "GLU") == 0  ||
                                strcmp(rnam, "GLH") == 0  ||
                                strcmp(rnam, "ASN") == 0  ||
                                strcmp(rnam, "GLN") == 0  ||
                                strcmp(rnam, "ARG") == 0  ||
                                strcmp(rnam, "LYS") == 0  ||
                                strcmp(rnam, "SER") == 0  ||
                                strcmp(rnam, "THR") == 0
                               ) aLast->setMolType('p');
                        // ligand from the command line
                        else if(strcmp(rnam, argv[2]) == 0) {
                            if(ligofi == rnum) aLast->setMolType('L');
                            else aLast->setMolType('l');
                        }
                        
                        // printf("%s %d %c x,y,z={%8.5f %8.5f %8.5f} q=%8.5f type=%c neighbors=%d N=%d\n",
                        //       rnam, rnum, chain, c.x, c.y, c.z, q, aLast->getMolType(), aLast->getNeighCount(), aLast->getElement());
		}
		printf("Read %d atoms.\n", i);
                
                a1 = aLast;
                while(a1) {
                    // if an atom of the ligand
                    if(a1->getMolType() == 'L') {
                        
                        // polar H: -OH, -NH, -SH
                        if(a1->getElement() == HYDROGEN) Xatom = index2atom(aLast, a1->getNeig(0));

                        // printf("%s\n", a1->getResName());
                        // loop on all the rest
                        a2 = aLast;
                        while(a2) {
                            if(a2->getMolType() != 'L') {
                                distance = giveDistance(a1->getCoords(), a2->getCoords());
                                if (distance <= 0.1) { printf("ERROR : distance of two atoms too short\n"); exit(1); }
                                
                                // LIPOPHILIC
                                if(isLipoAtom(a1) == YES  &&  isLipoAtom(a2) == YES) {
                                    e_hphob = give_e_hphob(a1, a2, distance);
                                         if(a2->getMolType() == 'p') tE_hphob_P = tE_hphob_P + e_hphob;
                                    else if(a2->getMolType() == 'l') tE_hphob_L = tE_hphob_L + e_hphob;
                                    else if(a2->getMolType() == 'w') tE_hphob_W = tE_hphob_W + e_hphob;
                                    else if(a2->getMolType() == 'm') tE_hphob_M = tE_hphob_M + e_hphob;
                                    else if(a2->getMolType() == 'h') tE_hphob_H = tE_hphob_H + e_hphob;
                                    tE_hphob = tE_hphob + e_hphob;
                                }

                                // ELECTROSTATIC
                                e_elstat = give_e_elstat(a1, a2, distance);
                                     if(a2->getMolType() == 'p') tE_elstat_P = tE_elstat_P + e_elstat;
                                else if(a2->getMolType() == 'l') tE_elstat_L = tE_elstat_L + e_elstat;
                                else if(a2->getMolType() == 'w') tE_elstat_W = tE_elstat_W + e_elstat;
                                else if(a2->getMolType() == 'm') tE_elstat_M = tE_elstat_M + e_elstat;
                                else if(a2->getMolType() == 'h') tE_elstat_H = tE_elstat_H + e_elstat;
                                tE_elstat = tE_elstat + e_elstat;
                                
                                // H_BONDING
                                // DONOR a1 is H and Xatom holds it (X-H...Y)
                                if(Xatom->getElement() == OXYGEN  ||  Xatom->getElement() == NITROGEN  ||  Xatom->getElement() == SULPHUR) {
                                    // if the other atom is an acceptor
                                    if(a2->getElement() == OXYGEN  ||  a2->getElement() == NITROGEN  ||  a2->getElement() == SULPHUR) {
                                        e_Hbond = giveHbondEnergy(Xatom->getCoords(), a1->getCoords(), a2->getCoords(), Xatom->getElement(), a2->getElement(), CorrectShortHB);
                                             if(a2->getMolType() == 'p') tE_Hbond_P = tE_Hbond_P + e_Hbond;
                                        else if(a2->getMolType() == 'l') tE_Hbond_L = tE_Hbond_L + e_Hbond;
                                        else if(a2->getMolType() == 'w') tE_Hbond_W = tE_Hbond_W + e_Hbond;
                                        else if(a2->getMolType() == 'm') tE_Hbond_M = tE_Hbond_M + e_Hbond;
                                        else if(a2->getMolType() == 'h') tE_Hbond_H = tE_Hbond_H + e_Hbond;
                                        tE_Hbond = tE_Hbond + e_Hbond;
                                        
                                        if(e_Hbond < MIN_EB_FOR_REPORTING) {
                                            printf("H-bond E=%5.2f r=%.2f found between ligand donor atom %d %s %d and %d %s %d\n",
                                                   e_Hbond, distance, a1->getNumber(), a1->getResName(), a1->getResNum(), a2->getNumber(), a2->getResName(), a2->getResNum());
                                        }
                                    }
                                }
                                // ACCEPTOR a1 is acceptor, a2 is H, Yatom holds a2 (H))
                                if(isHbondAcceptor(a1) == YES) {
                                    if(a2->getElement() == HYDROGEN) {
                                        Yatom = index2atom(aLast, a2->getNeig(0));
                                        if (Yatom != NULL) {
                                            if(Yatom->getElement() == OXYGEN  ||  Yatom->getElement() == NITROGEN  ||  Yatom->getElement() == SULPHUR) {
                                                e_Hbond = giveHbondEnergy(Yatom->getCoords(), a2->getCoords(), a1->getCoords(), Yatom->getElement(), a1->getElement(), CorrectShortHB);
                                                     if(a2->getMolType() == 'p') tE_Hbond_P = tE_Hbond_P + e_Hbond;
                                                else if(a2->getMolType() == 'l') tE_Hbond_L = tE_Hbond_L + e_Hbond;
                                                else if(a2->getMolType() == 'w') tE_Hbond_W = tE_Hbond_W + e_Hbond;
                                                else if(a2->getMolType() == 'm') tE_Hbond_M = tE_Hbond_M + e_Hbond;
                                                else if(a2->getMolType() == 'h') tE_Hbond_H = tE_Hbond_H + e_Hbond;
                                                tE_Hbond = tE_Hbond + e_Hbond;

                                                if(e_Hbond < MIN_EB_FOR_REPORTING) {
                                                    printf("H-bond E=%5.2f r=%.2f found between ligand acceptor atom %d %s %d and %d %s %d\n",
                                                           e_Hbond, distance, a1->getNumber(), a1->getResName(), a1->getResNum(), a2->getNumber(), a2->getResName(), a2->getResNum());
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            a2 = a2->getPrevious();
                        }
                    }
                    a1 = a1->getPrevious();
                }

                printf("All analyses done\n\n");

                // printf("totalE_hphob=%6.2f (P=%6.2f  L=%6.2f  W=%6.2f  M=%6.2f  H=%6.2f)\n",
                printf("%s totalE_hphob=%.2f P=%.2f  L=%.2f  W=%.2f  M=%.2f  H=%.2f\n",
                        argv[1], tE_hphob, tE_hphob_P, tE_hphob_L, tE_hphob_W, tE_hphob_M, tE_hphob_H);

                printf("%s totalE_elstat=%.2f P=%.2f  L=%.2f  W=%.2f  M=%.2f  H=%.2f\n",
                        argv[1], tE_elstat, tE_elstat_P, tE_elstat_L, tE_elstat_W, tE_elstat_M, tE_elstat_H);

                printf("%s totalE_Hbond=%.2f P=%.2f  L=%.2f  W=%.2f  M=%.2f  H=%.2f\n",
                        argv[1], tE_Hbond, tE_Hbond_P, tE_Hbond_L, tE_Hbond_W, tE_Hbond_M, tE_Hbond_H);

		printf("Done.\n");

		// C L O S E     F I L E S
		if(fclose(inDAT_file) != 0)  { printf("Error! Problem closing input .mae file.\n"); exit(1); }
	} else {
		printf("Program to calculate MM-GBSA contributions on a MacroModel .dat file.\n");
		printf("usage: myMMGBSA.x  input.dat  ligandName  ligandNumber  correctShortHB\n");
	}
	return 0;
}
