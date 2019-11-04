/* 
 * File:   points_vectors.h
 * Author: Martin Smiesko
 *
 * Created on June 19, 2014, 5:02 PM
 */

#ifndef POINTS_VECTORS_H
#define	POINTS_VECTORS_H

/**********************************************/
/*                                            */
/*    P O I N T S   A N  D   V E C T O R S    */
/*                                            */
/**********************************************/

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

double	giveVectorLength(Vector v1) { return sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z); }

double	giveVectorLengthSquared(Vector v1) { return (v1.x*v1.x + v1.y*v1.y + v1.z*v1.z); }

Vector	giveVector(Point A, Point B) {
	Vector	v;
	v.x = B.x - A.x;
	v.y = B.y - A.y;
	v.z = B.z - A.z;
	return v;
}

double	giveDistance(Point A, Point B) {
        A = A - B;
        return sqrt(A.x*A.x + A.y*A.y + A.z*A.z);
}

double	giveDistanceSquared(Point A, Point B) {
        A = A - B;
        return (A.x*A.x + A.y*A.y + A.z*A.z);
}

double giveAngle_SLOW(Point A, Point B, Point C) {
	double	a, b, c;
	a = giveVectorLength(point2vector(B - C));
	b = giveVectorLength(point2vector(A - C));
	c = giveVectorLength(point2vector(B - A));
	return acos((b * b - a * a - c * c) / (-2.0 * a * c));
}

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

double	giveVectorAngle(Vector v1, Vector v2) {
	double	asq, bsq, csq, arg;
	asq = giveVectorLengthSquared(v1 - v2);
        bsq = giveVectorLengthSquared(v1);
	csq = giveVectorLengthSquared(v2);
        arg = (asq - bsq - csq) / (-2.0 * sqrt(bsq * csq));
        if(arg < -1.0) arg = -1.0;
        else if(arg > 1.0) arg = 1.0;
        return acos(arg);
}

Vector	giveVectorProduct(Vector v1, Vector v2) {
	Vector v;
	v.x = v1.y * v2.z - v1.z * v2.y;
	v.y = v1.z * v2.x - v1.x * v2.z;
	v.z = v1.x * v2.y - v1.y * v2.x;
	return v;
}

double	giveDotProduct(Vector v1, Vector v2) { return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z); }

Vector	giveNormalVector(Point A, Point  B, Point C) { return giveVectorProduct(giveVector(B, A), giveVector(B, C)); }

Vector	multiplyVector(double factor, Vector v) {
	v.x = factor * v.x;
	v.y = factor * v.y;
	v.z = factor * v.z;
	return v;
}

Vector	giveUnitVector(Vector v1){ return multiplyVector(1 / giveVectorLength(v1), v1); }

double	giveTorsion(Point A, Point  B, Point C, Point D) {
	Vector	v1 = giveVector(A, B);
	Vector	v2 = giveVector(B, C);
	Vector	v3 = giveVector(C, D);
	
	Vector	term1 = multiplyVector(giveVectorLength(v2), v1);
	Vector	term2 = giveVectorProduct(v2, v3);
	Vector	term3 = giveVectorProduct(v1, v2);

	return atan2(giveDotProduct(term1, term2) , giveDotProduct(term3, term2));
}

Vector	giveVectorSum(Vector v1, Vector v2) {
	Vector v;
	v.x = v1.x + v2.x;
	v.y = v1.y + v2.y;
	v.z = v1.z + v2.z;
	return v;
}

/***************************************************/
/*                                                 */
/*       A N G L E S   &   R O T A T I O N S       */
/*                                                 */
/***************************************************/

double	giveAngleX(Point p){ return atan2(p.z, p.y); }

double	giveAngleY(Point p){ return atan2(p.z, p.x); }

double	giveAngleZ(Point p){ return atan2(p.y, p.x); }

Point	rotX(Point p, double alpha){
	double	Y, Z, ca, sa;
        
        sincos(alpha, &sa, &ca);
	
	Y = p.y * ca - p.z * sa;
	Z = p.y * sa + p.z * ca;
	p.y = Y;
	p.z = Z;		

	return p;	
}

Point	rotY(Point p, double alpha){
	double	X, Z, ca, sa;

        sincos(alpha, &sa, &ca);
	
	X = p.x * ca - p.z * sa;
	Z = p.x * sa + p.z * ca;
	p.x = X;
	p.z = Z;
		
	return p;
}

Point	rotZ(Point p, double alpha){
	double	X, Y, ca, sa;

        sincos(alpha, &sa, &ca);
	
	X = p.x * ca - p.y * sa;
	Y = p.x * sa + p.y * ca;
	p.x = X;
	p.y = Y;

	return p;
}

Point	generatePoint(Point A, Point B, Point C, double r, double alpha, double phi = 0){
	Point	D;
	double	x, y, z;
        double  sa, ca, sp, cp;
	
	// translate points by B
	A = A - B;
	C = C - B;

	// align BC to [ 1, 0, 0 ]
	// get angle along Z
	z = giveAngleZ(C);
	// rotate points
	A = rotZ(A, -z);
	C = rotZ(C, -z);
	
	// align BC to [ 1, 0, 0 ]
	// get angle along Y
	y = giveAngleY(C);
	// rotate points
	A = rotY(A, -y);
	C = rotY(C, -y);

	// align A to plane xy
	// get angle along X
	x = giveAngleX(A);
	// rotate points
	A = rotX(A, -x);
	// D O N E

	// initialize the point
	alpha = RAD_180 - alpha;
	
        sincos(alpha, &sa, &ca);
        sincos(phi, &sp, &cp);
        
        D.x = C.x + r * ca;
	D.y =       r * sa * cp;
	D.z =       r * sa * sp;
        
	// N O W   G O   B A C K W A R D S
	// rotate back by X
	D = rotX(D, x);
	// rotate back by Y
	D = rotY(D, y);
	// rotate back by Z
	D = rotZ(D, z);
	// translate the point back by B
	D = D + B;

	return D;
}

/***********************/
/*                     */
/*    P H Y S I C S    */
/*                     */
/***********************/

// linearity corrected (X-H...Y) H-bond energy function from teh Yeti force field
double	giveHbondEnergy(Point X, Point H, Point Y, unsigned char Nx, unsigned char Ny){
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


			C = -5 * wellDepth * pow(rOpt, 12);
			D = -6 * wellDepth * pow(rOpt, 10);

			cosine = cos(angleXHY);
			cosine = cosine * cosine;
			
			return ( cosine * (C/pow(radius, 12) - D/pow(radius, 10)) );
		}
	}
	return 0;
}

#endif	/* POINTS_VECTORS_H */

