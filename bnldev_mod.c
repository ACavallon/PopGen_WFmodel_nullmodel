#include <math.h>
#define PI 3.141592654


//___________________________________________________________________________
//_________________________NBLDEV_MOD: BINOMIALE_____________________________
//___________________________________________________________________________

double bnldev_mod(double pp, long long int n) // , long *idum)
{
	double gammln_mod(double xx);
	double genrand64_real1(void); //added by me
	int j;
	static long long int nold=(-1);
	double am,em,g,angle,p,bnl,sq,t,y;
	static double pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (genrand64_real1() < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= genrand64_real1();
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln_mod(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*genrand64_real1();
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln_mod(em+1.0)
				-gammln_mod(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (genrand64_real1() > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}
#undef PI
/* (C) Copr. 1986-92 Numerical Recipes Software "0j. */
