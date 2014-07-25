#include <math.h>
#define PI 3.141592654

double poidev_mod(double xm) // long *idum)
{
	float gammln(float xx);
	//float ran1(long *idum);
	double genrand64_real1(void);
	static float sq,alxm,g,oldm=(-1.0);
	double em,t,y;
	int em_int;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= genrand64_real1();
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*genrand64_real1());
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (genrand64_real1() > t);
	}
	em_int= em;
	return em_int;
}
#undef PI
/* (C) Copr. 1986-92 Numerical Recipes Software "0j. */
