//_____________________________________LIBRERIE______________________________________//
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <iomanip>
using namespace std;

#include "bnldev_mod.c"
#include "poidev_mod.c"
#include "gammln_mod.c"
#include "gammln.c"
#include "mt19937-64.c"

/*Valori dei parametri in SPK:
  N "infinito" (questo algoritmo regge per N<1E+09), s0=0.02, U=1E-05; si dovrebbero vedere a t=1900 un istogramma tra n=1300-1350, a t=1950 un istogramma tra n=1350-1400, a t=2000 un instogramma tra n=1400-1450*/

/*Valori dei parametri in DF:
  grafico <vk> (logN) --> U=1E-05, s=0.01
  grafico <vk> (logU) --> N=1E+06, s=0.01
  grafico <vk> (logs) --> N=1E+06, U=1E-05*/

/*Valori di riferimento nella tesi di Maria per vedere la stazionarieta' della velocita'
  {N=10^8, s=10^-3, U=10-4}*/

//__________________________________________________________________________________________//
//______________________________________*** MAIN ***________________________________________//
//__________________________________________________________________________________________//

int main() {

  init_genrand64 ( time(NULL) );
  ofstream wmeantotFile;
  ofstream kmeantotFile;
  ofstream advtotFile;
  ofstream vktotFile;
  ofstream vstotFile;
  ofstream histoFile;
  ofstream histowFile;
  ofstream numclassFile;
  wmeantotFile.open("wmeantot.dat", ios::trunc);
È bra  kmeantotFile.open("kmeantot.dat", ios::trunc); 
  advtotFile.open("advtot.dat", ios::trunc); 
  vktotFile.open("vktot.dat", ios::trunc);
  vstotFile.open("vstot.dat", ios::trunc);
  histoFile.open("histo.dat", ios::trunc);
  histowFile.open("histow.dat", ios::trunc);
  numclassFile.open("numclass.dat", ios::trunc);

  cout.setf(ios::fixed,ios::floatfield);
  cout.precision(8); 

  //Parametri dal fitting dei dati di Enoch
  /*  double N0 = 3E+10, N = N0;
  double s0 = 0.42, s = s0;
  double slog = 0.297;
  double sgeo = 0.25;
  double alpha = 0.19;
  double qq = 0.60;*/  

 //____Dichiarazione delle variabili globali (iterazioni ensemble)____// 
  int VarPar = 1;
  int iter = 30;
  double N0 = 1E+8, N = N0;
  double s0 = 0.01, s = s0;
  double U0 = 1E-3, U = U0;
  int t0 = 0;
  int tmax = 10001;
  int step = floor(tmax / 20);
  int cell = floor(tmax / step);
  int cc = 0;
  int counter0 = 2*step; int counter = counter0;
  
  //____Dichiarazione delle variabili di popolazione (iterazioni temporali)____//
  int k = 1;
  int kmax = tmax*2;
  int MaxMut;
  double wmean = 0.;
  double kmean = 0.;
  double adv = 0.;
  int numclass = 0;

   //____Dichiarazione delle variabili di classe (iterazioni sulle frequenze)____//
  double nr, Nres = N;
  double pr=0., qr=0., sumpr=0.;
  
  double Vwmeantot[tmax];
  double Vkmeantot[tmax];
  double Vadvtot[tmax];
  double Vvel_k[cell];
  double Vvel_s[cell];
  double Vnumclasstot[tmax];

  int *Vmutation = new int[kmax]; //semplice vettore del # di mutazioni {0,1,2,3...}
  double *Vprob = new double[kmax]; //vettore delle pr probabilita'  delle classi r {p0,p1,p2,p3...}
  double *Vqrob = new double[kmax]; //vettore delle qr prob ridotte delle classi r {q0,q1,q2,q3...}
  double *Vfreq = new double[kmax]; // vettore delle fr frequenze delle classi r {f0,f1,f2,f3...} 
  double *Vfitness = new double[kmax]; //vettore delle wr fitness delle classi r {w0,w1,w2,w3...}
  double *Vadvantage = new double[kmax]; //vettore delle s(k) delle classi r {s0(k),s1(k),s2(k)...}

  double *Vprob_new = new double[kmax]; //vettore delle pr delle classi r {p0,p1,p2,p3...}
  double *Vqrob_new = new double[kmax]; //vettore delle qr delle classi r {q0,q1,q2,q3...}
  double *Vfreq_new = new double[kmax]; //vettore delle fr delle classi r {f0,f1,f2,f3...}

 
  //_____________________________ ||||| FOR SUI PARAMETRI ||||| _____________________________//
  for (int m = 1; m<=VarPar; m++) {

    cout << "Variazione di parametro #" << m << "/" << VarPar << endl;
   
    //___________Inizializzazione della fitness_____________//
    for (int i=0; i<kmax; i++) {
      Vadvantage[i] = s*i;                                   //fitness a vantaggio fisso
      Vfitness[i] = exp(Vadvantage[i]); 
    }
    
       //___________________________ **** FOR SUGLI ENSEMBLE ****_____________________________//
       for (int a=1; a<=iter; a++) {
	 cout << "Iterazione #" << a << "/" << iter << endl;

	 //___________Inizializzazione dei vettori_____________//
	 
	 for (int i=0; i<kmax; i++) {
	   Vmutation[i]=i;
	   Vprob[i]=0.;
	   Vqrob[i]=0.;
	   Vfreq[i]=0.;
	   Vprob_new[i]=0.;
	   Vqrob_new[i]=0.;
	   Vfreq_new[i]=0.;
	 }
	 
	 Vfreq[0] = 1.;
      k = 1;
      counter = counter0; 
      
     //___________________________ °°° FOR TEMPORALE °°° ____________________________//
      for(int t=t0; t<tmax; t++) {
	k++;

	wmean = 0.;
	kmean = 0.;
	adv = 0.;
	for (int i=0; i<k; i++) { wmean += (Vfreq[i]*Vfitness[i]); }
	for (int i=0; i<k; i++) { kmean += (Vfreq[i]*Vmutation[i]); }
	for (int i=0; i<k; i++) { adv += (Vfreq[i]*Vadvantage[i]);}

	//_______________ ^^ FOR SULLE CLASSI DI FREQUENZA ^^ _________________//
	for (int r=0; r<k; r++) {

	  if (r!=0) {
	    pr  = ((Vfreq[r-1] * U * Vfitness[r-1]) + (Vfreq[r] * (1-U) * Vfitness[r]))/wmean;
	  }
	  else 
	    pr = (Vfreq[r] * (1-U) * Vfitness[r])/wmean;

	  sumpr += pr;
	  if (pr!=0.)
	    qr= pr/sumpr;
	  else 
	    qr = 0.;

	  if (pr==0) continue;

	  Vprob_new[r] = pr;
	  Vqrob_new[r] = qr;
	  qr = 0.;
	  MaxMut = k;

	  //cout << "t= " << t << "; k=" << Vmutation[r] << "; fr=" << Vfreq[r] <<  "; pr= " << Vprob_new[r]  << "; sumpr=" << sumpr << "; qr=" << Vqrob_new[r]  << "; wr=" << Vfitness[r] << endl;

	  //___________ + FOR SULLA BINOMIALE + ____________//
	  for (int i = MaxMut; i>=0; i--) {
	    qr = Vqrob_new[i];
		
	    if((qr*Nres) < 3){	
	      if((qr*Nres)< 1){
		if((qr*Nres)< 1E-4){ nr=0;
		}else{ 
		  if((qr*Nres)< 1E-2)
		    {
		      nr= floor(double(poidev_mod(Nres*qr*1000))/1000.);
		    }else{
		    nr= floor(double(poidev_mod(Nres*qr*10))/10.);
		  };
		};
	      }else{
		nr= poidev_mod(Nres*qr);
	      };
	    }else{
	      nr= bnldev_mod(qr, Nres);
	    };
	    /* if((qr*Nres)< 1E-2)
	       nr=0.; 
	       else
	       nr = bnldev_mod(qr, Nres);*/
	  
    	    double freq =nr/N;	  
	    Vfreq_new[i] = freq;
	    Nres -= nr;
	    nr = -1.;
	  } //__ + Chiusura del for sulla binomiale + __ 

	  Vprob[r] = Vprob_new[r];
	  Vqrob[r] = Vqrob_new[r];
	  Vfreq[r] = Vfreq_new[r];

	  //	  if (t==500 || t==5000 || t==25000 || t==50000 || t==100000 || t == 150000) {
	    if (t == counter) {
	    histoFile << r << " " << Vfreq[r] << endl;
	    histowFile << Vfitness[r] << " " << Vfreq[r] << endl;
	  }

	  //Vfreq[r] = pr; //Scommentare questa e commentare il for sulla binomiale per l'algoritmo deterministico
	  pr = 0.;
	  qr = 0.;
	  Nres = N;

	  if (Vfreq[r] != 0) {numclass++;}

	  if (r>0   &&   Vprob[r-1] != 0   &&   Vprob[r] == 0) break;

	}; // __ ^^ Chiusura del for sulle classi di frequenza ^^ __
 
	  if (t == counter) {
	    cout << t << endl;
	  histoFile << endl;
	  histowFile << endl;
	  // counter += counter0;
	}

	  if (t == counter)
	    {histoFile << endl;
	      counter += counter0;
	    }

	Vwmeantot[t] += wmean;
	Vkmeantot[t] += kmean;
	Vadvtot[t] += adv;
	Vnumclasstot[t] += (double) numclass;
	sumpr = 0.;
	numclass = 0; 

	//			if (t % 500 == 0) {
	//	  cout << t << " " << kmean << endl;
	//	}
	//cout << endl;

      }; // __ °°° Chiusura del for temporale °°° __
    
   
  
    };  // __ **** Chiusura della singola iterazione **** __

    cc = 0; 
   
    for (int n=0; n<tmax; n++) {
      Vwmeantot[n] /= iter;
      Vkmeantot[n] /= iter;
      Vadvtot[n] /= iter;
      Vnumclasstot[n] /= iter;
       
      if (n % step == 0) {
	wmeantotFile << n << " " << Vwmeantot[n] << endl;
	kmeantotFile << n << " " << Vkmeantot[n] << endl;
	advtotFile << Vkmeantot[n] << " " << Vadvtot[n] << endl;
	numclassFile << n << " " << Vnumclasstot[n]/2 << endl;

	if (n>0) {
	  Vvel_k[cc] = (Vkmeantot[n] - Vkmeantot[n-step])/step;
	  vktotFile << Vkmeantot[n] << " " << Vvel_k[cc] << endl;
	  Vvel_s[cc] = (Vadvtot[n] - Vadvtot[n-step])/step;
	  vstotFile << Vkmeantot[n] << " " << Vvel_s[cc] << endl;
	  if (cc == cell){
	    vktotFile << endl;
	    vstotFile << endl;
	  };
	};	 
	  cc++;
	}

      if (n == (tmax-1)) {
	wmeantotFile << endl;
	kmeantotFile << endl;
	advtotFile << endl;
      }
    }


  vktotFile << endl;

  //v teorica stochastic edge
  // for (int i=0; i < cell; i++){ 
  //vtheo = s*s*(2*log(N*s) - log(s/U))/(log(s/U)*log(s/U));
  //vktotFile << i*step << " " << vtheo << endl; 
  //}
  // vktotFile << endl;

// v teorica dalla gaussiana
  //    for (int i=0; i < cell; i++){ 
  //  vgauss = 2*s*s*log(N)/(log(U)*log(U));
  //vktotFile << i*step << " " << vgauss << endl; 
  //}


  //cout << U << " " << vtheo << endl;
  U *= 2;

  }; // __ ||||| Chiusura del for sui parametri ||||| __
 
  histoFile.close();
  histowFile.close();
  wmeantotFile.close();
  kmeantotFile.close();
  advtotFile.close();
  vktotFile.close();
  vstotFile.close();
  numclassFile.close();

  delete Vmutation;
  delete Vprob;
  delete Vqrob;
  delete Vfreq;
  delete Vfitness;
  delete Vadvantage;
  delete Vprob_new;
  delete Vqrob_new;
  delete Vfreq_new;
 
  return 4444;
  
}//Chiusura del main
