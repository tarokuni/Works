#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 15
#define M 5
#define KINL 10
#define LAMBDA 0.01
#define ALPHAI 1.1e132
double ALPHA = ALPHAI; //(M_Pl/M)^(2n+4m)
#define KIZAMIA 1.01
#define X2 1. //X^2 on/off
#define BETAII 2.7e161 //(M_Pl/Mtilde)^4(ell-1)
#define KIZAMIB 1.01
#define EPS 1.0e-7//
#define PHIINII 0.03
#define DPHIINI -4.0e-10
double PHIINI = PHIINII;
#define HINI 0.00001
#define HLIM 1.23e-8 //
#define OSC 30
#define EFOLD 60
int NUM = 1500;
double BETAI = BETAII;
double BETA = BETAII;

#define A2 1./5.
#define A3 3./10.
#define A4 3./5.
#define A5 1.
#define A6 7./8.

#define B21 1./5.
#define B31 3./40.
#define B32 9./40.
#define B41 3./10.
#define B42 (-9./10.)
#define B43 6./5.
#define B51 (-11./54.)
#define B52 5./2.
#define B53 (-70./27.)
#define B54 35./27.
#define B61 1631./55296.
#define B62 175./512.
#define B63 575./13824.
#define B64 44275./110592.
#define B65 253./4096.

#define C1 37./378.
#define C2 0.
#define C3 250./621.
#define C4 125./594.
#define C5 0.
#define C6 512./1771.

#define CS1 2825./27648.
#define CS2 0.
#define CS3 18575./48384.
#define CS4 13525./55296.
#define CS5 277./14336.
#define CS6 1./4.


FILE *out1;

double powint(double phi, int i);
double EOM(double phi, double dphi, double H);
double EOM5(double phi, double dphi, double H);
double Hcon(double phi, double dphi);
double HEOM(double phi, double dphi, double ddphi, double H);
double Fs(double phi, double dphi, double ddphi, double Hubble);
double Gs(double phi, double dphi, double ddphi, double H);
double cs2(double phi, double dphi, double ddphi, double Hubble);
int RK5(double phi, double dphi, double lna, double H);
double maxval(double x, double y, double z, double w);

int main(void) {
	int x;
	double phi, dphi, lna, H;

	out1 = fopen("./out28.dat", "w");

	lna = 0;
	phi = PHIINI;
	dphi = DPHIINI;
	H = Hcon(phi, dphi);

	while (1) {
		while (1) {
			x = RK5(phi, dphi, lna, H);
			if (x == 13) {
				phi *= 1.3;
				printf("phi=%e\n", phi);
				H = Hcon(phi, dphi);
			}
			if (x == 18) {
				phi *= 0.81;
				printf("phi=%e\n", phi);
				H = Hcon(phi, dphi);
			}
			if (x == 0) {
				BETA *= KIZAMIB;
			}
			if (x == 1) {
				//				break;
				return 0;
			}
			if (x == 2) {
				BETA = BETAI;
				ALPHA *= KIZAMIA;
				break;
			}
			if (x == 3) {
				BETAI *= 0.99;
				BETA = BETAI;
				break;
			}
			if (x == 4) {
				ALPHA/=KIZAMIA;
				BETAI *= 0.99;
				BETA = BETAI;
				break;
			}
		}

		//		printf("\n");
	}

	fclose(out1);

	return 0;
}


int RK5(double phi, double dphi, double lna, double H) {

	double k1A, k1B, k1C, k1D, k2A, k2B, k2C, k2D, k3A, k3B, k3C, k3D, k4A, k4B, k4C, k4D, k5A, k5B, k5C, k5D, k6A, k6B, k6C, k6D;
	double phis, dphis, lnas, Hs, Dphi, Ddphi, Dlna, DH, t, h, Delta, ts, lnaa, Pcalc, dlnass;
	int i, j, k, l, p, q, s, osc, n, m;

	double lnass[5001], css[5001], HHH[5001], Fss[5001], epss[5001], phiss[5001];

	j = 1;
	k = 0;
	l = 1;
	p = 1;
	osc = 0;

	phis = phi;
	dphis = dphi;
	lnas = lna;
	Hs = H;
	Pcalc = 1.;
	dlnass = 0;

	for (q = 0; q <= 5000; q++) {
		lnass[q] = lna;
		epss[q] = -HEOM(phi, dphi, EOM(phi, dphi, H), H) / H / H;
		css[q] = 1.;
		HHH[q] = H;
		Fss[q] = Fs(phi, dphi, EOM(phi, dphi, H), H);
		phiss[q] = phi;
	}

	//	printf("lnass=%e, epss=%e, css=%e, HHH=%e, Fss=%e, phiss=%e\n", lnass[NUM], epss[NUM], css[NUM], HHH[NUM], Fss[NUM],phiss[NUM]);

	lnaa = lna;

	Dphi = 0;
	Ddphi = 0;
	Dlna = 0;
	DH = 0;
	Delta = 0;

	h = HINI;
	t = 0;
	ts = 1.0;

	for (i = 0;; i++) {


		if (p) {

			if (lna >= lnaa + 0.02) {

				s = 1;

				while (s <= (lna - lnaa - 0.02)*50.) {
					s += 1;
				}

				for (q = 0; q <= 5000; q++) {
					if (q + s <= 5000) {
						lnass[q] = lnass[q + s];
						epss[q] = epss[q + s];
						css[q] = css[q + s];
						HHH[q] = HHH[q + s];
						Fss[q] = Fss[q + s];
						phiss[q] = phiss[q + s];
					}
					else {
						lnass[q] = lna;
						epss[q] = -HEOM(phi, dphi, EOM(phi, dphi, H), H) / H / H;
						css[q] = cs2(phi, dphi, EOM(phi, dphi, H), H);
						HHH[q] = H;
						Fss[q] = Fs(phi, dphi, EOM(phi, dphi, H), H);
						phiss[q] = phi;
					}

				}
				lnaa = lna;
			}



		}

		if (-HEOM(phi, dphi, EOM(phi, dphi, H), H) / H / H >= 1 && p) {
			lnaa = lna;
			p = 0;
			//				printf("N=%e, N2=%e, NUM=%d\n", lnaa - lnass[NUM], lnaa - lnass[NUM+1], NUM);

			n = 1;
			while (1) {
				dlnass = (lnass[NUM + n] - lnass[NUM - n]);
				if (dlnass != 0) { break; }
//					break;
				n++;
			}



			Pcalc = pow(2., (4 * epss[NUM] + 3 * (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] - (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]))
				/ (2 - 2 * epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM])))
				*(1. - epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] / 2. + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]) / 2.)
				*HHH[NUM] * HHH[NUM] / sqrt(css[NUM]) / Fss[NUM] / 8. / 3.141593 / 3.141593;
//				printf("Inflation end\n");
//				printf("N=%e, P=%e\n", lna, Pcalc);


			if (lna < 1.2*EFOLD) { return 13; }
			if (lna > 8. * EFOLD) { return 18; }
			if (Pcalc > 2.2e-9) {
					printf("N=%e, P=%e\n", lna, Pcalc);
//					printf("\n");
				break;
			}
		}

		for (;;) {

			k1A = h*dphi;
			k1B = h*EOM(phi, dphi, H);
			k1C = h*H;
			k1D = h*HEOM(phi, dphi, EOM(phi, dphi, H), H);

			k2A = h*(dphi + k1B*B21);
			k2B = h*EOM(phi + k1A*B21, dphi + k1B*B21, H + k1D*B21);
			k2C = h*(H + k1D*B21);
			k2D = h*HEOM(phi + k1A*B21, dphi + k1B*B21, EOM(phi + k1A*B21, dphi + k1B*B21, H + k1D*B21), H + k1D*B21);

			k3A = h*(dphi + k2B*B32 + k1B*B31);
			k3B = h*EOM(phi + k2A*B32 + k1A*B31, dphi + k2B*B32 + k1B*B31, H + k2D*B32 + k1D*B31);
			k3C = h*(H + k2D*B32 + k1D*B31);
			k3D = h*HEOM(phi + k2A*B32 + k1A*B31, dphi + k2B*B32 + k1B*B31, EOM(phi + k2A*B32 + k1A*B31, dphi + k2B*B32 + k1B*B31, H + k2D*B32 + k1D*B31), H + k2D*B32 + k1D*B31);

			k4A = h*(dphi + k3B*B43 + k2B*B42 + k1B*B41);
			k4B = h*EOM(phi + k3A*B43 + k2A*B42 + k1A*B41, dphi + k3B*B43 + k2B*B42 + k1B*B41, H + k3D*B43 + k2D*B42 + k1D*B41);
			k4C = h*(H + k3D*B43 + k2D*B42 + k1D*B41);
			k4D = h*HEOM(phi + k3A*B43 + k2A*B42 + k1A*B41, dphi + k3B*B43 + k2B*B42 + k1B*B41,
				EOM(phi + k3A*B43 + k2A*B42 + k1A*B41, dphi + k3B*B43 + k2B*B42 + k1B*B41, H + k3D*B43 + k2D*B42 + k1D*B41), H + k3D*B43 + k2D*B42 + k1D*B41);

			k5A = h*(dphi + k4B*B54 + k3B*B53 + k2B*B52 + k1B*B51);
			k5B = h*EOM(phi + k4A*B54 + k3A*B53 + k2A*B52 + k1A*B51, dphi + k4B*B54 + k3B*B53 + k2B*B52 + k1B*B51, H + k4D*B54 + k3D*B53 + k2D*B52 + k1D*B51);
			k5C = h*(H + k4D*B54 + k3D*B53 + k2D*B52 + k1D*B51);
			k5D = h*HEOM(phi + k4A*B54 + k3A*B53 + k2A*B52 + k1A*B51, dphi + k4B*B54 + k3B*B53 + k2B*B52 + k1B*B51,
				EOM(phi + k4A*B54 + k3A*B53 + k2A*B52 + k1A*B51, dphi + k4B*B54 + k3B*B53 + k2B*B52 + k1B*B51, H + k4D*B54 + k3D*B53 + k2D*B52 + k1D*B51), H + k4D*B54 + k3D*B53 + k2D*B52 + k1D*B51);

			k6A = h*(dphi + k5B*B65 + k4B*B64 + k3B*B63 + k2B*B62 + k1B*B61);
			k6B = h*EOM(phi + k5A*B65 + k4A*B64 + k3A*B63 + k2A*B62 + k1A*B61, dphi + k5B*B65 + k4B*B64 + k3B*B63 + k2B*B62 + k1B*B61, H + k5D*B65 + k4D*B64 + k3D*B63 + k2D*B62 + k1D*B61);
			k6C = h*(H + k5D*B65 + k4D*B64 + k3D*B63 + k2D*B62 + k1D*B61);
			k6D = h*HEOM(phi + k5A*B65 + k4A*B64 + k3A*B63 + k2A*B62 + k1A*B61, dphi + k5B*B65 + k4B*B64 + k3B*B63 + k2B*B62 + k1B*B61,
				EOM(phi + k5A*B65 + k4A*B64 + k3A*B63 + k2A*B62 + k1A*B61, dphi + k5B*B65 + k4B*B64 + k3B*B63 + k2B*B62 + k1B*B61, H + k5D*B65 + k4D*B64 + k3D*B63 + k2D*B62 + k1D*B61), H + k5D*B65 + k4D*B64 + k3D*B63 + k2D*B62 + k1D*B61);

			Dphi = (k1A*C1 + k2A*C2 + k3A*C3 + k4A*C4 + k5A*C5 + k6A*C6) - (k1A*CS1 + k2A*CS2 + k3A*CS3 + k4A*CS4 + k5A*CS5 + k6A*CS6);
			Ddphi = (k1B*C1 + k2B*C2 + k3B*C3 + k4B*C4 + k5B*C5 + k6B*C6) - (k1B*CS1 + k2B*CS2 + k3B*CS3 + k4B*CS4 + k5B*CS5 + k6B*CS6);
			Dlna = (k1C*C1 + k2C*C2 + k3C*C3 + k4C*C4 + k5C*C5 + k6C*C6) - (k1C*CS1 + k2C*CS2 + k3C*CS3 + k4C*CS4 + k5C*CS5 + k6C*CS6);
			DH = (k1D*C1 + k2D*C2 + k3D*C3 + k4D*C4 + k5D*C5 + k6D*C6) - (k1D*CS1 + k2D*CS2 + k3D*CS3 + k4D*CS4 + k5D*CS5 + k6D*CS6);

			Delta = maxval(fabs(Dphi / (fabs(k1A*C1) + fabs(k2A*C2) + fabs(k3A*C3) + fabs(k4A*C4) + fabs(k5A*C5) + fabs(k6A*C6) + fabs(k1A*CS1) + fabs(k2A*CS2) + fabs(k3A*CS3) + fabs(k4A*CS4) + fabs(k5A*CS5) + fabs(k6A*CS6))),
				fabs(Ddphi / (fabs(k1B*C1) + fabs(k2B*C2) + fabs(k3B*C3) + fabs(k4B*C4) + fabs(k5B*C5) + fabs(k6B*C6) + fabs(k1B*CS1) + fabs(k2B*CS2) + fabs(k3B*CS3) + fabs(k4B*CS4) + fabs(k5B*CS5) + fabs(k6B*CS6))),
				fabs(Dlna / (fabs(k1C*C1) + fabs(k2C*C2) + fabs(k3C*C3) + fabs(k4C*C4) + fabs(k5C*C5) + fabs(k6C*C6) + fabs(k1C*CS1) + fabs(k2C*CS2) + fabs(k3C*CS3) + fabs(k4C*CS4) + fabs(k5C*CS5) + fabs(k6C*CS6))),
				fabs(DH / (fabs(k1D*C1) + fabs(k2D*C2) + fabs(k3D*C3) + fabs(k4D*C4) + fabs(k5D*C5) + fabs(k6D*C6) + fabs(k1D*CS1) + fabs(k2D*CS2) + fabs(k3D*CS3) + fabs(k4D*CS4) + fabs(k5D*CS5) + fabs(k6D*CS6)))) / EPS;


			if (Delta >= 1) {
				h = ((0.95*h*pow(Delta, -0.2) >= 0.2*h) ? 0.95*h*pow(Delta, -0.2) : 0.2*h);
				if (h<HLIM) {
					h = HLIM;
					break;
				}
			}
			else {
				break;
			}

		}

		t += h;

		phi += k1A*C1 + k2A*C2 + k3A*C3 + k4A*C4 + k5A*C5 + k6A*C6;
		dphi += k1B*C1 + k2B*C2 + k3B*C3 + k4B*C4 + k5B*C5 + k6B*C6;
		lna += k1C*C1 + k2C*C2 + k3C*C3 + k4C*C4 + k5C*C5 + k6C*C6;
		H += k1D*C1 + k2D*C2 + k3D*C3 + k4D*C4 + k5D*C5 + k6D*C6;

		phis += k1A*CS1 + k2A*CS2 + k3A*CS3 + k4A*CS4 + k5A*CS5 + k6A*CS6;
		dphis += k1B*CS1 + k2B*CS2 + k3B*CS3 + k4B*CS4 + k5B*CS5 + k6B*CS6;
		lnas += k1C*CS1 + k2C*CS2 + k3C*CS3 + k4C*CS4 + k5C*CS5 + k6C*CS6;
		Hs += k1D*CS1 + k2D*CS2 + k3D*CS3 + k4D*CS4 + k5D*CS5 + k6D*CS6;

//		printf("%e", phi);

		if (Delta != 0.0) {
			h *= ((0.95*pow(Delta, -0.25) <= 5.) ? 0.95*pow(Delta, -0.25) : 5.);
		}

		if (Delta == 0.0) {
			h *= 5.;
		}


		if (j && (phi <= 0)) {
			j = 0;
		}

		if ((!j) && (phi >= 0)) {
			j = 1;
			osc += 1;
		}

		//			if(phi>2.*PHIINI|phi<-2.*PHIINI){
		//				printf("”­ŽUAphi=%e\n", phi);
		//				break;
		//			}

		if (p) {
			while (lnaa - lnass[NUM] > (EFOLD + 0.02)) { NUM += 1; }
			while (NUM >= 10 && lnaa - lnass[NUM] < (EFOLD - 0.02)) { NUM -= 1; } 
		}

		if (osc >= OSC) {
			break;
		}

		if (!p) {
			if (cs2(phi, dphi, EOM(phi, dphi, H), H) < 0) {
				//					fprintf(out1, "cs2<0\n");
				printf("cs2<0, LAMBDA=%e, ALPHA=%e, BETA=%e, N=%e\n", LAMBDA, ALPHA, BETA, lnaa);
				return 1;
			}

			if (EOM5(phi, dphi, H) < 0) {
				//					fprintf(out1, "EOM breakdown, ALPHA=%e, BETA=%e\n", ALPHA, BETA);
				printf("EOM breakdown, LAMBDA=%e, ALPHA=%e, BETA=%e\n", LAMBDA, ALPHA, BETA);
				return 1;
			}
		}





	}




	n = 1;
	while (1) {
		dlnass = (lnass[NUM + n] - lnass[NUM - n]);
		if (dlnass != 0.0) { break; }
		//			break;
		//			printf("%d, %e\n",n, dlnass);
		n++;
	}

	//		printf("%d, %e\n", n, dlnass);




	Pcalc = pow(2., (4 * epss[NUM] + 3 * (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] - (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]))
		/ (2 - 2 * epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM])))
		*(1. - epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] / 2. + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]) / 2.)
		*HHH[NUM] * HHH[NUM] / sqrt(css[NUM]) / Fss[NUM] / 8. / 3.141593 / 3.141593;

	if (Pcalc <= 2.2e-9) {



		if (BETA == BETAI) {
			//				printf("BETA too large\n");
			fprintf(out1, "%e, %e, %e, %f, %e, %e, %e, %e, %e, %e, %e\n", lna, ALPHA, BETA, lnaa - lnass[NUM],
				/*P*/(pow(2., (4 * epss[NUM] + 3 * (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] - (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]))
					/ (2 - 2 * epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM])))
					*(1. - epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] / 2. + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]) / 2.))
				*HHH[NUM] * HHH[NUM] / sqrt(css[NUM]) / Fss[NUM] / 8. / 3.141593 / 3.141593,
				/*ns*/1. - (4 * epss[NUM] + 3 * (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] - (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]))
				/ (2 - 2 * epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM])),
				/*r*/16 * Fss[NUM] * sqrt(css[NUM]), epss[NUM], (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM],
				(Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]), dlnass);
			printf("oscillation=%d, ALPHA=%e, BETA=%e, N=%f, P=%e, ns=%e, dlnass=%e\nr=%e, NUM=%d, lna=%e\n%e, %e, %e, %e\n", osc, ALPHA, BETA, lnaa - lnass[NUM],
				(pow(2., (4 * epss[NUM] + 3 * (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] - (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]))
					/ (2 - 2 * epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM])))
					*(1. - epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] / 2. + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]) / 2.))
				*HHH[NUM] * HHH[NUM] / sqrt(css[NUM]) / Fss[NUM] / 8. / 3.141593 / 3.141593,
				1. - (4 * epss[NUM] + 3 * (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] - (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]))
				/ (2 - 2 * epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM])), dlnass,
				16 * Fss[NUM] * sqrt(css[NUM]), NUM, lna, epss[NUM], (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM],
				(Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]), css[NUM]);
			if (Pcalc <= 2.19e-9) { return 4; }
			return 3;
		}



		fprintf(out1, "%e, %e, %e, %f, %e, %e, %e, %e, %e, %e, %e\n", lna, ALPHA, BETA, lnaa - lnass[NUM],
			/*P*/(pow(2., (4 * epss[NUM] + 3 * (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] - (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]))
				/ (2 - 2 * epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM])))
				*(1. - epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] / 2. + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]) / 2.))
			*HHH[NUM] * HHH[NUM] / sqrt(css[NUM]) / Fss[NUM] / 8. / 3.141593 / 3.141593,
			/*ns*/1. - (4 * epss[NUM] + 3 * (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] - (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]))
			/ (2 - 2 * epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM])),
			/*r*/16 * Fss[NUM] * sqrt(css[NUM]), epss[NUM], (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM],
			(Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]), dlnass);
		printf("oscillation=%d, ALPHA=%e, BETA=%e, N=%f, P=%e, ns=%e, dlnass=%e\nr=%e, NUM=%d, lna=%e\n%e, %e, %e, %e\n", osc, ALPHA, BETA, lnaa - lnass[NUM],
			(pow(2., (4 * epss[NUM] + 3 * (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] - (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]))
				/ (2 - 2 * epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM])))
				*(1. - epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] / 2. + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]) / 2.))
			*HHH[NUM] * HHH[NUM] / sqrt(css[NUM]) / Fss[NUM] / 8. / 3.141593 / 3.141593,
			1. - (4 * epss[NUM] + 3 * (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] - (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]))
			/ (2 - 2 * epss[NUM] - (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM] + (Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM])), dlnass,
			16 * Fss[NUM] * sqrt(css[NUM]), NUM, lna, epss[NUM], (Fss[NUM + n] - Fss[NUM - n]) / dlnass / Fss[NUM],
			(Fss[NUM + n] / css[NUM + n] - Fss[NUM - n] / css[NUM - n]) / dlnass / (Fss[NUM] / css[NUM]), css[NUM]);
		return 2;
	}


	return 0;
}


double HEOM(double phi, double dphi, double ddphi, double H) {
	return -1. / 2.*(dphi*dphi / 2. + X2*(BETA*powint(dphi, 2.*KINL) / KINL / powint(2., KINL)) - LAMBDA*phi*phi*phi*phi / 4. + (2 * N + 1)*ALPHA*powint(phi, 2 * N)*powint(dphi, 2 * M + 2) / powint(2., M)
		+ 2 * M*ALPHA*powint(phi, 2 * N + 1)*powint(dphi, 2 * M)*ddphi / powint(2., M) + 3 * H*H);
}


double Hcon(double phi, double dphi) {
	return (-6.*ALPHA*M*powint(phi, 2 * N + 1)*powint(dphi, 2 * M + 1) / powint(2., M) +
		sqrt(9 * ALPHA*ALPHA*M*M*powint(phi, 4 * N + 2)*powint(dphi, 4 * M + 2) / powint(2., 2 * M - 2)
			+ 6 * dphi*dphi + 3 * LAMBDA*phi*phi*phi*phi + X2*(BETA * 12 * (2 * KINL - 1)*powint(dphi, 2 * KINL) / powint(2., KINL) / KINL) + 6.*(2 * N + 1)*ALPHA*powint(phi, 2 * N)*powint(dphi, 2 * M + 2) / powint(2., M - 1))) / 6.;
}

double EOM(double phi, double dphi, double H) {
	return -1.*(3.*H*dphi + LAMBDA*phi*phi*phi - 9.*M*ALPHA*H*H*powint(phi, 2 * N + 1)*powint(dphi, 2 * M) / powint(2., M)
		- 3.*ALPHA*(M - 1)*(2 * N + 1)*H*powint(phi, 2 * N)*powint(dphi, 2 * M + 1) / powint(2., M - 1) + 3.*M*ALPHA*powint(phi, 2 * N + 1)*powint(dphi, 2 * M + 2) / powint(2., M + 1)
		- 3.*M*ALPHA*LAMBDA*powint(phi, 2 * N + 5)*powint(dphi, 2 * M) / powint(2., M + 2) + (2 * N + 1)*N*ALPHA*powint(phi, 2 * N - 1)*powint(dphi, 2 * M + 2) / powint(2., M - 1)
		+ 3.*M*(2 * N + 1)*ALPHA*ALPHA*powint(phi, 4 * N + 1)*powint(dphi, 4 * M + 2) / powint(2., 2 * M) + X2*(3.*BETA*H*powint(dphi, 2 * KINL - 1) / powint(2, KINL - 1) + 3.*BETA*M*ALPHA*powint(phi, 2 * N + 1)*powint(dphi, 2 * M + 2 * KINL) / powint(2., KINL + M) / KINL))
		/ (1 + 3.*M*M*ALPHA*ALPHA*powint(phi, 4 * N + 2)*powint(dphi, 4 * M) / powint(2., 2 * M - 1) - 6 * M*M*ALPHA*H*powint(phi, 2 * N + 1)*powint(dphi, 2 * M - 1) / powint(2., M - 1)
			+ (2 * N + 1)*(M + 1)*ALPHA*powint(phi, 2 * N)*powint(dphi, 2 * M) / powint(2., M - 1) + X2*(BETA*(2.*KINL - 1)*powint(dphi, 2 * (KINL - 1)) / powint(2., KINL - 1)));
}

double EOM5(double phi, double dphi, double H) {
	return 1 + 3.*M*M*ALPHA*ALPHA*powint(phi, 4 * N + 2)*powint(dphi, 4 * M) / powint(2., 2 * M - 1) - 6 * M*M*ALPHA*H*powint(phi, 2 * N + 1)*powint(dphi, 2 * M - 1) / powint(2., M - 1)
		+ (2 * N + 1)*(M + 1)*ALPHA*powint(phi, 2 * N)*powint(dphi, 2 * M) / powint(2., M - 1) + X2*(BETA*(2.*KINL - 1)*powint(dphi, 2 * (KINL - 1)) / powint(2., KINL - 1));
}



double Fs(double phi, double dphi, double ddphi, double H) {
	return (1 - (2 * N + 1)*(M - 1)*ALPHA*powint(phi, 2 * N)*powint(dphi, 2 * M) / powint(2., M - 1) - 4.*M*ALPHA*H*powint(phi, 2 * N + 1)*powint(dphi, 2 * M - 1) / powint(2., M - 1)
		- 2.*M*M*ALPHA*powint(phi, 2 * N + 1)*powint(dphi, 2 * M - 2)*ddphi / powint(2., M - 1) - M*M*ALPHA*ALPHA*powint(phi, 4 * N + 2)*powint(dphi, 4 * M) / powint(2., 2 * M - 1) + X2*(BETA*powint(dphi, 2 * (KINL - 1)) / powint(2., KINL - 1)))
		*dphi*dphi / 2. / (H + ALPHA*M*powint(phi, 2 * N + 1)*powint(dphi, 2 * M + 1) / powint(2., M)) / (H + ALPHA*M*powint(phi, 2 * N + 1)*powint(dphi, 2 * M + 1) / powint(2., M));
}

double Gs(double phi, double dphi, double ddphi, double H) {
	return (1 + (2 * N + 1)*(M + 1)*ALPHA*powint(phi, 2 * N)*powint(dphi, 2 * M) / powint(2., M - 1) - 6.*M*M*ALPHA*H*powint(phi, 2 * N + 1)*powint(dphi, 2 * M - 1) / powint(2., M - 1)
		+ 3.*M*M*ALPHA*ALPHA*powint(phi, 4 * N + 2)*powint(dphi, 4 * M) / powint(2., 2 * M - 1) + X2*(BETA*(2 * KINL - 1)*powint(dphi, 2 * (KINL - 1)) / powint(2., KINL - 1)))
		*dphi*dphi / 2. / (H + ALPHA*M*powint(phi, 2 * N + 1)*powint(dphi, 2 * M + 1) / powint(2., M)) / (H + ALPHA*M*powint(phi, 2 * N + 1)*powint(dphi, 2 * M + 1) / powint(2., M));
}

double cs2(double phi, double dphi, double ddphi, double H) {
	return ((1 - (2 * N + 1)*(M - 1)*ALPHA*powint(phi, 2 * N)*powint(dphi, 2 * M) / powint(2., M - 1) - 4.*M*ALPHA*H*powint(phi, 2 * N + 1)*powint(dphi, 2 * M - 1) / powint(2., M - 1)
		- 2.*M*M*ALPHA*powint(phi, 2 * N + 1)*powint(dphi, 2 * M - 2)*ddphi / powint(2., M - 1) - M*M*ALPHA*ALPHA*powint(phi, 4 * N + 2)*powint(dphi, 4 * M) / powint(2., 2 * M - 1) + X2*(BETA*powint(dphi, 2 * (KINL - 1)) / powint(2., KINL - 1)))
		/ (1 + (2 * N + 1)*(M + 1)*ALPHA*powint(phi, 2 * N)*powint(dphi, 2 * M) / powint(2., M - 1) - 6.*M*M*ALPHA*H*powint(phi, 2 * N + 1)*powint(dphi, 2 * M - 1) / powint(2., M - 1)
			+ 3.*M*M*ALPHA*ALPHA*powint(phi, 4 * N + 2)*powint(dphi, 4 * M) / powint(2., 2 * M - 1) + X2*(BETA*(2 * KINL - 1)*powint(dphi, 2 * (KINL - 1)) / powint(2., KINL - 1))));
}

double maxval(double x, double y, double z, double w) {
	double hoge = 0.0;

	hoge = ((0.0 >= x) ? 0.0 : x);
	hoge = ((hoge >= y) ? hoge : y);
	hoge = ((hoge >= z) ? hoge : z);
	hoge = ((hoge >= w) ? hoge : w);

	return hoge;
}

double powint(double phi, int i) {
	double hoge = 1.;
	int j;

	for (j = 0; j<i; j++) {
		hoge *= phi;
	}

	return hoge;
}