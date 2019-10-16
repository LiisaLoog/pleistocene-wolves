//-----------------------------------------------------------------------------------
// (c) Anders Eriksson 
//
// Expansion model
//
// Compile using:
// g++ -Wall -O3 -o wolf_expansion wolf_expansion.cpp ../Data/MersenneTwist.cpp -I ../Data
//
//-----------------------------------------------------------------------------------

#include <assert.h>

#include "Util.h"
#include "MersenneTwist.h"

//-----------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------

#include "../Data/neigh_inc.cpp"

#include "../Data/sample_inc.cpp"

//-----------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------

int lines[nDemes][nInds];
int nLines[nDemes];

const int nNodes  = 2*nInds-1;
struct Node {
	int left, right, deme;
	double t;
};

Node nodes[nNodes];

//-----------------------------------------------------------------------------------
// main
//-----------------------------------------------------------------------------------

int main(int ac, char *av[])
{
	//----------------------------------------------------------
	// Read parameters from command line
	//----------------------------------------------------------

	// first print general info
	printf("%% "); for (int i = 0; i < ac; i++) printf("%s ", av[i]); printf("\n\n");
	printf("%% seed = %d;\n", MT::SeedTime(RPC(seed,0)));

	Parameter(const int, nSamples, 2000);
	Parameter(const double, R2_rep_T, -100);

	Parameter(const double, m1, 0.01);
	Parameter(const double, K1, 0.01);
	Parameter(const double, xBottle1, 0);
	Parameter(const double, tStart1, 0);
	Parameter(const double, dT1, 0);
	
	Parameter(const double, m2, 100);
	Parameter(const double, K2, 100);
	Parameter(const double, xBottle2, 1);
	Parameter(const double, tStart2, 200);
	Parameter(const double, dT2, 100);
	
	Parameter(const bool, printT, false);
	
	Parameter(const int, firstDeme, 1); assert(1 <= firstDeme && firstDeme <= nDemes);
	
	const double inf = 1.0/0.0;
	
	//---------------------------------------------------------------------------
	// Init mc parameters
	//---------------------------------------------------------------------------
	
	const int nParam = 5;
	//                              0        1        2          3        4
	const double paramL[nParam] = { log(m1), log(K1), xBottle1, tStart1, dT1};
	const double paramH[nParam] = { log(m2), log(K2), xBottle2, tStart2, dT2};
	double param[nParam];
	
	printf("%% sample[1], R2[2], m[3], K[4], xBottle[5], tStart[6], dT[7], T[5:end]\n");
	fflush(stdout);
	
	//---------------------------------------------------------------------------
	// do uniform sweep
	//---------------------------------------------------------------------------

	double T[nInds][nInds];
	int nD[nNodes], Rpos[nNodes], ind_order[nInds];
	
	
	for (int sample = 1; sample <= nSamples; sample++) {

		//----------------------------------------------------------
		// Init spatial list
		//----------------------------------------------------------
		
		for (int i = 0; i < nParam; i++) {
			do {
				param[i] = paramL[i] + (paramH[i] - paramL[i])* MT::Uni();;
			} while (param[i] < paramL[i]  || param[i] > paramH[i]);
		}
		
		const double m         = exp(param[0]);      // 1
		const double K         = exp(param[1]);      // 2
		const double xBottle   = param[2];      // 3
		const double tStart    = param[3];      // 4
		const double dT        = param[4];      // 5
		
		const double D_bottle = -log(std::max(1e-20,1.0-xBottle));
		
		//----------------------------------------------------------
		// Init spatial list
		//----------------------------------------------------------
		
		int    colonDeme[nDemes];
		double colonTime[nDemes];
		
		int colon_order[nDemes];
		for (int i = 0; i < nDemes; i++) colon_order[i] = -1;
		colon_order[firstDeme-1] = 0;
		
		int n_colonised = 0;
		for (int i_colon = 0; n_colonised < nDemes; i_colon++) {
			const double t = tStart - i_colon*dT;
			if (t <= 0) break;
			
			for (int i = 0; i < nDemes; i++) {
				if (colon_order[i] >= 0 && colon_order[i] <= i_colon) {
					for (int j = neighStart[i]; j < neighStart[i+1]; j++) {
						if (colon_order[neigh[j]] < 0) {
							colonDeme[n_colonised] = neigh[j];
							colonTime[n_colonised] = t;
							colon_order[neigh[j]] = i_colon+1;
							n_colonised++;
						}
					}
				}
			}
		}
		n_colonised--;
		
		//----------------------------------------------------------
		// Init
		//----------------------------------------------------------

		for (int i = 0; i < nInds; i++) {
			nodes[i].left = nodes[i].right = -1;
			nodes[i].deme = ind_deme[i];
			nodes[i].t    = ind_t[i];
		}
		int i_node = nInds;

		for (int i = 0; i < nDemes; i++) nLines[i] = 0;
		int nLinesTot = 0;
		
		// assume individuals waiting are ordered with in increasing t
		double t = nodes[0].t;
		if (n_colonised >= 0) assert(t < colonTime[n_colonised]);

		int i_ind_waiting = 0;
		while (i_ind_waiting < nInds && nodes[i_ind_waiting].t == t) {
			// fprintf(stderr,"add %d %.12lg\n", i_ind_waiting, t);
			const int i = nodes[i_ind_waiting].deme;
			lines[i][nLines[i]] = i_ind_waiting;
			nLines[i]++;
			i_ind_waiting++;
			nLinesTot++;
		}

		double rate[nDemes];
		
		const double x = 1.0/K;
		
		while (i_node < nNodes) {

			double rate_tot = 0;
			for (int i = 0; i < nDemes; i++) {
				rate_tot += rate[i] = (x*0.5)*(nLines[i]*(nLines[i]-1)) + m*(nNeighs[i]*nLines[i]);
			}
			assert(rate_tot > 0);

			t += MT::Exp()/rate_tot;
			const double t_waiting = (i_ind_waiting < nInds) ? nodes[i_ind_waiting].t : inf;
			const double t_colon   = (n_colonised >= 0) ? colonTime[n_colonised] : inf;
			if (t > t_waiting || t > t_colon) {
				if (t_waiting < t_colon) {
					t = nodes[i_ind_waiting].t;
					while (i_ind_waiting < nInds && nodes[i_ind_waiting].t == t) {
						// fprintf(stderr,"add %d %.12lg\n", i_ind_waiting, t);
						const int i = nodes[i_ind_waiting].deme;
						lines[i][nLines[i]] = i_ind_waiting;
						nLines[i]++;
						i_ind_waiting++;
						nLinesTot++;
					}
				} else {
					t = colonTime[n_colonised];
					
					for (; n_colonised >= 0 && colonTime[n_colonised] == t; --n_colonised) {
						const int i_d = colonDeme[n_colonised];
						
						// do bottleneck
						for (double tt = D_bottle; nLines[i_d] > 1;) {
							
							tt -= MT::Exp()/(0.5*(nLines[i_d]*(nLines[i_d]-1)));
							if (tt < 0) break;
							
							// coalescent event
							nodes[i_node].t = t;
							nodes[i_node].deme = i_d;
							
							// assert(nLines[i_d]>1);
							int i = MT::Rand() % nLines[i_d];
							nodes[i_node].left = lines[i_d][i]; assert(lines[i_d][i] < i_node);
							--nLines[i_d]; assert(nLines[i_d]>0);
							lines[i_d][i] = lines[i_d][nLines[i_d]];
							
							i = MT::Rand() % nLines[i_d];
							nodes[i_node].right = lines[i_d][i];
							lines[i_d][i] = i_node;
							
							i_node++;
							--nLinesTot;
						}
						
						// move remaining lines to random neighbours
						while (nLines[i_d]>0) {
							--nLines[i_d];
							for (;;) {
								const int k = neigh[neighStart[i_d] + (MT::Rand() % nNeighs[i_d])];
								
								if (colon_order[k] < colon_order[i_d]) {
									lines[k][nLines[k]] = lines[i_d][nLines[i_d]];
									nLines[k]++;
									break;
								}
							}
						}
					}
				}
				
				continue;
			}
			
			double xi = MT::Uni()*rate_tot;
			for (int i_d = 0; i_d < nDemes; i_d++) {
				if (xi < rate[i_d]) {
					const double mm = m*nNeighs[i_d];
					
					if (xi < mm*nLines[i_d]) {
						const int k = (int)(xi/mm); // k = id of line
						
						// generate randomly chose active edge
						const int j_d = neigh[neighStart[i_d] + (MT::Rand() % nNeighs[i_d])];
						// move the chosen line from deme i_d to deme j_d
						lines[j_d][nLines[j_d]] = lines[i_d][k]; nLines[j_d]++;
						--nLines[i_d]; lines[i_d][k] = lines[i_d][nLines[i_d]];
						
					} else {

						// coalescent event
						nodes[i_node].t = t;
						nodes[i_node].deme = i_d;

						int i = MT::Rand() % nLines[i_d];
						nodes[i_node].left = lines[i_d][i]; assert(lines[i_d][i] < i_node);
						--nLines[i_d]; assert(nLines[i_d]>0);
						lines[i_d][i] = lines[i_d][nLines[i_d]];

						i = MT::Rand() % nLines[i_d];
						nodes[i_node].right = lines[i_d][i];
						lines[i_d][i] = i_node;
						
						i_node++;
						--nLinesTot;
					}
					break;
				}
				
				xi -= rate[i_d];
			} // end for i_d

		} // end while
		
		//----------------------------------------------------------
		// calc T matrix
		//----------------------------------------------------------
		
		for (int i = 0; i < nInds; i++) nD[i] = 1;
		for (int i = nInds; i < nNodes; i++) nD[i] = nD[nodes[i].left] + nD[nodes[i].right];
		
		Rpos[nNodes-1] = nD[nNodes-1];
		for (int i = nNodes-1; i >= nInds; --i) {
			Rpos[nodes[i].right] = Rpos[i];
			Rpos[nodes[i].left ] = Rpos[i] - nD[nodes[i].right];
		}
		
		for (int i = 0; i < nInds; i++) ind_order[Rpos[i]-1] = i;
		
		for (int k = nInds; k < nNodes; k++) {
			const int eL = Rpos[nodes[k].left], eR = Rpos[nodes[k].right];
			for (int i = eL - nD[nodes[k].left]; i < eL; ++i) {
				const int iL = ind_order[i];
				for (int j = eR - nD[nodes[k].right]; j < eR; ++j) {
					const int iR = ind_order[j];
					T[iL][iR] = T[iR][iL] = nodes[k].t;
				}
			}
		}
		
		for (int i = 0; i < nInds; i++) T[i][i] = nodes[i].t;
		
		//---------------------------------------------------------------------------
		// measure SSR
		//---------------------------------------------------------------------------
		
		double R2 = 0;
		
		double nSamples_full = 0;
		for (int i = 0; i < nInds; i++) {
			for (int j = 0; j <= i; j++) {
				if (usePop[i] == 1 && usePop[j] == 1) {
					const double dT = T[i][j] - data_T[i][j];
					R2 += dT*dT;
					nSamples_full++;
				}
			} // end loop j
		} // end loop i
		
		R2 = 1.0 - R2/(nSamples_full * data_T_var);
		
		//---------------------------------------------------------------------------
		// measure sumstat
		//---------------------------------------------------------------------------
		
		double T_sumstat[nSumStats], R2_sumstat = 0;
		for (int k = 0; k < nSumStats; k++) {
			const int * const p1 = idSup[idSumStats[k][0]];
			const int * const p2 = idSup[idSumStats[k][1]];
			const int n1 = nIdSup[idSumStats[k][0]], n2 = nIdSup[idSumStats[k][1]];
			
			double ts = 0;
			for (int i = 0; i < n1; i++) {
				for (int j = 0; j < n2; j++) ts += T[p1[i]][p2[j]];
			}
			ts /= n1*n2;
			T_sumstat[k] = ts;
			
			const double dT = ts - data_T_sumstat[k];
			R2_sumstat += dT*dT;
		}
		
		R2_sumstat = 1.0 - R2_sumstat/(nSumStats * data_T_var_sumstat);
		
		//---------------------------------------------------------------------------
		// print results
		//---------------------------------------------------------------------------

		if (R2_rep_T <= -100 || R2_sumstat > R2_rep_T)  {
			printf("%d\t%.6lg\t%.6lg\t%.6lg\t%.6lg\t%.6lg\t%.6lg\t%.6lg", sample, R2_sumstat, R2, m, K, xBottle, tStart, dT);
			for (int i = 0; i < nSumStats; i++) printf("\t%.4lf",T_sumstat[i]);
			if (printT) {
				for (int i = 0; i < nInds; i++) {
					for (int j = 0; j <= i; j++) printf("\t%.3lf", T[i][j]);
				}
			}
			printf("\n");
			fflush(stdout);
		}
	}

	return 0;
}
