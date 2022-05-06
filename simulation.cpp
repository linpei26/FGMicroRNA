#include <iostream>
#include <vector>
#include <chrono>
#include <stdlib.h>
#include <stack>
#include <array>
#include <queue>
#include <time.h>
#include <random>
#include <math.h>
#include <fstream>
using namespace std;

int main (int argc, char* argv[])
{
	
	int ngene = 100;
	int total_optim = 0;
	// 	
	int N = 10000;
	int totalrun = 100;
	double omega;
	double index;
	double strongeffectprop;
	double envchange_higherbound;
	double envchange_lowerbound;
	double m;
	double largemu;
	double smallmu;
	
							omega = atof(argv[1]);
							index = atof(argv[2]);
				 strongeffectprop = atof(argv[3]);
							    m = atof(argv[4]);
						  smallmu = atof(argv[5]);
						  largemu = atof(argv[6]);
		
	double m_in = 1/m;
	double mu1 = 1/smallmu;
	double mu2 = 1/largemu;
	double mind = 0.5 * smallmu;			
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	
	std::exponential_distribution<double> distribution_distance(m_in);
	std::exponential_distribution<double> distribution_effect_small(mu1);
	std::exponential_distribution<double> distribution_effect_large(mu2);	
	
	std::uniform_real_distribution<double> direction(0.0,1.0);
	std::uniform_real_distribution<double> category(0.0,1.0);
	std::uniform_real_distribution<double> distribution_prob_to_fix(0.0,1.0);
	std::uniform_real_distribution<double> distribution_prob_to_dice(0.0,1.0);
	
	std::poisson_distribution<int> num_of_optimum(1);
	
	int nresults = 1;
	
	int nfix = 0;
	int print = 0;
	
	int gene = 1;
	
	while(gene <= 10000) {
		int the_last_shoot_1 = 0;
		int the_last_shoot_2 = 0;
		double poss1 = 0;int posn1 = 0;double pospr1 = 0;int summ1 = 0;
		double poss2 = 0;int posn2 = 0;double pospr2 = 0;int summ2 = 0;
		
		// initialize, wild-type at the origin!!!
		double trait_0 = 0;
		double trait_1 = 0; // a temporary variable, for the mutant.
		double trait_optim;
		// for each gene;
		int Hi1 = 0;
		int Hi2 = 0;
		int Mi1 = 0;
		int Mi2 = 0;
		double effect1_sum = 0;
		double effect2_sum = 0;
				
		int nw = 0;
		int ns = 0;
		int n_fix = 0;

		int n_of_optim = num_of_optimum(generator);
		int nmut = 1;
		total_optim = total_optim + n_of_optim;
		
		int f10 = 0;
		int f11 = 0;
		int f20 = 0;
		int f21 = 0;
		
		if(n_of_optim > 0) {
			for(int noo = 1;noo <= n_of_optim;noo++) {
				
				int nfixed = 0;
				double distance;				
				distance = distribution_distance(generator);
					
				// the novel optimum is too close, relative to a small step
				while(distance < mind) {
					distance = distribution_distance(generator);
				}
				
				double direction_of_new_optim;
				direction_of_new_optim = direction(generator);				
				if(direction_of_new_optim > 0.5) {
					direction_of_new_optim = 1;
				} else {
					direction_of_new_optim = -1;
				}				
				distance = distance * direction_of_new_optim;				

				trait_optim = trait_0 + distance;
		
				double z = distance;
				// Lande 1976: wz = exp((-z^2)/(2*omega^2))
				// fitness for the wild-type
				double wz_0_optim = (-1*pow(z,2))/(2*pow(omega,index));
				wz_0_optim = exp(wz_0_optim);
				int n_mutation = 0;				
				int stop = 0;
				
				// this line works for the while loop.
				double abs_dist = abs(trait_0 - trait_optim);
				int the_last_shoot;
				while(abs_dist > mind) {
			
					double cate;
					double dire;
					cate = category(generator);
					dire = direction(generator);					
					double effect;					
					int effectsize;
								
					if(cate < strongeffectprop) {
						effectsize = 2;
						effect = distribution_effect_large(generator);					
						if(dire < 0.5) {
							effect = -1 * effect;
						}
					} else {
						effectsize = 1;
						effect = distribution_effect_small(generator);
						if(dire < 0.5) {
							effect = -1 * effect;
						}
					}
					
					trait_1 = trait_0 + effect;
					double z1;
					z1 = trait_optim - trait_1;
					double wz_1_optim = (-1*pow(z1,2))/(2*pow(omega,index));
					wz_1_optim = exp(wz_1_optim);					
					double selectioncoefficient = 1 - (wz_0_optim/wz_1_optim);
					double s = selectioncoefficient;
					double fixprob = (1 - exp(-2*s))/(1 - exp(-2*N*s));
					
					if(isnan(fixprob)) {
						fixprob = 0;
					}
					
					double dice;
					dice = distribution_prob_to_fix(generator);		
					int fix_or_not;
					if(dice < fixprob) {
						fix_or_not = 1;					
					} else {			
						fix_or_not = 0;						
					}								
					
					int rightmove = 0;
						if(trait_0 < trait_1 && trait_1 < trait_optim) { rightmove = 1; }
						if(trait_optim < trait_1 && trait_1 < trait_0) { rightmove = 1; }					
					/* 		
						if(effectsize == 1) {
							summ1 = summ1 + 1;
							if(s > 0) {
								posn1 = posn1 + 1;
								poss1 = poss1 + s;
								pospr1 = pospr1 + fixprob;
							}					
						}
						if(effectsize == 2) {
							summ2 = summ2 + 1;
							if(s > 0) {
								posn2 = posn2 + 1;
								poss2 = poss2 + s;
								pospr2 = pospr2 + fixprob;
							}					
						}
					 */
					if(fix_or_not == 1) {
						if(effectsize == 1 && rightmove == 1) {f11 = f11 + 1;}
						if(effectsize == 1 && rightmove == 0) {f10 = f10 + 1;}
						if(effectsize == 2 && rightmove == 1) {f21 = f21 + 1;}
						if(effectsize == 2 && rightmove == 0) {f20 = f20 + 1;}
					}
					
					if(fix_or_not == 1) {												
						// update the trait.
						n_fix = n_fix + 1;
						trait_0 = trait_1;						
						double deviation_1_raw = trait_0 - trait_optim;
						double deviation_1 = abs(deviation_1_raw);
																		
						// update the wild-type fitness 
						// wz = exp((-z^2)/(2*omega^2))			
						wz_0_optim = (-1*pow(deviation_1_raw,2))/(2*pow(omega,index));
						wz_0_optim = exp(wz_0_optim);
						
						the_last_shoot = effectsize;
						abs_dist = abs(trait_0 - trait_optim);	
						
						if(effectsize == 2) {
							effect2_sum = effect2_sum + effect;
							if(effect < 0) {Hi2 = Hi2 + 1;}
							if(effect > 0) {Mi2 = Mi2 + 1;}						
						} else {
							effect1_sum = effect1_sum + effect;
							if(effect < 0) {Hi1 = Hi1 + 1;}
							if(effect > 0) {Mi1 = Mi1 + 1;}							
						}						
					}
				}
				
				if(the_last_shoot == 1) {the_last_shoot_1 = the_last_shoot_1 + 1;}
				if(the_last_shoot == 2) {the_last_shoot_2 = the_last_shoot_2 + 1;}
			}

			if(n_fix > 0) {
				std::cout << gene << ":" << n_of_optim << "\t";
				std::cout << trait_optim << "\t";
				std::cout << trait_0 << "\t";
				std::cout << Hi1 << "\t" << Mi1 << "\t" << effect1_sum << "\t";
				std::cout << Hi2 << "\t" << Mi2 << "\t" << effect2_sum << "\t";				
				/* std::cout << summ1 << "\t";
				std::cout << summ2 << "\t";				
				std::cout << posn1 << "\t";
				std::cout << posn2 << "\t";				
				std::cout << poss1 << "\t";
				std::cout << poss2 << "\t";				
				std::cout << pospr1 << "\t";
				std::cout << pospr2 << "\t"; */
				// std::cout << the_last_shoot_1 << "\t";
				// std::cout << the_last_shoot_2 << "\t";
				std::cout << f11 << "\t";
				std::cout << f10 << "\t";
				std::cout << f21 << "\t";
				std::cout << f20 << "\n";
				gene = gene + 1;
			}
		}
	}

}
//
// Feb-04-2021
////////////////
////////////////
////////////////
////////////////
