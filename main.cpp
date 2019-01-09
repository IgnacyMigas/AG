#include <cstdlib>
#include <iostream>
#include <math.h>
#include <map>
#include <ga/ga.h>
#include <ctime>

#define INSTANTIATE_REAL_GENOME
#include <ga/GARealGenome.h>

using namespace std;

int popsize  = 100;
int ngen     = 1000;
float pmut   = 0.01;
float pcross = 0.9;

float objective(GAGenome &);

int main(int argc, char *argv[])
{
  GABin2DecPhenotype map;
  map.add(30, 0, 100000*M_PI);

  GABin2DecGenome genome(map, objective);

  GASimpleGA ga(genome);
  ga.populationSize(popsize);
  ga.pMutation(pmut);
//  ga.set("mutation_probability",pmut);
  ga.pCrossover(pcross);
  
//  GASigmaTruncationScaling scaling;
//  ga.scaling(scaling);

//  GARankSelector selector;
//  GARouletteWheelSelector selector;
//  GATournamentSelector selector;
//  GADSSelector selector;
//  GASRSSelector selector;
//  GAUniformSelector selector;
//  ga.selector(selector);  
  
  ga.elitist(gaFalse);
  
  ga.nGenerations(ngen);                            ga.terminator(GAGeneticAlgorithm::TerminateUponGeneration);
//  ga.pConvergence(0.7);  ga.nConvergence(10);       ga.terminator(GAGeneticAlgorithm::TerminateUponConvergence);
//  ga.pConvergence(0.1);                             ga.terminator(GAGeneticAlgorithm::TerminateUponPopConvergence);    
  ga.nBestGenomes(3);
  
  ga.scoreFilename("zbieznosc.dat");
  ga.scoreFrequency(10);
  ga.flushFrequency(50);
  ga.selectScores(GAStatistics::Mean | 
                  GAStatistics::Maximum | 
                  GAStatistics::Minimum | 
                  GAStatistics::Deviation | 
                  GAStatistics::Diversity);
  ga.recordDiversity(gaTrue);                
  ga.initialize((unsigned)time(0));
  while (!ga.done())
  {
 //   ++ga;
    ga.step();
    cout << ga.generation() << "\t conv=" << ga.convergence();
    genome = ga.statistics().bestIndividual();
    cout << "\t x=" << genome.phenotype(0) << "\t Fbest=" << objective(genome) << endl;
  }
  genome = ga.statistics().bestIndividual();
  cout << "Najlepsze rozwiazanie to F(x=";
  cout << genome.phenotype(0) << ")=";
  cout << objective(genome) << endl;

  getch();
  return EXIT_SUCCESS;
}

float objective(GAGenome & c)
{
  GABin2DecGenome & genome = (GABin2DecGenome &)c;

  float w,x;
  x=genome.phenotype(0);
  w=fabs(sin(M_PI*x)*exp(-x/(6*M_PI)));
  return w;
}

