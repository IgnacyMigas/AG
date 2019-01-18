#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <ga/ga.h>
#include <ctime>
#include <algorithm>

#define INSTANTIATE_REAL_GENOME
#include <ga/GARealGenome.h>

using namespace std;

int popsize  = 100;
int ngen     = 1500;
float pmut   = 1.0;
float pcross = 0.0;

const int MAX_NR = 1003;
float v_min=51, v_max=-1, p_max=1;
float prenty[MAX_NR];
int len=0;

float odchylenie_max=0;

float objective(GAGenome &);

float pole(float a, float b, float c);
void save_output(const GAGenome& g);
void licz(GAGenome & g, int& cnt, float& odchylenie);

float Objective(GAGenome&);
void  Initializer(GAGenome&);
int   Mutator(GAGenome&, float);
float Comparator(const GAGenome&, const GAGenome&); 

int main(int argc, char *argv[])
{
	float v;
	float m1=0, m2=0, m3=0;
	unsigned seed = time(0);
	
	ifstream we("prety.txt");
	
	while(!we.eof()){
		we>>v;
		we>>v;
		prenty[len++] = v;
		if (v>m1) {
			m3=m2;
			m2=m1;
			m1=v;
		}
		else if (v>m2){
			m3=m2;
			m2=v;
		}
		else if (v>m3){
			m3=v;
		}
		
		if (v<v_min) v_min = v;
		if (v>v_max) v_max = v;
	}
	
	
	we.close();

	
	p_max = pole(m1, m2, m3)/2.;

  GAListGenome<int> genome(Objective);
  genome.initializer(::Initializer);
  genome.mutator(::Mutator);
  genome.comparator(::Comparator);

  GASteadyStateGA ga(genome);
  ga.minimize();
  ga.pReplacement(1.0);
  ga.populationSize(popsize);
  ga.nGenerations(ngen);
  ga.pMutation(pmut);
  ga.pCrossover(pcross);
  ga.selectScores(GAStatistics::AllScores);
  ga.parameters(argc, argv);
  ga.initialize(seed);
  while(!ga.done()) {
    ga.step();
  }

  genome = ga.statistics().bestIndividual();
  float od = 0;
  int cnt = 0;
  licz(genome, cnt, od);
  cout<< "\n\nOdchylenie=" << od << ", ilosc=" << cnt << endl<<endl;
  save_output(genome);

  return 0;
}

float Objective(GAGenome & g)
{
	float odchylenie = 0;
	int cnt = 0;
	licz(g, cnt, odchylenie);
  
  float fitness = (len/3.-cnt)+(odchylenie/p_max*(len*1.5));
  if (fitness != fitness)
	  return len*100;
  return fitness;
}

void Initializer(GAGenome& g) {
  GAListGenome<int> &child=(GAListGenome<int> &)g;
  while(child.head()) child.destroy(); // destroy any pre-existing list
  
  child.insert(0,GAListBASE::HEAD);
  for (int i=1; i<len; ++i){
	  child.insert(i);
  }
}

int Mutator(GAGenome& g, float pmut) {
  GAListGenome<int> &child=(GAListGenome<int> &)g;
  register int n, i;
  if ((GARandomFloat() >= pmut) || (pmut <= 0)) return 0;

  n = child.size();
  
  if (GARandomFloat()<0.5) {
    child.swap(GARandomInt(0,n-1),GARandomInt(0,n-1)); // swap only one time
  }
  else {
    int nNodes = GARandomInt(1,((int)(n/2-1)));       // displace nNodes 
    child.warp(GARandomInt(0,n-1));                   // with or without
    GAList<int> TmpList;                              // inversion
    for(i=0;i<nNodes;i++) {
      int *iptr = child.remove();
      TmpList.insert(*iptr,GAListBASE::AFTER);
      delete iptr;
      child.next();
    }
    int invert;
    child.warp(GARandomInt(0,n-nNodes));
    invert = (GARandomFloat()<0.5) ? 0 : 1;
    if (invert) TmpList.head(); else TmpList.tail();

    for(i=0;i<nNodes;i++) {
      int *iptr = TmpList.remove();
      child.insert(*iptr,GAListBASE::AFTER);
      delete iptr;
      if (invert) TmpList.prev(); else TmpList.next();
    }
  }
  child.head();		// set iterator to root node

  return (1);
}

float Comparator(const GAGenome& g1, const GAGenome& g2) {
  GAListGenome<int> &a = (GAListGenome<int> &)g1;
  GAListGenome<int> &b = (GAListGenome<int> &)g2;

  int sum=0;
  for (int i=0; i<len; ++i){
	  if (*a.current() != *b.current())
		  sum++;
  }
  return sum;
}

template <> int GAListGenome<int>::write(ostream & os) const
{
  int *cur, *head;
  GAListIter<int> tmpiter(*this);
  if((head=tmpiter.head()) != 0) {
    os << *head+1 << " " << prenty[*head]<<endl;
    for(cur=tmpiter.next(); cur && cur != head; cur=tmpiter.next())
      os << *cur+1 << " " << prenty[*cur]<<endl;
  }

  return os.fail() ? 1 : 0;
}

float pole(float a, float b, float c){
	float p = (a+b+c)/2.0;
	return sqrt(p*(p-a)*(p-b)*(p-c));
}

void licz(GAGenome & g, int& cnt, float& odchylenie){
  GAListGenome<int> & genome = (GAListGenome<int> &)g;
  float sum=0;
  float tab[MAX_NR/3];
  if(genome.head()) {
    for(int i=0; i<len/3; i++) {
      float a = prenty[*genome.current()];
      float b = prenty[*genome.next()];
      float c = prenty[*genome.next()];
	  genome.next();
      if (a+b>c && a+c>b && b+c>a){
		  tab[cnt]=pole(a, b, c);
		  sum+=tab[cnt];
		  cnt++;
	  }
    }
	float avg = sum/cnt;
	sum = 0;
	for(int i=0; i<cnt; i++) {
		sum += pow((tab[i]-avg), 2);
	}
	odchylenie=sqrt(sum/cnt);
	if (odchylenie>odchylenie_max) {
		odchylenie_max=odchylenie;
	}
  }
}

void save_output(const GAGenome& g) {
	GAListGenome<int> & genome = (GAListGenome<int> &)g;

	ofstream wy("output.txt");
	
	wy << *genome.head()+1 << endl;
	for (int i=1; i<len; ++i){
		wy << *genome.next()+1 << endl;
	}
	
	wy.close();
}


