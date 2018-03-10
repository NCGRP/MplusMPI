#include <map>
#include <unordered_map>
#include <chrono>

//#include <future>
//#include <mutex>

/*
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <math.h>
#include <numeric>
#include <set>
#include <string>
#include <string.h>
#include <sstream>
#include <time.h>
#include <unistd.h>
#include <vector>
*/


using namespace std;

//code in rare.hpp is derived from https://github.com/hildebra/Rarefaction
//Paul Saary, Kristoffer Forslund, Peer Bork, Falk Hildebrand; RTK: efficient rarefaction 
//analysis of large datasets, Bioinformatics, Volume 33, Issue 16, 15 August 2017, Pages 
//2594–2595, https://doi.org/10.1093/bioinformatics/btx206

/***************TYPEDEFS*****************/
typedef double mat_fl;
typedef unsigned int uint;
typedef unordered_map <uint, uint> rare_map;
typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
               // e.g. keep one global instance (per thread)
typedef std::map<std::string, vector<string>> LvlUp;
typedef std::map<std::string, int> GeneIDidx;

/***************TEMPLATES*****************/
template<typename T> T getMedian(vector<T>& in){
	sort(in.begin(), in.end());
	size_t size = in.size();
	if (size == 0){ return (T)0; }
	if (size == 1){ return (in[0]) ; }
	if (size == 2){ return (in[0] + in[1]) / 2; }
	T median(in[size / 2]);
	if (size % 2 == 0)	{
		median = (in[size / 2 - 1] + in[size / 2]) / 2;
	}
	return median;
}

/***************STRUCT*****************/
struct options
{
public:
	//options(int argc, char** argv);
	//options(std::string, std::string , int repeats, std::vector<double> depth, int NoOfMatrices, bool verbose, unsigned int threads);
	options(void);
	void print_rare_details();
	//~options();

	//vars
  std::string input = "";
  std::string output = "";
  std::string mode  = "";
  std::string referenceDir = "";
  std::string referenceFile = "";
  std::string map = "";
  std::vector<double> depth;
  long depthMin;
  unsigned int repeats;
  unsigned int write = 0;
  unsigned int threads = 1;
  bool writeSwap = true;
  bool verbose = false;
  bool oldMapStyle = true;

    std::string modDB;
    int modRedund;
    float modEnzCompl;
    float modModCompl;
    bool modWrXtraInfo;
    bool modCollapse;
    bool calcCoverage;

  std::string modDescr;
  std::string modHiera;
  std::string xtra;
};


/***************CLASSES*****************/
class DivEsts{
public:
	DivEsts():richness(0),shannon(0),
		simpson(0),invsimpson(0),chao1(0),eve(0){}
	~DivEsts(){}
	//void print2file(const string);
	//data vectors
	vector<vector<long> > richness;
	vector<vector<double> > shannon,simpson,invsimpson,chao1,eve;
	vector<vector<vector<uint> > > cntsx; //contains the rarefied sampled frequency of each row
	string SampleName;
	int depth;
};

class smplVec{
public:
	//smplVec(const string, const int);
	smplVec(const vector<mat_fl>&, const int);
	~smplVec(){
		//delete[] arr;
	}
	void rarefy(vector<double> ,string o,int rep,DivEsts*, vector<vector<rare_map>>& RareSample,
		vector<string>& retCntsSampleName, string& skippedSample, vector<vector<vector<uint>>>* ,vector<vector<vector<uint>>>* , int=0,bool=false, bool=false);
	//long getRichness(rare_map& cnts);
	long getRichness(const vector<unsigned int>&);
	//int maxSiz(){return vector<unsigned short>::max_size();}
	vector < string > getRowNames(){ return(IDs); }

private:
	int binarySearch(vector<float>,const float x);
	//void shuffle();
	void shuffle_singl();

	//diversity indices
	//method: 1=shannon, 2=simpson, 3=invsimpson
	vector<double> calc_div(const vector<uint>& vec,int meth=1, float base=2.718282f);
	vector <double> calc_div(rare_map& , int meth=1, float base=2.718282f);
	double calc_chao1(const vector<uint> & vec,int corrBias=1);
	double calc_chao1(rare_map& , int corrBias=1); //corrBias: 0/1
	double calc_eveness(const vector<uint>& vec);
	double calc_eveness(rare_map& );

	void print2File(const vector<unsigned int>&,const string);
	//unsigned short * arr;
	vector<string> IDs;
	vector<unsigned int> arr;
	double totSum;
	vector<MyRNG> rng_P;
	MyRNG rng;
	int num_threads;
	long richness;
	double Shannon;
	int numFeatures;

	//vector<float> vec;
};

class column{
	public:
		double colsum;
		string id;

};

/*
class HMat
{
public:
	HMat(string L, vector<string> Samples, vector<string> Features);
	~HMat(){}
	//get
	//unsigned long operator [](int i) const    { return registers[i]; }
	//set
	void set(string kk, int j, mat_fl v);

	void print(ofstream&);

private:
	GeneIDidx Feat2mat;
	string LvlName;

	vector<string> FeatureNs, SampleNs;
	vector< mat_fl > empty;
	vector< vector< mat_fl > > mat;
};
*/

class Matrix
{
// convention: mat[smpl][feature]
public:
	Matrix(stringstream& in);
	//Matrix(const string inF);
	//read and write
	//Matrix(const string inF, const string, const string xtra, vector<string>& outFName, bool highLvl = false, bool NumericRowId = false, bool writeTmpFiles = true);
	//read to mem
	//Matrix(const string inF, const string xtra, bool highLvl = false); 
	//module abundance matrix
	//Matrix(const vector<string>& rnms, const vector<string>& cnms);
	//empty opbject
	//Matrix(void);
	//normalize on the fly on vector colSums
	//Matrix(const string inF, const string outF, vector< double> colsums, vector<string>colNmds);
	~Matrix(void);
	void addTtlSmpl(vector<mat_fl> x, int idx) { mat[idx] = x; }
	void splitOnHDD(string out_seed);
	void writeSums(string);
	void normalize();
	void transpose();
	void writeMatrix(const string ofile,bool onlyFilled=false);
	size_t smplNum(){ return colIDs.size(); }
	int rowNum(){ return rowIDs.size(); }

	smplVec* getSampleVec(uint which){ return new smplVec(mat[which],1); }
	string getSampleName(uint which){ return colIDs[which]; }

	int SmplNum() { return (int)mat.size(); }
	int FtNum() {
		if (mat.size() >= 1) { return (int)mat[0].size(); }
		else { return 0; }
	}
	void estimateModuleAbund(char ** args, int argc);
	void estimateModuleAbund(options*);
	void resizeMatRows(uint x,mat_fl def=(mat_fl)0);
	//for the R module, all used for rarefactions only
	void addRow(vector<mat_fl>);//idea is that a single row is added on to the matrix
	void setSampleNames(vector<string> in) { colIDs = in; }
	void setRowNames(vector<string> in) { rowIDs = in; }
	vector < string > getSampleNames(){ return(colIDs); }
	vector < string > getRowNames(){ return(rowIDs); }
	//void addCount(string, int, mat_fl);

	double getMinColSum();
	column getMinColumn(uint offset = 0);
	vector< pair <double, string>> getColSums(bool sorted = false);
	vector<double> getCSum() { return colSum; }
	void writeColSums(string outF);
protected:
	//subroutines
	//reads the number of columns and checks in first few lines
	//void readColNms(ifstream& in);
	void readColNms(stringstream& in); //PR
	//int iniCols(ifstream& in);
	int iniCols(stringstream& in); //PR
	void read_subset_genes(const string);
	void read_hierachy(const string );
	void addColumn(string);
	void readModuleFile(const string&);
	vector<mat_fl> getRowSums();
	void ini_mat();

	//storage
	vector< vector< mat_fl > > mat;
	vector< string > rowIDs,colIDs;
	unordered_map<string, int> colID_hash, rowID_hash;
	int maxCols;//number samples
	//vector<HMat*> HI;
	LvlUp LUp;
	int maxLvl;
	vector<string> LvlNms;
	string sampleNameSep;
	GeneIDidx subset;
	bool doSubsets, doHigh;
	vector<double> colSum;

	vector< pair <double, string>> colsums;
};

/***************MORESTRUCTS*****************/
struct rareStruct{
	int i;
	DivEsts* div;
	vector<string> cntsName;
	vector<vector<rare_map>> cnts;
	string skippedNames;
	vector<string> IDs;
	
};

/*
struct job {
  std::future <rareStruct*> fut;
  bool inUse = false;
};
*/

