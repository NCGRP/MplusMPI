#include "m+.hpp"
	
int MyDoRarify(int i, vector<vector<vector<int> > > AlleleList, std::set<int> AlleleSet, int CoreSize)
{
	std::stringstream inss;
	vector<int> CurrLoc;
	int d, j, M;
			
	//fill the ss header
	inss << "out";
	//for (j=0;j<CoreSize;++j) inss << "\tp" << j;
	for (j=0;j<CoreSize;++j) 
	{
		if (AlleleList[j][i].size() == 0) continue; //skip populations with no alleles (due to missing data)
		else inss << "\tp" << j;
	}
	inss << "\n";
	
	for (set<int>::iterator it=AlleleSet.begin(); it!=AlleleSet.end(); ++it)
	{
		d = *it; //dereference to get each unique allele from set
		inss << "a" << d; //add row label
		//for each accession in core
		for (j=0;j<CoreSize;++j)
		{
			if (AlleleList[j][i].size() == 0) continue; //skip populations with no alleles (due to missing data)
			else 
			{
				CurrLoc = AlleleList[j][i];	//all alleles at locus i in population j
				inss << "\t" << std::count(CurrLoc.begin(), CurrLoc.end(), d); //count number of allele d in pop j, add to stringstream
			}
		}
		inss << "\n";
	}
	

	/*
		for (set<int>::iterator it=AlleleSet.begin(); it!=AlleleSet.end(); ++it) cout << "AlleleSet[]=" << (*it) << "\n";
		cout << "AlleleList.size()=" << AlleleList.size() << "\n";
		for (j=0;j<CoreSize;++j)
		{
			cout << "AlleleList[" << j << "][" << i << "].size()=" << AlleleList[j][i].size() << "\n";
			for (unsigned int k=0;k<AlleleList[j][i].size();++k) cout << "AlleleList[j][i][]=" << AlleleList[j][i][k] << "\n";
			cout << inss.str() << "\n";
		}
	*/


	//send stringstream table to rtk codes for rarification	
	//rtkrare returns the median rarefied allele count, for the current locus, derived from all pops in core
	M = rtkrare(inss);
	return M;
}

int MyRarify(int i, vector<vector<vector<int> > > AlleleList, int sss)
{
	/*AlleleList structure = ActiveAlleleByPopList structure:
	  Pop1..r
		  locusarray1..n		
	Either 3D vector can be passed in as ActiveAlleleByPopList
	*/

	unsigned int j, k;
	int M;
	vector<int> CurrLoc;
	vector<vector<int> > CurrPop;  //contains vector of allelic states, for each population, for locus i only
	set<int> NewSet;
	
	//collect alleles at locus i, for all pops
	vector<vector<int> >().swap(CurrPop); //clear CurrPop
	for (j=0;j<AlleleList.size();j++)
	{
		CurrPop.push_back(AlleleList[j][i]); //extract alleles for all populations for this locus
	}
	
	//sort CurrPop by population/sample size, so that alleles can be counted from smallest pop to largest
	std::sort(CurrPop.begin(), CurrPop.end(), [](const vector<int> & a, const vector<int> & b){ return a.size() < b.size(); });
	
	//determine unique alleles across populations using a set
	//add a sample of sss alleles from each population to the set
	//examine populations in order from smallest to largest
	for (j=0;j<CurrPop.size();j++)
	{
		CurrLoc = CurrPop[j]; //all alleles for pop j, in order of increasing population size
		std::shuffle(std::begin(CurrLoc), std::end(CurrLoc), rng); //shuffle the order of the alleles so that a random sample of sss alleles can be taken for rarification

		for (k=0;k<sss;++k)
		{
			NewSet.insert(CurrLoc[k]); //place alleles into set to eliminate redundancies, returns pair<it,bool=true> if a new item is accepted
		}
	}
	
	if (NewSet.size() == 0) M=0;
	else M=NewSet.size();
	
	return M; //the set size after adding sss alleles from all populations for locus i is the rarified number of alleles

}

int MyCalculateDiversity(vector<vector<vector<int> > > AlleleList, vector<int> ActiveMaxAllelesList, std::string Standardize, std::string Rarify, double& RandomActiveDiversity, double& AltRandomActiveDiversity)
{
	/*AlleleList structure:
		  Pop1..r
			  locusarray1..n		*/
	int CoreSize = AlleleList.size();
	int NumLoci = AlleleList[0].size();
	int i, j, M;
	unsigned int k;
	vector<int> CurrLoc;
	vector<int> Mlist(NumLoci);
	set<int> AlleleSet;

	if (NumLoci == 0)
	{
		RandomActiveDiversity = -1; 
		AltRandomActiveDiversity = -1;
	}
	else
	{
		if (Rarify == "no")
		{
			//use set to eliminate redundancies
			for (i=0;i<NumLoci;i++)
			{
				//3. pass alleles from the same locus into a single set, for all populations in core, to remove redundancies
				AlleleSet.clear(); //clear AlleleSet
				for (j=0;j<CoreSize;j++)
				{
					CurrLoc = AlleleList[j][i];
					for (k=0;k<CurrLoc.size();++k)
					{
						AlleleSet.insert(CurrLoc[k]); //locus i for all population j
					}
				}
			
				if (AlleleSet.size() == 0) M=0;
				else M=AlleleSet.size();
				
				Mlist[i] = M;
			}
		}
		else if (Rarify == "yes")
		{
			for (i=0;i<NumLoci;i++)
			{
				
				/*** de novo coded rarification start ***/
				//sample no more than sss alleles from each accession in core, remove redundancies using set
				//add unique alleles to set in increasing order of population size

				M = MyRarify(i, AlleleList, sss);
				Mlist[i] = M;
				/*** de novo coded rarification end ***/

				
				/*** use rtk rarification start ***
					//determine unique alleles for current core using a set
					AlleleSet.clear(); //clear AlleleSet
					for (j=0;j<CoreSize;j++)
					{
						CurrLoc = AlleleList[j][i];
						for (k=0;k<CurrLoc.size();++k)
						{
							AlleleSet.insert(CurrLoc[k]); //unique alleles at locus i for all population j in core
						}
					}
			
					//perform the rarification, add the rarified allele count at this locus (i)
					if (AlleleSet.size() == 0) M = 0;
					else M = MyDoRarify(i, AlleleList, AlleleSet, CoreSize);
				
					Mlist[i] = M;
				*** use rtk rarification end ***/

			} //NumLoci
		}
		/*
				//print out the Mlist
				for (unsigned int i=0;i<Mlist.size();++i)
				{
					cout << "Mlist[" << i << "]=" << Mlist[i] << " ";
				}
				cout << "\n";
		*/




		
		//5. standardize the M values to the maximum possible number of alleles at that locus, 
		//and add them up to get final estimate of standardized allelic diversity in the core.
		//Then divide by the number of loci to get a number that is comparable across data sets.
		
		//calculate the standardized allelic diversity
		double SAD;
		double SADtemp = 0; //SADtemp is the summed value of standardized M across all loci for a single subcore
		for (i=0;i<NumLoci;i++)
		{
			SADtemp = SADtemp + ( (double) Mlist[i] / (double) ActiveMaxAllelesList[i] );
		}
		SAD = SADtemp / NumLoci;
		
		//calculate M by simply adding up the Mlist
		M = 0;
		for (k=0;k<Mlist.size();++k)
		{
			M = M + Mlist[k];
		}
		
		//6. Determine the way variables are updated so that the value for the desired optimality 
		//criterion is RandomActiveDiversity, while the other is AltRandomActiveDiversity.  
		//This way, the output can contain information on both M+ and M criteria, 
		//although only one will be used for optimization.
		if (Standardize == "yes")
		{
			RandomActiveDiversity = SAD;
			AltRandomActiveDiversity = M;
		}
		else if (Standardize == "no")
		{
			RandomActiveDiversity = M;
			AltRandomActiveDiversity = SAD;
		}
	}
	return 0;
}

void printProgBar( int percent)
{
	std::string bar;

	for(int i = 0; i < 50; i++){
		if( i < (percent/2)){
      		bar.replace(i,1,"=");
   		 }
    	else if( i == (percent/2)){
     		bar.replace(i,1,">");
    	}
    	else{
      		bar.replace(i,1," ");
   		 }
	}

  	std::cout<< "\r" "[" << bar << "] ";
  	std::cout.width( 3 );
  	std::cout<< percent << "% complete "<< std::flush;
}

// returns count of non-overlapping occurrences of 'sub' in 'str'
int countSubstring(const std::string& str, const std::string& sub)
{
    if (sub.length() == 0) return 0;
    int count = 0;
    for (size_t offset = str.find(sub); offset != std::string::npos;
	 offset = str.find(sub, offset + sub.length()))
    {
        ++count;
    }
    return count;
}


//M+	
void mp(
	int MinCoreSize,
	int MaxCoreSize,
	int SamplingFreq,
	int NumReplicates,
	char* OutFilePath,
	std::string Rarify,
	std::string Kernel,
	vector<int> KernelAccessionIndex,
	vector<int> AccessionNameList,
	vector<vector<vector<int> > > ActiveAlleleByPopList,
	vector<vector<vector<int> > > TargetAlleleByPopList,
	vector<int> ActiveMaxAllelesList,
	vector<int> TargetMaxAllelesList,
	vector<std::string> FullAccessionNameList
	)	
{

	//PERFORM INITIAL MPI STUFF
	MPI_Status status; //this struct contains three fields which will contain info about the sender of a received message
						 // MPI_SOURCE, MPI_TAG, MPI_ERROR
	
	//MPI::Init ();  //Initialize MPI.
	int nproc = MPI::COMM_WORLD.Get_size ( );  //Get the number of processes.
	int procid = MPI::COMM_WORLD.Get_rank ( );  //Get the individual process ID.
	

	//set up vectors to fill with results
	//below is a stupid way to calculate the number of rows in the output file, value l (which = V1) 
	//used to monitor progress and as the maximum vector index for shared output vectors
	int l=0;
	for (int i=MinCoreSize;i<MaxCoreSize+1;i=i+SamplingFreq)
	{
		for (int j=0;j<NumReplicates;j++)
		{
			l++;
		}
	}
	
	double V1 = (double)l; //(MaxCoreSize - MinCoreSize + 1)*NumReplicates; //number of rows in output vectors
	vector<vector<double> > Results(V1, vector<double>(9)); //will contain numerical results
	vector<vector<string> > Members(V1); //will contain core set members
	
	//***MPI:  RECEIVE RESULTS AT MASTER 0
	//receive values from any slave, in any order, exiting when the number of 'receives' = the top vector size
	if ( procid == 0 ) 
	{
		//set up variables for monitoring progress
		int percent; //percent of analysis completed
		int progindex = 0;  //index to monitor progress, percent = 100*(progindex/l)

		//receive and process results from slave processors
		unsigned int i = 0;
		while (i<2*(Results.size())) //two receives per row
		{
			//probe the incoming message to determine its tag
			int nchar; //will contain the length of the char array passed with tag=1
			int vchar; //will contain the length of the vector passed with tag=0
			int tag; //tag of message from sender
			int source; //procid of sender
			
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			//MPI_Get_count(&status, MPI_CHAR, &nchar); //probes the length of the message, saves it in nchar
			tag = status.MPI_TAG; //the tag defines which kind of comm it is, a vector of stats (0=resvec()) 
			                      //or a char array describing the members of the core (1=cc)
			source = status.MPI_SOURCE; //determine the source of the message so that you can define which sender to Recv from.  This will avoid an intervening message coming in after the MPI_Probe with a different length, causing a message truncated error.
			
			if (tag == 0)
			{
				//determine the length of the message tagged 0
				MPI_Get_count(&status, MPI_DOUBLE, &vchar);

				//cout <<" vchar="<<vchar<<" tag="<<tag<<" MPI_SOURCE="<<status.MPI_SOURCE<<" MPI_ERROR="<<status.MPI_ERROR<<"\n";

				//receive the vector of results, tagged 0, from:
				//MPI_Send(&resvec[0], resvec.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
				vector<double> t(10);
				MPI_Recv(&t[0], vchar, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
			
				//load data from vector received onto Results, row number is last item t[9]
				for (int j=0;j<9;++j)
				{
					Results[ t[9] ][j] = t[j];
				}
				t.clear();
			}

			else if (tag == 1)
			{
				//determine the length of the message tagged 1
				MPI_Get_count(&status, MPI_CHAR, &nchar); //probes the length of the message, saves it in nchar

				//cout <<" nchar="<<nchar<<" tag="<<tag<<" MPI_SOURCE="<<status.MPI_SOURCE<<" MPI_ERROR="<<status.MPI_ERROR<<"\n";
				
				//receive the vector<string> of the core set, tagged 1, from:
				//MPI_Send(&m[0], nchar, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
				//vector<string> m(nchar);
				char m[nchar];
				MPI_Recv(&m[0], nchar, MPI_CHAR, source, 1, MPI_COMM_WORLD, &status);
			
				//load core set onto Members
				//1. convert char array into a string
				string mstr(m);
			
				//2. split string on delimiter ',<!>,'
				string delim = ",<!>,";
				vector<string> mvec( countSubstring(mstr, delim) );
				unsigned int st = 0;
				std::size_t en = mstr.find(delim);
				int k = 0;
				while (en != std::string::npos)
				{
					mvec[k] = mstr.substr(st, en-st);
					st = en + delim.length();
					en = mstr.find(delim,st);
					++k;
				}
				string z = mstr.substr(st); //get row number as last item in mstr
				int zz = atoi(z.c_str()); //convert string to c-string then to int
			
				//3. load onto Members
				Members[zz] = mvec;
				
				//4. clean up
				memset(m, 0, nchar);; mstr=""; mvec.clear();
				

			}


			++i;
			
			//display progress
			progindex = progindex + 1;
			percent = 100*( progindex/(V1*2) ); //number of rows X 2 repeats needed to complete search
			printProgBar(percent); 
		}
	}//***MPI: END MASTER RECEIVE***/

	
	/***MPI:  SEND RESULTS FROM SLAVE PROCESSES***/
	else if ( procid != 0 )
	{
		unsigned int r; //r = core size, 
		//int nr, RandAcc, b, bsc, plateau; //nr = controller to repeat NumReplicates times
		int RandAcc, b, bsc, plateau; //nr = controller to repeat NumReplicates times
								//row = result vector row number, bsc = holds best sub core member, and other indexed accessions
									//plateau = index of the number of reps in optimization loop with same diversity value
		double RandomActiveDiversity;
		double AltRandomActiveDiversity;
		double StartingRandomActiveDiversity;
		double StartingAltRandomActiveDiversity;
		double RandomTargetDiversity;
		double AltRandomTargetDiversity;
		double StartingDiversity;
		double TempAltOptimizedActiveDiversity;
		double AltOptimizedActiveDiversity;
		double OptimizedTargetDiversity;
		double AltOptimizedTargetDiversity;
		double best;
		double nnew;
		vector<vector<vector<int> > > AlleleList;
		vector<vector<vector<int> > > CoreAlleles;
		vector<vector<vector<int> > > TdTempList;
		vector<vector<vector<int> > > BestSubCoreAlleles;
		std::string Standardize = "yes";  //a run that mimics the MSTRAT approach can be accomplished by setting Standardize="no", ensuring that Rarify="no", and setting up the var file so that each column in the .dat file is treated as a single locus, rather than two (or more) adjacent columns being treated as a single codominant locus.
		vector<int> AccessionsInCore;
		vector<int> AccessionsInSubCore;
		vector<int> BestSubCore;
		vector<int> BestSubCoreRevSorted;
		vector<int> TempList;
		vector<int> TempList2;
		vector<int> bestcore;
		vector<int> pbcore;
		vector<std::string> TempListStr;
	
		/*
		//seed the random number generator for each processor
		int tt;
		tt = (time(NULL));
		srand ( abs(((tt*181)*((procid-83)*359))%104729) );
		*/
		
		//seed the random number generator for each processor C++11
		unsigned long long seed = (unsigned long long)chrono::high_resolution_clock::now().time_since_epoch().count();
		//std::random_device rd;     // only used once to initialise (seed) engine
		std::mt19937_64 rng(seed*procid);    // random-number engine used (Mersenne-Twister x64 in this case), seeded with time*procid

		//do parallelization so that each rep by core size combo can be
		//handled by a distinct thread.  this involves figuring out the total
		//number of reps*coresizes taking into account the SamplingFreq
		int rsteps = 1 + floor( (MaxCoreSize - MinCoreSize) / SamplingFreq ); //number of steps from MinCoreSize to MaxCoreSize

		//***MPI: figure out where to start and stop loop for each processor
		int nreps = rsteps*NumReplicates;
		int count = nreps/(nproc-1); //p-1 assumes a master, i.e. one less processor than total
		int start = (procid-1) * count; //procid-1 makes you start at 0, assumes master is p0
		int stop;
		
		if (nreps % (nproc-1) > (procid-1))
		{
			start += procid - 1;
			stop = start + (count + 1); 
		}
		else
		{
			start += nreps % (nproc-1);
			stop = start + count;
		}
		
		//iterate thru the relevant rows
		for (int rnr=start;rnr<stop;++rnr)
		{
			r = MinCoreSize + ((rnr / NumReplicates) * SamplingFreq); //int rounds to floor
			
			//develop random starting core set
			//clear AccessionsInCore and set size
			AccessionsInCore.clear();
			AccessionsInCore.resize(r);
			
			//add kernel accessions to core, if necessary
			if (Kernel == "yes")
			{
				for (unsigned int i=0;i<KernelAccessionIndex.size();i++)
				{
					AccessionsInCore[i] = KernelAccessionIndex[i];
				}
			}

			//clear TempList and set size					
			TempList.clear();
			TempList.resize( AccessionNameList.size() );
			
			//set list of available accessions in TempList, by erasing those already in the core
			TempList = AccessionNameList;
			//expunge the kernel accessions, so they are not available for random addition below
			//KernelAccessionIndex has been reverse sorted so you don't go outside range after automatic resize by .erase
			for (unsigned int i=0;i<KernelAccessionIndex.size();i++)
			{
				b = KernelAccessionIndex[i];
				TempList.erase(TempList.begin()+b);
			}
		
			//randomly add accessions until r accessions are in the core. if there is a kernel, include those (done above)
			//plus additional, randomly selected accessions, until you get r accessions
			for (unsigned int i=KernelAccessionIndex.size();i<r;i++)
			{
				//choose an accession randomly from those available
				//RandAcc = rand() % TempList.size();
				std::uniform_int_distribution<int> uni(0,TempList.size()); // define random number generator interval to select items from TempList
				RandAcc = uni(rng); //rng is mt19937_64
				//add it to the list
				AccessionsInCore[i] = TempList[RandAcc];
			
				//remove it from the list of available accessions
				TempList.erase(TempList.begin()+RandAcc);
			}
	
			//assemble genotypes for random core and calculate diversity
			//1. put together initial list of active alleles
			CoreAlleles.clear();
			CoreAlleles.resize( AccessionsInCore.size() );
			for (unsigned int i=0;i<AccessionsInCore.size();i++)
			{
				b = AccessionsInCore[i];
				CoreAlleles[i] = ActiveAlleleByPopList[b];
			}

			//2. calculate diversity from random selection at active loci
			AlleleList.clear();
			AlleleList = CoreAlleles;
	
			MyCalculateDiversity(AlleleList, ActiveMaxAllelesList, Standardize, Rarify, RandomActiveDiversity, AltRandomActiveDiversity);
			//in MyCalculateDiversity, latter two variables are updated as references
			//save them away in non-updated variables
			StartingRandomActiveDiversity = RandomActiveDiversity;
			StartingAltRandomActiveDiversity = AltRandomActiveDiversity;

			//3. calculate diversity from random selection at target loci
			AlleleList.clear();
			AlleleList.resize( AccessionsInCore.size() );
			for (unsigned int j=0;j<AccessionsInCore.size();j++)
			{
				b = AccessionsInCore[j];
				AlleleList[j] = TargetAlleleByPopList[b];
			}
			MyCalculateDiversity(AlleleList, TargetMaxAllelesList, Standardize, Rarify, RandomTargetDiversity, AltRandomTargetDiversity);


			//BEGIN OPTIMIZATION
			StartingDiversity = 0; //this is the diversity recovered during the prior iteration.
			plateau = 0; //count of the number of times you have found the best value, evaluates when you are
						 //stuck on a plateau, assuming acceptance criterion allows downhill steps
			//this is the iterations step, now an indefinite loop that is broken when 
			//no improvement is made during the course of the optimization algorithm
			//If r = kernel size = MinCoreSize then do no optimization but still calculate all variables.
			if (KernelAccessionIndex.size() == r)
			{
				//assemble genotypes for core
				//1. put together initial list
				CoreAlleles.clear();
				CoreAlleles.resize(r);
				for (unsigned int i=0;i<r;i++)
				{
					b = AccessionsInCore[i];
					CoreAlleles[i] = ActiveAlleleByPopList[b];
				}
				
				AlleleList = CoreAlleles;
				
				MyCalculateDiversity(AlleleList, ActiveMaxAllelesList, Standardize, Rarify, RandomActiveDiversity, AltRandomActiveDiversity);
				best = RandomActiveDiversity; //best is equivalent to OptimizedActiveDiversity
				AltOptimizedActiveDiversity = AltRandomActiveDiversity;
			}
			else
			{
				//do optimization
				while ( true )
				{
					//assemble genotypes for core
					//1. put together initial list
					CoreAlleles.clear();
					CoreAlleles.resize(r);
					for (unsigned int i=0;i<r;i++)
					{
						b = AccessionsInCore[i];
						CoreAlleles[i] = ActiveAlleleByPopList[b];
					}
		
					//2. go through all possible subsets of size r-1, one at a time, noting which is best.
					//If there is a kernel, do not swap out any of those accessions (they are retained as the
					//first KernelAccessionIndex.size() items in CoreAlleles).  Accomplished by starting for loop
					//at KernelAccessionIndex.size().
					best=0;
					for (unsigned int i=KernelAccessionIndex.size();i<CoreAlleles.size();i++)
					{
						//remove each item consecutively from the list of all populations in the core
						AlleleList.clear();
						TdTempList.clear();
				
						TdTempList = CoreAlleles; //swap to temporary vector
						TdTempList.erase( TdTempList.begin() + i);
						AlleleList = TdTempList;
			
						TempList2.clear();
						TempList2 = AccessionsInCore;
						TempList2.erase(TempList2.begin() + i);
						AccessionsInSubCore = TempList2;

						/*Data structure for SubCoreAlleles:
						SubCore 1..r
							Population 1..(r-1)
								AlleleArray 1..NumLoci		
			
						--3. fuse alleles from the same locus into a single array, for all accessions, for the current subcore
						--4. assemble a list of diversity (M) for each locus separately
						--5. standardize the M values to the maximum possible number of alleles at that locus, and add them up to get final estimate of standardized allelic diversity in the core.  then divide by the number of loci to get a number that is comparable across data sets.
						--5.5. simultaneous to the calculation, keep track of which subcore is best
						*/
			
						MyCalculateDiversity(AlleleList, ActiveMaxAllelesList, Standardize, Rarify, RandomActiveDiversity, AltRandomActiveDiversity);
						nnew = RandomActiveDiversity;

						if (nnew >= best) // >= allows sideways movement during hill climbing
						{
							best = nnew;

							BestSubCore.clear();
							BestSubCore = AccessionsInSubCore;
							BestSubCoreAlleles.clear();
							BestSubCoreAlleles = AlleleList;
						}
					}  //for loop cycles thru all subcores

					//reverse sort BestSubCore to support easy assembly of pared TempList below
					BestSubCoreRevSorted = BestSubCore;
					std::sort(BestSubCoreRevSorted.begin(), BestSubCoreRevSorted.end(), std::greater<int>());
	
					/*
					6. take the subcore with greatest diversity and consecutively add each 
					possible additional accession from the base collection.  find the core of size r 
					(not r-1 subcore) that has the greatest diversity.

					suppress the IDs of those accessions found in the BestSubCore from the 
					list of all accessions to get a list of remaining accessions.*/
					TempList = AccessionNameList;
					for (unsigned int k=0;k<BestSubCoreRevSorted.size();k++)
					{
						bsc = BestSubCoreRevSorted[k];
						TempList.erase( TempList.begin() + bsc );
					}
			
					//shuffle the list of remaining accessions, so addition order is not predictable
					//std::random_shuffle (TempList.begin(), TempList.end());
					std::shuffle (std::begin(TempList), std::end(TempList), rng);
					
					//add each remaining accession consecutively, calculate diversity, test 
					//whether it is better than the prior one
					best = 0;
					for (unsigned int k=0;k<TempList.size();k++)
					{
						bsc = TempList[k];
			
						//define the core
						TempList2 = BestSubCore;
						TempList2.resize( TempList2.size() + 1 );
						//TempList2.push_back(i);
						TempList2[TempList2.size()-1] = bsc; //add new accession to last vector element
						AccessionsInCore = TempList2;
			
						//assemble the allelelist for the core
						TdTempList = BestSubCoreAlleles;
						TdTempList.resize( TdTempList.size() + 1 );
						TdTempList[TdTempList.size()-1] = ActiveAlleleByPopList[bsc];
						AlleleList = TdTempList;
			
						//calculate diversity
						MyCalculateDiversity(AlleleList, ActiveMaxAllelesList, Standardize, Rarify, nnew, TempAltOptimizedActiveDiversity); 
		
						
						


						
						for (unsigned int  z=0;z<AccessionsInCore.size();++z) cout << AccessionsInCore[z] << ",";
						cout << "  ";
						for (unsigned int  z=0;z<bestcore.size();++z) cout << bestcore[z] << ",";					
						cout << "  " << nnew << "  " << best << "  " << StartingDiversity << "  ";
						for (unsigned int  z=0;z<pbcore.size();++z) cout << pbcore[z] << ",";
						cout << "\n";					



						
						
						
						
						//test whether current diversity is higher than the best diversity found so far
						if (nnew >= best) // >= allows sideways movement during hill climbing
						{
							best = nnew;
							bestcore = AccessionsInCore;
							//save the alternative diversity value for the best core
							AltOptimizedActiveDiversity = TempAltOptimizedActiveDiversity;
						}
					}

					AccessionsInCore = bestcore; //define starting variable for next MSTRAT iteration
		
					
					
					
					cout << "best=" << best << " StartingDiversity=" << StartingDiversity << "\n";
					
					
					//when Rarify="yes", best may not equal StartingDiversity for the same core set
					//because diversity measure involves rarified sampling of populations.
					//In this case, allow loop to break when the same set of pops found in the
					//prior iteration is found again
					if (Rarify == "yes")
					{
						vector<int> v1 = bestcore;
						std::sort(v1.begin(), v1.end());
						vector<int> v2 = pbcore;
						std::sort(v2.begin(), v2.end());
						if ( v1 == v2 )
						{
							plateau++;
							if (plateau > 0) break;
						}
					}
					
					//if there has been no improvement from the prior iteration, you have reached
					// the plateau and should exit the repeat
					if (best == StartingDiversity) 
					{
						plateau++;
						if (plateau > 0) break;
					}
					else if (best > StartingDiversity) 
					{
						//update starting values and repeat
						StartingDiversity = best;
						pbcore = bestcore;
					}

				} //while(true) endless loop
			}

			//7. Calculate diversity at target loci
			//assemble the target loci allelelist for the accessions in the best core
			AlleleList.clear();
			AlleleList.resize( AccessionsInCore.size() );
			for (unsigned int j=0;j<AccessionsInCore.size();j++)
			{
				b = AccessionsInCore[j];
				AlleleList[j] = TargetAlleleByPopList[b];
			}
	
			//calculate diversity at target loci based upon the optimized core selection
			MyCalculateDiversity(AlleleList, TargetMaxAllelesList, Standardize, Rarify, OptimizedTargetDiversity, AltOptimizedTargetDiversity);

			//8. Assemble stats for optimized core and add to output vectors
			//create a list of accession names from the list of accession ID's in AccessionsInCore
			sort( AccessionsInCore.begin(), AccessionsInCore.end() );
			
			TempListStr.clear();
			TempListStr.resize(r);
			for (unsigned int i=0;i<AccessionsInCore.size();i++)
			{
				b = AccessionsInCore[i];
				TempListStr[i] = FullAccessionNameList[b];
			}

			/***MPI: BUILD & SEND RESULTS VECTOR***/
			//load the variables onto the results vectors
			
			//no need to calculate row number, it is the same as rnr, formula saved because it might be useful later
			//row = ((r - MinCoreSize)*NumReplicates) + nr - ( (NumReplicates*(SamplingFreq-1))*( (r-MinCoreSize)/SamplingFreq ) );
			// (r - MinCoreSize)*NumReplicates) + nr specifies row number if SamplingFreq=1
			// (NumReplicates*(SamplingFreq-1)) specifies a step value to correct when SamplingFreq>1
			// ( (r-MinCoreSize)/SamplingFreq ) specifies the replicate on core size, accounting for SamplingFreq
			// see file Calculation of row value.xlsx for development of the 'row' index
			
			//put results 0-8 into a vector, resvec, return row as last item
			vector<double> resvec(10);

			resvec[0] = double(r);
			resvec[1] = StartingRandomActiveDiversity;//RandomActiveDiversity;
			resvec[2] = best; //equivalent to OptimizedActiveDiversity
			resvec[3] = RandomTargetDiversity;
			resvec[4] = OptimizedTargetDiversity;
			resvec[5] = StartingAltRandomActiveDiversity;//AltRandomActiveDiversity;
			resvec[6] = AltOptimizedActiveDiversity;
			resvec[7] = AltRandomTargetDiversity;
			resvec[8] = AltOptimizedTargetDiversity;
			resvec[9] = double(rnr);
			
			
			//cout<<"MPI_Rank="<<MPI_Rank<<" 
			
			
			//send result vector to master 0, send row number, rnr, as last element.
			//message is tagged as 0
			//here you are pointing to the first element, then returning resvec.size() doubles-
			//worth of memory from that starting location.
			MPI_Send(&resvec[0], resvec.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			/***MPI: END BUILD & SEND RESULTS VECTOR***/


			/***MPI: BUILD & SEND MEMBERS VECTOR***/
			//add row number as last item in TempListStr
			TempListStr.resize(TempListStr.size()+1);
			stringstream ss;
			ss << rnr;	//convert int to stringstream to string			
			TempListStr[ TempListStr.size() - 1 ] = ss.str();
		
			//convert vector<string> to a single, ',<!>,' delimited, string
			string concat;
			for (unsigned int i=0;i<TempListStr.size();++i)
			{
				concat += TempListStr[i]; //add vector element
				if (i<TempListStr.size()-1) concat += ",<!>,"; //add delimiter, except for last item
			}
			//convert the string to a char array
			char cc[concat.size()+1];
			strcpy(cc, concat.c_str());
			
			//send the char array to master0 tagged as 1
			//tagged as 1 to distinguish from result vector send
			MPI_Send(&cc, sizeof(cc), MPI_CHAR, 0, 1, MPI_COMM_WORLD);

		} //end for loop over rows
	} //***MPI:  END SEND
	

	/*MPI: MASTER 0 WRITES OUTPUT*/
	if ( procid == 0 )
	{
		//set up file stream for output file
		ofstream output; 
		output.open(OutFilePath);
		output.close(); //quick open close done to clear any existing file each time program is run
		output.open(OutFilePath, ios::out | ios::app); //open file in append mode
		output << "core size	random reference diversity	optimized reference diversity	random target diversity	optimized target diversity	alt random reference diversity	alt optimized reference diversity	alt random target diversity	alt optimized target diversity	core members" << "\n";
		
		//write out results row by row
		for (int i=0;i<V1;i++)
		{
			//write variables
			output 	<< Results[i][0] 
					<< "	" << Results[i][1] 
					<< "	" << Results[i][2] 
					<< "	" << Results[i][3] 
					<< "	" << Results[i][4] 
					<< "	" << Results[i][5] 
					<< "	" << Results[i][6] 
					<< "	" << Results[i][7] 
					<< "	" << Results[i][8] 
					<< "	" << "(";
			//write Accessions retained
			for (unsigned int j=0;j<Members[i].size();j++)
			{
				if ( j==(Members[i].size() - 1) )
				{
					//add trailing parentheses and move to next row
					output << Members[i][j] << ")\n";
				}
				else
				{
					output << Members[i][j] << ",";
				}
			}
		}
	
		//wrap up write step
		output.close();
	} /***MPI: END MASTER WRITE***/
	
	//Terminate MPI.
	//MPI::Finalize ( );

}