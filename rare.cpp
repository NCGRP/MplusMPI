#include "m+.hpp"
#include "rare.hpp"

//code in rare.cpp is derived from https://github.com/hildebra/Rarefaction
//Paul Saary, Kristoffer Forslund, Peer Bork, Falk Hildebrand; RTK: efficient rarefaction 
//analysis of large datasets, Bioinformatics, Volume 33, Issue 16, 15 August 2017, Pages 
//2594–2595, https://doi.org/10.1093/bioinformatics/btx206

int medianDiv(vector<DivEsts*>& inD, bool printDIV, options* opts)
{
	int M; //rarified allele count

	vector<uint> cnts;
	vector<vector<uint> > rr(inD.size(),vector<uint>()); //holds the table reporting rarified allele count for each pop
	for (size_t i=0;i<inD.size();i++){
		for(uint j=0;j<inD[i]->cntsx.size();j++){
			for (uint k=0;k<inD[i]->cntsx[j].size();++k){
				//extract cnts, the listing of allele frequencies of occurrence from the vector of DivEsts divvs
				cnts = inD[i]->cntsx[j][k];
				rr[i] = cnts; //add current sample allele counts to table of all sample allele counts
			}
		}
	}

	/*
		//print out table of allele counts
		for (uint z=0;z<rr.size();++z)
		{
			cout << "inD.size()=" << inD.size() << " opts->depth.size()=" << opts->depth.size() 
				<< " rr.size()=" << rr.size() << " rr[" << z << "].size()=" << rr[z].size() << "\n";
			for (uint zi=0;zi<rr[z].size();++zi)
			{
				cout << rr[z][zi] << "\t";
			}
			cout << "\n";
		}
	*/
		
		//sum "columns" in rr, if column sum > 0, an allele of that column type was found, so index up M
		M=0;
		for (uint i=0;i<rr[0].size();++i)
		{
			for (uint j=0;j<rr.size();++j)
			{
				if (rr[j][i]>0) 
				{
					M++;
					break; //leave this column if an allele is found
				}
			}
		}

	return M;
}

int Matrix::iniCols(stringstream& in)
{
	int ini_ColPerRow = 0;
	int cnt = 0;
	string line;
	while (getline(in, line, '\n')) {
		if (line.substr(0, 1) == "#" || line.length() < 2) { continue; }
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss, segments, '\t')) {
			ColsPerRow++;
		}

		if (cnt == 0) {
			ini_ColPerRow = ColsPerRow;
		}
		else {
			if (ColsPerRow != ini_ColPerRow) {
				cerr << "C1: Number of columns on line " << cnt << " is " << ColsPerRow << ". Expected " << ini_ColPerRow << " columns.\n" << line << endl;
				std::exit(6);
			}
		}
		cnt++;
		if (cnt > 10) { break; } //checks only the first 10 rows for same column number
	}
	if (ini_ColPerRow == 0) {
		cerr << "Could not find valid columns in matrix.. exiting\n"; exit(432);
	}
	//reset input stream
	in.clear();
	in.seekg(0, ios::beg);
	colIDs.resize(ini_ColPerRow - 1, "");
	colSum.resize(ini_ColPerRow - 1, 0.0);

	return ini_ColPerRow;
}

std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();
    //from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();


    for (;;) {
        int c = sb->sbumpc();
        switch (c) {
            case '\n':
                return is;
            case '\r':
                if (sb->sgetc() == '\n')
                    sb->sbumpc();
                return is;
            case EOF:
                // Also handle the case when the last line has no line ending
                if (t.empty())
                    is.setstate(std::ios::eofbit);
                return is;
            default:
                t += (char)c;
        }
    }
}

void Matrix::readColNms(stringstream& in) {
	string segments; string line;
	safeGetline(in, line);
	while (line.substr(0, 1) == "#") {
		safeGetline(in, line);
	}
	stringstream sso;
	int cnt2(-2);
	sso << line;
	while (getline(sso, segments, '\t')) {
		cnt2++;
		if (segments.length() > 150) {
	#ifdef notRpackage
			cerr << segments << " error!\n"; std::exit(5);
	#endif
		}
		if (cnt2 == -1) { continue; }
		colIDs[cnt2] = segments;
	}
}


//Matrix::Matrix(const string inF, const string xtra, bool highLvl)
//	: rowIDs(0), colIDs(0), maxCols(0), maxLvl(0), sampleNameSep(""), doSubsets(false), doHigh(highLvl)
Matrix::Matrix(stringstream& in)
	: rowIDs(0), colIDs(0), maxCols(0), maxLvl(0), sampleNameSep(""), doSubsets(false)
{
	//reads matrix from HDD
	//into mem.. careful with big files!
/*	if (doHigh){
		read_hierachy(xtra);
	}
	else if (xtra.length() > 2){
		read_subset_genes(xtra);
	}
	string line;
	ifstream in(inF.c_str());
	if (!in){
#ifdef notRpackage
cerr << "Cant open file " << inF << endl; std::exit(11);
#endif
	}
*/	
	//make a demo input file as a stringstream, for further processing
/*	stringstream in;
	in << "out\tp1\tp2\tp3\tp4\n";
	in << "a1\t0\t1\t2\t4\n";
	in << "a2\t0\t1\t0\t10\n";
	in << "a3\t52\t1\t0\t10\n";
	in << "a4\t0\t0\t0\t20\n";
	in << "a5\t1\t1\t2\t40\n";
*/
	
	string line;
	int ini_ColPerRow = iniCols(in);
	//int ini_ColPerRow = iniCols(in);
	//cout << "ini_ColPerRow=" << ini_ColPerRow << "\n";
	
	in.clear(); //remove EOF flag
	in.seekg (0, ios::beg); //reset the stringstream to beginning


	
	int cnt(0);


	if (ini_ColPerRow == 0) {
		cerr << "Empty matrix provided\n";
		return;
	}
	cnt = -1;
		
	readColNms(in);

/*	if (doHigh){
		for (int i = 0; i < maxLvl; i++){
			HI.push_back(new HMat(LvlNms[i], colIDs, vector<string>(0)));
		}
	}
*/

	string segments;
	string rowID = "";
	int geneCnt(0);
//	int cntNA(0);
	//vector<mat_fl> emptyVec(ini_ColPerRow, (mat_fl)0);
	mat.resize(ini_ColPerRow -1, vector<mat_fl>(0));
	while (safeGetline(in, line)) {
		//while (getline(in, line, '\n')) {
		cnt++;
		if (line.substr(0, 1) == "#"){ continue; } //skip lines labeled with #
		//if (line.length()<5){ continue; }  //original code assumed line label like "OTU 1", hence <5
		if (line.length()<2){  continue; } //skip lines that contain no samples, assumes line label like "a1", hence <1
		int cnt2(-2);
		vector<string> taxa(0);
		stringstream ss;
		ss << line;
		bool breaker(false);
		std::map<std::string, vector<string>>::iterator fnd;
		while (getline(ss, segments, '\t')) {
			cnt2++;
			if (cnt2 == -1){
				rowID = segments;
				if (rowID == "mapped"){//to deal with mocat files
					breaker = true;
					break;
				}
				if (doSubsets && subset.find(rowID) == subset.end()){
					breaker = true;
					break;
				}
				//set index for rows..
				rowID_hash[rowID] = cnt;
				rowIDs.push_back(rowID);

/*				if (doHigh){
					fnd = LUp.find(rowID);
					if (fnd == LUp.end()){//needs to be added to HMat
						taxa = vector<string>(maxLvl, "-1");
						cntNA++;
						if (cntNA < 100) {
							#ifdef notRpackage
							std::cout << "Row ID " << rowID << " is not in hierachy.\n";// \nAborting..\n"; std::exit(24);
							if (cntNA == 99) { std::cout << " ..\n"; }
							#endif
						}
					}
					else {
						taxa = (*fnd).second;
					}
					
					string lngTax = "";
					for (int tt = 0; tt < maxLvl; tt++){
					lngTax += taxa[tt] ;
					taxa[tt] = lngTax;
					lngTax += +";";
					}
					
				}
*/				geneCnt++;
				continue;
			}
			mat_fl tmp = atof(segments.c_str());

				//these two moved from if clause below
				mat[cnt2].push_back(tmp);
				colSum[cnt2] += (double)tmp;
/*			if (doHigh){//1:finds relevant rowID, extracts taxa; 2:add on all HighLvl mats
				for (int tt = 0; tt< maxLvl; tt++){
					HI[tt]->set(taxa[tt], cnt2, tmp);
				}
				colSum[cnt2] += (double)tmp;
			}
			else {//would make it a sparse matrix: if (tmp>0) (but looses indexing)
				mat[cnt2].push_back(tmp);
				colSum[cnt2] += (double)tmp;
			}
			
*/
		}

		if (breaker){
			continue;
		}
		if (cnt2 + 2 != ini_ColPerRow){
			cerr << "C2: Number of columns on line " << cnt << " is " << cnt2 + 2 << ". Expected " << ini_ColPerRow << " columns.\n";
			std::exit(64);
		}

	}
//	in.close();
	maxCols = (int)mat.size();
//	std::cout << "Read " << geneCnt << " genes" << endl;

	if (geneCnt == 0) {
		cerr << "No genes read.. aborting\n";
		exit(0);
	}
}

rareStruct* calcDivRar(int i, Matrix* Mo, DivEsts* div, options* opts,
        vector<vector<vector<uint>>>* abundInRow, vector<vector<vector<uint>>>* occuencesInRow){

    smplVec* cur        = Mo->getSampleVec(i);
    string curS         = Mo->getSampleName(i);
    div->SampleName     = curS;
    std::vector<vector<uint>> cnts;
    vector<vector<rare_map>> cntsMap(opts->depth.size());
    vector<string> cntsName(opts->depth.size());
    string skippedNames;
    bool wrAtAll(opts->write > 0);

    cur->rarefy(opts->depth, opts->output, opts->repeats,
            div, cntsMap, cntsName, skippedNames, abundInRow, occuencesInRow,
            opts->write, false, wrAtAll);

    // put everything in our nice return container
    rareStruct* tmpRS       = new rareStruct();
    tmpRS->div              = div;
    tmpRS->cnts             = cntsMap;
    tmpRS->cntsName         = cntsName;
    tmpRS->skippedNames     = skippedNames;
    tmpRS->i                = i;

    delete cur;
    return tmpRS;
}

/*
string printSimpleMap(const rare_map & vec, string outF, string id, vector<string> rowNames){
    // takes a map from the rarefaction function and writes the vector
    // to the disk.
    // this way we dont need memory to do
    ofstream out(outF.c_str(),  ios::binary);
    if (!out){
        cerr << "Couldn't open tmpvec file " << outF << endl; std::exit(99);
    }
    for(uint i = 0; i < rowNames.size(); i++){
        uint value = 0;
        auto fnd = vec.find(i);
        if(fnd != vec.end()){
            value = fnd->second;
        }
        out.write((char*) &value, sizeof(uint));

    }
    out.close();
    return outF;
}
*/
/*
void binaryStoreSample(options* opts, vector<vector<vector< string >>>& tmpMatFiles, rareStruct* tmpRS, vector<string>& rowNames, string outF,  vector<vector<string>>& cntsNames, bool reshapeMap){
    // store vectors of rarefied matrix on hard disk for memory reduction
    if(reshapeMap){
        vector < string > rowIDs = tmpRS->IDs;
        vector < uint > nrowIDs(rowIDs.size());
        // convert ids into integer vector, for remapping the values
        for(uint i = 0; i < rowIDs.size(); i++){
            nrowIDs[i] = std::stoi(rowIDs[i]);
        }
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            // reshape each vector, as some are zero, and we need to rematch values and rows
            // we use the row Ids which we created correctly when splitting the vector from the input file
            rare_map tmpVec;
            for (auto const& x : tmpRS->cnts[i][ii]){
                tmpVec[nrowIDs[x.first]] = x.second;
            }
            string vecLocation = printSimpleMap(tmpVec,	outF + "tmp_"  + std::to_string(ii)+"_" + std::to_string(i) + tmpRS->cntsName[i] + ".binary",	tmpRS->cntsName[i], rowNames);

            tmpMatFiles[i][ii].push_back(vecLocation);
        }
        }
    }else{
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            string vecLocation = printSimpleMap(tmpRS->cnts[i][ii],	outF + "tmp_" + std::to_string(ii) + "_" + std::to_string(i) + tmpRS->cntsName[i] + ".binary",	tmpRS->cntsName[i], rowNames);
            tmpMatFiles[i][ii].push_back(vecLocation);
        } }
    }
    // save sample name
    for(uint di = 0; di < opts->depth.size(); di++){
        if(tmpRS->cntsName[di].size() != 0){
            cntsNames[di].push_back(tmpRS->cntsName[di]);

        }
    }
}


void memoryStoreSample(options* opts, rareStruct* tmpRS, vector< vector< vector< rare_map >> >& MaRare,  vector<vector<string>>& cntsNames, bool reshapeMap){
    if(reshapeMap){
        vector < string > rowIDs = tmpRS->IDs;
        vector < uint > nrowIDs(rowIDs.size());
        // convert ids into integer vector, for remapping the values
        for(uint i = 0; i < rowIDs.size(); i++){
            nrowIDs[i] = std::stoi(rowIDs[i]);
        }
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            // reshape each vector, as some are zero, and we need to rematch values and rows
            // we use the row Ids which we created correctly when splitting the vector from the input file
            rare_map tmpVec;
            for (auto const& x : tmpRS->cnts[i][ii]){
                tmpVec[nrowIDs[x.first]] = x.second;
            }
            MaRare[i][ii].push_back(tmpVec);
        }
    }
    }else{
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            MaRare[i][ii].push_back(tmpRS->cnts[i][ii]);
            }
        }
    }
    // save sample name
    for(uint di = 0; di < opts->depth.size(); di++){
        if(tmpRS->cntsName[di].size() != 0){
            cntsNames[di].push_back(tmpRS->cntsName[di]);
        }
    }
}
*/

//minimal sum per column, used for rarefaction min raredep
double Matrix::getMinColSum() {
	if (colSum.size() > 0) {
		double minE = colSum[0];
		for (uint i = 0; i < colSum.size(); i++) {
			if (minE > colSum[i]) {
				minE = colSum[i];
			}
		}
		return minE;
	}
	else {
		return 0;
	}
}

vector<double> parseDepths(string a){
    std::vector<double> vect;
    std::stringstream ss(a);

    float i;

    while (ss >> i)
    {
        vect.push_back(i);
        cout << i << " ";
        if (ss.peek() == ',')
            ss.ignore();
    }


    return vect;
}

options::options(void) :input(""), output(""), mode(""),
    referenceDir(""), referenceFile(""),
    depth(), repeats(1), write(0), threads(1), writeSwap(true), verbose(false), oldMapStyle(false),
    modDB(""), modRedund(5), modEnzCompl(0.5f), modModCompl(0.5f), modWrXtraInfo(false), 
    modCollapse(false), calcCoverage(false),
	modDescr(""), modHiera(""), xtra("") {


 /*       bool hasErr = false;
        if (argc == 0) { return; }//could be the pseudo mod opt object

        //bool newIDthrs = false; string newIDthStr("");
        mode = argv[1];

        for (int i = 1; i < argc; i++)
        {
            if (!strcmp(argv[i], "-i"))
                input = argv[++i];
            else if (!strcmp(argv[i], "-o"))
                output = argv[++i];
            ///else if (!strcmp(argv[i], "-m"))
            else if (!strcmp(argv[i], "-d")){
                // split by any given komma
                depth = parseDepths(argv[++i]);
                depthMin = (long) *std::min_element(depth.begin(), depth.end());
            }else if (!strcmp(argv[i], "-r"))
                repeats = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-w"))
                write = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-t"))
                threads = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-ns"))   // no swap
                writeSwap = false;
            else if (!strcmp(argv[i], "-v"))
                verbose = true;
            //else if (!strcmp(argv[i], "-h"))
            //   helpMsg();
            //geneMat specific args
            else if (!strcmp(argv[i], "-map"))
                map = argv[++i];
            else if (!strcmp(argv[i], "-refD"))
                referenceDir = argv[++i];
            else if (!strcmp(argv[i], "-reference"))
                referenceFile = argv[++i];
            //module specific args
            else if (!strcmp(argv[i], "-refMods"))
                modDB = (argv[++i]);
            else if (!strcmp(argv[i], "-redundancy"))
                modRedund = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-enzymeCompl"))
                modEnzCompl = (float)atof(argv[++i]);
            else if (!strcmp(argv[i], "-moduleCompl"))
				modModCompl = (float)atof(argv[++i]);
			else if (!strcmp(argv[i], "-oldMapStyle"))
				oldMapStyle = true;
			else if (!strcmp(argv[i], "-writeExtraModEstimates"))
				modWrXtraInfo = true;
			else if (!strcmp(argv[i], "-collapseDblModules"))
                modCollapse = true;
            else if (!strcmp(argv[i], "-useCoverage"))//for gene catalog, default is counts
                calcCoverage = true;
			else if (!strcmp(argv[i], "-description")) //description of single modules
				modDescr = (argv[++i]);
			else if (!strcmp(argv[i], "-hiera")) // hierachy for modules
				modHiera = (argv[++i]);
			else if (!strcmp(argv[i], "-xtra"))
                xtra = (argv[++i]);
        }

        // sanity checks
        // we need input
        if (input == "") {//just set some defaults
            cerr << "Input must be specified\n";
            hasErr = true;
        }

        if (output == "") {//just set some defaults
            cerr << "Output must be specified\n";
            hasErr = true;
        }
*/
        // default to min*0.95
        if(depth.size() == 0){
            depth.push_back(0.95);
        }
/*
        if (hasErr) {
            cerr << "Use \"rtk -h\" to get full help.\n";
            exit(98);
        }
        if (mode == "rarefaction") {
            if (writeSwap) { mode = "swap"; 
            } else { mode = "memory"; }
        }
        if( write > repeats){
            cerr << "Can not create more tables than repeats of rarefaction\n";
            cerr << "Writes is now equal to repeats ("<< repeats <<")\n";
            write = repeats;
        }
*/
    }


smplVec::smplVec(const vector<mat_fl>& vec, const int nt) :IDs(0),totSum(0),
    num_threads(nt), richness(-1), Shannon(-1.f){
        double cumSum(0.f);
        for (uint i = 0; i < vec.size(); i++){
            cumSum += vec[i];
        }
//        if (verbose){
//            cerr << (long)cumSum << " allocating ";
//        }
        //arr = (int*) malloc((int)cumSum * sizeof(int));
        //arr = new unsigned short[(int)cumSum];
        arr.resize((long)cumSum);
//        if (verbose){
//            cerr << "memory";
//      }
        totSum = cumSum;
        long k(0); uint posInVec(-1);
        //numFeatures = 0;
        for (size_t i = 0; i< vec.size(); i++){

            long maxG = (long)vec[i];
            IDs.push_back( std::to_string(i));

            posInVec++;
            if (maxG == 0){ continue; }//not really a feature, doesnt need ot be counted as cat

            maxG += k; //bring up to current level
            for (; k<maxG; k++){
                arr[k] = posInVec;
            }
            //numFeatures++;
            // some simple numeric id for refernce in chao2 abundance calculations, so it behaves as the swap mode
        }
        posInVec++;
        numFeatures = posInVec;
//        if (verbose){
//            cerr << "..\n";
//        }
    }

/*
smplVec::smplVec(const string inF, const int nt) :IDs(0),totSum(0), num_threads(nt),
    richness(-1),Shannon(-1.f) {

        vector<double> vec;
        ifstream in(inF.c_str());
        string line; double cumSum(0.f);
        while(getline(in,line,'\n')) {
            string ID; float num;
            stringstream ss(line);
            ss>>ID; ss>>num;
            cumSum += num;
            vec.push_back(num); IDs.push_back(ID);
        }
        in.close();

//        if (verbose){
//            cerr<<(long)cumSum<<" allocating ";
//        }
        arr.resize((long)cumSum);
//        if (verbose){
//            cerr<<"memory";
//        }
        totSum = cumSum;
        long k(0); uint posInVec(0);
        for (size_t i = 0; i< vec.size(); i++){

            long maxG = (long)vec[i];
            maxG += k;
            if (maxG == 0){ continue; }//not really a feature, doesnt need ot be counted as cat
            for (; k<maxG; k++){
                arr[k] = posInVec;
            }
            posInVec++;
        }
        numFeatures = posInVec;
//       if (verbose){
//            cerr<<"..\n";
//        }
    }
*/

void smplVec::rarefy(vector<double> depts, string ofile, int rep,
        DivEsts* divs, std::vector<vector<rare_map>> & RareSample,
        vector<string>& retCntsSampleName, string& skippedSample,
        vector<vector<vector<uint>>>* abundInRow, vector<vector<vector<uint>>>* occuencesInRow,
        int writes,bool write, bool fillret){
    bool doShuffle = true;
    long curIdx = 0;
    long dep;
    // resize divvs
    divs->richness.resize(depts.size());
    divs->cntsx.resize(depts.size());
//    divs->shannon.resize(depts.size());
//    divs->simpson.resize(depts.size());
//    divs->invsimpson.resize(depts.size());
//    divs->chao1.resize(depts.size());
//    divs->eve.resize(depts.size());
    
    for(uint i = 0; i < depts.size(); i++){
        dep = (long)depts[i];

        if (dep > totSum){
            skippedSample = divs->SampleName;
//            if (verbose){cout<<"skipped sample, because rowSums < depth \n";}
            return;
        }
        //long curIdx=(long)totSum+1;
        for (int curRep=0;curRep<rep;curRep++){
            if(curIdx+dep >= (long) totSum || doShuffle == true){
                shuffle_singl();	
                curIdx = 0;
                doShuffle = false;
            }


            //count up
            vector<uint> cnts(numFeatures, 0);
            for (long i=(0+curIdx);i<(dep+curIdx);i++){
                cnts[arr[i]]++;
            }

		/*
            curIdx += dep;
            string t_out = ofile;
            if (rep!=1){
                std::ostringstream oss;
                oss<<curRep;
                t_out += "_" +oss.str();
            }

            if (curRep < writes && fillret) {
                rare_map cntsMap;
                // fill map:
                for(uint i = 0; i < cnts.size(); i++){
                    if(cnts[i] != 0){
                        cntsMap.insert( std::make_pair(i, cnts[i]) );
                    }
                }
                RareSample[i].push_back(cntsMap);
                if(curRep == 0){
                    retCntsSampleName[i] = divs->SampleName; // safe the sample name as well
                }
            }
       */     
            
            
            richness = 0;
            //below line is critical to normal functioning of rtk, not needed for M+
            //divs->richness[i].push_back(this->getRichness(cnts));

			divs->cntsx[i].push_back(cnts);
//            vector<double> three = this->calc_div(cnts, 4);
//            divs->shannon[i].push_back(three[0]);
//            divs->simpson[i].push_back(three[1]);
//            divs->invsimpson[i].push_back(three[2]);
//            divs->chao1[i].push_back(this->calc_chao1(cnts,1));
 //           divs->eve[i].push_back(this->calc_eveness(cnts));
            richness = 0;

            // save abundance for chao2 calculations later
/*            rarefyMutex.lock();
            for(uint im = 0; im < IDs.size(); im++){
                //sparse convertions in swap mode
                int id = std::stoi(IDs[im]);
                if(cnts[im] != 0){
                    abundInRow->at(i)[curRep][id]++;	
                    occuencesInRow->at(i)[curRep][id] += cnts[im];
                }
            }
            rarefyMutex.unlock();
*/
        }
    }
}

//below function is critical to normal functioning of rtk, not needed for M+
// vector version of diversity fuctions:
/*
long smplVec::getRichness(const vector<unsigned int>& cnts){
    richness = 0;
    for (size_t i = 0; i<cnts.size(); i++){
        if (cnts[i]>0){
            richness++;
        }
    }
    return richness;
}
*/

void smplVec::shuffle_singl() {
    //auto engine = std::default_random_engine{};
//    std::random_device rd;
//    auto engine = std::mt19937_64{rd()};
	unsigned long long seed = (unsigned long long)chrono::high_resolution_clock::now().time_since_epoch().count();
	auto engine = std::mt19937_64{ seed };
    std::shuffle(std::begin(arr), std::end(arr), engine);	
}

//int main(int argc, char* argv[])
int rtkrare(stringstream& in)
{

    options* opts = new options();

/*  string inF = opts->input;
    string outF = opts->output;
    string mode = opts->mode;
    uint numThr = opts->threads;
    string arg4 = "";//std::to_string(opts->depth[0]);
    string map = opts->map;
    string refD = opts->referenceDir;
*/
	//make a demo input file as a stringstream, for further processing
/*	stringstream in;
	in << "out\tp1\tp2\tp3\tp4\n";
	in << "a1\t0\t1\t2\t4\n";
	in << "a2\t0\t1\t0\t10\n";
	in << "a3\t52\t1\t0\t10\n";
	in << "a4\t0\t0\t0\t20\n";
	in << "a5\t1\t1\t2\t40\n";
*/

	Matrix* Mo = new Matrix(in);//no arg for outfile &  hierachy | gene subset
    vector<DivEsts*> divvs(Mo->smplNum(), NULL);
    vector< string > rowNames = Mo->getRowNames();
    // transform all percentage sizes into correct values
    for(uint i = 0; i < opts->depth.size(); i++){
        if (opts->depth[i] < 1.) {
            // rarefy to smallest colSum
            opts->depth[i] = (uint)round(opts->depth[i] * Mo->getMinColSum());
            if (opts->depth[i] == 0.0) {
                cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
                exit(1);
            }
        } 
    }

    // hold rarefied matrices
    // stores : repeats - sampels eg rows - vectors of columns
    //int NoOfMatrices = opts->write;
    //vector< vector< vector< rare_map > >> MaRare(opts->depth.size(), vector< vector< rare_map> > (opts->write));
    //std::vector<vector<string>> cntsNames(opts->depth.size(), vector<string>());


    // declare and size abundance vectors to hold the number of occurences of genes per row
    // this will be used for Chao2 estimation
    vector<vector<vector<uint>>> occuencesInRow(opts->depth.size(), vector<vector<uint>>(opts->repeats, vector<uint>(Mo->rowNum(),0)));
    vector<vector<vector<uint>>> abundInRow(opts->depth.size(), vector<vector<uint>>(opts->repeats, vector<uint>(Mo->rowNum(),0)));

	/*	
		for (unsigned int i=0;i<occuencesInRow.size();++i)
		{
			for (unsigned int j=0;j<occuencesInRow[i].size();++j)
			{
				for (unsigned int k=0;k<occuencesInRow[i][j].size();++k)
				{
					cout <<"i:j:k=" << i << ":" << j << ":" << k << " " << occuencesInRow[i][j][k] << "\n";
				}
			}
		}
	
		for (unsigned int i=0;i<abundInRow.size();++i)
		{
			for (unsigned int j=0;j<abundInRow[i].size();++j)
			{
				for (unsigned int k=0;k<abundInRow[i][j].size();++k)
				{
					cout <<"i:j:k=" << i << ":" << j << ":" << k << " " << abundInRow[i][j][k] << "\n";
				}
			}
		}
	*/







    //object to keep matrices
    //vector < vector < vector < string >> > tmpMatFiles(opts->depth.size(), vector<vector <string>>(opts->write));
    //cerr << "TH";
    // vector keeping all the slots
    vector < job > slots(opts->threads);

    // now start a async in each slot
    uint i          = 0; 

    size_t smpls = Mo->smplNum();
    bool breakpoint(true);
    while (breakpoint) {
        // check for any finished jobs
        for( uint j = 0; j < slots.size(); j++ ){
            if( i >= smpls){
                breakpoint = false;
                // break in case we have more slots than work
                break;
            }

            if(slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(20)) == std::future_status::ready){

                // move the information
                rareStruct* tmpRS;
                tmpRS               = slots[j].fut.get();
                divvs[tmpRS->i]     = tmpRS->div;
                string curS         = Mo->getSampleName(tmpRS->i);

                // add the matrices to the container
/*                if (opts->write > 0) {
                    if (opts->writeSwap) {
                       binaryStoreSample(opts, tmpMatFiles, tmpRS, rowNames, outF, cntsNames, false);
                    }
                    else {
                        memoryStoreSample(opts, tmpRS, MaRare, cntsNames, false);
                    }
                }
*/
                delete tmpRS;
                // free slot
                slots[j].inUse = false;
            }

            // open new slots
            // the parallel mechanism is: a sample (column), governed by i, is sent to independent threads
            // to calculate diversity measures using calcDivRar(i,...)
            if( slots[j].inUse == false){

                slots[j].inUse = true;
                // launch an async task
                DivEsts * div   = new DivEsts();
                slots[j].fut    = async(std::launch::async, calcDivRar, i, Mo, div, opts, &abundInRow, &occuencesInRow);

/*                // send user some feedback
                int thirds = 1;
                if(smpls > 5){
                    thirds = (int) ceil((smpls-3)/3);
                }
                if(i < 3 || i % thirds == 0  ){
                    cout << "At Sample " << i+1 << " of " << smpls << " Samples" << std::endl  ;
                    if(i > 3 && i % thirds == 0 ){
                        cout << "..." << std::endl ;
                    }
                }else if( i == 3){
                    cout << "..." << std::endl ;
                }
*/
                i++;

            }

        }


    }

    // wait for the threads to finish up.
    for(uint j = 0; j < slots.size(); j++){
        if(slots[j].inUse == false ){
            // only copy if there is work to be done
            continue;
        }
        slots[j].fut.wait();
        // move the information
        rareStruct* tmpRS;
        tmpRS               = slots[j].fut.get();
        divvs[tmpRS->i]     = tmpRS->div;
        string curS         = Mo->getSampleName(tmpRS->i);

/*        // add the matrices to the container
        if (opts->write > 0) {
            if (opts->writeSwap) {
               binaryStoreSample(opts, tmpMatFiles, tmpRS, rowNames, outF, cntsNames, false);
            }
            else {
                memoryStoreSample(opts, tmpRS, MaRare, cntsNames, false);
            }
        }
*/
                
        delete tmpRS;
        // free slot
        slots[j].inUse = false;
    }

	/*
		// print out some stats on the vector holding diversity values, divvs
		for (size_t i=0;i<divvs.size();i++){
			for(uint j=0;j<divvs[i]->cntsx.size();j++){
				for (uint k=0;k<divvs[i]->cntsx[j].size();++k){
					for (uint l=0;l<divvs[i]->cntsx[j][k].size();++l){
						cout << "divvs[" << i <<"]->cntsx[" << j << "][" << k << "][" << l << "]=" << divvs[i]->cntsx[j][k][l] << "\n";
						//cout << "rr[" << i << "] = " << rr[i] << "\n";
					}
				}
			}
		}
    */
    
    //calculate median allelic diversity from options::repeats=1 rarifications
    //only 1 rarification is used because replication is controlled elsewhere in M+ (so "median" is meaningless)
    //this value is the number of unique alleles found among all members of the subset
   	int M = medianDiv(divvs, true, opts);

	return M;
}
