// Console.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "../Includes/Includes.h"


using namespace std;
const int MAX_LINE = 128;
//list<Point> points;
//Point req[3];
DWORD WINAPI ThreadProc( LPVOID lpParameter );
Cell LoadXDATCAR(vector<Elipsoid> & pList);
void OutIns(Elipsoid & el);
int _sign(flo a);
void Analize_symmety(nsShelxFile::ShelxData & shelx, vector<Elipsoid> & pList);

CRITICAL_SECTION CriticalSection; 

OldParam param;
int main(int argn, char * argv[]) {
	vector<string> OUT_st;
	param.AddNoname("NThreads");
	param.TakeAgrs(argn,argv);
	param.ReadNoname(OUT_st);
	int NThread = 1;
	if(OUT_st.size()>0) NThread = atoi(OUT_st[0].data());
	cout << "Program parametrs: NThreads" <<endl;
	cout <<"\tNThreads: "<< NThread << endl; 
	vector<DWORD>   dwThreadIdArray;
    vector<HANDLE>  hThreadArray; 
	dwThreadIdArray.resize(NThread-1);
	hThreadArray.resize(NThread-1);
	vector<list<Elipsoid*>> ElList;
	if (!InitializeCriticalSectionAndSpinCount(&CriticalSection, 0x00000400) ) 
		return -1; 
	vector<Elipsoid> El;
	Cell cell = LoadXDATCAR(El);

	ElList.resize(NThread);
	int j = 0;
	int size = El.size();
	for(int i = 0; i < size;) {
		for(int j = 0; (j < NThread) && (i < size); j++, i++) {
			El[i].DefineCell(cell);
			ElList[j].push_back(&El[i]);
		}
	}
	cout << "Prefunction complete." << endl;

	for(int i = 0; i<(NThread-1); i++) {
		hThreadArray[i] = CreateThread( 
            NULL,                   // default security attributes
            0,                      // use default stack size  
            ThreadProc,       // thread function name
            &ElList[i],          // argument to thread function 
            0,                      // use default creation flags 
            &dwThreadIdArray[i]);   // returns the thread identifier 
        if (hThreadArray[i] == NULL) 
        {
           ExitProcess(3);
        }
    }

	ThreadProc(&ElList[NThread-1]);
    WaitForMultipleObjects(NThread-1, hThreadArray.data(), TRUE, INFINITE);
	{
		ifstream old("old.res");
		if (old.is_open() == true) {
			nsShelxFile::ShelxData shelx(old);
			cout << "Symmetry found." << endl;
			Analize_symmety(shelx, El);
			ofstream nfile("new.res");
			shelx.OutIns(nfile);
			nfile.close();

		}
		old.close();
	}
	for (int i = 0; i < (NThread); i++) {
		auto iter = ElList[i].begin();
		while (iter != ElList[i].end()) {
			OutIns(**iter);
			iter++;
		}
	}
	for(int i=0; i<NThread-1; i++)
    {
        CloseHandle(hThreadArray[i]);
    }
    DeleteCriticalSection(&CriticalSection);
	return 0;
}
Cell LoadXDATCAR(vector<Elipsoid> & pList)
{

	vector<string> labbuf;
	int NAtoms = 0;
	char * fName = "XDATCAR";
	ifstream file(fName);
	if (file.is_open()) {
		cout << "XDATCAR opened." << endl;
	}
	else {
		cout << "Error! XDATCAR not opened." << endl;
		system("pause");
		exit(10);
	}
	char buf[MAX_LINE];

	file.getline(buf, MAX_LINE);
	file.getline(buf, MAX_LINE);
	flo U[9];
	file >> U[0] >> U[1] >> U[2];
	file >> U[3] >> U[4] >> U[5];
	file >> U[6] >> U[7] >> U[8];
	Matrix mat(&U[0], 3, 3);
	Cell cell(mat, true);

	file.getline(buf, MAX_LINE);


	file.getline(buf, MAX_LINE);
	std::stringstream str(buf);
	while (!str.eof()) {
		string temp;
		str >> temp;
		labbuf.push_back(temp);
	}
	labbuf.pop_back();
	labbuf.shrink_to_fit();

	size_t size = labbuf.size();

	vector<int> counts(size);

	for (int i = 0; i<size; i++)
		file >> counts[i];
	file.getline(buf, MAX_LINE);
	{
		ofstream temp("a.ins", ios::trunc);
		temp << "CELL 0.71073 " << cell.Lat_dir(0) << " " << cell.Lat_dir(1) << " " << cell.Lat_dir(2) << " "
			<< cell.Angle_grad(0) << " " << cell.Angle_grad(1) << " " << cell.Angle_grad(2) << " " << endl;
		temp << "ZERR 4 0.001 0.001 0.001 0 0 0" << endl;
		temp << "SFAC";
		for (int i = 0; i<size; i++)
			temp << " " << labbuf[i];
		temp << endl << "UNIT";
		for (int i = 0; i<size; i++)
			temp << " " << counts[i];
		temp << endl;
		temp << endl;
		temp.close();
	}


	for (int i = 0; i < size; i++) {
		NAtoms += counts[i];
	}
	pList.resize(NAtoms);
	for (int i = 0, k = 0; i < size; i++) {
		for (int j = 1; j <= counts[i]; j++, k++) {
			pList[k].DefineLabel(labbuf[i], j, i+1);
		}
	}

	int Counter = 0;
	int AllSteps = 0;
	vector<vector<Point> > tempVec(NAtoms);
	bool cleared = false;
	while (!file.eof()) {
		file.getline(buf, MAX_LINE);

		for (int i = 0; i<NAtoms; i++) {
			Point p;
			file >> buf;
			p.a[0] = (flo)atof(buf);
			file >> buf;
			p.a[1] = (flo)atof(buf);
			file >> buf;
			p.a[2] = (flo)atof(buf);
			tempVec[i].push_back(p);
		}

		if (AllSteps != 0) {
			for (int i = 0; i<NAtoms; i++) {
				Point dp = tempVec[i][AllSteps] - tempVec[i][AllSteps - 1];
				for (int j = 0; j < 3; j++) {
					if (abs(dp.a[j]) > 0.5)
						tempVec[i][AllSteps].a[j] -= _sign(dp.a[j]);
					else dp.a[j] = 0;
				}
			}
		}
		AllSteps++;
		if (AllSteps == 2001 && cleared == false) { //!
			vector<vector<Point>> tv2; 
			tv2.resize(NAtoms);
			tv2.swap(tempVec);
			for (int i = 0; i < NAtoms; i++) {
				tempVec[i].push_back(tv2[i].back());
			}
			AllSteps = 1;
			cleared = true;
		}
		file.getline(buf, MAX_LINE);
		Counter++;
		if (Counter%100 == 0) {
			cout << "Step " << Counter << " readed." << endl;
		}
	}
	cout << "Last step " << Counter-1 << " readed." << endl;

	for (int i = 0; i < NAtoms; i++) {
		tempVec[i].pop_back();
		tempVec[i].shrink_to_fit();
		pList[i].AddVecPoints(move(tempVec[i]), true);
	}
	cout << "VecPoints added." << endl;

	return cell;
}

void OutIns(Elipsoid & el)
{
	ofstream file("a.ins",ios::app);
	if (file.is_open()) {
	} else {
		cout << "Error! File not opened." << endl;
	}

	file << el.OutShelxString() << endl;

	file.close();

}
DWORD WINAPI ThreadProc( LPVOID lpParameter )
{
	list<Elipsoid*> * El = reinterpret_cast<list<Elipsoid*> *>(lpParameter);
	auto iter = El->begin();
	while(iter!= El->end()) {
		Elipsoid * pEl = *iter;
		//iter->Correlation();
		//iter->Process();
		//iter->UProc();
		ELIPSOID_RESULT res = iter.operator*()->CalculateDinmat();
		if (res != ELIPSOID_RESULT::OK) {
			cout << "Error" << endl;
		}
		iter++;
	}
	iter = El->begin();

    EnterCriticalSection(&CriticalSection); 
	while(iter!= El->end()) {
		//OutIns(*iter);
		//Out(*iter);
		iter++;
	}
    LeaveCriticalSection(&CriticalSection);

return 0;
}
int _sign(flo a) {
	if (a < 0) return -1;
	else return 1;
}
void Analize_symmety(nsShelxFile::ShelxData & shelx, vector<Elipsoid> & pList) {
	size_t size_a = shelx.atom.size();
	size_t size_s = shelx.symm.size();
	size_t size_el = pList.size();

	for (int i = 0; i < size_el; i++) {
		char buf[5];
		strcpy(buf, pList[i].label.c_str());
		_strupr_s(buf);
		for (int j = 1; j < shelx.sfac.size(); j++) {
			if (shelx.sfac[j].compare(buf) != 0) continue;
			pList[i].type = j;
			break;
		}
	}

	for (int i = 0; i < size_a; i++) {
		for (int j = 0; j < 6; j++) {
			shelx.atom[i].dinmat.U[j] = 0;
		}
	}
	vector<nsShelxFile::Atom> atsh(std::move(shelx.GenerateSymmAtom()));
	size_t size_ats = atsh.size();
	for (int j = -1; j <= 1; j++) {
		for (int k = -1; k <= 1; k++) {
			for (int l = -1; l <= 1; l++) {
				if (j == 0 && k == 0 && l == 0) continue;
				for (int i = 0; i < size_ats; i++) {
					atsh.push_back(atsh[i]);
					(atsh.back()).point += Point(j, k, l);
				}
			}
		}
	}
	size_ats = atsh.size();
	for (int i = 0; i < size_el; i++) {
		nsShelxFile::Atom * ap = &atsh[0];
		flo dp = (pList[i].vecPoints[0] - shelx.cell.FracToCart() * (atsh[0].point)).r();
		for (int j = 1; j < size_ats; j++) {
			if (atsh[j].type != pList[i].type) continue;
			flo tempdp = (pList[i].vecPoints[0] - shelx.cell.FracToCart() * (atsh[j].point)).r();
			if (tempdp >= dp) continue;
			dp = tempdp;
			ap = &(atsh[j]);
		}
		ap->point = pList[i].fracCenter;
		ap->dinmat = pList[i].dinmat;
		ap->type += 256;
	} 
	{
		vector<nsShelxFile::Atom> tempA;
		tempA.swap(atsh);
		for (int i = 0; i < size_ats; i++) {
			if ((tempA[i].type & 256) != 0) {
				tempA[i].type -= 256;
				atsh.push_back(std::move(tempA[i]));
			}
		}
	}
	size_ats = atsh.size();
	vector<nsShelxFile::Atom> acpy(shelx.atom);
	for (int i = 0; i < size_a; i++) {
		acpy[i].dinmat = Dinmat();
		acpy[i].point = Point();
	}
	for (int i = 0; i < size_a; i++) {
		size_t count = 0;
		for (int j = 0; j < size_ats; j++) {
			if (strcmp(atsh[j].label, shelx.atom[i].label) != 0) continue;
			vector<nsShelxFile::Atom> atsh2(std::move(shelx.GenerateSymmAtom(atsh[j])));
			nsShelxFile::Atom * ap = &atsh2[0];
			flo dp = (shelx.cell.FracToCart() * (shelx.atom[i].point - (atsh2[0].point))).r();
			for (int k = 0; k < atsh2.size(); k++) {
				flo tempdp = (shelx.cell.FracToCart() * (shelx.atom[i].point - (atsh2[k].point))).r();
				if (tempdp >= dp) continue;
				dp = tempdp;
				ap = &(atsh2[k]);
			}
			acpy[i].point += ap->point;
			for (int k = 0; k < 6; k++)
				acpy[i].dinmat.U[k] += ap->dinmat.U[k];
			count++;
		}
		cout << "Atom " << i << " connected to " << count << " atoms." << endl;
		acpy[i].point /= count;
		for (int k = 0; k < 6; k++)
			acpy[i].dinmat.U[k] /= count;
	}

	shelx.atom.swap(acpy);
}
