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

CRITICAL_SECTION CriticalSection; 

Param param;
int main(int argn, char * argv[]) {
	vector<string> OUT_st;
	param.AddNoname("NThreads");
	param.TakeAgrs(argn,argv);
	param.ReadNoname(OUT_st);
	int NThread = 1;
	if(OUT_st.size()>0) NThread = atoi(OUT_st[0].data());
	cout << "Program parametrs: NThreads" <<endl;
	cout <<"\tNThreads: "<< NThread << endl; 
	{
		ofstream temp("a.txt",ios::trunc);
		temp.close();

	}
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
		for(int j = 0; j<NThread; j++, i++) {
			El[i].DefineCell(cell);
			ElList[j].push_back(&El[i]);
			if(i >= size)
				break;
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

	int size = labbuf.size();

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
	while (!file.eof()) {
		file.getline(buf, MAX_LINE);

		for (int i = 0; i<NAtoms; i++) {
			Point p;
			file >> buf;
			p.a[0] = atof(buf);
			file >> buf;
			p.a[1] = atof(buf);
			file >> buf;
			p.a[2] = atof(buf);
			tempVec[i].push_back(p);
		}

		for (int i = 0; i<NAtoms; i++) {
			if (AllSteps != 0) {
				Point dp = tempVec[i][AllSteps] - tempVec[i][AllSteps - 1];
				for (int j = 0; j < 3; j++) {
					if (abs(dp.a[j]) > 0.5)
						tempVec[i][AllSteps].a[j] -= _sign(dp.a[j]);
					else dp.a[j] = 0;
				}
			}
		}
		AllSteps++;
		file.getline(buf, MAX_LINE);
		Counter++;
		if (Counter == 100) {
			Counter = 0;
			cout << "Step " << AllSteps << " readed." << endl;
		}
	}
	cout << "Last step " << AllSteps << " readed." << endl;

	for (int i = 0; i < NAtoms; i++) {
		//tempVec[i].pop_back();
		tempVec[i].shrink_to_fit();
		pList[i].AddVecPoints(tempVec[i], true);
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
	auto iter = (*El).begin();
	while(iter!= (*El).end()) {
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