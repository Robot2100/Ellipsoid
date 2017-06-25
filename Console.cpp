// Console.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "Elipsoid.h"
#include <sstream>


using namespace std;

//list<Point> points;
//Point req[3];
void LoadLog(char * fName, list<list<Point>> & pList);
DWORD WINAPI ThreadProc( LPVOID lpParameter );
void LoadXDATCAR(list<list<Point>> & pList);
CRITICAL_SECTION CriticalSection; 

Matrix FractToDec, DecToFract;

vector<string> Labels;
vector<string> labbuf;
Param param;
int main(int argn, char * argv[]) {
	vector<string> OUT_st;
	param.TakeAgrs(argn,argv);
	param.ReadNoname(OUT_st);
	int NThread = atoi(OUT_st[0].data());
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
	vector<list<Elipsoid>> ElList;
	if (!InitializeCriticalSectionAndSpinCount(&CriticalSection, 0x00000400) ) 
		return -1; 
	list<list<Point>> points;
	//LoadLog("opt2.log",points);
	LoadXDATCAR(points);
	ElList.resize(NThread);
	auto iter = points.begin();
	int j = 0;
	while(iter != points.end()) {
		for(int i = 0; i<NThread; i++) {
			ElList[i].push_back(Elipsoid(*iter,Labels[j],j+1));
			iter++;
			j++;
			if(iter == points.end())
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
	for(int i=0; i<NThread-1; i++)
    {
        CloseHandle(hThreadArray[i]);
    }
    DeleteCriticalSection(&CriticalSection);
	//system("pause");
	return 0;
}
void Load(char * fName, list<Point> & points)
{
	
	ifstream file(fName);
	if (file.is_open()) {
		cout << "File opened." << endl;
	} else {
		cout << "Error! File not opened." << endl;
	}
	char buf[MAX_LINE];
	//auto iter = points.begin();
	while (!file.eof())
	{
		Point p;
		file>> buf;
		p.x = atof(buf);
		file>> buf;
		p.y = atof(buf);//!!!
		file>> buf;
		p.z = atof(buf);//!!!
		points.push_back(p);
	}
}
void LoadLog(char * fName, list<list<Point>> & pList)
{
	char strNAtoms[] = "NAtoms=";
	int NAtoms = 0 ;
	char strStop[] = "                          Input orientation:                          ";
	char ZMatrix[] = "                            Z-MATRIX (ANGSTROMS AND DEGREES)";
	ifstream file(fName);
	if (file.is_open()) {
		cout << "File opened." << endl;
	} else {
		cout << "Error! File not opened." << endl;
	}

	char buf[MAX_LINE];
	while (!file.eof()) {
		file >> buf;
		if(strcmp(buf,strNAtoms)==0) {
			file >> buf;
			NAtoms = atoi(buf);
			break;
		}
	}
	if(NAtoms==0)
		cout << "Error! No info about NAtoms." << endl;
	else
		cout << "NAtoms = " << NAtoms << endl;
	

	//auto iter = points.begin();
	pList.resize(NAtoms);
	while (!file.eof()) {
		file.getline(buf,MAX_LINE);
		if(strcmp(buf,ZMatrix)==0) {
			break;
		}
	}
	file.getline(buf,MAX_LINE);
	file.getline(buf,MAX_LINE);
	for(int i = 0; i<NAtoms; i++) {
		file >> buf >>buf>> buf;
		Labels.push_back(buf);
		file.getline(buf,MAX_LINE);
	}
	Labels.shrink_to_fit();

	while (!file.eof())
	{
		while (!file.eof()) {
			file.getline(buf,MAX_LINE);
			if(strcmp(buf,strStop)==0) {
				break;
			}
		}
		if(file.eof()) {
			cout << "Log loaded." << endl;
			return;
		}
		file.getline(buf,MAX_LINE);
		file.getline(buf,MAX_LINE);
		file.getline(buf,MAX_LINE);
		file.getline(buf,MAX_LINE);
		auto iter = pList.begin();
		for(int i = 0 ; i<NAtoms;i++) {
			file >> buf >> buf >> buf;
			Point p;
			file>> buf;
			p.x = atof(buf);
			file>> buf;
			p.y = atof(buf);
			file>> buf;
			p.z = atof(buf);
			iter->push_back(p);
			iter++;
		}
		//Point p;
		//file>> buf;
		//p.x = atof(buf);
		//file>> buf;
		//p.y = atof(buf);//!!!
		//file>> buf;
		//p.z = atof(buf);//!!!
		//points.push_back(p);
	}
}
void LoadXDATCAR(list<list<Point>> & pList)
{
	int NAtoms = 0 ;
	char * fName = "XDATCAR";
	ifstream file(fName);
	if (file.is_open()) {
		cout << "XDATCAR opened." << endl;
	} else {
		cout << "Error! XDATCAR not opened." << endl;
		system("pause");
		exit(10);
	}
	char buf[MAX_LINE];
	
	file.getline(buf,MAX_LINE);
	file.getline(buf,MAX_LINE);
	
	file >> FractToDec.U[0] >> FractToDec.U[3] >> FractToDec.U[4];
	file >> FractToDec.U[1] >> FractToDec.U[1] >> FractToDec.U[5];
	file >> FractToDec.U[2] >> FractToDec.U[2] >> FractToDec.U[2];
	DecToFract = FractToDec.Mirror();
	Cell cell(FractToDec);
	
	file.getline(buf,MAX_LINE);

	
	file.getline(buf,MAX_LINE);
	std::stringstream str(buf);
	while(!str.eof()) {
		string temp;
		str >> temp;
		labbuf.push_back(temp);

	}
	labbuf.pop_back();
	labbuf.shrink_to_fit();
	int size = labbuf.size();
	vector<int> counts;
	counts.resize(size);
	for(int i = 0; i<size; i++)
		file >> counts[i];
	file.getline(buf,MAX_LINE);
	{
		ofstream temp("a.ins",ios::trunc);
		temp << "CELL 0.71073 " << cell.a << " " << cell.b << " " << cell.c << " " 
			<< cell.alpha << " " << cell.beta << " " << cell.gamma << " " << endl;
		temp << "ZERR 4 0.015 0.015 0.007 0 0 0" << endl;
		//temp << "LATT 2" << endl;
		temp << "SFAC";
		for(int i = 0; i<size;i++)
			temp << " " << labbuf[i];
		temp << endl << "UNIT";
		for(int i = 0; i<size;i++)
			temp << " " << counts[i];
		temp <<  endl;
		temp <<  endl;
		temp.close();
	}

	//Point t = FractToDec.Transform(Point(15,4,6));

	for(int i = 0; i <size; i++) {
		NAtoms+=counts[i];
		Labels.reserve(NAtoms);
		for(int j = 0; j < counts[i]; j++) {
			Labels.push_back(labbuf[i]);
		}
	}
	Labels.shrink_to_fit();
	pList.resize(NAtoms);
	//file.getline(buf,MAX_LINE);
	bool reverse = false;
	Point lastcenter(0,0,0);
	int Counter = 0;
	int AllSteps = 0;
	while(!file.eof()) {
		list<Point> tempList;
		Point center;
		file.getline(buf,MAX_LINE);
		auto iter = tempList.begin();
		
		for(int i =0 ; i<NAtoms;i++) {
			Point p;
			file>> buf;
			p.x = atof(buf);
			file>> buf;
			p.y = atof(buf);
			file>> buf;
			p.z = atof(buf);
			tempList.push_back(p);
			center=p+center;
			//iter++;
		}
		Elipsoid El(tempList,string("H"),1);
		El.ShiftCenter(lastcenter);
		//El.req[1] = *(++El.points.begin());
		//El.req[2] = *(++++(El.points.begin()));
		//El.ReturnToStartPos();

		//center.x/=NAtoms;
		//center.y/=NAtoms;
		//center.z/=NAtoms;
		//center=center-Point(0.5,0.5,0.5);


		//tempList.swap(El.points);


		auto iter2 = pList.begin();
		iter = El.points.begin();
		for(int i =0 ; i<NAtoms;i++) {
			if(reverse) {
				if(abs(iter2->begin()->x - iter->x) > 0.5)
					iter->x+=_sign(iter2->begin()->x - iter->x);
				if(abs(iter2->begin()->y - iter->y) > 0.5)
					iter->y+=_sign(iter2->begin()->y - iter->y);
				if(abs(iter2->begin()->z - iter->z) > 0.5)
					iter->z+=_sign(iter2->begin()->z - iter->z);
			}
			iter2->push_back(*iter);
			iter++;
			iter2++;
		}
		center = El.SearchCenter()-Point(0.5,0.5,0.5);
		El.ShiftCenter(center);
		tempList.swap(El.points);
		lastcenter = lastcenter+center;

		iter2 = pList.begin();
		iter = tempList.begin();
		for(int i =0 ; i<NAtoms;i++) {
			iter2->push_front(*iter);
			iter++;
			iter2++;
		}
		
		reverse = true;

		if(true);
		file.getline(buf,MAX_LINE);
		Counter++;
		if(Counter == 100) {
			Counter=0;
			AllSteps+=100;
			cout << "Step " << AllSteps << " readed." <<endl;
		}
	}
	cout << "Last step " << AllSteps+Counter-1 << " readed." <<endl;
	auto iter = pList.begin();
	for(int i =0 ; i<NAtoms;i++) {
		iter->pop_front();
		iter++;
	}
}
flo _quad(flo a) 
{
	return a*a;
}
int _sign(flo a)
{
	return a>0?1:-1;
}
void Out(Elipsoid & el)
{
	ofstream file("a.txt",ios::app);
	if (file.is_open()) {
		//cout << "File opened." << endl;
	} else {
		cout << "Error! File not opened." << endl;
	}
	//auto iter = points.begin();
	//file << "%chk=d4.chk" << endl
	//	<< "%mem=16GB" <<endl
	//	<<"%nprocshared=8" <<endl
	//	<< "#p pbe1pbe/6-311g ADMP( MaxPoints=100)" <<endl
	//	<<endl
	//	<< "d4" <<endl
	//	<<endl
	//	<<"0 1"<<endl;
	auto iter = el.points.begin();
	int i=0;
	while(iter!= el.points.end()) {
		if(i!=0) {
			iter++;
			i++;
			i=(i)%10;
			continue;

		}
		Point & p = *iter;
		file /*<< "Si "*/ << iter->x<< " " <<
			iter->y<< " " <<
			iter->z<< endl;
		iter++;
		i++;
		i=(i)%10;
	}
	file.close();
	//for(int i=0; i<18; i++) {

	//	file /*<< "O "*/ << el._pts[i].x*10+0.00000001 << " " <<
	//		el._pts[i].y*10+0.00000001 << " " <<
	//		el._pts[i].z*10+0.00000001 << endl;
	//}
}
void OutIns(Elipsoid & el)
{
	ofstream file("a.ins",ios::app);
	if (file.is_open()) {
		//cout << "File opened." << endl;
	} else {
		cout << "Error! File not opened." << endl;
	}
	//auto iter = points.begin();
	//file << "%chk=d4.chk" << endl
	//	<< "%mem=16GB" <<endl
	//	<<"%nprocshared=8" <<endl
	//	<< "#p pbe1pbe/6-311g ADMP( MaxPoints=100)" <<endl
	//	<<endl
	//	<< "d4" <<endl
	//	<<endl
	//	<<"0 1"<<endl;
	auto iter = el.points.begin();
	file << fixed;

	file.precision(5);
	//while(iter!= el.points.end()) {
	//	Point p = *iter;
	file << el.label<< el.number<< " ";
	int atype = 0;
	for(;atype < labbuf.size();atype++)
		if(labbuf[atype].compare(el.label)==0) break;
	file << atype+1 <<" " << el.center.x<< " " <<
			el.center.y<< " " <<
			el.center.z << " 11.00000 " <<
			el.mat.U[0] << " " <<
			el.mat.U[1] << " " <<
			el.mat.U[2] << " =" << endl << " " <<
			el.mat.U[5] << " " <<
			el.mat.U[4] << " " <<
			el.mat.U[3] << " " <<endl;


	file.close();

}
DWORD WINAPI ThreadProc( LPVOID lpParameter )
{
	list<Elipsoid> * El = reinterpret_cast<list<Elipsoid> *>(lpParameter);
	auto iter = El->begin();
	while(iter!= El->end()) {
		//iter->Correlation();
		//iter->Process();
		iter->UProc();
		iter++;
	}
	iter = El->begin();

    EnterCriticalSection(&CriticalSection); 
	while(iter!= El->end()) {
		OutIns(*iter);
		//Out(*iter);
		iter++;
	}
    LeaveCriticalSection(&CriticalSection);

return 0;
}
