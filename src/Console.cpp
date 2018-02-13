#include "stdafx.h"
#include "../Includes/Includes.h"


using namespace std;
constexpr int MAX_LINE = 128;
constexpr size_t cutoff = 2000;
string filenamein;

std::vector<Point> fPos;
void LoadXDATCAR(nsShelxFile::ShelxData & shelx, vector<vector<Point> > & pList);
int _sign(flo a);
Point TakePoint(istream & file, const size_t NAtoms);
void Analize_symmety(nsShelxFile::ShelxData & shelx, vector<vector<Point> > & pList);
void ffunc(const int l, std::vector<string> & in) {
	bool help = false;
	auto size = in.size();
	switch(l) {
	case 0:
		if (size != 0)
			filenamein = std::move(in[0]);
		break;
	case 1:
		help = true;
		break;
	case -1: 
		if (help) {
			for (int i = 0; i < size; i++)
				cout << in[i] << endl;
			exit(0);
		}
	}
}
int main(int argn, char * argv[]) {
	{
		constexpr BaseParam bp[2]{ {"","Input shelx file (optional)"},
									{"-help", "View help information"} };
		constexpr ConstParam<2> cp(bp);
		Param<2> param(&cp);
		try {
			param.TakeAgrs(argn, argv, ffunc);
		}
		catch (invalid_argument inv) {
			cout << "Error! Program termination. Reason:" << endl << inv.what() << endl;
			return 1;
		}
	}
	bool is_SYMM = false;
	nsShelxFile::ShelxData shelx;
	if (filenamein.length() != 0) {
		cout << "Inputed <Shelx_File> name: " << filenamein << endl;
		ifstream old(filenamein);
		if (old.is_open() == true) {
			shelx = nsShelxFile::ShelxData(old);
			shelx.atom.clear();
			cout << "Symmetry found." << endl;
			is_SYMM = true;
		}
		else {
			cout << "Cannot open <Shelx_File>. Continue without symmetry." << endl;
		}
		old.close();
	}
	vector<vector<Point> > El;
	LoadXDATCAR(shelx,El);

	if (is_SYMM == true) {
		Analize_symmety(shelx, El);
	}
	size_t Elsize = El.size();
	{
		vector<nsShelxFile::Atom> atombuf;
		for (size_t i = 0; i < Elsize; i++)
		{
			atombuf.push_back(nsShelxFile::Atom(shelx.atom[i].label, shelx.atom[i].type, shelx.atom[i].occup, shelx.cell, El[i], true));
		}
		shelx.atom = move(atombuf);
	}
	ofstream out("a.ins");
	shelx.OutIns(out);
	return 0;
}
void LoadXDATCAR(nsShelxFile::ShelxData & shelx, vector<vector<Point> > & pList)
{
	shelx.sfac.clear();
	size_t NAtoms = 0;
	constexpr char fName[] = "XDATCAR";
	ifstream file(fName);
	if (file.is_open()) {
		cout << "XDATCAR opened." << endl;
	}
	else {
		cout << "Error! XDATCAR not opened." << endl;
		exit(10);
	}
	char buf[MAX_LINE];

	file.getline(buf, MAX_LINE);
	file.getline(buf, MAX_LINE);
	flo U[9];
	
	file >> U[0] >> U[1] >> U[2];
	file >> U[3] >> U[4] >> U[5];
	file >> U[6] >> U[7] >> U[8];
	shelx.cell = Cell(Matrix(&U[0], 3, 3), true);

	file.getline(buf, MAX_LINE);


	file.getline(buf, MAX_LINE);
	std::stringstream str(buf);
	while (!str.eof()) {
		string temp;
		str >> temp;
		shelx.sfac.push_back(move(temp));
	}
	shelx.sfac.pop_back();
	shelx.sfac.shrink_to_fit();
	size_t size = shelx.sfac.size();

	shelx.unit.reserve(size);
	shelx.unit.resize(size);

	for (size_t i = 0; i < size; i++) {
		file >> buf;
		shelx.unit[i] = atof(buf);
	}
	file.getline(buf, MAX_LINE);

	for (size_t i = 0; i < size; i++) {
		NAtoms += shelx.unit[i];
	}
	pList.resize(NAtoms);
	for (int i = 0, k = 0; i < size; i++) {
		for (int j = 1; j <= shelx.unit[i]; j++, k++) {
			char str[128];
			sprintf(str, "%s%d", shelx.sfac[i].c_str(), j);
			shelx.atom.push_back(nsShelxFile::Atom(str, i + 1, Point(), flo(1.0), Dinmat()));
		}
	}
	file.getline(buf, MAX_LINE);
	fPos.reserve(NAtoms);
	for (int i = 0; i<NAtoms; i++) {
		Point p(TakePoint(file, NAtoms));
		pList[i].push_back(p);
		fPos.push_back(p);
	}

	int Counter = 1;
	int AllSteps = 1;
	for (int k = 1; k < cutoff && !file.eof(); k++) {
		file.getline(buf,MAX_LINE);
		for (size_t i = 0; i < NAtoms; i++) {
			file.getline(buf, MAX_LINE);
		}
		Counter++;
		if (Counter % 100 == 0) {
			cout  << Counter << " steps skipped." << endl;
		}
	}
	while (!file.eof()) {
		file.getline(buf, MAX_LINE);

		for (int i = 0; i<NAtoms; i++) {
			pList[i].push_back(TakePoint(file, NAtoms));
		}
		for (int i = 0; i<NAtoms; i++) {
			Point dp = pList[i][AllSteps] - pList[i][AllSteps - 1];
			for (int j = 0; j < 3; j++) {
				if (abs(dp.a[j]) > 0.5)
					pList[i][AllSteps].a[j] -= _sign(dp.a[j]);
				else dp.a[j] = 0;
			}
		}
		AllSteps++;
		Counter++;
		if (Counter%100 == 0) {
			cout << Counter << " steps readed. Step " << AllSteps-1 << " uses." << endl;
		}
	}
	cout << "Last step " << Counter-1 << " readed." << endl;

	for (int i = 0; i < NAtoms; i++) {
		pList[i].pop_back();
		pList[i].shrink_to_fit();
	}

}

int _sign(flo a) {
	if (a < 0) return -1;
	else return 1;
}

void Analize_symmety(nsShelxFile::ShelxData & shelx, vector<vector<Point> > & pList) {
	using namespace nsShelxFile;
	size_t size_s = shelx.symm.size();
	size_t size_el = pList.size();
	vector<vector<vector<SYMM*> > > table(size_el, vector<vector<SYMM*> >(size_el));
	vector<vector<Point> > addtable(size_el, vector<Point>(size_el));
	constexpr int _p = 1; 
	constexpr size_t _d = 2*_p+1;
	constexpr size_t _sizemod = (_d*_d*_d);
	size_t size_b = size_el*_sizemod;
	vector<Point> basis(size_b);


	for (int j = -_p, iter = 0; j <= _p; j++) {
		for (int k = -_p; k <= _p; k++) {
			for (int l = -_p; l <= _p; l++) {
				for (int i = 0; i < size_el; i++,iter++) {
					basis[iter] = shelx.cell.FracToCart() *( fPos[i] + Point(j, k, l));
				}
			}
		}
	}
	vector<nsShelxFile::SYMM> mirror;
	for (int i = 0; i < size_s; i++) {
		mirror.push_back((shelx.symm[i].MirrorSymm()));
	}



	for (int s = 0; s < size_s; s++) {
		vector<Point> pvec(size_el);
		for (int i = 0; i < size_el; i++) {
			pvec[i] = shelx.cell.FracToCart() * (shelx.symm[s].GenSymm(fPos[i]));
		}
		for (size_t i = 0; i < size_el; i++) {
			flo drm = (flo)1.0;
			size_t j1 = 0;
			for (size_t j = 0; j < size_b; j++) {
				flo dr = (pvec[i] - basis[j]).r();
				if (dr < drm) {
					drm = dr;
					j1 = j;
				}
			}
			size_t j2 = j1%size_el;
			if (i == j2 || table[i][j2].size() != 0) 
				continue;
			table[i][j2].push_back(&mirror[s]);
			addtable[i][j2] = shelx.cell.CartToFrac()*basis[j1] - fPos[j2];
		}
	}
	{




		bool changed = false;
		do {
			changed = false;
			for (size_t i = 0; i < size_el; i++) {
				if (table[i].size() == 0) continue;
				for (size_t j = 0; j < size_el; j++) {
					if (table[i][j].size() == 0 || i==j || table[j].size() == 0) continue;
					for (size_t j1 = 0; j1 < size_el; j1++) {
						if (table[j][j1].size() == 0 || table[i][j1].size() != 0 || i==j1) continue;
						table[i][j1].reserve(table[i][j].size() + table[j][j1].size());
						table[i][j1].insert(table[i][j1].end(), table[j][j1].begin(), table[j][j1].end());
						table[i][j1].insert(table[i][j1].end(), table[i][j].begin(), table[i][j].end());
						addtable[i][j1] = addtable[i][j] + addtable[j][j1];
					}
					table[j].clear();
					changed = true;
				}
			}
		} while (changed == true);
	}




	vector<Point> dcheck;
	{
		size_t size_l = pList[0].size();
		for (int i = 0; i < size_el; i++) {
			if (table[i].size() == 0) continue;
			for (int j = 0; j < size_el; j++) {
				if (table[i][j].size() == 0 || pList[j].size() == 0) continue;
				size_t k = pList[i].size();
				pList[i].reserve(k + pList[j].size());
				pList[i].insert(pList[i].end(),
					std::make_move_iterator(pList[j].begin()),
					std::make_move_iterator(pList[j].end()));
				pList[j].clear();
				size_t it = pList[i].size();

				for (int k2 = 0; k2 < table[i][j].size(); k2++) {
					pList[i][k] = (table[i][j][k2]->RetroGenSymm(pList[i][k]));
				}
				Point dfix = pList[i][k] + Point(100.5, 100.5, 100.5) - fPos[i];
				dfix = (Point(100, 100, 100) - Point(int(dfix.a[0]), int(dfix.a[1]), int(dfix.a[2])));
				pList[i][k] += dfix;

				for (++k; k < it; k++) {
					for (int k2 = 0; k2 < table[i][j].size(); k2++) {
						pList[i][k] = (table[i][j][k2]->RetroGenSymm(pList[i][k]));
					}
					pList[i][k] += dfix;
				}
				dcheck.push_back(pList[i].back());
			}
		}
	}
	std::remove_reference<decltype(pList)>::type temp;
	decltype(shelx.atom) tempAtom;
	for (int i = 0; i < size_el; i++) {
		if (pList[i].empty() == false) {
			temp.push_back(move(pList[i]));
			tempAtom.push_back(move(shelx.atom[i]));
		}
	}
	pList = move(temp);
	shelx.atom = move(tempAtom);
}
Point TakePoint(istream & file, const size_t NAtoms) {
	char buf[MAX_LINE];
	file.getline(buf, (MAX_LINE-1));
	char * end = buf;
	Point out;
	out.a[0] = strtod(end, &end);
	out.a[1] = strtod(end, &end);
	out.a[2] = strtod(end, NULL);
	return out;
}
