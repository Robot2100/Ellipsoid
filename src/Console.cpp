#include "stdafx.h"


using namespace std;
constexpr int MAX_LINE = 128;
size_t cutoff = 2000;
string filenamein;

std::vector<Point> fPos;
int _sign(flo a);
void Analize_symmety(nsShelxFile::ShelxData & shelx, vector<vector<Point> > & pList);
void ffunc(const int l, std::vector<string> & in) {
	bool help = false;
	auto size = in.size();
	switch(l) {
	case 0:
		if (size == 1)
			filenamein = std::move(in[0]);
		else
			throw invalid_argument("Wrong number of parameters.");
		break;
	case 1:
		if (size == 1)
			cutoff = static_cast<size_t>(stoi(in[0]));
		else
			throw invalid_argument("Wrong number of parameters of '-cut'.");
		break;
	}
}




int main(int argn, char * argv[]) {
	ios::sync_with_stdio(false);
	{
		constexpr BaseParam bp[] {	
			{"",	"",		"<Filename>",	"Take symmetry from shelx file [optional]"},
			{"c",	"cut",	"<N>",			"Ignore first N steps [default=2000]" } };
		constexpr Param<2> param(bp);
		try {
			param.TakeAgrs(argn, argv, ffunc);
		}
		catch (invalid_argument & inv) {
			cout << "Error! Program termination. Reason:\n" << inv.what() << endl;
			return 1;
		}
		catch (IncExceptions::ParamException & inv) {
			cout << "Error! Unknown parameter: " << inv.what()
				<< "\nUse -h or --help parameter for more information." << endl;
			return 1;
		}
	}
	cout << "Program Ellipsoid. Version 1.1.1\n" << endl;
	cout << "Ignore first " << cutoff << " steps." << endl;
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
	vector<vector<Point> > El = shelx.LoadXDATCAR(cutoff, &fPos);


	if (is_SYMM == true) {
		Analize_symmety(shelx, El);
	}
	cout << "Symmetry analise complited." << endl;
	size_t Elsize = El.size();
	{
		vector<nsShelxFile::Atom> atombuf;
		for (size_t i = 0; i < Elsize; i++)
		{
			atombuf.push_back(nsShelxFile::Atom(shelx.atom[i].label, shelx.atom[i].type, shelx.atom[i].occup, shelx.cell, El[i], true));
		}
		shelx.atom = move(atombuf);
	}
	cout << "Writing to file 'a.ins'." << endl;
	ofstream out("a.ins");
	shelx.OutIns(out);
	cout << "Program normal termination. " << endl;
	return 0;
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
				for (size_t i = 0; i < size_el; i++,iter++) {
					basis[iter] = shelx.cell.FracToCart() *( fPos[i] + Point(j, k, l));
				}
			}
		}
	}

	vector<nsShelxFile::SYMM> mirror;
	for (size_t i = 0; i < size_s; i++) {
		mirror.push_back((shelx.symm[i].MirrorSymm()));
	}

	for (size_t s = 0; s < size_s; s++) {
		vector<Point> pvec(size_el);
		for (size_t i = 0; i < size_el; i++) {
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
		for (size_t i = 0; i < size_el; i++) {
			if (table[i].size() == 0) continue;
			for (size_t j = 0; j < size_el; j++) {
				if (table[i][j].size() == 0 || pList[j].size() == 0) continue;
				size_t k = pList[i].size();
				pList[i].reserve(k + pList[j].size());
				pList[i].insert(pList[i].end(),
					std::make_move_iterator(pList[j].begin()),
					std::make_move_iterator(pList[j].end()));
				pList[j].clear();
				size_t it = pList[i].size();

				for (size_t k2 = 0; k2 < table[i][j].size(); k2++) {
					pList[i][k] = (table[i][j][k2]->RetroGenSymm(pList[i][k]));
				}
				Point dfix = pList[i][k] + Point(100.5, 100.5, 100.5) - fPos[i];
				dfix = (Point(100, 100, 100) - Point(int(dfix.a[0]), int(dfix.a[1]), int(dfix.a[2])));
				pList[i][k] += dfix;

				for (++k; k < it; k++) {
					for (size_t k2 = 0; k2 < table[i][j].size(); k2++) {
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
	for (size_t i = 0; i < size_el; i++) {
		if (pList[i].empty() == false) {
			temp.push_back(move(pList[i]));
			tempAtom.push_back(move(shelx.atom[i]));
		}
	}
	pList = move(temp);
	shelx.atom = move(tempAtom);
}