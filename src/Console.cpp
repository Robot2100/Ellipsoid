#include "stdafx.h"


using namespace std;

constexpr int MAX_LINE = 128;
size_t cutoff = 2000;
string filenamein;

std::vector<Point> fPos;
void Analize_symmety(nsShelxFile::ShelxData & shelx, vector<vector<Point> > & pList);
void ffunc(const int l, std::vector<string> & in) {
	bool help = false;
	auto size = in.size();
	switch(l) {
	case 0:
		if (size == 1)
			filenamein = std::move(in[0]);
		else
			throw invalid_argument("Wrong number of parameters. Must be equal 1.");
		break;
	case 1:
		if (size == 1)
			cutoff = static_cast<size_t>(stoi(in[0]));
		else
			throw invalid_argument("Wrong number of parameters of '-c' or '--cut'. Must be equal 1.");
		break;
	case 2:
		if (size == 0)
			cout.rdbuf(NULL);
		else
			throw invalid_argument("Wrong number of parameters of '-q' or '--quiet'. Must be equal 0.");
		break;
	}
}

int main(int argn, char * argv[]) {
	ios::sync_with_stdio(false);
	{
		constexpr BaseParam bp[]{
			{"",	"",		"<Filename>",	"Take symmetry from shelx file [optional]"},
			{ "c",	"cut",	"<N>",			"Ignore first N steps [default=2000]" },
			{ "q", "quiet", "", "Output only error messages" } };
		constexpr Param<3> param(bp);
		try {
			param.TakeAgrs(argn, argv, ffunc);
		}
		catch (invalid_argument & inv) {
			cerr << "Error! Program termination. Reason:\n" << inv.what() << endl;
			return 1;
		}
		catch (IncExceptions::ParamException & inv) {
			cerr << "Error! Unknown parameter: " << inv.what()
				<< "\nUse -h or --help parameter for more information." << endl;
			return 1;
		}
		catch (...) {
			cerr << "Unknown error during parsing parameters." << endl;
			return 1;
		}
	}
	cout << "Program Ellipsoid. Version 1.1.2\n" << endl;
	cout << "Ignore first " << cutoff << " steps." << endl;
	bool is_SYMM = false;
	nsShelxFile::ShelxData shelx;
	if (filenamein.length() != 0) {
		cout << "Inputed <Shelx_File> name: " << filenamein << endl;
		ifstream old(filenamein);
		if (old.is_open() == true) {
			shelx = nsShelxFile::ShelxData(old);
			cout << "Symmetry found." << endl;
			is_SYMM = true;
		}
		else {
			cout << "Cannot open <Shelx_File>. Continue without symmetry." << endl;
		}
		old.close();
	}
	vector<vector<Point> > El;
	try {
		El = shelx.LoadXDATCAR(cutoff, &fPos);
	}
	catch (IncExceptions::OpenXDATCAR_Exception & ex) {
		cerr << ex.what() << endl;
		return 1;
	}
	catch (IncExceptions::ReadXDATCAR_Exception & ex) {
		cerr << ex.what() << endl;
		return 1;
	}
	catch (...) {
		cerr << "Unknown error during loading XDATCAR file." << endl;
		return 1;
	}

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
void Analize_symmety(nsShelxFile::ShelxData & shelx, vector<vector<Point> > & pList) {
	using namespace nsShelxFile;
	size_t size_s = shelx.symm.size();
	size_t size_el = pList.size();
	auto size_shelx_atom = shelx.atom.size();
	
	vector<Matrix> tables(Matrix({1,0,0,0,1,0,0,0,1},3,3);
	constexpr int _p = 1; 
	constexpr size_t _d = 2*_p+1;
	constexpr size_t sizemod = (_d*_d*_d);
	size_t size_b = size_el*sizemod;
	vector<Point> basis(size_b);
	vector<vector<size_t> > AtomList(size_shelx_atom);


	// Create basis
	for (int j = -_p, iter = 0; j <= _p; j++) {
		for (int k = -_p; k <= _p; k++) {
			for (int l = -_p; l <= _p; l++) {
				for (size_t i = 0; i < size_el; i++,iter++) {
					basis[iter] = shelx.cell.FracToCart() * ( fPos[i] + Point(j, k, l));
				}
			}
		}
	}
	// Create AtomList[][0]
	for (size_t i = 0; i < size_shelx_atom; i++)
	{
		size_t n = size_el;
		flo d = 1.0;	
		for (size_t j = 0; j < size_el; j++)
		{
			flo nd = (basis[j] - shelx.atom[i].point).r();
			if (nd < d) {
				d = nd;
				n = j / sizemod;
			}
		}
		if (n == size_el)
			throw invalid_argument("Bad shelx file.");
		AtomList[i].push_back(n);
	}
	// Create AtomList[][i]
	for (size_t s = 0; s < size_s; s++) {

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

	{
		vector<Point> dcheck;
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
Matrix EqualMatrix(const size_t n) {
	using namespace std;
	vector<vector<flo> > vec(n, vector<flo>(n, static_cast<flo>(0.0)));
	for (size_t i = 0; i < n; i++)
	{
		vec[i][i] = static_cast<flo>(1.0);
	}
	return Matrix(std::move(vec));
}
