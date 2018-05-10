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
	cout << "Program Ellipsoid. Version 1.2.0\n" << endl;
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
			cerr << "Cannot open <Shelx_File>. Continue without symmetry." << endl;
		}
		old.close();
	}
	vector<vector<Point> > El;
	try {
		nsShelxFile::ShelxData sheltemp(nsShelxFile::XDATCAR);
		shelx.cell = move(sheltemp.cell);
		El = nsShelxFile::ShelxData::LoadXDATCAR(cutoff, &fPos);
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
		cout << "Symmetry analise started." << endl;
		try {
			Analize_symmety(shelx, El);
		}
		catch (invalid_argument & ex) {
			cerr << ex.what() << endl;
			return 1;
		}
		catch (...) {
			cerr << "Unknown error during symmetry analising." << endl;
			return 1;
		}
		cout << "Symmetry analise complited." << endl;
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
	auto eqMat = Matrix::EqualMatrix(3);
	vector<Matrix> Tables(size_el, eqMat);
	vector<Point> Shift(size_el);
	vector<int> To_n(size_el, -1);
	//vector<bool> Used(size_el, false);
	constexpr int _p = 1;
	constexpr size_t _d = 2*_p+1;
	constexpr size_t sizemod = (_d*_d*_d);
	size_t size_b = size_el*sizemod;
	vector<Point> basis(size_b);


	// Create basis
	for (int j = -_p, iter = 0; j <= _p; j++) {
		for (int k = -_p; k <= _p; k++) {
			for (int l = -_p; l <= _p; l++) {
				for (size_t i = 0; i < size_el; i++,iter++) {
					basis[iter] =( fPos[i] + Point(j, k, l));
				}
			}
		}
	}
	// Create first atoms
	for (size_t i = 0; i < size_shelx_atom; i++)
	{
		size_t n = size_el;
		flo d = 2;	
		for (size_t j = 0; j < size_b; j++)
		{
			flo nd = (shelx.cell.FracToCart() * (basis[j] - shelx.atom[i].point)).r();
			if (nd < d) {
				d = nd;
				n = j % size_el;
			}
		}
		if (n == size_el)
			throw invalid_argument("Bad shelx file.");
		To_n[n] = n;
		for (size_t p = 0; p < 3; p++)
		{
			flo s = shelx.atom[i].point.a[p] - fPos[n].a[p];
			if (s > static_cast<flo>(0.5)) Shift[n].a[p]+=1;
			else if (s <= static_cast<flo>(-0.5)) Shift[n].a[p] -= 1;
		}
	}
	// T = sizemod * size_s * size_el^2
	{
		vector<int> tom;
		tom.reserve(size_el);
		for (size_t i = 0; i < size_el; i++)
		{
			if (To_n[i] != -1) {
				tom.push_back(i);
			}
		}
		vector<Matrix> Invsymm;
		Invsymm.reserve(size_s);
		for (size_t s = 0; s < size_s; s++) {
			Invsymm.push_back(shelx.symm[s].mat.Invert());
		}
		for (size_t i = 0; i < size_el; i++) {
			if (i >= tom.size()) {
				throw invalid_argument("Need more atoms in shelx file.");
			}
			for (size_t s = 0; s < size_s; s++) {
				Point ps = shelx.symm[s].GenSymmNorm(fPos[tom[i]]);
				size_t n = size_el;
				flo d = 2;
				for (size_t j = 0; j < size_b; j++)
				{
					flo nd = (shelx.cell.FracToCart() * (basis[j] - ps)).r();
					if (nd < d) {
						d = nd;
						n = j % size_el;
					}
				}
				if (n == size_el)
					throw invalid_argument("Bad symmetry in shelx file.");
				if (To_n[n] != -1) continue;
				To_n[n] = To_n[tom[i]];
				Tables[n] = Invsymm[s] * Tables[tom[i]];
				Shift[n] = (Tables[n].Invert() * fPos[To_n[n]]) - fPos[n];
				tom.push_back(n);
				const size_t nsize = pList[n].size();
				pList[To_n[n]].reserve(pList[To_n[n]].size() + nsize);
				for (size_t l = 0; l < nsize; l++)
				{
					pList[To_n[n]].push_back(Tables[n] * (Shift[n] + pList[n][l]));
				}
				pList[n].clear();
			}
		}
	}
	{
		vector<vector<Point> > temp;
		temp.reserve(size_shelx_atom);
		for (size_t j = 0; j < size_shelx_atom; j++) {
			for (size_t i = 0; i < size_el; i++) {
				if (!(pList[i].empty())) {
					if(Point(fmod(fPos[i].a[0] - shelx.atom[j].point.a[0] + _d, 1),
						fmod(fPos[i].a[0] - shelx.atom[j].point.a[0] + _d, 1),
						fmod(fPos[i].a[0] - shelx.atom[j].point.a[0] + _d, 1)).r() > 0.01) continue;
					temp.push_back(move(pList[i]));
				}
			}
		}
		pList.swap(temp);
		if (size_shelx_atom != pList.size()) {
			throw invalid_argument("Bad shelx file.");
		}
	}
}
