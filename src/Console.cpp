#include "stdafx.h"


using namespace std;
constexpr int MAX_LINE = 128;
size_t cutoff = 2000;
string filenamein;

std::vector<Point> fPos;
void Analize_symmety(nsShelxFile::ShelxData & shelx, nsShelxFile::ShelxData & xdat, vector<vector<Point> > & pList);
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
	cout << "Program Ellipsoid. Version 1.2.0\n" << endl;
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
	cout << "Ignore first " << cutoff << " steps." << endl;
	bool is_SYMM = false;
	nsShelxFile::ShelxData shelx;
	nsShelxFile::ShelxData xdata;
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
			is_SYMM = false;
		}
		old.close();
	}
	vector<vector<Point> > El;
	try {
		xdata = move(nsShelxFile::ShelxData(nsShelxFile::XDATCAR));
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
			Analize_symmety(shelx, xdata, El);
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
	else {
		size_t Elsize = El.size();
		vector<nsShelxFile::Atom> atombuf;
		for (size_t i = 0; i < Elsize; i++)
		{
			atombuf.push_back(nsShelxFile::Atom(xdata.atom[i].label, xdata.atom[i].type, 1.0, xdata.cell, move(El[i]), true));
		}
		xdata.atom = move(atombuf);
	}
	cout << "Writing to file 'a.ins'." << endl;
	ofstream out("a.ins");
	xdata.OutIns(out);
	cout << "Program normal termination. " << endl;
	return 0;
}
void Analize_symmety(nsShelxFile::ShelxData & shelx, nsShelxFile::ShelxData & xdat, vector<vector<Point> > & pList) {
	using namespace nsShelxFile;
	const size_t size_s = shelx.symm.size();
	const size_t size_el = pList.size();
	const auto size_shelx_atom = shelx.atom.size();
	const auto eqMat = Matrix::EqualMatrix(3);
	vector<Matrix> Tables(size_el, eqMat);
	vector<Point> Shift(size_el), checkPoint(size_el);
	vector<int> To_n(size_el, -1);
	constexpr int _p = 1;
	constexpr size_t _d = 2*_p+1;
	constexpr size_t sizemod = (_d*_d*_d);
	const size_t size_b = size_el*sizemod;
	vector<Point> basis(size_b);
	vector<size_t> from_to(xdat.sfac.size());
	auto Ntypes = from_to.size();
	vector<size_t> shelxToXdat(Ntypes), xdatToShelx(Ntypes);
	vector<bool> rotor(size_el, false);
	from_to[0] = 0;
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
	// Creating list of types
	for (size_t i = 1; i < Ntypes; i++)
	{
		from_to[i] = from_to[i-1] + xdat.unit[i];

		for (size_t j = 1; j < Ntypes; j++)
		{
			if (_strcmpi(shelx.sfac[i].c_str(), xdat.sfac[j].c_str()) == 0)
			{
				shelxToXdat[i] = j;
				break;
			}
		}
		for (size_t j = 1; j < Ntypes; j++)
		{
			if (strcmpi(shelx.sfac[j].c_str(), xdat.sfac[i].c_str()) == 0)
			{
				xdatToShelx[j] = i;
				break;
			}
		}
	}

	// Connect first atoms
	for (size_t i = 0; i < size_shelx_atom; i++)
	{
		size_t from = from_to[shelxToXdat[shelx.atom[i].type] - 1];
		size_t to = from_to[shelxToXdat[shelx.atom[i].type]];
		size_t n = size_el;
		flo d = 2;
		Point tcheckpoint;
		for (size_t j = 0; j < size_b; j++)
		{
			auto temp = j % size_el;
			if (temp < from || temp >= to)
				continue;
			flo nd = (xdat.cell.FracToCart() * (basis[j] - shelx.atom[i].point)).r();
			if (nd < d) {
				d = nd;
				n = temp;
				tcheckpoint = basis[j];
				Shift[n] = fPos[n] - basis[j];
			}
		}
		if (n == size_el)
			throw invalid_argument("Bad shelx file.");
		To_n[n] = n;
		Tables[n] = eqMat;
		checkPoint[n] = tcheckpoint;
		rotor[n] = true;
		strcpy_s(xdat.atom[n].label, shelx.atom[i].label);
		for (size_t p = 0; p < 3; p++)
		{
			flo s = shelx.atom[i].point.a[p] - fPos[n].a[p];
			if (s > static_cast<flo>(0.5)) Shift[n].a[p] += 1;
			else if (s <= static_cast<flo>(-0.5)) Shift[n].a[p] -= 1;
		}
	}
	bool changed = true;
	while (changed == true) {
		changed = false;
		for (size_t i = 0; i < size_el; i++)
		{
			if (rotor[i] == false) continue;
			for (size_t s = 0; s < size_s; s++)
			{
				Point tp = shelx.symm[s].GenSymm(fPos[i]);

				size_t from = from_to[xdat.atom[i].type - 1];
				size_t to = from_to[xdat.atom[i].type];
				size_t n = size_el;
				flo d = 2;

				for (size_t j = 0; j < size_b; j++)
				{
					auto temp = j % size_el;
					if (temp < from || temp >= to)
						continue;
					flo nd = (xdat.cell.FracToCart() * (basis[j] - tp)).r();
					if (nd < d) {
						d = nd;
						n = temp;
					}
				}
				if (n == size_el)
					throw invalid_argument("Bad shelx file.");
				if (To_n[n] != -1)
					continue;
				To_n[n] = To_n[i];
				Tables[n] = shelx.symm[s].mat * Tables[i];
				checkPoint[n] = shelx.symm[s].GenSymmNorm(checkPoint[i]);
				rotor[n] = true;
				changed = true;
			}
			rotor[i] = false;
		}
	}
	
	for (size_t i = 0; i < size_el; i++)
	{

		if (i == To_n[i]) {
			const size_t p_size = pList[i].size();
			Shift[i] = checkPoint[i] - fPos[i];
			for (size_t j = 0; j < p_size; j++)
			{
				pList[i][j] += Shift[i];
			}
			continue;
		}

		Tables[i] = std::move(Tables[i].Invert());
		Shift[i] = checkPoint[To_n[i]] - Tables[i] * checkPoint[i];

		// Coping data
		const size_t p_size = pList[i].size();
		pList[To_n[i]].reserve(p_size + pList[To_n[i]].size());
		for (size_t j = 0; j < p_size; j++)
		{
			pList[To_n[i]].push_back(Tables[i] * pList[i][j] + Shift[i]);
		}
		pList[i].clear();
	}
	vector<nsShelxFile::Atom> atombuf;

	for (size_t i = 0; i < size_el; i++)
	{
		if (pList[i].empty()) 
			continue;
		atombuf.push_back(nsShelxFile::Atom(xdat.atom[i].label, xdat.atom[i].type, 1, xdat.cell, move(pList[i]), true));
		pList[i].clear();
	}
	xdat.atom = move(atombuf);
	xdat.LATT = shelx.LATT;
	xdat.symm = shelx.symm;
}