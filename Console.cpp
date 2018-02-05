// Console.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "../Includes/Includes.h"


using namespace std;
const int MAX_LINE = 128;

std::vector<Point> fPos;

Cell LoadXDATCAR(vector<Elipsoid> & pList);
void OutIns(Elipsoid & el);
int _sign(flo a);
void Analize_symmety(nsShelxFile::ShelxData & shelx, vector<Elipsoid> & pList);


int main(int argn, char * argv[]) {
	vector<string> OUT_st;
	vector<Elipsoid> El;
	Cell cell = LoadXDATCAR(El);

	for (int i = 0; i < El.size();i++) {
		El[i].DefineCell(cell);
	}

	ifstream old("old.res");
	if (old.is_open() == true) {
		nsShelxFile::ShelxData shelx(old);
		cout << "Symmetry found." << endl;
		Analize_symmety(shelx, El);
		ofstream temp("a.ins", ios::app);
		temp << "LATT " << shelx.LATT << endl;
		for (int i = 0; i < shelx.symm.size(); i++) {
			if(shelx.symm[i].LATT == false) temp << "SYMM " << shelx.symm[i].sstr << endl;
		}
		temp.close();
	}
	old.close();

	int j = 0;
	int size = El.size();
	for(int i = 0; i < size;i++) {
		ELIPSOID_RESULT res = El[i].CalculateDinmat();
		if (res != ELIPSOID_RESULT::OK) {
			cout << "Error" << endl;
			return 1;
		}
		OutIns(El[i]);
	}
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
	Cell cell(Matrix(&U[0], 3, 3), true);

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
			fPos.resize(NAtoms);
			tv2.swap(tempVec);
			for (int i = 0; i < NAtoms; i++) {
				fPos[i] = tv2[i].front();
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
int _sign(flo a) {
	if (a < 0) return -1;
	else return 1;
}
void Analize_symmety(nsShelxFile::ShelxData & shelx, vector<Elipsoid> & pList) {
	using namespace nsShelxFile;
	shelx.cell = (*(pList[0].pCell));
	size_t size_s = shelx.symm.size();
	size_t size_el = pList.size();
	vector<vector<vector<SYMM*> > > table(size_el, vector<vector<SYMM*> >(size_el));
	vector<vector<Point> > addtable(size_el, vector<Point>(size_el));
	constexpr int _p = 1; 
	constexpr size_t _d = 2*_p+1;
	constexpr size_t _sizemod = (_d*_d*_d);
	size_t size_b = size_el*_sizemod;
	// Неподвижная Декартова суперячейка
	vector<Point> basis(size_b);
	// Заполнение Базиса копией начальных положений атомов суперячейки
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
	// Составление первого круга таблицы связности
	for (int s = 0; s < size_s; s++) {
		// Декартовы точки после операции симметрии
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
		// Были ли изменения
		bool changed = false;
		// Полное сведение таблицы связности 
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
	// Перевод точек по схеме
	vector<Point> dcheck;
	{
		size_t size_l = pList[0].vecPoints.size();
		for (int i = 0; i < size_el; i++) {
			if (table[i].size() == 0) continue;
			for (int j = 0; j < size_el; j++) {
				if (table[i][j].size() == 0 || pList[j].vecPoints.size() == 0) continue;
				size_t k = pList[i].vecPoints.size();
				pList[i].vecPoints.reserve(k + pList[j].vecPoints.size());
				pList[i].vecPoints.insert(pList[i].vecPoints.end(),
					std::make_move_iterator(pList[j].vecPoints.begin()),
					std::make_move_iterator(pList[j].vecPoints.end()));
				pList[j].vecPoints.clear();
				pList[j].useVecPoints = false;
				size_t it = pList[i].vecPoints.size();

				for (int k2 = 0; k2 < table[i][j].size(); k2++) {
					pList[i].vecPoints[k] = (table[i][j][k2]->RetroGenSymm(pList[i].vecPoints[k]));
				}
				Point dfix = pList[i].vecPoints[k] + Point(100.5, 100.5, 100.5) - fPos[i];
				dfix = (Point(100, 100, 100) - Point(int(dfix.a[0]), int(dfix.a[1]), int(dfix.a[2])));
				pList[i].vecPoints[k] += dfix;

				for (++k; k < it; k++) {
					for (int k2 = 0; k2 < table[i][j].size(); k2++) {
						pList[i].vecPoints[k] = (table[i][j][k2]->RetroGenSymm(pList[i].vecPoints[k]));
					}
					pList[i].vecPoints[k] += dfix;
				}
				dcheck.push_back(pList[i].vecPoints.back());
			}
		}
	}
	vector<Elipsoid> temp;
	for (int i = 0; i < size_el; i++) {
		if(pList[i].useVecPoints == true)
			temp.push_back(move(pList[i]));
	}
	temp.swap(pList);
}
//void Analize_symmety_old(nsShelxFile::ShelxData & shelx, vector<Elipsoid> & pList) {
//	size_t size_a = shelx.atom.size();
//	size_t size_s = shelx.symm.size();
//	size_t size_el = pList.size();
//
//	for (int i = 0; i < size_el; i++) {
//		char buf[5];
//		strcpy(buf, pList[i].label.c_str());
//		_strupr_s(buf);
//		for (int j = 1; j < shelx.sfac.size(); j++) {
//			if (shelx.sfac[j].compare(buf) != 0) continue;
//			pList[i].type = j;
//			break;
//		}
//	}
//	vector<nsShelxFile::Atom> atsh(std::move(shelx.GenerateSymmAtom()));
//	size_t size_ats = atsh.size();
//	for (int j = -1; j <= 1; j++) {
//		for (int k = -1; k <= 1; k++) {
//			for (int l = -1; l <= 1; l++) {
//				if (j == 0 && k == 0 && l == 0) continue;
//				for (int i = 0; i < size_ats; i++) {
//					atsh.push_back(atsh[i]);
//					(atsh.back()).point += Point(j, k, l);
//				}
//			}
//		}
//	}
//	size_ats = atsh.size();
//	for (int i = 0; i < size_el; i++) {
//		nsShelxFile::Atom * ap = &atsh[0];
//		flo dp = (pList[i].vecPoints[0] - shelx.cell.FracToCart() * (atsh[0].point)).r();
//		for (int j = 1; j < size_ats; j++) {
//			if (atsh[j].type != pList[i].type) continue;
//			flo tempdp = (pList[i].vecPoints[0] - shelx.cell.FracToCart() * (atsh[j].point)).r();
//			if (tempdp >= dp) continue;
//			dp = tempdp;
//			ap = &(atsh[j]);
//		}
//		ap->point = pList[i].fracCenter;
//		ap->dinmat = pList[i].dinmat;
//		ap->type += 256;
//	}
//	{
//		vector<nsShelxFile::Atom> tempA;
//		tempA.swap(atsh);
//		for (int i = 0; i < size_ats; i++) {
//			if ((tempA[i].type & 256) != 0) {
//				tempA[i].type -= 256;
//				atsh.push_back(std::move(tempA[i]));
//			}
//		}
//	}
//	size_ats = atsh.size();
//	vector<nsShelxFile::Atom> acpy(shelx.atom);
//	for (int i = 0; i < size_a; i++) {
//		acpy[i].dinmat = Dinmat();
//		acpy[i].point = Point();
//	}
//	for (int i = 0; i < size_a; i++) {
//		size_t count = 0;
//		for (int j = 0; j < size_ats; j++) {
//			if (strcmp(atsh[j].label, shelx.atom[i].label) != 0) continue;
//			vector<nsShelxFile::Atom> atsh2(std::move(shelx.GenerateSymmAtom(atsh[j])));
//			nsShelxFile::Atom * ap = &atsh2[0];
//			flo dp = (shelx.cell.FracToCart() * (shelx.atom[i].point - (atsh2[0].point))).r();
//			for (int k = 0; k < atsh2.size(); k++) {
//				flo tempdp = (shelx.cell.FracToCart() * (shelx.atom[i].point - (atsh2[k].point))).r();
//				if (tempdp >= dp) continue;
//				dp = tempdp;
//				ap = &(atsh2[k]);
//			}
//			//acpy[i].point += ap->point;
//			//for (int k = 0; k < 6; k++)
//			//acpy[i].dinmat.U[k] += ap->dinmat.U[k];
//			count++;
//		}
//		cout << "Atom " << i << " connected to " << count << " atoms." << endl;
//		acpy[i].point /= count;
//		for (int k = 0; k < 6; k++)
//			acpy[i].dinmat.U[k] /= count;
//	}
//
//	shelx.atom.swap(acpy);
//}
