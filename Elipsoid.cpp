#include "stdafx.h"
#include "Elipsoid.h"

extern Matrix FractToDec, DecToFract;

Elipsoid::Elipsoid()
{
	a=b=c=0;
}

Elipsoid::Elipsoid(list <Point> & l, string lab, int num)
{
	label = lab;
	number = num;
	points.swap(l);
}
Elipsoid::~Elipsoid()
{
}

void Elipsoid::Process()
{
	req[1]=OX;
	req[2]=OZ;	
	//{
	//	auto iter = points.begin();
	//	while(iter!=points.end()) {
	//		Point & p = *iter;
	//		p = FractToDec.Transform(*iter);
	//		iter++;
	//	}
	//}
	ShiftCenter(SearchCenter());

	//Point last;
	//for(int i = 1 ; i<=50; i++) {
	//	last = SearchCenterQuad(i);
	//	ShiftCenter(last);
	//}
	bool q = false;
	Line line;
	line = SearchXYLine(50, q);
	RotateZ(line, OY);
	line = SearchYZLine(50, q);
	RotateX(line,OZ);
	line = SearchXYLine(50, q);
	RotateZ(line, OX);
	PlainFirstElipsoid();
	//cout << "kV = " << a*b*c<<endl;
	RotationAndRefinding((int)25, (int)12);
	ListReload();

	Draw_pts();
	cout << "kV = " << a*b*c<<endl;
		
	pbuf.splice(pbuf.begin(),points);
	points.push_back(*pbuf.begin());
	//AddEllipsoidToPoints();
	ReturnToStartPos();	
	//{
	//	auto iter = points.begin();
	//	while(iter!=points.end()) {
	//		Point & p = *iter;
	//		p = DecToFract.Transform(*iter);
	//		iter++;
	//	}
	//	Point p = DecToFract.Transform(Point(a,b,c));
	//}


	//Out(*this);
	if(true);
}
void Elipsoid::UProc()
{
	int size = points.size();
	vector<flo> Xi[3];
	flo avXi[3]={0,0,0};
	for(int i=0;i<3;i++) 
		Xi[i].reserve(size);
	auto iter = points.begin();
	while(iter!= points.end()) {
		Xi[0].push_back(iter->x);
		Xi[1].push_back(iter->y);
		Xi[2].push_back(iter->z);
		iter++;
	}
	for(int i=0;i<3;i++) {
		for(int j=0;j<size;j++) {
			avXi[i] += Xi[i][j];
		}
		avXi[i]/=size;
			
	}
	center=Point(avXi[0],avXi[1],avXi[2]);
	for(int i=0;i<6;i++)
		mat.U[i] = 0;
	for(int i=0;i<size;i++) {
		
		mat.U[0] += (Xi[0][i]-avXi[0])*(Xi[0][i]-avXi[0]);
		mat.U[1] += (Xi[1][i]-avXi[1])*(Xi[1][i]-avXi[1]);
		mat.U[2] += (Xi[2][i]-avXi[2])*(Xi[2][i]-avXi[2]);

		mat.U[3] += (Xi[0][i]-avXi[0])*(Xi[1][i]-avXi[1]);
		mat.U[4] += (Xi[0][i]-avXi[0])*(Xi[2][i]-avXi[2]);
		mat.U[5] += (Xi[1][i]-avXi[1])*(Xi[2][i]-avXi[2]);

	}
	for(int i=0;i<6;i++) {
		
		mat.U[i] /= (flo)size/200;

	}
	//Out(*this);
	if(true);
}