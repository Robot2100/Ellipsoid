#pragma once
#include "stdafx.h"

using namespace std;
typedef float flo;
flo _quad(flo);
int _sign(flo);

	const flo RadtoGrad = 90/1.57079632679489661923;
	const flo GradtoRad = 1.57079632679489661923/90;

void Load(char *, list<Point> & );
typedef Point Line;

const Point OX(1,0,0), OY(0,1,0), OZ(0,0,1);

const int MAX_LINE = 255;

class Elipsoid
{
public:
	Point req[3];
	list <Point> points;
	list <Point> pbuf;
	Point center;
	//Point & X, x, Y, y, Z, z, XY, Xy, xY, xy, XZ, Xz, xZ, xz, YZ, Yz, yZ, yz;
	flo a,b,c;
	Matrix mat;
	string label;
	int number;
	Point _pts[18];


	Elipsoid();
	Elipsoid(list <Point> &, string lab, int num);
	~Elipsoid();
	
	Point SearchCenter()
	{
		Point sum;
		auto iter = points.begin();
		int i = 0;
		while(iter!= points.end()) {
			sum.x+=iter->x;
			sum.y+=iter->y;
			sum.z+=iter->z;
			i++;
			iter++;
		}
		sum.x/=i;
		sum.y/=i;
		sum.z/=i;
		//cout << "SearchCenter: " << sum.x << " "
		//	<< sum.y << " "
		//	<< sum.z << endl;
		return sum;
	}
	Point SearchCenterQuad(unsigned int it)
	{
		//it*=1.5;
		Point sum;
		auto iter = points.begin();
		int i = points.size();
		while(iter!= points.end()) {
			sum.x+=_quad(iter->x)*_sign(iter->x);
			sum.y+=_quad(iter->y)*_sign(iter->y);
			sum.z+=_quad(iter->z)*_sign(iter->z);
			iter++;
		}
		sum.x=_sign(sum.x)*sqrt(abs(sum.x)/i)/it;
		sum.y=_sign(sum.y)*sqrt(abs(sum.y)/i)/it;
		sum.z=_sign(sum.z)*sqrt(abs(sum.z)/i)/it;
		//cout << "SearchCenterQuad: " << sum.x << " "
		//	<< sum.y << " "
		//	<< sum.z << endl;
		
		return sum;
	}
	void ShiftCenter(Point p)
	{
		auto iter = points.begin();
		//int i = 0;
		while(iter!= points.end()) {
			iter->x-=p.x;
			iter->y-=p.y;
			iter->z-=p.z;
			iter++;
		}
		for(int j = 0; j<3;j++)
			req[j]=req[j]-p;
	}
	Point dLine(Line l, Point p)
	{
		Point aa;
		flo r = sqrt(_quad(l.x)+_quad(l.y)+_quad(l.z));
		aa.x = (p.y*l.z - p.z*l.y)/r;
		aa.y = (p.z*l.x - p.x*l.z)/r;
		aa.z = (p.x*l.y - p.y*l.x)/r;
		return aa;
	}
	Line SearchXYLine(int it, bool q)
	{
		flo delta = 1.57079632679489661923;
		Line l1,l2,k;
		l1.x=l2.x=k.x=1;
		flo R = 0;
		Point sum;
		{
			auto iter = points.begin();
			while(iter!= points.end()) {
				Point p = dLine(k,*iter);
				if(q)
					R+=_quad(p.r());
				else
					R+=p.r();
				iter++;
			}
		}

		for(int i =0; i<it;i++) {
			l1.x=cos(delta)*k.x-sin(delta)*k.y;
			l1.y=cos(delta)*k.y+sin(delta)*k.x;		
			l2.x=cos(-delta)*k.x-sin(-delta)*k.y;
			l2.y=cos(-delta)*k.y+sin(-delta)*k.x;
			flo R1=0,R2=0;
			{
				auto iter = points.begin();
				while(iter!= points.end()) {
					Point p = dLine(l1,*iter);
					if(q)
						R1+=_quad(p.r());
					else
						R1+=p.r();
					iter++;
				}		
			
				sum.x = 0;
				sum.y = 0;
				iter = points.begin();
				while(iter!= points.end()) {
					Point p = dLine(l2,*iter);
					if(q)
						R2+=_quad(p.r());
					else
						R2+=p.r();
					iter++;
				}
			}
			if(R>R1 || R>R2) {
				if(R1<R2) {
					k=l1;
					R=R1;
				} else {
					k=l2;
					R=R2;
				}
			}

			delta/=2;
		}
		//k.x/=r;
		//k.y/=r;
		return k;
	
	}
	Line SearchYZLine(int it, bool q)
	{
		flo delta = 1.57079632679489661923;
		Line l1,l2,k;
		l1.y=l2.y=k.y=1;
		flo R = 0;
		Point sum;
		{
			auto iter = points.begin();
			while(iter!= points.end()) {
				Point p = dLine(k,*iter);
				if(q)
					R+=_quad(p.r());
				else
					R+=p.r();
				iter++;
			}
		}

		for(int i =0; i<it;i++) {
			l1.y=cos(delta)*k.y-sin(delta)*k.z;
			l1.z=cos(delta)*k.z+sin(delta)*k.y;		
			l2.y=cos(-delta)*k.y-sin(-delta)*k.z;
			l2.z=cos(-delta)*k.z+sin(-delta)*k.y;
			flo R1=0,R2=0;
			{
				auto iter = points.begin();
				while(iter!= points.end()) {
					Point p = dLine(l1,*iter);
					if(q)
						R1+=_quad(p.r());
					else
						R1+=p.r();
					iter++;
				}		
			
				sum.y = 0;
				sum.z = 0;
				iter = points.begin();
				while(iter!= points.end()) {
					Point p = dLine(l2,*iter);
					if(q)
						R2+=_quad(p.r());
					else
						R2+=p.r();
					iter++;
				}
			}
			if(R>R1 || R>R2) {
				if(R1<R2) {
					k=l1;
					R=R1;
				} else {
					k=l2;
					R=R2;
				}
			}

			delta/=2;
		}
		//k.x/=r;
		//k.y/=r;
		return k;
	
	}
	flo RotateZ(Line line, Point c)
	{
		if(line.x==0&&line.y==0) return 0;
		if(c.x==0&&c.y==0) return 0;
		c.z=0;
		c.normalize();
		line.z=0;
		line.normalize();
		flo angle1 = acos(line.y);
		if(line.x<0) angle1 = (3.14159265358979323846*2)-angle1;
		flo angle2 = acos(c.y);
		if(c.x<0) angle2 = (3.14159265358979323846*2)-angle2;
		flo angle = -angle1+angle2;
		//if(Point(c+line).r()<sqrt((flo)2.0)) angle +=3.14159265358979323846;
		//if((c.x*line.x+c.y*line.y)>0) angle +=3.14159265358979323846;
		//flo transangle = angle * (180/3.14159265358979323846);
		//cout << "RotateZ: rotation " << angle*(90/1.57079632679489661923) << " grad." << endl;
		auto iter = points.begin();
		while(iter!= points.end()) {
			Point p = *iter;
			iter->x=cos(angle)*p.x+sin(angle)*p.y;
			iter->y=cos(angle)*p.y-sin(angle)*p.x;	
			//transangle = (atan(iter->y/iter->x)-atan(p.y/p.x));
				
			//if((iter->x*p.x+iter->y*p.y)>0) angle +=3.14159265358979323846;
			//transangle *=(180/3.14159265358979323846);
			iter++;
		}
		for(int j = 0; j<3;j++) {
			Point p = req[j];
			req[j].x=cos(angle)*p.x+sin(angle)*p.y;
			req[j].y=cos(angle)*p.y-sin(angle)*p.x;
		}
		return angle;
	}
	flo RotateX(Line line, Point c)
	{
		
		if(line.z==0&&line.y==0) return 0;
		if(c.z==0&&c.y==0) return 0;
		c.x=0;
		c.normalize();
		line.x=0;
		line.normalize();
		//flo angle = -atan(c.y/c.x)+atan(line.y/line.x);

		
		flo angle1 = acos(line.z);
		if(line.y<0) angle1 = (3.14159265358979323846*2)-angle1;
		flo angle2 = acos(c.z);
		if(c.y<0) angle2 = (3.14159265358979323846*2)-angle2;
		flo angle = -angle1+angle2;

		//if(Point(c+line).r()<sqrt((flo)2.0)) angle +=3.14159265358979323846;
		//if((c.x*line.x+c.y*line.y)>0) angle +=3.14159265358979323846;
		//cout << "RotateX: rotation " << angle*(90/1.57079632679489661923) << " grad." << endl;
		auto iter = points.begin();
		while(iter!= points.end()) {
			Point p = *iter;
			iter->y=cos(angle)*p.y+sin(angle)*p.z;
			iter->z=cos(angle)*p.z-sin(angle)*p.y;	
			iter++;
		}
		for(int j = 0; j<3;j++) {
			Point p = req[j];
			req[j].y=cos(angle)*p.y+sin(angle)*p.z;
			req[j].z=cos(angle)*p.z-sin(angle)*p.y;
		}
		return angle;
	}
	Point Checking()
	{
		Point res;
		flo R=1;
		auto iter = points.begin();
		int i=0;
		int size = points.size();
		while(i<size) {
			flo r1 = iter->InEllipsoid(a,b,c);
			if(r1<=1) {
				pbuf.push_back(*iter);
				points.erase(iter++);
			} else {
				if(r1>R) {
					R=r1;
					res=*iter;
				}
				iter++;
			}
			i++;
		}
		return res;
	}
	void ListReload()
	{
		points.splice(points.begin(),pbuf);
	}
	void Grow(Point p)
	{
		flo da,db,dc;
		da=sqrt(abs(p.x));
		db=sqrt(abs(p.y));
		dc=sqrt(abs(p.z));

		flo grid = sqrt(p.InEllipsoid(da*a,db*b,dc*c));
		da=abs(da)*grid;
		db=abs(db)*grid;
		dc=abs(dc)*grid;
		a*=da>1?da:1;
		b*=db>1?db:1;
		c*=dc>1?dc:1;
	}
	void PlainFirstElipsoid()
	{
		a=b=c=0;
		Point am,bm,cm;
		auto iter = points.begin();
		ListReload();
		while(iter!= points.end()) {
			if(abs(iter->x)>a) {
				a = abs(iter->x);
				am = *iter;
			}
			if(abs(iter->y)>b) {
				b = abs(iter->y);
				bm = *iter;
			}
			if(abs(iter->z)>c) {
				c = abs(iter->z);
				cm = *iter;
			}
			iter++;
		}
		Point t = Checking();


		while(!points.empty())
		{
			Grow(t);
			t=Checking();
		}
		ListReload();
	}
	void Draw_pts()
	{
		_pts[0]=Point(a,0,0);
		_pts[1]=Point(-a,0,0);
		_pts[2]=Point(0,b,0);
		_pts[3]=Point(0,-b,0);
		_pts[4]=Point(0,0,c);
		_pts[5]=Point(0,0,-c);
		_pts[6]=Point(a*sqrt(0.5f),b*sqrt(0.5f),0);
		_pts[7]=Point(a*sqrt(0.5f),-b*sqrt(0.5f),0);
		_pts[8]=Point(-a*sqrt(0.5f),b*sqrt(0.5f),0);
		_pts[9]=Point(-a*sqrt(0.5f),-b*sqrt(0.5f),0);
		_pts[10]=Point(a*sqrt(0.5f),0,c*sqrt(0.5f));
		_pts[11]=Point(a*sqrt(0.5f),0,-c*sqrt(0.5f));
		_pts[12]=Point(-a*sqrt(0.5f),0,c*sqrt(0.5f));
		_pts[13]=Point(-a*sqrt(0.5f),0,-c*sqrt(0.5f));
		_pts[14]=Point(0,b*sqrt(0.5f),c*sqrt(0.5f));
		_pts[15]=Point(0,b*sqrt(0.5f),-c*sqrt(0.5f));
		_pts[16]=Point(0,-b*sqrt(0.5f),c*sqrt(0.5f));
		_pts[17]=Point(0,-b*sqrt(0.5f),-c*sqrt(0.5f));
	}
	Point RefindingZero(int Ntry)
	{
		Point ret(0,0,0);
		ListReload();
		PlainFirstElipsoid();
		Point Qreq(req[0]);
		flo exA = a;
		flo exB = b;
		flo exC = c;
		flo exV = a*b*c;
		flo separ = 0.5;
		Point lastlake;
		for(int i=0; i<Ntry;i++) {
			PlainFirstElipsoid();
			Point lake(0,0,0);
			auto iter = points.begin();
			while(iter!= points.end()) {
				if(iter->InEllipsoid(a,b,c)==1) {
					lake = lake + *iter;
				}
				iter++;
			}
			if(lastlake!=lake) {
				separ = 0.5;
				i=0;
				lastlake = lake;
			}
			lake.x*=separ;
			lake.y*=separ;
			lake.z*=separ;

			ListReload();
			ShiftCenter(lake);
			PlainFirstElipsoid();
			separ*=0.5;
			if(flo(a*b*c)<exV) {
				exV=a*b*c;
				//exA = a;
				//exB = b;
				//exC = c;
				ret=ret-lake;
			} else {
				lake.x*=-1;
				lake.y*=-1;
				lake.z*=-1;
				//a=exA;
				//b=exB;
				//c=exC;
				ShiftCenter(lake);
			}
		}
		a=b=c=0;
		
		ListReload();
		PlainFirstElipsoid();
		//cout << "RefindingZero:\tkV = " << a*b*c<<endl; 
		return ret;
	}
	void RotationAndRefinding(int Ntry,int Nangle)
	{
		int total=0;
		//cout << "RotationAndRefinding Start: kV = " << a*b*c<<endl; 
		Ntry=25;
		Nangle=12;
		RefindingZero(Ntry);
		Point Qreq(0,0,0);
		flo exV = a*b*c;
		flo angle = 3.14159265358979323846/3;
		ListReload();
		for(int i=1; i<=Nangle;i++) {
			if(total>100) break;
			//exV = a*b*c;
			//cout << "Circle " << i <<"\tkV = " << a*b*c <<endl; 
			RotateZ(Line(cos(angle),sin(angle),0),OX);
			a=b=c=0;
			Qreq = RefindingZero(Ntry);
			if(flo(a*b*c)<exV) {
				exV=a*b*c;
				//RefindingZero(Ntry);
				i--;
				total++;
				continue;
			} else {
				ShiftCenter(Qreq);
				RotateZ(OX,Point(cos(angle),sin(angle),0));
			}
			RotateZ(OX,Point(cos(angle),sin(angle),0));
			Qreq = RefindingZero(Ntry);
			if(flo(a*b*c)<exV) {
				exV=a*b*c;
				//RefindingZero(Ntry);
				i--;
				total++;
				continue;
			} else {
				ShiftCenter(Qreq);
				RotateZ(Line(cos(angle),sin(angle),0),OX);
			}

			RotateX(Line(0,cos(angle),sin(angle)),OY);
			Qreq = RefindingZero(Ntry);
			if(flo(a*b*c)<exV) {
				exV=a*b*c;
				//RefindingZero(Ntry);
				i--;
				total++;
				continue;
			} else {
				ShiftCenter(Qreq);
				RotateX(OY,Point(0,cos(angle),sin(angle)));
			}
			RotateX(OY,Point(0,cos(angle),sin(angle)));
			a=b=c=0;
			Qreq = RefindingZero(Ntry);
			if(flo(a*b*c)<exV) {
				exV=a*b*c;
				i--;
				total++;
				continue;
			} else {
				ShiftCenter(Qreq);
				RotateX(Line(0,cos(angle),sin(angle)),OY);
			}
			Qreq = RefindingZero(Ntry);
			if(flo(a*b*c)<exV) {
				exV=a*b*c;
			} else {
				ShiftCenter(Qreq);
			}
			angle *= 0.5;
		}

	}
	void AddEllipsoidToPoints() {
		for(int i = 0;i<18;i++)
			points.push_back(_pts[i]);
	}
	void ReturnToStartPos()
	{
		mat=Matrix();
		mat.U[0]=1/_quad(a);
		mat.U[1]=1/_quad(b);
		mat.U[2]=1/_quad(c);
		center = Point()-req[0];
		ShiftCenter(req[0]);
		mat.RotateZ(RotateZ(req[2],OY));
		mat.RotateX(RotateX(req[2],OZ));
		mat.RotateZ(RotateZ(req[1],OX));
	}
	void Process();
	void UProc();


	void FFTAnalysis(double *AVal, double *FTvl, int Nvl, int Nft) {
		const double TwoPi = 6.283185307179586;
		int i, j, n, m, Mmax, Istp;
		double Tmpr, Tmpi, Wtmp, Theta;
		double Wpr, Wpi, Wr, Wi;
		double *Tmvl;

		n = Nvl * 2; Tmvl = new double[n];

		for (i = 0; i < n; i+=2) {
			Tmvl[i] = 0;
			Tmvl[i+1] = AVal[i/2];
		}

		i = 1; j = 1;
		while (i < n) {
			if (j > i) {
				Tmpr = Tmvl[i]; Tmvl[i] = Tmvl[j]; Tmvl[j] = Tmpr;
				Tmpr = Tmvl[i+1]; Tmvl[i+1] = Tmvl[j+1]; Tmvl[j+1] = Tmpr;
			}
			i = i + 2; m = Nvl;
			while ((m >= 2) && (j > m)) {
				j = j - m; m = m >> 1;
			}
			j = j + m;
		}

		Mmax = 2;
		while (n > Mmax) {
			Theta = -TwoPi / Mmax; Wpi = sin(Theta);
			Wtmp = sin(Theta / 2); Wpr = Wtmp * Wtmp * 2;
			Istp = Mmax * 2; Wr = 1; Wi = 0; m = 1;

			while (m < Mmax) {
				i = m; m = m + 2; Tmpr = Wr; Tmpi = Wi;
				Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
				Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

				while (i < n) {
					j = i + Mmax;
					Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j-1];
					Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j-1];

					Tmvl[j] = Tmvl[i] - Tmpr; Tmvl[j-1] = Tmvl[i-1] - Tmpi;
					Tmvl[i] = Tmvl[i] + Tmpr; Tmvl[i-1] = Tmvl[i-1] + Tmpi;
					i = i + Istp;
				}
			}

			Mmax = Istp;
		}

		for (i = 0; i < Nft; i++) {
			j = i * 2; FTvl[i] = 2*sqrt(pow(Tmvl[j],2) + pow(Tmvl[j+1],2))/Nvl;
		}

		delete []Tmvl;
	}
};



void Out(Elipsoid &);
void OutIns(Elipsoid & el);
