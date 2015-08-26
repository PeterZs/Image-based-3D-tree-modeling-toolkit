#include <iostream>	
#include <string>   
#include <iomanip>  
#include <sstream>  
#include <stdio.h>
#include <math.h>
#include<ctime>
#include<cstdlib>
#include <string>


#include "opencv\cv.h"
#include "opencv\highgui.h"
#include "opencv\cvaux.h"
#include "opencv\cxcore.h"

#include <windows.h>
#include "glut.h"
#include "bitmap.h"

#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtx/transform.hpp>

#include <fstream>

using namespace std;
using namespace cv;
#define CV_EVENT_RBUTTONDOWN 2

#define maxW 480
#define maxH 640
#define LeastbranchLength 20

const int picname = 15;
const string picstr = "tree15.jpg";

Mat image;
IplImage* img;
//-----------OpenGL---------------------
int eye = 4;
int angleX = 0,angleY = 0,angleZ = 0;
GLfloat scaleX = 0.0f,scaleY = 0.0f,scaleZ = 0.0f;
int mouseX = -1,mouseY = -1;
#define GLUT_WHEEL_UP 3
#define GLUT_WHEEL_DOWN 4
#define projectWidth 5.0
BITMAPINFO* TexInfo;
GLubyte* TexBits;
int number = 1;
//-----------MousePoint variants--------
struct point
{
	Point p;
	int parent;
	int num;
};

struct tdpoint{
	tdpoint() {}
	tdpoint(double _x,double _y){ x = _x; y = _y;}
	double x,y;
};

Point prev_pt = Point(-1,-1);
Point cur_pt = Point(-1,-1);
bool yellow = true;
bool first =  true;
vector<point>skltnpts;  
int pos = -1;

//-----------Crown
vector<Point> crown;

//-----------GetOriginBran variants-------
struct D3Point{
	D3Point() {}
	D3Point(double _x,double _y,double _z){ x = _x; y = _y; z = _z;}
	void set(double _x,double _y,double _z) { x = _x; y = _y; z = _z;}
	double x,y,z;
};
struct branch{
	int left,right;
	double l,r,i,a;
	double x,y;
	double z;
	bool iftip;
	int p;
};
vector<branch> vbr;
int tableI[100];
int tableA[100];
int tableN[100];

D3Point Nv(0,1,0);

vector<double> b2;
double LratioLower = 0.5;
double LratioUpper = 0.8;
const double Rratio = 0.8;

int rnSize,riSize,raSize;
vector<int> rnTable;
vector<int> raTable;
vector<int> riTable;
vector<bool> ifexp;
int Origineva;
int Originbran;
vector<double> incTable;

double rmap[maxH][maxW*2];
int eval[maxH][maxW*2][2];
double ub[maxH],lb[maxH];
double newmap[maxH][maxW*2][2];
int tmpt[maxH];
int times;

glm::vec3 YAxis(0,1,0);
glm::vec3 ZAxis(0,0,1);

static GLuint m_Texture[2];
GLfloat Tcb = 47;
GLfloat Tcg = 49;
GLfloat Tcr = 35;
GLfloat Lcb = 114;
GLfloat Lcg = 147;
GLfloat Lcr = 76;

bool FirstGenerate = true;
int curbranch = 1;
int ifleaf = -1;


int bskpnt[20];
int curbsk;

const double PI = 3.14159265;
const double AVar = 15.0;
const double IVar = 10.0;
const double EvaImprove = 0.0;
double sumInc = 0.0;
double sumVar = 0.0;
double sumAzi = 0.0;
double sumAziVar = 0.0;
int sumNumber = 0;
int TotMapNumber = 0;

#define LEAF_STEP 15
#define LEAF_SIZE 2

//--Record the Point on the 2D tree crown
void CrownPoint( int event, int x, int y, int flags, void* zhang)
{
    if( event == CV_EVENT_LBUTTONUP || !(flags & CV_EVENT_FLAG_LBUTTON) )
	{
        //prev_pt = cvPoint(-1,-1);
		prev_pt = Point(-1,-1);
	}
    else if( event == CV_EVENT_LBUTTONDOWN )
	{
        prev_pt = cvPoint(x,y);
		crown.push_back(prev_pt);
		line( image, prev_pt, prev_pt, cvScalar(0,0,255), 3, 8, 0 );
		imshow("Display window",image);
	}
    else if( event == CV_EVENT_MOUSEMOVE && (flags & CV_EVENT_FLAG_LBUTTON) )
    {
        CvPoint pt = cvPoint(x,y);
		if( prev_pt.x < 0 || prev_pt.y < 0)
            prev_pt = pt;
        //cvLine( inpaint_mask, prev_pt, pt, cvScalar(255,0,0), 5, 8, 0 );
		line( image, prev_pt, pt, cvScalar(0,255,255), 3, 8, 0 );
		crown.push_back(pt);
		prev_pt = pt;
        imshow("Display window",image);
    }
}

//Record the branch
void MousePoint(int event, int x, int y, int flags, void* zhang)
{
    if( event == CV_EVENT_LBUTTONDOWN )
	{
		cur_pt = Point(x,y);
		point* p = new point;
		p->p = cur_pt;
		p->num = 0;
		p->parent = pos;
		if(pos != -1)
		{
			skltnpts[pos].num++;
		}
		skltnpts.push_back(*p);
		skltnpts[skltnpts.size()-1].parent = pos;
		if(first)
		{
			prev_pt = Point(x,y);
			first = false;
			//p->parent = -1;
		}
        if(yellow)
		    line( image, prev_pt, cur_pt, Scalar(255,0,0), 2, 8, 0 );
		else
            line( image, prev_pt, cur_pt, Scalar(0,255,255), 2, 8, 0 );

		imshow("Display window",image);
	}
	else if( event = CV_EVENT_RBUTTONDOWN && (flags & CV_EVENT_FLAG_RBUTTON) )
	{
		//find the nearest point in the skeleton as the current parent nodeW
		if(skltnpts.size() != 0)
		{
			Point cursor = Point(x,y);
			unsigned int dis = 999999;
			for(int i=0;i<(int)skltnpts.size();i++)
			{
				int d = (cursor.x - skltnpts[i].p.x) * (cursor.x - skltnpts[i].p.x) 
					  + (cursor.y - skltnpts[i].p.y) * (cursor.y - skltnpts[i].p.y);
				if(d < dis)
				{
					dis = d;
					pos = i;
				}
			}
			prev_pt = skltnpts[pos].p;
		}
		else
		{
			pos = -1;
		}
	}
}

//Get the Inclination of the two braches
double GetI(D3Point x,D3Point y){
	return acos((x.x*y.x+x.y*y.y+x.z*y.z)/sqrt(x.x*x.x+x.y*x.y+x.z*x.z)/sqrt(y.x*y.x+y.y*y.y+y.z*y.z));
}

//Get the angle of rotation
double GetA(D3Point x,D3Point y,D3Point nv,D3Point& ev){
	ev.x = x.y*y.z - x.z*y.y;
	ev.y = x.z*y.x - x.x*y.z;
	ev.z = x.x*y.y - x.y*y.x;
	return GetI(nv,ev);
}

//calculate the inclination and angle of rotation of the main branches
void GetVbr(int index,D3Point v,int aha){
	D3Point ev(0,0,0);
	D3Point tv(0,0,0);
	int ti,tti;
	int i,j,tmp,t1,t2,t3,t4,tmp1,tmp2,tmp3,tmp4;
	double idealL,x1,x2,y1,y2,t,dtmp,new1;
	int Idata[5] = {10,8,6,4,2};
	if (index == 1){
		vbr[index].a = 0;
		vbr[index].i = 0;
		vbr[index].r = 7.5;
		vbr[index].l = b2[index];
		vbr[index].y = 0;
		ev.x = 0;
		ev.y = 1;
		ev.z = 0;
	}
	else{
		ti = skltnpts[index].parent;
		idealL = vbr[ti].l * LratioLower;
		if (idealL < b2[index]) idealL = b2[index];
		if (aha % 2 == 1) vbr[index].y = 10 + vbr[ti].y;
		else vbr[index].y = -10 + vbr[ti].y;
		idealL = sqrt(b2[index] * b2[index] + 100);
		tti = skltnpts[ti].parent;
		
		D3Point *p1 = new D3Point(vbr[ti].x-vbr[tti].x,vbr[ti].y-vbr[tti].y,vbr[ti].z-vbr[tti].z);
		D3Point *p2 = new D3Point(vbr[index].x-vbr[ti].x,vbr[index].y-vbr[ti].y,vbr[index].z-vbr[ti].z);
		
		vbr[index].i = GetI(*p1,*p2) / PI * 180.0;
		if (vbr[index].i > 90) vbr[index].i = 180.0 - vbr[index].i;
		vbr[index].a = GetA(*p1,*p2,v,ev) / PI * 180.0;
		//if (vbr[index].a > 180) vbr[index].a = 360.0 - vbr[index].a;
		//if (index == 2 || index == 5) vbr[index].i -= 20;
		tmp1 = (int)vbr[index].i;
		tmp2 = (int)vbr[index].a;
		tmp2 = tmp2 % 90;
		cout << index << " " << skltnpts[index].parent << " " << vbr[index].y << " " << vbr[skltnpts[index].parent].y << endl;

		vbr[index].r = vbr[ti].r * Rratio;
		vbr[index].l = idealL;
		
		int ty;
		idealL = vbr[ti].l * LratioUpper;	
		if (idealL < b2[index]) idealL = b2[index];
		ty = (int)sqrt(idealL * idealL - b2[index] * b2[index]) + vbr[ti].y;
		ty = 0;
		p1 = new D3Point(vbr[ti].x-vbr[tti].x,vbr[ti].y-vbr[tti].y,vbr[ti].z-vbr[tti].z);
		p2 = new D3Point(vbr[index].x-vbr[ti].x,ty-vbr[ti].y,vbr[index].z-vbr[ti].z);
		dtmp = GetI(*p1,*p2) / PI * 180;
		if (dtmp > 90) dtmp = 180.0 - dtmp;
		tmp3 = (int)dtmp;
		dtmp = GetA(*p1,*p2,v,tv) / PI * 180;
		tmp4 = (int)dtmp;
		tmp4 = tmp4 % 90;
		//cout << tmp2 << " " << tmp4 << endl;
		for ( j = max(min(tmp1,tmp3),0); j <= min(90,max(tmp1,tmp3)); ++ j)
			tableI[j] += 10;
		for ( j = max(0,min(tmp2,tmp4)); j <= min(90,max(tmp2,tmp4)); ++ j)
			tableA[j] += 10;
		t1 = min(tmp1,tmp3);
		t2 = max(tmp1,tmp3);
		t3 = min(tmp2,tmp4);
		t4 = max(tmp2,tmp4);
		for ( j = max(t1-4,0); j <= max(t1-1,0); ++ j)
			tableI[j] += Idata[abs(t1-j)];
		for ( j = min(t2+1,90); j <= min(t2+4,90); ++ j)
			tableI[j] += Idata[abs(j-t2)];
		for ( j = max(t3-4,0); j <= max(t3-1,0); ++ j)
			tableA[j] += Idata[abs(t3-j)];
		for ( j = min(t4+1,90); j <= min(t4+4,90); ++ j)
			tableA[j] += Idata[abs(j-t4)];
	}
	int cur;
	cur = vbr[index].left;
	aha = 1;
	while (cur != -1){
		GetVbr(cur,ev,aha);
		cur = vbr[cur].right;
		++ aha;
	}
}

//Get the data of main branches and calculate the probability distribution
void GetOriginBran(){
	int i,j,ti,tti,tmp,t1,t2,t3,t4,tmp1,tmp2,tmp3,tmp4;
	double idealL,x1,x2,y1,y2,t,dtmp,new1;
	int Idata[5] = {10,8,6,4,2};
	int cur;
	rnTable.clear();
	riTable.clear();
	raTable.clear();
	b2.clear();

	memset(tableN,0,sizeof(tableN));
	memset(tableI,0,sizeof(tableI));
	memset(tableA,0,sizeof(tableA));
	for ( i = 0; i < skltnpts.size(); ++ i)
		++ tableN[skltnpts[i].num];
	for ( i = 0; i < skltnpts.size(); ++ i){
		branch* pb = new branch();
		pb->left = -1;
		pb->right = -1;
		if (i > 0) pb->p = skltnpts[i].parent;
		vbr.push_back(*pb);
		ifexp.push_back(true);
	}
	for ( i = 0; i < skltnpts.size(); ++ i){
		ti = skltnpts[i].parent;
		b2.push_back(i == 0? 0 : sqrt((double)(skltnpts[ti].p.x-skltnpts[i].p.x) * (skltnpts[ti].p.x-skltnpts[i].p.x) + (skltnpts[ti].p.y-skltnpts[i].p.y) * (skltnpts[ti].p.y-skltnpts[i].p.y)));
	}
	Point ori;
	ori.x = skltnpts[0].p.x;
	ori.y = skltnpts[0].p.y;
	for ( i = 0; i < skltnpts.size(); ++ i){
		vbr[i].x = skltnpts[i].p.x-ori.x;
		vbr[i].z = ori.y-skltnpts[i].p.y;
	}
	for ( i = 1; i < vbr.size(); ++ i){
		ti = skltnpts[i].parent;
		vbr[i].right = vbr[ti].left;
		vbr[ti].left = i;
	}
	GetVbr(1,Nv,1);
	Originbran = vbr.size()-1;
	//cout << "ok" << endl;
//---------------------random table--------------------
	for ( i = 0; i <= 60; ++ i)
		tableI[i] += i / 4;
	rnSize = 0;
	for ( i = 1; i <= 4; ++ i){
		tableN[i] = tableN[i] * 4 + 1;
		if (i == 1) -- tableN[i];
		if (i == 3) tableN[i] += 3;
		tableN[i] *= i;
		rnSize += tableN[i];
	}
	//cout << rnSize << endl;
	cur = 1;
	for ( i = 1; i <= rnSize; ++ i){
		while (tableN[cur] == 0) ++ cur;
		rnTable.push_back(cur);
		-- tableN[cur];
	}
	//cout << "ok" << endl;
	//cout << "hehe" << endl;
	riSize = 0;
	for ( i = 0; i <= 90; ++ i){
		riSize += tableI[i];
	}
	//cout << riSize << endl;
	cur = 0;
	for ( i = 1; i <= riSize; ++ i){
		while (tableI[cur] == 0 && cur <= 90) ++ cur;
		riTable.push_back(cur);
		-- tableI[cur];
	}
	//cout << "ok" << endl;
	raSize = 0;
	for ( i = 0; i <= 90; ++ i){
		tableA[i] += 3;
		raSize += tableA[i];
	}
	//cout << raSize << endl;
	cur = 0;
	for ( i = 1; i <= raSize; ++ i){
		while (tableA[cur] == 0 && cur <= 90) ++ cur;
		raTable.push_back(cur);
		-- tableA[cur];
	}
	//cout << "ok" << endl;
	if (FirstGenerate){
		FirstGenerate = false;
		int da[maxH];
		int xiao[maxH];
		for ( i = 0; i < maxH; ++ i){
			da[i] = -100000;
			xiao[i] = 100000;
		}
		memset(rmap,0,sizeof(rmap));
		for ( i = 0; i < crown.size(); ++ i){
			int x = crown[i].x - ori.x;
			int y = ori.y - crown[i].y;
			if (x < xiao[y]) xiao[y] = x;
			if (x > da[y]) da[y] = x;
		}
		for ( i = 0; i < maxH; ++ i)
			if (da[i] != -100000 && xiao[i] != 100000) break;
		int head = i;
		for ( i = maxH-1; i >= 0; -- i)
			if (da[i] != -100000 && xiao[i] != 100000) break;
		int tail = i;
		
		for ( i = head+1; i <= tail; ++ i){
			if (da[i] == -100000 || da[i-1]-da[i] > 20) da[i] = da[i-1];
			if (xiao[i] == 100000 || xiao[i]-xiao[i-1] > 20) xiao[i] = xiao[i-1];
		}
		
		for ( i = 0; i < maxH; ++ i){
			double r = (double)(da[i]-xiao[i]) / 2.0;
			if (da[i] == -100000 || xiao[i] == 100000) continue;
			for ( j = xiao[i]+1; j < da[i]; ++ j){
				
				int t = j-xiao[i];
				rmap[i][j+maxW] = sqrt((double)r*r-(r-t)*(r-t))+500;
				++ TotMapNumber;
				
			}
		}
	}
}

branch* GenerateBran(int index){
	branch *pb = new branch();
	pb->i = riTable[rand() % riSize];
	pb->a = raTable[rand() % raSize];
	int ttt = rand() % 11;
	double tratio = LratioLower + (double)ttt / 100.0;
	if (vbr[index].l < 100) tratio = 0.8;
	if (vbr[index].l < 50) tratio = 0.85;
	//if (vbr[index].l < 30) tratio = 0.95;
	pb->l = vbr[index].l * tratio;
	pb->r = vbr[index].r * Rratio;
	pb->left = -1;
	pb->right = -1;
	return pb;
}

double max(double x,double y){
	if (x > y) return x;
	return y;
}

double min(double x,double y){
	if (x < y) return x;
	return y;
}

//pending whether the branch is in the tree crown
bool inborder(branch* pb,glm::mat4 modelview,int index){
	glm::vec4 test(0,0,0,1);
	glm::vec4 tmpvec(0,0,0,1);
	tmpvec = modelview * tmpvec;
	modelview = glm::rotate(modelview,(float)pb->a,ZAxis);
	modelview = glm::rotate(modelview,(float)pb->i,YAxis);
	modelview = glm::translate(modelview,glm::vec3(0.f,0.f,(float)pb->l));
	test = modelview * test;
	double tx = test.x;
	double ty = test.z;
	pb->x = test.x;
	pb->y = test.y;
	pb->z = test.z;
	//cout << "inborder:" << index << " " << pb->x << " " << pb->y << " " << pb->z << " " << abs(pb->y-tmpvec.y) << " " << rmap[(int)pb->z][(int)pb->x+maxW] << endl;
	if (abs(pb->y-tmpvec.y) <= rmap[(int)pb->z][(int)pb->x+maxW]) return true;
	else{
		return false;
	}
}


//generate random branches
void GetRandomBran(int index,glm::mat4 modelview){
	int cur,SonNum,tmp;
	branch* pbtest;
	if (vbr[index].l <= LeastbranchLength) return;
	//cout << index << endl;
	modelview = glm::rotate(modelview,(float)vbr[index].a,ZAxis);
	modelview = glm::rotate(modelview,(float)vbr[index].i,YAxis);
	modelview = glm::translate(modelview,glm::vec3(0.f,0.f,(float)vbr[index].l));
	
	
	if (vbr[index].left != -1){
		cur = vbr[index].left;
		while (cur != -1){
			if (ifexp[cur] && vbr[cur].l >= LeastbranchLength) GetRandomBran(cur,modelview);
			cur = vbr[cur].right;
		}
	}
	else{
		SonNum = rnTable[rand() % rnSize];
		if (SonNum == 1){
			pbtest = GenerateBran(index);
			pbtest->a += (rand() % 4) * 90;
			pbtest->p = index;
			ifexp.push_back(inborder(pbtest,modelview,index));
				vbr.push_back(*pbtest);
				tmp = vbr.size() - 1;
				vbr[index].left = tmp;
			
		}
		if (SonNum == 2){
			pbtest = GenerateBran(index);
			pbtest->p = index;
			ifexp.push_back(inborder(pbtest,modelview,index));
				vbr.push_back(*pbtest);
				tmp = vbr.size() - 1;
				vbr[index].left = tmp;
			
			pbtest = GenerateBran(index);
			pbtest->a += 180;
			pbtest->p = index;
			ifexp.push_back(inborder(pbtest,modelview,index));
				vbr.push_back(*pbtest);
				tmp = vbr.size() - 1;
				vbr[tmp].right = vbr[index].left;
				vbr[index].left = tmp;
			
		}
		if (SonNum == 3){
			pbtest = GenerateBran(index);
			pbtest->p = index;
			ifexp.push_back(inborder(pbtest,modelview,index));
				vbr.push_back(*pbtest);
				tmp = vbr.size() - 1;
				vbr[index].left = tmp;
			
			pbtest = GenerateBran(index);
			pbtest->a += 120;
			pbtest->p = index;
			ifexp.push_back(inborder(pbtest,modelview,index));
				vbr.push_back(*pbtest);
				tmp = vbr.size() - 1;
				vbr[tmp].right = vbr[index].left;
				vbr[index].left = tmp;
			
			pbtest = GenerateBran(index);
			pbtest->a += 240;
			pbtest->p = index;
			ifexp.push_back(inborder(pbtest,modelview,index));
				vbr.push_back(*pbtest);
				tmp = vbr.size() - 1;
				vbr[tmp].right = vbr[index].left;
				vbr[index].left = tmp;
				
			
		}
		if (SonNum == 4){
			pbtest = GenerateBran(index);
			pbtest->p = index;
			ifexp.push_back(inborder(pbtest,modelview,index));
				vbr.push_back(*pbtest);
				tmp = vbr.size() - 1;
				vbr[index].left = tmp;
			
			pbtest = GenerateBran(index);
			pbtest->a += 90;
			pbtest->p = index;
			ifexp.push_back(inborder(pbtest,modelview,index));
				vbr.push_back(*pbtest);
				tmp = vbr.size() - 1;
				vbr[tmp].right = vbr[index].left;
				vbr[index].left = tmp;
			
			pbtest = GenerateBran(index);
			pbtest->a += 180;
			pbtest->p = index;
			ifexp.push_back(inborder(pbtest,modelview,index));
				vbr.push_back(*pbtest);
				tmp = vbr.size() - 1;
				vbr[tmp].right = vbr[index].left;
				vbr[index].left = tmp;
			
			pbtest = GenerateBran(index);
			pbtest->a += 270;
			pbtest->p = index;
			ifexp.push_back(inborder(pbtest,modelview,index));
				vbr.push_back(*pbtest);
				tmp = vbr.size() - 1;
				vbr[tmp].right = vbr[index].left;
				vbr[index].left = tmp;
			
		}
		cur = vbr[index].left;
		while (cur != -1){
			if (ifexp[cur] && vbr[cur].l >= LeastbranchLength) GetRandomBran(cur,modelview);
			cur = vbr[cur].right;
		}
	}
}



// find tree tips
void Findtips(int index){
	int cur;
	cur = vbr[index].left;
	while (cur != -1){
		if (ifexp[cur]){
			vbr[index].iftip = false;
			Findtips(cur);
		}
		cur = vbr[cur].right;
	}
		
}




void RenderLeaves( int len )
{
	
	glPushMatrix();
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	double i = 0; 
	while(i < len){
		glColor3f(Lcr/(float)255.0, Lcg/(float)255.0, Lcb/(float)255.0);		
		glPushMatrix();
		glRotatef(45, 1.0, 0.0, 0.0);
		glTranslatef(LEAF_STEP, 0, 0);
		glBegin(GL_QUADS);
		glVertex3f(-LEAF_SIZE, LEAF_SIZE, 0.0f);	
		glVertex3f( LEAF_SIZE, LEAF_SIZE, 0.0f);	
		glVertex3f( LEAF_SIZE,-LEAF_SIZE, 0.0f);	
		glVertex3f(-LEAF_SIZE,-LEAF_SIZE, 0.0f);
		glEnd();	
		glPopMatrix();
		glRotatef(60, 0.0, 0, 1.0);
		glTranslatef(0,0,(LEAF_STEP));
		i += LEAF_STEP;
	}
	glPopAttrib();
	glPopMatrix();

}

void OnDrawCylinder(int index,int cen)
{
   int cur;
   if (cen > 25) return;
   //if (index > 17 || index == 3) return;
   GLUquadricObj* pQuadObj = gluNewQuadric();
   gluQuadricNormals(pQuadObj, GLU_SMOOTH);	
   //glEnable(GL_TEXTURE_2D);
   if (index == curbranch){
	   glColor3f(1.0,0,0);
   }
   else
	   glColor3f(Tcr/(float)255.0, Tcg/(float)255.0, Tcb/(float)255.0);
   gluQuadricDrawStyle(pQuadObj, GLU_FILL);
   gluQuadricTexture(pQuadObj, GL_FALSE);
	//glBindTexture(GL_TEXTURE_2D, m_Texture[0]);
   glRotatef(vbr[index].a,0.0f,0.0f,1.0f);
   glRotatef(vbr[index].i,0.0f,1.0f,0.0f);
   gluCylinder(pQuadObj, vbr[index].r, vbr[index].r*0.6, vbr[index].l, 40, 5);
   glTranslatef(0,0,vbr[index].l);
   cur = vbr[index].left;
   if (cur == -1){
	   glTranslatef(0,0,-vbr[index].l);
	   if (ifleaf == 1) RenderLeaves(vbr[index].l);
	   glTranslatef(0,0,vbr[index].l);
   }
   while (cur != -1){
		glPushMatrix();
			OnDrawCylinder(cur,cen+1);
		glPopMatrix();
		cur = vbr[cur].right;
   }
}



    

void RenderScene()
{
	glEnable(GL_DEPTH_TEST);
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);	

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0,eye,eye,0.0,eye,0.0,0.0,1.0,0.0);

	//glTranslatef(-pBranSkltn->node.x,-pBranSkltn->node.y,0.0f);
	glScalef(0.015f+scaleX,0.015f+scaleY,0.015f+scaleX);
	glRotatef(angleX,0.0f,1.0f,0.0f);
	glRotatef(angleY,1.0f,0.0f,0.0f);
	glRotatef(angleZ,0.0f,0.0f,1.0f);
	glLineWidth(3.0f);
	//distance field
	glTranslatef(0,-80,0);
	glBegin(GL_POLYGON);
		glColor3f(1.0f,1.0f,1.0f);
		glVertex3f(-64.0f,0.0f,-64.0f);
		glVertex3f(64.0f,0.0f,-64.0f);
		glVertex3f(64.0f,0.0f,64.0f);
		glVertex3f(-64.0f,0.0f,64.0f);
	glEnd();
	//draw coordinate system
	glBegin(GL_LINES);
		//x-axis
		glColor3f(1.0f,0.0f,0.0f);
		glVertex3f(0.0f,0.0f,0.0f);
		glVertex3f(50.0f,0.0f,0.0f);
		//y-axis
		glColor3f(0.0f,1.0f,0.0f);
		glVertex3f(0.0f,0.0f,0.0f);
		glVertex3f(0.0f,50.0f,0.0f);
		//z-axis
		glColor3f(0.0f,0.0f,1.0f);
		glVertex3f(0.0f,0.0f,0.0f);
		glVertex3f(0.0f,0.0f,50.0f);
	glEnd();
	
    //draw the branch
    /*
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);	
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glTexImage2D(GL_TEXTURE_2D, 0, 3, TexInfo->bmiHeader.biWidth, TexInfo->bmiHeader.biHeight,
                0, GL_RGB, GL_UNSIGNED_BYTE, TexBits);
    glEnable(GL_TEXTURE_2D);
    */
    //Load_Texture();
    //glColor3f(1.0f,1.0f,1.0f);
	//glFlush();
	glRotatef(-90.0,1.0f,0.0f,0.0f);
	OnDrawCylinder(1,1);
	glutPostRedisplay();
	glutSwapBuffers();
}

void SetupRC()
{
	glClearColor(0.0f,0.0f,0.0f,1.0f);
	glLineWidth(3.0f);
}
void ChangeSize(GLsizei w,GLsizei h)
{
	if(h == 0)
		h = 1;
	//visible area
	glViewport(0,0,w,h);
	//projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if(w <= h)
	{
		glOrtho(-projectWidth,projectWidth,-projectWidth*h/w,projectWidth*h/w,1,2000);
	}
	else
	{
		glOrtho(-projectWidth*w/h,projectWidth*w/h,-projectWidth,projectWidth,1,2000);
	}
	glMatrixMode(GL_MODELVIEW);
}
void MouseHandler(int button,int state,int x,int y)
{
	if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		//MessageBeep(-1);
		mouseX = x;
		mouseY = y;
	}
	//the wheel of mouse to zoom in and out the tree
	if(button == GLUT_WHEEL_UP)
	{
		printf("the wheel of mouse is up.\n");
		scaleX += 0.0025;
		scaleY += 0.0025;
		scaleZ += 0.0025;
	}
	else if(button == GLUT_WHEEL_DOWN)
	{
		//printf("the wheel of mouse is down.\n");
		scaleX -= 0.0025;
		scaleY -= 0.0025;
		scaleZ -= 0.0025;
	}
}

void MouseMoveHandler(int x,int y)
{
	//MessageBeep(-1);
	angleX += x - mouseX;
	angleY += y - mouseY;
	angleZ += y - mouseY;
//	printf("angleX=%d angleY=%d\n",angleX,angleY);
	mouseX = x;
	mouseY = y;
}
				

double GetAngle(tdpoint aindex){
	double answ = asin(aindex.y / sqrt(aindex.x*aindex.x+aindex.y*aindex.y));
	if (aindex.x <= 0 && aindex.y >= 0) answ = 180.0 - answ;
	if (aindex.x <= 0 && aindex.y <= 0) answ = 180.0 + answ;
	if (aindex.x >= 0 && aindex.y <= 0) answ = 360.0 - answ;
	return answ;
}

bool direc(tdpoint x,tdpoint y){
	double tx = GetAngle(x);
	double ty = GetAngle(y);
	return tx > ty;
}

// main branch location adjustment
void DoLocation(int index,glm::mat4 mm){
	glm::mat4 tmm = mm;
	
	tmm = glm::rotate(tmm,(float)vbr[index].a,ZAxis);
	tmm = glm::rotate(tmm,(float)vbr[index].i,YAxis);
	tmm = glm::translate(tmm,glm::vec3(0.f,0.f,(float)vbr[index].l));

	glm::vec4 test(0,0,0,1);
	test = tmm * test;
	if (index > 1){
		int itmp = skltnpts[index].parent;
		int fitmp = skltnpts[itmp].parent;
		if (direc(tdpoint(vbr[index].x - vbr[itmp].x,vbr[index].z - vbr[itmp].z),
				  tdpoint(vbr[itmp].x - vbr[fitmp].x,vbr[itmp].z - vbr[fitmp].z)) != 
			direc(tdpoint(test.x - vbr[itmp].x,test.z - vbr[itmp].z),
				  tdpoint(vbr[itmp].x - vbr[fitmp].x,vbr[itmp].z - vbr[fitmp].z))) vbr[index].a = 180.0 - vbr[index].a;
		if (vbr[index].a < 0) vbr[index].a += 360.0;
		if ((vbr[index].y-vbr[itmp].y) * (test.y - vbr[itmp].y) < 0) vbr[index].a = 360.0 - vbr[index].a;
	}
	mm = glm::rotate(mm,(float)vbr[index].a,ZAxis);
	mm = glm::rotate(mm,(float)vbr[index].i,YAxis);
	mm = glm::translate(mm,glm::vec3(0.f,0.f,(float)vbr[index].l));
	cout << index << " " << test.x << " " << test.y << " " << test.z << " " << vbr[index].x << " " << vbr[index].y << " " << vbr[index].z << endl;
	int cur = vbr[index].left;
	while (cur != -1){
		DoLocation(cur,mm);
		cur = vbr[cur].right;
	}
}

void GenerateLocation(){
	glm::mat4 mm = glm::scale(1.f,1.f,1.f);
	DoLocation(1,mm);
}



void KeyBoardHandler(unsigned char key, int x, int y)
{
	if (key == 'l'){
		ifleaf = -ifleaf;
	}
	if (key == 'b'){
		curbranch = 0;
	}
	if (key == 'n') curbranch = 1;
	if (key == 'a'){
		if (curbranch != 1){
			if (vbr[vbr[curbranch].p].left != curbranch){
				int ttttt = vbr[vbr[curbranch].p].left;
				while (vbr[ttttt].right != curbranch){
					ttttt = vbr[ttttt].right;
				}
				curbranch = ttttt;
			}
		}
	}
	if (key == 'd'){
		if (vbr[curbranch].right != -1) curbranch = vbr[curbranch].right;
	}
	if (key == 'w'){
		if (curbranch != 1) curbranch = vbr[curbranch].p;
	}
	if (key == 's'){
		if (vbr[curbranch].left != -1) curbranch = vbr[curbranch].left;
	}
	if (key == 'z'){
		vbr[curbranch].i += 5;
		if (vbr[curbranch].i > 180) vbr[curbranch].i = 180;
		glm::mat4 mmm = glm::scale(1.f,1.f,1.f);
		
		mmm = glm::scale(1.f,1.f,1.f);
		GetRandomBran(1,mmm);
	}
	if (key == 'x'){
		vbr[curbranch].i -= 5;
		if (vbr[curbranch].i < 0) vbr[curbranch].i = 0;
		glm::mat4 mmm = glm::scale(1.f,1.f,1.f);
		
		mmm = glm::scale(1.f,1.f,1.f);
		GetRandomBran(1,mmm);
		
	}
	if (key == 'c'){
		vbr[curbranch].a += 10;
		if (vbr[curbranch].a >= 360) vbr[curbranch].a -= 360;
		glm::mat4 mmm = glm::scale(1.f,1.f,1.f);
		
		mmm = glm::scale(1.f,1.f,1.f);
		GetRandomBran(1,mmm);
		
	}
	if (key == 'v'){
		vbr[curbranch].a -= 10;
		if (vbr[curbranch].a < 0) vbr[curbranch].a += 360;
		glm::mat4 mmm = glm::scale(1.f,1.f,1.f);
		
		mmm = glm::scale(1.f,1.f,1.f);
		GetRandomBran(1,mmm);
		
	}
	
	if (key == 'p'){
		angleX = 0;
		angleY = 0;
		angleZ = 0;
	}
	
	
	if (key == 'q'){
		ofstream fout("t.txt");
		fout << Tcr << " " << Tcg << " " << Tcb << endl;
		fout << Lcr << " " << Lcg << " " << Lcb << endl;
		fout << vbr.size() << endl;
		for (int i = 0; i < vbr.size(); ++ i){
			fout << vbr[i].a << " " << vbr[i].i << " " << vbr[i].iftip << " " << vbr[i].l << " " << vbr[i].left << " " << vbr[i].p << " " << vbr[i].r << " "
			<< vbr[i].right << " " << vbr[i].x << " " << vbr[i].y << " " << vbr[i].z << endl;
		}
		for ( int i = 0; i < maxH; ++ i){
			for ( int j = 0; j < maxW; ++ j)
				fout << rmap[i][j] << " ";
			fout << endl;
		}
		fout.close();
	}
	if (key == 'r'){
		ifstream fin("t.txt");
		int branchNum;
		fin >> Tcr >> Tcg >> Tcb >> Lcr >> Lcg >> Lcb >> branchNum;
		vbr.clear();
		for ( int i = 0; i < branchNum; ++ i){
			branch pb;
			fin >> pb.a >> pb.i >> pb.iftip >> pb.l >> pb.left >> pb.p >> pb.r
			>> pb.right >> pb.x >> pb.y >> pb.z;
			vbr.push_back(pb);
		}
		for ( int i = 0; i < maxH; ++ i){
			for ( int j = 0; j < maxW; ++ j)
				fin >> rmap[i][j];
		}
		fin.close();
	}
	if (key == 'o'){
		freopen("c.txt","w",stdout);
		cout << vbr.size() << endl;
		for ( int i = 0; i < vbr.size(); ++ i)
			cout << i << " " << vbr[i].l << " " << vbr[i].a << " " << vbr[i].i << " " << vbr[i].r << " " << vbr[i].left << " " << vbr[i].right << endl;
	}
	
	if (key == 'm'){
		angleX += 45;
	}
	if (key == 't'){
		vbr[curbranch].a += 180;
	}
	if (key == 'g'){
		glm::mat4 mmm = glm::scale(1.f,1.f,1.f);
		GetRandomBran(1,mmm);
	}
	
}



int main(int argc, char *argv[]){
	image = imread(picstr,CV_LOAD_IMAGE_COLOR);
	img = cvLoadImage(picstr.c_str());
	namedWindow("Display window",CV_WINDOW_AUTOSIZE);
	imshow("Display window",image);
	setMouseCallback("Display window",MousePoint,0);
	srand( (unsigned)time(NULL));
	for (;;){
		char c = waitKey(0);
		if (c == 'o'){
			freopen("a.txt","w",stdout);
			/*
			cout << crown.size() << endl;
			for ( int i = 0; i < vbr.size(); ++ i)
				cout << i << " " << vbr[i].l << " " << vbr[i].a << " " << vbr[i].i << " " << vbr[i].r << " " << vbr[i].left << " " << vbr[i].right << endl;
			*/
			//cout << "Evaluation: " << counter[0] * 10000 + counter[1] * 100 + counter[2] << " Base Rank: " << MinRank << endl;
			
		}
		
		if (c == 's'){
			
			vbr.clear();
			ifexp.clear();
			glm::mat4 modelview = glm::scale(1.f,1.f,1.f);
			memset(eval,-1,sizeof(eval));
			LratioLower = 0.7;
			LratioUpper = 0.8;
			GetOriginBran();
			LratioLower = 0.75;
			LratioUpper = 0.85;
			//GetRandomBran(1,modelview);
			cout << "branch Number:" << vbr.size() << endl;
			GenerateLocation();
			GenerateLocation();
			
			TexBits = LoadDIBitmap("tree.bmp",&TexInfo);
			glutInit(&argc,argv);
			glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
			glutInitWindowSize(780,840);
			glutInitWindowPosition(200,200);
			glutCreateWindow("3D tree branch!");
			//printf("Current Version: %s.\n",glGetString(GL_VERSION));
			glutDisplayFunc(RenderScene);
			glutReshapeFunc(ChangeSize);
			glutMouseFunc(MouseHandler);
			glutMotionFunc(MouseMoveHandler);
			glutKeyboardFunc(KeyBoardHandler);
			SetupRC();
			glutMainLoop();
		}
		if (c == 'e'){
			setMouseCallback("Display window",CrownPoint,0);
		}
		
            
	}
	return 0;
}
