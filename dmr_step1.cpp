#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <algorithm>
#include <fcntl.h>
#include <map>
#include <math.h>

#define maxN 20000000//28000
#define sWinN 1000 //300 //1000000 //500 //249 //1000000
#define maxLimit 300
#define chromoNo 7
#define maxLarge 56000 //60000 ,number of proteins
#define maxNumFea 30000 //600000
#define Half 8000
#define trainPairs 14
#define M0 4  //4
#define M1 10  //10  positions
#define M2 10  //10 intervals, jumping distance

//#define PNclass "+1"
using namespace std;
//int pairRedund[maxLarge][2];
//int pairNoR[maxLarge][2];

//unsigned m,way,actualway;
string feaArray[maxNumFea];
bool flagDomainA[maxLarge];
bool flagFeature[maxNumFea];
bool flagFeature2[maxNumFea];

long int totalF=0;
long int totalF_ys=0;
long int totalLine=0,totalE1=0,totalE2=0;
long double threshV=0.01, eValueOri=1.0, nLE_eValueOri;
vector <long int>arrayMethy[chromoNo];
vector <long int>arrayMethy_2[chromoNo];

//ofstream fileOut("irregular format.txt");
//ofstream fileOutE("Evalue larger than threshV.txt");
//ofstream fileOutE2("Evalue is NA.txt");
//ofstream coutTest("test signT.txt"); //the sum of lines of all files for one species
//string pfamFeature[maxNumFea];
//string proteinPfam[maxN][3];

int TOUPPER(int c)   
{  
    return   toupper(c);  
}
struct feaNode
{
	int No;
	long double negaLogE;
	int realN;
	int biN;
//}interProArray[maxLarge][maxNumFea];
}tmpNode;
vector <feaNode> interProArray[maxLarge];

struct positionNoNode
{
	long int pos;
	long int end;
	long double pV;
	long int num;// hyper or hypo
	long int num2;//λ������ 
	long int numCom;
	int markR;
	int DMR_S;
	int DMR_E;
	vector <long int> posV;
}tmpNode2,tmpNode3,tmpNode4;
vector <positionNoNode> arrayMethy1[chromoNo];
vector <positionNoNode> arrayMethy2[chromoNo];
vector <positionNoNode> arrayMethy3[chromoNo];

bool ascendSort(const feaNode & f1, const feaNode & f2) 
{
    if (f1.No < f2.No) return true;
	else return false;
}
bool cmpV(const positionNoNode & f1, const positionNoNode & f2) 
{
    if (f1.pos < f2.pos) return true;
	else return false;
}
bool cmpM(long int a, long int b)
{
	if (a<b) return true;
	else return false;
}
bool cmpCom(const positionNoNode & f1, const positionNoNode & f2) 
{
    if (f1.numCom > f2.numCom) return true;
	else return false;
}
int lineNo(ifstream &file)
{
    string line;
    int count = 0;
    while (getline(file, line))
    {
        count++;  //���������⣬���ܶ��һ��
    }

    file.clear();
    file.seekg(0, ios::beg);

    return count;
}
/*
struct pairNode
{
	vector<string> LV;
	vector<string> RV;
}pairArray[maxLarge];
*/
char *removelast(char *p)
{
	if(p==NULL) return NULL;
	if(strlen(p)==0) return "";
	char *pnew = (char*)malloc((strlen(p)+1)*sizeof(char));
	strcpy(pnew,p);
	pnew[strlen(pnew)-1]=0;
	return pnew;
}
/*
long double changeM(string a)
{
	int posE=0, flagS=0;
	string partOne,partTwo; partThree;
	flag=0;
	for (int i=0; i<a.length(); i++)
	{
		if (a[i]=='e' || a[i]=='E')
		{
			if (a[i+1]=='+')
			{
				posE=i;
				flagS=1;
				break;
			}
			else if (a[i+1]=='-')
			{
				posE=i;
				flagS=-1;
				break;
			}
			//cout<<"xxxe-xxx error!"<<endl;
		}		
	}
	if (posE) 
	{
		partOne=a.substr(0,posE);
		partTwo=a.substr(posE+2, a.length()-posE-1);
		flag=1;  // 1 means 'e' is included
	}
	else partOne=a;
	
}
*/
long double changeR(string a)//�ַ���תΪ˫���ȸ����� 
{
    long double b=0.0;
    int i,f,s,s1,expN,flag,posE;  //f=sign   
    f=1;
	s=-1;  //s=-1:integer part   s=1:decimal part
	s1=0;  //s1  '+' or '-' after 'e'
	expN=0;
	flag=0;
	i=0;
	if(a[0]=='-')
	{
		f=-1;
		i=1;
	}
	for ( ; i<a.length(); i++)
    {
		if((a[i]=='e') || (a[i]=='E')) {flag=1; posE=i; break;}
		else
		{
			if(a[i]=='.') 
			{
				s=1;
				continue;
			}
			//else if(s<0) b=b*10+(*p-48);
			if(s<0) b=b*10.0+(long double)(a[i]-48);
			else
			{
				//if(*p=='.')p++;
				s=s*10;
				//cout<<(float)(a[i]-48)<<endl;
				//cout<<1.0/4.0<<endl;
				//cout<<((long double)(a[i]-48))/(long double)s<<endl;
				b=b+(long double)(a[i]-48)/(long double)s;
				//cout<<b<<endl;
			}
		}
	}
	i++;
    b=b*f;
	if (flag)
	{
		if (a[i]=='-') s1=-1;
		else if (a[i]=='+') s1=1;
		else 
		{
			cout<<"no sign after 'e'!"<<endl;
		}		
		
		for (int k=posE+2; k<a.length(); k++)
		{
			expN=expN*10+(a[k]-48);
		}
		if (s1>0)
		{
			for (int r1=0; r1<expN; r1++)
				b=b*(long double)(10.0);
		}
		else if (s1<0)
		{
			for (int r1=0; r1<expN; r1++)
				b=b/(long double)(10.0);
		}
	}
    return b;
} 
string subSoyName(string str)
{
	string str1="";
	int strLen,cp=0,pos1=0,pos2=0;
	strLen=str.length();
	
	for (int i=0; i<strLen; i++)
	{
		if (str[i]=='|') 
		{
			cp++;
			if (cp==1) pos1=i;
			else if (cp==2) 
			{
				pos2=i;
				break;
			}
		}
	}
	str1=str.substr(pos1+1,pos2-pos1-1);
	return str1;
}
//int main()
int main(int argc, char * argv[])
{
    /*
	string str,str1,str2, str3,str4,str5,amiAci[20],secStr[3],tmpFea,tmp,tmp1,typeInt,firstP,secondP,strTmp,strOri;
	string str6,str7,str8,str9,str10,str11,leftP,rightP,strNewTmp,strTmp2,strN1,strN3,strN1Tmp;
	//double leftArray[Half+1], rightArray[Half+1];
	int number,number1,number2,number3,number4,number5, numberApp, leftNo,rightNo,flag, flag1,flag2,flag3,flagNew,flagNew1,flagNew2,countColu;
	int sign, No, cir, cirNew4, afterFirst, afterSecond,cirNew6,lenStr, lenStr2, lenThree, lenFour,cirPA;
	double p;
	long double valueTmp;
	long int totalTest;
	char ch, dic[80];
	int Te,Tr,cir1,cirNew,cirNew2,cirTmp,SGDstart,positionPair1,positionPair2,leftPosition, rightPosition, leftPositionNew, rightPositionNew, cirNew3;
	int NoTmp,flagDnum,cirInt,cirRealPro,cirSoyF,position3, position4, position5, position6,featureP, numberF, pNumber, featureN, sgdNumber, NoYP;
	string strIni, threeStr, fourStr, str001, str002 ,str003, str004, str005, str006, str007, str008, eValueTmp, strYfea;
	string strm1,strm2,strm3,strm4,strm5,strGiven,strPath,strName,strName1,strName2;
	bool flagExtraFeaA[maxLarge];
	*/
	//cir1=0;
	/*
	if (0.0<1e-323) cout<<"OK!"<<endl;
	else cout<<"No!"<<endl;
	*/
	//string noToNameArray[maxN],newStrArray[maxN];
	string strX,strX1,strX0,str,str0,str1,str2,str3,str4,str5,str6,str7,str8,str10,strEnd,strName1,strName2,strName3,strName4;
	string str01,str02,str03,str20,str21,str30,str31,str32,str33,str34,str35,str36,strFile1,strFile2,strFile3,strFile4,strFile5,strFile6,strFile7,strFile8,strFile9;
	int flagX1,cir1,cir2,Order,NoP1,strN1,strN2,flag,orderNo,flag0,flag1,flag2,flag3,flag4,flag5,flag6,flag7,num10,cir20;
	long int startP,endP,cirI1,cirI2,cirI3,numS1,numS2,numS3,cirN,cirN0,cirN1,valueM,methyP, numTmp,numTmp1,cir,cirTmp,posTmp,posTmp1,posTmp2,numTmp2,numTmp3,numTmp4,tmpV1,tmpV2;
	long int cir10,num0,num1,num2,num3,num4,firstP,lastP,maxCom,num01,num02,num03,num04,num05,num06,num07,num08,num09,numTmp01,numTmp02,numTmp03,numTmp5,numTmp09;
	typedef map<string, int> strToNumDefine;
	strToNumDefine speciesA;
	map<int, int> filterP;
	map<long int, int>posVmap;
	char dic[50],charS1,charS2;
	long int arrayCir[7], arrayCir_2[7], arrayCir_3[7];
	long double pValue,pValueTmp,ratioI;
	vector<long int>posV;
	//int M0;
	memset(dic,'0',sizeof(dic));
	memset(arrayCir, 0, sizeof(arrayCir));
	memset(arrayCir_2, 0, sizeof(arrayCir_2));
	memset(arrayCir_3, 0, sizeof(arrayCir_3));
	
	
	
	//str01="msh1_col_FET_smaller0.01.txt";
	//ifstream fin1(str01.c_str());
	//ofstream cout2("new output_msh1_col_FET_smaller0.01.txt");
	//ofstream cout3("slidingW_msh1_col_FET_smaller0.01.txt");
	
	//ifstream fin01("slidingW_msh1_col_FET_smaller0.01.txt");
	//ofstream cout02("allDMCs_new_Standardized_slidingW_msh1_col_FET_smaller0.01.txt");
	//ofstream cout03("noTitle_allDMCs_new_Standardized_slidingW_msh1_col_FET_smaller0.01.txt");
	
	str01=argv[1];
	//str01="file1.txt";
	
	ifstream fin1(str01.c_str());
	
	strFile1="newOutput_"+str01;
	strFile2="slidingW_"+str01;
	strFile3="allDMCs_new_Standardized_slidingW_"+str01;
	strFile4="noTitle_allDMCs_new_Standardized_slidingW_"+str01;
	strFile5="runningOutput_"+str01;
	strFile6="DMR_"+str01;
	strFile7="DMR_questionNoOverlap_"+str01;
	strFile8="unitWin_data_"+str01;
	strFile9="DMR_noOverlap_"+str01;
	
	ofstream cout2(strFile1.c_str());//�����λ������ 
	ofstream cout3(strFile2.c_str());
	
	ifstream fin01(strFile2.c_str());
	ofstream cout02(strFile3.c_str());
	ofstream cout03(strFile4.c_str());
	ofstream cout04(strFile5.c_str());//������־ 
	ofstream cout05(strFile6.c_str());
	ofstream cout06(strFile7.c_str());
	ofstream cout07(strFile8.c_str());
	ofstream cout08(strFile9.c_str());
	
	ratioI=0.1;
	//M2=((double)(1)/(double)ratioI);
	//M2=10/1;
	cout<<"jumping gap is:"<<M2<<endl;
	
	//section1
	//start  reading in the file, str01, sorting, and outputing
	getline(fin1,str02);               
	cirN1=0;
	numS1=0;
	numS2=0;
	numS3=0;
	//cir2=0;
	str="";
	while(getline(fin1,str))
	{
		//cout<<str<<endl;
		//cir1++;
		istringstream instr(str);
			
		instr>>num1>>pValue>>num2;
		tmpNode2.pos=num1;
		tmpNode2.pV=pValue;
		tmpNode2.num=num2; //.num is the mark of methy or unmethy
		arrayMethy1[0].push_back(tmpNode2);	
	}
	
	
	if (arrayMethy1[0].empty())  
	{
		cout04<<"the arrayMethy1[0] is empty, check!"<<endl;
	}
	sort(arrayMethy1[0].begin(),arrayMethy1[0].end(),cmpV);//����λ���������� 
	
	vector<positionNoNode>::iterator it0;
	vector<positionNoNode>::iterator it2;
	it0=arrayMethy1[0].end();
	cout04<<"the end of arrayMethy1 +1	"<<it0->pos<<"	"<<it0->num<<"	"<<it0->num2<<endl;
	it0--;
	cout04<<"the end of arrayMethy1	"<<it0->pos<<"	"<<it0->num<<"	"<<it0->num2<<endl;
	flag0=1;
		
	if (!arrayMethy1[0].empty()) it0=arrayMethy1[0].begin();
	else 
	{
		flag0=0;
		cout04<<"the arrayMethy1[0] is empty, check!"<<endl;
	}
	
	cir=0;
	while (flag0 && (it0!=arrayMethy1[0].end()))
	{
		cir++;
		numTmp1=it0->pos;
		pValueTmp=it0->pV;
		numTmp2=it0->num;
		cout2<<numTmp1<<"	"<<pValueTmp<<"	"<<numTmp2<<endl;
		it0++;
	}
	
	if (it0==arrayMethy1[0].end())
	{
		it0--;
		numTmp1=it0->pos;
		lastP=numTmp1; //lastP is the maximal position in the current file, fin1
		cout04<<"lastP="<<lastP<<endl;
	}
	else cout<<"Wait, the lastP hasn't been decided, Check!"<<endl;
	
	it0=arrayMethy1[0].begin();
	numTmp1=it0->pos;
	firstP=-1;
	firstP=numTmp1;
	if (firstP>0) cout04<<"firstP="<<firstP<<endl;
	else cout<<"Wait, the firstP hasn't been decided, Check!"<<endl;
	//end of read, sort, ...
	
	
	
	
	//section2
	//start calculating the number of positions within each of the unitSlidingWindows
	
	//ratioI=0.0001;
	//ratioI=0.001;
	//ratioI=0.1;
	
	it0=arrayMethy1[0].begin();
	num1=it0->pos;
	cirN0=(num1-1)/sWinN;//�������ڴ�С��100
	cir2=0;
	vector<long int>::iterator it30;
	while (1)
	{
		numS1=0;
		numS2=0;
		startP=cirN0*sWinN+sWinN*ratioI*cir2+1;//ratio:0.1 
		endP=cirN0*sWinN+sWinN*(ratioI*cir2+1);
		
		//cout<<startP<<"	"<<endP<<endl;
		
		if (startP>lastP) break;
		it0=arrayMethy1[0].begin();
		num1=it0->pos;
		num2=it0->num;
		//cout<<num1<<"	"<<num2<<endl;		
		flag=0;
		
		posV.clear();
		
		while (1)
		{
			if ((num1>=startP) && (num1<=endP))
			{
				flag=1;
				if (num2) numS1++;
				else if (!num2) numS2++;
				posV.push_back(num1);
			}
			it0++;
			if (it0==arrayMethy1[0].end()) break;
			num1=it0->pos;
			num2=it0->num;
			if (num1>endP) break;
		}
		if (flag)
		{
			tmpNode2.pos=startP;
			tmpNode2.end=endP;
			tmpNode2.num=numS1; //methy
			tmpNode2.num2=numS2; //unMethy
			tmpNode2.posV=posV;
			//cout<<tmpNode2.pos<<"	"<<tmpNode2.end<<"	"<<tmpNode2.num<<"	"<<tmpNode2.num2<<endl;
			arrayMethy1[1].push_back(tmpNode2);	
		}
		else
		{
			tmpNode2.pos=startP;
			tmpNode2.end=endP;
			tmpNode2.num=0;
			tmpNode2.num2=0;
			tmpNode2.posV=posV;
			//cout<<tmpNode2.pos<<"	"<<tmpNode2.end<<"	"<<tmpNode2.num<<"	"<<tmpNode2.num2<<endl;
			arrayMethy1[1].push_back(tmpNode2);	
		}
		cir2++;
	}
	cout04<<endl<<endl<<"the total number of intervals: "<<cir2<<endl; //or cir2-1?
	
	vector<positionNoNode>::iterator it1;
	flag0=1;
	if (!arrayMethy1[1].empty()) it1=arrayMethy1[1].begin();
	else 
	{
		flag0=0;
		cout04<<"the arrayMethy1[1] is empty, check!"<<endl;
	}
		
	cir=0;
	while (flag0 && (it1!=arrayMethy1[1].end()))
	{
		cir++;
		numTmp1=it1->pos;
		numTmp2=it1->end;
		//pValueTmp=it1->pV;
		numTmp3=it1->num;
		numTmp4=it1->num2;
		cout3<<numTmp1<<"	"<<numTmp2<<"	"<<numTmp3<<"	"<<numTmp4<<"	";
		for (it30=it1->posV.begin(); it30!=it1->posV.end(); it30++)
		{
			numTmp09=*it30;
			cout3<<numTmp09<<"| ";
		}
		cout3<<endl;
		it1++;
	}
	cout3.close(); //important;
	//end calculating the number of positions within the unitSlidingWindows

	

	
	
	//section3
	//start copy arrayMethy1[1] to arrayMethy1[0] and arrayMethy1[2] with adding numCom
	//standardized:
	    
	arrayMethy1[0].clear();
	if (arrayMethy1[0].empty()) cout<<"Yes, correct, the vector has been cleared to be empty!"<<endl;
	else cout<<"No, incorrect, the vector has not been cleared to be empty, check!"<<endl;
	
	arrayMethy1[2].clear();
	while(getline(fin01,str))
	{
		//cout<<str<<endl;
		//cir1++;
		istringstream instr1(str);
			
		instr1>>num1>>num2>>num3>>num4;
		tmpNode2.pos=num1;
		tmpNode2.end=num2;
		tmpNode2.num=num3;
		tmpNode2.num2=num4;
		tmpNode2.numCom=num3+num4;
		tmpNode2.markR=0;
		arrayMethy1[0].push_back(tmpNode2);	
		//arrayMethy1[1].push_back(tmpNode2);
		arrayMethy1[2].push_back(tmpNode2);			
	}
	sort(arrayMethy1[0].begin(),arrayMethy1[0].end(),cmpCom);//����numcom�������� 
	//sort(arrayMethy1[1].begin(),arrayMethy1[1].end(),cmpU);
	
	vector<positionNoNode>::iterator it00;
	vector<positionNoNode>::iterator it02;
	
	it00=arrayMethy1[0].begin();
	//it1=arrayMethy1[1].begin();
	//maxO=it0->num;
	//maxU=it1->num2;
	maxCom=it00->numCom;
	cout04<<"max number of methy+unMethy in each interval SW: "<<maxCom<<endl;
	
	it02=arrayMethy1[2].begin();
	//str1="start of window";
	//str2="end of window";
	//str3="num of allDMCs";
	//str4="num of under-methy";
	//str5="standardized num_DMCs";
	//str6="standardized Under";
	
	while (it02!=arrayMethy1[2].end())
	{
		//cir++;
		numTmp1=it02->pos;
		numTmp2=it02->end;
		numTmp3=it02->numCom;
		// Tail-window fix v2: keep windows whose left boundary is still within the DMP range.
			// The original condition (numTmp2 <= lastP) dropped tail windows such as
			// 30408101-30409100 when lastP is 30409031, even though the window still
			// contains DMPs. This matches the Python behavior: stop only after startP > lastP.
			if (numTmp1<=lastP) cout03<<setiosflags(ios::left)<<setw(20)<<numTmp1<<setw(20)<<numTmp2<<setw(20)<<numTmp3<<setw(30)<<(double)numTmp3/(double)maxCom<<endl;
		it02++;
	}
	//end copy arrayMethy1[1] to arrayMethy1[0] and arrayMethy1[2] with adding numCom





    //start from here Oct25 2012, why is the output file empty?
	//actually starting at Nov07 2012, noon

	//section4
	//start DMR extension, mark the important window intervals
	vector<positionNoNode>::iterator it03;
	vector<positionNoNode>::iterator it04;
	num04=0;
	//num05=0;
	//M2=(int)((double)sWinN/((double)sWinN*ratioI));
	
	for (it02=arrayMethy1[2].begin(); it02!=arrayMethy1[2].end(); it02++)
	{
		it03=it02;
		it04=it02;
		
		num1=-1;
		num2=-1;
		num3=-1;
		
		num1=it02->pos;
		num2=it02->end;
		num3=it02->numCom;
		
		num01=-1;
		num02=-1;
		num03=-1;
			
		num01=num1;
		num02=num2;
		num03=num3;
		
		num04=0;
		//num05=0;
		
		num06=-1;
		num07=-1;
		
		num09=-1;
		num10=-1;
		
		flag1=0;
		flag2=0;
		flag3=0;
		flag4=0;
		flag5=0;
		flag6=0;
		
		while( (num03>=M0) && (num01>=firstP) )//M0Ϊ4
		{
			if (!flag3)//�������չ 
			{
				flag1=1;
				flag3=1;
				num07=num02;
				
				num09=num01;//DMR����ʼλ�� 
				num10=num02;//DMR�����λ�� 
			}
			num04=num04+num03;//�ۼ���λ���� 
			num06=num01;
			
			it03->markR=1; //???? in this part? yes
			
			if (it03!=arrayMethy1[2].begin()) 
			{
				//it03--;
				//if (it03!=arrayMethy1[2].begin()) 
				flag5=0;
				for (int k=0; k<(M2-1); k++)//��Ծ���룺m2 
				{
					//it03=it03-M2;
					it03--;
					if (it03==arrayMethy1[2].begin()) 
					{
						flag5=1;
						break;
					}
				}
				
				if (flag5)
				{
					break;
				}
				else
				{
					it03--;
					num01=it03->pos;
					//num02=it03->end;
					num03=it03->numCom;
				}
			}
			else break;
		}
		
		if (flag1)//���Ҳ���չ 
		{
			for (int k=0; k<M2; k++)
			{
				//it03=it03-M2;
				it04++;
				if (it04==arrayMethy1[2].end()) 
				{
					flag6=1;
					break;
				}
			}
			// Tail-window fix v3: if there are fewer than M2 windows to the right,
			// do NOT break the outer seed loop. The current seed/left-extension
			// can still define a valid tail boundary, matching the Python behavior.
			if (!flag6)
			{			
				while (it04!=arrayMethy1[2].end())
				{
					num02=it04->end;
					num03=it04->numCom;
					//if (it04==arrayMethy1[2].end()) break;
					
					// Tail-window fix v2: extend while the skipped-to window still has enough DMPs;
						// do not stop only because its right boundary exceeds lastP.
						if (num03>=M0)
					{
						if (!flag4)
						{
							flag2=1;
							flag4=1;
						}
						num07=num02;
						num04=num04+num03;
					
						it04->markR=1;
					
						//it04--;
						//it04=it04+M2;
						flag7=0;
						for (int k=0; k<M2; k++)
						{
							//it03=it03-M2;
							it04++;
							if (it04==arrayMethy1[2].end()) 
							{
								flag7=1;
								break;
							}
						}
						if (flag7) break;
					
						//num01=it04->pos;
						//num02=it04->end;
						//num03=it04->numCom;
					}
					else break;
				}
			}
		}
		
		if (flag1 || flag2)
		{
			tmpNode2.pos=num06;
			tmpNode2.end=num07;
			tmpNode2.numCom=num04;
			tmpNode2.DMR_S=num09;
			tmpNode2.DMR_E=num10;
			num08=num07-num06+1;
			
			//if ( (num04>=M1) && (num08>=(2*sWinN)) ) arrayMethy1[3].push_back(tmpNode2);  //do not use this condition for looking for boundaries
			arrayMethy1[3].push_back(tmpNode2); 
		}
	}
	
	for (it02=arrayMethy1[3].begin(); it02!=arrayMethy1[3].end(); it02++)
	{
		numTmp1=it02->pos;//dmr��ʼλ�� 
		numTmp2=it02->end;//dmr����λ�� 
		numTmp3=it02->numCom;
		numTmp4=it02->DMR_S;//dmr������ʼλ�� 
		numTmp5=it02->DMR_E;//dmr�������λ�� 
		cout05<<setiosflags(ios::left)<<setw(20)<<numTmp1<<setw(20)<<numTmp2<<setw(20)<<numTmp3<<"	["<<numTmp4<<" "<<numTmp5<<"]"<<endl;
	}
	//end DMR extension, mark the important window intervals
	
	
	/*
	while (it02!=arrayMethy1[2].end())
	{
		//cir++;
		numTmp1=it02->pos;
		numTmp2=it02->end;
		numTmp3=it02->numCom;
		if (numTmp2<=lastP) cout03<<setiosflags(ios::left)<<setw(20)<<numTmp1<<setw(20)<<numTmp2<<setw(20)<<numTmp3<<setw(30)<<(double)numTmp3/(double)maxCom<<endl;
		it02++;
	}
	*/
	
	
}	
	
	
