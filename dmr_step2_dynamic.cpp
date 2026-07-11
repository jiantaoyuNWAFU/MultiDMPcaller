#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <algorithm>
#include <fcntl.h>
#include <map>
#include <cmath>

#define maxN 20000000
#define sWinN 1000
#define maxLimit 300
#define chromoNo 7
#define maxLarge 56000
#define maxNumFea 30000
#define Half 8000
#define trainPairs 14
#define M0 4
#define M1 10
#define M2 10

using namespace std;

string feaArray[maxNumFea];
bool flagDomainA[maxLarge];
bool flagFeature[maxNumFea];
bool flagFeature2[maxNumFea];

long int totalF = 0;
long int totalF_ys = 0;
long int totalLine = 0, totalE1 = 0, totalE2 = 0;
long double threshV = 0.01, eValueOri = 1.0, nLE_eValueOri;

vector<long int> arrayMethy[chromoNo];
vector<long int> arrayMethy_2[chromoNo];

int TOUPPER(int c)
{
    return toupper(c);
}

struct feaNode
{
    int No;
    long double negaLogE;
    int realN;
    int biN;
} tmpNode;

vector<feaNode> interProArray[maxLarge];

struct positionNoNode
{
    long int pos;
    long int end;
    long double pV;
    long double ratio;
    long int num;
    long int num2;
    long int numCom;
    int markR;
    long int DMR_S;
    long int DMR_E;
    vector<long int> posV;
    vector<long double> logPV;
    vector<int> meUnV;
} tmpNode2, tmpNode3, tmpNode4;

vector<positionNoNode> arrayMethy1[chromoNo];
vector<positionNoNode> arrayMethy2[chromoNo];
vector<positionNoNode> arrayMethy3[chromoNo];

bool ascendSort(const feaNode &f1, const feaNode &f2)
{
    return f1.No < f2.No;
}

bool cmpV(const positionNoNode &f1, const positionNoNode &f2)
{
    return f1.pos < f2.pos;
}

bool cmpM(long int a, long int b)
{
    return a < b;
}

bool cmpCom(const positionNoNode &f1, const positionNoNode &f2)
{
    return f1.numCom > f2.numCom;
}

bool cmpCom2(const positionNoNode &f1, const positionNoNode &f2)
{
    return f1.ratio > f2.ratio;
}

int lineNo(ifstream &file)
{
    string line;
    int count = 0;
    while (getline(file, line))
    {
        count++;
    }

    file.clear();
    file.seekg(0, ios::beg);

    return count;
}

char *removelast(char *p)
{
    if (p == NULL) return NULL;
    if (strlen(p) == 0) return (char *)"";
    char *pnew = (char *)malloc((strlen(p) + 1) * sizeof(char));
    strcpy(pnew, p);
    pnew[strlen(pnew) - 1] = 0;
    return pnew;
}

long double changeR(string a)
{
    long double b = 0.0;
    int i, f, s, s1, expN, flag, posE;
    f = 1;
    s = -1;
    s1 = 0;
    expN = 0;
    flag = 0;
    posE = 0;
    i = 0;

    if (a.empty()) return 0.0;

    if (a[0] == '-')
    {
        f = -1;
        i = 1;
    }

    for (; i < (int)a.length(); i++)
    {
        if ((a[i] == 'e') || (a[i] == 'E'))
        {
            flag = 1;
            posE = i;
            break;
        }
        else
        {
            if (a[i] == '.')
            {
                s = 1;
                continue;
            }

            if (s < 0)
                b = b * 10.0 + (long double)(a[i] - 48);
            else
            {
                s = s * 10;
                b = b + (long double)(a[i] - 48) / (long double)s;
            }
        }
    }

    i++;
    b = b * f;

    if (flag)
    {
        if (i < (int)a.length() && a[i] == '-') s1 = -1;
        else if (i < (int)a.length() && a[i] == '+') s1 = 1;
        else
        {
            cout << "no sign after 'e'!" << endl;
        }

        for (int k = posE + 2; k < (int)a.length(); k++)
        {
            expN = expN * 10 + (a[k] - 48);
        }

        if (s1 > 0)
        {
            for (int r1 = 0; r1 < expN; r1++)
                b = b * (long double)(10.0);
        }
        else if (s1 < 0)
        {
            for (int r1 = 0; r1 < expN; r1++)
                b = b / (long double)(10.0);
        }
    }

    return b;
}

string subSoyName(string str)
{
    string str1 = "";
    int strLen, cp = 0, pos1 = 0, pos2 = 0;
    strLen = str.length();

    for (int i = 0; i < strLen; i++)
    {
        if (str[i] == '|')
        {
            cp++;
            if (cp == 1) pos1 = i;
            else if (cp == 2)
            {
                pos2 = i;
                break;
            }
        }
    }

    str1 = str.substr(pos1 + 1, pos2 - pos1 - 1);
    return str1;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <boundary_file> <dmp_file>" << endl;
        return 1;
    }

    string strX, strX1, strX0, str, str0, str1, str2, str3, str4, str5, str6, str7, str8, str10, strEnd, strName1, strName2, strName3, strName4, str04, str05;
    string str01, str02, str03, str20, str21, str30, str31, str32, str33, str34, str35, str36, strFile1, strFile2, strFile3, strFile4, strFile5, strFile6, strFile7, strFile8, strFile9;

    int flagX1, cir1, cir2, Order, NoP1, strN1, strN2, flag, orderNo, flag0, flag1, flag2, flag3, flag4, flag5, flag6, flag7, num10, cir20;
    long int pos1, pos2, startP, endP, cirI1, cirI2, cirI3, numS1, numS2, numS3, cirN, cirN0, cirN1, valueM, methyP, numTmp, numTmp1, cir, cirTmp, posTmp, posTmp1, posTmp2, numTmp2, numTmp3, numTmp4, tmpV1, tmpV2;
    long int maxV, caseN, cir10, num0, num1, num2, num3, num4, num8, firstP, lastP, maxCom, num01, num02, num03, num04, num05, num06, num07, num08, num09, numTmp01, numTmp02, numTmp03, numTmp5, numTmp09;

    typedef map<string, int> strToNumDefine;
    strToNumDefine speciesA;
    map<int, int> filterP;
    map<long int, int> posVmap;

    char dic[50], charS1, charS2;
    long int arrayCir[7], arrayCir_2[7], arrayCir_3[7];
    long double pValue, pValueTmp, ratioI, ratioTmp;
    vector<long int> posV;
    vector<long double> logPV;
    vector<int> meUnV;

    memset(dic, '0', sizeof(dic));
    memset(arrayCir, 0, sizeof(arrayCir));
    memset(arrayCir_2, 0, sizeof(arrayCir_2));
    memset(arrayCir_3, 0, sizeof(arrayCir_3));

    str01 = argv[1];
    str02 = argv[2];

    str03 = "error_" + str01;
    str04 = "boundaries_noOverlapping_" + str01;
    str05 = "DMR_list_" + str01;

    ofstream coutE(str03.c_str());
    ofstream cout1(str04.c_str());
    ofstream cout3(str05.c_str());

    if (!coutE)
    {
        cerr << "Error: cannot open output file: " << str03 << endl;
        return 1;
    }
    if (!cout1)
    {
        cerr << "Error: cannot open output file: " << str04 << endl;
        return 1;
    }
    if (!cout3)
    {
        cerr << "Error: cannot open output file: " << str05 << endl;
        return 1;
    }

    /*
      动态 chromoL 修改点：
      原版为：
          #define chromoL 44000000
          int chromoArray[chromoL];
      这里改为先扫描 argv[1]，找到最大 boundary end，再动态分配。
      额外 +2 是为了保留 maxV+1 位置的 0 作为连续区间扫描的哨兵。
    */
    maxV = 0;
    {
        ifstream fin_scan(str01.c_str());
        if (!fin_scan)
        {
            cerr << "Error: cannot open boundary file: " << str01 << endl;
            return 1;
        }

        string line_scan;
        while (getline(fin_scan, line_scan))
        {
            long int s = -1;
            long int e = -1;
            istringstream iss(line_scan);
            if (!(iss >> s >> e)) continue;

            if (e <= s) coutE << "error, end <= start" << endl;

            if (s > maxV) maxV = s;
            if (e > maxV) maxV = e;
        }
    }

    long int chromoL = maxV + 2;
    if (chromoL < 2) chromoL = 2;

    vector<unsigned char> chromoArray((size_t)chromoL, 0);

    ifstream fin1(str01.c_str());
    ifstream fin2(str02.c_str());

    if (!fin1)
    {
        cerr << "Error: cannot open boundary file: " << str01 << endl;
        return 1;
    }
    if (!fin2)
    {
        cerr << "Error: cannot open DMP file: " << str02 << endl;
        return 1;
    }

    ratioI = 0.1;

    pos1 = 0;
    pos2 = 0;
    str = "";

    while (getline(fin1, str))
    {
        num1 = -1;
        num2 = -1;
        caseN = 0;

        istringstream instr(str);
        if (!(instr >> num1 >> num2)) continue;

        if (num2 <= num1) coutE << "error, end <= start" << endl;

        if (num1 < 0 || num2 < 0)
        {
            coutE << "error, negative boundary: " << num1 << "\t" << num2 << endl;
            continue;
        }

        if (num2 >= chromoL)
        {
            coutE << "error, boundary exceeds dynamic chromoL: " << num1 << "\t" << num2 << "\t" << chromoL << endl;
            continue;
        }

        for (long int i = num1; i <= num2; i++)
            chromoArray[(size_t)i] = 1;
    }

    cout << maxV << "\t" << chromoL << endl;

    for (long int i = 0; i < chromoL; i++)
    {
        pos1 = 0;

        for (long int j = i; j < chromoL; j++)
        {
            if (chromoArray[(size_t)j] == 1)
            {
                pos1 = j;
                break;
            }
        }

        if (pos1 == 0) break;

        flag = 0;
        num1 = pos1;
        i = pos1 + 1;

        while (1)
        {
            if (i >= chromoL)
            {
                num2 = chromoL - 1;
                break;
            }

            if (chromoArray[(size_t)i] == 1)
            {
                flag = 1;
                i++;
            }
            else if (chromoArray[(size_t)i] == 0)
            {
                num2 = i - 1;
                break;
            }
        }

        if (flag)
        {
            tmpNode2.DMR_S = num1;
            tmpNode2.DMR_E = num2;
            arrayMethy1[0].push_back(tmpNode2);
        }
        else
        {
            coutE << "no continuous regions found for this starting point:" << pos1 << endl;
        }
    }

    vector<positionNoNode>::iterator it0;

    for (it0 = arrayMethy1[0].begin(); it0 != arrayMethy1[0].end(); it0++)
    {
        num1 = it0->DMR_S;
        num2 = it0->DMR_E;
        cout1 << num1 << "\t" << num2 << endl;
    }

    cirN1 = 0;
    str = "";

    while (getline(fin2, str))
    {
        istringstream instr001(str);

        if (!(instr001 >> num1 >> pValue >> num2))
        {
            continue;
        }

        tmpNode2.pos = num1;
        tmpNode2.pV = pValue;
        tmpNode2.num = num2;
        arrayMethy1[5].push_back(tmpNode2);
    }

    if (arrayMethy1[0].empty())
    {
        coutE << "arrayMethy1[0] is empty; no merged boundary regions." << endl;
        cout << endl << endl << "the total number of intervals: 0" << endl;
        cout << 0 << endl;
        return 0;
    }

    if (arrayMethy1[5].empty())
    {
        coutE << "arrayMethy1[5] is empty; no valid DMP records." << endl;
        cout << endl << endl << "the total number of intervals: 0" << endl;
        cout << 0 << endl;
        return 0;
    }

    cir2 = 0;

    vector<positionNoNode>::iterator it40;
    it40 = arrayMethy1[0].begin();

    while (1)
    {
        numS1 = 0;
        numS2 = 0;

        startP = it40->DMR_S;
        endP = it40->DMR_E;

        it0 = arrayMethy1[5].begin();
        num1 = it0->pos;
        pValue = it0->pV;
        pValue = -log(pValue);
        num2 = it0->num;

        flag = 0;

        posV.clear();
        logPV.clear();
        meUnV.clear();

        while (1)
        {
            if ((num1 >= startP) && (num1 <= endP))
            {
                flag = 1;

                if (num2) numS1++;
                else if (!num2) numS2++;

                posV.push_back(num1);
                logPV.push_back(pValue);
                meUnV.push_back(num2);
            }

            it0++;

            if (it0 == arrayMethy1[5].end()) break;

            num1 = it0->pos;
            pValue = it0->pV;
            pValue = -log(pValue);
            num2 = it0->num;

            if (num1 > endP) break;
        }

        if (flag)
        {
            tmpNode2.pos = startP;
            tmpNode2.end = endP;
            tmpNode2.num = numS1;
            tmpNode2.num2 = numS2;
            tmpNode2.numCom = numS1 + numS2;
            tmpNode2.posV = posV;
            tmpNode2.logPV = logPV;
            tmpNode2.meUnV = meUnV;

            num8 = endP - startP + 1;
            numS3 = 0;
            numS3 = numS1 + numS2;

            ratioTmp = (long double)numS3 / (long double)num8;
            tmpNode2.ratio = ratioTmp;

            if ((numS3 >= M1) && (num8 >= (1 * sWinN)))
                arrayMethy1[6].push_back(tmpNode2);
        }

        cir2++;
        it40++;

        if (it40 == arrayMethy1[0].end()) break;
    }

    cout << endl << endl << "the total number of intervals: " << cir2 << endl;

    sort(arrayMethy1[6].begin(), arrayMethy1[6].end(), cmpV);

    cir = 0;

    vector<positionNoNode>::iterator it1;

    for (it1 = arrayMethy1[6].begin(); it1 != arrayMethy1[6].end(); it1++)
    {
        numTmp1 = it1->pos;
        numTmp2 = it1->end;
        numTmp3 = it1->num;
        numTmp4 = it1->num2;
        numTmp5 = it1->numCom;

        cout3 << setiosflags(ios::left)
              << setw(20) << numTmp1
              << setw(20) << numTmp2
              << setw(20) << numTmp3
              << setw(20) << numTmp4
              << setw(20) << numTmp5
              << endl;

        cir++;
    }

    cout << cir << endl;

    return 0;
}
