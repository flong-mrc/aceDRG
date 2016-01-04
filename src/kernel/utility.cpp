/*
 * File:   Utility.h
 * Author: flong
 *
 * Created on Oct 5, 2011, 7:27 PM
 */

#include "utility.h"

namespace LIBMOL
{
    /* This part contains utility function for 
     * string processing
     */
    
    bool isInt(std::string tS)
    {
        for (int i =0; i < (int)tS.size(); i++)
        {
            if( !std::isdigit(tS[i]))
            {
               return false;
            }
        }
        return true;
    }
    
    int StrToInt(std::string tS)
    {        
        std::stringstream convert(tS);
        int   tNum;
        convert >> tNum;
        return tNum;
    }
    
    std::string IntToStr(int tN)
    {
        std::stringstream convert;
        convert << tN;
        return convert.str();
    }
    
    REAL StrToReal(std::string tS)
    {
        
        std:: stringstream convert(tS);
        REAL tNum;
        convert >> tNum;
        return tNum;
    }

    std::string RealToStr(REAL tN)
    {
        std::stringstream convert;
        convert << tN;
        return convert.str();
    }
    
    REAL StrFractToReal(std::string tS)
    {
        REAL tR=0.0;
        std::vector<std::string> tBuf;
        if (tS.find("/") !=std::string::npos)
        {
            StrTokenize(tS, tBuf, '/');
            if ((int)tBuf.size()==2)
            {
                REAL t1=StrToReal(tBuf[0]);
                REAL t2=StrToReal(tBuf[1]);
                if (std::fabs(t2) > 0.99)
                {
                   tR=t1/t2; 
                }
                else
                {
                    std::cout << "can not transfer "  << tS 
                              << " to a real number " << std::endl;
                }
            }
            else
            {
                std::cout << "can not transfer "  << tS 
                          << " to a real number " << std::endl;
            }
        } 
        else
        {
            StrToReal(tS);
        }
        
        return tR;
    }
    
    void getDigitSec(std::string & tStrIn, std::string & tStrOut)
    {
        for (unsigned i=0; i < tStrIn.size(); i++)
        {
            if (std::isdigit(tStrIn[i]))
            {
                tStrOut.push_back(tStrIn[i]);
            }
        }
    }
    std::string TrimSpaces( std::string  tStr)
    {
        // Trim Both leading and trailing spaces
       std::size_t startPos = tStr.find_first_not_of(" \t"); 
       std::size_t endPos   = tStr.find_last_not_of(" \t");
       
       // if all spaces or empty return an empty string
       if(( std::string::npos == startPos) 
           || (std::string::npos == endPos))
       {
           tStr = "";
       }
       else
       {
           tStr = tStr.substr( startPos, endPos-startPos+1);   
       }
       
       return tStr;
    }
    
    void cleanChar(std::string & tStr,
                   char          tChar)
    {
        size_t found =tStr.find_first_of(tChar);
        while (found!=std::string::npos)
        {
            tStr.erase(found,1);
            found=tStr.find_first_of(tChar, found+1);
        }
    }
    
    void  StrUpper(std::string & tStr)
    {
        std::transform(tStr.begin(), tStr.end(), 
                       tStr.begin(), ::toupper);
    }
    
    void  StrLower(std::string & tStr)
    {
        std::transform(tStr.begin(), tStr.end(), 
                       tStr.begin(), ::tolower);
    }
     
    void  StrTokenize(std::string     tStr, 
            std::vector<std::string>  & tokens)
    {
        std::string buf;
        std::stringstream ss(tStr);
        
        while (ss >> buf)
        {
            tokens.push_back(buf);
        }
         
    }
    
    std::vector<std::string> StrTokenize(const std::string &tStr,  char delim)
    {
        std::stringstream sS(tStr);
        std::string tToken;
        std::vector<std::string> tV;
        
        while (std::getline(sS, tToken, delim))
        {
            tV.push_back(tToken);
        }
         
        return tV;
    }
    
    void  StrTokenize(const std::string &tStr, std::vector<std::string>  & tV, char delim)
    {
        std::stringstream sS(tStr);
        std::string tToken;
        
        while (std::getline(sS, tToken, delim))
        {
            tV.push_back(tToken);
        }
    }
    
    bool compareNoCase (std::string first, std::string second)
    {
        unsigned int i=0;
        while ( (i<first.length()) && (i<second.length()) )
        {
            if (std::toupper(first[i])<std::toupper(second[i])) return true;
            else if (std::toupper(first[i])>std::toupper(second[i])) return false;
            ++i;
        }
        return (first.length() > second.length() ?  true : false);
    }
    
    bool compareNoCase2 (std::string first, std::string second)
    {
        if (first.length() > second.length())
        {
             return  true;
        }
        else if (first.length() < second.length())
        {
            return false;
        }
        
        unsigned int i=0;
        while ( (i<first.length()) && (i<second.length()) )
        {
            if (std::toupper(first[i])<std::toupper(second[i])) return true;
            else if (std::toupper(first[i])>std::toupper(second[i])) return false;
            ++i;
        }
        return true;
    }
    
    bool compareNoCase3 (std::string first, std::string second)
    {
        if (first.length() < second.length())
        {
             return  true;
        }
        else if (first.length() > second.length())
        {
            return false;
        }
        
        unsigned int i=0;
        while ( (i<first.length()) && (i<second.length()) )
        {
            if (std::toupper(first[i])<std::toupper(second[i])) return true;
            else if (std::toupper(first[i])>std::toupper(second[i])) return false;
            ++i;
        }
        return true;
    }
    
    bool sortStrByLen(std::string first, std::string second)
    {
        return (int)first.size() > (int)second.size();
    }
    
    
    bool desSortMapKey(const sortMap& a ,const sortMap & b)
    {
        return a.key.length() > b.key.length();
    }
    
    bool desSortMapKey1(const sortMap& a ,const sortMap & b)
    {
        return compareNoCase3(a.key, b.key);
    }
    
    bool desSortMapKey2(const sortMap2 & a ,const sortMap2 & b)
    {
        if (a.key.length() > b.key.length())
        {
            return true;
        }
        else if (a.key.length() == b.key.length())
        {
            if(a.nNB > b.nNB)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    
    bool desSortMapKey3(const sortMap3 & a ,const sortMap3 & b)
    {
        return compareNoCase2(a.key, b.key);
    }
    
    bool sortMapkey4(const sortMap4 & a ,const sortMap4 & b)
    {
        return compareNoCase2(a.key, b.key);
    }
    
    bool compareNoCaseClass(const sortLine & first, const sortLine & second)
    {
        unsigned int i=0;
        while ( (i<first.key.length()) && (i<second.key.length()) )
        {
            if (std::toupper(first.key[i])<std::toupper(second.key[i])) return true;
            else if (std::toupper(first.key[i])>std::toupper(second.key[i])) return false;
            ++i;
        }
        return (first.key.length() > second.key.length() ?  true : false);
    }
    
    void cleanSymbol(std::string & tStr, std::string symb)
    {
        std::size_t found;
        found=tStr.find_first_of(symb);
        while (found!=std::string::npos)
        {
            tStr.erase(found,1);
            found=tStr.find_first_of(symb, found+1);
        }
    }

    int getKeyWordPos(std::string  tKey, std::vector<std::string> & tKeyVect)
    {
        int aP=-1;
        for (int i=0; i < (int)tKeyVect.size(); i++)
        {
            if (tKeyVect[i].compare(tKey)==0)
            {
                aP=i;
                break;
            }
        }
        
        return aP;
    }
    
    /*  end of string processing parts */
    
    // zero-> the file does not exist
    int isFileExist(FileName aF)
    {
        struct stat buffer;
        
        if ( stat( aF, &buffer ) ) return 0 ;
        
        return 1 ;
    }
    
    // trigonometrical functions
    extern REAL degreeToRadians(REAL tDeg)
    {
        return tDeg*PI180;
    }
    extern REAL RadiansToDegree(REAL tRad)
    {
        return tRad*PID180;
    }
    
    /* The following are some math functions 
     *  acting on REAL type (change them to template later on)
     */

    extern bool inVect(std::vector<REAL> & tVect,
                       REAL tVal, REAL tErr)
    {
        for (std::vector<REAL>::iterator iV=tVect.begin();
                iV !=tVect.end(); iV++)
        {
            if (*iV-tVal <tErr)
            {
                return true;
            }
        }
        return false;
    }
    
    extern bool inVectABS(std::vector<REAL> & tVect,
                       REAL tVal, REAL tErr)
    {
        for (std::vector<REAL>::iterator iV=tVect.begin();
                iV !=tVect.end(); iV++)
        {
            if (fabs(*iV-tVal) <tErr)
            {
                return true;
            }
        }
        
        return false;
    }
    
    extern bool outVectAbsDiff(std::vector<REAL> & tVect,
                              REAL tVal, REAL tErr)
    {
        REAL x1=fabs(tVal);
        for (std::vector<REAL>::iterator iV=tVect.begin();
                iV !=tVect.end(); iV++)
        {
            REAL x2 = fabs(*iV);
            // std::cout << " diff is " << 2.0*fabs(x1-x2)/(x1+x2) << std::endl;
            if (x1 >= 0.0000001 && x2 >=0.0000001)
            {
                if (2.0*fabs(x1-x2)/(x1+x2) > tErr)
                {
                    return true;
                }
            }
            else if (fabs(x1-x2) > tErr)
            {
                return true;
            }
        }
        
        return false;
    }
    
    extern void getFracReal(REAL tR, REAL & tFrac, REAL tTar)
    {
        // select a fraction of real number tR which is small than tTar,
        // give fraction number;
        int aN = 10;
        
        tFrac = 1.0;
        while (tFrac*tR > tTar && aN > 1)
        {
            tFrac -=0.1;
            aN--;
        }
    }
    
    extern REAL lengthV(std::vector<REAL> & tV1)
    {
        REAL tBuf = 0.0;
        for (int i = 0; i < (int)tV1.size(); i++)
        {
            tBuf +=(pow((tV1[i]), 2));
        }
        
        return sqrt(tBuf);         
    }
    
    extern REAL distanceV(std::vector<REAL> & tV1, std::vector<REAL> & tV2)
    {
        REAL tBuf = 0.0;
        
        if (tV1.size() == tV2.size())
        {
            for (int i = 0; i < (int)tV1.size(); i++)
            {
                tBuf +=(pow((tV1[i]-tV2[i]), 2));
            }
        }
        
        return  std::sqrt(tBuf);
    }
    
    extern void normalizeV(std::vector<REAL> & tV)
    {
        REAL tMax=0.0;
        for (int i=0; i < (int)tV.size(); i++)
        {
            REAL tF=fabs(tV[i]);
            if (tF > tMax)
            {
                tMax=tF;
            }
        }
        
        if (tMax > 1)
        {
            for(int i=0; i < (int)tV.size(); i++)
            {
                tV[i]=tV[i]/tMax;
            }
        }
    }
    
    extern REAL dotP(std::vector<REAL> & tV1, std::vector<REAL> & tV2)
    {
        REAL tBuf = 0.0;
        
        if (tV1.size() == tV2.size())
        {
            for (int i = 0; i < (int)tV1.size(); i++)
            {
                tBuf += tV1[i]*tV2[i];
            }
        }
        
        return tBuf;
        
    }
    
    void  crossP2V(std::vector<REAL> & tV1, 
                   std::vector<REAL> & tV2,
                   std::vector<REAL> & tV3)
    {   
       
        if(tV1.size()==3 && tV2.size()==3)
        {
            tV3.push_back(tV1[1]*tV2[2]-tV1[2]*tV2[1]);
            tV3.push_back(tV1[2]*tV2[0]-tV1[0]*tV2[2]);
            tV3.push_back(tV1[0]*tV2[1]-tV1[1]*tV2[0]);
        }
    }
    
    REAL getAngle2V(std::vector<REAL> & tV1, std::vector<REAL> & tV2)
    {
        REAL tA = 0.0;
        
        REAL tMOD1 =0.0, tMOD2=0.0, tMOD=0.0;

        tMOD1 = lengthV(tV1);
        if(tMOD1 <=1.0e-16)
        {
            return tA;
        }
        
        
        tMOD2 = lengthV(tV2);
        if(tMOD2 <=1.0e-16)
        {
            return tA;
        }       
        
        tMOD = dotP(tV1, tV2);
        
        return acos(tMOD/(tMOD1*tMOD2));
    }
    
    extern REAL getTorsion3V(std::vector<REAL> & tV1, std::vector<REAL> & tV2,
                      std::vector<REAL> & tV3)
    {
        std::vector<REAL> tV2V3;
        crossP2V(tV2, tV3, tV2V3);
        
        //std::cout << "tV2V3[0]: " << tV2V3[0] << std::endl;
        //std::cout << "tV2V3[1]: " << tV2V3[1] << std::endl;
        //std::cout << "tV2V3[2]: " << tV2V3[2] << std::endl;
        REAL tA = -dotP(tV2,tV2)*dotP(tV1,tV3)+dotP(tV1,tV2)*dotP(tV2,tV3);
        //std::cout << "tA: " << tA << std::endl;
        
        REAL tB = std::sqrt(dotP(tV2,tV2))*dotP(tV1,tV2V3);
        //std::cout << "tB: " << tB << std::endl;
        REAL tAB2 = std::sqrt(pow(tA,2.0)+pow(tB,2.0));
        //std::cout << "tAB2: " << tAB2 << std::endl;
        
        if (fabs(tB) <=1.0e-16 ||  tAB2 <=1.0e-16 )
        {
            return 0.0;
        }
       
        return Sign(1.0,tB)*acos(tA/tAB2)*PID180;
    }
    
    extern REAL getDet(std::vector<REAL> & v1,
                       std::vector<REAL> & v2,
                       std::vector<REAL> & v3)
    {
        REAL t1, t2, t3;
        t1 =0.0;
        t2 =0.0;
        t3 =0.0;
        
        if(v1.size() == 3 && v2.size() ==3 
               && v3.size() == 3)
        {
            t1 = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1]);
            t2 = v1[1]*(v2[2]*v3[0]-v2[0]*v3[2]);
            t2 = v1[2]*(v2[0]*v3[1]-v2[1]*v3[0]); 
            return  (t1 + t2 + t3);
        }
        else
        {
            return 0.0;
        }
    }
    
    extern void matCopyInt(std::vector<std::vector<int> > & tM1,
                           std::vector<std::vector<int> > & tM2)
    {
        if (not tM2.empty())
        {
            tM2.clear();
        }
        
        for (unsigned i=0; i < tM1.size(); i++ )
        {
            std::vector<int> aRow;
            for (unsigned j=0; j < tM1[i].size(); j++)
            {
                aRow.push_back(tM1[i][j]);
            }
            tM2.push_back(aRow);
        }
    }
    extern void matMultVec(std::vector<std::vector<REAL> >  & tMat,
                std::vector<REAL> & tInitV, std::vector<REAL> & tFinV)
    {
        tFinV.clear();
        for (std::vector<std::vector<REAL> >::iterator iRow=tMat.begin();
                iRow != tMat.end(); iRow++)
        {
            if (iRow->size() == tInitV.size())
            {
                REAL tSum =0.0;
                for (unsigned i=0; i<iRow->size(); i++)
                {
                    tSum+=((*iRow)[i]*tInitV[i]);
                }
                tFinV.push_back(tSum);       
            }
            else
            {
                std::cout << "dimension error in a matrix times a vector" << std::endl;
                exit(1);
            }
        }
    }
    
    extern void matMultMatInt(std::vector<std::vector<int> >    & tMat1,
                              std::vector<std::vector<int> >    & tMat2,
                              std::vector<std::vector<int> >    & tOutMat)
    {
        unsigned nRow1 = tMat1.size();
        unsigned nRow2 = tMat2.size();
        
        if (nRow1 > 0 && nRow2 > 0)
        {
            unsigned nCol1 = tMat1[0].size();
            unsigned nCol2 = tMat2[0].size();
            
            if (nCol1 == nRow2 && nCol2 > 0)
            {
                if (!tOutMat.empty())
                {
                    tOutMat.clear();
                }
                for (unsigned i=0; i < nRow1; i++)
                {
                    std::vector<int> aRow;
                    for (unsigned j=0; j < nCol2; j++)
                    {
                        aRow.push_back(0.0);
                    }
                    tOutMat.push_back(aRow);
                }
                
                for (unsigned i=0; i < nRow1; i++)
                {
                    for (unsigned j=0; j < nCol2; j++)
                    {
                        // std::cout << i << " " << j << std::endl;
                        for (unsigned k=0; k < nCol1; k++)
                        {
                            // std::cout << "k  " << tMat1[i][k] << "  " 
                            //           <<  tMat2[k][j] << std::endl;
                            
                            tOutMat[i][j]+=(tMat1[i][k]*tMat2[k][j]);
                        }
                    }
                }
            }
        }
    }
    
    extern void printMatrix(std::vector<std::vector<int> >    & tMat)
    {
        unsigned nR=tMat.size();
        if (nR > 0)
        {
            unsigned nC=tMat[0].size();
            std::cout << std::setw(6) << "C-R";
            for (unsigned i=0; i < nR; i++)
            {
                 std::cout << std::setw(6) << i;
            }
            std::cout << std::endl;
            
            for (unsigned i=0; i < nR; i++)
            {
                std::cout << std::setw(6) << i;
                for (unsigned j=0; j < nC; j++)
                {
                    std::cout << std::setw(6) << tMat[i][j];
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // 
    //   The following are functions imported directly from
    //   a old global minimization package
    //
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    // The following are methods on the operation of
    // vector, matrix, bond angle and torsion angle
 
    // 1.  Calculate the length of a vector

    extern REAL length_v(REAL *v, int tDim)
    {
        int i;
        REAL len;

        len = 0.0;
 
         for (i =0; i < tDim; i++)
         {
             len += (v[i]*v[i]);
         }

         len = sqrt(len);

         return len;
    }
    
    extern REAL DotP(REAL *v1,  REAL *v2)
    {

        int dim =3;
        REAL product=0.0;
        
        for (int i =0; i < dim; i++) 
        {
            product += v1[i]*v2[i];
        }
        
        return product;
    }
    
    extern REAL DotP(int n_size, REAL *v1,  REAL *v2)
    {
        REAL product=0;
  
        for (int i =0; i < n_size; i++) 
        {
            product += v1[i]*v2[i];
        }
        
        return product;
    }
    
    extern REAL DotP(std::vector<REAL> & v1, REAL *v2)
    {
        REAL p=0.0;
        
        for (int i=0; i < (int)v1.size(); i++)
        {
            p += (v1[i]*v2[i]);
        }
        
        return p;
    }
    
    extern void CrossP(REAL *v1, REAL *v2, REAL *v3)
    {
        // consider three dimension case only at this stage
        int      i;
   
        int  dim = 3;   
        REAL v_tm[dim];

        for (i =0; i < dim; i++)
        {
            v3[i] =0.0;
        }

        v_tm[0] = v1[1]*v2[2]-v1[2]*v2[1];
        v_tm[1] = v1[2]*v2[0]-v1[0]*v2[2];
        v_tm[2] = v1[0]*v2[1]-v1[1]*v2[0];

        for ( i =0; i < dim; i++)
        {
            v3[i] = v_tm[i];
        }
    }
    
    extern REAL GetAngle(REAL *v1, REAL *v2)
    {
        REAL MOD1, MOD2, MOD;
        REAL  Angle;
        /*
        for (int i =0; i < 3; i++)
        {
            std::cout << "v1["<<i<<"]= " << v1[i] << std::endl;
            std::cout << "v2["<<i<<"]= " << v1[i] << std::endl;
        }
         */

        MOD1 = sqrt(DotP(v1,v1));
        MOD2 = sqrt(DotP(v2,v2));
        MOD  = DotP(v1,v2);

        if (MOD1 <=1.0e-8 || MOD2 <= 1.0e-8 )
        {
            std::cout << " one of vectors constructing the angle is zero. " << std::endl;
            exit(1);
        }
        else 
        {  
            REAL MOLT= MOD/(MOD1*MOD2);
            if(MOLT>0.9999)
            {
                Angle = PI;
            }
            else
            {
                Angle = acos(MOD/(MOD1*MOD2));
            }
            
            //std::cout << "MOD1 " << MOD1 << std::endl;
            //std::cout << "MOD2 " << MOD2 << std::endl;
            //std::cout << "MOD  " << MOD  << std::endl;
            return Angle;
            
        }
        
    }

    
    // matrix operation 
    
    // Matrix copy
    extern void MatCopy(int          n_size,
                        REAL      ** A_old,
                        REAL      ** A_new)
    {
        for (int i =0; i < n_size; i++)
        {
            for (int j =0; j < n_size; j++)
            {
                A_new[i][j] = A_old[i][j];
            }
        }
    }

    // Vector operation. convert the vector V_old to the vector V_new
    //   [V_new] = [AM][V_old]
    extern void MatMultip(int         n_dim,
                          REAL     ** AM,
                          REAL     *  V_old,
                          REAL     *  V_new)
    {
        for (int i =0; i < n_dim; i++)
        {
            V_new[i] =0.0;
            for (int j =0; j < n_dim; j++)
            {
                V_new [i] += AM[i][j]*V_old[j];
                //  std::cout << " AM[" << i<<"][" << j <<" ] " << AM[i][j]
                //      << " V_old["<< j << "] " << V_old[j] 
                //      << " V_new+=  " << V_new [i] ;
            }
            //   std::cout << endl;
            //   std::cout << "V_new [" << i << "] =" << V_new [i] << endl;
        }
    }

    
    extern REAL CalcDet(REAL ** Matr)
    {
        
        REAL t1=0.0, t2=0.0, t3=0.0;
        REAL T =0.0;

        t1 = Matr[0][0]*(Matr[1][1]*Matr[2][2]-Matr[1][2]*Matr[2][1]);
        t2 = Matr[0][1]*(Matr[1][2]*Matr[2][0]-Matr[1][0]*Matr[2][2]);
        t3 = Matr[0][2]*(Matr[1][0]*Matr[2][1]-Matr[1][1]*Matr[2][0]);

        T = t1+t2+t3;
 
        return T;
    }
    
    extern REAL CalcDet(REAL  * v1,
                        REAL  * v2,
                        REAL  * v3)
    {
        REAL t1=0.0, t2=0.0, t3=0.0;
        REAL T =0.0;
        
        t1 = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1]);
        t2 = v1[1]*(v2[2]*v3[0]-v2[0]*v3[2]);
        t2 = v1[2]*(v2[0]*v3[1]-v2[1]*v3[0]); 

        T = t1+t2+t3;

        return T;
        
    }
    
    // For sparse matrix format 
    extern void Prod_Mat_Vec( int  n_dim,   int  n_block,
                              REAL    * secDerivCartSparse_row,
                              REAL    * secDerivCartSparse_col,
                              REAL  *** secDerivCartSparse,
                              REAL    * vect_old, 
                              REAL    * vect_new)                          
    {
        int i,j;
        int i_atom, i_block, i_block_fin; 
        int i_row, j_col;

        for (i_atom = 0; i_atom < n_block; i_atom++)  // block rows
        {

            i_block = secDerivCartSparse_row[i_atom];  // starting and finishing
            i_block_fin = secDerivCartSparse_row[i_atom+1]; // blocks in a row

            for( i = 0; i < n_dim; i++)
            {
                i_row = i_atom*n_dim+i;
                vect_new [i_row] = 0.0;
            }
            
            do
            {
                for ( i = 0; i < n_dim; i++)
                {
                    i_row = i_atom*n_dim+i;
                    for (j = 0; j < n_dim; j++)
                    {
                        j_col = secDerivCartSparse_col[i_block] + j;
                        vect_new [i_row]+=(secDerivCartSparse [i_block][i][j] 
                                           *vect_old[j_col]);
                    }
                }
                i_block++;   
            }while(i_block < i_block_fin);
        
        }
    }
    
    // Normalization of the matrix, A, and its associate
    // vector b, where A*X=b.
    // Overlaoded version 1.
    extern void MatNormolization(int      n_size,
                                 REAL  *  elemA_dia,
                                 REAL  ** A,
                                 REAL  *  b)
    {
        // Normalize A and b 
        
        
        REAL  t_dia = 0;
       
        //std::ofstream Matrix_File("Matrix_B.data");
        for (int i =0; i < n_size; i++)
        {
            b[i] = b[i]/elemA_dia[i];
            for (int j =0; j < n_size; j++)
            {
                t_dia = elemA_dia[i]*elemA_dia[j];
                //Matrix_File << "A["<<i+1 << "," << j+1 << "]= " 
                //     << A[i][j] << std::endl;
                //Matrix_File << "t_dia=" << t_dia << std::endl;
                A[i][j] = A[i][j]/t_dia;
                //Matrix_File << "A["<<i+1 << "," << j+1 << "]= " 
                //     << A[i][j] << std::endl;
            }
        }
        //Matrix_File.close();
   
    }
    
    // Normalization of matrix A only.
    // Overload version 2
    extern void MatNormolization(int      n_size,
                                 REAL  ** A)
    {
        int i,j;

        REAL ** T = new REAL * [n_size];
 
        for(i = 0; i < n_size; i++)
        {
            T[i] = new REAL [n_size];
            for (j = 0; j < n_size; j++)
            {
                T[i][j] = 0.0;
            }
        }

        for (i = 0; i < n_size; i++)
        {
            if(fabs(A[i][i]) < 1.0e-16)
            {
                A[i][i] = 1.0;
            }
        }

        REAL t_dia = 0;
        for (i =0; i < n_size; i++)
        {
            t_dia = A[i][i]; 
            for (j =0; j < n_size; j++)
            {
                // t_dia = sqrt(A[i][i]*A[j][j]);
                T[i][j] = A[i][j]/t_dia;         
            }
        }

        MatCopy(n_size, T, A);
  
        // ofstream Matrix_File("Matrix_A.data");
        // for (i = 0; i < n_size; i++)
        //  {
        //    for (j =0; j < n_size; j+=3)
        //    {
        //         Matrix_File << "T["<<i+1 << "," << j+1 << "]= " 
        //             << T[i][j] << "   " 
        //             << "T["<<i+1 << "," << j+2 << "]= " 
        //              << T[i][j+1] << "   " 
        //              << "T[" << i+1 << "," << j+3 << "]= " 
        //             << T[i][j+2] << endl;
        //    }                   
        //  }

        // Matrix_File.close();

        for (i =0; i < n_size; i++)
        {
            delete [] T[i];
            T[i] = 0;
        }
        delete [] T;
        T = 0;

    }  
    
    // Normalization of the matrix, A, and its associate
    // vector b, where A*X=b. Sparse matrix format BCRS 
    // is used.
    // Overlaoded version 3. 
    
    extern void MatNormolization(int        n_size,
                                 int        n_dim,
                                 int        n_block,
                                 REAL    *  elemA_dia,
                                 REAL    *  b,
                                 REAL    *  secDerivCartSparse_row,
                                 REAL    *  secDerivCartSparse_col,
                                 REAL  ***  secDerivCartSparse)                                
    {
        int i,j,k;
        int idx_atom, idx_block;

        REAL  *** T_spa;
        T_spa = new REAL ** [n_size];
        for (i=0; i < n_size; i++)
        {
            T_spa[i] = new REAL * [n_dim];
            for(j =0; j < n_dim; j++)
            {
                T_spa[i][j] = new REAL [n_dim];
                for (k =0; k < n_dim; k++)
                {
                    T_spa[i][j][k] = 0.0;
                }
            }
        }

        // Normalize the matrix
  
        int i_row_s, i_row_f, j_col, j_col_s;
        REAL product_ij;

        for (idx_atom = 0; idx_atom < n_block; idx_atom++)
        {
            i_row_s = secDerivCartSparse_row[idx_atom];
            i_row_f = secDerivCartSparse_row[idx_atom+1];
  
            //      cout << "for block row " << idx_atom+1 << endl
            //     << "The row start at " << i_row_s << endl
            //     << " and end at " << i_row_f-1 << " in A " << endl;
      
            for (idx_block =i_row_s; idx_block < i_row_f; idx_block++)
            {  
                j_col   = secDerivCartSparse_col[idx_block]/n_dim;
                j_col_s = secDerivCartSparse_row[j_col];
          
                //  cout << "j_col = " << j_col << endl;
                //   cout << "j_col_s = " << j_col_s << endl;
          
                REAL aa, bb;

                for (i = 0; i < n_dim; i++)
                {
                    for(j = 0; j < n_dim; j++)
                    {
                        aa = secDerivCartSparse[i_row_s][i][i];
                        if(fabs(aa) < 1.0e-8)
                        {
                            aa = 1.0;
                        }

                        bb = secDerivCartSparse[j_col_s][j][j];
                        if(fabs(bb) < 1.0e-8)
                        {
                            bb = 1.0;
                        }
                        product_ij = aa*bb;

                        //  product_ij = secDerivCartSparse[i_row_s][i][i]
                        //             *secDerivCartSparse[j_col_s][j][j];
                        product_ij = sqrt(product_ij);
                        if(fabs(product_ij) < 1.0e-16)
                        {
                            std::cout << "zero diagonal element at block "
                                      << idx_block << " position " << i
                                      << " or " << j << std::endl;
                            exit(1);
                        }
                        else
                        {
                            T_spa[idx_block][i][j] = 
                            (secDerivCartSparse[idx_block][i][j]/product_ij);
                        }

                    }
                }
            }
        }

        // Normalize the vector 

        int n_dim_all = n_dim*n_block;

        for(i = 0; i < n_dim_all; i++)
        {
            b[i] /= elemA_dia[i];
            // std::cout << " elemA_dia[" << i << "]= " << elemA_dia[i] << "\t";
            // std::cout << "b[" << i << "] = " << b[i] << endl;
        }
  
        //  copy back

        for (idx_atom = 0; idx_atom < n_block; idx_atom++)
        {
            i_row_s = secDerivCartSparse_row[idx_atom];
            i_row_f = secDerivCartSparse_row[idx_atom+1]; 
      
            for (idx_block =i_row_s; idx_block < i_row_f; idx_block++)
            { 
                for(i =0; i < n_dim; i++)
                {
                    for(j =0; j < n_dim; j++)
                    {
                        secDerivCartSparse[idx_block][i][j] 
                                   =  T_spa[idx_block][i][j];
                    }
                }
            }

        }
 
       // Release memory

      for (i=0; i < n_size; i++)
      {
          for(j =0; j < n_dim; j++)
          {
              delete [] T_spa[i][j];
              T_spa[i][j] = 0;
          }
          delete [] T_spa[i];
          T_spa[i] = 0;
      }
      delete [] T_spa;
      T_spa = 0;
  
    }

    extern void UnNormalization (int n_size, REAL *  A_mat_dia,
                                 REAL *  X_vect)
    {
        for (int i = 0; i < n_size; i++)
        {
            X_vect[i] /= A_mat_dia[i];
        }
    }
    
    // Full matrix version 
    extern REAL ConjugateP(int n_size,  REAL * Vect_1,
                           REAL ** A,   REAL *  Vect_2)
    {
        REAL * T1 = new REAL [n_size];
        REAL   T2 = 0.0;

        for (int i =0; i < n_size; i++)
        {
            T1[i] = 0.0;
        }

        MatMultip(n_size, A, Vect_1, T1);

        T2 = DotP(n_size, T1,Vect_2);

        delete [] T1;
        T1 = 0;

        return T2;
    }
    
    // Sparse matrix version
    extern REAL ConjugateP(int n_size,
                           int n_dim,
                           int n_block,
                           REAL    * secDerivCartSparse_row,
                           REAL    * secDerivCartSparse_col,
                           REAL  *** secDerivCartSparse,
                           REAL *    Vect_1,
                           REAL *    Vect_2)
    {
        
          REAL * T1= new REAL [n_size];
          REAL   T2;

          for (int i =0; i < n_size; i++)
          {
              T1[i] = 0.0;
          }

          Prod_Mat_Vec(n_dim, n_block, secDerivCartSparse_row,
                       secDerivCartSparse_col, secDerivCartSparse,
                       Vect_1, T1);
 
          // for(int i = 0; i < n_size; i++)
          //  {
          //    std::cout << "Vect[" << i<<"] " << Vect_1[i] << "\t"
          //         << "T[" << i<<"] " << T1[i] << std::endl;
          //  }
 

          T2 = 0.0;

          T2 = DotP(n_size, T1,Vect_2);
 

          delete [] T1;
          T1 = 0;

          return T2;   
    }
    
    
    // Solve an eigenvalue equation using TNT lib 
    extern void EigenSolve(int         size_matrix,
                           REAL    **  in_matrix,
                           REAL    *   e_value,
                           REAL    **  e_vect)
    {
        
        int i,j;
        
        Array2D<REAL> A_matrix(size_matrix, size_matrix);
        Array2D<REAL> EigenVec_matrix(size_matrix, size_matrix);
        Array1D<REAL> EigenVal(size_matrix);
        
        //  realtype product_tm[size_matrix][size_matrix];
        
        for(i =0; i < size_matrix; i++)
        {
            for(j =0; j < size_matrix;  j++)
            {
                A_matrix[i][j] = in_matrix[i][j];
            }
        }
        
        // Instanciate the EigenValue class and solve the eigenvalue problem

        JAMA::Eigenvalue<REAL> Eigen_prob(A_matrix);
        
        // Obtain the eigenvalues and eigenvectors

        Eigen_prob.getRealEigenvalues(EigenVal);

        Eigen_prob.getV(EigenVec_matrix);

        for(i = 0; i < size_matrix; i++)
        {
            e_value[i] = EigenVal[i];
            for (j =0; j < size_matrix; j++)
            {
                e_vect[i][j] = EigenVec_matrix[j][i];
            }
        }
    }
    
    
    // Construct a transformation matrix given the components of
    // the new coordinate axis in the old (fixed) coordinate
 
    extern void InitAMat2(REAL  * X_n,
                          REAL  * Y_n,
                          REAL  * Z_n,
                          REAL ** a_mat,
                          int   t_dim)
    {
        int i, j;
        
        for (i =0; i < t_dim; i++)
        {
            for (j =0; j < t_dim; j++)
            {
                a_mat[i][j] = 0.0;
            }
        }
        
        for ( i = 0; i < t_dim; i++)
        {
            a_mat[i][0] = X_n [i];
            a_mat[i][1] = Y_n [i];
            a_mat[i][2] = Z_n [i];
        }

    }
    
    extern void  NBD_MCOPY(REAL ** A, REAL ** B, int t_dim)
    {
        
        int i, j;
        for (i =0; i < t_dim; i++)
        {
            for (j =0; j < t_dim; j++)
            {
                B[i][j] = 0.0;
            }
        }
  
        for (i =0; i < t_dim; i++)
        {
            for (j =0; j < t_dim; j++)
            {
                B[i][j] = A[i][j];
            }
        }
    }

    extern void  NBD_MATMLT(REAL ** a_mat, 
                            REAL ** Y,
                            REAL ** Z, 
                            int     t_dim)
    {
        int i, j, k; 
        
        for (i =0; i < t_dim; i++)
        {
            for (j =0; j < t_dim; j++)
            {
                Z[i][j] = 0.0;
                for(k =0; k < t_dim; k++)
                {
                    Z[i][j] = Z[i][j] + a_mat[i][k]*Y[k][j];
                }
            }
        }

    }
    
    // Geometry related 
    extern void Plane_Find(int        n_points, 
                           REAL   **  X_in,
                           REAL   *   co)
    {
        
        int i1, i2, i3;

        // Setup a 3 x 3 matrix
        int dim = 3;
        REAL  ** A = new REAL * [dim];
  
        for (i1 =0; i1 < dim; i1++)
        {
            A[i1] = new REAL [dim];
            for(i2 = 0; i2 < dim; i2++)
            {
                A[i1][i2] = 0.0;
                for(i3 = 0; i3 < n_points; i3++)
                {
                    A[i1][i2] += (X_in[i3][i1]*X_in[i3][i2]);
                }
            }
        }
   
        // Solve eigenvalue equation of [A][X] = lamda*[X]
        
        REAL *  eigenValue_A = new REAL [dim];
        REAL ** eigenVect_A  = new REAL * [dim];
        for(i1 =0; i1 < dim; i1++)
        {
            eigenVect_A[i1] = new REAL [dim];
        }
        
        EigenSolve(dim, A, eigenValue_A, eigenVect_A);
        
        // Get the coefficients. The eigen vectors corresponding
        // to the smallest eigenvalue
        REAL min_va = eigenValue_A[0];
        int      i_min  = 0;
        for(i1 = 1; i1 < dim; i1++)
        {
            if(eigenValue_A[i1] < min_va)
            {
                i_min = i1;
            }
        }
        
        for (i1 =0; i1 < dim; i1++)
        {
            co[i1] = eigenVect_A[i_min][i1];
        }
        
        for (i1 =0; i1 < dim; i1++)
        {
            delete [] A[i1];
            A[i1] = 0;
        }
        delete [] A;
        A = 0;

        delete []  eigenValue_A;
        eigenValue_A = 0;

        for(i1 =0; i1 < dim; i1++)
        {
            delete [] eigenVect_A[i1];
            eigenVect_A[i1] =0;
        }
        delete []  eigenVect_A;
        eigenVect_A = 0;

    }
    
    // Find coefficients that describe a plane which is fitting,
    // via least square, to a group of points, and the derivatives
    // of these coefficients wrt the atomic coordinates  
    extern void  Plane_and_Deriv_Find(int         n_points,
                                      REAL     ** X_in,
                                      REAL      * co,
                                      REAL     ** d_co)
    {
        int dim = 3;
        int i_dim, i_point;

        
        // 1. Find the coefficients of the plane equation

        int i1, i2, i3;
             
        // Setup a t_dim x t_dim matrix

        REAL  ** A = new REAL * [dim];
        for (i1 =0; i1 < dim; i1++)
        {
            A[i1] = new REAL [dim];
            for(i2 = 0; i2 < dim; i2++)
            { 
                A[i1][i2] = 0.0;
                for(i3 = 0; i3 < n_points; i3++)
                {
                    A[i1][i2] += (X_in[i3][i1]*X_in[i3][i2]);
                }
            }
        }
        
        
        // Solve eigenvalue equation of [A][X] = lamda*[X]

        REAL *  eigenValue_A = new REAL   [dim];
        REAL ** eigenVect_A  = new REAL * [dim];
        for(i1 =0; i1 < dim; i1++)
        {
            eigenVect_A[i1] = new REAL [dim];
        }
        
        EigenSolve(dim, A, eigenValue_A, eigenVect_A);

        // Get the coefficients. The eigen vectors corresponding
        // to the smallest eigenvalue

        REAL min_va = eigenValue_A[0];
        
        // std::cout << "1st EigenValue " <<   eigenValue_A[0] << std::endl;
        
        int      i_min  = 0;
        for(i1 = 1; i1 < dim; i1++)
        {
            if(eigenValue_A[i1] < min_va)
            {
                i_min = i1;
            }
            //std::cout << i1+1 << " eigenValue " << eigenValue_A[i1] << std::endl;
        }
        
        for (i1 =0; i1 < dim; i1++)
        {
            co[i1] = eigenVect_A[i_min][i1];
        }
        
        // 2. Find the relative eigenvalues and its inverses

        REAL * rel_eigenValue_A = new REAL [dim];
        REAL * inv_rel_eigenValue_A = new REAL [dim];
        for (i1 =0; i1 < dim; i1++)
        {
            rel_eigenValue_A[i1] = eigenValue_A[i1]-eigenValue_A[i_min];
            if(rel_eigenValue_A[i1] > 1.0e-6)
            {
                inv_rel_eigenValue_A[i1] = 1.0/rel_eigenValue_A[i1];
            }
            else
            {
                inv_rel_eigenValue_A[i1] = 0.0;
            }   
        }
        
        // 3. Find the inverse of the Matrix A via product of the 
        //    inverse of the relative eigenvalues and the matrix
        //    of eigenvectors.
        
        REAL ** inv_A;
        inv_A = new REAL * [dim];
        for (i1 =0; i1 < dim; i1++)
        {
            inv_A[i1] = new REAL [dim];
            for (i2 =0; i2 < dim; i2++)
            {
                inv_A[i1][i2] =0.0;
                for (i3 =0; i3 < dim; i3++)
                {
                    inv_A[i1][i2] += (eigenVect_A[i3][i1]                
                                     *inv_rel_eigenValue_A[i3]
                                     *eigenVect_A[i3][i2]);
                }
            }
        }   
        
        // 4. setup vectors g(x), the righthand side
        //    of plane derivative equation 

        REAL  d_lm_x, d_lm_y, d_lm_z;
        
        int total_dim = n_points*dim;
        REAL  ** g = new REAL * [total_dim];
        for(i1 =0; i1 < total_dim; i1++)
        {
            g[i1] = new REAL [3];
        }
        
        i_dim = 0;
  
        for (i_point = 0; i_point < n_points; i_point++)
        {
            // derivative g1, g2, g3 rwt x of atom i_point
            d_lm_x = 2.0*co[0]*(X_in[i_point][0]*co[0]
                    +X_in[i_point][1]*co[1]
                    +X_in[i_point][2]*co[2]);

            g[i_dim][0] = -(2.0*co[0]*X_in[i_point][0]
                          +co[1]*X_in[i_point][1]
                          +co[2]*X_in[i_point][2])       
                          +d_lm_x*co[0];
    
            g[i_dim][1] = -co[0]*X_in[i_point][1]+d_lm_x*co[1];

            g[i_dim][2] = -co[0]*X_in[i_point][2]+d_lm_x*co[2];
            
            // derivative g1, g2, g3 rwt y of atom i_point
   
            i_dim ++;

            d_lm_y = 2.0*co[1]*(X_in[i_point][0]*co[0]
                    +X_in[i_point][1]*co[1]
                    +X_in[i_point][2]*co[2]);

            g[i_dim][0] = -co[1]*X_in[i_point][0]       
                    +d_lm_y*co[0];

            g[i_dim][1] = -(co[0]*X_in[i_point][0]
                          +2.0*co[1]*X_in[i_point][1]
                          +co[2]*X_in[i_point][2])
                          +d_lm_y*co[1];

            g[i_dim][2] = -co[1]*X_in[i_point][2]+d_lm_y*co[2];

    
            // derivative g1, g2, g3 rwt z of of atom i_point
   
            i_dim ++;

            d_lm_z = 2.0*co[2]*(X_in[i_point][0]*co[0]
                    +X_in[i_point][1]*co[1]
                    +X_in[i_point][2]*co[2]);

            g[i_dim][0] = -co[2]*X_in[i_point][0]       
                          +d_lm_z*co[0];

            g[i_dim][1] = - co[2]*X_in[i_point][1]
                         +d_lm_z*co[1];
    
            g[i_dim][2] = -(co[0]*X_in[i_point][0]
                          +co[1]*X_in[i_point][1]
                          +2.0*co[2]*X_in[i_point][2])
                          +d_lm_z*co[2];
            
            i_dim++;
        }
        
        
        // 5. Finally, get the derivatives of coefficients wrt 
        //   atomic coordinates

        for (i1 = 0; i1 < total_dim; i1++)
        {    
            MatMultip(dim, inv_A, g[i1], d_co[i1]);
        }
        
        
        // 6. Release memory
        
        for (i1 =0; i1 < dim; i1++)
        {
            delete [] A[i1];
            A[i1] = 0;
        }
        delete [] A;
        A = 0;
        
        for (i1 =0; i1 < dim; i1++)
        {
            delete [] inv_A[i1];
            inv_A[i1] = 0;
        }   
        delete [] inv_A;
        inv_A = 0;

        delete []  eigenValue_A;
        eigenValue_A = 0;

        for(i1 =0; i1 < dim; i1++)
        {
            delete [] eigenVect_A[i1];
            eigenVect_A[i1] = 0;
        }
        delete [] eigenVect_A;
        eigenVect_A = 0;
  
        delete [] rel_eigenValue_A;
        rel_eigenValue_A = 0;
        delete [] inv_rel_eigenValue_A;
        inv_rel_eigenValue_A = 0;

        for(i1 =0; i1 < total_dim; i1++)
        {
            delete [] g[i1];
            g[i1] = 0;
        }   
        delete [] g;
        g = 0;
     
    }
    
  
    
    
    // Random number generator
    // Generate a random number between 0.0 to 1.0
    extern REAL GetRand()
    {
        const int     IMUL   = 314159269;
        const int     IIMUL  = 14156;
        const int     IADD   = 453806245;
        const int     MASK   = 2147483647;
        const double  SCALE  = 0.46566112873E-9;
        unsigned long int mySeed;

        mySeed = (unsigned)std::time(0)+rand()*IIMUL;

        mySeed = (mySeed * IMUL + IADD) & MASK;

        return (mySeed*SCALE);

    }
    
    // For hashing code of COD atom types
    extern void initPrimeTab(std::vector<int> & tPrimeTab,
                             std::string tLibmolTabDir)
    {
        // should use a env variable related path
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string primName = clibMonDir + "primes.table";
        //std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        std::string primName = tLibmolTabDir + "/primes.table";
        // std::cout << primName << std::endl;
        
        std::ifstream primFile;
        primFile.open(primName.c_str(), std::ios::in);
        
        std::string tRecord="";
        
        if (primFile.is_open())
        {
            while(!primFile.eof())
            {
                std::getline(primFile, tRecord);
                tPrimeTab.push_back(StrToInt(TrimSpaces(tRecord)));
            }
            
            primFile.close();
        }
        else
        {
            std::cout << "can not open " << primName << " for hashing table "
                      << std::endl;
        }
    }
    
    // Chemistry-related 

    extern void initOrgTable(std::vector<std::string> & tOrgTab)
    {
        ID orgSet[] = {"AT", "At", "at", "B", "b", "BR", "Br", "br", 
                       "C", "c", "CL", "Cl", "cl", "F", "f", "H", "h",
                       "I", "i", "N","n",  "O", "o", "P", "p", "S", "s", 
                       "SE", "Se", "se"};
        tOrgTab.assign(orgSet, orgSet+30);
        
    }
    
    extern bool isOrganc(std::vector<ID> & tOrgTab, ID tID)
    {
        std::vector<ID>::iterator iFind = std::find(tOrgTab.begin(), tOrgTab.end(), tID);
        if (iFind != tOrgTab.end())
        {
            return true;
        }
        
        return false;
        
    }
    
    extern void initAminoAcidTab(std::vector<ID> & tAATab)
    {
        ID aminoAcids [] = {"ALA", "ARG", "ASN", "ASP", "CYS",
                            "GLN", "GLU", "GLY", "HIS", "ILE",
                            "LEU", "LYS", "MET", "PHE", "PRO",
                            "SER", "THR", "TRP", "TYR", "VAL"};
        
        tAATab.assign(aminoAcids, aminoAcids+19);
    }
    extern bool isAminoAcid(std::vector<ID> & tAATab, ID tID)
    {
        std::vector<ID>::iterator iFind = std::find(tAATab.begin(), tAATab.end(), tID);
        if (iFind != tAATab.end())
        {
            return true;
        }
        
        return false;
    }
    
    extern void initMetalTab(std::vector<ID> & tMeTab)
    {
        ID metals[] = {"Li", "li", "Na", "na", "K",  "k",  "Rb", "rb", "Cs", "cs", "Fr", "fr",
                       "Be", "be", "Mg", "mg", "Ca", "ca", "Sr", "sr", "Ba", "ba", "Ra", "ra",
                       "Sc", "sc", "Y",  "y",
                       "Si", "si", "Ge", "ge", "As", "as", "Sb", "sb", "Te", "te", "Po", "po", 
                       "Ti", "ti", "Zr", "zr", "Hf", "hf", "Rf", "rf",
                       "V",  "v"   "Nb", "nb", "Ta", "ta", "Db", "db", 
                       "Cr", "cr", "Mo", "mo", "W",  "w",  "Sg", "sg", 
                       "Mn", "mn", "Tc", "tc", "Re", "re", "Bh", "bh",  
                       "Fe", "fe", "Ru", "ru", "Os", "os", "Hs", "hs",   
                       "Co", "co", "Rh", "rh", "Ir", "ir", "Mt", "mt",  
                       "Ni", "ni", "Pd", "pd", "Pt", "pt", "Ds", "ds",  
                       "Cu", "cu", "Ag", "ag", "Au", "au", "Rg", "rg",   
                       "Zn", "zn", "Cd", "cd", "Hg", "hg",   
                       "Al", "al", "Ga", "ga", "In", "in", "Ti", "ti", 
                       "Sn", "sn", "Pb", "pb", "Bi", "bi", "Pu", "pu", "Nd", "nd", "Ce", "ce",
                       "La", "la","Pr", "pr", "Pm", "pm", "Sm", "sm", "Eu", "eu", "Gd", "gd", 
                       "Tb", "tb", "Dy", "dy", "Ho", "ho", "Er", "er", "Tm", "tm", "Yb", "yb"
                       "Lu", "lu", "Ac", "ac", "Th", "th", "Pa", "pa", "U", "u", "Np", "np",
                       "Am", "am", "Cm", "cm", "Bk", "bk", "Cf", "cf", "Es", "es", "Fm", "fm",
                       "Md", "md", "No", "no", "Lr", "lr"};
        
        tMeTab.assign(metals, metals+181);
        /*
        for (std::vector<ID>::iterator iM=tMeTab.begin();
                iM !=tMeTab.end(); iM++)
        {
            std::cout << "Metal " << *iM << std::endl;
        }
         */
    }
    
    
    extern bool isMetal(std::vector<ID> & tMeTab, ID tID)
    {
         
        std::vector<ID>::iterator iFind = std::find(tMeTab.begin(), tMeTab.end(), tID);
        if (iFind != tMeTab.end())
        {
            return true;
        }
        
        return false;
    }
    
    extern void fromIdToChemType(ID tId, ID & tChemType)
    {
        unsigned j=0;
        for (unsigned i=0; i < tId.size(); i++)
        {
            char tCh=tId[i];
            if (!std::isdigit(tCh))
            {
                std::string tStr(tId.substr(i, 1));
                if(j==0)
                {
                    StrUpper(tStr);
                }
                else
                {
                    StrLower(tStr);
                }
                tChemType.append(TrimSpaces(tStr));
                j++;
            }
        }
    }
    
    extern REAL StrToOrder(std::string  & tStrOrder)
    {
        REAL tOrder = -1.0;
        
        StrUpper(tStrOrder);
        
        std::string a4 = tStrOrder.substr(0,4);
        if (a4.find("SING") !=std::string::npos )
        {
            tOrder = 1.0;
        }
        else if (a4.find("DOUB") !=std::string::npos)
        {
            tOrder = 2.0;
        }
        else if (a4.find("TRIP") !=std::string::npos)
        {
            tOrder = 3.0;
        }        
        else if (a4.find("AROM") !=std::string::npos)
        {
            //tOrder = 4.0;
            tOrder = 1.5;
        }
        else if (a4.find("DELO") !=std::string::npos)
        {
            tOrder = 1.5;
        }
        else if (a4.find("BOTH") !=std::string::npos)
        {
            tOrder = 5.0;
        }
        else if (a4.find("META") !=std::string::npos)
        {
            tOrder = 10.0;
        }
        else if (a4.find(".") !=std::string::npos)
        {
            tOrder = 100.0;
        }
        // what about "single or aromatic"->6, "double or aromatic"->7
        // and "any" -> 8
        return tOrder;
        
    }
    
    extern REAL StrToOrder2(std::string  & tStrOrder)
    {
        REAL tOrder = -1.0;
        
        StrUpper(tStrOrder);
        std::string aStr;
        if (tStrOrder.size() >=4)
        {
            aStr = tStrOrder.substr(0,4);
        }
        else
        {
            aStr = tStrOrder;
        }
        
        
        if (aStr.find("SING") !=std::string::npos || aStr.find("1") !=std::string::npos )
        {
            tOrder = 1.0;
        }
        else if (aStr.find("DOUB") !=std::string::npos || aStr.find("2") !=std::string::npos)
        {
            tOrder = 2.0;
        }
        else if (aStr.find("TRIP") !=std::string::npos || aStr.find("3") !=std::string::npos)
        {
            tOrder = 3.0;
        }        
        else if (aStr.find("AROM") !=std::string::npos || aStr.find("AR") !=std::string::npos)
        {
            //tOrder = 4.0;
            tOrder = 1.5;
        }
        else if (aStr.find("DELO") !=std::string::npos)
        {
            tOrder = 1.5;
        }
        else if (aStr.find("BOTH") !=std::string::npos)
        {
            tOrder = 5.0;
        }
        else if (aStr.find("META") !=std::string::npos)
        {
            tOrder = 10.0;
        }
        else if (aStr.find(".") !=std::string::npos || aStr.find("UN") !=std::string::npos)
        {
            tOrder = 100.0;
        }
        // what about "single or aromatic"->6, "double or aromatic"->7
        // and "any" -> 8
        return tOrder;
        
    }
    
    extern void OrderToStr(REAL tOrder, std::string  & sOrder)
    {
        if (tOrder == 1)
        {
            sOrder = "single";   
        }
        else if (tOrder == 2)
        {
            sOrder = "double";
        }
        else if (tOrder == 3)
        {
            sOrder = "trip";
        }
        else if (tOrder == 4 || tOrder==1.5)
        {
            sOrder = "aromatic";
        }      
        else if (tOrder == 5)
        {
            sOrder = "both";
        } 
        else if (tOrder == 9)
        {
            sOrder = "deloc";
        } 
        else if (tOrder == 10)
        {
            sOrder = "metal";
        }
        else if (tOrder == 100)
        {
            sOrder = ".";
        }
        
    }
    
    extern void OrderStrToStr(std::string & tOrder, std::string  & sOrder)
    {
        if (tOrder == "1")
        {
            sOrder = "single";   
        }
        else if (tOrder == "2")
        {
            sOrder = "double";
        }
        else if (tOrder == "3")
        {
            sOrder = "trip";
        }
        else if (tOrder == "4" || tOrder=="1.5" )
        {
            sOrder = "aromatic";
        }      
        else if (tOrder == "5")
        {
            sOrder = "both";
        } 
        else if (tOrder == "9")
        {
            sOrder = "deloc";
        } 
        else if (tOrder == "10")
        {
            sOrder = "metal";
        }
        else
        {
            sOrder = ".";
        }
        
    }
    
    extern void unifyStrForOrder(std::string & tOrder)
    {
        StrLower(tOrder);
        
        if (tOrder.size() ==4)
        {
            if (tOrder.find("sing") !=std::string::npos)
            {
                tOrder = "single";
            }
            else if (tOrder.find("doub") !=std::string::npos)
            {
                tOrder = "double";
            }
        }
    }
    
    // Chiral center conversion (mol/sdf file)
    extern void ChiToStr(int & tCIdx, std::string & tCStr)
    {
        if (tCIdx==2)
        {
            tCStr="negative";
        }
        else if (tCIdx==1)
        {
            tCStr="positive";
        }
        else
        {
            tCStr="both";
        }
    }
    
    // formal charge conversion (mol/sdf file)
    extern REAL strToCharge(std::string & tStr)
    {
        REAL aCharge =0.0;
        
        if (tStr.compare("1")==0)
        {
            aCharge = 3.0;
        }
        else if (tStr.compare("2")==0)
        {
            aCharge = 2.0;
        }
        else if (tStr.compare("3")==0)
        {
            aCharge = 1.0;
        }
        else if (tStr.compare("4")==0)
        {
            // doublet radical ?
            aCharge = 0.0;
        }
        else if (tStr.compare("5")==0)
        {
            aCharge = -1.0;
        }
        else if (tStr.compare("6")==0)
        {
            aCharge = -2.0;
        }
        else if (tStr.compare("7")==0)
        {
            aCharge = -3.0;
        }
        
        return aCharge;
    }
    
    // Symmetry-related functions
    
    extern void StrToSymmOps(std::vector<std::string>           & tStrs, 
                             std::vector<std::vector<REAL> >    & tMat)
    {
        for (std::vector<std::string>::iterator iS=tStrs.begin();
                iS !=tStrs.end(); iS++)
        {
            std::vector<REAL> tV;
            tV.push_back(0.0);
            tV.push_back(0.0);
            tV.push_back(0.0);
            tV.push_back(0.0);
            StrToSymmOneRow2(*iS, tV);
            tMat.push_back(tV);
        }
        
        std::vector<REAL> tVe;
        tVe.push_back(0.0);
        tVe.push_back(0.0);
        tVe.push_back(0.0);
        tVe.push_back(1.0);
        tMat.push_back(tVe);
    }
    
    void StrToSymmOneRow(std::string       & tStr,
                                std::vector<REAL> & tRow)
    {
        StrLower(tStr);
        std::string tB="";
        for (unsigned i=0; i < tStr.size(); i++)
        {
            if (tStr[i]=='-')
            {
                if ((int(tB.size()==0)))
                {
                    tB.append("-");
                }
            }
            else if (tStr[i]=='x')
            {
                tB.append("1");
                tRow[0]=StrToReal(tB);
                tB.clear();
            }
            else if (tStr[i]=='y')
            {
                tB.append("1");
                tRow[1]=StrToReal(tB);
                tB.clear();
            }
            else if (tStr[i]=='z')
            {
                tB.append("1");
                tRow[2]=StrToReal(tB);
                tB.clear();
            }
            else if(std::isdigit(tStr[i]))
            {
                if (tB.find("/") != std::string::npos)
                {
                    tB +=tStr[i];
                    tRow[3] = StrFractToReal(tB);
                    tB.clear();
                }
                else 
                {
                    tB+=tStr[i];
                }
            }
            else if (tStr[i]=='/')
            {
                tB.append("/");
            }
        }
    }
    
    extern void StrToSymmOneRow2(std::string         & tStr,
                          std::vector<double> & tRow)
    { 
        StrLower(tStr);
        std::vector<std::string> tBuf;
        if (tStr.find("+") != std::string::npos)
        {
           StrTokenize(tStr, tBuf, '+');
           for (unsigned j=0; j < tBuf.size(); j++)
           {
              if(tBuf[j].find("/") != std::string::npos)
              {
                 tRow[3] = StrFractToReal(tBuf[j]);
              }
              else if (tBuf[j].find(".") != std::string::npos)
              {
                 tRow[3] = StrToReal(tBuf[j]);
              }
              else if (tBuf[j].find("x")!= std::string::npos)
              {
                 if (tBuf[j].find("-")!= std::string::npos)
                 {
                    tRow[0]=-1.0;
                 }
                 else
                 {
                    tRow[0]=1.0;
                 }
              }
              else if (tBuf[j].find("y")!= std::string::npos)
              {
                 if (tBuf[j].find("-")!= std::string::npos)
                 {
                    tRow[1]=-1.0;
                 }
                 else
                 {
                    tRow[1]=1.0;
                 }
              }
              else if (tBuf[j].find("z")!= std::string::npos)
              {
                 if (tBuf[j].find("-")!= std::string::npos)
                 {
                    tRow[2]=-1.0;
                 }
                 else
                 {
                    tRow[2]=1.0;
                 }
              }
           }
        }
        else if (tStr.find("-") != std::string::npos)
        {
           std::string  tB="";
           for (unsigned i=0; i < tStr.size(); i++)
           {
              if (tStr[i]=='-')
              {
                  if (tB.size()==0)
                  {
                    tB.append("-");
                  }
                  else
                  {
                     if (tB.find("/") != std::string::npos)
                     {
                         tRow[3] = StrFractToReal(tB);
                         tB.clear();
                     }
                     else if (tB.find(".") != std::string::npos)
                     {
                         tRow[3] = StrToReal(tB);
                         tB.clear();
                     }
                     tB.append("-");
                  }
               }
               else if (tStr[i]=='x')
               {
                   tB.append("1");
                   tRow[0]=StrToReal(tB);
                   tB.clear();
               }
               else if (tStr[i]=='y')
               {
                   tB.append("1");
                   tRow[1]=StrToReal(tB);
                   tB.clear();
               }
               else if (tStr[i]=='z')
               { 
                   tB.append("1");
                   tRow[2]=StrToReal(tB);
                   tB.clear();
               }
               else if (i==(tStr.size()-1))
               {
                  tB.push_back(tStr[i]);
                  if (tB.find("/") != std::string::npos)
                  {
                     tRow[3] = StrFractToReal(tB);
                     tB.clear();
                  }
                  else if (tB.find(".") != std::string::npos)
                  {
                     tRow[3] = StrToReal(tB);
                     tB.clear();
                  }
               }
               else
               {
                  tB.push_back(tStr[i]);
               }
           }
        }
        else
        {
           if (tStr.find("x") != std::string::npos)
           {
              tRow[0]=1.0;
           }
           else if (tStr.find("y") != std::string::npos)
           {
              tRow[1]=1.0;
           }
           else if (tStr.find("z") != std::string::npos)
           {
              tRow[2]=1.0;
           }
        }
    }
              
    extern void TranslateIntoUnitCell(std::vector<REAL> & tInitFracX,
                                    std::vector<REAL> & tFinFracX)
    {
        tFinFracX.clear();
        
        for (std::vector<REAL>::iterator iIX=tInitFracX.begin();
                iIX !=tInitFracX.end(); iIX++)
        {
            if ((*iIX) >=1.0 || (*iIX) <=0.0 )
            {
                tFinFracX.push_back(*iIX-floor((*iIX)));
            }
            else
            {
                tFinFracX.push_back((*iIX));
            }
        }
    }
    
    extern void TranslateIntoUnitCell(std::vector<REAL> & tFracX)
    {
        std::vector<REAL> tmpFracX;
        for (std::vector<REAL>::iterator iIX=tFracX.begin();
                iIX !=tFracX.end(); iIX++)
        {
            if ((*iIX) >=1.0 || (*iIX) <=0.0 )
            {
                tmpFracX.push_back(*iIX-floor((*iIX)));
            }
            else
            {
                tmpFracX.push_back((*iIX));
            }
        }
        
        tFracX.clear();
        
        for (std::vector<REAL>::iterator iX=tmpFracX.begin();
                iX !=tmpFracX.end(); iX++)
        {
            tFracX.push_back(*iX);
        }
        
    }
    
    extern void FractToOrtho(std::vector<REAL> & tFractCoords,
                             std::vector<REAL> & tOrthoCoords,
                             REAL a, REAL b, REAL c,
                             REAL alpha, REAL beta, REAL gamma)
    {
        if (tOrthoCoords.size()==0)
        {
            for (unsigned i=0; i < tFractCoords.size(); i++)
            {
                tOrthoCoords.push_back(0.0);
            }
        }
        
        REAL coA   = cos(alpha*PI180);
        REAL coB   = cos(beta*PI180);
        REAL coG   = cos(gamma*PI180);
        REAL siG   = sin(gamma*PI180);
        
        
        tOrthoCoords[0] = tFractCoords[0]*a + tFractCoords[1]*b*coG + tFractCoords[2]*c*coB;
        tOrthoCoords[1] = tFractCoords[1]*b*siG + tFractCoords[2]*c*(coA-coB*coG)/siG;
        tOrthoCoords[2] = tFractCoords[2]*c*sqrt(pow(siG,2.0) - pow(coB,2.0) - pow(coA,2.0) + 2*coA*coB*coG)/siG;
        
        
    }
    
    extern void OrthoToFract(std::vector<REAL> tOrthoCoords,
                             std::vector<REAL> tFractCoords,
                             REAL a, REAL b, REAL c,
                             REAL alpha, REAL beta, REAL gamma)
    {
        if (tFractCoords.size()==0)
        {
            for (unsigned i=0; i < tOrthoCoords.size(); i++)
            {
                tFractCoords.push_back(0.0);
            }
        }
        
        REAL coA   = cos(alpha*PI180);
        REAL coB   = cos(beta*PI180);
        REAL coG   = cos(gamma*PI180);
        REAL siG   = sin(gamma*PI180);
        REAL tgG   = coG/siG;
        REAL reV   = sqrt(pow(siG,2.0) - pow(coB,2.0) - pow(coA, 2.0) + 2*b*c*coA*coB*coG);
        
        tFractCoords[0] = tOrthoCoords[0]/a - tOrthoCoords[1]*tgG/a + tOrthoCoords[2]*tgG*(coA-coB*coG)/(a*reV);
        tFractCoords[1] = tOrthoCoords[1]/(b*siG) - tOrthoCoords[2]*(coA - coB*coG)/(b*reV*siG);
        tFractCoords[2] = tOrthoCoords[2]*siG/(c*reV);
        
    }
    
    // executing external program using command line arguments
    // a simple version 
    extern void cmdExecute(std::string & tCom)
    {
        /*
        pid_t  pid=fork();
        int    status;

        if (pid < 0) 
        {     
            // fork a child process           
            std::cout << "*** ERROR: forking child process failed for command lines: " 
                      << std::endl << tCom << std::endl;
            exit(1);
        }
        else if (pid == 0)
        {
            // for the child process:         
            system(tCom.c_str());              // execute the command  
        }
        else
        {                                        // for the parent:      
            while (wait(&status) != pid);        // wait for completion  
        }
        
        */

    }
    
    extern std::string getLibmolDir()
    {
        return NullString; 
    }
    
    // The class for graph isomorphism 
    isomorGraph::isomorGraph()
    {
    }
    
    isomorGraph::~isomorGraph()
    {
    }
    
    void isomorGraph::isomorSubGraph(Graph & tSubGraph, 
                                     Graph & tGraph, 
                                     int     tMode, 
                                     std::vector<std::map<int,int> > & tOutMatchs)
    {  
        
        std::vector<std::vector<int> > adjA, adjB, M0, tM;
        
        setInitMatrixs(adjA, adjB, M0, tMode, tSubGraph, tGraph);        
        
        std::vector<int> usedCols;
        
        matCopyInt(M0, tM);
        
        unsigned curRow = 0;
        
        bool tDone = false;
        recurseMandIsomor(adjA, adjB, M0, usedCols, curRow, tM, 
                          tOutMatchs, tMode, tDone);
        
    }
    
    void isomorGraph::copyNodes(std::map<ID,ID> & tNodeI, 
                                std::map<ID,ID> & tNodeD)
    {
        for(std::map<ID, ID>::iterator iMI=tNodeI.begin();
                iMI !=tNodeI.end(); iMI++)
        {
            tNodeD[iMI->first] = iMI->second; 
        }
    }
    
    void isomorGraph::copyAdjacencies(std::vector<int> & tConnI, 
                                      std::vector<int> & tConnD)
    {
        for (std::vector<int>::iterator iC=tConnI.begin();
                iC !=tConnI.end(); iC++)
        {
            tConnD.push_back(*iC);
        }
    }
    
    void isomorGraph::setOneGraph(FileName tFName, Graph& tG)
    {  
        std::ifstream FofG(tFName);
        
        if (FofG.is_open())
        {
            bool lA = false;
            bool lC = false;
            
            std::string tRecord="";
            
            while(!FofG.eof())
            {   
                std::getline(FofG, tRecord);
                tRecord = TrimSpaces(tRecord);
                
                std::vector<std::string> tBuf;
                
                if (tRecord.find("ATOMS:") !=std::string::npos )
                {
                    lA=true;
                    lC=false;
                }
                else if (tRecord.find("CONNECTIONS:") !=std::string::npos )
                {
                    lA=false;
                    lC=true;
                }
                else if(lA)
                {
                    StrTokenize(tRecord, tBuf);
                    if (tBuf.size()==4)
                    {
                        int a = StrToInt(tBuf[0]);
                        tG.nodes[a]["elem"]=TrimSpaces(tBuf[1]);
                        tG.nodes[a]["name"]=TrimSpaces(tBuf[2]);
                        tG.nodes[a]["type"]=TrimSpaces(tBuf[3]);
                    }
                }
                else if(lC)
                {
                    StrTokenize(tRecord, tBuf);
                    if (tBuf.size()==2)
                    {
                        int c1=StrToInt(tBuf[0]);
                        int c2=StrToInt(tBuf[1]);
                        if(std::find(tG.adjacencies[c1].begin(), 
                                     tG.adjacencies[c1].end(), c2)==tG.adjacencies[c1].end())
                        {
                            tG.adjacencies[c1].push_back(c2);
                        }
                        
                        if(std::find(tG.adjacencies[c2].begin(), 
                                     tG.adjacencies[c2].end(), c1)==tG.adjacencies[c2].end())
                        {
                            tG.adjacencies[c2].push_back(c1);
                        }
                    }
                }
            }
            FofG.close();  
        }
        
        // Check
        /*
        if (tG.nodes.size() !=0)
        {
            for (std::map<int, std::map<ID,ID> >::iterator iN=tG.nodes.begin();
                    iN !=tG.nodes.end(); iN++)
            {
                std::cout << "For Node " << iN->first << std::endl
                          << "    its name: "    << iN->second["name"] << std::endl
                          << "    its element: " << iN->second["elem"] << std::endl
                          << "    its type: "    << iN->second["type"] << std::endl;
                std::cout << "The Node connects to the following Node: " << std::endl;
                for (std::vector<int>::iterator iC=tG.adjacencies[iN->first].begin();
                        iC !=tG.adjacencies[iN->first].end(); iC++)
                {
                    std::cout << "Node " << *iC << " of " 
                              << tG.nodes[*iC]["name"] << std::endl;
                } 
            }
        }
         */
    }
    
    void isomorGraph::reducedGraph(Graph& tFullG, Graph& tReducedG,
                                   std::map<int, int> & tNonHMapR2F,
                                   std::map<int, int> & tNonHMapF2R,
                                   std::map<int, std::vector<int> > & tLinkedH )
    {
        if (!tReducedG.nodes.empty())
        {
            tReducedG.nodes.clear();
        }
        
        if (!tReducedG.adjacencies.empty())
        {
            tReducedG.adjacencies.clear();
        }
        
        if (!tNonHMapR2F.empty())
        {
            tNonHMapR2F.clear();
        }
        
        if (!tNonHMapF2R.empty())
        {
            tNonHMapF2R.clear();
        }
        
        
        if (!tLinkedH.empty())
        {
            tLinkedH.clear();
        }
        
        unsigned nNonH=0;
        
        for (unsigned i=0; i < tFullG.nodes.size(); i++)
        {
            if (tFullG.nodes[i]["elem"].find("H")==std::string::npos)
            {
                copyNodes(tFullG.nodes[i], tReducedG.nodes[nNonH]);
                tNonHMapR2F[nNonH]=i;
                tNonHMapF2R[i] = nNonH;
                nNonH++;
            }
            else
            {
                if (tFullG.adjacencies[i].size()==1)
                {
                    tLinkedH[tFullG.adjacencies[i][0]].push_back(i);
                }
                else
                {
                    std::cout << "Error : H node " << tFullG.nodes[i]["name"]
                              << " connects " << tFullG.adjacencies[i].size()
                              << " nodes " << std::endl;
                    exit(1);
                   
                }
            }
        }
        
        // copy connections, be careful not to let H atoms into them 
        for (unsigned i=0; i < tFullG.adjacencies.size(); i++)
        {
            if (tFullG.nodes[i]["elem"].find("H")==std::string::npos)
            {
                for (std::vector<int>::iterator iC=tFullG.adjacencies[i].begin();
                        iC !=tFullG.adjacencies[i].end(); iC++)
                {
                    if (tFullG.nodes[*iC]["elem"].find("H")==std::string::npos)
                    {
                        tReducedG.adjacencies[tNonHMapF2R[i]].push_back(tNonHMapF2R[*iC]);
                    }
                }
            }
        }
        
        // Check 
        
        std::cout << "Original graph : " << std::endl;
        for (unsigned i=0; i < tFullG.nodes.size(); i++)
        {
            std::cout << "Node " << i << std::endl
                      << "name : " << tFullG.nodes[i]["name"] << std::endl
                      << "type : " << tFullG.nodes[i]["type"] << std::endl;
            
            std::cout << "This node connects to : " << std::endl;
            for (unsigned j=0; j < tFullG.adjacencies[i].size(); j++)
            {
                std::cout << "Node " << tFullG.nodes[tFullG.adjacencies[i][j]]["name"]
                          << std::endl;
            }
            std::cout << std::endl;
        }
        
        std::cout << "\nReduced graph now : " << std::endl;
        
        for (unsigned i=0; i < tReducedG.nodes.size(); i++)
        {
            std::cout << "Node " << i << std::endl
                      <<" name : " << tReducedG.nodes[i]["name"] << std::endl
                      <<" type : " << tReducedG.nodes[i]["type"] << std::endl
                      << "It corresponds to Node  " << tFullG.nodes[tNonHMapR2F[i]]["name"]
                      << " in original graph " << std::endl; 
            std::cout << "This node connects to : " << std::endl;
            for (unsigned j=0; j < tReducedG.adjacencies[i].size(); j++)
            {
                std::cout << "Node " << tReducedG.nodes[tReducedG.adjacencies[i][j]]["name"]
                          << std::endl;
            }
        }
        
        std::cout << "Nodes with H links " << std::endl;
        
        for (std::map<int, std::vector<int> >::iterator iLH=tLinkedH.begin();
                iLH !=tLinkedH.end(); iLH++)
        {
            std::cout << "Node " << tFullG.nodes[iLH->first]["name"] 
                      << "Linked with following H nodes: " << std::endl;
            
            for (std::vector<int>::iterator iH=iLH->second.begin();
                    iH!=iLH->second.end(); iH++)
            {
                std::cout << "Node " << tFullG.nodes[*iH]["name"] << std::endl;
            }
            
        }
        
    }
    
    void isomorGraph::setHLinks(Graph& tFullG1, 
                                std::map<int,int>& tNonHMapF2R_1, 
                                std::map<int,std::vector<int> > & tLinkedH1, 
                                Graph& tFullG2, 
                                std::map<int,int>& tNonHMapR2F_2, 
                                std::map<int,std::vector<int> >& tLinkedH2, 
                                std::map<int,int>& tOneSetR2RMap,
                                std::map<int,int>& tOneSetH2HMap)
    {
        // 
        for (std::map<int, std::vector<int> >::iterator iLH1=tLinkedH1.begin();
                iLH1 !=tLinkedH1.end(); iLH1++)
        {
            int nF2R_1 = tNonHMapF2R_1[iLH1->first];
            if (tOneSetR2RMap.find(nF2R_1) != tOneSetR2RMap.end())
            {
                int nR2F_2=tNonHMapR2F_2[tOneSetR2RMap[nF2R_1]];
                
                if (iLH1->second.size() ==tLinkedH2[nR2F_2].size())
                {
                    for (unsigned iH=0; iH<iLH1->second.size(); iH++)
                    {
                        tOneSetH2HMap[iLH1->second[iH]]= tLinkedH2[nR2F_2][iH];
                    }
                }
                else
                {
                    std::cout << "Numbers of H atoms are not consistent " << std::endl;
                    exit(1);
                }
            }
        }
    }
    
    void isomorGraph::recoverFullGraphMatch(Graph& tGraph1,
                                            std::map<int,int>& tNonHMapF2R_1,
                                            std::map<int,int>& tNonHMapR2F_1,
                                            std::map<int,std::vector<int> >& tLinkedH1,
                                            Graph& tGraph2,
                                            std::map<int,int>& tNonHMapF2R_2,
                                            std::map<int,int>& tNonHMapR2F_2,
                                            std::map<int,std::vector<int> >& tLinkedH2,
                                            std::map<int,int>    & tOneSetReducedOutMatch,
                                            std::map<int,int>    & tOneSetOutMatch)
    {
        std::map<int,int> aSetH2HMap;
        
        setHLinks(tGraph1, tNonHMapF2R_1, tLinkedH1, 
                  tGraph2, tNonHMapR2F_2, tLinkedH2, 
                  tOneSetReducedOutMatch, aSetH2HMap);
        
        // Non H
        for (std::map<int, int>::iterator iNonH=tOneSetReducedOutMatch.begin();
                iNonH !=tOneSetReducedOutMatch.end(); iNonH++)
        {
           tOneSetOutMatch[tNonHMapR2F_1[iNonH->first]] =  tNonHMapR2F_2[iNonH->second];
        }
        
        // H
        for (std::map<int, int>::iterator iH=aSetH2HMap.begin();
                iH !=aSetH2HMap.end(); iH++)
        {
            tOneSetOutMatch[iH->first] =  iH->second;
        }
        
        // Check
        /*
        std::cout << "Number of matched nodes : " << tOneSetOutMatch.size()
                  << std::endl;
        for (std::map<int, int>::iterator iM=tOneSetOutMatch.begin();
                iM != tOneSetOutMatch.end(); iM++)
        {
            std::cout << "Node " << tGraph1.nodes[iM->first]["name"] 
                      << " corresponds to Node " 
                      << tGraph2.nodes[iM->second]["name"] << std::endl;
        }
         */
        
        
    }
    
    void isomorGraph::setInitMatrixs(std::vector<std::vector<int> >& tAdjA, 
                                     std::vector<std::vector<int> >& tAdjB, 
                                     std::vector<std::vector<int> >& tM0, 
                                     int   tMode,
                                     Graph& tGraph1, Graph& tGraph2)
    {
        unsigned nN1 = tGraph1.nodes.size();
        unsigned nN2 = tGraph2.nodes.size();
        
        std::cout << "For Node A: " << std::endl;
        
        for (unsigned i=0; i < nN1; i++)
        {
            std::vector<int> tN;
            for (unsigned j=0; j < nN1; j++)
            {
                if (std::find(tGraph1.adjacencies[i].begin(), 
                              tGraph1.adjacencies[i].end(), j)
                        !=tGraph1.adjacencies[i].end())
                {
                    tN.push_back(1);
                    std::cout << "Node " << i << " connects Node " << j << std::endl;
                }
                else
                {
                    tN.push_back(0);
                }
            }
            tAdjA.push_back(tN);   
        }
        std::cout << "\n\n";
        
        
        std::cout << "For Node B: " << std::endl;
        for (unsigned i=0; i < nN2; i++)
        {
            std::vector<int> tN;
            for (unsigned j=0; j < nN2; j++)
            {
                if (std::find(tGraph2.adjacencies[i].begin(), 
                              tGraph2.adjacencies[i].end(), j)
                        !=tGraph2.adjacencies[i].end())
                {
                    tN.push_back(1);
                    std::cout << "Node " << i << " connects Node " << j << std::endl;
                }
                else
                {
                    tN.push_back(0);
                }
            }
            tAdjB.push_back(tN);
        }
        std::cout << "\n\n";
        
        if (tMode==1)
        {
            setExactMatch_M0(tM0, tGraph1, tGraph2);
            std::cout << "Number of rows in M0 " << tM0.size() << std::endl;
            std::cout << "Number of cols in M0 " << tM0[0].size() << std::endl;
        }
    }
    
    void isomorGraph::setExactMatch_M0(std::vector<std::vector<int> > & tM0,
                          Graph  & tGraph1, Graph & tGraph2)
    {
        unsigned nN1 = tGraph1.nodes.size();
        unsigned nN2 = tGraph2.nodes.size();
        
        for (unsigned i=0; i < nN1; i++)
        {
            std::vector<int> tN;
            for (unsigned j=0; j < nN2; j++)
            {
                if (tGraph1.nodes[i]["type"] == tGraph2.nodes[j]["type"])
                {
                    tN.push_back(1);
                    std::cout << "Node " << i << " and Node " << j 
                              << " have the same types "  << std::endl;
                    std::cout << "confirm : " << std::endl
                              << "Graph 1: " << tGraph1.nodes[i]["type"] << std::endl
                              << "Graph 2: " << tGraph2.nodes[j]["type"] << std::endl;
                    std::cout << "M0 " << i << " and " << tN.size()-1 
                              << "==" << tN[tN.size()-1] << std::endl;
                    std::cout << std::endl;
                }
                else
                {
                    tN.push_back(0);
                }
            }
            tM0.push_back(tN);
        }
        
        /*
        std::cout << "Matrix M0 is: " << std::endl;
        printMatrix(tM0);
        std::cout << std::endl;
        */
        
        int nMax=1;
        for (unsigned i=0; i < tM0.size(); i++)
        {
            int sumR=0;
            for (unsigned j=0; j < tM0[i].size(); j++)
            {
                sumR+=tM0[i][j];
            }
            if (sumR !=0)
            {  
                nMax = nMax*sumR;
                std::cout << "row " << i << " has " << sumR << " 1s " << std::endl; 
                std::cout << "nMax now " << nMax << std::endl;
            }
            else
            {
                std::cout << "Error : " << std::endl
                          << "row " << i << " does not has any match candidate" << std::endl; 
            }
        }
        std::cout << "Max number of tries : " <<  nMax << std::endl;
       
    }
    
    void isomorGraph::recurseMandIsomor(std::vector<std::vector<int> >  & tAdjA, 
                                        std::vector<std::vector<int> >  & tAdjB, 
                                        std::vector<std::vector<int> >  & tM0, 
                                        std::vector<int> & usedCols, int curRow, 
                                        std::vector<std::vector<int> >  & tM,
                                        std::vector<std::map<int,int> > & tOutMatch,
                                        int                               tMode,
                                        bool                            & tDone)
    {
        
        if (curRow==tM0.size())
        {
            //std::cout << "None Zero elements in One M are : " << std::endl;
            
            for (unsigned i=0; i < tM.size(); i++)
            {
                unsigned nNR=0;
                for (unsigned j=0; j < tM[i].size(); j++)
                {
                    if (tM[i][j]==1)
                    {
                        //std::cout << "None zero element at row " << i 
                        //          << " and col " << j << std::endl;
                        nNR++;
                    }
                    if (nNR > 1)
                    {
                        std::cout << "Error : row " << i 
                                  << " has more than one none zero elements "
                                  << std::endl;
                        exit(1);
                    }
                }
            }        
            // printMatrix(tM);
            checkIsomor(tAdjA, tAdjB, tM, tOutMatch);
            std::cout << "Number of match set " << tOutMatch.size() << std::endl;
            if (tMode==1 && tOutMatch.size() !=0)
            {
                tDone = true;
            }
            
        }
        else
        {
            //std::cout << "\ntM0 row size " << tM0.size() << std::endl 
            //          << "tM0 col size " << tM0[curRow].size() << std::endl;
            std::cout << "\nCurrent row Now " << curRow << std::endl;
            
            
            unsigned iCol=0;
            do 
            {
                //std::cout << "current col " << iCol << std::endl;
                //std::cout << "M0 [" << curRow << "][" << iCol 
                //          << "]=" << tM0[curRow][iCol] << std::endl;
                while(std::find(usedCols.begin(), usedCols.end(),iCol)==usedCols.end()
                      && iCol < tM0[curRow].size() && !tDone)
                {
                    //std::cout << "current row inside " << curRow << std::endl;
                    //std::cout << "current col inside " << iCol << std::endl;
                    if (tM0[curRow][iCol]==1)
                    {
                        tM[curRow][iCol] =1;
                        usedCols.push_back(iCol);
                        // std::cout << "tM " << curRow << "  " << iCol << " ==1 " << std::endl;
                        /*
                        std::cout << "Number of used cols " << usedCols.size() << std::endl;
                        for (unsigned iC=0; iC < usedCols.size(); iC++)
                        {
                            std::cout << "Col: " << usedCols[iC] << std::endl;
                        }
                        */    
                        for (unsigned j=0; j < tM[curRow].size(); j++)
                        {
                            if (j!=iCol)
                            {
                                tM[curRow][j]=0;
                                //std::cout << "set col " << j << " in row " << curRow 
                                //          << " to be zero " << std::endl;
                            }
                        }
                        curRow = curRow+1;
                        recurseMandIsomor(tAdjA, tAdjB, tM0, usedCols, 
                                          curRow, tM,  tOutMatch, tMode, tDone);
                        curRow = curRow-1;
                        usedCols.pop_back();
                    }
                    //std::cout << "row " << curRow << " and col " << iCol 
                    //          << " not used (zero M0 element) " << std::endl;
                    iCol++;
                }
           
                //std::cout << "row " << curRow << " and col " << iCol 
                //          << " not used (in used_list) " << std::endl;
                iCol++;
            }while(iCol <tM0[curRow].size() && !tDone);
        }
    }
    
    void isomorGraph::checkIsomor(std::vector<std::vector<int> >  & tAdjA, 
                                  std::vector<std::vector<int> >  & tAdjB, 
                                  std::vector<std::vector<int> >  & tM,
                                  std::vector<std::map<int,int> > & tOutMatch)
    {
        std::vector<std::vector<int> >  aM, transM, cM;
        
        //MxtAdjB
        matMultMatInt(tM, tAdjB, aM);
        //std::cout << "Matrix B: " << std::endl;
        //printMatrix(tAdjB);
        //std::cout << "Matrix (M X B) " << std::endl;
        //printMatrix(aM);
        
        //(MxtAdjB)^T
        unsigned nAM=aM.size();
        if (nAM >0)
        {
            unsigned mAM = aM[0].size();  
            for (unsigned i=0; i < mAM; i++ )
            {
                std::vector<int> aT;
                for (unsigned j=0; j < nAM; j++)
                {
                    aT.push_back(aM[j][i]);
                }
                transM.push_back(aT);
            }
        }
        //std::cout << "Matrix (M X B)^T " << std::endl;
        //printMatrix(transM);
        
        
        // M x (MxtAdjB)^T
        matMultMatInt(tM, transM, cM);
        
        //std::cout << "Matrix C: " << std::endl;
        //printMatrix(cM);
        
        //std::cout << "Matrix A: " << std::endl;
        //printMatrix(tAdjA);
        // Final check
        
        bool lIso=true;
        unsigned nRow = tAdjA.size();
        unsigned nRowCM= cM.size();
        //std::cout << "A with " << nRow   << " rows " << std::endl;
        //std::cout << "C with " << nRowCM << " rows " << std::endl;
        if (nRow > 0 && nRowCM==nRow)
        {
            unsigned nCol   = tAdjA[0].size();
            unsigned nColCM = cM[0].size();
            if (nCol==nColCM)
            {
                for (unsigned i=0; i < nRow; i++)
                {
                    for (unsigned j=0; j < nCol; j++)
                    {
                        if (tAdjA[i][j]!=cM[i][j])
                        {
                            //std::cout << "A[" << i << "][" << j << "]" 
                            //          << tAdjA[i][j] << std::endl;
                            //std::cout << "C[" << i << "][" << j << "]" 
                            //          << cM[i][j] << std::endl;
                            lIso = false;
                            break;
                        }
                    }
                }
            }
            else
            {
                lIso = false;
            }
        }
        else
        {
            lIso = false;
        }
        
        if (lIso)
        {
            std::map<int, int> aMatch;
            for (unsigned i=0; i < tM.size(); i++)
            {
                for (unsigned j=0; j < tM[0].size(); j++)
                {
                    if (tM[i][j]==1)
                    {
                        aMatch[i]=j;
                    }
                }
            }
            tOutMatch.push_back(aMatch);   
        }
    }
    
    void isomorGraph::graphMatch(FileName tSubGraphFName, 
                                 Graph &  tSubGraph,
                                 FileName tGraphFName,
                                 Graph &  tGraph,
                                 std::vector<std::map<int,int> > & tOutMatch,
                                 int tMode)
    {     
        
        Graph reducedSubG, reducedG;
        std::map<int, int> nonHMapR2F_SubG, nonHMapF2R_SubG;
        std::map<int, int> nonHMapR2F_G, nonHMapF2R_G;
        std::map<int, std::vector<int> > linkedHMap_subG, linkedHMap_G;
        std::vector<std::map<int,int> >  tOutMatchOnreducedG;
        
        setOneGraph(tSubGraphFName, tSubGraph);
        reducedGraph(tSubGraph, reducedSubG, nonHMapR2F_SubG, 
                     nonHMapF2R_SubG, linkedHMap_subG);
        
        setOneGraph(tGraphFName, tGraph);
        reducedGraph(tGraph, reducedG, 
                     nonHMapR2F_G, nonHMapF2R_G, linkedHMap_G);
        
        if (reducedSubG.nodes.size() !=0 && reducedSubG.adjacencies.size() !=0
                && reducedG.nodes.size() !=0 && reducedG.adjacencies.size() !=0)
        {
            isomorSubGraph(reducedSubG, reducedG, tMode, tOutMatchOnreducedG);
        }
        else
        {
            if (reducedSubG.nodes.size() ==0 || reducedSubG.adjacencies.size() ==0)
            {
                std::cout << "Check : " << tSubGraphFName
                          << " does not provide the correct info to build a graph"
                          << std::endl;
            }
            else if (reducedSubG.nodes.size() ==0 || reducedSubG.adjacencies.size() ==0)
            {
                std::cout << "Check : " << tGraphFName
                          << " does not provide the correct info to build a graph"
                          << std::endl;
            }
            else
            {
                std::cout << "Check: File format errors in  " << tSubGraphFName
                          << " or " << tGraphFName << std::endl;
            }
        }
        
        unsigned nMSets=tOutMatchOnreducedG.size();
        std::cout << "Number of matched sets of nodes (reduced set): " << nMSets << std::endl;
        if (nMSets !=0)
        {
            for (unsigned i=0; i < nMSets; i++ )
            {
                std::cout << "For matched node set " << i+1 << std::endl;
                for (std::map<int, int>::iterator iM=tOutMatchOnreducedG[i].begin(); 
                        iM !=tOutMatchOnreducedG[i].end(); iM++)
                {  
                    std::cout << "Atom " << reducedSubG.nodes[iM->first]["name"]  
                              << " corresponds to atom " 
                              << reducedG.nodes[iM->second]["name"] << std::endl;
                }
            }
            
            for (unsigned i=0; i < nMSets; i++ )
            {
                std::map<int, int> aSetMatchFull;
                recoverFullGraphMatch(tSubGraph, nonHMapF2R_SubG, 
                                      nonHMapR2F_SubG, linkedHMap_subG,
                                      tGraph, nonHMapF2R_G, 
                                      nonHMapR2F_G, linkedHMap_G,
                                      tOutMatchOnreducedG[i],
                                      aSetMatchFull);
                if (aSetMatchFull.size() !=0)
                {
                    tOutMatch.push_back(aSetMatchFull);
                }
            }
        }
        
        
        
        /*
        if (tSubGraph.nodes.size() !=0 && tSubGraph.adjacencies.size() !=0
                && tGraph.nodes.size() !=0 && tGraph.adjacencies.size() !=0)
        {
            isomorSubGraph(tSubGraph, tGraph, tMode, tOutMatch);
        }
        else
        {
            if (tSubGraph.nodes.size() ==0 || tSubGraph.adjacencies.size() ==0)
            {
                std::cout << "Check : " << tSubGraphFName
                          << " does not provide the correct info to build a graph"
                          << std::endl;
            }
            else if (tSubGraph.nodes.size() ==0 || tSubGraph.adjacencies.size() ==0)
            {
                std::cout << "Check : " << tGraphFName
                          << " does not provide the correct info to build a graph"
                          << std::endl;
            }
            else
            {
                std::cout << "Check: File format errors in  " << tSubGraphFName
                          << " or " << tGraphFName << std::endl;
            }
        }
        */ 
    }
    
    void isomorGraph::outputMatchedGraphs(Graph & tSubGraph, 
                                          Graph & tGraph, 
                                          std::vector<std::map<int,int> > & tOutMatch,
                                          int                               tMode,
                                          FileName    tOutMatchFileName)
    {
        std::ofstream outFile(tOutMatchFileName);
        
        if (outFile.is_open())
        {
            for (unsigned i=0; i < tOutMatch.size(); i++)
            {
                if (tMode==1)
                {
                    outFile << "#Exact match between two graphs" << std::endl;
                }
                else
                {
                    outFile << "#Match between two graphs, set " << i+1 << std::endl;
                }
                
                for (std::map<int,int>::iterator iAt=tOutMatch[i].begin();
                        iAt != tOutMatch[i].end(); iAt++)
                {
                    outFile << std::setw(10) << tSubGraph.nodes[iAt->first]["name"] 
                            << std::setw(10) << tGraph.nodes[iAt->second]["name"] 
                            << std::endl;
                }
            }
            outFile.close();
            
        }
        
    }
    
}
