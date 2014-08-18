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
    
    bool sortStrByLen(std::string first, std::string second)
    {
        return (int)first.size() > (int)second.size();
    }
    
    
    bool desSortMapKey(const sortMap& a ,const sortMap & b)
    {
        return a.key.length() > b.key.length();
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
    extern void initPrimeTab(std::vector<int> & tPrimeTab)
    {
        // should use a env variable related path
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string primName = clibMonDir + "primes.table";
        std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        std::string primName = clibMonDir + "/tables/primes.table";
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
        if (a4.find("SING") !=std::string::npos)
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
        else if (a4.find("BOTH") !=std::string::npos)
        {
            tOrder = 5.0;
        }
        else if (a4.find("DELO") !=std::string::npos)
        {
            tOrder = 9.0;
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
        else if (tOrder == 4)
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
        else if (tOrder == "4")
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
    
    
        
}
