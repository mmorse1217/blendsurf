#include "ccsubmatlib.hpp"

// ---------------------------------------------------------------------- 
int CCSubMatLib::setup(istream& fin)
{
  
    assert(fin.good());
    _interior.resize(MAXVALENCE);
    _boundary.resize(MAXVALENCE);
    _cornerconvex.resize(MAXVALENCE);
    _cornerconcave.resize(MAXVALENCE);
  
    char tmp[100];
    fin>>tmp;
    while(fin.eof()==false) {
      //        cout<<tmp<<endl;
        assert(strcmp(tmp, "interior")==0 || strcmp(tmp, "boundary")== 0 || 
               strcmp(tmp, "cornerconvex") == 0 || strcmp(tmp, "cornerconcave") == 0);
        if (strcmp(tmp, "interior") == 0){
            vector<Entry>& cur = _interior;
            //get valence
            int k;
            fin>>k;	 assert(k>=3 && k<=MAXVALENCE        );
            Entry& entry = cur[k];
            //get size
            int size;
            fin>>size;	 assert(size == 6*k+1);
            //read in SM : subdivision matrix
            NumMatrix& sm = entry.SM();
            sm.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>sm(i,j); //read row by row
            //read in IV: inverse of eigen vectors
            NumMatrix& iv = entry.IV();
            iv.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>iv(i,j); //read row by row
            //read in V : eigen vector matrix
            NumMatrix& v = entry.V();
            v.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>v(i,j);  //read row by row
            //read in D : eigen values
            NumVector& d = entry.D();
            d.resize(size);
            for(int i=0; i<size; i++)
                fin>>d(i);
      

        }
        else if (strcmp(tmp, "boundary") == 0){
            vector<Entry>& cur = _boundary;
            //get valence
            int k;
            fin>>k;	 assert(k>=2 && k<MAXVALENCE);
            Entry& entry = cur[k];
            //get size
            int size;
            fin>>size;	 assert(size == 6*k+3);
            //read in SM : subdivision matrix
            NumMatrix& sm = entry.SM();
            sm.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>sm(i,j); //read row by row
            //read in IV: inverse of eigen vectors
            NumMatrix& iv = entry.IV();
            iv.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>iv(i,j); //read row by row
            //read in V : eigen vector matrix
            NumMatrix& v = entry.V();
            v.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>v(i,j);  //read row by row
            //read in D : eigen values
            NumVector& d = entry.D();
            d.resize(size);
            for(int i=0; i<size; i++)
                fin>>d(i);
      
        }
        else if (strcmp(tmp, "cornerconvex") == 0){
            vector<Entry>& cur = _cornerconvex;
            //get valence
            int k;
            fin>>k;	 assert(k>=1 && k<MAXVALENCE);
            Entry& entry = cur[k];
            //get size
            int size;
            fin>>size;	 assert(size == 6*k+3);
            //read in SM : subdivision matrix
            NumMatrix& sm = entry.SM();
            sm.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>sm(i,j); //read row by row
            //read in IV: inverse of eigen vectors
            NumMatrix& iv = entry.IV();
            iv.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>iv(i,j); //read row by row
            //read in V : eigen vector matrix
            NumMatrix& v = entry.V();
            v.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>v(i,j);  //read row by row
            //read in D : eigen values
            NumVector& d = entry.D();
            d.resize(size);
            for(int i=0; i<size; i++)
                fin>>d(i);
      
        }
        else if (strcmp(tmp, "cornerconcave") == 0){
            vector<Entry>& cur = _cornerconcave;
            //get valence
            int k;
            fin>>k;	 assert(k>=2 && k<MAXVALENCE);
            Entry& entry = cur[k];
            //get size
            int size;
            fin>>size;	 assert(size == 6*k+3);
            //read in SM : subdivision matrix
            NumMatrix& sm = entry.SM();
            sm.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>sm(i,j); //read row by row
            //read in IV: inverse of eigen vectors
            NumMatrix& iv = entry.IV();
            iv.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>iv(i,j); //read row by row
            //read in V : eigen vector matrix
            NumMatrix& v = entry.V();
            v.resize(size,size);
            for(int i=0; i<size; i++)
                for(int j=0; j<size; j++)
                    fin>>v(i,j);  //read row by row
            //read in D : eigen values
            NumVector& d = entry.D();
            d.resize(size);
            for(int i=0; i<size; i++)
                fin>>d(i);
        }

        fin>>tmp;
    }
    return 0;
}


CCSubMatLib::Entry& CCSubMatLib::chooseMatrix(int flag, int k)
{
    if(flag==INTERIOR) 
        return _interior[k];
  
    else if(flag==BOUNDARY) 
        return _boundary[k];
  
    else if (flag==CORNERCONVEX)
        return _cornerconvex[k];
  
    else if (flag == CORNERCONCAVE)
        return _cornerconcave[k];
    else assert(0);
  
    //control can not reach here anyway
    return _interior[k];
}

