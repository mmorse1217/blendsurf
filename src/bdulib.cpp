#include "bdulib.hpp"

// ---------------------------------------------------------------------- 
int BdULib::setup(istream& fin)
{
    assert(fin.good());
    _evec.resize(MAXVALENCE);
    char tmp[100];
    fin>>tmp;
    while(fin.eof()==false) {
        assert(strcmp(tmp, "valence")==0);
        int k;	 fin>>k;	 assert(k>=3 && k<MAXVALENCE);
        int m;	 fin>>m;
        int n;	 fin>>n;
        Entry& e = _evec[k];
        NumMatrix& U = e.U();
        U.resize(m,n);
        for(int i=0; i<m; i++)
            for(int j=0; j<n; j++)
                fin>>U(i,j);
        NumMatrix& IU = e.IU();
        IU.resize(n,m);
        for(int i=0; i<n; i++)
            for(int j=0; j<m; j++)
                fin>>IU(i,j);
        fin>>tmp;
    }
    return 0;
}
