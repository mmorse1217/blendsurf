#include "ccsubmatlib.hpp"
#include "ccsurf.hpp"
#include "bdsurf.hpp"

#include "ballviewer.hpp"
#include "ccsurfobj.hpp"
#include "bdsurfobj.hpp"
#include "evalfromw.hpp"

int optionsCreate(const char* optfile, map<string,string>& options)
{
    options.clear();
    ifstream fin(optfile); assert(fin.good());
    string name;  fin>>name;
    while(fin.good()) {
        char cont[100];   fin.getline(cont, 99);
        options[name] = string(cont);
        fin>>name;
    }
    fin.close();
    return 0;
}


// ---------------------------------------------------------------------- 
int main(int argc, char** argv)
{
  
  
    assert( argc == 2 );

    map<string,string> opts;
    optionsCreate(argv[1], opts);
  
    map<string,string>::iterator mi;

    int objtype;
    mi = opts.find("-objtype"); assert(mi!=opts.end());  { istringstream ss((*mi).second); ss>>objtype; }
    GeoObject* obj = NULL;

    if(objtype==0) { //dealing with catmull-clark surface
        //------------------------

        //submatlib stores a subdivision matrix library for a two ring around a vertex with any valance.
        char filea[100];	 mi = opts.find("-ccsurf_submatlibfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filea; }
        ifstream tmpa(filea);
        CCSubMatLib* submatlib = new CCSubMatLib();
        cerr<<"submatlib setup"<<endl;
        submatlib->setup(tmpa);
    
        //setup the surface(read the mesh, share, etc)
        char fileb[100];	 mi = opts.find("-ccsurf_meshfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>fileb; }
        ifstream tmpb(fileb);

        CCSurf* ccsurf = new CCSurf( );
        double flats;  mi = opts.find("-ccsurf_flats"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>flats;} 
        cerr<<"ccsurf setup"<<endl;
        ccsurf->setup(flats, tmpb,0,Point3(1.0), (void*)submatlib);
  
        int lvl;  
        mi = opts.find("-ccsurf_renderlvl"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>lvl; }

        int tg = lvl;
        if( ccsurf->validpos(tg)==false ) {
            //level of subdivion not enough for renderer.need tos subdivide more
            cout<<"valid pos is false"<<endl;
            int lv = tg;
            while(ccsurf->validpos(lv)==false && lv>0)
                lv--; //find the nearest valid one
            assert(ccsurf->validpos(lv)==true);
      
            for(int i=lv; i<tg; i++) 
                ccsurf->subdivide(i, *submatlib); //subdivide if necessary      
        }
        cerr<<"ccsurfobj setup begin"<<endl;
        //surfobj takes care of rendering.
        obj = new CCSurfObj(ccsurf, lvl);
        cerr<<"ccsurfobj setup end"<<endl;
    } 

    else if(objtype==1) { //the blended surface
        //------------------------
        char filea[100];	 mi = opts.find("-bdsurf_submatlibfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filea; }
        ifstream tmpa(filea);
        CCSubMatLib* submatlib = new CCSubMatLib();
        cerr<<"submatlib setup"<<endl;
        submatlib->setup(tmpa);
        char fileb[100];	 mi = opts.find("-bdsurf_bdulibfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>fileb; }
        ifstream tmpb(fileb);
        BdULib* bdulib = new BdULib();
        cerr<<"bdulib setup"<<endl;
        bdulib->setup(tmpb);
        char filec[100];	 mi = opts.find("-bdsurf_meshfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filec; }
        ifstream tmpc(filec);
    
        int cl;     mi = opts.find("-bdsurf_ctrllvl");   assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>cl; } assert(cl>=0 && cl<=2);
        int pc;	    mi = opts.find("-bdsurf_pouctrl");   assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>pc; } assert(pc==0 || pc==1 || pc==2);

        int ct;     mi = opts.find("-bdsurf_chttyp");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>ct; } 
        int bt;     mi = opts.find("-bdsurf_bsstyp");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>bt; }
        int pp;	    mi = opts.find("-bdsurf_stpppg");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>pp; }
        double lb;	mi = opts.find("-bdsurf_lb");        assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>lb; }
        double ub;	mi = opts.find("-bdsurf_ub");        assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>ub; }
        int gl;	    mi = opts.find("-bdsurf_indepboun"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>gl; }
        double fs;	mi = opts.find("-bdsurf_flats");     assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>fs; }
        int cot;    mi = opts.find("-bdsurf_concat");    assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>cot; }
        int sd;	    mi = opts.find("-bdsurf_spdeg");   assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>sd; }
	int poudeg; mi = opts.find("-bdsurf_poubsdeg");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>poudeg; }

	int loadW;  mi = opts.find("-bdsurf_loadw");  if(mi==opts.end()) loadW = 0; else 
	  { istringstream ss((*mi).second); ss>>loadW; }

	mi = opts.find("-bdsurf_matdir"); if(mi!=opts.end()) { istringstream ss((*mi).second); ss>> MatricesDir; }


        BdSurf* bdsurf = new BdSurf(  );
        bdsurf->setParams(cl,pc,ct,bt,pp,lb,ub, gl, fs, cot, sd ,poudeg, loadW );
        cerr<<"bdsurf setup"<<endl;
        bdsurf->setup(tmpc,0,Point3(1.0),(void*)submatlib,(void*)bdulib);
        //ofstream fout("bddump.dat"); bdsurf->dump(fout); fout.close();
    
        int lvl; mi = opts.find("-bdsurf_renderlvl");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>lvl; }
        int gen; mi = opts.find("-bdsurf_rendergen");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>gen; }
        int alt; mi = opts.find("-bdsurf_alt");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>alt; }

	cerr<<"bdsurfobj setup"<<endl;
        char filetex[100];  mi = opts.find("-bdsurfobj_texfile"); assert(mi!=opts.end());   { istringstream ss((*mi).second); ss>>filetex; }
           obj = new BdSurfObj(bdsurf, lvl, gen, alt, filetex);
	   //bdsurf->evalall(lvl, gen);
    }
    //2. viewer
    
    char vrname[100];  mi = opts.find("-vrname"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>vrname; }

    Viewer::initGL(&argc, argv);
    BallViewer *ballviewer = new BallViewer(vrname,600,600);
    ballviewer->setObject(obj);
    glutMainLoop();
    glCheck();
    
}

