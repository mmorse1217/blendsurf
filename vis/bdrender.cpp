#include "ccsubmatlib.hpp"
#include "ccsurf.hpp"
#include "bdsurf.hpp"

#include "camera.hpp"
#include "ballviewer.hpp"
#include "ccsurfobj.hpp"
#include "bdsurfobj.hpp"
#include "evalfromw.hpp"

int optionsCreate(int argc,  char* argv[], map<string,string>& options)
{
    options.clear();
    ifstream fin(argv[1]); assert(fin.good());
    string name;  fin>>name;
    while(fin.good()) {
        char cont[100];   fin.getline(cont, 99);
        options[name] = string(cont);
        fin>>name;
    }    
    fin.close();
    int arg_cnt = 2;
    while (arg_cnt < argc) { 
      options[string(argv[arg_cnt])] = string(argv[arg_cnt+1]);
      arg_cnt += 2;
    }
    return 0;
}

int WriteAndQuit = 0;

// ---------------------------------------------------------------------- 
int main(int argc, char** argv)
{ 
    assert( argc %2 == 0 );

    map<string,string> opts;
    optionsCreate(argc, argv, opts);

    map<string,string>::iterator mi;

    int objtype;
    mi = opts.find("-objtype"); assert(mi!=opts.end());  { istringstream ss((*mi).second); ss>>objtype; }
    GeoObject* obj = NULL;
    Camera* camera = 0; 

    int renderctrl=1;
    mi = opts.find("-rendertype"); if(mi!=opts.end()) { istringstream ss((*mi).second); ss>>renderctrl; }
    int surfctrl = 0;
    mi = opts.find("-surfrendertype"); if(mi!=opts.end()) { istringstream ss((*mi).second); ss>>surfctrl; }

        char fileb[100];	 mi = opts.find("-meshfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>fileb; }
        ifstream tmpb(fileb);
    

    if(objtype==0) { //dealing with catmull-clark surface
        //------------------------

        //submatlib stores a subdivision matrix library for a two ring around a vertex with any valance.
        char filea[100];	 mi = opts.find("-ccsurf_submatlibfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filea; }
        ifstream tmpa(filea);
        CCSubMatLib* submatlib = new CCSubMatLib();
        cerr<<"submatlib setup"<<endl;
        submatlib->setup(tmpa);
    
        //setup the surface(read the mesh, share, etc)

        CCSurf* ccsurf = new CCSurf( );
        double flats;  mi = opts.find("-ccsurf_flats"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>flats;} 
        cerr<<"ccsurf setup"<<endl;
        ccsurf->setup(flats, tmpb,0,Point3(1.0), (void*)submatlib);
  
        int lvl;  
        mi = opts.find("-ccsurf_renderlvl"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>lvl; }

        int tg = lvl;
        if( ccsurf->validpos(tg)==false ) {
            //level of subdivion not enough for renderer.need tos subdivide more
	  cerr <<"valid pos is false"<<endl;
            int lv = tg;
            while(ccsurf->validpos(lv)==false && lv>0)
                lv--; //find the nearest valid one
            assert(ccsurf->validpos(lv)==true);
      
            for(int i=lv; i<tg; i++) 
                ccsurf->subdivide(i, *submatlib); //subdivide if necessary      
        }
        cerr<<"ccsurfobj setup begin"<<endl;
        //surfobj takes care of rendering.
        obj = new CCSurfObj(ccsurf, lvl,renderctrl,surfctrl);
        cerr<<"ccsurfobj setup end"<<endl;
    } 

    else if(objtype==1) { //the blended surface
        //------------------------
        char filea[100];	 mi = opts.find("-bdsurf_submatlibfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filea; }
        ifstream tmpa(filea);
        CCSubMatLib* submatlib = new CCSubMatLib();
        cerr<<"submatlib setup"<<endl;
        submatlib->setup(tmpa);
        char filebd[100];	 mi = opts.find("-bdsurf_bdulibfile"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>filebd; }
        ifstream tmpbd(filebd);
        BdULib* bdulib = new BdULib();
        cerr<<"bdulib setup"<<endl;
        bdulib->setup(tmpbd);
    
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

        extern int NUM_PATCHES;	    mi = opts.find("-bdsurf_g");   if(mi!=opts.end()) { istringstream ss((*mi).second); ss>> NUM_PATCHES; }


	mi = opts.find("-bdsurf_matdir"); if(mi!=opts.end()) { istringstream ss((*mi).second); ss>> MatricesDir; }

    int activevert = -1;  
    mi = opts.find("-bdsurf_activevert"); if(mi!=opts.end()) { istringstream ss((*mi).second); ss>>activevert; }

        BdSurf* bdsurf = new BdSurf(  );
        bdsurf->setParams(cl,pc,ct,bt,pp,lb,ub, gl, fs, cot, sd ,poudeg, loadW );
        cerr<<"bdsurf setup"<<endl;
        bdsurf->setup(tmpb,0,Point3(1.0),(void*)submatlib,(void*)bdulib);
        //ofstream fout("bddump.dat"); bdsurf->dump(fout); fout.close();
    
        int lvl; mi = opts.find("-bdsurf_renderlvl");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>lvl; }
        int gen; mi = opts.find("-bdsurf_rendergen");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>gen; }
        int alt; mi = opts.find("-bdsurf_alt");  assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>alt; }


	cerr<<"bdsurfobj setup"<<endl;
        char filetex[100];  mi = opts.find("-bdsurfobj_texfile"); assert(mi!=opts.end());   { istringstream ss((*mi).second); ss>>filetex; }
	obj = new BdSurfObj(bdsurf, lvl, gen, alt, filetex,renderctrl,surfctrl,activevert);
	   //bdsurf->evalall(lvl, gen);
    }
    //2. viewer
    
    char vrname[100];  mi = opts.find("-vrname"); assert(mi!=opts.end()); { istringstream ss((*mi).second); ss>>vrname; }


    mi = opts.find("-writeandquit");  if(mi!=opts.end())  
	   { istringstream ss((*mi).second); ss>>WriteAndQuit; }

    mi = opts.find("-camfile"); if(mi!=opts.end()) {  
      string camfname;
      istringstream ss((*mi).second); ss>>camfname;	  
      cerr << "loading camera from "  << camfname.c_str() << endl;
      camera = new Camera(); 
      ifstream camif( camfname.c_str() );
      assert(camif.good());
      camera->read(camif); 
      camif.close();
    }

    Viewer::initGL(&argc, argv);
    BallViewer *ballviewer = new BallViewer(vrname,600,600);
    ballviewer->setObject(obj);
    if (camera) ballviewer->setCamera(camera);
    glutMainLoop();
    glCheck();
    
}

