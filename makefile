#------------------------------------------------------
# include machine dependent parameters
#------------------------------------------------------
include makefile.in

#------------------------------------------------------
# source files to be compiled
#------------------------------------------------------
LIB_SRC =	src/bdsurf.cpp  	src/bdulib.cpp  	src/ccsubmatlib.cpp \
			src/ccsurf.cpp  	src/ccsurfop.cpp  	src/gpmesh.cpp \
			src/vecmatop.cpp    src/quaternion.cpp  src/svdrep.cpp \
			vis/viewer3d.cpp  	vis/bdsurfobj.cpp 	vis/reflmap.cpp \
			vis/ccsurfobj.cpp  	vis/ballviewer.cpp 	vis/viewer.cpp \
			vis/psOpenGL.cpp  	vis/camera.cpp  	vis/arcball.cpp \
			vis/texture.cpp 	src/evalfromw.cpp 	src/utils.cpp

LIB_SRC_C = src/splinebasis.c
#------------------------------------------------------
# library object
#------------------------------------------------------
LIB_OBJ = 	$(LIB_SRC:.cpp=.o) $(LIB_SRC_C:.c=.o)

#------------------------------------------------------
# executable
#------------------------------------------------------
TST_SRC = 	vistt0.cpp

DEP     = 	$(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

#------------------------------------------------------
# build library archive
#------------------------------------------------------
libblend.a: $(LIB_OBJ)
	$(AR) $(ARFLAGS) $@ $(LIB_OBJ)
	$(RANLIB) $@

#------------------------------------------------------
# build executable
#------------------------------------------------------
vistt0:		vistt0.o libblend.a
	${CXX} -o vistt0 vistt0.o libblend.a ${LDFLAGS}


sampleout:		vis/vistt0_sampleout.o lib/libblend.a
	${CXX} -o sampleout vistt0_sampleout.o libblend.a ${LDFLAGS}


bdrender:		vis/bdrender.o lib/libblend.a
	${CXX} -o bdrender bdrender.o libblend.a ${LDFLAGS}


#------------------------------------------------------
# dependency information for make
#------------------------------------------------------
-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ *.d *.a *.o

tags:
	etags *hpp *cpp

dox:
	doxygen
