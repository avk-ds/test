anal_basic_c217_newformat:
	g++ -g -pthread -m64 -Wno-deprecated -I${ROOTSYS}/include -o anal_rpcdata_c217_v0 EveTree.C StraightLineFit.cc anal_rpcdata_c217_v0.cc `root-config --cflags --libs` -lMinuit

clean:
	$(RM) *.o *~ 
