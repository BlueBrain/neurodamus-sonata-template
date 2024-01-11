/*********************************************************
Model Name      : lookupTableV2
Filename        : lookupTableV2.mod
NMODL Version   : 6.2.0
Vectorized      : true
Threadsafe      : true
Created         : Thu Jan 11 15:15:20 2024
Backend         : C++ (api-compatibility)
NMODL Compiler  : 0.6 [6f6db3b6 2023-10-20 15:10:33 +0200]
*********************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <coreneuron/gpu/nrn_acc_manager.hpp>
#include <coreneuron/mechanism/mech/mod2c_core_thread.hpp>
#include <coreneuron/mechanism/register_mech.hpp>
#include <coreneuron/nrnconf.h>
#include <coreneuron/nrniv/nrniv_decl.h>
#include <coreneuron/sim/multicore.hpp>
#include <coreneuron/sim/scopmath/newton_thread.hpp>
#include <coreneuron/utils/ivocvect.hpp>
#include <coreneuron/utils/nrnoc_aux.hpp>
#include <coreneuron/utils/randoms/nrnran123.h>


namespace coreneuron {
    #ifndef NRN_PRCELLSTATE
    #define NRN_PRCELLSTATE 0
    #endif


    /** channel information */
    static const char *mechanism[] = {
        "6.2.0",
        "lookupTableV2",
        0,
        0,
        0,
        "ptr",
        0
    };


    /** all global variables */
    struct lookupTableV2_Store {
        int point_type{};
        int reset{};
        int mech_type{};
    };
    static_assert(std::is_trivially_copy_constructible_v<lookupTableV2_Store>);
    static_assert(std::is_trivially_move_constructible_v<lookupTableV2_Store>);
    static_assert(std::is_trivially_copy_assignable_v<lookupTableV2_Store>);
    static_assert(std::is_trivially_move_assignable_v<lookupTableV2_Store>);
    static_assert(std::is_trivially_destructible_v<lookupTableV2_Store>);
    lookupTableV2_Store lookupTableV2_global;


    /** all mechanism instance variables and global variables */
    struct lookupTableV2_Instance  {
        double* v_unused{};
        const double* node_area{};
        void** point_process{};
        double* ptr{};
        lookupTableV2_Store* global{&lookupTableV2_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {nullptr, nullptr}
    };


    /** connect global (array) variables to hoc -- */
    static DoubVec hoc_vector_double[] = {
        {nullptr, nullptr, 0}
    };


    static inline int first_pointer_var_index() {
        return 2;
    }


    static inline int float_variables_size() {
        return 1;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return lookupTableV2_global.mech_type;
    }


    static inline Memb_list* get_memb_list(NrnThread* nt) {
        if (!nt->_ml_list) {
            return nullptr;
        }
        return nt->_ml_list[get_mech_type()];
    }


    static inline void* mem_alloc(size_t num, size_t size, size_t alignment = 16) {
        void* ptr;
        posix_memalign(&ptr, alignment, num*size);
        memset(ptr, 0, size);
        return ptr;
    }


    static inline void mem_free(void* ptr) {
        free(ptr);
    }


    static inline void coreneuron_abort() {
        abort();
    }

    // Allocate instance structure
    static void nrn_private_constructor_lookupTableV2(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new lookupTableV2_Instance{};
        assert(inst->global == &lookupTableV2_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(lookupTableV2_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_lookupTableV2(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<lookupTableV2_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &lookupTableV2_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(lookupTableV2_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<lookupTableV2_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &lookupTableV2_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(lookupTableV2_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->v_unused = ml->data+0*pnodecount;
        inst->node_area = nt->_data;
        inst->point_process = nt->_vdata;
        inst->ptr = nt->_data;
    }



    static void nrn_alloc_lookupTableV2(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_lookupTableV2(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<lookupTableV2_Instance*>(ml->instance);

        #ifndef CORENEURON_BUILD
        	
        	if(sizeof(float)!=4 || sizeof(int)!=4){
        		printf("sizeof does not match. Need to specify data type sizes more explicitly\n");
        		return;
        	}
        	tableStruct** tempTable = (tableStruct**)(&(inst->ptr[indexes[2*pnodecount + id]]));
        	tableStruct* tbl = 0;
        	tbl = (tableStruct*)hoc_Emalloc(sizeof(tableStruct));
        	tbl->numToRead = 0;
        	if(ifarg(1)&&hoc_is_str_arg(1)){
        		char tableName[128];
        		sprintf(tableName,"%s",gargstr(1));
        		FILE *file;
        		
        		if((file = fopen(tableName, "r"))==NULL) {
        			
        			
        			*tempTable = tbl;   			 
        			return;
        		}
        		float retVal = 0;
        		
        		int i,j,k;
        		int gidNum;
        		char header[128] = {0x0};
        		char readChar = 0x1;
        		i=0;
        		int iAmAFloat;
        		while(readChar!=0xA){
        			fread(&readChar,1,1,file);
        			header[i++]=readChar;
        		}
        		if(strncmp(header,"ExtracellularElectrodeLookupTable",33)!=0){
        			*tempTable = tbl;
        			printf("Header does not match: \n");
        			printf("%s", header);
        			return;
        		}
        		fread(&(tbl->vInfo_lookupTableV2),sizeof(char),1,file);		
        		char circuitPath[512];
        		readChar = 0x1;
        		i=0;		
        		while(readChar!=0xA){
        			fread(&readChar,1,1,file);
        			circuitPath[i++]=readChar;
        		}
        		
        		fread(&(tbl->numToRead),sizeof(int),1,file);        
        		tbl->numToRead = (int)htonl(tbl->numToRead);
        		
        		tbl->gids = (int**)hoc_Emalloc((tbl->numToRead)*sizeof(int *));
        		tbl->factors = (float**)hoc_Emalloc((tbl->numToRead)*sizeof(float *));
        		tbl->secData = (sectionEntry**)hoc_Emalloc((tbl->numToRead)*sizeof(sectionEntry *));
        		tbl->tableSizes = (int*)hoc_Emalloc((tbl->numToRead)*sizeof(int));
        		tbl->xPos = (float*)hoc_Emalloc((tbl->numToRead)*sizeof(float));
        		tbl->yPos = (float*)hoc_Emalloc((tbl->numToRead)*sizeof(float));
        		tbl->zPos = (float*)hoc_Emalloc((tbl->numToRead)*sizeof(float));
        		
        		for(i = 0; i < tbl->numToRead; ++i){				
        			fread((&iAmAFloat),sizeof(float),1,file);
        			tbl->xPos[i] = swapForFloat(iAmAFloat);
        			fread((&iAmAFloat),sizeof(float),1,file);
        			tbl->yPos[i] = swapForFloat(iAmAFloat);
        			fread((&iAmAFloat),sizeof(float),1,file);								
        			tbl->zPos[i] = swapForFloat(iAmAFloat);				
        			int tableSize;
        			fread((&tableSize),sizeof(int),1,file);
        			tableSize = (int)htonl(tableSize);
        			
        			(tbl->tableSizes)[i]=tableSize;
        			(tbl->gids)[i] = (int*)hoc_Emalloc(tableSize*sizeof(int));
        			(tbl->factors)[i] = (float*)hoc_Emalloc(tableSize*sizeof(float));
        			(tbl->secData)[i] = (sectionEntry*)hoc_Emalloc(tableSize*sizeof(sectionEntry));
        			if((tbl->gids)[i]==0 || (tbl->factors)[i]==0 || (tbl->secData)[i]==0){
        				printf("Problem allocating memory for factor tables\n");
        				return;
        			}						
        			int index;				
        			float factor;
        			
        			
        			unsigned int somaSize, axonSize, dendriteSize, apicalSize;
        			for (j = 0; j < tableSize; ++j){
        				fread(&index,4,1,file);
        				(tbl->gids)[i][j]=(int)htonl(index);
        				
        				fread((&iAmAFloat),sizeof(iAmAFloat),1,file);	
        				
        				
        				
        				(tbl->factors)[i][j]=swapForFloat(iAmAFloat);
        				fread(&somaSize,4,1,file);
        				somaSize = (unsigned int)htonl(somaSize);
        				fread(&axonSize,4,1,file);
        				axonSize = (unsigned int)htonl(axonSize);
        				(tbl->secData)[i][j].axonSz = axonSize;
        				fread(&dendriteSize,4,1,file);
        				dendriteSize = (unsigned int)htonl(dendriteSize);
        				(tbl->secData)[i][j].basalSz = dendriteSize;
        				fread(&apicalSize,4,1,file);
        				apicalSize = (unsigned int)htonl(apicalSize);
        				(tbl->secData)[i][j].apicalSz = apicalSize;
        				if(somaSize!=1){
        					printf("Need exactly one value for the soma. Got %d",somaSize);
        					return;
        				}
        				
        				(tbl->secData)[i][j].axonsSz = (int*)hoc_Emalloc(axonSize*sizeof(int));
        				(tbl->secData)[i][j].axon = (float**)hoc_Emalloc(axonSize*sizeof(float*));				
        				(tbl->secData)[i][j].basalsSz = (int*)hoc_Emalloc(dendriteSize*sizeof(int));
        				(tbl->secData)[i][j].basals = (float**)hoc_Emalloc(dendriteSize*sizeof(float*));
        				(tbl->secData)[i][j].apicalsSz = (int*)hoc_Emalloc(apicalSize*sizeof(int));
        				(tbl->secData)[i][j].apicals = (float**)hoc_Emalloc(apicalSize*sizeof(float*));
        				
        				fread(&somaSize,4,1,file);
        				somaSize = (unsigned int)htonl(somaSize);
        				if(somaSize!=1){
        					printf("Need exactly one value for the soma. Got %d",somaSize);
        					return;
        				}				
        				readAndSwapInt((tbl->secData)[i][j].axonsSz,axonSize,file);
        				readAndSwapInt((tbl->secData)[i][j].basalsSz,dendriteSize,file);
        				readAndSwapInt((tbl->secData)[i][j].apicalsSz,apicalSize,file);
        				
        				fread((&iAmAFloat),sizeof(iAmAFloat),1,file);
        				(tbl->secData)[i][j].soma = swapForFloat(iAmAFloat);		
        				for(k = 0; k < axonSize; ++k){
        					(tbl->secData)[i][j].axon[k] = (float*)hoc_Emalloc((tbl->secData)[i][j].axonsSz[k]*sizeof(float));
        					readAndSwapFloat((tbl->secData)[i][j].axon[k],(tbl->secData)[i][j].axonsSz[k],file);					
        				}
        				for(k = 0; k < dendriteSize; ++k){
        					(tbl->secData)[i][j].basals[k] = (float*)hoc_Emalloc((tbl->secData)[i][j].basalsSz[k]*sizeof(float));
        					readAndSwapFloat((tbl->secData)[i][j].basals[k],(tbl->secData)[i][j].basalsSz[k],file);					
        				}
        				for(k = 0; k < apicalSize; ++k){
        					(tbl->secData)[i][j].apicals[k] = (float*)hoc_Emalloc((tbl->secData)[i][j].apicalsSz[k]*sizeof(float));
        					readAndSwapFloat((tbl->secData)[i][j].apicals[k],(tbl->secData)[i][j].apicalsSz[k],file);					
        				}								
        			}
        		}   		
        	}
        	*tempTable = tbl;
        #endif  // CORENEURON_BUILD
        	

        #endif
    }


    void nrn_destructor_lookupTableV2(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<lookupTableV2_Instance*>(ml->instance);

        #endif
    }


    inline double vInfo_lookupTableV2(int id, int pnodecount, lookupTableV2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double getValueForGid_lookupTableV2(int id, int pnodecount, lookupTableV2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double getValueForSection_lookupTableV2(int id, int pnodecount, lookupTableV2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
}


using namespace coreneuron;


#ifndef CORENEURON_BUILD
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netinet/in.h>
typedef struct{
	float **axon;
	float **apicals;
	float **basals;
	int *axonsSz;
	int *apicalsSz;
	int *basalsSz;
	float soma;
	unsigned int axonSz;
	unsigned int apicalSz;
	unsigned int basalSz;
} sectionEntry;
typedef struct{
	int **gids;
	sectionEntry **secData;
	float **factors;
	int numToRead;
	int *tableSizes;
	char vInfo_lookupTableV2;
	float *xPos;
	float *yPos;
	float *zPos;
}tableStruct;
float swapForFloat(int iAmAFloat){
	float factor;
	unsigned char *inPoint = (unsigned char *) &iAmAFloat;
	unsigned char *outPoint = (unsigned char *) &factor;
	iAmAFloat = (int)htonl(iAmAFloat);
	int k;
	for(k = 0; k < 4; ++k)
		outPoint[k]=inPoint[k];
	return factor;
};
void readAndSwapInt(int *readInto, int number, FILE *file){
	int i = 0;
	fread(readInto,sizeof(int),number,file);
	for(i = 0; i < number; ++i){
		readInto[i] = (int)htonl(readInto[i]);
	}        	
};
void readAndSwapFloat(float *readInto, int number, FILE *file){
	int i = 0;
	int tempReadArray[number];
	fread(&tempReadArray,sizeof(int),number,file);
	for(i = 0; i < number; ++i){
		readInto[i] = swapForFloat(tempReadArray[i]);
	}        	
};
#endif  // CORENEURON_BUILD


namespace coreneuron {


    inline double vInfo_lookupTableV2(int id, int pnodecount, lookupTableV2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_vInfo = 0.0;
        #ifndef CORENEURON_BUILD
        	tableStruct **tempData = (tableStruct**)(&inst->ptr[indexes[2*pnodecount + id]]);
        	tableStruct *tbl = (tableStruct*) *tempData;
        	return tbl->ret_vInfo;
        #endif  // CORENEURON_BUILD
        	

        return ret_vInfo;
    }


    inline double getValueForGid_lookupTableV2(int id, int pnodecount, lookupTableV2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_getValueForGid = 0.0;
        #ifndef CORENEURON_BUILD
        	tableStruct **tempData = (tableStruct**)(&inst->ptr[indexes[2*pnodecount + id]]);
        	tableStruct *tbl = (tableStruct*) *tempData;
        	if((tbl->numToRead)==0)
        		return 1;
        	float retVal = 0;
        	int i,j;
        	if(ifarg(1)){
        		int targetGid = (int) *getarg(1);
        		for(i = 0; i < (tbl->numToRead); ++i){
        			for (j = 0; j < (tbl->tableSizes)[i]; ++j){
        				if((tbl->gids)[i][j]==targetGid){
        					retVal+=(tbl->factors)[i][j];
        					break; 
        				}
        			}
        		}
        	}
        	return retVal;   
        #endif  // CORENEURON_BUILD
        	

        return ret_getValueForGid;
    }


    inline double getValueForSection_lookupTableV2(int id, int pnodecount, lookupTableV2_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_getValueForSection = 0.0;
        #ifndef CORENEURON_BUILD
        	tableStruct **tempData = (tableStruct**)(&inst->ptr[indexes[2*pnodecount + id]]);
        	tableStruct *tbl = (tableStruct*) *tempData;
        	if((tbl->numToRead)==0)
        		return 1;
        	float retVal = 0;
        	int i,j;
        	if(ifarg(4)){
        		int targetGid = (int) *getarg(1);
        		int targetType = (int) *getarg(2);
        		int targetSec = (int) *getarg(3);
        		int usedCompIndex;
        		float compDist = (float) *getarg(4);
        		for(i = 0; i < (tbl->numToRead); ++i){	
        			for (j = 0; j < (tbl->tableSizes)[i]; ++j){
        				if((tbl->gids)[i][j]==targetGid){ 
        					sectionEntry tEntry = (tbl->secData)[i][j];
        					if(targetType == 3)		
        						if(targetSec<tEntry.basalSz){
        							usedCompIndex = (int)(tEntry.basalsSz[targetSec]*compDist);
        							retVal+=tEntry.basals[targetSec][usedCompIndex];	
        						} else
        							printf("Target Section out of bounds: %i, %i, %i",targetGid,targetType,targetSec);
        					else if(targetType == 4)
        						if(targetSec<tEntry.apicalSz){
        							usedCompIndex = (int)(tEntry.apicalsSz[targetSec]*compDist);
        							retVal+=tEntry.apicals[targetSec][usedCompIndex];
        						} else
        							printf("Target Section out of bounds!");
        					else if(targetType == 2)
        						if(targetSec<tEntry.axonSz){
        							usedCompIndex = (int)(tEntry.axonsSz[targetSec]*compDist);
        							retVal+=tEntry.axon[targetSec][usedCompIndex];
        						} else
        							printf("Target Section out of bounds!");
        					else if(targetType == 1)
        						retVal+=tEntry.soma;
        					break; 
        				}
        			}
        		}
        	}
        	return retVal;
        #endif  // CORENEURON_BUILD
        	

        return ret_getValueForSection;
    }


    /** initialize channel */
    void nrn_init_lookupTableV2(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<lookupTableV2_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma omp simd
            #pragma ivdep
            for (int id = 0; id < nodecount; id++) {
                double v = 0.0;
            }
        }
    }


    /** register channel with the simulator */
    void _lookupTableV2_reg() {

        int mech_type = nrn_get_mechtype("lookupTableV2");
        lookupTableV2_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        point_register_mech(mechanism, nrn_alloc_lookupTableV2, nullptr, nullptr, nullptr, nrn_init_lookupTableV2, nrn_private_constructor_lookupTableV2, nrn_private_destructor_lookupTableV2, first_pointer_var_index(), nrn_constructor_lookupTableV2, nullptr, 1);

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "area");
        hoc_register_dparam_semantics(mech_type, 1, "pntproc");
        hoc_register_dparam_semantics(mech_type, 2, "pointer");
        add_nrn_artcell(mech_type, 0);
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
