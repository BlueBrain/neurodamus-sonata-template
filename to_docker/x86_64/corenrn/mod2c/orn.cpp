/*********************************************************
Model Name      : orn
Filename        : orn.mod
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
        "orn",
        "g_e_max",
        "cc_peak",
        "g_e_baseline",
        "std_e",
        "tau_e",
        0,
        "i",
        "g_e",
        "D_e",
        0,
        "O",
        "C",
        "D",
        0,
        "donotuse",
        0
    };


    /** all global variables */
    struct orn_Store {
        int point_type{};
        double O0{};
        double C0{};
        double D0{};
        int reset{};
        int mech_type{};
        double net_receive_on{0};
        double E_e{0};
        int slist1[3]{8, 9, 10};
        int dlist1[3]{14, 15, 16};
        int slist2[3]{8, 9, 10};
        ThreadDatum ext_call_thread[3]{};
    };
    static_assert(std::is_trivially_copy_constructible_v<orn_Store>);
    static_assert(std::is_trivially_move_constructible_v<orn_Store>);
    static_assert(std::is_trivially_copy_assignable_v<orn_Store>);
    static_assert(std::is_trivially_move_assignable_v<orn_Store>);
    static_assert(std::is_trivially_destructible_v<orn_Store>);
    orn_Store orn_global;


    /** all mechanism instance variables and global variables */
    struct orn_Instance  {
        double* g_e_max{};
        double* cc_peak{};
        double* g_e_baseline{};
        double* std_e{};
        const double* tau_e{};
        double* i{};
        double* g_e{};
        double* D_e{};
        double* O{};
        double* C{};
        double* D{};
        double* g_e1{};
        double* exp_e{};
        double* amp_e{};
        double* DO{};
        double* DC{};
        double* DD{};
        double* v_unused{};
        double* g_unused{};
        double* tsave{};
        const double* node_area{};
        const int* point_process{};
        void** donotuse{};
        orn_Store* global{&orn_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"net_receive_on_orn", &orn_global.net_receive_on},
        {"E_e_orn", &orn_global.E_e},
        {nullptr, nullptr}
    };


    /** connect global (array) variables to hoc -- */
    static DoubVec hoc_vector_double[] = {
        {nullptr, nullptr, 0}
    };


    static inline int first_pointer_var_index() {
        return 2;
    }


    static inline int num_net_receive_args() {
        return 4;
    }


    /** thread specific helper routines for derivimplicit */

    static inline int* deriv1_advance(ThreadDatum* thread) {
        return &(thread[0].i);
    }

    static inline int dith1() {
        return 1;
    }

    static inline void** newtonspace1(ThreadDatum* thread) {
        return &(thread[2]._pvoid);
    }


    static inline int float_variables_size() {
        return 20;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return orn_global.mech_type;
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


    /** thread memory allocation callback */
    static void thread_mem_init(ThreadDatum* thread)  {
        thread[dith1()].pval = nullptr;
    }


    /** thread memory cleanup callback */
    static void thread_mem_cleanup(ThreadDatum* thread)  {
        free(thread[dith1()].pval);
        nrn_destroy_newtonspace(static_cast<NewtonSpace*>(*newtonspace1(thread)));
    }

    // Allocate instance structure
    static void nrn_private_constructor_orn(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new orn_Instance{};
        assert(inst->global == &orn_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(orn_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_orn(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<orn_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &orn_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(orn_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<orn_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &orn_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(orn_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->g_e_max = ml->data+0*pnodecount;
        inst->cc_peak = ml->data+1*pnodecount;
        inst->g_e_baseline = ml->data+2*pnodecount;
        inst->std_e = ml->data+3*pnodecount;
        inst->tau_e = ml->data+4*pnodecount;
        inst->i = ml->data+5*pnodecount;
        inst->g_e = ml->data+6*pnodecount;
        inst->D_e = ml->data+7*pnodecount;
        inst->O = ml->data+8*pnodecount;
        inst->C = ml->data+9*pnodecount;
        inst->D = ml->data+10*pnodecount;
        inst->g_e1 = ml->data+11*pnodecount;
        inst->exp_e = ml->data+12*pnodecount;
        inst->amp_e = ml->data+13*pnodecount;
        inst->DO = ml->data+14*pnodecount;
        inst->DC = ml->data+15*pnodecount;
        inst->DD = ml->data+16*pnodecount;
        inst->v_unused = ml->data+17*pnodecount;
        inst->g_unused = ml->data+18*pnodecount;
        inst->tsave = ml->data+19*pnodecount;
        inst->node_area = nt->_data;
        inst->point_process = ml->pdata;
        inst->donotuse = nt->_vdata;
    }



    static void nrn_alloc_orn(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_orn(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<orn_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_orn(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<orn_Instance*>(ml->instance);

        #endif
    }


    inline double normrand123_orn(int id, int pnodecount, orn_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int oup_orn(int id, int pnodecount, orn_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int initstream_orn(int id, int pnodecount, orn_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int noiseFromRandom_orn(int id, int pnodecount, orn_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int noiseFromRandom123_orn(int id, int pnodecount, orn_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
}


using namespace coreneuron;


#if !NRNBBCORE /* running in NEURON */
/*
   1 means noiseFromRandom was called when _ran_compat was previously 0 .
   2 means noiseFromRandom123 was called when _ran_compat was previously 0.
*/
static int _ran_compat; /* specifies the noise style for all instances */
#endif /* running in NEURON */


static void bbcore_write(double* x, int* d, int* xx, int *offset, int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
#if !NRNBBCORE
	/* error if using the legacy normrand */
	if (!nt->_vdata[indexes[2*pnodecount + id]]) {
		fprintf(stderr, "orn: cannot use the legacy normrand generator for the random stream.\n");
		assert(0);
	}
	if (d) {
		uint32_t* di = ((uint32_t*)d) + *offset;
		if (_ran_compat == 1) {
			Rand** pv = (Rand**)(&nt->_vdata[indexes[2*pnodecount + id]]);
			/* error if not using Random123 generator */
			if (!nrn_random_isran123(*pv, di, di+1, di+2)) {
				fprintf(stderr, "orn: Random123 generator is required\n");
				assert(0);
			}
		}else{
			nrnran123_State** pv = (nrnran123_State**)(&nt->_vdata[indexes[2*pnodecount + id]]);
			nrnran123_getids3(*pv, di, di+1, di+2);
		}
		/*printf("orn bbcore_write %d %d %d\n", di[0], di[1], di[3]);*/
	}
#endif
	*offset += 3;
}
static void bbcore_read(double* x, int* d, int* xx, int* offset, int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
	uint32_t* di = ((uint32_t*)d) + *offset;
	nrnran123_State** pv = (nrnran123_State**)(&nt->_vdata[indexes[2*pnodecount + id]]);
#if !NRNBBCORE
    if(*pv) {
        nrnran123_deletestream(*pv);
    }
#endif
	*pv = nrnran123_newstream3(di[0], di[1], di[2]);
	*offset += 3;
}


namespace coreneuron {


    inline int oup_orn(int id, int pnodecount, orn_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_oup = 0;
        if (inst->tau_e[id] != 0.0) {
            double normrand123_in_0;
            {
                  if (inst->donotuse[indexes[2*pnodecount + id]]) {
                    /*
                      :Supports separate independent but reproducible streams for
                      : each instance. However, the corresponding hoc Random
                      : distribution MUST be set to Random.negexp(1)
                    */
                    #if !NRNBBCORE
                      if(_ran_compat == 1) {
                        normrand123_in_0 = nrn_random_pick((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
                      }else{
                        normrand123_in_0 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
                      }
                    #else
                      #pragma acc routine(nrnran123_normal) seq
                      normrand123_in_0 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
                    #endif
                  }else{
                    /* only use Random123 */
                    assert(0);
                  }

            }
            inst->g_e1[id] = inst->exp_e[id] * inst->g_e1[id] + inst->amp_e[id] * normrand123_in_0;
        }
        return ret_oup;
    }


    inline int initstream_orn(int id, int pnodecount, orn_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_initstream = 0;
          if (inst->donotuse[indexes[2*pnodecount + id]]) {
            uint32_t id1, id2, id3;
            #if NRNBBCORE
              nrnran123_setseq((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]], 0, 0);
            #else
              if (_ran_compat == 1) {
                nrn_random_reset((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
              }else{
                nrnran123_setseq((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]], 0, 0);
              }
            #endif
          }	

        return ret_initstream;
    }


    inline int noiseFromRandom_orn(int id, int pnodecount, orn_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_noiseFromRandom = 0;
        #if !NRNBBCORE
         {
        	void** pv = (void**)(&inst->donotuse[indexes[2*pnodecount + id]]);
        	if (_ran_compat == 2) {
        		fprintf(stderr, "orn.noiseFromRandom123 was previously called\n");
        		assert(0);
        	}
        	_ran_compat = 1;
        	if (ifarg(1)) {
        		*pv = nrn_random_arg(1);
        	}else{
        		*pv = (void*)0;
        	}
         }
        #endif

        return ret_noiseFromRandom;
    }


    inline int noiseFromRandom123_orn(int id, int pnodecount, orn_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_noiseFromRandom123 = 0;
        #if !NRNBBCORE
         {
                nrnran123_State** pv = (nrnran123_State**)(&inst->donotuse[indexes[2*pnodecount + id]]);
                if (_ran_compat == 1) {
                  fprintf(stderr, "orn.noiseFromRandom was previously called\n");
                  assert(0);
                }
                _ran_compat = 2;
                if (*pv) {
                  nrnran123_deletestream(*pv);
                  *pv = (nrnran123_State*)0;
                }
                if (ifarg(3)) {
        	      *pv = nrnran123_newstream3((uint32_t)*getarg(1), (uint32_t)*getarg(2), (uint32_t)*getarg(3));
                }else if (ifarg(2)) {
        	      *pv = nrnran123_newstream((uint32_t)*getarg(1), (uint32_t)*getarg(2));
                }
         }
        #endif

        return ret_noiseFromRandom123;
    }


    inline double normrand123_orn(int id, int pnodecount, orn_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_normrand123 = 0.0;
          if (inst->donotuse[indexes[2*pnodecount + id]]) {
            /*
              :Supports separate independent but reproducible streams for
              : each instance. However, the corresponding hoc Random
              : distribution MUST be set to Random.negexp(1)
            */
            #if !NRNBBCORE
              if(_ran_compat == 1) {
                ret_normrand123 = nrn_random_pick((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
              }else{
                ret_normrand123 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
              }
            #else
              #pragma acc routine(nrnran123_normal) seq
              ret_normrand123 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
            #endif
          }else{
            /* only use Random123 */
            assert(0);
          }

        return ret_normrand123;
    }


    namespace {
        struct _newton_states_orn {
            int operator()(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) const {
                auto* const inst = static_cast<orn_Instance*>(ml->instance);
                double* savstate1 = static_cast<double*>(thread[dith1()].pval);
                auto const& slist1 = inst->global->slist1;
                auto const& dlist1 = inst->global->dlist1;
                double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (3*pnodecount);
                double KO, KC1, KC2, KD1, KD2;
                KO = 1.0 / 100.0;
                KC1 = 1.0 / 100.0;
                KC2 = 1e-4;
                KD1 = 1.0 / 6000.0;
                KD2 = 1.0 / 100.0;
                inst->DO[id] = KO * (1.0 - inst->C[id] - inst->O[id]);
                inst->DC[id] = KC1 * (1.0 - inst->C[id]) * inst->C[id] + KC2 * (1.0 - inst->C[id]);
                inst->DD[id] = KD1 * inst->O[id] * (1.0 - inst->D[id]) - KD2 * inst->D[id] * (1.0 - inst->O[id]);
                int counter = -1;
                for (int i=0; i<3; i++) {
                    if (*deriv1_advance(thread)) {
                        dlist2[(++counter)*pnodecount+id] = data[dlist1[i]*pnodecount+id]-(data[slist1[i]*pnodecount+id]-savstate1[i*pnodecount+id])/nt->_dt;
                    } else {
                        dlist2[(++counter)*pnodecount+id] = data[slist1[i]*pnodecount+id]-savstate1[i*pnodecount+id];
                    }
                }
                return 0;
            }
        };
    }

    int states_orn(int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
        auto* const inst = static_cast<orn_Instance*>(ml->instance);
        double* savstate1 = (double*) thread[dith1()].pval;
        auto const& slist1 = inst->global->slist1;
        auto& slist2 = inst->global->slist2;
        double* dlist2 = static_cast<double*>(thread[dith1()].pval) + (3*pnodecount);
        for (int i=0; i<3; i++) {
            savstate1[i*pnodecount+id] = data[slist1[i]*pnodecount+id];
        }
        int reset = nrn_newton_thread(static_cast<NewtonSpace*>(*newtonspace1(thread)), 3, slist2, _newton_states_orn{}, dlist2, id, pnodecount, data, indexes, thread, nt, ml, v);
        return reset;
    }




    /** initialize block for net receive */
    static void net_init(Point_process* pnt, int weight_index, double flag) {
        // do nothing
    }


    static inline void net_receive_kernel_orn(double t, Point_process* pnt, orn_Instance* inst, NrnThread* nt, Memb_list* ml, int weight_index, double flag) {
        int tid = pnt->_tid;
        int id = pnt->_i_instance;
        double v = 0;
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        double* data = ml->data;
        double* weights = nt->weights;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        double* peak = weights + weight_index + 0;
        double* base = weights + weight_index + 1;
        double* gemax = weights + weight_index + 2;
        double* stde = weights + weight_index + 3;
        inst->tsave[id] = t;
        {
            inst->C[id] = 0.0;
            if (inst->global->net_receive_on) {
                inst->cc_peak[id] = (*peak);
                inst->g_e_baseline[id] = (*base);
                inst->g_e_max[id] = (*gemax);
                inst->std_e[id] = (*stde);
            }
        }
    }


    static void net_receive_orn(Point_process* pnt, int weight_index, double flag) {
        NrnThread* nt = nrn_threads + pnt->_tid;
        Memb_list* ml = get_memb_list(nt);
        NetReceiveBuffer_t* nrb = ml->_net_receive_buffer;
        if (nrb->_cnt >= nrb->_size) {
            realloc_net_receive_buffer(nt, ml);
        }
        int id = nrb->_cnt;
        nrb->_pnt_index[id] = pnt-nt->pntprocs;
        nrb->_weight_index[id] = weight_index;
        nrb->_nrb_t[id] = nt->_t;
        nrb->_nrb_flag[id] = flag;
        nrb->_cnt++;
    }


    void net_buf_receive_orn(NrnThread* nt) {
        Memb_list* ml = get_memb_list(nt);
        if (!ml) {
            return;
        }

        NetReceiveBuffer_t* nrb = ml->_net_receive_buffer;
        auto* const inst = static_cast<orn_Instance*>(ml->instance);
        int count = nrb->_displ_cnt;
        #pragma omp simd
        #pragma ivdep
        for (int i = 0; i < count; i++) {
            int start = nrb->_displ[i];
            int end = nrb->_displ[i+1];
            for (int j = start; j < end; j++) {
                int index = nrb->_nrb_index[j];
                int offset = nrb->_pnt_index[index];
                double t = nrb->_nrb_t[index];
                int weight_index = nrb->_weight_index[index];
                double flag = nrb->_nrb_flag[index];
                Point_process* point_process = nt->pntprocs + offset;
                net_receive_kernel_orn(t, point_process, inst, nt, ml, weight_index, flag);
            }
        }
        nrb->_displ_cnt = 0;
        nrb->_cnt = 0;
    }


    /** initialize channel */
    void nrn_init_orn(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<orn_Instance*>(ml->instance);


        int& deriv_advance_flag = *deriv1_advance(thread);
        deriv_advance_flag = 0;
        auto ns = newtonspace1(thread);
        auto& th = thread[dith1()];
        if (*ns == nullptr) {
            int vec_size = 2*3*pnodecount*sizeof(double);
            double* vec = makevector(vec_size);
            th.pval = vec;
            *ns = nrn_cons_newtonspace(3, pnodecount);
        }
        if (_nrn_skip_initmodel == 0) {
            #pragma omp simd
            #pragma ivdep
            for (int id = 0; id < nodecount; id++) {
                inst->tsave[id] = -1e20;
                int node_id = node_index[id];
                double v = voltage[node_id];
                #if NRN_PRCELLSTATE
                inst->v_unused[id] = v;
                #endif
                inst->O[id] = inst->global->O0;
                inst->C[id] = inst->global->C0;
                inst->D[id] = inst->global->D0;
                {
                      if (inst->donotuse[indexes[2*pnodecount + id]]) {
                        uint32_t id1, id2, id3;
                        #if NRNBBCORE
                          nrnran123_setseq((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]], 0, 0);
                        #else
                          if (_ran_compat == 1) {
                            nrn_random_reset((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
                          }else{
                            nrnran123_setseq((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]], 0, 0);
                          }
                        #endif
                      }	

                }
                inst->g_e1[id] = 0.0;
                if (inst->tau_e[id] != 0.0) {
                    inst->D_e[id] = 2.0 * inst->std_e[id] * inst->std_e[id] / inst->tau_e[id];
                    inst->exp_e[id] = exp( -nt->_dt / inst->tau_e[id]);
                    inst->amp_e[id] = inst->std_e[id] * sqrt((1.0 - exp( -2.0 * nt->_dt / inst->tau_e[id])));
                }
                inst->O[id] = 0.0;
                inst->C[id] = 1.0;
                inst->D[id] = 0.0;
            }
        }
        deriv_advance_flag = 1;
    }


    inline double nrn_current_orn(int id, int pnodecount, orn_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        double SORN;
        SORN = inst->O[id] * (1.0 - inst->D[id]);
        inst->g_e[id] = inst->g_e1[id] + SORN * inst->cc_peak[id] * inst->g_e_max[id] + inst->g_e_baseline[id];
        if (inst->g_e[id] < 0.0) {
            inst->g_e[id] = 0.0;
        }
        inst->i[id] = inst->g_e[id] * (v - inst->global->E_e);
        current += inst->i[id];
        return current;
    }


    /** update current */
    void nrn_cur_orn(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        double* vec_rhs = nt->_actual_rhs;
        double* vec_d = nt->_actual_d;
        double* shadow_rhs = nt->_shadow_rhs;
        double* shadow_d = nt->_shadow_d;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<orn_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            double g = nrn_current_orn(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double rhs = nrn_current_orn(id, pnodecount, inst, data, indexes, thread, nt, v);
            g = (g-rhs)/0.001;
            double mfactor = 1.e2/inst->node_area[indexes[0*pnodecount + id]];
            g = g*mfactor;
            rhs = rhs*mfactor;
            #if NRN_PRCELLSTATE
            inst->g_unused[id] = g;
            #endif
            shadow_rhs[id] = rhs;
            shadow_d[id] = g;
        }
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            vec_rhs[node_id] -= shadow_rhs[id];
            vec_d[node_id] += shadow_d[id];
        }
    }


    /** update state */
    void nrn_state_orn(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<orn_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            if (inst->tau_e[id] != 0.0) {
                double normrand123_in_0;
                {
                      if (inst->donotuse[indexes[2*pnodecount + id]]) {
                        /*
                          :Supports separate independent but reproducible streams for
                          : each instance. However, the corresponding hoc Random
                          : distribution MUST be set to Random.negexp(1)
                        */
                        #if !NRNBBCORE
                          if(_ran_compat == 1) {
                            normrand123_in_0 = nrn_random_pick((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
                          }else{
                            normrand123_in_0 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
                          }
                        #else
                          #pragma acc routine(nrnran123_normal) seq
                          normrand123_in_0 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
                        #endif
                      }else{
                        /* only use Random123 */
                        assert(0);
                      }

                }
                inst->g_e1[id] = inst->exp_e[id] * inst->g_e1[id] + inst->amp_e[id] * normrand123_in_0;
            }
            states_orn(id, pnodecount, data, indexes, thread, nt, ml, v);
        }
    }


    /** register channel with the simulator */
    void _orn_reg() {

        int mech_type = nrn_get_mechtype("orn");
        orn_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        point_register_mech(mechanism, nrn_alloc_orn, nrn_cur_orn, nullptr, nrn_state_orn, nrn_init_orn, nrn_private_constructor_orn, nrn_private_destructor_orn, first_pointer_var_index(), nullptr, nullptr, 4);

        thread_mem_init(orn_global.ext_call_thread);
        _nrn_thread_reg0(mech_type, thread_mem_cleanup);
        _nrn_thread_reg1(mech_type, thread_mem_init);
        hoc_reg_bbcore_read(mech_type, bbcore_read);
        hoc_reg_bbcore_write(mech_type, bbcore_write);
        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "area");
        hoc_register_dparam_semantics(mech_type, 1, "pntproc");
        hoc_register_dparam_semantics(mech_type, 2, "bbcorepointer");
        hoc_register_net_receive_buffering(net_buf_receive_orn, mech_type);
        set_pnt_receive(mech_type, net_receive_orn, net_init, num_net_receive_args());
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
