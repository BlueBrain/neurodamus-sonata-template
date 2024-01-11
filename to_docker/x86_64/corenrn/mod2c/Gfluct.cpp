/*********************************************************
Model Name      : Gfluct
Filename        : Gfluct.mod
NMODL Version   : 6.2.0
Vectorized      : true
Threadsafe      : true
Created         : Thu Jan 11 15:15:19 2024
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
        "Gfluct",
        "E_e",
        "E_i",
        "g_e0",
        "g_i0",
        "std_e",
        "std_i",
        "tau_e",
        "tau_i",
        0,
        "i",
        "g_e",
        "g_i",
        "g_e1",
        "g_i1",
        "D_e",
        "D_i",
        0,
        0,
        "donotuse",
        0
    };


    /** all global variables */
    struct Gfluct_Store {
        int point_type{};
        int reset{};
        int mech_type{};
    };
    static_assert(std::is_trivially_copy_constructible_v<Gfluct_Store>);
    static_assert(std::is_trivially_move_constructible_v<Gfluct_Store>);
    static_assert(std::is_trivially_copy_assignable_v<Gfluct_Store>);
    static_assert(std::is_trivially_move_assignable_v<Gfluct_Store>);
    static_assert(std::is_trivially_destructible_v<Gfluct_Store>);
    Gfluct_Store Gfluct_global;


    /** all mechanism instance variables and global variables */
    struct Gfluct_Instance  {
        const double* E_e{};
        const double* E_i{};
        const double* g_e0{};
        const double* g_i0{};
        const double* std_e{};
        const double* std_i{};
        const double* tau_e{};
        const double* tau_i{};
        double* i{};
        double* g_e{};
        double* g_i{};
        double* g_e1{};
        double* g_i1{};
        double* D_e{};
        double* D_i{};
        double* exp_e{};
        double* exp_i{};
        double* amp_e{};
        double* amp_i{};
        double* v_unused{};
        double* g_unused{};
        const double* node_area{};
        const int* point_process{};
        void** donotuse{};
        Gfluct_Store* global{&Gfluct_global};
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
        return 21;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return Gfluct_global.mech_type;
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
    static void nrn_private_constructor_Gfluct(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new Gfluct_Instance{};
        assert(inst->global == &Gfluct_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(Gfluct_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_Gfluct(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<Gfluct_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Gfluct_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Gfluct_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<Gfluct_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &Gfluct_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(Gfluct_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->E_e = ml->data+0*pnodecount;
        inst->E_i = ml->data+1*pnodecount;
        inst->g_e0 = ml->data+2*pnodecount;
        inst->g_i0 = ml->data+3*pnodecount;
        inst->std_e = ml->data+4*pnodecount;
        inst->std_i = ml->data+5*pnodecount;
        inst->tau_e = ml->data+6*pnodecount;
        inst->tau_i = ml->data+7*pnodecount;
        inst->i = ml->data+8*pnodecount;
        inst->g_e = ml->data+9*pnodecount;
        inst->g_i = ml->data+10*pnodecount;
        inst->g_e1 = ml->data+11*pnodecount;
        inst->g_i1 = ml->data+12*pnodecount;
        inst->D_e = ml->data+13*pnodecount;
        inst->D_i = ml->data+14*pnodecount;
        inst->exp_e = ml->data+15*pnodecount;
        inst->exp_i = ml->data+16*pnodecount;
        inst->amp_e = ml->data+17*pnodecount;
        inst->amp_i = ml->data+18*pnodecount;
        inst->v_unused = ml->data+19*pnodecount;
        inst->g_unused = ml->data+20*pnodecount;
        inst->node_area = nt->_data;
        inst->point_process = ml->pdata;
        inst->donotuse = nt->_vdata;
    }



    static void nrn_alloc_Gfluct(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_Gfluct(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Gfluct_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_Gfluct(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Gfluct_Instance*>(ml->instance);

        #endif
    }


    inline double normrand123_Gfluct(int id, int pnodecount, Gfluct_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int oup_Gfluct(int id, int pnodecount, Gfluct_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int noiseFromRandom_Gfluct(int id, int pnodecount, Gfluct_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
}


using namespace coreneuron;


static void bbcore_write(double* x, int* d, int* xx, int *offset, int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
#if !NRNBBCORE
	/* error if using the legacy normrand */
	if (!nt->_vdata[indexes[2*pnodecount + id]]) {
		fprintf(stderr, "Gfluct: cannot use the legacy normrand generator for the random stream.\n");
		assert(0);
	}
	if (d) {
		uint32_t* di = ((uint32_t*)d) + *offset;
		Rand** pv = (Rand**)(&nt->_vdata[indexes[2*pnodecount + id]]);
		/* error if not using Random123 generator */
		if (!nrn_random_isran123(*pv, di, di+1, di+2)) {
			fprintf(stderr, "Gfluct: Random123 generator is required\n");
			assert(0);
		}
		/*printf("Gfluct bbcore_write %d %d %d\n", di[0], di[1], di[3]);*/
	}
	*offset += 3;
#endif
}
static void bbcore_read(double* x, int* d, int* xx, int* offset, int id, int pnodecount, double* data, Datum* indexes, ThreadDatum* thread, NrnThread* nt, Memb_list* ml, double v) {
	assert(!nt->_vdata[indexes[2*pnodecount + id]]);
	uint32_t* di = ((uint32_t*)d) + *offset;
	nrnran123_State** pv = (nrnran123_State**)(&nt->_vdata[indexes[2*pnodecount + id]]);
	*pv = nrnran123_newstream3(di[0], di[1], di[2]);
	*offset += 3;
}


namespace coreneuron {


    inline int oup_Gfluct(int id, int pnodecount, Gfluct_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_oup = 0;
        if (inst->tau_e[id] != 0.0) {
            double normrand123_in_2;
            {
                	if (inst->donotuse[indexes[2*pnodecount + id]]) {
                		/*
                		:Supports separate independent but reproducible streams for
                		: each instance. However, the corresponding hoc Random
                		: distribution MUST be set to Random.negexp(1)
                		*/
                        #if !NRNBBCORE
                		normrand123_in_2 = nrn_random_pick((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
                        #else
                        #pragma acc routine(nrnran123_normal) seq
                        normrand123_in_2 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
                        #endif
                	}else{
                		/* only use Random123 */
                        assert(0);
                	}

            }
            inst->g_e1[id] = inst->exp_e[id] * inst->g_e1[id] + inst->amp_e[id] * normrand123_in_2;
        }
        if (inst->tau_i[id] != 0.0) {
            double normrand123_in_3;
            {
                	if (inst->donotuse[indexes[2*pnodecount + id]]) {
                		/*
                		:Supports separate independent but reproducible streams for
                		: each instance. However, the corresponding hoc Random
                		: distribution MUST be set to Random.negexp(1)
                		*/
                        #if !NRNBBCORE
                		normrand123_in_3 = nrn_random_pick((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
                        #else
                        #pragma acc routine(nrnran123_normal) seq
                        normrand123_in_3 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
                        #endif
                	}else{
                		/* only use Random123 */
                        assert(0);
                	}

            }
            inst->g_i1[id] = inst->exp_i[id] * inst->g_i1[id] + inst->amp_i[id] * normrand123_in_3;
        }
        return ret_oup;
    }


    inline int noiseFromRandom_Gfluct(int id, int pnodecount, Gfluct_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_noiseFromRandom = 0;
        #if !NRNBBCORE
         {
        	void** pv = (void**)(&inst->donotuse[indexes[2*pnodecount + id]]);
        	if (ifarg(1)) {
        		*pv = nrn_random_arg(1);
        	}else{
        		*pv = (void*)0;
        	}
         }
        #endif

        return ret_noiseFromRandom;
    }


    inline double normrand123_Gfluct(int id, int pnodecount, Gfluct_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_normrand123 = 0.0;
        	if (inst->donotuse[indexes[2*pnodecount + id]]) {
        		/*
        		:Supports separate independent but reproducible streams for
        		: each instance. However, the corresponding hoc Random
        		: distribution MUST be set to Random.negexp(1)
        		*/
                #if !NRNBBCORE
        		ret_normrand123 = nrn_random_pick((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
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


    /** initialize channel */
    void nrn_init_Gfluct(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<Gfluct_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma omp simd
            #pragma ivdep
            for (int id = 0; id < nodecount; id++) {
                int node_id = node_index[id];
                double v = voltage[node_id];
                #if NRN_PRCELLSTATE
                inst->v_unused[id] = v;
                #endif
                inst->g_e1[id] = 0.0;
                inst->g_i1[id] = 0.0;
                if (inst->tau_e[id] != 0.0) {
                    inst->D_e[id] = 2.0 * inst->std_e[id] * inst->std_e[id] / inst->tau_e[id];
                    inst->exp_e[id] = exp( -nt->_dt / inst->tau_e[id]);
                    inst->amp_e[id] = inst->std_e[id] * sqrt((1.0 - exp( -2.0 * nt->_dt / inst->tau_e[id])));
                }
                if (inst->tau_i[id] != 0.0) {
                    inst->D_i[id] = 2.0 * inst->std_i[id] * inst->std_i[id] / inst->tau_i[id];
                    inst->exp_i[id] = exp( -nt->_dt / inst->tau_i[id]);
                    inst->amp_i[id] = inst->std_i[id] * sqrt((1.0 - exp( -2.0 * nt->_dt / inst->tau_i[id])));
                }
            }
        }
    }


    inline double nrn_current_Gfluct(int id, int pnodecount, Gfluct_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        if (inst->tau_e[id] == 0.0) {
            double normrand123_in_0;
            {
                	if (inst->donotuse[indexes[2*pnodecount + id]]) {
                		/*
                		:Supports separate independent but reproducible streams for
                		: each instance. However, the corresponding hoc Random
                		: distribution MUST be set to Random.negexp(1)
                		*/
                        #if !NRNBBCORE
                		normrand123_in_0 = nrn_random_pick((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
                        #else
                        #pragma acc routine(nrnran123_normal) seq
                        normrand123_in_0 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
                        #endif
                	}else{
                		/* only use Random123 */
                        assert(0);
                	}

            }
            inst->g_e[id] = inst->std_e[id] * normrand123_in_0;
        }
        if (inst->tau_i[id] == 0.0) {
            double normrand123_in_1;
            {
                	if (inst->donotuse[indexes[2*pnodecount + id]]) {
                		/*
                		:Supports separate independent but reproducible streams for
                		: each instance. However, the corresponding hoc Random
                		: distribution MUST be set to Random.negexp(1)
                		*/
                        #if !NRNBBCORE
                		normrand123_in_1 = nrn_random_pick((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
                        #else
                        #pragma acc routine(nrnran123_normal) seq
                        normrand123_in_1 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
                        #endif
                	}else{
                		/* only use Random123 */
                        assert(0);
                	}

            }
            inst->g_i[id] = inst->std_i[id] * normrand123_in_1;
        }
        inst->g_e[id] = inst->g_e0[id] + inst->g_e1[id];
        if (inst->g_e[id] < 0.0) {
            inst->g_e[id] = 0.0;
        }
        inst->g_i[id] = inst->g_i0[id] + inst->g_i1[id];
        if (inst->g_i[id] < 0.0) {
            inst->g_i[id] = 0.0;
        }
        inst->i[id] = inst->g_e[id] * (v - inst->E_e[id]) + inst->g_i[id] * (v - inst->E_i[id]);
        current += inst->i[id];
        return current;
    }


    /** update current */
    void nrn_cur_Gfluct(NrnThread* nt, Memb_list* ml, int type) {
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
        auto* const inst = static_cast<Gfluct_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            double g = nrn_current_Gfluct(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double rhs = nrn_current_Gfluct(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_Gfluct(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<Gfluct_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            if (inst->tau_e[id] != 0.0) {
                double normrand123_in_2;
                {
                    	if (inst->donotuse[indexes[2*pnodecount + id]]) {
                    		/*
                    		:Supports separate independent but reproducible streams for
                    		: each instance. However, the corresponding hoc Random
                    		: distribution MUST be set to Random.negexp(1)
                    		*/
                            #if !NRNBBCORE
                    		normrand123_in_2 = nrn_random_pick((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
                            #else
                            #pragma acc routine(nrnran123_normal) seq
                            normrand123_in_2 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
                            #endif
                    	}else{
                    		/* only use Random123 */
                            assert(0);
                    	}

                }
                inst->g_e1[id] = inst->exp_e[id] * inst->g_e1[id] + inst->amp_e[id] * normrand123_in_2;
            }
            if (inst->tau_i[id] != 0.0) {
                double normrand123_in_3;
                {
                    	if (inst->donotuse[indexes[2*pnodecount + id]]) {
                    		/*
                    		:Supports separate independent but reproducible streams for
                    		: each instance. However, the corresponding hoc Random
                    		: distribution MUST be set to Random.negexp(1)
                    		*/
                            #if !NRNBBCORE
                    		normrand123_in_3 = nrn_random_pick((Rand*)inst->donotuse[indexes[2*pnodecount + id]]);
                            #else
                            #pragma acc routine(nrnran123_normal) seq
                            normrand123_in_3 = nrnran123_normal((nrnran123_State*)inst->donotuse[indexes[2*pnodecount + id]]);
                            #endif
                    	}else{
                    		/* only use Random123 */
                            assert(0);
                    	}

                }
                inst->g_i1[id] = inst->exp_i[id] * inst->g_i1[id] + inst->amp_i[id] * normrand123_in_3;
            }
        }
    }


    /** register channel with the simulator */
    void _Gfluct_reg() {

        int mech_type = nrn_get_mechtype("Gfluct");
        Gfluct_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        point_register_mech(mechanism, nrn_alloc_Gfluct, nrn_cur_Gfluct, nullptr, nrn_state_Gfluct, nrn_init_Gfluct, nrn_private_constructor_Gfluct, nrn_private_destructor_Gfluct, first_pointer_var_index(), nullptr, nullptr, 1);

        hoc_reg_bbcore_read(mech_type, bbcore_read);
        hoc_reg_bbcore_write(mech_type, bbcore_write);
        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "area");
        hoc_register_dparam_semantics(mech_type, 1, "pntproc");
        hoc_register_dparam_semantics(mech_type, 2, "bbcorepointer");
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
