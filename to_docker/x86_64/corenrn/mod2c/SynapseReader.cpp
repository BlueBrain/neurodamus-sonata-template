/*********************************************************
Model Name      : 
Filename        : SynapseReader.mod
NMODL Version   : 6.2.0
Vectorized      : true
Threadsafe      : true
Created         : Wed Sep 27 13:48:15 2023
Backend         : C (api-compatibility)
NMODL Compiler  : 0.6 [2ce4a2b9 2023-06-05 16:57:21 +0200]
*********************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <coreneuron/nrnconf.h>
#include <coreneuron/sim/multicore.hpp>
#include <coreneuron/mechanism/register_mech.hpp>
#include <coreneuron/gpu/nrn_acc_manager.hpp>
#include <coreneuron/utils/randoms/nrnran123.h>
#include <coreneuron/nrniv/nrniv_decl.h>
#include <coreneuron/utils/ivocvect.hpp>
#include <coreneuron/utils/nrnoc_aux.hpp>
#include <coreneuron/mechanism/mech/mod2c_core_thread.hpp>
#include <coreneuron/sim/scopmath/newton_thread.hpp>


namespace coreneuron {
    #ifndef NRN_PRCELLSTATE
    #define NRN_PRCELLSTATE 0
    #endif


    /** channel information */
    static const char *mechanism[] = {
        "6.2.0",
        "",
        0,
        0,
        0,
        0
    };


    /** all global variables */
    struct _Store {
        int reset{};
        int mech_type{};
    };
    static_assert(std::is_trivially_copy_constructible_v<_Store>);
    static_assert(std::is_trivially_move_constructible_v<_Store>);
    static_assert(std::is_trivially_copy_assignable_v<_Store>);
    static_assert(std::is_trivially_move_assignable_v<_Store>);
    static_assert(std::is_trivially_destructible_v<_Store>);
    _Store _global;


    /** all mechanism instance variables and global variables */
    struct _Instance  {
        double* v_unused{};
        _Store* global{&_global};
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
        return -1;
    }


    static inline int float_variables_size() {
        return 1;
    }


    static inline int int_variables_size() {
        return 0;
    }


    static inline int get_mech_type() {
        return _global.mech_type;
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
    static void nrn_private_constructor_(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new _Instance{};
        assert(inst->global == &_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->v_unused = ml->data+0*pnodecount;
    }



    static void nrn_alloc_(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<_Instance*>(ml->instance);

        #endif
    }


    /** initialize channel */
    void nrn_init_(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma ivdep
            #pragma omp simd
            for (int id = 0; id < nodecount; id++) {
                int node_id = node_index[id];
                double v = voltage[node_id];
                #if NRN_PRCELLSTATE
                inst->v_unused[id] = v;
                #endif
            }
        }
    }


    /** register channel with the simulator */
    void _SynapseReader_reg() {

        int mech_type = nrn_get_mechtype("");
        _global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        register_mech(mechanism, nrn_alloc_, nullptr, nullptr, nullptr, nrn_init_, nrn_private_constructor_, nrn_private_destructor_, first_pointer_var_index(), 1);

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
