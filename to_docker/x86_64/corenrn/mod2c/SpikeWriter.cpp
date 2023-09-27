/*********************************************************
Model Name      : SpikeWriter
Filename        : SpikeWriter.mod
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
        "SpikeWriter",
        0,
        0,
        0,
        0
    };


    /** all global variables */
    struct SpikeWriter_Store {
        int point_type{};
        int reset{};
        int mech_type{};
    };
    static_assert(std::is_trivially_copy_constructible_v<SpikeWriter_Store>);
    static_assert(std::is_trivially_move_constructible_v<SpikeWriter_Store>);
    static_assert(std::is_trivially_copy_assignable_v<SpikeWriter_Store>);
    static_assert(std::is_trivially_move_assignable_v<SpikeWriter_Store>);
    static_assert(std::is_trivially_destructible_v<SpikeWriter_Store>);
    SpikeWriter_Store SpikeWriter_global;


    /** all mechanism instance variables and global variables */
    struct SpikeWriter_Instance  {
        double* v_unused{};
        const double* node_area{};
        void** point_process{};
        SpikeWriter_Store* global{&SpikeWriter_global};
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
        return 2;
    }


    static inline int get_mech_type() {
        return SpikeWriter_global.mech_type;
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
    static void nrn_private_constructor_SpikeWriter(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new SpikeWriter_Instance{};
        assert(inst->global == &SpikeWriter_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(SpikeWriter_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_SpikeWriter(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<SpikeWriter_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &SpikeWriter_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(SpikeWriter_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<SpikeWriter_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &SpikeWriter_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(SpikeWriter_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->v_unused = ml->data+0*pnodecount;
        inst->node_area = nt->_data;
        inst->point_process = nt->_vdata;
    }



    static void nrn_alloc_SpikeWriter(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_SpikeWriter(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<SpikeWriter_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_SpikeWriter(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<SpikeWriter_Instance*>(ml->instance);

        #endif
    }


    inline int write_SpikeWriter(int id, int pnodecount, SpikeWriter_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
}


using namespace coreneuron;


#ifndef CORENEURON_BUILD
#include <stdlib.h>
#ifndef DISABLE_MPI
#include <mpi.h>
#endif
#ifndef NRN_VERSION_GTEQ_8_2_0
extern double* vector_vec();
extern int vector_capacity();
extern void* vector_arg();
#endif
extern int nrnmpi_myid;
#endif  // CORENEURON_BUILD


namespace coreneuron {


    inline int write_SpikeWriter(int id, int pnodecount, SpikeWriter_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_write = 0;
        #ifndef CORENEURON_BUILD
                double *time = NULL, *gid = NULL;
                int num_spikes = 0;
                char *filePath = NULL;
                if (ifarg(1)) {
                    IvocVect *v = vector_arg(1);
                    time = vector_vec(v);
                    num_spikes = vector_capacity(v);
                }
                if (ifarg(2)) {
                    IvocVect *v = vector_arg(2);
                    gid = vector_vec(v);
                }
                if(ifarg(3)) {
                    filePath = hoc_gargstr(3);
                } else  {
        #ifndef DISABLE_MPI
                    if(nrnmpi_myid == 0) {
                        fprintf(stderr, " Error : No spike file path provided, can't write spikes! \n");
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
        #else
                    fprintf(stderr, " Error : No spike file path provided, can't write spikes! \n");
                    exit(-1);
        #endif
                }
                unsigned num_entries = nrnmpi_myid == 0 ? (num_spikes + 1) : num_spikes;
                const int spike_record_length = 48;
                unsigned num_bytes = (sizeof(char) * num_entries * spike_record_length) + 1;
                char *spike_data = (char *) malloc(num_bytes);
                if(spike_data == NULL) {
                    fprintf(stderr, "Error : Memory allocation failed for spike buffer I/O!\n");
        #ifndef DISABLE_MPI
                    MPI_Abort(MPI_COMM_WORLD, 1);
        #else
                    exit(-1);
        #endif
                }
                strcpy(spike_data, "");
                if(nrnmpi_myid == 0) {
                    strcat(spike_data, "/scatter\n");
                }
                int i;
                for(i = 0; i < num_spikes; i++) {
                    char str[spike_record_length];
                    int nstr = snprintf(str, spike_record_length, "%.3f\t%d\n", time[i], (int)gid[i]);
                    if (nstr >= spike_record_length) {
                        fprintf(stderr, "Error : Record written is larger than spike record buffer\n");
                        free(spike_data);
        #ifndef DISABLE_MPI
                        MPI_Abort(MPI_COMM_WORLD, 1);
        #else
                        exit(-1);
        #endif
                    }
                    strcat(spike_data, str);
                }
                unsigned long num_chars = strlen(spike_data);
        #ifndef DISABLE_MPI
                unsigned long offset = 0;
                MPI_Exscan(&num_chars, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
                if (nrnmpi_myid == 0) {
                    offset = 0;
                }
                MPI_File fh;
                MPI_Status status;
                MPI_File_delete(filePath, MPI_INFO_NULL);
                MPI_File_open(MPI_COMM_WORLD, filePath, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
                MPI_File_write_at_all(fh, offset, spike_data, num_chars, MPI_BYTE, &status);
                MPI_File_close(&fh);
        #else
                FILE *spike_file = fopen(filePath, "w");
                fprintf(spike_file, spike_data);
                fclose(spike_file);
        #endif
                free(spike_data);
        #endif  // CORENEURON_BUILD

        return ret_write;
    }


    /** initialize channel */
    void nrn_init_SpikeWriter(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<SpikeWriter_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma ivdep
            #pragma omp simd
            for (int id = 0; id < nodecount; id++) {
                double v = 0.0;
            }
        }
    }


    /** register channel with the simulator */
    void _SpikeWriter_reg() {

        int mech_type = nrn_get_mechtype("SpikeWriter");
        SpikeWriter_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        point_register_mech(mechanism, nrn_alloc_SpikeWriter, nullptr, nullptr, nullptr, nrn_init_SpikeWriter, nrn_private_constructor_SpikeWriter, nrn_private_destructor_SpikeWriter, first_pointer_var_index(), nullptr, nullptr, 1);

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "area");
        hoc_register_dparam_semantics(mech_type, 1, "pntproc");
        add_nrn_artcell(mech_type, 0);
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
