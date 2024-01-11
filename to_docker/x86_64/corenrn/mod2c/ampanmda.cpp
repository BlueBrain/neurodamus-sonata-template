/*********************************************************
Model Name      : AmpaNmda
Filename        : ampanmda.mod
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
        "AmpaNmda",
        "mg",
        "gmax",
        "ltdinvl",
        "ltpinvl",
        "sighalf",
        "sigslope",
        "sigexp",
        "x",
        "mgid",
        "ggid",
        "srcgid",
        "training",
        0,
        "i",
        "inmda",
        "iampa",
        "gnmda",
        "Rinf",
        "Rtau",
        0,
        "Ron",
        "Roff",
        "gampa",
        0,
        0
    };


    /** all global variables */
    struct AmpaNmda_Store {
        int point_type{};
        double Ron0{};
        double Roff0{};
        double gampa0{};
        int reset{};
        int mech_type{};
        double Cdur{1};
        double Alpha{0.35};
        double Beta{0.012};
        double E{0};
        double ampatau{3};
        double gampafactor{0.001};
        double nmdafactor{0.07};
        int slist1[3]{18, 19, 20};
        int dlist1[3]{22, 23, 24};
    };
    static_assert(std::is_trivially_copy_constructible_v<AmpaNmda_Store>);
    static_assert(std::is_trivially_move_constructible_v<AmpaNmda_Store>);
    static_assert(std::is_trivially_copy_assignable_v<AmpaNmda_Store>);
    static_assert(std::is_trivially_move_assignable_v<AmpaNmda_Store>);
    static_assert(std::is_trivially_destructible_v<AmpaNmda_Store>);
    AmpaNmda_Store AmpaNmda_global;


    /** all mechanism instance variables and global variables */
    struct AmpaNmda_Instance  {
        const double* mg{};
        const double* gmax{};
        const double* ltdinvl{};
        const double* ltpinvl{};
        const double* sighalf{};
        const double* sigslope{};
        const double* sigexp{};
        const double* x{};
        const double* mgid{};
        const double* ggid{};
        const double* srcgid{};
        const double* training{};
        double* i{};
        double* inmda{};
        double* iampa{};
        double* gnmda{};
        double* Rinf{};
        double* Rtau{};
        double* Ron{};
        double* Roff{};
        double* gampa{};
        double* synon{};
        double* DRon{};
        double* DRoff{};
        double* Dgampa{};
        double* v_unused{};
        double* g_unused{};
        double* tsave{};
        const double* node_area{};
        const int* point_process{};
        const int* tqitem{};
        AmpaNmda_Store* global{&AmpaNmda_global};
    };


    /** connect global (scalar) variables to hoc -- */
    static DoubScal hoc_scalar_double[] = {
        {"Cdur_AmpaNmda", &AmpaNmda_global.Cdur},
        {"Alpha_AmpaNmda", &AmpaNmda_global.Alpha},
        {"Beta_AmpaNmda", &AmpaNmda_global.Beta},
        {"E_AmpaNmda", &AmpaNmda_global.E},
        {"ampatau_AmpaNmda", &AmpaNmda_global.ampatau},
        {"gampafactor_AmpaNmda", &AmpaNmda_global.gampafactor},
        {"nmdafactor_AmpaNmda", &AmpaNmda_global.nmdafactor},
        {nullptr, nullptr}
    };


    /** connect global (array) variables to hoc -- */
    static DoubVec hoc_vector_double[] = {
        {nullptr, nullptr, 0}
    };


    static inline int first_pointer_var_index() {
        return -1;
    }


    static inline int num_net_receive_args() {
        return 6;
    }


    static inline int float_variables_size() {
        return 28;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return AmpaNmda_global.mech_type;
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
    static void nrn_private_constructor_AmpaNmda(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new AmpaNmda_Instance{};
        assert(inst->global == &AmpaNmda_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(AmpaNmda_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_AmpaNmda(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<AmpaNmda_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &AmpaNmda_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(AmpaNmda_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<AmpaNmda_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &AmpaNmda_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(AmpaNmda_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->mg = ml->data+0*pnodecount;
        inst->gmax = ml->data+1*pnodecount;
        inst->ltdinvl = ml->data+2*pnodecount;
        inst->ltpinvl = ml->data+3*pnodecount;
        inst->sighalf = ml->data+4*pnodecount;
        inst->sigslope = ml->data+5*pnodecount;
        inst->sigexp = ml->data+6*pnodecount;
        inst->x = ml->data+7*pnodecount;
        inst->mgid = ml->data+8*pnodecount;
        inst->ggid = ml->data+9*pnodecount;
        inst->srcgid = ml->data+10*pnodecount;
        inst->training = ml->data+11*pnodecount;
        inst->i = ml->data+12*pnodecount;
        inst->inmda = ml->data+13*pnodecount;
        inst->iampa = ml->data+14*pnodecount;
        inst->gnmda = ml->data+15*pnodecount;
        inst->Rinf = ml->data+16*pnodecount;
        inst->Rtau = ml->data+17*pnodecount;
        inst->Ron = ml->data+18*pnodecount;
        inst->Roff = ml->data+19*pnodecount;
        inst->gampa = ml->data+20*pnodecount;
        inst->synon = ml->data+21*pnodecount;
        inst->DRon = ml->data+22*pnodecount;
        inst->DRoff = ml->data+23*pnodecount;
        inst->Dgampa = ml->data+24*pnodecount;
        inst->v_unused = ml->data+25*pnodecount;
        inst->g_unused = ml->data+26*pnodecount;
        inst->tsave = ml->data+27*pnodecount;
        inst->node_area = nt->_data;
        inst->point_process = ml->pdata;
        inst->tqitem = ml->pdata;
    }



    static void nrn_alloc_AmpaNmda(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_AmpaNmda(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<AmpaNmda_Instance*>(ml->instance);

        #endif
    }


    void nrn_destructor_AmpaNmda(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<AmpaNmda_Instance*>(ml->instance);

        #endif
    }


    inline double mgblock_AmpaNmda(int id, int pnodecount, AmpaNmda_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v);
    inline double plast_AmpaNmda(int id, int pnodecount, AmpaNmda_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double step);


    inline double mgblock_AmpaNmda(int id, int pnodecount, AmpaNmda_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double arg_v) {
        double ret_mgblock = 0.0;
        ret_mgblock = 1.0 / (1.0 + exp(0.062 *  -arg_v) * (inst->mg[id] / 3.57));
        return ret_mgblock;
    }


    inline double plast_AmpaNmda(int id, int pnodecount, AmpaNmda_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v, double step) {
        double ret_plast = 0.0;
        ret_plast = pow((1.0 - 1.0 / (1.0 + exp((step - inst->sighalf[id]) / inst->sigslope[id]))), inst->sigexp[id]);
        return ret_plast;
    }


    static inline void net_send_buffering(const NrnThread* nt, NetSendBuffer_t* nsb, int type, int vdata_index, int weight_index, int point_index, double t, double flag) {
        int i = 0;
        i = nsb->_cnt++;
        if (i >= nsb->_size) {
            nsb->grow();
        }
        if (i < nsb->_size) {
            nsb->_sendtype[i] = type;
            nsb->_vdata_index[i] = vdata_index;
            nsb->_weight_index[i] = weight_index;
            nsb->_pnt_index[i] = point_index;
            nsb->_nsb_t[i] = t;
            nsb->_nsb_flag[i] = flag;
        }
    }


    /** initialize block for net receive */
    static void net_init(Point_process* pnt, int weight_index, double flag) {
        int tid = pnt->_tid;
        int id = pnt->_i_instance;
        double v = 0;
        NrnThread* nt = nrn_threads + tid;
        Memb_list* ml = nt->_ml_list[pnt->_type];
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        double* data = ml->data;
        double* weights = nt->weights;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<AmpaNmda_Instance*>(ml->instance);

        double* weight = weights + weight_index + 0;
        double* s = weights + weight_index + 1;
        double* w = weights + weight_index + 2;
        double* tlast = weights + weight_index + 3;
        double* r0 = weights + weight_index + 4;
        double* t0 = weights + weight_index + 5;
        if ((*s) == 0.0) {
            (*w) = 0.0;
        } else {
            double plast_in_0;
            {
                double step_in_0;
                step_in_0 = (*s);
                plast_in_0 = pow((1.0 - 1.0 / (1.0 + exp((step_in_0 - inst->sighalf[id]) / inst->sigslope[id]))), inst->sigexp[id]);
            }
            (*w) = (*weight) * plast_in_0;
        }
        (*tlast) =  -1e9;
        (*r0) = 0.0;
        (*t0) =  -1e9;
        auto& nsb = ml->_net_send_buffer;
    }


    static inline void net_receive_kernel_AmpaNmda(double t, Point_process* pnt, AmpaNmda_Instance* inst, NrnThread* nt, Memb_list* ml, int weight_index, double flag) {
        int tid = pnt->_tid;
        int id = pnt->_i_instance;
        double v = 0;
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        double* data = ml->data;
        double* weights = nt->weights;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        double* weight = weights + weight_index + 0;
        double* s = weights + weight_index + 1;
        double* w = weights + weight_index + 2;
        double* tlast = weights + weight_index + 3;
        double* r0 = weights + weight_index + 4;
        double* t0 = weights + weight_index + 5;
        inst->tsave[id] = t;
        {
            if (flag == 0.0) {
                if (t - (*tlast) < inst->ltpinvl[id]) {
                    if (inst->training[id]) {
                        (*s) = (*s) + 1.0;
                        if ((*s) > 75.0) {
                            (*s) = 75.0;
                        }
                    }
                } else if (t - (*tlast) > inst->ltdinvl[id]) {
                } else {
                    if (inst->training[id]) {
                        (*s) = (*s) - 1.0;
                        if ((*s) < 0.0) {
                            (*s) = 0.0;
                        }
                    }
                }
                (*tlast) = t;
                if ((*s) == 0.0) {
                    (*w) = 0.0;
                } else {
                    double plast_in_1;
                    {
                        double step_in_1;
                        step_in_1 = (*s);
                        plast_in_1 = pow((1.0 - 1.0 / (1.0 + exp((step_in_1 - inst->sighalf[id]) / inst->sigslope[id]))), inst->sigexp[id]);
                    }
                    (*w) = (*weight) * plast_in_1;
                }
                inst->gampa[id] = inst->gampa[id] + (*w) * inst->gmax[id] * inst->global->gampafactor;
                (*r0) = (*r0) * exp( -inst->global->Beta * (t - (*t0)));
                (*t0) = t;
                inst->synon[id] = inst->synon[id] + (*w);
                inst->Ron[id] = inst->Ron[id] + (*r0);
                inst->Roff[id] = inst->Roff[id] - (*r0);
                net_send_buffering(nt, ml->_net_send_buffer, 0, inst->tqitem[2*pnodecount+id], weight_index, inst->point_process[1*pnodecount+id], t+inst->global->Cdur, (*w) + 1.0);
            } else {
                (*r0) = (flag - 1.0) * inst->Rinf[id] + ((*r0) - (flag - 1.0) * inst->Rinf[id]) * exp( -(t - (*t0)) / inst->Rtau[id]);
                (*t0) = t;
                inst->synon[id] = inst->synon[id] - (flag - 1.0);
                inst->Ron[id] = inst->Ron[id] - (*r0);
                inst->Roff[id] = inst->Roff[id] + (*r0);
            }
        }
    }


    static void net_receive_AmpaNmda(Point_process* pnt, int weight_index, double flag) {
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


    void net_buf_receive_AmpaNmda(NrnThread* nt) {
        Memb_list* ml = get_memb_list(nt);
        if (!ml) {
            return;
        }

        NetReceiveBuffer_t* nrb = ml->_net_receive_buffer;
        auto* const inst = static_cast<AmpaNmda_Instance*>(ml->instance);
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
                net_receive_kernel_AmpaNmda(t, point_process, inst, nt, ml, weight_index, flag);
            }
        }
        nrb->_displ_cnt = 0;
        nrb->_cnt = 0;

        NetSendBuffer_t* nsb = ml->_net_send_buffer;
        for (int i=0; i < nsb->_cnt; i++) {
            int type = nsb->_sendtype[i];
            int tid = nt->id;
            double t = nsb->_nsb_t[i];
            double flag = nsb->_nsb_flag[i];
            int vdata_index = nsb->_vdata_index[i];
            int weight_index = nsb->_weight_index[i];
            int point_index = nsb->_pnt_index[i];
            net_sem_from_gpu(type, vdata_index, weight_index, tid, point_index, t, flag);
        }
        nsb->_cnt = 0;
    }


    /** initialize channel */
    void nrn_init_AmpaNmda(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<AmpaNmda_Instance*>(ml->instance);

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
                inst->Ron[id] = inst->global->Ron0;
                inst->Roff[id] = inst->global->Roff0;
                inst->gampa[id] = inst->global->gampa0;
                inst->Rinf[id] = inst->global->Alpha / (inst->global->Alpha + inst->global->Beta);
                inst->Rtau[id] = 1.0 / (inst->global->Alpha + inst->global->Beta);
                inst->synon[id] = 0.0;
                inst->gampa[id] = 0.0;
            }
        }

        NetSendBuffer_t* nsb = ml->_net_send_buffer;
        for (int i=0; i < nsb->_cnt; i++) {
            int type = nsb->_sendtype[i];
            int tid = nt->id;
            double t = nsb->_nsb_t[i];
            double flag = nsb->_nsb_flag[i];
            int vdata_index = nsb->_vdata_index[i];
            int weight_index = nsb->_weight_index[i];
            int point_index = nsb->_pnt_index[i];
            net_sem_from_gpu(type, vdata_index, weight_index, tid, point_index, t, flag);
        }
        nsb->_cnt = 0;
    }


    inline double nrn_current_AmpaNmda(int id, int pnodecount, AmpaNmda_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double current = 0.0;
        double mgblock_in_0;
        {
            double v_in_0;
            v_in_0 = v;
            mgblock_in_0 = 1.0 / (1.0 + exp(0.062 *  -v_in_0) * (inst->mg[id] / 3.57));
        }
        inst->gnmda[id] = mgblock_in_0 * (inst->Ron[id] + inst->Roff[id]) * inst->gmax[id] * inst->global->nmdafactor;
        inst->inmda[id] = inst->gnmda[id] * (v - inst->global->E);
        inst->iampa[id] = inst->gampa[id] * (v - inst->global->E);
        inst->i[id] = inst->iampa[id] + inst->inmda[id];
        current += inst->i[id];
        return current;
    }


    /** update current */
    void nrn_cur_AmpaNmda(NrnThread* nt, Memb_list* ml, int type) {
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
        auto* const inst = static_cast<AmpaNmda_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            double g = nrn_current_AmpaNmda(id, pnodecount, inst, data, indexes, thread, nt, v+0.001);
            double rhs = nrn_current_AmpaNmda(id, pnodecount, inst, data, indexes, thread, nt, v);
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
    void nrn_state_AmpaNmda(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<AmpaNmda_Instance*>(ml->instance);

        #pragma omp simd
        #pragma ivdep
        for (int id = 0; id < nodecount; id++) {
            int node_id = node_index[id];
            double v = voltage[node_id];
            #if NRN_PRCELLSTATE
            inst->v_unused[id] = v;
            #endif
            inst->Ron[id] = inst->Ron[id] + (1.0 - exp(nt->_dt * (((( -1.0))) / inst->Rtau[id]))) * ( -((((inst->synon[id]) * (inst->Rinf[id]))) / inst->Rtau[id]) / (((( -1.0))) / inst->Rtau[id]) - inst->Ron[id]);
            inst->Roff[id] = inst->Roff[id] + (1.0 - exp(nt->_dt * (( -inst->global->Beta) * (1.0)))) * ( -(0.0) / (( -inst->global->Beta) * (1.0)) - inst->Roff[id]);
            inst->gampa[id] = inst->gampa[id] + (1.0 - exp(nt->_dt * (( -1.0) / inst->global->ampatau))) * ( -(0.0) / (( -1.0) / inst->global->ampatau) - inst->gampa[id]);
        }
    }


    /** register channel with the simulator */
    void _ampanmda_reg() {

        int mech_type = nrn_get_mechtype("AmpaNmda");
        AmpaNmda_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        point_register_mech(mechanism, nrn_alloc_AmpaNmda, nrn_cur_AmpaNmda, nullptr, nrn_state_AmpaNmda, nrn_init_AmpaNmda, nrn_private_constructor_AmpaNmda, nrn_private_destructor_AmpaNmda, first_pointer_var_index(), nullptr, nullptr, 1);

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "area");
        hoc_register_dparam_semantics(mech_type, 1, "pntproc");
        hoc_register_dparam_semantics(mech_type, 2, "netsend");
        hoc_register_net_receive_buffering(net_buf_receive_AmpaNmda, mech_type);
        set_pnt_receive(mech_type, net_receive_AmpaNmda, net_init, num_net_receive_args());
        hoc_register_net_send_buffering(mech_type);
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
