/*********************************************************
Model Name      : HDF5Reader
Filename        : HDF5reader.mod
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
        "HDF5Reader",
        0,
        0,
        0,
        "ptr",
        0
    };


    /** all global variables */
    struct HDF5Reader_Store {
        int point_type{};
        int reset{};
        int mech_type{};
    };
    static_assert(std::is_trivially_copy_constructible_v<HDF5Reader_Store>);
    static_assert(std::is_trivially_move_constructible_v<HDF5Reader_Store>);
    static_assert(std::is_trivially_copy_assignable_v<HDF5Reader_Store>);
    static_assert(std::is_trivially_move_assignable_v<HDF5Reader_Store>);
    static_assert(std::is_trivially_destructible_v<HDF5Reader_Store>);
    HDF5Reader_Store HDF5Reader_global;


    /** all mechanism instance variables and global variables */
    struct HDF5Reader_Instance  {
        double* v_unused{};
        double* tsave{};
        const double* node_area{};
        void** point_process{};
        double* ptr{};
        HDF5Reader_Store* global{&HDF5Reader_global};
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


    static inline int num_net_receive_args() {
        return 1;
    }


    static inline int float_variables_size() {
        return 2;
    }


    static inline int int_variables_size() {
        return 3;
    }


    static inline int get_mech_type() {
        return HDF5Reader_global.mech_type;
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
    static void nrn_private_constructor_HDF5Reader(NrnThread* nt, Memb_list* ml, int type) {
        assert(!ml->instance);
        assert(!ml->global_variables);
        assert(ml->global_variables_size == 0);
        auto* const inst = new HDF5Reader_Instance{};
        assert(inst->global == &HDF5Reader_global);
        ml->instance = inst;
        ml->global_variables = inst->global;
        ml->global_variables_size = sizeof(HDF5Reader_Store);
    }

    // Deallocate the instance structure
    static void nrn_private_destructor_HDF5Reader(NrnThread* nt, Memb_list* ml, int type) {
        auto* const inst = static_cast<HDF5Reader_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &HDF5Reader_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(HDF5Reader_Store));
        delete inst;
        ml->instance = nullptr;
        ml->global_variables = nullptr;
        ml->global_variables_size = 0;
    }

    /** initialize mechanism instance variables */
    static inline void setup_instance(NrnThread* nt, Memb_list* ml) {
        auto* const inst = static_cast<HDF5Reader_Instance*>(ml->instance);
        assert(inst);
        assert(inst->global);
        assert(inst->global == &HDF5Reader_global);
        assert(inst->global == ml->global_variables);
        assert(ml->global_variables_size == sizeof(HDF5Reader_Store));
        int pnodecount = ml->_nodecount_padded;
        Datum* indexes = ml->pdata;
        inst->v_unused = ml->data+0*pnodecount;
        inst->tsave = ml->data+1*pnodecount;
        inst->node_area = nt->_data;
        inst->point_process = nt->_vdata;
        inst->ptr = nt->_data;
    }



    static void nrn_alloc_HDF5Reader(double* data, Datum* indexes, int type) {
        // do nothing
    }


    void nrn_constructor_HDF5Reader(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<HDF5Reader_Instance*>(ml->instance);

         {
        #ifndef CORENEURON_BUILD
        #ifdef DISABLE_HDF5
            if(ifarg(1)) {
                fprintf(stderr, "HDF5 support is not available\n");
                exit(-1);
            }
        #else
            char nameoffile[512];
            int nFiles = 1;
            if( ifarg(2) ) {
                nFiles = *getarg(2);
            }
            if(ifarg(1) && hoc_is_str_arg(1)) {
                INFOCAST;
                Info* info = 0;
                strncpy(nameoffile, gargstr(1),512);
                info = (Info*) hoc_Emalloc(sizeof(Info)); hoc_malchk();
                initInfo( info );
                if( ifarg(3) ) {
                    info->verboseLevel = *getarg(3);
                }
                *ip = info;
                if( nFiles == 1 ) {
                    hid_t plist_id = 0;
                    if( ifarg(2) ) {
                        if( nrnmpi_myid == 0 ) { fprintf( stderr, "using parallel hdf5 is disabled\n" ); }
                    }
                    openFile( info, nameoffile, -1, 0, 0, 0 );
                }
                else {
                    info->synapseCatalog.rootName = strdup( nameoffile );
                    int mpi_size=1, mpi_rank=0;
        #ifndef DISABLE_MPI
                    MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
                    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
        #endif
                    int plength = strlen(nameoffile);
                    char *nrnPath = (char*) malloc( plength + 128 );
                    strncpy( nrnPath, nameoffile, plength+1 );
                    char *term = strrchr( nrnPath, '/' ); *term = 0;
                    strcat( nrnPath, "/mergeAllH5.sh" );
                    FILE *fin = fopen( nrnPath, "r" );
                    if( !fin ) {
                        *term = 0;
                        strcat( nrnPath, "/merge_nrn_positions.sh" );
                        fin = fopen( nrnPath, "r" );
                    }
                    free(nrnPath);
                    if( fin ) {
                        info->synapseCatalog.directFile = fin;
                        return;
                    }
                    if( mpi_size > nFiles ) { 
                        int nRanksPerFile = (int) mpi_size / nFiles;
                        if( mpi_size % nFiles != 0 )
                            nRanksPerFile++;
                        if( nRanksPerFile * nFiles > mpi_size ) { 
                            info->file_ = -1;
                            return;
                        }
                        int fileIndex = (int) mpi_rank / nRanksPerFile;  
                        int startRank = fileIndex * nRanksPerFile;
                        openFile( info, nameoffile, fileIndex, nRanksPerFile, startRank, mpi_rank );
                    } else {
                        int nFilesPerRank = nFiles / mpi_size;
                        if( nFiles % mpi_size != 0 ) {
                            nFilesPerRank++;
                        }
                        int startFile = mpi_rank * nFilesPerRank;
                        if( startFile+nFilesPerRank > nFiles ) {
                            nFilesPerRank = nFiles-startFile;
                            if( nFilesPerRank <= 0 ) {
                                info->file_ = -1;
                                return;
                            }
                        }
                        int fileIndex=0;
                        for( fileIndex=0; fileIndex<nFilesPerRank; fileIndex++ ) {
                            openFile( info, nameoffile, startFile+fileIndex, 1, 0, 0 );
                        }
                    }
                }
            }
        #endif  // DISABLE_HDF5
        #endif  // CORENEURON_BUILD
        }

        #endif
    }


    void nrn_destructor_HDF5Reader(NrnThread* nt, Memb_list* ml, int type) {
        #ifndef CORENEURON_BUILD
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;
        auto* const inst = static_cast<HDF5Reader_Instance*>(ml->instance);

         {
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST; Info* info = *ip;
            if(info->file_>=0)
            {
                if( info->acc_tpl1 != -1 ) {
                    if( nrnmpi_myid == 0 ) { fprintf( stderr, "terminating parallel h5 access\n" ); }
                    H5Pclose(info->acc_tpl1);
                }
                H5Fclose(info->file_);
                info->file_ = -1;
            }
            if(info->datamatrix_ != NULL)
            {
                free(info->datamatrix_);
                info->datamatrix_ = NULL;
            }
            if( info->datavector_ != NULL )
            {
                free( info->datavector_ );
                info->datavector_ = NULL;
            }
        #endif
        #endif  // CORENEURON_BUILD
        }

        #endif
    }


    inline double redirect_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double checkVersion_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double loadData_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double getNoOfColumns_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double numberofrows_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double getData_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double getDataInt_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double getColumnDataRange_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double getColumnData_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double getAttributeValue_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double closeFile_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline double exchangeSynapseLocations_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int getDimensions_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
    inline int getDataString_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v);
}


using namespace coreneuron;


#ifndef CORENEURON_BUILD
#if defined(NRN_VERSION_GTEQ)
#if NRN_VERSION_GTEQ(9,0,0)
#define NRN_VERSION_GTEQ_9_0_0
#endif
#endif
#ifndef DISABLE_HDF5
#undef ptr
#define H5_USE_16_API 1
#undef dt
#include "hdf5.h"
#define nt->_dt nrn_threads->_dt
#ifndef DISABLE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
#ifndef NRN_VERSION_GTEQ_8_2_0
extern double* hoc_pgetarg(int iarg);
extern double* getarg(int iarg);
extern char* gargstr(int iarg);
extern int hoc_is_str_arg(int iarg);
extern int ifarg(int iarg);
extern double chkarg(int iarg, double low, double high);
extern double* vector_vec(void* vv);
extern int vector_capacity(void* vv);
extern void* vector_arg(int);
#endif
extern int nrnmpi_numprocs;
extern int nrnmpi_myid;
/**
 * During synapse loading initialization, the h5 files with synapse data are catalogged
 * such that each cpu looks at a subset of what's available and gets ready to report to
 * other cpus where to find postsynaptic cells
 */
struct SynapseCatalog
{
    char* rootName;
    FILE *directFile;
    int fileID;
    int nFiles;
    int *fileIDs;
    int *availablegidCount;
    int **availablegids;
    int gidIndex;
};
typedef struct SynapseCatalog SynapseCatalog;
/**
 * After loading, the cpus will exchange requests for info about which files contain synapse
 * data for their local gids.  This is used to store the received info provided those gids receive
 */
struct ConfirmedCells
{
    int count;
    int *gids;
    int *fileIDs;
};
typedef struct ConfirmedCells ConfirmedCells;
#define NONE 0
#define FLOAT_MATRIX 1
#define LONG_VECTOR 2
/**
 * Hold persistent HDF5 data such as file handle and info about latest dataset loaded
 */
struct Info {
    hid_t file_;
    float * datamatrix_;
    long long *datavector_;
    char name_group[256];
    hsize_t rowsize_;
    hsize_t columnsize_;
    hid_t acc_tpl1;
    int mode;
    int verboseLevel;
    SynapseCatalog synapseCatalog;
    ConfirmedCells confirmedCells;
};
typedef struct Info Info;
#define INFOCAST Info** ip = (Info**)(&(nt->_data[indexes[2*pnodecount + id]]))
#define dp double*
/**
 * Utility function to ensure that all members of an Info struct are initialized.
 */
void initInfo( Info *info )
{
    info->file_ = -1;
    info->datamatrix_ = NULL;
    info->datavector_ = NULL;
    info->mode = NONE;
    info->name_group[0] = '\0';
    info->rowsize_ = 0;
    info->columnsize_ = 0;
    info->acc_tpl1 = -1;
    info->verboseLevel = 0;
    info->synapseCatalog.rootName = NULL;
    info->synapseCatalog.directFile = NULL;
    info->synapseCatalog.fileID = -1;
    info->synapseCatalog.fileIDs = NULL;
    info->synapseCatalog.nFiles = 0;
    info->synapseCatalog.availablegidCount = NULL;
    info->synapseCatalog.availablegids = NULL;
    info->synapseCatalog.gidIndex = 0;
    info->confirmedCells.count = 0;
    info->confirmedCells.gids = NULL;
    info->confirmedCells.fileIDs = NULL;
}
/**
 * Use to sort gids for look up during call of readDirectMapping func
 */
int gidcompare( const void *a, const void *b ) {
    return ( *(double*)a - *(double*)b );
}
/**
 * Use to confirm a gid is present on local cpu to help readDirectMapping func.  use binary search on sorted gidList
 */
int gidexists( double searchgid, int ngids, double *gidList ) {
    int first = 0;
    int last = ngids-1;
    int middle = (first+last)/2;
    while( first <= last ) {
        if( gidList[middle] < searchgid ) {
            first = middle + 1;
        } else if( gidList[middle] == searchgid ) {
            return 1;
        } else {
            last = middle - 1;
        }
        middle = (first+last)/2;
    }
    return 0;
}
/**
 * Use Circuit Builder output script to know which files hold this cpu's gid synapse data.  Should be
 * short term solution until better loading method is available.
 */
void readDirectMapping( Info *info, int ngids, double *gidList ) {
    double *tmpList = (double*) malloc( ngids * sizeof(double) );
    memcpy( tmpList, gidList, ngids*sizeof(double) );
    qsort( tmpList, ngids, sizeof(double), gidcompare );
    info->confirmedCells.count = ngids;
    info->confirmedCells.gids = (int*) malloc(ngids*sizeof(int));
    info->confirmedCells.fileIDs = (int*) malloc(ngids*sizeof(int));
    FILE *fin = info->synapseCatalog.directFile;
    info->synapseCatalog.directFile = NULL;
    char line[1024];
    char *res = fgets( line, 1024, fin );
    int gid = -1;
    int gidIndex = 0;
    while( res != NULL ) {
        if( strncmp( line, "$CMD", 4 ) == 0 ) {
            char *gidloc = strstr( line, "/a" );
            if( gidloc != NULL ) {
                gid = atoi(&gidloc[2]);
            }
            if( gidexists( (double) gid, ngids, tmpList ) ) {
                int fileID = atoi(&line[12]);
                info->confirmedCells.gids[gidIndex] = gid;
                info->confirmedCells.fileIDs[gidIndex] = fileID;
                gidIndex++;
            }
        }
        res = fgets( line, 1024, fin );
    }
    fclose( fin );
    free(tmpList);
}
/**
 * Callback function for H5Giterate - if the dataset opened corresponds to a gid, it is catalogged so the
 * local cpu can inform other cpus the whereabouts of that gid
 *
 * @param loc_id hdf5 handle to the open file
 * @param name name of the dataset to be accessed during this iteration step
 * @param opdata not used since we have global Info object
 */
herr_t loadShareData( hid_t loc_id, const char *name, void *opdata )
{
    INFOCAST;
    Info* info = *ip;
    assert( info->file_ >= 0 );
    int gid = atoi( name+1 );
    char rebuild[32];
    snprintf( rebuild, 32, "a%d", gid );
    if( strcmp( rebuild, name ) != 0 ) {
        return 0;
    }
    int fileIndex = info->synapseCatalog.nFiles-1;
    info->synapseCatalog.availablegids[ fileIndex ][ info->synapseCatalog.gidIndex++ ] = gid;
    return 1;
}
/**
 * Open an HDF5 file for reading.  In the event of synapse data, the datasets of the file may be iterated in order to
 * build a catalog of available gids and their file locations.
 *
 * @param info Structure that manages hdf5 info
 * @param filename File to open
 * @param fileID Integer to identify this file (attached as suffix to filename)
 * @param nNodesPerFile 0: open file, but don't load data; 1: open file for catalogging; N: read portion of file for catalogging
 * @param startRank used to help calculate data range to load when file subportion is loaded
 * @param myRank used to help calculate data range to load when file subportion is loaded
 */
int openFile( Info* info, const char *filename, int fileID, int nRanksPerFile, int startRank, int myRank )
{
    if( info->file_ != -1 ) {
        H5Fclose(info->file_);
    }
    char nameoffile[512];
    if( fileID != -1 ) {
        snprintf( nameoffile, 512, "%s.%d", filename, fileID );
    } else {
        strncpy( nameoffile, filename, 512 );
    }
    info->name_group[0]='\0';
    hid_t file_driver = (info->acc_tpl1 != -1)? info->acc_tpl1 : H5P_DEFAULT;
    info->file_ = H5Fopen( nameoffile, H5F_ACC_RDONLY, file_driver);
    int result = (info->file_ < 0);
    int failed = result;
#ifndef DISABLE_MPI
    if( info->acc_tpl1 != -1 ) {
        MPI_Allreduce( &result, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
        if( failed ) {
            int canreport = (result<0)?nrnmpi_myid:nrnmpi_numprocs, willreport = 0;
            MPI_Allreduce( &canreport, &willreport, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
            if( willreport == nrnmpi_myid ) {
                fprintf(stderr, "%d ERROR: %d ranks failed collective open of synapse file: %s\n", nrnmpi_myid, failed, nameoffile );
            }
            info->file_ = -1;
            H5Eprint(stderr);
            return -1;
        }
    } else  
#endif
    if( failed ) {
        info->file_ = -1;
        fprintf(stderr, "ERROR: Failed to open synapse file: %s\n", nameoffile );
        H5Eprint(stderr);
        return -1;
    }
    if( nRanksPerFile == 0 ) {
        return 0;
    }
    info->synapseCatalog.fileID = fileID;
    int nDatasetsToImport=0, startIndex=0;
    hsize_t nObjects;
    H5Gget_num_objs( info->file_, &nObjects );
    if( nRanksPerFile == 1 ) {
        nDatasetsToImport = (int) nObjects;
    }
    else {
        nDatasetsToImport = (int) nObjects / nRanksPerFile;
        if( nObjects%nRanksPerFile != 0 )
            nDatasetsToImport++;
        startIndex = (myRank-startRank)*nDatasetsToImport;
        if( startIndex + nDatasetsToImport > (int) nObjects ) {
            nDatasetsToImport = (int) nObjects - startIndex;
            if( nDatasetsToImport <= 0 ) {
                return 0;
            }
        }
    }
    int nFiles = ++info->synapseCatalog.nFiles;
    info->synapseCatalog.fileIDs = (int*) realloc ( info->synapseCatalog.fileIDs, sizeof(int)*nFiles );
    info->synapseCatalog.fileIDs[nFiles-1] = fileID;
    info->synapseCatalog.availablegidCount = (int*) realloc ( info->synapseCatalog.availablegidCount, sizeof(int)*nFiles );
    info->synapseCatalog.availablegids = (int**) realloc ( info->synapseCatalog.availablegids, sizeof(int*)*nFiles );
    info->synapseCatalog.availablegidCount[nFiles-1] = nDatasetsToImport;
    info->synapseCatalog.availablegids[nFiles-1] = (int*) calloc ( nDatasetsToImport, sizeof(int) );
    info->synapseCatalog.gidIndex=0;
    int i, verify=startIndex;
    for( i=startIndex; i<startIndex+nDatasetsToImport && i<nObjects; i++ ) {
        assert( verify == i );
        result = H5Giterate( info->file_, "/", &verify, loadShareData, NULL );
        if( result != 1 )
            continue;
    }
    return 0;
}
/**
 * Load a dataset so that the dimensions are available, but don't retrieve any data
 *
 * @param info Structure that manages hdf5 info, its datamatrix_ variable is populated with hdf5 data on success
 * @param name The name of the dataset to access
 * @return 0 on success, < 0 on error
 */
int loadDimensions( Info *info, char* name )
{
    int isCurrentlyLoaded = strncmp( info->name_group, name, 256 ) == 0;
    if( isCurrentlyLoaded )
        return 0;
    hsize_t dims[2] = {0}, offset[2] = {0};
    hid_t dataset_id, dataspace;
    if( H5Lexists(info->file_, name, H5P_DEFAULT) == 0)
    {
        fprintf(stderr, "Error accessing to dataset %s in synapse file\n", name);
        return -1;
    }
    dataset_id = H5Dopen(info->file_, name);
    strncpy(info->name_group, name, 256);
    dataspace = H5Dget_space(dataset_id);
    int dimensions = H5Sget_simple_extent_ndims(dataspace);
    H5Sget_simple_extent_dims(dataspace,dims,NULL);
    info->rowsize_ = (unsigned long)dims[0];
    if( dimensions > 1 )
        info->columnsize_ = dims[1];
    else
        info->columnsize_ = 1;
    H5Sclose(dataspace);
    H5Dclose(dataset_id);
    return 0;
}
/**
 * Given the name of a dataset, load it from the current hdf5 file into the matrix pointer
 *
 * @param info Structure that manages hdf5 info, its datamatrix_ variable is populated with hdf5 data on success
 * @param name The name of the dataset to access and load in the hdf5 file
 */
int loadDataMatrix( Info *info, char* name )
{
    int isCurrentlyLoaded = strncmp( info->name_group, name, 256 ) == 0;
    if( isCurrentlyLoaded )
        return 0;
    hsize_t dims[2] = {0}, offset[2] = {0};
    hid_t dataset_id, dataspace;
    if( H5Lexists(info->file_, name, H5P_DEFAULT) == 0) {
        return -1;
    }
    dataset_id = H5Dopen(info->file_, name);
    strncpy(info->name_group, name, 256);
    dataspace = H5Dget_space(dataset_id);
    int dimensions = H5Sget_simple_extent_ndims(dataspace);
    H5Sget_simple_extent_dims(dataspace,dims,NULL);
    info->rowsize_ = (unsigned long)dims[0];
    if( dimensions > 1 )
        info->columnsize_ = dims[1];
    else
        info->columnsize_ = 1;
    if(info->datamatrix_ != NULL)
    {
        free(info->datamatrix_);
    }
    info->datamatrix_ = (float *) malloc(sizeof(float) *(info->rowsize_*info->columnsize_));
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
    hid_t dataspacetogetdata=H5Screate_simple(dimensions,dims,NULL);
    H5Dread(dataset_id,H5T_NATIVE_FLOAT,dataspacetogetdata,H5S_ALL,H5P_DEFAULT,info->datamatrix_);
    H5Sclose(dataspace);
    H5Sclose(dataspacetogetdata);
    H5Dclose(dataset_id);
    return 0;
}
/**
 * Given the name of a dataset with id values, load it from the current hdf5 file into the matrix pointer
 *
 * @param info Structure that manages hdf5 info, its datavector_ variable is populated with hdf5 data on success
 * @param name The name of the dataset to access and load in the hdf5 file
 */
int loadDataVector( Info *info, char* name )
{
    hsize_t dims[2] = {0}, offset[2] = {0};
    hid_t dataset_id, dataspace;
    if( H5Lexists(info->file_, name, H5P_DEFAULT) == 0)
    {
        fprintf(stderr, "Error accessing to dataset %s in synapse file\n", name);
        return -1;
    }
    dataset_id = H5Dopen(info->file_, name);
    dataspace = H5Dget_space(dataset_id);
    strncpy(info->name_group, name, 256);
    int dimensions = H5Sget_simple_extent_ndims(dataspace);
    H5Sget_simple_extent_dims(dataspace,dims,NULL);
    info->rowsize_ = (unsigned long)dims[0];
    if( dimensions > 1 )
        info->columnsize_ = dims[1];
    else
        info->columnsize_ = 1;
    if(info->datavector_ != NULL) {
        free(info->datavector_);
    }
    info->datavector_ = (long long *) malloc(sizeof(long long)*(info->rowsize_*info->columnsize_));
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
    hid_t dataspacetogetdata=H5Screate_simple(dimensions,dims,NULL);
    H5Dread(dataset_id,H5T_NATIVE_ULONG,dataspacetogetdata,H5S_ALL,H5P_DEFAULT,info->datavector_);
    H5Sclose(dataspace);
    H5Sclose(dataspacetogetdata);
    H5Dclose(dataset_id);
    return 0;
}
/**
 * Load an individual value from a dataset
 *
 * @param info Shared data
 * @param name dataset to open and read from
 * @param row  data item to retrieve
 * @param dest int address where data is to be stored
 */
int loadDataInt( Info* info, char* name, hid_t row, int *dest )
{
    #if defined(NRN_VERSION_GTEQ_9_0_0)
    hsize_t dims[1] = {1}, offset[1] = {static_cast<hsize_t>(row)}, offset_out[1] = {0}, count[1] = {1};
    #else
    hsize_t dims[1] = {1}, offset[1] = {row}, offset_out[1] = {0}, count[1] = {1};
    #endif
    hid_t dataset_id, dataspace, memspace, space; 
    herr_t status;
    int ndims = 0;
    long long temp;
    if( H5Lexists(info->file_, name, H5P_DEFAULT) == 0) {
        fprintf(stderr, "Error accessing to dataset %s in h5 file\n", name);
        return -1;
    }
    dataset_id = H5Dopen(info->file_, name);
    dataspace = H5Dget_space(dataset_id);
    space = H5Dget_space (dataset_id);
    ndims = H5Sget_simple_extent_dims (space, dims, NULL);
    status = H5Sselect_hyperslab( space, H5S_SELECT_SET, offset, NULL, count, NULL );
    memspace = H5Screate_simple( 1, dims, NULL );
    status = H5Sselect_hyperslab( memspace, H5S_SELECT_SET, offset_out, NULL, count, NULL );
    status = H5Dread (dataset_id, H5T_NATIVE_ULONG, memspace, space, H5P_DEFAULT, &temp );
    *dest = temp;
    status = H5Sclose (space);
    status = H5Sclose(dataspace);
    status = H5Dclose(dataset_id);
    return 0;
}
/**
 * Given the name of a dataset that should contain variable length string data, load those strings
 */
int loadDataString( Info* info, char* name, hid_t row, char **hoc_dest )
{
    #if defined(NRN_VERSION_GTEQ_9_0_0)
    hsize_t dims[1] = {1}, offset[1] = {static_cast<hsize_t>(row)}, offset_out[1] = {0}, count[1] = {1};
    #else
    hsize_t dims[1] = {1}, offset[1] = {row}, offset_out[1] = {0}, count[1] = {1};
    #endif
    hid_t dataset_id, dataspace, memspace, space, filetype, memtype;
    herr_t status;
    char** rdata;
    int ndims = 0;
    if( H5Lexists(info->file_, name, H5P_DEFAULT) == 0)
    {
        fprintf(stderr, "Error accessing to dataset %s in h5 file\n", name);
        return -1;
    }
    dataset_id = H5Dopen(info->file_, name);
    dataspace = H5Dget_space(dataset_id);
    filetype = H5Dget_type (dataset_id);
    space = H5Dget_space (dataset_id);
    ndims = H5Sget_simple_extent_dims (space, dims, NULL);
    rdata = (char **) malloc (sizeof (char *));
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, H5T_VARIABLE);
    H5Tset_cset(memtype, H5T_CSET_UTF8);
    status = H5Sselect_hyperslab( space, H5S_SELECT_SET, offset, NULL, count, NULL );
    memspace = H5Screate_simple( 1, dims, NULL );
    status = H5Sselect_hyperslab( memspace, H5S_SELECT_SET, offset_out, NULL, count, NULL );
    status = H5Dread (dataset_id, memtype, memspace, space, H5P_DEFAULT, rdata);
    hoc_assign_str( hoc_dest, rdata[0] );
    free(rdata);
    status = H5Sclose (space);
    status = H5Sclose (dataspace);
    status = H5Tclose (filetype);
    status = H5Tclose (memtype);
    status = H5Dclose (dataset_id);
    return 0;
}
#endif  // DISABLE_HDF5
#endif


namespace coreneuron {


    inline int getDimensions_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_getDimensions = 0;
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            if( info->file_ >= 0 && ifarg(1) && hoc_is_str_arg(1) ) {
                loadDimensions( info, gargstr(1) );
            }
        #endif
        #endif  // CORENEURON_BUILD

        return ret_getDimensions;
    }


    inline int getDataString_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        int ret_getDataString = 0;
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            if( info->file_ >= 0 && ifarg(1) && hoc_is_str_arg(1) && ifarg(2) && ifarg(3) )
            {
                loadDataString( info, gargstr(1), *getarg(2), hoc_pgargstr(3) );
            }
        #endif
        #endif  // CORENEURON_BUILD

        return ret_getDataString;
    }


    inline double redirect_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_redirect = 0.0;
         {
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
        #ifndef DISABLE_MPI
            FILE *fout;
            char fname[128];
            int mpi_size, mpi_rank;
            MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
            MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
            if( mpi_rank != 0 ) {
                sprintf( fname, "NodeFiles/%d.%dnode.out", mpi_rank, mpi_size );
                fout = freopen( fname, "w", stdout );
                if( !fout ) {
                    fprintf( stderr, "failed to redirect.  Terminating\n" );
                    exit(0);
                }
                sprintf( fname, "NodeFiles/%d.%dnode.err", mpi_rank, mpi_size );
                fout = freopen( fname, "w", stderr );
                setbuf( fout, NULL );
            }
        #endif  // DISABLE_MPI
        #endif  // DISABLE_HDF5
        #endif  // CORENEURON_BUILD
        }

        return ret_redirect;
    }


    inline double checkVersion_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_checkVersion = 0.0;
         {
        #ifndef CORENEURON_BUILD
            int versionNumber = 0;
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            int mpi_size=1, mpi_rank=0;
        #ifndef DISABLE_MPI
            MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
            MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
        #endif
            if( mpi_rank == 0 )
            {
                if( H5Lexists(info->file_, "info", H5P_DEFAULT) == 0) {
                    versionNumber = 0;
                } else {
                    hid_t dataset_id = H5Dopen( info->file_, "info" );
                    hid_t attr_id = H5Aopen_name( dataset_id, "version" );
                    H5Aread( attr_id, H5T_NATIVE_INT, &versionNumber );
                    H5Aclose(attr_id);
                    H5Dclose(dataset_id);
                }
            }
        #ifndef DISABLE_MPI
            MPI_Bcast( &versionNumber, 1, MPI_INT, 0, MPI_COMM_WORLD );
        #endif
        #endif  // DISABLE_HDF5
            return versionNumber;
        #endif  // CORENEURON_BUILD
        }

        return ret_checkVersion;
    }


    inline double loadData_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_loadData = 0.0;
         {
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            if(info->file_>=0 && ifarg(1) && hoc_is_str_arg(1))
            {
                if( ifarg(2) ) {
                    if( *getarg(2) == 1 ) { 
                        info->mode = LONG_VECTOR;
                        return loadDataVector( info, gargstr(1) );
                    }
                }
                info->mode = FLOAT_MATRIX;
                return loadDataMatrix( info, gargstr(1) );
            }
            else if( ifarg(1) )
            {
                int gid = *getarg(1);
                int gidIndex=0;
                for( gidIndex=0; gidIndex<info->confirmedCells.count; gidIndex++ ) {
                    if( info->confirmedCells.gids[gidIndex] == gid ) {
                        openFile( info, info->synapseCatalog.rootName, info->confirmedCells.fileIDs[gidIndex], 0, 0, 0 );
                        char cellname[256];
                        sprintf( cellname, "a%d", gid );
                        info->mode = FLOAT_MATRIX;
                        return loadDataMatrix( info, cellname );
                    }
                }
                if( info->verboseLevel > 0 )
                    fprintf( stderr, "Warning: failed to find data for gid %d\n", gid );
            }
            return -1;
        #endif  // DISABLE_HDF5
        #endif  // CORENEURON_BUILD
        }

        return ret_loadData;
    }


    inline double getNoOfColumns_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_getNoOfColumns = 0.0;
         {
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            if(info->file_>=0 && ifarg(1) && hoc_is_str_arg(1))
            {
                char name[256];
                strncpy(name,gargstr(1),256);
                if(strncmp(info->name_group,name,256) == 0)
                {
                    int res = (int) info->columnsize_;
                    return res;
                }
                return 0;
            }
            else
            {
                return 0;
            }
        #endif
        #endif  // CORENEURON_BUILD
        }

        return ret_getNoOfColumns;
    }


    inline double numberofrows_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_numberofrows = 0.0;
         {
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            if(info->file_>=0 && ifarg(1) && hoc_is_str_arg(1))
            {
                char name[256];
                strncpy(name,gargstr(1),256);
                if(strncmp(info->name_group,name,256) == 0)
                {
                    int res = (int) info->rowsize_;
                    return res;
                }
                return 0;
            }
            else
            {
                fprintf(stderr, "general error: bad file handle, missing arg?");
                return 0;
            }
        #endif
        #endif  // CORENEURON_BUILD
        }

        return ret_numberofrows;
    }


    inline double getData_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_getData = 0.0;
         {
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            if(info->file_>=0&& ifarg(1) && hoc_is_str_arg(1) && ifarg(2) && ifarg(3))
            {
                char name[256];
                strncpy(name,gargstr(1),256);
                if(strncmp(info->name_group,name,256) == 0)
                {
                    hsize_t row,column;
                    row = (hsize_t) *getarg(2);
                    column = (hsize_t) *getarg(3);
                    if(row<0 || row >=info->rowsize_ || column < 0 || column>=info->columnsize_)
                    {
                        fprintf(stderr, "ERROR: trying to access to a row and column erroneus on %s, size: %lld,%lld accessing to %lld,%lld\n",
                                name, info->rowsize_, info->columnsize_, row, column);
                        return 0;
                    }
                    if( info->mode == FLOAT_MATRIX ) {
                        return info->datamatrix_[row*info->columnsize_ + column];
                    } else if( info->mode == LONG_VECTOR ) {
                        return (double) info->datavector_[row];
                    } else {
                        fprintf( stderr, "unexpected mode: %d\n", info->mode );
                    }
                }
                fprintf(stderr, "(Getting data)Error on the name of last loaded data: access:%s loaded:%s\n",name,info->name_group);
                return 0;
            }
            else
            {
                fprintf( stderr, "ERROR:Error on number of rows of %s\n", gargstr(1) );
                return 0;
            }
        #endif
        #endif  // CORENEURON_BUILD
        }

        return ret_getData;
    }


    inline double getDataInt_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_getDataInt = 0.0;
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            int value = -1;
            ret_getDataInt = 0;
            if( info->file_ >= 0 && ifarg(1) && hoc_is_str_arg(1) && ifarg(2) )
            {
                if( loadDataInt( info, gargstr(1), *getarg(2), &value ) == 0 ) {
                    ret_getDataInt = value;
                }
            }
        #endif
        #endif  // CORENEURON_BUILD

        return ret_getDataInt;
    }


    inline double getColumnDataRange_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_getColumnDataRange = 0.0;
         {
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            IvocVect* pdVec = NULL;
            double* pd  = NULL;
            int i = 0;
            int nStart, nEnd, count;
            if(info->file_>=0&& ifarg(1) && hoc_is_str_arg(1) && ifarg(2) )
            {
                char name[256];
                strncpy(name,gargstr(1),256);
                if(strncmp(info->name_group,name,256) == 0)
                {
                    hsize_t column;
                    column  = (hsize_t) *getarg(2);
                    if(column<0 || column >=info->columnsize_ )
                    {
                        fprintf(stderr, "ERROR: trying to access to a column erroneus on %s, size: %lld,%lld accessing to column %lld\n ",name,info->rowsize_,info->columnsize_,column);
                        return 0;
                    }
                    pdVec = vector_arg(3);
                    nStart = (int)*getarg(4);
                    nEnd  = (int)*getarg(5);
                    vector_resize(pdVec, nEnd-nStart+1);
                    pd = vector_vec(pdVec);
                    count =0;
                    for( i=nStart; i<=nEnd; i++){
                        pd[count] = info->datamatrix_[i*info->columnsize_ + column];
                        count = count +1;
                    }
                    return 1;
                }
                fprintf(stderr, "(Getting data)Error on the name of last loaded data: access:%s loaded:%s\n",name,info->name_group);
                return 0;
            }
            else
            {
                return 0;
            }
        #endif
        #endif  // CORENEURON_BUILD
        }

        return ret_getColumnDataRange;
    }


    inline double getColumnData_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_getColumnData = 0.0;
         {
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            IvocVect* pdVec = NULL;
            double* pd  = NULL;
            int i = 0;
            if(info->file_>=0&& ifarg(1) && hoc_is_str_arg(1) && ifarg(2) )
            {
                char name[256];
                strncpy(name,gargstr(1),256);
                if(strncmp(info->name_group,name,256) == 0)
                {
                    hsize_t column;
                    column  = (hsize_t) *getarg(2);
                    if(column<0 || column >=info->columnsize_ )
                    {
                        fprintf(stderr, "ERROR: trying to access to a column erroneus on %s, size: %lld,%lld accessing to column %lld\n ",name,info->rowsize_,info->columnsize_,column);
                        return 0;
                    }
                    pdVec = vector_arg(3);
                    vector_resize(pdVec, (int) info->rowsize_);
                    pd = vector_vec(pdVec);
                    for( i=0; i<info->rowsize_; i++){
                        pd[i] = info->datamatrix_[i*info->columnsize_ + column];
                    }
                    return 1;
                }
                fprintf(stderr, "(Getting data)Error on the name of last loaded data: access:%s loaded:%s\n",name,info->name_group);
                return 0;
            }
            else
            {
                return 0;
            }
        #endif
        #endif  // CORENEURON_BUILD
        }

        return ret_getColumnData;
    }


    inline double getAttributeValue_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_getAttributeValue = 0.0;
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            if( info->file_ >= 0 && ifarg(1) && hoc_is_str_arg(1) && ifarg(2) && hoc_is_str_arg(2) )
            {
                if( H5Lexists(info->file_, gargstr(1), H5P_DEFAULT) == 0)
                {
                    fprintf( stderr, "Error: no dataset with name %s available.\n", gargstr(1) );
                    return 0;
                }
                hid_t dataset_id = H5Dopen( info->file_, gargstr(1) );
                double soughtValue;
                hid_t attr_id = H5Aopen_name( dataset_id, gargstr(2) );
                if( attr_id < 0 ) {
                    fprintf( stderr, "Error: failed to open attribute %s\n", gargstr(2) );
                    return 0;
                }
                H5Aread( attr_id, H5T_NATIVE_DOUBLE, &soughtValue );
                H5Dclose(dataset_id);
                return soughtValue;
            }
            return 0;
        #endif
        #endif  // CORENEURON_BUILD

        return ret_getAttributeValue;
    }


    inline double closeFile_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_closeFile = 0.0;
         {
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
            INFOCAST;
            Info* info = *ip;
            if(info->file_ >=0)
            {
                H5Fclose(info->file_);
                info->file_ = -1;
            }
            if(info->datamatrix_ != NULL)
            {
                free(info->datamatrix_);
                info->datamatrix_ = NULL;
            }
        #endif
        #endif  // CORENEURON_BUILD
        }

        return ret_closeFile;
    }


    inline double exchangeSynapseLocations_HDF5Reader(int id, int pnodecount, HDF5Reader_Instance* inst, double* data, const Datum* indexes, ThreadDatum* thread, NrnThread* nt, double v) {
        double ret_exchangeSynapseLocations = 0.0;
        #ifndef CORENEURON_BUILD
        #ifndef DISABLE_HDF5
        #ifndef DISABLE_MPI
            INFOCAST;
            Info* info = *ip;
            IvocVect* vv = vector_arg(1);
            int gidCount = vector_capacity(vv);
            double *gidList = vector_vec(vv);
            if( info->synapseCatalog.directFile != NULL ) {
                readDirectMapping( info, gidCount, gidList );
                return 0;
            }
            int mpi_size, mpi_rank;
            MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
            MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
            int *foundCountsAcrossCPUs = (int*) malloc( sizeof(int)*mpi_size );
            int *foundDispls = (int*) malloc( sizeof(int)*mpi_size );
            int bufferSize;
            MPI_Allreduce( &gidCount, &bufferSize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
            double* gidsRequested = (double*) malloc ( sizeof(double)*bufferSize );
            int* gidsFound = (int*) malloc( sizeof(int)*bufferSize );
            int* fileIDsFound = (int*) malloc( sizeof(int)*bufferSize );
            double *tempRequest;  
            int activeRank, requestCount;
            for( activeRank=0; activeRank<mpi_size; activeRank++ ) {
                if( activeRank == mpi_rank ) {
                    requestCount = gidCount;
                    tempRequest = gidsRequested;
                    gidsRequested = gidList;
                }
                MPI_Bcast( &requestCount, 1, MPI_INT, activeRank, MPI_COMM_WORLD );
                MPI_Bcast( gidsRequested, requestCount, MPI_DOUBLE, activeRank, MPI_COMM_WORLD );
                int nFiles = 0;
                int fileIndex, gidIndex, requestIndex;
                for( fileIndex=0; fileIndex < info->synapseCatalog.nFiles; fileIndex++ ) {
                    for( gidIndex=0; gidIndex < info->synapseCatalog.availablegidCount[fileIndex]; gidIndex++ ) {
                        for( requestIndex=0; requestIndex < requestCount; requestIndex++ ) {
                            if( info->synapseCatalog.availablegids[fileIndex][gidIndex] == gidsRequested[requestIndex] ) {
                                gidsFound[nFiles] = gidsRequested[requestIndex];
                                fileIDsFound[nFiles++] = info->synapseCatalog.fileIDs[fileIndex];
                            }
                        }
                    }
                }
                MPI_Gather( &nFiles, 1, MPI_INT, foundCountsAcrossCPUs, 1, MPI_INT, activeRank, MPI_COMM_WORLD );
                if( activeRank == mpi_rank ) {
                    info->confirmedCells.count = 0;
                    int nodeIndex;
                    for( nodeIndex=0; nodeIndex<mpi_size; nodeIndex++ ) {
                        foundDispls[nodeIndex] = info->confirmedCells.count;
                        info->confirmedCells.count += foundCountsAcrossCPUs[nodeIndex];
                    }
                    info->confirmedCells.gids = (int*) malloc ( sizeof(int)*info->confirmedCells.count );
                    info->confirmedCells.fileIDs = (int*) malloc ( sizeof(int)*info->confirmedCells.count );
                }
                MPI_Gatherv( gidsFound, nFiles, MPI_INT, info->confirmedCells.gids, foundCountsAcrossCPUs, foundDispls, MPI_INT, activeRank, MPI_COMM_WORLD );
                MPI_Gatherv( fileIDsFound, nFiles, MPI_INT, info->confirmedCells.fileIDs, foundCountsAcrossCPUs, foundDispls, MPI_INT, activeRank, MPI_COMM_WORLD );
                if( activeRank == mpi_rank ) {
                    gidsRequested = tempRequest;
                }
            }
            free(gidsRequested);
            free(gidsFound);
            free(fileIDsFound);
            free(foundCountsAcrossCPUs);
            free(foundDispls);
        #endif  // DISABLE_MPI
        #endif  // DISABLE_HDF5
        #endif  // CORENEURON_BUILD

        return ret_exchangeSynapseLocations;
    }


    static inline void net_receive_HDF5Reader(Point_process* pnt, int weight_index, double flag) {
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
        auto* const inst = static_cast<HDF5Reader_Instance*>(ml->instance);

        double t = nt->_t;
        inst->tsave[id] = t;
        {
        }
    }


    /** initialize channel */
    void nrn_init_HDF5Reader(NrnThread* nt, Memb_list* ml, int type) {
        int nodecount = ml->nodecount;
        int pnodecount = ml->_nodecount_padded;
        const int* node_index = ml->nodeindices;
        double* data = ml->data;
        const double* voltage = nt->_actual_v;
        Datum* indexes = ml->pdata;
        ThreadDatum* thread = ml->_thread;

        setup_instance(nt, ml);
        auto* const inst = static_cast<HDF5Reader_Instance*>(ml->instance);

        if (_nrn_skip_initmodel == 0) {
            #pragma omp simd
            #pragma ivdep
            for (int id = 0; id < nodecount; id++) {
                inst->tsave[id] = -1e20;
                double v = 0.0;
            }
        }
    }


    /** register channel with the simulator */
    void _HDF5reader_reg() {

        int mech_type = nrn_get_mechtype("HDF5Reader");
        HDF5Reader_global.mech_type = mech_type;
        if (mech_type == -1) {
            return;
        }

        _nrn_layout_reg(mech_type, 0);
        point_register_mech(mechanism, nrn_alloc_HDF5Reader, nullptr, nullptr, nullptr, nrn_init_HDF5Reader, nrn_private_constructor_HDF5Reader, nrn_private_destructor_HDF5Reader, first_pointer_var_index(), nrn_constructor_HDF5Reader, nrn_destructor_HDF5Reader, 1);

        hoc_register_prop_size(mech_type, float_variables_size(), int_variables_size());
        hoc_register_dparam_semantics(mech_type, 0, "area");
        hoc_register_dparam_semantics(mech_type, 1, "pntproc");
        hoc_register_dparam_semantics(mech_type, 2, "pointer");
        add_nrn_artcell(mech_type, 0);
        set_pnt_receive(mech_type, net_receive_HDF5Reader, nullptr, num_net_receive_args());
        hoc_register_var(hoc_scalar_double, hoc_vector_double, NULL);
    }
}
