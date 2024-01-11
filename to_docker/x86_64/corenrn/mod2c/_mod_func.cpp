#include <cstdio>
namespace coreneuron {
extern int nrnmpi_myid;
extern int nrn_nobanner_;
extern int _ALU_reg(void),
  _CoreNEURONArtificialCell_reg(void),
  _GRC_CA_reg(void),
  _GRC_KM_reg(void),
  _Gfluct_reg(void),
  _GluSynapse_reg(void),
  _HDF5reader_reg(void),
  _Kca11_reg(void),
  _Kca22_reg(void),
  _Kir23_reg(void),
  _Kv11_reg(void),
  _Kv15_reg(void),
  _Kv22_reg(void),
  _Kv34_reg(void),
  _Kv43_reg(void),
  _Leak_reg(void),
  _ProbAMPANMDA_EMS_reg(void),
  _ProbGABAAB_EMS_reg(void),
  _SonataReportHelper_reg(void),
  _SonataReports_reg(void),
  _SpikeWriter_reg(void),
  _SynapseReader_reg(void),
  _ThreshDetect_reg(void),
  _VecStim_reg(void),
  _ampanmda_reg(void),
  _distrt_reg(void),
  _exp2syn_reg(void),
  _expsyn_reg(void),
  _fi_reg(void),
  _fi_stdp_reg(void),
  _gap_reg(void),
  _hh_reg(void),
  _kamt_reg(void),
  _kdrmt_reg(void),
  _ks_reg(void),
  _lookupTableV2_reg(void),
  _naxn_reg(void),
  _netstim_reg(void),
  _netstim_inhpoisson_reg(void),
  _orn_reg(void),
  _ostimhelper_reg(void),
  _passive_reg(void),
  _pattern_reg(void),
  _stim_reg(void),
  _svclmp_reg(void);

void modl_reg() {
    if (!nrn_nobanner_ && nrnmpi_myid < 1) {
        fprintf(stderr, " Additional mechanisms from files\n");
        fprintf(stderr, " ALU.mod");
        fprintf(stderr, " CoreNEURONArtificialCell.mod");
        fprintf(stderr, " GRC_CA.mod");
        fprintf(stderr, " GRC_KM.mod");
        fprintf(stderr, " Gfluct.mod");
        fprintf(stderr, " GluSynapse.mod");
        fprintf(stderr, " HDF5reader.mod");
        fprintf(stderr, " Kca11.mod");
        fprintf(stderr, " Kca22.mod");
        fprintf(stderr, " Kir23.mod");
        fprintf(stderr, " Kv11.mod");
        fprintf(stderr, " Kv15.mod");
        fprintf(stderr, " Kv22.mod");
        fprintf(stderr, " Kv34.mod");
        fprintf(stderr, " Kv43.mod");
        fprintf(stderr, " Leak.mod");
        fprintf(stderr, " ProbAMPANMDA_EMS.mod");
        fprintf(stderr, " ProbGABAAB_EMS.mod");
        fprintf(stderr, " SonataReportHelper.mod");
        fprintf(stderr, " SonataReports.mod");
        fprintf(stderr, " SpikeWriter.mod");
        fprintf(stderr, " SynapseReader.mod");
        fprintf(stderr, " ThreshDetect.mod");
        fprintf(stderr, " VecStim.mod");
        fprintf(stderr, " ampanmda.mod");
        fprintf(stderr, " distrt.mod");
        fprintf(stderr, " exp2syn.mod");
        fprintf(stderr, " expsyn.mod");
        fprintf(stderr, " fi.mod");
        fprintf(stderr, " fi_stdp.mod");
        fprintf(stderr, " gap.mod");
        fprintf(stderr, " hh.mod");
        fprintf(stderr, " kamt.mod");
        fprintf(stderr, " kdrmt.mod");
        fprintf(stderr, " ks.mod");
        fprintf(stderr, " lookupTableV2.mod");
        fprintf(stderr, " naxn.mod");
        fprintf(stderr, " netstim.mod");
        fprintf(stderr, " netstim_inhpoisson.mod");
        fprintf(stderr, " orn.mod");
        fprintf(stderr, " ostimhelper.mod");
        fprintf(stderr, " passive.mod");
        fprintf(stderr, " pattern.mod");
        fprintf(stderr, " stim.mod");
        fprintf(stderr, " svclmp.mod");
        fprintf(stderr, "\n\n");
    }

    _ALU_reg();
    _CoreNEURONArtificialCell_reg();
    _GRC_CA_reg();
    _GRC_KM_reg();
    _Gfluct_reg();
    _GluSynapse_reg();
    _HDF5reader_reg();
    _Kca11_reg();
    _Kca22_reg();
    _Kir23_reg();
    _Kv11_reg();
    _Kv15_reg();
    _Kv22_reg();
    _Kv34_reg();
    _Kv43_reg();
    _Leak_reg();
    _ProbAMPANMDA_EMS_reg();
    _ProbGABAAB_EMS_reg();
    _SonataReportHelper_reg();
    _SonataReports_reg();
    _SpikeWriter_reg();
    _SynapseReader_reg();
    _ThreshDetect_reg();
    _VecStim_reg();
    _ampanmda_reg();
    _distrt_reg();
    _exp2syn_reg();
    _expsyn_reg();
    _fi_reg();
    _fi_stdp_reg();
    _gap_reg();
    _hh_reg();
    _kamt_reg();
    _kdrmt_reg();
    _ks_reg();
    _lookupTableV2_reg();
    _naxn_reg();
    _netstim_reg();
    _netstim_inhpoisson_reg();
    _orn_reg();
    _ostimhelper_reg();
    _passive_reg();
    _pattern_reg();
    _stim_reg();
    _svclmp_reg();
}
} //namespace coreneuron
