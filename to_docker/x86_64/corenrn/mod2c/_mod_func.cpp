#include <cstdio>
namespace coreneuron {
extern int nrnmpi_myid;
extern int nrn_nobanner_;
extern int _ALU_reg(void),
  _CoreNEURONArtificialCell_reg(void),
  _DetAMPANMDA_reg(void),
  _DetGABAAB_reg(void),
  _GluSynapse_reg(void),
  _HDF5reader_reg(void),
  _Im_ms_reg(void),
  _NO_reg(void),
  _ProbAMPANMDA_EMS_reg(void),
  _ProbGABAAB_EMS_reg(void),
  _SonataReportHelper_reg(void),
  _SonataReports_reg(void),
  _SpikeWriter_reg(void),
  _SynapseReader_reg(void),
  _TTXDynamicsSwitch_reg(void),
  _VecStim_reg(void),
  _bk_ch_reg(void),
  _bk_fs_reg(void),
  _bk_ms_reg(void),
  _ca_ch_reg(void),
  _cacumm_reg(void),
  _cacummb_reg(void),
  _cadyn_fs_reg(void),
  _cadyn_ms_reg(void),
  _cagk_reg(void),
  _cal12_ms_reg(void),
  _cal13_ms_reg(void),
  _cal2_reg(void),
  _cal_ch_reg(void),
  _caldyn_ms_reg(void),
  _can2_reg(void),
  _can_fs_reg(void),
  _can_ms_reg(void),
  _cap_ch_reg(void),
  _caq_fs_reg(void),
  _caq_ms_reg(void),
  _car_fs_reg(void),
  _car_ms_reg(void),
  _cat_reg(void),
  _cat32_ms_reg(void),
  _cat33_ms_reg(void),
  _exp2syn_reg(void),
  _expsyn_reg(void),
  _gap_reg(void),
  _h_reg(void),
  _hcn12_ch_reg(void),
  _hd_lts_reg(void),
  _hh_reg(void),
  _im_lts_reg(void),
  _it_lts_reg(void),
  _kadist_reg(void),
  _kaf_fs_reg(void),
  _kaf_lts_reg(void),
  _kaf_ms_reg(void),
  _kaprox_reg(void),
  _kas_fs_reg(void),
  _kas_ms_reg(void),
  _kca_reg(void),
  _kcnq_ch_reg(void),
  _kd_reg(void),
  _kd2_reg(void),
  _kdb_reg(void),
  _kdb_lts_reg(void),
  _kdr_fs_reg(void),
  _kdr_lts_reg(void),
  _kdr_ms_reg(void),
  _kdrb_lts_reg(void),
  _kdrbca1_reg(void),
  _kdrca1_reg(void),
  _kir23_ch_reg(void),
  _kir23_lts_reg(void),
  _kir2_ch_reg(void),
  _kir_fs_reg(void),
  _kir_ms_reg(void),
  _km_reg(void),
  _kmb_reg(void),
  _kv2_ch_reg(void),
  _kv4_ch_reg(void),
  _lookupTableV2_reg(void),
  _na2_ch_reg(void),
  _na3_lts_reg(void),
  _na3n_reg(void),
  _na_ch_reg(void),
  _naf_fs_reg(void),
  _naf_lts_reg(void),
  _naf_ms_reg(void),
  _naxn_reg(void),
  _netstim_reg(void),
  _netstim_inhpoisson_reg(void),
  _par_ggap_reg(void),
  _passive_reg(void),
  _pattern_reg(void),
  _sk_ch_reg(void),
  _sk_fs_reg(void),
  _sk_ms_reg(void),
  _stim_reg(void),
  _svclmp_reg(void),
  _tmgabaa_reg(void),
  _tmglut_reg(void),
  _tmglut_M1RH_D1_reg(void),
  _tmglut_double_reg(void),
  _vecevent_reg(void);

void modl_reg() {
    if (!nrn_nobanner_ && nrnmpi_myid < 1) {
        fprintf(stderr, " Additional mechanisms from files\n");
        fprintf(stderr, " ALU.mod");
        fprintf(stderr, " CoreNEURONArtificialCell.mod");
        fprintf(stderr, " DetAMPANMDA.mod");
        fprintf(stderr, " DetGABAAB.mod");
        fprintf(stderr, " GluSynapse.mod");
        fprintf(stderr, " HDF5reader.mod");
        fprintf(stderr, " Im_ms.mod");
        fprintf(stderr, " NO.mod");
        fprintf(stderr, " ProbAMPANMDA_EMS.mod");
        fprintf(stderr, " ProbGABAAB_EMS.mod");
        fprintf(stderr, " SonataReportHelper.mod");
        fprintf(stderr, " SonataReports.mod");
        fprintf(stderr, " SpikeWriter.mod");
        fprintf(stderr, " SynapseReader.mod");
        fprintf(stderr, " TTXDynamicsSwitch.mod");
        fprintf(stderr, " VecStim.mod");
        fprintf(stderr, " bk_ch.mod");
        fprintf(stderr, " bk_fs.mod");
        fprintf(stderr, " bk_ms.mod");
        fprintf(stderr, " ca_ch.mod");
        fprintf(stderr, " cacumm.mod");
        fprintf(stderr, " cacummb.mod");
        fprintf(stderr, " cadyn_fs.mod");
        fprintf(stderr, " cadyn_ms.mod");
        fprintf(stderr, " cagk.mod");
        fprintf(stderr, " cal12_ms.mod");
        fprintf(stderr, " cal13_ms.mod");
        fprintf(stderr, " cal2.mod");
        fprintf(stderr, " cal_ch.mod");
        fprintf(stderr, " caldyn_ms.mod");
        fprintf(stderr, " can2.mod");
        fprintf(stderr, " can_fs.mod");
        fprintf(stderr, " can_ms.mod");
        fprintf(stderr, " cap_ch.mod");
        fprintf(stderr, " caq_fs.mod");
        fprintf(stderr, " caq_ms.mod");
        fprintf(stderr, " car_fs.mod");
        fprintf(stderr, " car_ms.mod");
        fprintf(stderr, " cat.mod");
        fprintf(stderr, " cat32_ms.mod");
        fprintf(stderr, " cat33_ms.mod");
        fprintf(stderr, " exp2syn.mod");
        fprintf(stderr, " expsyn.mod");
        fprintf(stderr, " gap.mod");
        fprintf(stderr, " h.mod");
        fprintf(stderr, " hcn12_ch.mod");
        fprintf(stderr, " hd_lts.mod");
        fprintf(stderr, " hh.mod");
        fprintf(stderr, " im_lts.mod");
        fprintf(stderr, " it_lts.mod");
        fprintf(stderr, " kadist.mod");
        fprintf(stderr, " kaf_fs.mod");
        fprintf(stderr, " kaf_lts.mod");
        fprintf(stderr, " kaf_ms.mod");
        fprintf(stderr, " kaprox.mod");
        fprintf(stderr, " kas_fs.mod");
        fprintf(stderr, " kas_ms.mod");
        fprintf(stderr, " kca.mod");
        fprintf(stderr, " kcnq_ch.mod");
        fprintf(stderr, " kd.mod");
        fprintf(stderr, " kd2.mod");
        fprintf(stderr, " kdb.mod");
        fprintf(stderr, " kdb_lts.mod");
        fprintf(stderr, " kdr_fs.mod");
        fprintf(stderr, " kdr_lts.mod");
        fprintf(stderr, " kdr_ms.mod");
        fprintf(stderr, " kdrb_lts.mod");
        fprintf(stderr, " kdrbca1.mod");
        fprintf(stderr, " kdrca1.mod");
        fprintf(stderr, " kir23_ch.mod");
        fprintf(stderr, " kir23_lts.mod");
        fprintf(stderr, " kir2_ch.mod");
        fprintf(stderr, " kir_fs.mod");
        fprintf(stderr, " kir_ms.mod");
        fprintf(stderr, " km.mod");
        fprintf(stderr, " kmb.mod");
        fprintf(stderr, " kv2_ch.mod");
        fprintf(stderr, " kv4_ch.mod");
        fprintf(stderr, " lookupTableV2.mod");
        fprintf(stderr, " na2_ch.mod");
        fprintf(stderr, " na3_lts.mod");
        fprintf(stderr, " na3n.mod");
        fprintf(stderr, " na_ch.mod");
        fprintf(stderr, " naf_fs.mod");
        fprintf(stderr, " naf_lts.mod");
        fprintf(stderr, " naf_ms.mod");
        fprintf(stderr, " naxn.mod");
        fprintf(stderr, " netstim.mod");
        fprintf(stderr, " netstim_inhpoisson.mod");
        fprintf(stderr, " par_ggap.mod");
        fprintf(stderr, " passive.mod");
        fprintf(stderr, " pattern.mod");
        fprintf(stderr, " sk_ch.mod");
        fprintf(stderr, " sk_fs.mod");
        fprintf(stderr, " sk_ms.mod");
        fprintf(stderr, " stim.mod");
        fprintf(stderr, " svclmp.mod");
        fprintf(stderr, " tmgabaa.mod");
        fprintf(stderr, " tmglut.mod");
        fprintf(stderr, " tmglut_M1RH_D1.mod");
        fprintf(stderr, " tmglut_double.mod");
        fprintf(stderr, " vecevent.mod");
        fprintf(stderr, "\n\n");
    }

    _ALU_reg();
    _CoreNEURONArtificialCell_reg();
    _DetAMPANMDA_reg();
    _DetGABAAB_reg();
    _GluSynapse_reg();
    _HDF5reader_reg();
    _Im_ms_reg();
    _NO_reg();
    _ProbAMPANMDA_EMS_reg();
    _ProbGABAAB_EMS_reg();
    _SonataReportHelper_reg();
    _SonataReports_reg();
    _SpikeWriter_reg();
    _SynapseReader_reg();
    _TTXDynamicsSwitch_reg();
    _VecStim_reg();
    _bk_ch_reg();
    _bk_fs_reg();
    _bk_ms_reg();
    _ca_ch_reg();
    _cacumm_reg();
    _cacummb_reg();
    _cadyn_fs_reg();
    _cadyn_ms_reg();
    _cagk_reg();
    _cal12_ms_reg();
    _cal13_ms_reg();
    _cal2_reg();
    _cal_ch_reg();
    _caldyn_ms_reg();
    _can2_reg();
    _can_fs_reg();
    _can_ms_reg();
    _cap_ch_reg();
    _caq_fs_reg();
    _caq_ms_reg();
    _car_fs_reg();
    _car_ms_reg();
    _cat_reg();
    _cat32_ms_reg();
    _cat33_ms_reg();
    _exp2syn_reg();
    _expsyn_reg();
    _gap_reg();
    _h_reg();
    _hcn12_ch_reg();
    _hd_lts_reg();
    _hh_reg();
    _im_lts_reg();
    _it_lts_reg();
    _kadist_reg();
    _kaf_fs_reg();
    _kaf_lts_reg();
    _kaf_ms_reg();
    _kaprox_reg();
    _kas_fs_reg();
    _kas_ms_reg();
    _kca_reg();
    _kcnq_ch_reg();
    _kd_reg();
    _kd2_reg();
    _kdb_reg();
    _kdb_lts_reg();
    _kdr_fs_reg();
    _kdr_lts_reg();
    _kdr_ms_reg();
    _kdrb_lts_reg();
    _kdrbca1_reg();
    _kdrca1_reg();
    _kir23_ch_reg();
    _kir23_lts_reg();
    _kir2_ch_reg();
    _kir_fs_reg();
    _kir_ms_reg();
    _km_reg();
    _kmb_reg();
    _kv2_ch_reg();
    _kv4_ch_reg();
    _lookupTableV2_reg();
    _na2_ch_reg();
    _na3_lts_reg();
    _na3n_reg();
    _na_ch_reg();
    _naf_fs_reg();
    _naf_lts_reg();
    _naf_ms_reg();
    _naxn_reg();
    _netstim_reg();
    _netstim_inhpoisson_reg();
    _par_ggap_reg();
    _passive_reg();
    _pattern_reg();
    _sk_ch_reg();
    _sk_fs_reg();
    _sk_ms_reg();
    _stim_reg();
    _svclmp_reg();
    _tmgabaa_reg();
    _tmglut_reg();
    _tmglut_M1RH_D1_reg();
    _tmglut_double_reg();
    _vecevent_reg();
}
} //namespace coreneuron
