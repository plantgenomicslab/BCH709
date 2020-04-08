---
layout: page
title: Frequently Asked Questions
---

Thank you for your interest in BCH709. Below you will find answers to some frequently asked questions about this curriculum. If the answer to your question doesn’t appear, please contact [wyim@unr.edu](mailto:wyim@unr.edu). 
In 2019 Fall, this is first time of BCH709. We don't have that many FAQ.


### What does this workshop cover? 

This class teaches data management and analysis for genomics research including: best practices for organization of bioinformatics projects and data, use of command line utilities, use of command line tools to analyze sequence quality and perform variant calling, and connecting to and using cloud computing. 

### What experience do learners need to have before this workshop? What will they be able to do by the end of the workshop? 

Learners are expected to have some familiarity with biological concepts, including the concept of genomic variation within a population. By the end of the workshop, learners will be able to: 

- structure their metadata, organize and document their genomics data and bioinformatics workflow, and access data on the NCBI sequence read archive (SRA) database,
- navigate their file systems, create, copy, move, and remove files and directories, and automate repetitive tasks using scripts and wildcards,
- use command-line tools to perform quality control, align reads to a reference genome, and identify and visualize between-sample variation,
- work with Amazon AWS cloud computing and transfer data between their local computer and cloud resources.

### What are the software, hardware, and connectivity needs for this workshop?
Learners will need to bring a laptop (not a tablet) with any spreadsheet program installed (e.g. LibreOffice, Microsoft Excel). Learners using laptop need to check setup pages to prepare the computer. There are no other hardware or software requirements. Learners will need a stable, strong internet connection in order to work on the remote computing system used for this workshop.

### My lab has its own compute cluster, or our research group uses a different cloud computing resource. Can we use that system?
Yes you can use what ever you want but need to know what it is.

### What experience do learner need to have for this class?
Anyone who has some experience using the Bash shell can be an effective, but we provide free prerequisite course for any learners that doesn't have any experience.



### What is this error ?

```bash
Traceback (most recent call last):
  File "/data/gpfs/home/wyim/scratch/bin/miniconda3/envs/preprocessing/bin/multiqc", line 6, in <module>
    from multiqc.__main__ import multiqc
  File "/data/gpfs/home/wyim/scratch/bin/miniconda3/envs/preprocessing/lib/python3.6/site-packages/multiqc/__main__.py", line 44, in <module>
    multiqc.run_cli(prog_name='multiqc')
  File "/data/gpfs/home/wyim/scratch/bin/miniconda3/envs/preprocessing/lib/python3.6/site-packages/click/core.py", line 829, in __call__
    return self.main(*args, **kwargs)
  File "/data/gpfs/home/wyim/scratch/bin/miniconda3/envs/preprocessing/lib/python3.6/site-packages/click/core.py", line 760, in main
    _verify_python3_env()
  File "/data/gpfs/home/wyim/scratch/bin/miniconda3/envs/preprocessing/lib/python3.6/site-packages/click/_unicodefun.py", line 130, in _verify_python3_env
    " mitigation steps.{}".format(extra)
RuntimeError: Click will abort further execution because Python 3 was configured to use ASCII as encoding for the environment. Consult https://click.palletsprojects.com/python3/ for mitigation steps.

This system lists a couple of UTF-8 supporting locales that you can pick from. The following suitable locales were discovered: aa_DJ.utf8, aa_ER.utf8, aa_ET.utf8, af_ZA.utf8, am_ET.utf8, an_ES.utf8, ar_AE.utf8, ar_BH.utf8, ar_DZ.utf8, ar_EG.utf8, ar_IN.utf8, ar_IQ.utf8, ar_JO.utf8, ar_KW.utf8, ar_LB.utf8, ar_LY.utf8, ar_MA.utf8, ar_OM.utf8, ar_QA.utf8, ar_SA.utf8, ar_SD.utf8, ar_SY.utf8, ar_TN.utf8, ar_YE.utf8, as_IN.utf8, ast_ES.utf8, ayc_PE.utf8, az_AZ.utf8, be_BY.utf8, bem_ZM.utf8, ber_DZ.utf8, ber_MA.utf8, bg_BG.utf8, bho_IN.utf8, bn_BD.utf8, bn_IN.utf8, bo_CN.utf8, bo_IN.utf8, br_FR.utf8, brx_IN.utf8, bs_BA.utf8, byn_ER.utf8, ca_AD.utf8, ca_ES.utf8, ca_FR.utf8, ca_IT.utf8, crh_UA.utf8, cs_CZ.utf8, csb_PL.utf8, cv_RU.utf8, cy_GB.utf8, da_DK.utf8, de_AT.utf8, de_BE.utf8, de_CH.utf8, de_DE.utf8, de_LU.utf8, doi_IN.utf8, dv_MV.utf8, dz_BT.utf8, el_CY.utf8, el_GR.utf8, en_AG.utf8, en_AU.utf8, en_BW.utf8, en_CA.utf8, en_DK.utf8, en_GB.utf8, en_HK.utf8, en_IE.utf8, en_IN.utf8, en_NG.utf8, en_NZ.utf8, en_PH.utf8, en_SG.utf8, en_US.utf8, en_ZA.utf8, en_ZM.utf8, en_ZW.utf8, es_AR.utf8, es_BO.utf8, es_CL.utf8, es_CO.utf8, es_CR.utf8, es_CU.utf8, es_DO.utf8, es_EC.utf8, es_ES.utf8, es_GT.utf8, es_HN.utf8, es_MX.utf8, es_NI.utf8, es_PA.utf8, es_PE.utf8, es_PR.utf8, es_PY.utf8, es_SV.utf8, es_US.utf8, es_UY.utf8, es_VE.utf8, et_EE.utf8, eu_ES.utf8, fa_IR.utf8, ff_SN.utf8, fi_FI.utf8, fil_PH.utf8, fo_FO.utf8, fr_BE.utf8, fr_CA.utf8, fr_CH.utf8, fr_FR.utf8, fr_LU.utf8, fur_IT.utf8, fy_DE.utf8, fy_NL.utf8, ga_IE.utf8, gd_GB.utf8, gez_ER.utf8, gez_ET.utf8, gl_ES.utf8, gu_IN.utf8, gv_GB.utf8, ha_NG.utf8, he_IL.utf8, hi_IN.utf8, hne_IN.utf8, hr_HR.utf8, hsb_DE.utf8, ht_HT.utf8, hu_HU.utf8, hy_AM.utf8, ia_FR.utf8, id_ID.utf8, ig_NG.utf8, ik_CA.utf8, is_IS.utf8, it_CH.utf8, it_IT.utf8, iu_CA.utf8, iw_IL.utf8, ja_JP.utf8, ka_GE.utf8, kk_KZ.utf8, kl_GL.utf8, km_KH.utf8, kn_IN.utf8, ko_KR.utf8, kok_IN.utf8, ks_IN.utf8, ku_TR.utf8, kw_GB.utf8, ky_KG.utf8, lb_LU.utf8, lg_UG.utf8, li_BE.utf8, li_NL.utf8, lij_IT.utf8, lo_LA.utf8, lt_LT.utf8, lv_LV.utf8, mag_IN.utf8, mai_IN.utf8, mg_MG.utf8, mhr_RU.utf8, mi_NZ.utf8, mk_MK.utf8, ml_IN.utf8, mn_MN.utf8, mni_IN.utf8, mr_IN.utf8, ms_MY.utf8, mt_MT.utf8, my_MM.utf8, nb_NO.utf8, nds_DE.utf8, nds_NL.utf8, ne_NP.utf8, nhn_MX.utf8, niu_NU.utf8, niu_NZ.utf8, nl_AW.utf8, nl_BE.utf8, nl_NL.utf8, nn_NO.utf8, nr_ZA.utf8, nso_ZA.utf8, oc_FR.utf8, om_ET.utf8, om_KE.utf8, or_IN.utf8, os_RU.utf8, pa_IN.utf8, pa_PK.utf8, pap_AN.utf8, pl_PL.utf8, ps_AF.utf8, pt_BR.utf8, pt_PT.utf8, ro_RO.utf8, ru_RU.utf8, ru_UA.utf8, rw_RW.utf8, sa_IN.utf8, sat_IN.utf8, sc_IT.utf8, sd_IN.utf8, se_NO.utf8, shs_CA.utf8, si_LK.utf8, sid_ET.utf8, sk_SK.utf8, sl_SI.utf8, so_DJ.utf8, so_ET.utf8, so_KE.utf8, so_SO.utf8, sq_AL.utf8, sq_MK.utf8, sr_ME.utf8, sr_RS.utf8, ss_ZA.utf8, st_ZA.utf8, sv_FI.utf8, sv_SE.utf8, sw_KE.utf8, sw_TZ.utf8, szl_PL.utf8, ta_IN.utf8, ta_LK.utf8, te_IN.utf8, tg_TJ.utf8, th_TH.utf8, ti_ER.utf8, ti_ET.utf8, tig_ER.utf8, tk_TM.utf8, tl_PH.utf8, tn_ZA.utf8, tr_CY.utf8, tr_TR.utf8, ts_ZA.utf8, tt_RU.utf8, ug_CN.utf8, uk_UA.utf8, unm_US.utf8, ur_IN.utf8, ur_PK.utf8, ve_ZA.utf8, vi_VN.utf8, wa_BE.utf8, wae_CH.utf8, wal_ET.utf8, wo_SN.utf8, xh_ZA.utf8, yi_US.utf8, yo_NG.utf8, yue_HK.utf8, zh_CN.utf8, zh_HK.utf8, zh_SG.utf8, zh_TW.utf8, zu_ZA.utf8

Click discovered that you exported a UTF-8 locale but the locale system could not pick up from it because it does not exist. The exported locale is 'C.UTF-8' but it is not supported

```

 - solve with below commands

```bash
export LC_ALL=en_US.utf-8

export LANG=en_US.utf-8
```


### Unknown terminal type
```bash
test@localhost:~# top
'xterm-256color': unknown terminal type.
test@localhost
```
## Solution
```bash
export TERM=xterm

```

### Improve nanorc
```bash
nano ~/.nanorc
```
```
set nowrap
set softwrap
set const
## Nanorc files
include "/usr/share/nano/nanorc.nanorc"

## C/C++
include "/usr/share/nano/c.nanorc"

## HTML
include "/usr/share/nano/html.nanorc"

## TeX
include "/usr/share/nano/tex.nanorc"

## Quoted emails (under e.g. mutt)
include "/usr/share/nano/mutt.nanorc"

## Patch files
include "/usr/share/nano/patch.nanorc"

## Manpages
include "/usr/share/nano/man.nanorc"

## Groff
include "/usr/share/nano/groff.nanorc"

## Perl
include "/usr/share/nano/perl.nanorc"

## Python
include "/usr/share/nano/python.nanorc"

## Ruby
include "/usr/share/nano/ruby.nanorc"

## Java
include "/usr/share/nano/java.nanorc"

## Assembler
include "/usr/share/nano/asm.nanorc"

## Bourne shell scripts
include "/usr/share/nano/sh.nanorc"

## POV-Ray
include "/usr/share/nano/pov.nanorc"
```
