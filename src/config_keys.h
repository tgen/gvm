#ifndef GVM_CONFIG_KEYS_H
#define GVM_CONFIG_KEYS_H

#define GVM_CONFIG_bam_list         "bamList"
#define GVM_CONFIG_ref_file         "refGenome"
#define GVM_CONFIG_bed_file         "regionsFile"
#define GVM_CONFIG_snp_vcf_path     "snpVCFpath"
#define GVM_CONFIG_snp_vcf_name     "snpVCFname"
#define GVM_CONFIG_cosm_vcf_path    "cosmicVCF"
#define GVM_CONFIG_normal_base_path "NormalBase"
#define GVM_CONFIG_normal_base_out  "outfile"
#define GVM_CONFIG_out_name         "outName"

#define GVM_CONFIG_CHECK_STROPT(N) \
    if (strcmp(option, GVM_CONFIG_##N) == 0) { \
        strncpy(s->N, value, sizeof(s->N)); \
        return 1; \
    }

#endif
