#ifndef GVM_CONFIG_KEYS_H
#define GVM_CONFIG_KEYS_H

#include <stdint.h>

#define GVM_CONFIG_bam_list         "bamList"
#define GVM_CONFIG_ref_file         "refGenome"
#define GVM_CONFIG_bed_file         "regionsFile"
#define GVM_CONFIG_snp_vcf_path     "snpVCFpath"
#define GVM_CONFIG_snp_vcf_name     "snpVCFname"
#define GVM_CONFIG_cosm_vcf_path    "cosmicVCF"
#define GVM_CONFIG_normal_base_path "NormalBase"
#define GVM_CONFIG_normal_base_out  "outfile"
#define GVM_CONFIG_out_name         "outName"
#define GVM_CONFIG_min_bq           "minBQ"
#define GVM_CONFIG_min_mq           "minMQ"
#define GVM_CONFIG_default_bq       "defaultBQ"
#define GVM_CONFIG_default_mq       "defaultMQ"
#define GVM_CONFIG_min_b_count      "minBCount"
#define GVM_CONFIG_pv_freq          "pvFreq"
#define GVM_CONFIG_prior_map_error  "priorMapError"

#define GVM_CONFIG_CHECK_STROPT(N) \
    if (strcmp(option, GVM_CONFIG_##N) == 0) { \
        strncpy(s->N, value, sizeof(s->N)); \
        return 1; \
    }

/* I had to go to a confessional after writing this */
#define GVM_CONFIG_CHECK_NUMOPT(N) \
	if (strcmp(option, GVM_CONFIG_##N) == 0) { \
		N = _Generic((N), /* This is how to DRY, right? */ \
                        unsigned int: strtol(value, &tmp, 10), \
                        double: strtod(value, &tmp)); \
		if (*tmp != '\0') { \
			err_printf("Invalid value for " GVM_CONFIG_##N ": %s\n", value); \
			return 0; \
		} \
		s->N = (N); \
		return 1; \
	}

#endif
