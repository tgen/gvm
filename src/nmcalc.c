#include <assert.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include "nmcalc.h"

void get_ab(struct variant_table *v, struct variant_counts **a, struct variant_counts **b);
int check_ab(struct variant_counts *a, struct variant_counts *b);

// Convenience function to mimic MATLAB function
static inline double binopdf(unsigned int x, unsigned int n, double p) {
	return gsl_ran_binomial_pdf(x, p, n);
}

#define NONZERO_OR_NAN(X, NAN_CODE) ((X) == (X) ? (X) : nan(#NAN_CODE))

void nmcalc(struct context *context,
            struct variant_table *v,
            double prior_map_error,
            int ploidy,
            struct nm_entry *out)
{
	(void) context;

	double a_read_depth, b_read_depth;
	double a_mean_mq, b_mean_mq;

	double prior_hom, prior_het;
	double p_data_map_error, p_data_hom, p_data_het, p_data;
	double p_map_error;

	double per_read_pass, ab_frac;

	struct variant_counts *a, *b;

	get_ab(v, &a, &b);
	assert(check_ab(a, b));

	// Initialization of useful variables

	a_read_depth = (a->count_f + a->count_r) / 2.0;
	b_read_depth = (b->count_f + b->count_r) / 2.0;

	a_mean_mq = a->total_mq / NONZERO_OR_NAN(a_read_depth, 200);
	b_mean_mq = b->total_mq / NONZERO_OR_NAN(b_read_depth, 200);

	// Normal metrics calculation starts here

	// priors
	if (ploidy == 2) {
		prior_het = 2 * a->pop_af * b->pop_af;
	} else {
		prior_het = 0;
	}

	prior_hom = a->pop_af * a->pop_af;

	// likelihoods
	p_data_map_error = pow(10, -fmin(a_mean_mq, b_mean_mq)/10);
	p_data_hom = binopdf(b_read_depth, v->read_count_pass, pow(10, (-b_mean_mq / 10)));
	p_data_het = binopdf(b_read_depth, v->read_count_pass,  0.5);

	p_data = prior_map_error * p_data_map_error + prior_het * p_data_het + prior_hom * p_data_hom;

	// posteriors
	if (v->read_count_pass > 0) {
		p_map_error = prior_map_error * p_data_map_error / p_data;
	} else {
		p_map_error = nan("300");
	}

	// other metrics
	per_read_pass = v->read_count_pass / NONZERO_OR_NAN(v->read_count, 200);

	ab_frac = (a_read_depth + b_read_depth) / NONZERO_OR_NAN(v->read_count_pass, 200);

	out->pos = v->offset;

	if (ploidy == 1 || ploidy == 2) {
		out->norm_read_depth = (3-ploidy) * v->read_count_pass; // sneak
		out->prob_map_err = p_map_error;
		out->read_pass = per_read_pass;
		out->a_or_b = ab_frac;
	} else  {
		out->norm_read_depth = nan("101");
		out->prob_map_err = nan("102");
		out->read_pass = nan("103");
		out->a_or_b = nan("104");
	}
}
