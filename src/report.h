#ifndef GVM_REPORT_H
#define GVM_REPORT_H


static const char *alignment_type_strings[] = {
	"SINGLE", "INSERTION", "DELETION", "SOFTCLIP"
};


enum alignment_type { at_single = 0, at_ins, at_del, at_sclip };
static inline const char *at_tostr(enum alignment_type x)
	{ return alignment_type_strings[x]; }

struct alignment_report {
	enum alignment_type align_type;

	/* position in the entire genome */
	uint32_t pos;
	/* position in this single read */
	uint32_t spos;

	/* data associated with the variant
	   (can be a base or many bases) */
	uint32_t data;

	/* only used for indels */
	uint32_t size;
};

struct context;

typedef void (*reporter_func)(struct context *, struct alignment_report, void *);

#endif
