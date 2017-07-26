#ifndef GVM_UTILS_H
#define GVM_UTILS_H

#define err_printf_helper(fmtstr, ...) fprintf(stderr, "%s:%d: " fmtstr "%c", __FILE__, __LINE__, __VA_ARGS__)
#define err_printf(...) err_printf_helper(__VA_ARGS__, 0)

#ifdef VERBOSE
#define verbose_fprintf fprintf
#else
#define verbose_fprintf(...) do{} while(0);
#endif

#define MAX(X,Y) ( (X) > (Y) ? (X) : (Y) )

#endif
