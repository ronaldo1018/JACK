#include    <unistd.h>
#include    "jack.h"

#if PAR_DEBUG
static const char *ParamNames[] = {
    "P_INTEGER", "P_BOOLEAN", "P_REAL", "P_STRING", "P_LIST", "P_FLAG", "P_PATH", "P_READ", "P_WRIT"
};
#endif

static void ZeroValue(enum ParamType type, void *dst) {
    assert(dst);
    switch (type) {
	case P_FLAG:
	case P_BOOLEAN:
	case P_INTEGER: {
	    int	*tmp = (int *) dst;
	    *tmp = 0;
	    break;
	}
	case P_REAL: {
	    double *tmp = (double *) dst;
	    *tmp = 0.0;
	    break;
	}
	case P_WRIT:
	case P_READ:
	case P_PATH:
	case P_LIST:
	case P_STRING: {
	    char **tmp = (char **) dst;
	    *tmp = NULL;
	    break;
	}
    }
}

static int SetValue(enum ParamType type, char *src, void *dst, int quiet) {
    char *ptr;
    if (type != P_FLAG) assert(src);
    assert(dst);
    switch (type) {
	case P_INTEGER: {
	    int	*tmp, junk;
	    tmp = (int *) dst;
	    junk = strtol(src, &ptr, 0);
	    if (ptr == src) {
		if (! quiet) fprintf(stderr, "could not interpret \"%s\" as an integer\n", src);
		return(-1);
	    }
	    *tmp = junk;
	    break;
	}
	case P_FLAG: {
	    int	*tmp;
	    tmp = (int *) dst;
	    *tmp = 1;
	    break;
	}
	case P_BOOLEAN: {
	    int	*tmp, junk;
	    tmp = (int *) dst;
	    junk = strtol(src, &ptr, 0);
	    if (ptr == src) {
		if (! quiet) fprintf(stderr, "could not interpret \"%s\" as an integer\n", src);
		return(-1);
	    }
	    *tmp = junk;
	    break;
	}
	case P_REAL: {
	    double *tmp, junk;
	    tmp = (double *) dst;
	    junk = strtod(src, &ptr);
	    if (ptr == src) {
		if (! quiet) fprintf(stderr, "could not interpret \"%s\" as a real number\n", src);
		return(-1);
	    }
	    *tmp = junk;
	    break;
	}
	case P_PATH: {
	    if (access(src, R_OK)) {
		if (! quiet) fprintf(stderr, "could not open \"%s\" as a path\n", src);
		return(-1);
	    }
	    char **tmp;
	    tmp = (char **) dst;
	    *tmp = strdup(src);
	    break;
	}
	case P_READ: {
	    char **tmp;
	    tmp = (char **) dst;
	    if (src) {
		if (access(src, R_OK)) {
		    if (! quiet) fprintf(stderr, "could not open \"%s\" as a file to read\n", src);
		    return(-1);
		}
		*tmp = strdup(src);
	    } else {
		*tmp = NULL;
	    }
	    break;
	}
	case P_WRIT: {
	    char **tmp;
	    tmp = (char **) dst;
	    if (src) {
		if (! access(src, F_OK)) {
		    if (access(src, W_OK)) {
			if (! quiet) fprintf(stderr, "could not open \"%s\" as a file to write\n", src);
			return(-1);
		    }
		}
		*tmp = strdup(src);
	    } else {
		*tmp = NULL;
	    }
	    break;
	}
	case P_LIST:
	case P_STRING: {
	    char **tmp;
	    tmp = (char **) dst;
	    *tmp = strdup(src);
	    break;
	}
    }
    return(0);
}

int ParseParams(int argc, char **argv, ParamDescriptor *table, void **destinations) {
    static int set[128];
    int i, j, ntab;
    char *name, *value;
/* first set defaults if they were specified and mark all fields as not being set */
    for (ntab=0; table[ntab].key; ntab++) {
	if (table[ntab].deflt && table[ntab].type != P_FLAG) {
#if PAR_DEBUG
	    fprintf(stderr, "SetValue %d: type=%s src=%s\n", ntab, ParamNames[table[ntab].type], table[ntab].deflt);
#endif
	    if (SetValue(table[ntab].type, (char*) table[ntab].deflt, destinations[ntab], 1) != 0) {
/* default appears to be incorrect, so we mark this variable as one that we need a setting for */
		table[ntab].required = 1;
        printf("%d!!!!!!!!!!!!!!!!!!!!!!!\n",table[ntab].required);
	    }
	} else {
#if PAR_DEBUG
	    fprintf(stderr, "ZeroValue %d: type=%s src=%s\n", ntab, ParamNames[table[ntab].type], table[ntab].deflt);
#endif
	    ZeroValue(table[ntab].type, destinations[ntab]);
	}
	set[ntab] = 0;
    }
    assert(ntab<128);
/* check to see if the arguments are positional or via -<flags> */
    if (argc && argv[0][0] == '-') {
/* go through the arguments */
	for (i=0; i<argc; i++) {
	    name = argv[i];
	    if (name[0] != '-') {
		fprintf(stderr, "flag %s missing initial \"-\"\n", name);
		return(-1);
	    }
	    ++name;
	    for (j=0; j<ntab; j++) {
		if (! strcmp(name, table[j].key)) break;
	    }
	    if (! table[j].key) {
		fprintf(stderr, "flag %s is not valid\n", name);
		fprintf(stderr, "valid flags are:\n");
		for (j=0; table[j].key; j++) {
		    fprintf(stderr, "\t%s\n", table[j].key);
		}
		return(-1);
	    }
/* flags do not get values */
	    if (table[j].type == P_FLAG) {
		value = NULL;
	    } else {
		if (++i < argc) {
		    value = argv[i];
		} else {
		    fprintf(stderr, "flag %s requires a value\n", name);
		    return(-1);
		}
	    }
#if PAR_DEBUG
	    fprintf(stderr, "SetValue %d: type=%s src=%s new value=%s\n", j, ParamNames[table[j].type], table[j].deflt, value);
#endif
	    if (SetValue(table[j].type, value, destinations[j], 0) != 0) {
		return(-1);
	    }
	    set[j] = 1;
	}
    }
    j = 0;
/* make sure those parameters without defaults were specified */
    for (i=0; i<ntab; i++) {
	if (table[i].required && set[i] == 0) {
	    fprintf(stderr, "parameter %s must be specified\n", table[i].key);
	    j = -1;
	}
    }
    return(j);
}
