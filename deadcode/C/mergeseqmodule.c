#include <Python.h>
#include <string.h>

#define min(a,b) ((a) < (b) ? (a) : (b))

// Declarations.
char *merge(const char *, const char *, const char *,
      const char *, double);
int find_frame(const char *, const char*, const double);

// Python interface.
static PyObject *
mergeseq_merge(PyObject *self, PyObject *args)
{

    double errors;
    const char *seq1;
    const char *seq2;
    const char *qual1;
    const char *qual2;

    if (!PyArg_ParseTuple(args, "ssssd",
             &seq1, &seq2, &qual1, &qual2, &errors))
        return NULL;

    char *consensus = merge(seq1, seq2, qual1, qual2, errors);

    PyObject * value = Py_BuildValue("s", consensus);
    free(consensus);

    return value;

}


static PyMethodDef MergeSeqMethods[] = {
    {"merge",  mergeseq_merge, METH_VARARGS, "Merge sequences."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
initmergeseq(void)
{
    (void) Py_InitModule("mergeseq", MergeSeqMethods);
}

// C functions.
char *
merge(
   const char * seq1,
   const char * seq2,
   const char * qual1,
   const char * qual2,
         double errors
)
{

   const int ov = find_frame(seq1, seq2, errors);
   if (!ov) return NULL;

   size_t sz1 = strlen(seq1);
   size_t sz2 = strlen(seq2);

   // Initialize consensus.
   char *cons = malloc(sz1 + sz2 - ov + 1);
   strncpy(cons, seq1, sz1-ov);
   strncpy(cons + (sz1-ov), seq2, sz2 + 1);

   int i;

   for (i = 0 ; i < ov ; i++) {
      // In case of conflict use read quality.
      // Remember that the nucleotide by default if seq2.
      if (seq1[sz1-ov+i] != seq2[i] && qual1[sz1-ov+i] > qual2[i]) {
         // DEBUG
         cons[sz1-ov+i] = seq1[sz1-ov+i];
      }
   }

   return cons;

}


int
find_frame
(
   const char * seq1,
   const char * seq2,
   const double errors
)
{

   size_t sz1 = strlen(seq1);
   size_t sz2 = strlen(seq2);

   const int ovmin = 20;
   const int ovmax = min(sz1, sz2);

   int i;
   int ov;

   for (ov = ovmin ; ov < ovmax ; ov++) {
      int nerrors = 0;
      for (i = 0 ; i < ov ; i++) {
         nerrors += (seq1[sz1-ov+i] != seq2[i] || seq2[i] == 'N');
      }
      // If the number of errors is less than tolerance
      // we found the frame: return the overlap length.
      if (nerrors < (ov * errors)) return ov;
   }

   // Did not find the frame
   return 0;

}


