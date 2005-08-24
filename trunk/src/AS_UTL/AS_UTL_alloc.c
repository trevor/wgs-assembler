#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "AS_UTL_alloc.h"


// Allocate and return a pointer to an array of  num  elements of
// len  bytes each.  All are set to 0.  Exit if fai.
//
void *
safe_calloc(size_t num, size_t len) {
  void  *p;

   p = calloc(num, len);
   if  (p == NULL) {
     fprintf(stderr, "Could not calloc memory (%d number * %d bytes) \n",
             (int)num,
             (int)len);
     assert(p != NULL);
   }
   return(p);
}



// Allocate and return a pointer to len bytes of memory.
// Len  bytes each.  Exit if fail.
//
void *
safe_malloc(size_t len) {
  void  *p;

  p = malloc(len);
  if (p == NULL) {
    fprintf(stderr, "Could not malloc memory (%d bytes) \n",(int) len);
    assert(p != NULL);
  }
  
  return(p);
}




// Reallocate memory for q to len  bytes and return a pointer
// to the new memory.  Exit if fail.
//
void *
safe_realloc(void *q, size_t len) {
  void  *p;

  p = realloc(q, len);
  if (p == NULL) {
    fprintf(stderr, "Could not malloc memory (%d bytes) \n",(int) len);
    assert(p != NULL);
  }
  
  return(p);
}
