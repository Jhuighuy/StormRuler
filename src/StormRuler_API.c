#include <stdio.h>

void SR_PrintPointer(const void* p) {
  fprintf(stdout, "PRINT_PTR=%p\n", p);
  fflush(stdout);
} // SR_PrintPointer
