#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* line;
int max_line_len;

char* readline(FILE* input) {
  int len;
  
  if (fgets(line, max_line_len, input) == NULL) {
    return NULL;
  }

  while(strrchr(line, '\n') == NULL) {
    max_line_len *= 2;
    line = (char*) realloc(line, max_line_len);
    len = (int) strlen(line);
    if(fgets(line + len, max_line_len - len, input) == NULL)
      break;
  }
  return line;
}
