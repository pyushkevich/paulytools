#include <stdio.h>
#include <string.h>

int usage()
{
  printf("tempfile - create temporary file\n");
  printf("usage: \n");
  printf("  tempfile pattern\n");
  printf("notes: \n");
  printf("  The pattern must end with 6 symbols XXXXXX, e.g., /tmp/sock_XXXXXX\n");
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 2)
    return usage();

  char tempy[L_tmpnam];
  strcpy(tempy, argv[1]);
  int handle = mkstemp(tempy);
  if(handle >= 0)
    printf("%s\n",tempy);
}
