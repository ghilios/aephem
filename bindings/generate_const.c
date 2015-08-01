#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <aephem.h>

const char *blacklist[] = {"AEPHEM_H", 
                           "AE_FROM_J2000", "AE_TO_J2000",
                           "AE_ADD_ABERRATION", "AE_REMOVE_ABERRATION",
                           "AE_PHYSICAL_USE_D", "AE_PHYSICAL_USE_T",
                           NULL};

void do_python(char **const_name, char **enum_name) {
  int i;

  printf("// This code automatically generated. Changes made here will be "
         "lost.\n");
  printf("#include \"aephem_py.h\"\n");
  printf("\n");
  printf("const struct aepy_const_t aepy_const_list[] = {\n");
  for (i = 0; const_name[i] != NULL; i++)
    printf("  {\"%s\", %s},\n", const_name[i] + 3, const_name[i]);
  printf("{NULL, 0}\n};\n\n");

  printf("const struct aepy_const_t aepy_enum_list[] = {\n");
  for (i = 0; enum_name[i] != NULL; i++)
    printf("  {\"%s\", %s},\n", enum_name[i] + 3, enum_name[i]);
  printf("{NULL, 0}\n};\n");
}

int in_blacklist(const char *name) {
  int i;
  
  for (i = 0; blacklist[i] != NULL; i++) {
    if (!strcmp(name, blacklist[i]))
      return 1;
  }

  if (strncmp(name, "AE_", 3))
    return 1;

  return 0;
}

void skip_whitespace(char **ptr) {
  while (**ptr == ' ' || **ptr == '\t')
    (*ptr)++;

  return;
}

int main(int argc, char *argv[]) {
  FILE *fp;
  int i, n_const, n_enum, do_enum;
  char line[128], *line_ptr, tmp_name[128], **const_name, **enum_name;
  
  if (argc != 3) {
    fprintf(stderr, "Usage: generate_const <relative path to aephem.h> p\n");
    exit(1);
  }

  if ((fp = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "Could not open %s for reading.\n", argv[1]);
    exit(1);
  }

  n_const = 0;
  const_name = NULL;
  n_enum = 0;
  enum_name = NULL;
  do_enum = 0;
  while (fgets(line, 127, fp) != NULL) {
    line_ptr = line;
    skip_whitespace(&line_ptr);

    if (!strncmp(line_ptr, "enum", 3))
      do_enum = 1;
    else if (do_enum) {
      if (!strncmp(line_ptr, "};", 2))
        do_enum = 0;
      else if (strncmp(line_ptr, "//", 2)) {
        sscanf(line_ptr, "%s\n", tmp_name);
        if (tmp_name[strlen(tmp_name) - 1] == ',')
          tmp_name[strlen(tmp_name) - 1] = '\0';
        if (!in_blacklist(tmp_name)) {
          enum_name = (char **)realloc(enum_name,
                                       (n_enum + 1) * sizeof(char *));
          enum_name[n_enum++] = strdup(tmp_name);
        }
      }
    }
    else if (!strncmp(line_ptr, "#define", 7)) {
      line_ptr += 7;
      skip_whitespace(&line_ptr);
      sscanf(line_ptr, "%s", tmp_name);
      for (i = 0; i < strlen(tmp_name); i++) {
        if (islower(tmp_name[i]))
          break;
      }
      if (i == strlen(tmp_name) && !in_blacklist(tmp_name)) {
        const_name = (char **)realloc(const_name, 
                                      (n_const + 1) * sizeof(char *));
        const_name[n_const++] = strdup(tmp_name);
      }
    }
  }

  const_name = (char **)realloc(const_name, (n_const + 1) * sizeof(char *));
  const_name[n_const] = NULL;
  enum_name = (char **)realloc(enum_name, (n_enum + 1) * sizeof(char *));
  enum_name[n_enum] = NULL;

  if (!strcmp(argv[2], "p"))
    do_python(const_name, enum_name);

  return 0;
}
