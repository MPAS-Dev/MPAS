
#define TABLESIZE 271

struct dnode {
   char key[1024];
   struct dnode * next;
};

struct dtable {
   int size;
   struct dnode * table[TABLESIZE];
};

void dict_alloc(struct dtable **);
void dict_insert(struct dtable *, char *);
void dict_remove(struct dtable *, char *);
int dict_search(struct dtable *, char *);
int dict_size(struct dtable *);
void dict_free(struct dtable **);
