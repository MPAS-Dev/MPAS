void gen_namelists(struct namelist *);
void gen_history_attributes(char * modelname, char * corename, char * version);
void gen_field_defs(struct group_list * groups, struct variable *, struct dimension *);
void gen_reads(struct group_list * groups, struct variable *, struct dimension *);
void gen_writes(struct group_list * groups, struct variable *, struct dimension *, struct namelist *);
