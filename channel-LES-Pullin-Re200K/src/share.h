void get_periodic_data(variables *var, neighbors *shared, parameters *par);

void get_neignboring_data(variables *var, neighbors *shared, parameters *par);

void finalize_neignboring(neighbors *shared, parameters *par);

void share_u(int my_rank, variables *var, double *u_in_r, double *u_in_l, parameters *par, int tag);

void share_v(int my_rank, variables *var, double *v_in_r, double *v_in_l, parameters *par, int tag);

void share_w(int my_rank, variables *var, double *w_in_r, double *w_in_l, parameters *par, int tag);

void share_cn(int my_rank, double *Txx, double *in_r, double *in_l, parameters *par, int tag);
void share_cn_periodic(int my_rank, double *Txx, double *in_r, double *in_l, parameters *par, int tag);
void share_ct_periodic(int my_rank, double *Tij_ct, double *in_r, double *in_l, parameters *par, int tag);
void share_les_periodic__2(int my_rank, double *les, double *in_l, parameters *par, int tag);
