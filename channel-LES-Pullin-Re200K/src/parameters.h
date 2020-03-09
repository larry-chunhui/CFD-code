extern void init_par(
                     parameters *par);

extern void init_var(
                     int flag,
                     variables *var,
                     parameters *par);

extern void init_vel_bc(
                        variables *var,
                        parameters *par);

extern void finalize_par(
                         parameters *par);

extern void finalize_var(
                         int flag,
                         variables *var,
                         parameters *par);

extern void finalize_vel_bc(
                            variables *var,
                            parameters *par);
