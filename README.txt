Various notes/comments taken from the files. Need to be cleaned up.

/* factor to make the units work...
    energy (joule) = 1/2 mass(kg) velocity(m/s)^2
   we have:
    energy in eV = 1.60217653E-19 Joules
    mass in au = 1.66E-27 Kg
    distance in mm = 1E-3 m
    time in ns = 1E-9 sec
   so this works out to:
    energy(eV) = 5180.45 mass(au) velocity(mm/ns)^2
*/
#define	HALF	    5180.45

/* what is a mass */
#define	MASS	    double

/* particle masses, in atomic mass units, one au = 1.66E-27 KG */
MASS particle_mass[NPARTICLE] = { 1.007276466812, 5.4857990946E-4, 1.00866491600, 4.001506179125 };

/* particle masses, in atomic mass units, one au = 1.66E-27 KG */
extern REAL particle_mass[NPARTICLE];
extern char *particle_name[NPARTICLE];
