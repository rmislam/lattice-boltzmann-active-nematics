struct point {
    int i;
    int j;    
};

//Line number in the x direction. Range is from 0 to I-1.
int i_vr(int l) {
    return l % I; 
}

//Line number in the y direction
int j_vr(int l) {
    return (l % (I * J)) / I; 
}

//Characteristic velocity vectors of the LB model.
void initialiseE()
{
    E[0][0] = 0;		E[0][1] = 0;
    E[1][0] = 1;		E[1][1] = 0;
    E[2][0] = 0;		E[2][1] = 1;
    E[3][0] = -1;		E[3][1] = 0;
    E[4][0] = 0; 		E[4][1] = -1;
    E[5][0] = 1;		E[5][1] = 1;
    E[6][0] = -1;		E[6][1] = 1;
    E[7][0] = -1;		E[7][1] = -1;
    E[8][0] = 1;		E[8][1] = -1;
   
    //Constants for the model.
    int i;
    WKONST[0] = 4.f / 9.f;
    for (i = 1; i <= 4; i++) WKONST[i] = 1. / 9.;
    for (i = 5; i <= 8; i++) WKONST[i] = 1. / 36.;
}

//Write the density and velocity field to a file
void write_velocity(int t) {
    char buf[120];
    sprintf(buf, "./%s_velocity_%d.dat", FILE_NAME, t);
    FILE *file;
    file = fopen(buf, "w");
    
    for (int l = 0; l < NMAX; l++) {
        fprintf(file, "%d %d %lf %lf %lf\n", i_vr(l), j_vr(l), U[0][l], U[1][l], U[2][l]);
    }
    
    fclose(file);
}

//Write the degree of order and director angle field to a file
void write_orientation(int t) {
    char buf[120];
    sprintf(buf, "./%s_orientation_%d.dat", FILE_NAME, t);
    FILE *file;
    file = fopen(buf, "w");
    
    for (int l = 0; l < NMAX; l++) {
        double degree_of_order = sqrt(4 * (Q[0][l] * Q[0][l] + Q[1][l] * Q[1][l]));
        double angle = 0.5 * atan2(Q[1][l], Q[0][l]);
        fprintf(file, "%d %d %lf %lf\n", i_vr(l), j_vr(l), degree_of_order, angle);
    }
    
    fclose(file);
}

struct point rotate_point(int l, int i_center, int j_center, double angle) {
    struct point point_rot;
    int i = i_vr(l);
    int j = j_vr(l);
    int rel_i = i - i_center;
    //int rel_j = j - j_center;
    int rel_j = j_center - j;
    double cosTheta = cos(angle);
    double sinTheta = sin(angle);
    point_rot.i = round(i_center + cosTheta * rel_i - sinTheta * rel_j);
    point_rot.j = round(j_center - (sinTheta * rel_i + cosTheta * rel_j));
    return point_rot;
}