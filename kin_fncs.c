/*
Author: Fotios Lygerakis
email: fotios.lygerakis@mavs.uta.edu
*/

#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415927
#endif
const int M = 4;
const int N = 4;
const double l_0 = 0.25;
const double l_1 = 0.2;
const double l_2 = 0.2;
const double l_3 = 0.15;
const double d_1 = -0.04;
const double d_2 = 0.04;
const double d_3 = -0.04;
const double d_4 = -0.04;

void display(double result[][N]);
void display_array(double result[], int n);
void copy_array(double array_1[], double array_2[]);

// function to multiply two matrices
void mat_mul(double first[][4],
                      double second[][4],
                      double result[][4]) {

   // Multiplying first and second matrices and storing it in result
   for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
      	 result[i][j] = 0;
         for (int k = 0; k < N; ++k) {
            result[i][j] += first[i][k] * second[k][j];
         }
      }
   }
}

void array_mul(double array1D[], double array2D[][4], double final_result[]) {
  for (int i = 0; i < 4; i++) { // use nested loops to multiply the two arrays together
      final_result[i] = 0;
      for (int j = 0; j < 4; j++) {
          final_result[i] += array2D[i][j] * array1D[j]; // sum the products of each set of elements
      }
   }
} // end of getOutputArray method

void display_array(double result[], int n) {

   printf("\ array:\n");
   for (int i = 0; i < n; ++i) {
	printf("%lf ", result[i]);
   }
   printf("\n");
}
void display(double result[M][N]) {

   printf("\Matrix:\n");
   for (int i = 0; i < M; ++i) {
      for (int j = 0; j < M; ++j) {
         printf("%lf  ", result[i][j]);
         if (j == 4 - 1)
            printf("\n");
      }
   }
}

void copy_array(double array_1[], double array_2[]){
// copies the elements of array 1 to array 2
	for (int i=0; i<M-1; i++){
		array_2[i] = array_1[i];

	}
}
void copy_matrix(double array_1[][M], double array_2[][M-1]){
// copies the elements of array 1 to array 2
	for (int i=0; i<M-1; i++){
		for (int j=0; j<M-1; j++){
			array_2[i][j] = array_1[i][j];
		}
	}
}

double* fwd_kin(theta, x)
double *theta; //6 joint angles theta (the last angle corresponds to the gripper andis of no interest here)
double x[3]; // 3 dimensional position of the tool frame (x[0] = x, x[1]= y, x[2] = z)

{
double x_augmented[4];
double D_z_l_0[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
double D_y_d_1[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
double D_x_l_1[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
double D_y_d_2[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
double D_x_l_2[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
double D_y_d_3[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
double D_z_d_4[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
double D_x_l_3[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
double R_z_theta_0[][4] = {
			{cos( theta[0]), -sin( theta[0]),0,0},
			{sin( theta[0]),cos( theta[0]),0,0},
			{0,0,1,0},
			{0,0,0,1}};
double R_y_theta_1[][4] = {
			{cos(theta[1]), 0, sin(theta[1]),0},
			{0, 1, 0,0},
			{-1*sin( theta[1]), 0, cos( theta[1]),0},
			{0,0,0,1}};
double R_y_theta_2[][4] = {
			{cos( theta[2]), 0, sin( theta[2] ),0},
			{0, 1, 0,0},
			{-1*sin( theta[2]), 0, cos( theta[2]),0},
			{0,0,0,1}};
double R_y_theta_3[][4] = {
			{cos( theta[3] ), 0, sin( theta[3]),0},
			{0, 1, 0,0},
			{-1*sin( theta[3]), 0, cos( theta[3]),0},
			{0,0,0,1}};
double R_x_theta_4[][4] = {
			{1,0,0,0},
			{0,cos( theta[4]) ,-1*sin( theta[4]) ,0},
			{0,sin( theta[4]),cos( theta[4]) ,0},
			{0,0,0,1}};


// prepare the disposition matrices
int i_x = 0;
int i_y = 1;
int i_z = 2;
D_z_l_0[i_z][3] = l_0; // z
D_y_d_1[i_y][3] = d_1;	// y
D_x_l_1[i_x][3] = l_1;	// x
D_y_d_2[i_y][3] = d_2;	// y
D_x_l_2[i_x][3] = l_2;	// x
D_y_d_3[i_y][3] = d_3;	// y
D_z_d_4[i_z][3] = d_4;	// z
D_x_l_3[i_x][3] = l_3;	// x
// augment x vector from 3 to 4 dimensions
for (int i=0;i<4;i++){
	x_augmented[i] = 0;
	if (i == 3){
	x_augmented[i] = 1;
	}
}
// do the multiplications
double result[4][4];
mat_mul(D_z_l_0, R_z_theta_0, result);
double result_1[4][4];
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);
mat_mul(result_1, R_y_theta_1, result);
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);
mat_mul(result_1, D_y_d_1, result);
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);
mat_mul(result_1, D_x_l_1, result);
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);
mat_mul(result_1, R_y_theta_2, result);
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);
mat_mul(result_1, D_y_d_2, result);
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);
mat_mul(result_1, D_x_l_2, result);
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);
mat_mul(result_1, R_y_theta_3, result);
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);
mat_mul(result_1, D_z_d_4, result);
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);
mat_mul(result_1, D_y_d_3, result);
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);
mat_mul(result_1, D_x_l_3, result);
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);
mat_mul(result_1, R_x_theta_4, result);
memcpy (result_1, result, M*N*sizeof(double));
//display(result_1);

//display_array(result);

double final[4]; 
array_mul(x_augmented, result_1, final);
//display_array(final);
memcpy(x, &final[0], 3*sizeof(*final));
//display_array(x);
return x;
}




double* inv_kin(x, theta)
double *x;
double theta[6];
{
// calcualte helping angles for theta_0 (see notes)
double theta_prime = atan2(x[1], x[0]);
double theta_prime2 =  asin(d_3/ sqrt( pow(x[0], 2) + pow(x[1], 2)));
theta[0] = theta_prime - theta_prime2;

double X_arc[3] = {0 , 0, -1};
double x_prime[3];
// move the ee point towards l3 on the z axis
for (int i=0; i<3; i++){
	x_prime[i] = x[i]- l_3*X_arc[i];
} 
// translate point on the frame of the rotated robot around theta_0
double x_prime2[3];
x_prime2[0] = sqrt( pow(x_prime[0],2) + pow(x_prime[1],2) - pow(d_4,2)) - d_4; 
x_prime2[1] = 0;
x_prime2[2] = x_prime[2]-l_0;

// calculate alpha and gamma angles (see on the notes)
double alpha = acos((pow(x_prime2[0],2) + pow(x_prime2[2],2)) 
		/ (2 * l_1 * sqrt(pow(x_prime2[0],2) + pow(x_prime2[2],2))));
double gamma = atan2(x_prime2[2], x_prime2[0]);
// choose positive alpha
double theta_small = gamma + alpha;
theta[1] = -theta_small;
theta[2] = acos((pow(x_prime2[0],2) + pow(x_prime2[2],2) - pow(l_1,2) - pow(l_2,2))/(2*l_1*l_2));
theta[3] = M_PI/2 - theta[1] - theta[2];
theta[4] = 0;
theta[5] = 0;
return theta;
}








