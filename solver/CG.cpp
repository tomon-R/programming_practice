//反復解法　GS法
// static const int MAX_ITERATIONS = 10000;
// static const double EPSILON = 0.0000001;
// void solve_mat_GS(double* sol, double** mat, double* rhs, int length) {
//     for (int l = 0; l < MAX_ITERATIONS; l++) {
//         for (int i = 0; i < length; i++) {
//             sol[i] = rhs[i];
//             for (int j = 0; j < length; j++) {
//                 if (i != j) {
//                     sol[i] -= mat[i][j] * sol[j];
//                 }
//             }
//             sol[i] /= mat[i][i];
//         }

//         double residual = 0.0;
//         for (int i = 0; i < length; i++) {
//             double r_i = -rhs[i];
//             for (int j = 0; j < length; j++) {
//                 r_i += mat[i][j] * sol[j];
//             }
//             residual += r_i * r_i;
//         }
//         residual = sqrt(residual);
//         printf("GS_loop %d: %e\n", l, residual);
//         if (residual < EPSILON) return;
//     }
// }
