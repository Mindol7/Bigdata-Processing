#include <mpi.h>
#include <stdio.h>

// Function to integrate
double f(double x) {
    return x * x; // f(x) = x^2
}

double f2(double x){
    return f(x) * f(x) - 3*f(x) + x + 4; // f(x) = x^4-3x^2+x+4
}

// Compute local trapezoidal approximation
// Trap function about f(x) = x^2
double Trap(double left_endpt, double right_endpt, int trap_count, double base_len) {
    double estimate, x;
    int i;

    estimate = (f(left_endpt) + f(right_endpt)) / 2.0;
    for (i = 1; i <= trap_count - 1; i++) {
        x = left_endpt + i * base_len;
        estimate += f(x);
    }
    estimate = estimate * base_len;
    
    return estimate;
}

// Trap function about f(x) = x^4-3x^2+x+4
double Trap_2(double left_endpt, double right_endpt, int trap_count, double base_len) {
    double estimate_2, x;
    int i;

    estimate_2 = (f2(left_endpt) + f2(right_endpt)) / 2.0;
    for (i = 1; i <= trap_count - 1; i++) {
        x = left_endpt + i * base_len;
        estimate_2 += f2(x);
    }
    estimate_2 = estimate_2 * base_len;
    
    return estimate_2;
}

int main(int argc, char** argv) {
    int my_rank, comm_sz;
    int n[4] = {256, 1024, 4096, 16348}; // Array of number of Trapezoid
    double a = 0.0, b = 2.0;
    double h[4] = {0.0, }, local_a[4] = {0.0, }, local_b[4] = {0.0, };
    int local_n[4] = {0, };

    double local_int[4] = {0.0, }, total_int[4] = {0.0, };
    double start[4] = {0.0, }, finish[4] = {0.0, }, elapsed[4] = {0.0, };

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);


    for(int i = 0; i < 4; i++){
        // Start timing the entire program
        start[i] = MPI_Wtime();

        h[i] = (b - a) / n[i]; // width of each trapezoid
        local_n[i] = n[i] / comm_sz; // number of trapezoids for each process

        // Each process calculates its local interval
        local_a[i] = a + my_rank * local_n[i] * h[i];
        local_b[i] = local_a[i] + local_n[i] * h[i];

        // Calculate the local integral
        
        local_int[i] = Trap(local_a[i], local_b[i], local_n[i], h[i]);
        // local_int[i] = Trap_2(local_a[i], local_b[i], local_n[i], h[i]);
            
        // Sum up the integrals from all processes using MPI_Reduce
        MPI_Reduce(&local_int[i], &total_int[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        // End timing the entire program
        finish[i] = MPI_Wtime();
    } 


    // Output the result and total time
    if (my_rank == 0) {
        for(int i = 0; i < 4; i++){    
            printf("==================================================================================================\n\n");
            printf("Part 2: Scalability Test\n");
            printf("%d Processes used\n", comm_sz);
            printf("With n = %d trapezoids, our estimate of the integral from %f to %f = %.15f\n", n[i], a, b, total_int[i]);
            printf("Total elapsed time = %f seconds\n\n", elapsed[i]);
            printf("==================================================================================================\n\n");
        }
    }

    MPI_Finalize();
    return 0;
}
