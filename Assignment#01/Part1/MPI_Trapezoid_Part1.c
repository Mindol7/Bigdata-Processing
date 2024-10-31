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
    int my_rank, comm_sz, n = 4096; // 4096 is number of Trapezoid
    double a = 0.0, b = 2.0, h, local_a, local_b;
    int local_n;
    double local_int, total_int;
    double start, finish, elapsed;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    // Start timing the entire program
    start = MPI_Wtime();

    h = (b - a) / n;  // width of each trapezoid
    local_n = n / comm_sz;  // number of trapezoids for each process

    // Each process calculates its local interval
    local_a = a + my_rank * local_n * h;
    local_b = local_a + local_n * h;
    
    // Calculate the local integral
    // local_int = Trap(local_a, local_b, local_n, h);
    local_int = Trap_2(local_a, local_b, local_n, h);

    //PI Sum up the integrals from all processes using MPI_Reduce
    MPI_Reduce(&local_int, &total_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // End timing the entire program
    finish = MPI_Wtime();
    elapsed = finish - start;

    // Output the result and total time
    if (my_rank == 0) {
        printf("==================================================================================================\n\n");
        printf("Part 1: Serial vs Parallel Execution Time Comparasion\n");
        printf("%d Processes used\n", comm_sz);
        printf("With n = %d trapezoids, our estimate of the integral from %f to %f = %.15e\n", n, a, b, total_int);
        printf("Total elapsed time = %f seconds\n\n", elapsed);
        printf("==================================================================================================\n\n");
    }

    MPI_Finalize();
    return 0;
}
