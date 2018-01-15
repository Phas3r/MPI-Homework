# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <mpi.h>
# include <string>

using namespace std;

int main(int argc, char *argv[]);
double f(double x);
double integrate(double a, double b, int n, int method = 1);

int main(int argc, char *argv[])
{
    double a, b;
    double error;
    double my_a, my_b; // intergrating limits for each process
    int my_n;
    double total, my_total, wtime;
    int p;// total number of processes
    int source;
    MPI_Status status;
    int tag, target, id;
    double x;
    int master = 0;
    int method;


    a = 0.0;
    b = 1.0;

    int n = 10000000;
    double exact = 3.1415926535897932384626433832795028841971693993751;
    //
    //  Initialize MPI.
    //
    MPI_Init(&argc, &argv);
    //
    //  Get this processor's ID.
    //
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    //
    //  Get the number of processors.
    //
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (id == 0)
    {
        my_n = n / (p);
        n = (p) * my_n;
        wtime = MPI_Wtime();

        if (argc > 1)
        {
            try
            {
                method = std::stoul(argv[1]);
            }
            catch (std::invalid_argument e) {
                method = 1;
            }
            if (argc > 2)
            {
                try
                {
                    b = std::stod(argv[2]);
                }
                catch (std::invalid_argument e) {
                    b = 1.0;
                }
            }
        }
    }

    source = 0;
    MPI_Bcast(&my_n, 1, MPI_INT, source, MPI_COMM_WORLD);
    //
    //  Process 0 assigns each process a subinterval of [A,B].
    //
    if (id == 0)
    {
        for (int q = 1; q <= p - 1; q++)
        {
            my_a = a + ((b - a) / double(p))*double(q);

            target = q;//process to which we send arguments
            tag = 1;
            MPI_Send(&my_a, 1, MPI_DOUBLE, target, tag, MPI_COMM_WORLD);

            my_b = a + ((b - a) / double(p))*double(q+1);

            target = q;
            tag = 2;
            MPI_Send(&my_b, 1, MPI_DOUBLE, target, tag, MPI_COMM_WORLD);
        }
        my_a = a;
        my_b = a + (b - a) / p;
    }
    //
    //  Processes receive MY_A, MY_B, and compute their part of the integral.
    //
    else
    {
        source = 0;
        tag = 1;
        MPI_Recv(&my_a, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);

        source = 0;
        tag = 2;
        MPI_Recv(&my_b, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
    }

    my_total = integrate(my_a, my_b, my_n, method);
    
    //
    //  Each process sends its value to the master process.
    //
    MPI_Reduce(&my_total, &total, 1, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    //
    //  Compute the weighted estimate.
    //
    if (id == 0)
    {
        error = fabs(total - exact);
        wtime = MPI_Wtime() - wtime;

        /*cout << "\n";
        cout << "  Estimate = " << setw(24) << setprecision(16) << total << "\n";
        cout << "  Error = " << error << "\n";*/
        cout << p << ", " << wtime << "\n";
    }
    //
    //  Terminate MPI.
    //
    MPI_Finalize();
    //
    //  Terminate.
    //
    /*if (id == 0)
    {
        cout << "\n";
        cout << "QUAD_MPI:\n";
        cout << "  Normal end of execution.\n";
        cout << "\n";
    }*/
    return 0;
}
//****************************************************************************80

double f(double x)
{
    double value;

    value = 4.0 / (x * x + 1.0);

    return value;
}

double rectangular(double a, double b, int n)
{
    double h = (b - a) / n;
    double x = a + h/2, total = 0.0;
    for (int i = 1; i <= n; i++)
    {
        total += f(x);
        x += h;
    }
    total = (b - a) * total / (double)(n);
    return total;
}

double trapezoidal(double a, double b, int n)
{
    double h = (b - a) / n;
    double x = a, total = 0.0;
    double value = f(x);
    for (int i = 1; i <= n; i++)
    {
        double prev_value = value;
        x += h;
        value = f(x);
        total += (prev_value + value)/2;
    }
    total = (b - a) * total / (double)(n);
    return total;
}

double simpson(double a, double b, int n)
{
    double h = (b - a) / n;
    double h_2 = h / 2;
    double x = a, total = 0.0;
    double value = f(x);
    for (int i = 1; i <= n; i++)
    {
        double prev_value = value;
        x += h;
        value = f(x);
        total += (prev_value + f(x - h_2) * 4 + value) / 6;
    }
    total = (b - a) * total / (double)(n);
    return total;
}

double gaussian(double a, double b, int n)
{
    double h = (b - a) / n;
    double h_2 = h / 2;
    double const1 = h_2 - h_2 / sqrt(3);
    double const2 = h_2 + h_2 / sqrt(3);
    double x = a, total = 0.0;;
    for (int i = 1; i <= n; i++)
    {
        total += (f(x + const1) + f(x+const2));
        x += h;
    }
    total /= 2;
    total = (b - a) * total / (double)(n);
    return total;
}


double integrate(double a, double b, int n, int method)
{
    switch (method)
    {
    case 1:
        return rectangular(a, b, n);
    case 2:
        return trapezoidal(a, b, n);
    case 3:
        return simpson(a, b, n);
    case 4:
        return gaussian(a, b, n);
    default:
        return rectangular(a, b, n);
    }
}