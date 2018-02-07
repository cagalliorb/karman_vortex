#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

const double Re = 70.0;   // Reynolds Number
const double cfl = 0.02;  // CFL Number

/*SOR Pamameters*/
const double omegap = 1.00;
const int maxitp = 50;
const double errorp = 0.0001;

/* No. of Time Steps*/
const int nlast = 50000;

/* set x-grid parameters*/
const int mx = 401;   // x軸格子点数(1~401),x=30.0がmx=401に対応
const int i_1 = 96;   // x軸格子点数基準での,角柱の左端
const int i_2 = 106;  // x軸格子点数基準での,角柱の右端

/* set y-grid parameters*/
const int my = 201;   // y軸格子点数(1~201)
const int j_1 = 96;   // y軸格子点数基準での,角柱の下端
const int j_2 = 106;  // y軸格子点数基準での,角柱の上端

/* set delta x,y,t*/
const double dx = 1.0 / (i_2 - i_1);
const double dy = 1.0 / (j_2 - j_1);
const double dt = cfl * fmin(dx, dy);
const double u_inf = cfl * dx / dt;
/*配列の定義*/
double x[mx + 1];
double y[my + 1];
double omega[mx + 2][my + 2], phi[mx + 2][my + 2];
double icent = (i_1 + i_2) / 2;
double jcent = (j_1 + j_2) / 2;

void make_xygrid(double x[mx + 1], double y[my + 1]) {
    for (int i = 1; i <= mx; i++) {
	x[i] = dx * double(i - icent);
    }
    for (int j = 1; j <= my; j++) {
	y[j] = dy * double(j - jcent);
    }
}

void init_condition(double phi[mx + 2][my + 2], double omega[mx + 2][my + 2]) {
    for (int i = 1; i <= mx; i++) {
	for (int j = 1; j <= my; j++) {
	    omega[i][j] = 0.0;
	    phi[i][j] = 0.0;
	}
    }
}

void bcfor_phi(double phi[mx + 2][my + 2]) {
    for (int j = 1; j <= my; j++) {
	phi[1][j] = u_inf * y[j];
	phi[mx][j] = phi[mx - 1][j];
    }
    for (int i = 1; i <= mx; i++) {
	phi[i][1] = phi[1][1];
	phi[i][my] = phi[1][my];
    }

    for (int j = j_1; j <= j_2; j++) {
	phi[i_1][j] = 0;
	phi[i_2][j] = 0;
    }
    for (int i = i_1; i <= i_2; i++) {
	phi[i][j_1] = 0;
	phi[i][j_2] = 0;
    }
}

void bcfor_omega(double omega[mx + 2][my + 2]) {
    for (int j = 1; j <= my; j++) {
	omega[1][j] = 0;
	omega[mx][j] = omega[mx - 1][j];
    }
    for (int i = 1; i <= mx; i++) {
	omega[i][1] = -2 * (phi[i][2] - phi[i][1] - u_inf * dy) / (dy * dy);
	omega[i][my] =
	    -2 * (phi[i][my - 1] - phi[i][my] + u_inf * dy) / (dy * dy);
    }

    for (int i = i_1; i <= i_2; i++) {
	omega[i][j_1] = -2 * phi[i][j_1 - 1] / (dy * dy);
	omega[i][j_2] = -2 * phi[i][j_2 + 1] / (dy * dy);
    }
    for (int j = j_1; j <= j_2; j++) {
	omega[i_1][j] = -2 * phi[i_1 - 1][j] / (dx * dx);
	omega[i_2][j] = -2 * phi[i_2 + 1][j] / (dx * dx);
    }
}
void poisson_eq(double phi[mx + 2][my + 2], double omega[mx + 2][my + 2]) {
    for (int itr = 1; itr <= maxitp; itr++) {
	double res = 0.0;
	for (int i = 2; i <= mx - 1; i++) {
	    for (int j = 2; j <= my - 1; j++) {
		if (i_1 < i && i < i_2 && j_1 < j && j < j_2) {
		    continue;
		}
		/*double dphi = (phi[i + 1][j] + phi[i - 1][j]) / (dx * dx) +
			      (phi[i][j + 1] + phi[i][j - 1]) / (dy * dy) +
			      omega[i][j];*/
		double phi_new = (phi[i + 1][j] + phi[i - 1][j] +
				  phi[i][j + 1] + phi[i][j - 1]) /
				     4 +
				 dx * dx / 4 * omega[i][j];

		double dphi = phi_new - phi[i][j];
		res = res + dphi * dphi;
		phi[i][j] = phi_new;
	    }
	}
	bcfor_phi(phi);
	res = sqrt(res / double(mx * my));
	if (res < errorp) break;
    }
}

void voriticity_eq(double phi[mx + 2][my + 2], double omega[mx + 2][my + 2]) {
    double omega_new[mx + 2][my + 2];
    for (int i = 2; i <= mx - 1; i++) {
	for (int j = 2; j <= my - 1; j++) {
	    if (i_1 < i && i < i_2 && j_1 < j && j < j_2) {
		continue;
	    }
	    double temp1 = phi[i][j + 1] - phi[i][j - 1];
	    double temp2 = omega[i + 1][j] - omega[i - 1][j];
	    double rhs1 = -temp1 * temp2 / (4 * dx * dy);

	    double temp3 = phi[i + 1][j] - phi[i - 1][j];
	    double temp4 = omega[i][j + 1] - omega[i][j - 1];
	    double rhs2 = temp3 * temp4 / (4 * dx * dy);

	    double temp5 = omega[i + 1][j] + omega[i - 1][j] - 2 * omega[i][j];
	    double temp6 = omega[i][j + 1] + omega[i][j - 1] - 2 * omega[i][j];
	    double rhs3 = (temp5 / (dx * dx) + temp6 / (dy * dy)) / Re;

	    omega_new[i][j] = omega[i][j] + dt * (rhs1 + rhs2 + rhs3);
	}
    }
    for (int i = 2; i <= mx - 1; i++) {
	for (int j = 2; j <= my - 1; j++) {
	    omega[i][j] = omega_new[i][j];
	}
    }
}

int file_write_val(double val[mx + 2][my + 2], const char *file_name) {
    FILE *fp;
    fp = fopen(file_name, "w");
    for (int i = my; i >= 1; i--) {
	for (int j = 1; j <= mx; j++) {
	    fprintf(fp, "%f", val[j][i]);
	    if (j == mx) {
		fprintf(fp, "\n");
	    } else {
		fprintf(fp, ",");
	    }
	}
    }
    fclose(fp);
    return 0;
}
int main() {
    make_xygrid(x, y);
    init_condition(phi, omega);
    bcfor_phi(phi);
    bcfor_omega(omega);
    std::chrono::system_clock::time_point start_time,
	end_time;  //時間計測用変数を確保
    start_time = std::chrono::system_clock::now();  //計測スタート:
    for (int n = 1; n <= nlast; n++) {
	voriticity_eq(phi, omega);
	bcfor_omega(omega);
	poisson_eq(phi, omega);
	bcfor_phi(phi);
	end_time = std::chrono::system_clock::now();  //
	double elapsed_time =
	    std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
								  start_time)
		.count() /
	    1000.0;
	std::cout << "step:" << n << "/" << nlast
		  << " elapsed_time:" << elapsed_time << "[s]" << std::endl;
	std::cout << "omega[200][200]= " << omega[200][100] << std::endl;
	std::string s = "./data30/" + std::to_string(n) + "data.csv";
	const char *cs = s.data();
	file_write_val(omega, cs);
    }
    return 0;
}
