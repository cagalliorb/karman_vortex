#define _USE_MATH_DEFINES
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

const double Re = 100.0;  // Reynolds Number
const double cfl = 0.2;   // CFL Number

/*SOR Pamameters*/
const double omegap = 1.00;
const int maxitp = 100;
const double errorp = 0.0001;

/* No. of Time Steps*/
const int nlast = 200000;

/* set x-grid parameters*/
const int mx = 401;  // x軸格子点数(1~401),x=30.0がmx=401に対応
const double r = 60;

/* set y-grid parameters*/
const int my = 201;   // y軸格子点数(1~201)
const int j_1 = 96;   // y軸格子点数基準での,角柱の下端
const int j_2 = 106;  // y軸格子点数基準での,角柱の上端

/* set delta xi,theta,t*/
const double dxi = std::log(r) / (mx - 1);
const double dtheta = 2 * M_PI / (my - 1);
const double dt = cfl * fmin(dxi, dtheta);
const double u_inf = cfl * fmin(dxi, dtheta) / dt;
/*配列の定義*/
double xi[mx + 1], theta[my + 1];
double omega[mx + 2][my + 2], phi[mx + 2][my + 2];
double jcent = (j_1 + j_2) / 2;

void make_xygrid(double xi[mx + 1], double theta[my + 1]) {
    for (int i = 1; i <= mx; i++) {
	xi[i] = dxi * double(i);
    }
    for (int j = 1; j <= my; j++) {
	theta[j] = dtheta * double(j - jcent);
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
	if (j <= float(my) / 4 || j >= float(my) / 4 * 3) {
	    phi[mx][j] = u_inf * r * std::sin(j * dtheta - M_PI);
	} else {
	    phi[mx][j] = phi[mx - 1][j];
	}
	phi[1][j] = 0;
    }
    /*
	for (int i = 1; i <= mx; i++) {
	    phi[i][1] = phi[i][2];
	    phi[i][my] = phi[i][my - 1];
	}*/
}
void bcfor_omega(double omega[mx + 2][my + 2]) {
    for (int j = 1; j <= my; j++) {
	if (j <= float(my) / 4 || j >= float(my) / 4 * 3) {
	    omega[mx][j] = 0;
	} else {
	    omega[mx][j] = omega[mx - 1][j];
	}
	omega[1][j] = -2 * phi[2][j] * std::exp(-2 * xi[1]) / (dxi * dxi);
    }
    /*
for (int i = 1; i <= mx; i++) {
    omega[i][1] = omega[i][2];
    omega[i][my] = omega[i][my - 1];
}*/
}
void poisson_eq(double phi[mx + 2][my + 2], double omega[mx + 2][my + 2]) {
    for (int itr = 1; itr <= maxitp; itr++) {
	double res = 0.0;
	for (int i = 2; i <= mx - 1; i++) {
	    for (int j = 2; j <= my - 1; j++) {
		double phi_new =
		    (phi[i + 1][j] + phi[i - 1][j]) / (dxi * dxi) +
		    (phi[i][j + 1] + phi[i][j - 1]) / (dtheta * dtheta) +
		    std::exp(2 * xi[i]) * omega[i][j];
		phi_new = phi_new / (2 / (dxi * dxi) + 2 / (dtheta * dtheta));
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
	    double temp1 = phi[i][j + 1] - phi[i][j - 1];
	    double temp2 = omega[i + 1][j] - omega[i - 1][j];
	    double rhs1 = -temp1 * temp2 / (4 * dxi * dtheta);

	    double temp3 = phi[i + 1][j] - phi[i - 1][j];
	    double temp4 = omega[i][j + 1] - omega[i][j - 1];
	    double rhs2 = temp3 * temp4 / (4 * dxi * dtheta);

	    double temp5 = omega[i + 1][j] + omega[i - 1][j] - 2 * omega[i][j];
	    double temp6 = omega[i][j + 1] + omega[i][j - 1] - 2 * omega[i][j];
	    double rhs3 =
		(temp5 / (dxi * dxi) + temp6 / (dtheta * dtheta)) / Re;

	    omega_new[i][j] =
		omega[i][j] + dt * std::exp(-2 * xi[i]) * (rhs1 + rhs2 + rhs3);
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

int file_write_xitheta(double xi[mx + 1], double theta[mx + 1]) {
    FILE *fp;
    fp = fopen("xitheta.csv", "a");
    for (int i = 1; i <= mx; i++) {
	for (int j = 1; j <= my; j++) {
	    fprintf(fp, "%f,%f\n", xi[i], theta[j]);
	}
    }

    fclose(fp);
    return 0;
}

int main() {
    make_xygrid(xi, theta);
    // file_write_xitheta(xi, theta);
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
	if (n % 1000 == 0) {
	    std::string s = "./data80/" + std::to_string(n) + "data.csv";
	    const char *cs = s.data();
	    file_write_val(omega, cs);
	}
    }
    return 0;
}
