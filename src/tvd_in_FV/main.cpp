#include <iostream>
#include <vector>
#include <cmath>

#define Nx 41
#define Ny 41
#define Nt 400

#define POSITIVE
//#undef POSITIVE
std::vector<double> dxy = {1., 1.};
std::vector<std::vector<std::vector<double>>> vel(2, std::vector<std::vector<double>>(Nx, std::vector<double>(Ny, 0.)));
std::vector<std::vector<double>> vel_edg_x(Nx - 1, std::vector<double>(Ny    , 0.));
std::vector<std::vector<double>> vel_edg_y(Nx    , std::vector<double>(Ny - 1, 0.));
#ifdef POSITIVE
std::vector<double> constv = {2., 0.};
#else
std::vector<double> constv = {-2., 0.};
#endif
int ix_center = Nx / 2, iy_center = Ny / 2;
std::vector<std::vector<double>> field    (Nx, std::vector<double>(Ny, 0.));
std::vector<std::vector<double>> field_new(Nx, std::vector<double>(Ny, 0.));
std::vector<std::vector<double>> R_plus (Nx, std::vector<double>(Ny, 0.));
std::vector<std::vector<double>> R_minus(Nx, std::vector<double>(Ny, 0.));
std::vector<std::vector<double>> Phi_x(Nx - 1, std::vector<double>(Ny    , 0.));
std::vector<std::vector<double>> Phi_y(Nx    , std::vector<double>(Ny - 1, 0.));
std::vector<std::vector<std::vector<double>>> K (Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(5, 0.)));
std::vector<std::vector<std::vector<double>>> D (Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(5, 0.)));
std::vector<std::vector<std::vector<double>>> F (Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(5, 0.)));
std::vector<std::vector<std::vector<double>>> Mx(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(5, 0.)));
std::vector<int> biases_x = {0, -1, 0, 1, 0};
std::vector<int> biases_y = {-1, 0, 0, 0, 1};

double dt = 0.1;

//#ifdef POSITIVE
//int ix_start = 4;
//int ix_end = 7;
//#else
//int ix_start = Nx - 1 - 7;
//int ix_end = Nx - 1 - 4;
//#endif
int iy_start = 4;
int iy_end   = 7;

int ix_start = Nx / 2 - 3;
int ix_end   = Nx / 2 + 3;
//int iy_start = 0;
//int iy_end   = Ny - 1;

int it = 0;

double tvd_limit(const double r){
	
	if (r < 0) return 0.;
	return fmax(fmin(1.0, 2 * r), fmin(2.0, r));
}

double get_r_i(int ix, int iy, bool positive_stream, int ort);
double get_r_i_biased(int ix, int iy, int ort);
void fill_R(const int ort);
void fill_F(void);
void fill_Phi_xy(void);

int main(){
	
	for(int ort = 0; ort < 2; ort++)
	for(int ix = 0; ix < Nx; ix++)
	for(int iy = 0; iy < Ny; iy++){
		vel[ort][ix][iy] = constv[ort];
		if(ix < Nx / 2)
			vel[ort][ix][iy] = -constv[ort];
		else if(ix == Nx / 2)
			vel[ort][ix][iy] = 0.;
		else
			vel[ort][ix][iy] = constv[ort];
//#ifdef POSITIVE
//		vel[ort][ix][iy] = constv[ort] * fmin(3., 3. * (double) (Nx - 1 + 5 - ix) / (double) Nx);
////		vel[ort][ix][iy] = constv[ort] * fmin(3., 3. * (double) (5 + ix) / (double) Nx);
//#else
////		vel[ort][ix][iy] = constv[ort] * fmin(3., 3. * (double) (Nx - 1 + 5 - ix) / (double) Nx);
//		vel[ort][ix][iy] = constv[ort] * fmin(3., 3. * (double) (5 + ix) / (double) Nx);
//#endif
	}
// 	//Круговое поле скоростей
// 	for(int ix = 0; ix < Nx; ix++)
// 	for(int iy = 0; iy < Ny; iy++){
// 		
// 		vel[0][ix][iy] = - constv[0] * (iy - iy_center) / (double)Ny;
// 		vel[1][ix][iy] =   constv[1] * (ix - ix_center) / (double)Nx;
// 	}
	
	for(int ix = 0; ix < Nx - 1; ix++)
	for(int iy = 0; iy < Ny; iy++){
		vel_edg_x[ix][iy] = (vel[0][ix][iy] + vel[0][ix + 1][iy]) * 0.5;
	}
	
	for(int ix = 0; ix < Nx; ix++)
	for(int iy = 0; iy < Ny - 1; iy++){
		vel_edg_y[ix][iy] = (vel[1][ix][iy] + vel[1][ix][iy + 1]) * 0.5;
	}
	
	for(int ix = 0; ix < Nx; ix++)
	for(int iy = 0; iy < Ny; iy++){
		
		if(ix >= ix_start && ix <= ix_end &&
		   iy >= iy_start && iy <= iy_end)
			field[ix][iy] = 1.;
		else
			field[ix][iy] = 0.;
	}
	
	{
		char str[20]; sprintf(str, "vel_field.txt");
		FILE *res = fopen(str, "w");
		for(int iy = 0; iy < Ny - 1; iy++){
			for(int ix = 0; ix < Nx - 1; ix++){
				double x = (double) ix * dxy[0] + dxy[0] * 0.5;
				double y = (double) iy * dxy[1] + dxy[1] * 0.5;
				fprintf(res, "%12.6e %12.6e %12.6e %12.6e\n", x, y, vel_edg_x[ix][iy], vel_edg_y[ix][iy]);
			}
		}
		fclose(res);
	}
	
	for(it = 0; it < Nt; it++){
		
		// Fill matrix Mx
		fill_F();
		for(int ix = 0; ix < Nx; ix++)
		for(int iy = 0; iy < Ny; iy++)
		for(int adj = 0; adj < 5; adj++)
			Mx[ix][iy][adj] = F[ix][iy][adj];
		
		// Solve
//		for(int ix = 0; ix < Nx; ix++)
//		for(int iy = 0; iy < Ny; iy++)
//			if(iy == iy_start && ix >= ix_start && ix <= ix_end)
//				field[ix][iy] += 1 * dt;
		
		for(int ix = 0; ix < Nx; ix++)
		for(int iy = 0; iy < Ny; iy++){
			field_new[ix][iy] = field[ix][iy];
			for(int adj = 0; adj < 5; adj++){
				if(ix + biases_x[adj] < 0 || ix + biases_x[adj] > Nx - 1 ||
				   iy + biases_y[adj] < 0 || iy + biases_y[adj] > Ny - 1)
					continue;
				field_new[ix][iy] -= dt * Mx[ix][iy][adj] * field[ix + biases_x[adj]][iy + biases_y[adj]];
			}
		}
		
		double field_sum = 0.;
		for(int ix = 0; ix < Nx; ix++)
		for(int iy = 0; iy < Ny; iy++){
			field[ix][iy] = field_new[ix][iy];
			field_sum += field[ix][iy];
		}
		printf("%d: %12.6e\n", it, field_sum);
		
		char str[20]; sprintf(str, "%04d.txt", it);
		FILE *res = fopen(str, "w");
		for(int iy = 0; iy < Ny; iy++){
			for(int ix = 0; ix < Nx; ix++){
				fprintf(res, "%04d %04d %12.6e\n", ix, iy, field_new[ix][iy]);
			}
			fprintf(res, "\n");
		}
		fclose(res);
		
		{
			sprintf(str, "colsum_%04d.txt", it);
			res = fopen(str, "w");
			
			std::vector<std::vector<double> > colsums(Nx, std::vector<double>(Ny, 0.));
			for(int iy = 0; iy < Ny; iy++){
				for(int ix = 0; ix < Nx; ix++){
					for(int adj = 0; adj < 5; adj++){
						int ix_biased = ix + biases_x[adj];
						int iy_biased = iy + biases_y[adj];
						if(ix_biased < 0 || ix_biased > Nx - 1 || iy_biased < 0 || iy_biased > Ny - 1)
							continue;
						colsums[ix_biased][iy_biased] += F[ix][iy][adj];
					}
				}
			}
			for(int iy = 0; iy < Ny; iy++){
				for(int ix = 0; ix < Nx; ix++)
					fprintf(res, "%04d %04d %12.6e\n", ix, iy, colsums[ix][iy]);
			}
			fclose(res);
		}
		
		{
			sprintf(str, "xline_%04d.txt", it);
			res = fopen(str, "w");
			int iy = Ny / 2;
			for(int ix = 0; ix < Nx; ix++){
				fprintf(res, "%04d %04d %12.6e %12.6e %12.6e\n", ix, iy, field_new[ix][iy],
				        tvd_limit(R_plus[ix][iy]), tvd_limit(R_minus[ix][iy]));
			}
			fclose(res);
		}
		
		{
			sprintf(str, "yline_%04d.txt", it);
			res = fopen(str, "w");
			int ix = Nx / 2;
			for(int iy = 0; iy < Ny; iy++){
				fprintf(res, "%04d %04d %12.6e %12.6e %12.6e\n", ix, iy, field_new[ix][iy],
				        tvd_limit(R_plus[ix][iy]), tvd_limit(R_minus[ix][iy]));
			}
			fclose(res);
		}
	}
	
	return 0;
}

void fill_R(const int ort){
	
	for(int ix = 0; ix < Nx; ix++){
		for(int iy = 0; iy < Ny; iy++){
			R_plus [ix][iy] = 0;
			R_minus[ix][iy] = 0.;
		}
	}
	
	for(int ix = 0; ix < Nx; ix++){
		for(int iy = 0; iy < Ny; iy++){
			
			R_plus [ix][iy] = get_r_i(ix, iy, true,  ort);
			R_minus[ix][iy] = get_r_i(ix, iy, false, ort);
		}
	}
	
	return;
}

void fill_F(void){
	
	for(int ix = 0; ix < Nx; ix++)
		for(int iy = 0; iy < Ny; iy++)
			for(int adj = 0; adj < 5; adj++)
				F[ix][iy][adj] = 0;
	
	fill_Phi_xy();
	
	std::vector<int> Nxy = {Nx, Ny};
	
	for(int ort = 0; ort < 2; ort++){
		// смотрим потоки вдоль смежного направления
		int tro = (ort + 1) % 2;
		const double yxd = 1 / dxy[tro];
		const double yxd2 = 0.5 * yxd;
		
		char str[100];
		if(ort == 0)
			sprintf(str, "yvals_%04d.txt", it);
		else
			sprintf(str, "xvals_%04d.txt", it);
		FILE *vals = fopen(str, "w");
		fprintf(vals, "#%13s %14s %14s %14s\n", "ix", "iy", "ri", "phi_term");
//		fill_R(tro);
		// поле скоростей вдоль нужного напрвления
		std::vector<std::vector<double> > *v1 = (tro == 0) ? &vel_edg_x : &vel_edg_y;
		
		std::vector<std::vector<double> > *Phi1 = (tro == 0) ? &Phi_x : &Phi_y;
		// смещение вдоль выбранного направления
		int adj_prev = (tro == 0) ? 1 : 0;
		int adj_next = (tro == 0) ? 3 : 4;
		
		for(int i1 = 0; i1 < Nxy[ort]; i1++){
			
			// расчет переноса через ребра выбранного направления
			for(int i2 = 0; i2 < Nxy[tro] - 1; i2++){
				
				// рассматривается ребро (i1, i2) между ячейками (i1, i2) и (i1, i2+1)
				std::vector<int> ixy_prev(2); ixy_prev[ort] = i1; ixy_prev[tro] = i2;
				std::vector<int> ixy_next(2); ixy_next[ort] = i1; ixy_next[tro] = i2 + 1;
				
				double vel_xy = (*v1)[ixy_prev[0]][ixy_prev[1]] / dxy[tro];
				
				double phi_term = 0.5 * fabs(vel_xy) * (*Phi1)[ixy_prev[0]][ixy_prev[1]];
				
				if(vel_xy < 0){
					F[ixy_prev[0]][ixy_prev[1]][2]        += - phi_term;
					F[ixy_prev[0]][ixy_prev[1]][adj_next] +=   phi_term + vel_xy;
					
					F[ixy_next[0]][ixy_next[1]][adj_prev] +=   phi_term;
					F[ixy_next[0]][ixy_next[1]][2]        += - phi_term - vel_xy;
				} else{
					F[ixy_prev[0]][ixy_prev[1]][2]        += - phi_term + vel_xy;
					F[ixy_prev[0]][ixy_prev[1]][adj_next] +=   phi_term;
					
					F[ixy_next[0]][ixy_next[1]][adj_prev] +=   phi_term - vel_xy;
					F[ixy_next[0]][ixy_next[1]][2]        += - phi_term;
				}
				fprintf(vals, "%14d %14d %14e %14e\n", ixy_prev[0], ixy_prev[1],
				        get_r_i_biased(ixy_prev[0], ixy_prev[1], tro), phi_term);
			}
			fprintf(vals, "\n");
		}
		fclose(vals);
	}
}

double get_r_i(int ix, int iy, bool positive_stream, int ort){
	
	double ri_curr;
	double field_up, field_dn, field_cr;
	int ixy_prev[2] = {ix, iy}; ixy_prev[ort]--;
	int ixy_next[2] = {ix, iy}; ixy_next[ort]++;
	
	int Nxy[2] = {Nx, Ny};
	if(ixy_prev[ort] < 0)
		ixy_prev[ort] = 0;
	if(ixy_next[ort] > Nxy[ort] - 1)
		ixy_next[ort] =  Nxy[ort] - 1;
	
	if(positive_stream)
	{
		field_up = field[ixy_next[0]][ixy_next[1]];
		field_dn = field[ixy_prev[0]][ixy_prev[1]];
	} else
	{
		field_up = field[ixy_prev[0]][ixy_prev[1]];
		field_dn = field[ixy_next[0]][ixy_next[1]];;
	}
	
	field_cr = field[ix][iy];
	
	if(field_cr == field_up){
		if(field_cr - field_dn < 0.)
			ri_curr = -1e6;
		else if(field_cr - field_dn > 0.)
			ri_curr = 1e6;
		else
			ri_curr = 1.;
	}
	else{
		ri_curr = (field_cr - field_dn) / (field_up - field_cr);
		if(std::isinf(ri_curr) || std::isnan(ri_curr))
			ri_curr = 1;
	}
	return ri_curr;
}

double get_r_i_biased(int ix, int iy, int ort){
	
	int tro = (ort + 1) % 2;
	double ri_curr;
	double field_up, field_dn, field_cr;
	int Nxy[2] = {Nx, Ny};
	int ixy[2] = {ix, iy};
	
	double vel = (ort == 0) ? vel_edg_x[ix][iy] : vel_edg_y[ix][iy];
	
	if((vel > 0 && ixy[ort] != 0) || ixy[ort] == Nxy[ort] - 2)
	{
		field_cr = field[ix      ][iy      ];
		field_up = field[ix + tro][iy + ort];
		field_dn = field[ix - tro][iy - ort];
	} else
	{
		field_cr = field[ix + tro    ][iy + ort    ];
		field_up = field[ix          ][iy          ];
		field_dn = field[ix + 2 * tro][iy + 2 * ort];
	}
	
	// TODO: сравнивать до заданной точности
	if(field_cr == field_up){
		if(field_cr - field_dn < 0.)
			ri_curr = -1e6;
		else if(field_cr - field_dn > 0.)
			ri_curr = 1e6;
		else
			ri_curr = 1.;
	}
	else{
		ri_curr = (field_cr - field_dn) / (field_up - field_cr);
		if(std::isinf(ri_curr) || std::isnan(ri_curr))
			ri_curr = 1;
	}
	return ri_curr;
}

void fill_Phi_xy(void){
	
	///////////////////////////////
	// x
	///////////////////////////////
	for(unsigned int iy = 0; iy < Ny; iy++){
		for(unsigned int ix = 0; ix < Nx - 1; ix++){
			double r = get_r_i_biased(ix, iy, 0);
			Phi_x[ix][iy] = tvd_limit(r);
		}
	}
	
	///////////////////////////////
	// y
	///////////////////////////////
	for(unsigned int iy = 0; iy < Ny - 1; iy++){
		for(unsigned int ix = 0; ix < Nx; ix++){
			double r = get_r_i_biased(ix, iy, 1);
			Phi_y[ix][iy] = tvd_limit(r);
		}
	}
}
