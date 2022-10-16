#ifndef DEF
#define DEF

#include <iostream>
#include <math.h>
#include <SFML\\Graphics.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <omp.h>
#include <time.h> 
#include <chrono>
#include <amp.h>


using namespace std;
using namespace concurrency;


double radians(double a) {
	return (a / 180 * 3.141592653);
}
double to_double(int a) {
	return a;
}
int to_int(double a) {
	return a;
}

int sfUint8_to_int(sf::Uint8 uint8_num) {
	return uint8_num * 1;
}

class Engine3D {

public:
	int rx, ry, RESOLUTION_X, RESOLUTION_Y;
	double x, y, angle;
	vector<vector<vector<char>>> render_matrix;
	vector<vector<vector<float>>> rays_steps;
	vector<vector<int>> map;
	vector<vector<vector<vector<char>>>> textures;
	vector<vector<double>> heights;
	int d_angle;
	int integer_x, integer_y, code, page, map_width, map_height;
	
	sf::Uint8* bg_main;
	sf::Texture bg_texture;
	sf::Image bg_image;
	double x_past, y_past;
	double ray_step;
	double current_x, current_y, dx, dy, ox, oy, oz; // текущие x y
	sf::Uint8* pixels;


	void init(double x_, double y_, double angle_, string bg_main_filename, double ray_step_k, int rx_, int ry_) {
		RESOLUTION_X = rx_;
		RESOLUTION_Y = ry_;
		rx = rx_;
		ry = ry_;
		//render_matrix = render_matrix(360, vector<vector<int>>(rx, vector<int>(2)));
		//rays_steps = rays_steps(360, vector<vector<double>>(rx, vector<double>(2)));
		const sf::Uint8* bg_const = new sf::Uint8[RESOLUTION_X * RESOLUTION_Y * 4 + 4];
		pixels = new sf::Uint8[RESOLUTION_X * RESOLUTION_Y * 4 + 4];
		bg_main = new sf::Uint8[RESOLUTION_X * RESOLUTION_Y * 4 + 4];
		page = 1;
		heights.resize(rx, vector<double>(10, 0));
		render_matrix.resize(rx, vector<vector<char>>(ry, vector<char>(3, 0)));
		rays_steps.resize(360, vector<vector<float>>(rx, vector<float>(2, 0)));
		textures.resize(15, vector<vector<vector<char>>>(ry, vector<vector<char>>(ry, vector<char>(3, 0))));
		bg_image.loadFromFile(bg_main_filename);
		bg_const = bg_image.getPixelsPtr();
		int i1;
		for (i1 = 0; i1 < rx * ry * 4; i1++) {
			bg_main[i1] = bg_const[i1];
		}
		x = x_;
		y = y_;
		angle = angle_;
		d_angle = 90.0;
		int i, j, k;
		double rs;
		ray_step = 1.0 / to_double(rx) *  ray_step_k;
		std::cout << "inizialisation rays_steps" << endl;
		//rays_steps[359][127][1] = 2.4;
		std::cout << to_double(abs(rx / 2 - 0)) / to_double(rx) * 10 << endl;
		for (i = 0; i < 360; i++) {
			for (j = 0; j < rx; j++) {
				k = 0;
				rays_steps[i][j][0] = cos(radians(i - (-rx / 2.0 + j) / rx * d_angle)) * ray_step;
				rays_steps[i][j][1] = sin(radians(i - (-rx / 2.0 + j) / rx * d_angle)) * ray_step;
			}
		}
		//cout << rays_steps[0][0][0] << endl;

	}
	void load_texture(int texture_code, string filename) {
		sf::Image texture_img;
		texture_img.loadFromFile(filename);
		const sf::Uint8* pixels1 = new sf::Uint8[RESOLUTION_X * RESOLUTION_Y * 4 + 4];
		sf::Uint8* pixels2 = new sf::Uint8[RESOLUTION_X * RESOLUTION_Y * 4 + 4];

		pixels1 = texture_img.getPixelsPtr();
		int i, j, k;
		for (i = 0;i < ry; i++) {
			for (j = 0; j < ry; j++) {
				textures[texture_code][i][j][0] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 4]);
				textures[texture_code][i][j][1] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 3]);
				textures[texture_code][i][j][2] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 2]);
			}
		}


	}
	void set_map(string filename) {
		std::string line;
		std::ifstream in(filename);
		int str_len, i, k;
		k = 0;
		if (in.is_open())
		{
			cout << "File with MAP is loaded" << endl;
			while (getline(in, line))
			{	
				str_len = line.length();
				k++;
			}
		}
		in.close();
		map.resize(k, vector<int>(str_len, 0));
		int k2;
		k2 = k;
		cout << "width: " << str_len << " height: " << k << endl;
		map_width = str_len;
		map_height = k;
		std::ifstream in2(filename);
		k = 0;
		if (in2.is_open())
		{
			cout << "File with MAP is loaded 2" << endl;
			while (getline(in2, line))
			{
				for (i = 0; i < str_len; i++) {
					//cout << line[i] << " " << static_cast<int>(line[i]) - 48 << endl;
					map[k][i] = static_cast<int>(line[i]) - 48;
					//cout << map[k][i];
				}
				str_len = line.length();
				k++;
			}
		}
		in2.close();
		cout << "File with MAP is loaded 3" << endl;


	}
	void stage_heights() {
		if (page == 1) {
			while (angle >= 360) {
				angle -= 360.0;
			}
			while (angle < 0) {
				angle += 360.0;
			}
			//cout << "stage-heights" << endl; // debug info
			int i, j, k; // итераторы
			int int_angle, int_x, int_y; // целые
			double x_step, y_step, rs;
			int_angle = round(angle);
			
			for (i = 0; i < rx; i++) {
				int_x = trunc(x); // приводим к целому чтобы был индекс массива
				int_y = trunc(y); // приводим к целому чтобы был индекс массива
				k = 0;
				//cout << "x:" << current_x << " y:" << current_y << " a:" << int_angle << " i:" << i << endl;
				current_x = x;
				current_y = y;
					x_step = rays_steps[int_angle][i][0];
					y_step = rays_steps[int_angle][i][1];

				//cout << x << " " << y << ":" << map[0][2] << " " << int_x <<  " " << int_y << endl;

				while ((map[int_x][int_y] == 0) and (k < 10 / ray_step)) { // Проверка не врезалось ли
					current_x += x_step; // делаем шаг x
					current_y += y_step; // делаем шаг y
					int_x = trunc(current_x); // приводим к целому чтобы был индекс массива
					int_y = trunc(current_y); // приводим к целому чтобы был индекс массива
					k++;
				}

				dx = current_x - x;
				dy = current_y - y;
				//cout << cos(radians(abs(int_angle - d_angle / to_double(rx) * i))) << " " << i << " " << d_angle / to_double(rx) << endl;

				heights[i][3] = dx; // x расстояние
				heights[i][4] = dy; // y расстояние
				heights[i][0] = abs((1 / (sqrt(dx*dx + dy*dy) * cos(radians(to_double(i) / to_double(rx) * 90.0 - 45.0)))) * ry / 2.0); // высота
				//cout << heights[i][0] << endl;
				code = map[int_x][int_y]; // код текстуры
				if (k >= 10 / ray_step) {
					code = 10;
				}
				ox = current_x - trunc(current_x); // смещение текстуры при пересчении y - оси
				oy = current_y - trunc(current_y); // смещение текстуры при пересчении x - оси
				if (abs(ox - 0.5) < abs(oy - 0.5)) {
					oz = ox;
				}
				else if (abs(ox - 0.5) > abs(oy - 0.5)) {
					oz = oy;
				}
				heights[i][1] = oz; // нужное смещение
				heights[i][2] = to_double(code); // код текстуры

			}
		}
	}
	sf::Uint8* stage_render(sf::Texture texture_target, sf::Sprite sprite_target) {
		//cout << "stage_render" << endl;

		
		int i, j, k, code, pix0, pix1, pix2, pix02, pix12, pix22, arxf, aryf;
		double oz, i1, i2, i3, i4, distance_for_floor, dx_f, dy_f, ox_f, oy_f;
		if (page == 1) {
#pragma omp parallel for 
			for (i = 0; i < rx; i++) {
				//array_view <const int, 2> avA(ry, 10, heights);
				//array_view <const int, 3> avB(ry, ry, 3, textures);
				//array_view <int, 2> avC(m, n, C);
				//avC.discard_data();
				for (j = 0; j < ry / 2; j++) {

					if ((ry / 2 - j) < heights[i][0]) {
						oz = heights[i][1];
						code = to_int(heights[i][2]);
						i1 = trunc(oz * ry);
						i2 = ry / 2 - ry / 2 / heights[i][0] * (ry / 2 - j);
						i3 = ry / 2 - ry / 2 / heights[i][0] * (j - ry / 2);
						//cout << heights[i][0] / to_double(ry) * 2.0 * (ry / 2 - j - 1) << endl;
						pix0 = textures[code][i2][i1][0];
						pix1 = textures[code][i2][i1][1];
						pix2 = textures[code][i2][i1][2];
						//cout << ry / 2 / heights[i][0] * (ry / 2 + j - 1) << endl;
						pix02 = textures[code][i3][i1][0];
						pix12 = textures[code][i3][i1][1];
						pix22 = textures[code][i3][i1][2];
						pixels[j * rx * 4 + i * 4] = pix0;
						pixels[j * rx * 4 + i * 4 + 1] = pix1;
						pixels[j * rx * 4 + i * 4 + 2] = pix2;
						pixels[(ry - j - 1) * rx * 4 + i * 4] = pix02;
						pixels[(ry - j - 1) * rx * 4 + i * 4 + 1] = pix12;
						pixels[(ry - j - 1) * rx * 4 + i * 4 + 2] = pix22;
						//cout << pix0 << " " <<  pix1 << endl;
					}
					else {
						distance_for_floor = 1.0/(ry/2 - j)*ry/2/ cos(radians(to_double(i) / rx * 90.0 - 45.0));
						dx_f = distance_for_floor * rays_steps[angle][i][0] / ray_step + x;
						dy_f = distance_for_floor * rays_steps[angle][i][1] / ray_step + y;
						ox_f = dx_f - trunc(dx_f);
						oy_f = dy_f - trunc(dy_f);
						arxf = abs(ry * ox_f);
						aryf = abs(ry * oy_f);
						pixels[j * rx * 4 + i * 4] = textures[11][abs(ry * ox_f)][abs(ry * oy_f)][0];
						pixels[j * rx * 4 + i * 4 + 1] = textures[11][abs(ry * ox_f)][abs(ry * oy_f)][1];
						pixels[j * rx * 4 + i * 4 + 2] = textures[11][abs(ry * ox_f)][abs(ry * oy_f)][2];
						
						//cout << angle << endl;
						pixels[(ry - j - 1) * rx * 4 + i * 4] = textures[0][abs(ry*ox_f)][abs(ry * oy_f)][0];
						pixels[(ry - j - 1) * rx * 4 + i * 4 + 1] = textures[0][abs(ry * ox_f)][abs(ry * oy_f)][1];
						pixels[(ry - j - 1) * rx * 4 + i * 4 + 2] = textures[0][abs(ry * ox_f)][abs(ry * oy_f)][2];
					}
					pixels[j * rx * 4 + i * 4 + 3] = 255;
					pixels[(ry - j - 1) * rx * 4 + i * 4 + 3] = 255;
				}
			}
		}
		else if (page == 0) {
			return bg_main;
		}
		return pixels;
	}

	void step_left(double step) {
		x_past = x;
		y_past = y;
		x -= sin(radians(angle)) * step;
		y += cos(radians(angle)) * step;
		integer_x = trunc(x);
		integer_y = trunc(y);
		if (map[integer_x][integer_y] > 0) {
			x = x_past;
			y = y_past;
		}
	}
	void step_right(double step) {
		x_past = x;
		y_past = y;
		x += sin(radians(angle)) * step;
		y -= cos(radians(angle)) * step;
		integer_x = trunc(x);
		integer_y = trunc(y);
		if (map[integer_x][integer_y] > 0) {
			x = x_past;
			y = y_past;
		}
	}
	void step_forward(double step) {
		x_past = x;
		y_past = y;
		x += cos(radians(angle)) * step;
		y += sin(radians(angle)) * step;
		integer_x = trunc(x);
		integer_y = trunc(y);
		if (map[integer_x][integer_y] > 0) {
			x = x_past;
			y = y_past;
		}
	}
	void step_backward(double step) {
		x_past = x;
		y_past = y;
		x -= cos(radians(angle)) * step;
		y -= sin(radians(angle)) * step;
		integer_x = trunc(x);
		integer_y = trunc(y);
		if (map[integer_x][integer_y] > 0) {
			x = x_past;
			y = y_past;
		}
	}

};

	//void renderTo(window_){
	//	int i, j, k;
	//	for (i = 0; i < rx; i++) {
	//		for (j = 0; j < ry; j++) {
	//			sf::Vertex point(sf::Vector2f(i, j), sf::Color{128,128,128});
	//			window_.draw(&point, 1, sf::Points);
	//		}
	//	}
	//}






#endif