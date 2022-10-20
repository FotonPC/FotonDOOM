#ifndef DEF
#define DEF

#include <iostream>
#include <math.h>
#include <SFML\\Graphics.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
//#include <omp.h>
#include <time.h> 
#include <chrono>
//#include <amp.h>
#include <thread>
#include "sprites.h"
#include <algorithm>

using namespace std;
//using namespace concurrency;

int index_global;

bool comp_fs(vector<double> a, vector<double> b) {
	return a[1] < b[1];
}

struct OZ_Code_struct {
	double oz;
	int part;
	int a;
	int b;
	int c;
};

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

class Point
{
public:
	float x, y;
}; 
class LineSegment
{
public:
	Point top;
	Point bottom;
};

Point* getIntersectionPoint(LineSegment line1, LineSegment line2)
{
	float s1_x, s1_y, s2_x, s2_y;
	float p1_x = line1.bottom.x, p1_y = line1.bottom.y;
	float p0_x = line1.top.x, p0_y = line1.top.y;
	float p3_x = line2.bottom.x, p3_y = line2.bottom.y;
	float p2_x = line2.top.x, p2_y = line2.top.y;
	s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
	s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

	float s, t;
	s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
	t = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

	if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
	{
		Point* intsPoint = (Point*)malloc(sizeof(Point));
		intsPoint->x = p0_x + (t * s1_x);
		intsPoint->y = p0_y + (t * s1_y);
		return intsPoint;
	}

	return NULL;
}

class RectSprite3D {
public:
	double x, y;
	double x_size;
	double y_size;
	bool orient_x;
	int res_x_texture, res_y_texture, ry, tc1, tc2, tc3;
	vector<vector<vector<vector<char>>>> textures_x;
	vector<vector<vector<vector<char>>>> textures_y;
	long long id;
	void init(double x_, double y_, double x_size_, double y_size_, int res_x_texture_, int res_y_texture_, int ry_, long long id_, int text_code_1, int text_code_2, int text_code_3) {
		id = id_;
		x = x_;
		y = y_;
		tc1 = text_code_1;
		tc2 = text_code_2;
		tc3 = text_code_3;
		x_size = x_size_;
		y_size = y_size_;
		res_x_texture = res_x_texture_;
		res_y_texture = res_y_texture_;
		ry = ry_;
		orient_x = true;
		//_engine = ENGINE;
	}
	void load_texture_x(string filename, int texture_code) {
		sf::Image texture_img;
		texture_img.loadFromFile(filename);
		const sf::Uint8* pixels1 = new sf::Uint8[res_x_texture * ry * 4 + 4];
		sf::Uint8* pixels2 = new sf::Uint8[res_x_texture * ry * 4 + 4];

pixels1 = texture_img.getPixelsPtr();
int i, j, k;
for (i = 0; i < res_x_texture; i++) {
	for (j = 0; j < ry; j++) {
		textures_x[texture_code][i][j][0] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 4]);
		textures_x[texture_code][i][j][1] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 3]);
		textures_x[texture_code][i][j][2] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 2]);
	}
}
	}
	void load_texture_y(string filename, int texture_code) {
		sf::Image texture_img;
		texture_img.loadFromFile(filename);
		const sf::Uint8* pixels1 = new sf::Uint8[res_y_texture * ry * 4 + 4];
		sf::Uint8* pixels2 = new sf::Uint8[res_y_texture * ry * 4 + 4];

		pixels1 = texture_img.getPixelsPtr();
		int i, j, k;
		for (i = 0; i < res_y_texture; i++) {
			for (j = 0; j < ry; j++) {
				textures_y[texture_code][i][j][0] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 4]);
				textures_y[texture_code][i][j][1] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 3]);
				textures_y[texture_code][i][j][2] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 2]);
			}
		}
	}
	void set_orient_x(bool is23) {

	}
	bool is_collision(double x_, double y_) {
		if (abs(x - x_) <= x_size / 2 && abs(y - y_) <= y_size / 2) {
			return true;
		}
		else return false;
	}
	OZ_Code_struct get_oz(double x_, double y_, double angle) {
		double ox, oy;
		OZ_Code_struct res;
		//cout << orient_x;

		if (orient_x) {
			ox = (x_ - (x - x_size / 2)) / x_size;
			oy = (y_ - (y - y_size / 2)) / y_size;
			if (abs(ox - 0.5) < abs(oy - 0.5)) {
				res = OZ_Code_struct();
				res.oz = ox;
				if (id == 1) {
					if (angle > 180) {
						res.part = tc2;
					}
					else {
						res.part = tc1;
					}
				}
				else {
					res.part = tc1;
				}
				return res;
			}
			else {
				res = OZ_Code_struct();
				res.oz = 1 - oy;
				if (id == 0) {
					res.part = tc2;
				}
				else if (id == 1) {
					res.part = tc3;
				}
				return res;
			}
		}
		else {
			ox = (x_ - (x - x_size / 2)) / x_size;
			oy = (y_ - (y - y_size / 2)) / y_size;
			if (abs(ox - 0.5) < abs(oy - 0.5)) {
				res = OZ_Code_struct();
				res.oz = ox;
				if (id == 0) {
					res.part = tc2;
				}
				else if (id == 1) {
					res.part = tc3;
				}
				return res;
			}
			else {
				res = OZ_Code_struct();
				res.oz = oy;
				res.part = tc1;
				return res;
			}
		}
	}
};

class FlatSprite3D {
public:
	double x, y, size;
	bool orient_x;
	long long id;
	int ry, texture_size_x;
	vector<vector<vector<char>>> image;
	vector<vector<vector<double>>> triangles;
	void (*get_capture)(double, double, long long, bool, vector<vector<vector<char>>>&, vector<vector<vector<double>>>&) {

	};


	void init(double x_, double y_, bool orient_x_, double size_, long long id_, int ry_, int texture_size_x_) {
		//void (*get_capture_)(double, double, long long, bool, vector<vector<vector<char>>>&, vector<vector<vector<double>>>&), int ry_, int texture_size_x_) {
		//get_capture = get_capture_;
		x = x_;
		y = y_;
		orient_x = orient_x_;
		size = size_;
		id = id_;
		ry = ry_;
		texture_size_x = texture_size_x_;
		image.resize(texture_size_x, vector<vector<char>>(ry, vector<char>(4, 0)));
	}
	void render_stage(double angle, double x_, double y_) {

	}

	void load_texture(int texture_code, string filename) {
		sf::Image texture_img;
		texture_img.loadFromFile(filename);
		const sf::Uint8* pixels1 = new sf::Uint8[texture_size_x * ry * 4 + 4];
		sf::Uint8* pixels2 = new sf::Uint8[texture_size_x * ry * 4 + 4];

		pixels1 = texture_img.getPixelsPtr();
		int i, j, k;
		for (i = 0; i < ry; i++) {
			for (j = 0; j < ry; j++) {
				image[i][j][0] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 4]);
				image[i][j][1] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 3]);
				image[i][j][2] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 2]);
				image[i][j][3] = sfUint8_to_int(pixels1[i * ry * 4 + j * 4 - 1]);
			}
		}


	}
	bool is_collision(double x1, double y1, double x2, double y2) {
		double x3, x4, y3, y4;
		if (orient_x) {
			x3 = x - size / 2;
			x4 = x + size / 2;
			y3 = y - 0.01;
			y4 = y + 0.01;
		}
		else {
			y3 = y - size / 2;
			y4 = y + size / 2;
			x3 = x - 0.01;
			x4 = x + 0.01;
		}
		double dx0 = x2 - x1;
		double dx1 = x4 - x3;
		double dy0 = y2 - y1;
		double dy1 = y4 - y3;
		double p0 = dy1 * (x4 - x1) - dx1 * (y4 - y1);
		double p1 = dy1 * (x4 - x2) - dx1 * (y4 - y2);
		double p2 = dy0 * (x2 - x3) - dx0 * (y2 - y3);
		double p3 = dy0 * (x2 - x4) - dx0 * (y2 - y4);
		return (p0 * p1 <= 0) && (p2 * p3 <= 0);
	}
	double get_oz(double c_x, double c_y) {
		if (orient_x) {
			return (c_x - (x - size / 2)) / size;
		}
		else {
			return abs(c_y - (y - size / 2)) / size;
		}
	}
};

class Engine3D {

public:
	int rx, ry, RESOLUTION_X, RESOLUTION_Y;
	double x, y, angle;
	vector<vector<vector<char>>> render_matrix;
	vector<vector<vector<float>>> rays_steps;
	vector<vector<int>> map;
	vector<vector<vector<double>>> sprites_flat_heights;
	vector<vector<vector<double>>> sorted_i_sf;
	vector<vector<vector<vector<char>>>> textures;
	vector<vector<double>> heights;
	vector<vector<vector<vector<float>>>> floor_ceil_pre_counted_values;
	vector<RectSprite3D> sprites_objects;
	vector<FlatSprite3D> flat_sprites;
	float* floor_ceil_pcv_arr;
	thread* threads_;
	int d_angle;
	int threads_num;
	int integer_x, integer_y, page, map_width, map_height, n_of_objects, n_flat_sprites;
	
	
	sf::Uint8* bg_main;
	sf::Texture bg_texture;
	sf::Image bg_image;
	double x_past, y_past;
	double ray_step;
	
	sf::Uint8* pixels;


	void init(double x_, double y_, double angle_, string bg_main_filename, double ray_step_k, int rx_, int ry_, int n_of_objects_, int n_flat_sprites_) {
		RESOLUTION_X = rx_;
		RESOLUTION_Y = ry_;
		index_global = 0;
		rx = rx_;
		ry = ry_;
		n_of_objects = n_of_objects_;
		threads_num = 8;
		threads_ = new thread[threads_num];
		floor_ceil_pcv_arr = new float[180 * rx * ry * 2];
		sprites_objects.resize(n_of_objects);
		n_flat_sprites = n_flat_sprites_;
		//floor_ceil_pre_counted_values.resize(360, vector<vector<vector<float>>>( rx, vector<vector<float>>(ry/2, vector<float>(2, 0))));
		//render_matrix = render_matrix(360, vector<vector<int>>(rx, vector<int>(2)));
		//rays_steps = rays_steps(360, vector<vector<double>>(rx, vector<double>(2)));
		const sf::Uint8* bg_const = new sf::Uint8[RESOLUTION_X * RESOLUTION_Y * 4 + 4];
		pixels = new sf::Uint8[RESOLUTION_X * RESOLUTION_Y * 4 + 4];
		bg_main = new sf::Uint8[RESOLUTION_X * RESOLUTION_Y * 4 + 4];
		page = 1;
		flat_sprites.resize(n_flat_sprites);
		heights.resize(rx, vector<double>(10, 0));
		sprites_flat_heights.resize(n_flat_sprites, vector<vector<double>>(rx, vector<double>(10, 0)));
		int gg, g;
		sorted_i_sf.resize(rx, vector<vector<double>>(n_flat_sprites, vector<double>(2,0)));
		for (g = 0; g < rx; g++) {
			for (gg = 0; gg < n_flat_sprites; gg++) {
				sorted_i_sf[g][gg][0] = 0;
				sorted_i_sf[g][gg][1] = 0;
			}
		}
		render_matrix.resize(rx, vector<vector<char>>(ry, vector<char>(3, 0)));
		rays_steps.resize(360, vector<vector<float>>(rx, vector<float>(2, 0)));
		textures.resize(25, vector<vector<vector<char>>>(ry, vector<vector<char>>(ry, vector<char>(3, 0))));
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
		ray_step = 1.0 / to_double(rx) *  ray_step_k;
		std::cout << "inizialisation rays_steps" << endl;
		//rays_steps[359][127][1] = 2.4;
		//std::cout << to_double(abs(rx / 2 - 0)) / to_double(rx) * 10 << endl;
		for (i = 0; i < 360; i++) {
			for (j = 0; j < rx; j++) {
				k = 0;
				rays_steps[i][j][0] = cos(radians(i - (-rx / 2.0 + j) / rx * d_angle)) * ray_step;
				rays_steps[i][j][1] = sin(radians(i - (-rx / 2.0 + j) / rx * d_angle)) * ray_step;
			}
		}
		for (k = 0; k < 360; k++) {
			for (i = 0; i < rx; i++) {
				for (j = 0; j < ry / 2; j++) {
					floor_ceil_pcv_arr[k * rx * ry + i * ry + j * 2] = 1.0 / (ry / 2 - j) * ry / 2 / cos(radians(to_double(i) / rx * 90.0 - 45.0)) * rays_steps[k][i][0] / ray_step;;
					floor_ceil_pcv_arr[k * rx * ry + i * ry + j * 2 + 1] = 1.0 / (ry / 2 - j) * ry / 2 / cos(radians(to_double(i) / rx * 90.0 - 45.0)) * rays_steps[k][i][1] / ray_step;
					//floor_ceil_pre_counted_values[k][i][j][0] = 1.0 / (ry / 2 - j) * ry / 2 / cos(radians(to_double(i) / rx * 90.0 - 45.0)) * rays_steps[k][i][0] / ray_step;
					//floor_ceil_pre_counted_values[k][i][j][1] = 1.0 / (ry / 2 - j) * ry / 2 / cos(radians(to_double(i) / rx * 90.0 - 45.0)) * rays_steps[k][i][1] / ray_step;
				}
			}
		}
		//cout << rays_steps[0][0][0] << endl;

	}
	/*void sort_fs_by_H(int ind) {
		int i, j, k, ij, id_max;
		bool in_s;
		double h_max;
		for (i = 0; i < n_flat_sprites; i++) {
			sorted_i_sf[i] = -1;
		}
		for (i = 0; i < n_flat_sprites; i++) {
			id_max = 0;
			h_max = 0;
			for (j = 0; j < n_flat_sprites; j++) {
				in_s = false;
				for (k)
				if (not j in sorted_i_sf) {
					if (sprites_flat_heights[j][ind][2] > h_max) {
						h_max = sprites_flat_heights[j][ind][2];
						id_max = j;
					}
				}
			}
			sorted_i_sf[n_flat_sprites - ij] = id_max;
			ij += 1;
		}
	}*/
	void add_rect_sprite(RectSprite3D sprite, int index_) {
		sprites_objects[index_] = sprite;
	}
	void add_flat_sprite(FlatSprite3D sprite, int index_) {
		flat_sprites[index_] = sprite;
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
	void thread_stage_heights(int i_start, int i_end, int num_threads) {
		
		//cout << "stage-heights" << endl; // debug info
		int i, j, k, ijk; // ���������
		int int_angle, int_x, int_y; // �����
		int code;
		double x_step, y_step;
		double current_x, current_y, dx, dy, ox, oy, oz; // ������� x y
		bool sprite_is;
		int_angle = round(angle);
		oz = 0;

		for (i = 0; i < rx; i++) {
			//int_x = trunc(x); // �������� � ������ ����� ��� ������ �������
			//int_y = trunc(y); // �������� � ������ ����� ��� ������ �������
			k = 0;
			//cout << "x:" << current_x << " y:" << current_y << " a:" << int_angle << " i:" << i << endl;
			current_x = x;
			current_y = y;
			x_step = rays_steps[int_angle][i][0];
			y_step = rays_steps[int_angle][i][1];
			for (ijk = 0; ijk < n_flat_sprites; ijk++) {
				sprites_flat_heights[ijk][i][1] = 0;
			}
			//cout << x << " " << y << ":" << map[0][2] << " " << int_x <<  " " << int_y << endl;
			sprite_is = false;
			while ((map[current_x][current_y] == 0) and (k < 10 / ray_step)) { // �������� �� ��������� ��
				current_x += x_step; // ������ ��� x
				current_y += y_step; // ������ ��� y
				//int_x = trunc(current_x); // �������� � ������ ����� ��� ������ �������
				//int_y = trunc(current_y); // �������� � ������ ����� ��� ������ �������
				for (ijk = 0; ijk < n_of_objects; ijk++) {
					
					if (sprites_objects[ijk].is_collision(current_x, current_y)) {
						//cout << "gg" << endl;
						code = sprites_objects[ijk].get_oz(current_x, current_y, angle).part;
						oz = sprites_objects[ijk].get_oz(current_x, current_y, angle).oz;
						k = to_int(10 / ray_step);
						sprite_is = true;
					}
				}
				for (ijk = 0; ijk < n_flat_sprites; ijk++) {
					if (flat_sprites[ijk].is_collision(current_x, current_y, current_x + x_step, current_y + y_step)) {
						sprites_flat_heights[ijk][i][0] = flat_sprites[ijk].get_oz(current_x, current_y);
						sprites_flat_heights[ijk][i][1] = 1;
						dx = current_x - x;
						dy = current_y - y;
						
						sprites_flat_heights[ijk][i][2] = abs((1 / (sqrt(dx * dx + dy * dy) * cos(radians(to_double(i) / to_double(rx) * 90.0 - 45.0)))) * ry / 2.0);
					}
				}
				k++;
			}

			
			//cout << cos(radians(abs(int_angle - d_angle / to_double(rx) * i))) << " " << i << " " << d_angle / to_double(rx) << endl;
			dx = current_x - x;
			dy = current_y - y;
			heights[i][3] = dx; // x ����������
			heights[i][4] = dy; // y ����������
			heights[i][0] = abs((1 / (sqrt(dx * dx + dy * dy) * cos(radians(to_double(i) / to_double(rx) * 90.0 - 45.0)))) * ry / 2.0); // ������
			//cout << heights[i][0] << endl;
			if (not sprite_is) {
				code = map[current_x][current_y]; // ��� ��������
				if (k >= 10 / ray_step) {
					code = 10;
				}
				ox = current_x - trunc(current_x); // �������� �������� ��� ���������� y - ���
				oy = current_y - trunc(current_y); // �������� �������� ��� ���������� x - ���
				if (abs(ox - 0.5) < abs(oy - 0.5)) {
					oz = ox;
				}
				else if (abs(ox - 0.5) > abs(oy - 0.5)) {
					oz = oy;
				}
			}
			heights[i][1] = oz; // ������ ��������
			heights[i][2] = code; // ��� ��������
		}
	}
	void stage_heights() {
		int i;
		threads_ = new thread[threads_num];
		if (page == 1) {
			while (angle >= 360) {
				angle -= 360.0;
			}
			while (angle < 0) {
				angle += 360.0;
			}
			for (i = 0; i < threads_num; i++) {
				threads_[i] = thread(&Engine3D::thread_stage_heights, this, rx / threads_num * i, rx / threads_num * (i + 1), threads_num);
			}
			for (i = 0; i < threads_num; i++) {
				threads_[i].join();
			}
			
		}
	}
	void thread_stage_render(int i_start, int i_end, int threads_number) {
		int i, j, k, code, pix0, pix1, pix2, pix02, pix12, pix22, arxf, aryf, up_code,down_code, ijk, ij;
		double oz, i1, i2, i3, i4, distance_for_floor, dx_f, dy_f, ox_f, oy_f, light1, max_h, i2_s, i3_s;
		int pix0_s, pix1_s, pix2_s, pix00_s, pix11_s, pix22_s;
		bool sprite_is;
		sprite_is = true;
		if (page == 1) {
#pragma omp parallel for 
			for (i = i_start; i < i_start + rx / threads_number; i++) {
				//array_view <const int, 2> avA(ry, 10, heights);
				//array_view <const int, 3> avB(ry, ry, 3, textures);
				//array_view <int, 2> avC(m, n, C);
				//avC.discard_data();
				max_h = 0;
				for (ij = 0; ij < n_flat_sprites; ij++) {
					sorted_i_sf[i][ij][0] = ij;
					sorted_i_sf[i][ij][1] = sprites_flat_heights[ij][i][2];
				}
				//std::sort(sorted_i_sf.begin(), sorted_i_sf.end(), comp_fs);
				std::sort(sorted_i_sf[i].begin(), sorted_i_sf[i].end(), [](auto& left, auto& right) {
					return left[1] > right[1];
					});
				for (j = 0; j < ry / 2; j++) {
					if ((ry / 2 - j) < heights[i][0]) {
						light1 = sqrt(sqrt(sqrt(sqrt(1 / heights[i][0] * 10000))));
						oz = heights[i][1];
						code = to_int(heights[i][2]);
						i1 = trunc(oz * ry);
						i2 = ry / 2 - ry / 2 / heights[i][0] * (ry / 2 - j);
						i3 = ry / 2 - ry / 2 / heights[i][0] * (j - ry / 2);
						//cout << heights[i][0] / to_double(ry) * 2.0 * (ry / 2 - j - 1) << endl;
						pix0_s = 0;
						pix1_s = 0;
						pix2_s = 0;
						pix00_s = 0;
						pix11_s = 0;
						pix22_s = 0;
						index_global = i;
						//bool (*ptr) (vector<vector<double>>, vector<vector<double>>) = nullptr;
						//ptr = &comp_fs;
						
						for (ijk = 0; ijk < n_flat_sprites; ijk++) {
							//cout << sprites_flat_heights[ijk][i][1] << endl;
							//cout << to_int(sorted_i_sf[ijk][0]) << "g" << i << "h" << endl;
							if (sprites_flat_heights[to_int(sorted_i_sf[i][ijk][0])][i][1] == 1) {
								//cout << "g";
								if (((pix0_s == 0) && ((pix1_s == 0) && (pix2_s == 0)))) {
									//cout << sorted_i_sf[ijk][0] << endl;
									i2_s = ry / 2 - ry / 2 / sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][2] * (ry / 2 - j);
									i3_s = ry / 2 - ry / 2 / sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][2] * (j - ry / 2);
									pix0_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i2_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][0];
									pix1_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i2_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][1];
									pix2_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i2_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][2];
								}
								if (pix00_s == 0 && pix11_s == 0 and pix22_s == 0) {
									pix00_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i3_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][0];
									pix11_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i3_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][1];
									pix22_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i3_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][2];
								}
							}
						}

						if (pix0_s == 0 && pix1_s == 0 && pix2_s == 0) {
							pix0 = textures[code][i2][i1][0];
							pix1 = textures[code][i2][i1][1];
							pix2 = textures[code][i2][i1][2];
							//cout << ry / 2 / heights[i][0] * (ry / 2 + j - 1) << endl;

						}
						else {
							pix0 = pix0_s;
							pix1 = pix1_s;
							pix2 = pix2_s;
						}
						if (pix00_s == 0 && pix11_s == 0 && pix22_s == 0) {
							pix02 = textures[code][i3][i1][0];
							pix12 = textures[code][i3][i1][1];
							pix22 = textures[code][i3][i1][2];
							//cout << ry / 2 / heights[i][0] * (ry / 2 + j - 1) << endl;

						}
						else {
							pix02 = pix00_s;
							pix12 = pix11_s;
							pix22 = pix22_s;
						}
						pixels[j * rx * 4 + i * 4] = pix0;
						pixels[j * rx * 4 + i * 4 + 1] = pix1;
						pixels[j * rx * 4 + i * 4 + 2] = pix2;
						pixels[(ry - j - 1) * rx * 4 + i * 4] = pix02;
						pixels[(ry - j - 1) * rx * 4 + i * 4 + 1] = pix12;
						pixels[(ry - j - 1) * rx * 4 + i * 4 + 2] = pix22;
						//cout << pix0 << " " <<  pix1 << endl;
					}
					else {
						sprite_is = false;
						for (ijk = 0; ijk < n_flat_sprites; ijk++) {
							//cout << ijk << " " << i << endl;
							if ((ry / 2 - j) < sprites_flat_heights[ijk][i][2] && sprites_flat_heights[ijk][i][2] > heights[i][0])
								sprite_is = true;
						}
						pix0_s = 0;
						pix1_s = 0;
						pix2_s = 0;
						pix00_s = 0;
						pix11_s = 0;
						pix22_s = 0;
						if (sprite_is) {
							
							//cout << heights[i][0] / to_double(ry) * 2.0 * (ry / 2 - j - 1) << endl;

							for (ijk = 0; ijk < n_flat_sprites; ijk++) {
								//cout << sprites_flat_heights[ijk][i][1] << endl;
								if (sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][1] == 1 && (ry / 2 - j) < sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][2]) {
									//cout << "g";
									if (pix0_s == 0 && pix1_s == 0 and pix2_s == 0) {
										
										i2_s = ry / 2 - ry / 2 / sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][2] * (ry / 2 - j);
										//cout << min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x) << endl;
										pix0_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i2_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][0];
										pix1_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i2_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][1];
										pix2_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i2_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][2];
									}
									if (pix00_s == 0 && pix11_s == 0 and pix22_s == 0) {
										i3_s = ry / 2 - ry / 2 / sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][2] * (j - ry / 2);
										pix00_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i3_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][0];
										pix11_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i3_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][1];
										pix22_s = flat_sprites[sorted_i_sf[i][ijk][0]].image[i3_s][min((double)ry - 1, sprites_flat_heights[sorted_i_sf[i][ijk][0]][i][0] * flat_sprites[sorted_i_sf[i][ijk][0]].texture_size_x)][2];
									}
									
								}
							}


						}




						//distance_for_floor = 1.0 / (ry / 2 - j) * ry / 2 / cos(radians(to_double(i) / rx * 90.0 - 45.0));
						//dx_f = floor_ceil_pre_counted_values[(int)angle][i][j][0] + x;
						//dy_f = floor_ceil_pre_counted_values[(int)angle][i][j][1] + y;
						dx_f = floor_ceil_pcv_arr[(int)angle * rx * ry + i * ry + j * 2] + x;
						dy_f = floor_ceil_pcv_arr[(int)angle * rx * ry + i * ry + j * 2 + 1] + y;
						//dx_f = distance_for_floor * rays_steps[angle][i][0]/ray_step + x;
						//dy_f = distance_for_floor * rays_steps[angle][i][1] / ray_step + y;

						ox_f = dx_f - trunc(dx_f);
						oy_f = dy_f - trunc(dy_f);
						arxf = abs(ry * ox_f);
						aryf = abs(ry * oy_f);
						//cout << arxf << endl;
						if (dy_f > 9) {
							up_code = 4;
						}
						else {
							up_code = 11;
						}
						if (dy_f > 9) {
							down_code = 15;
						}
						else {
							down_code = 0;
						}

						if (pix0_s == 0 && pix1_s == 0 && pix2_s == 0) {
							pixels[j * rx * 4 + i * 4] = textures[up_code][arxf][aryf][0] * 0.95;
							pixels[j * rx * 4 + i * 4 + 1] = textures[up_code][arxf][aryf][1] * 0.95;
							pixels[j * rx * 4 + i * 4 + 2] = textures[up_code][arxf][aryf][2] * 0.95;
							//cout << ry / 2 / heights[i][0] * (ry / 2 + j - 1) << endl;

						}
						else {
							pixels[j * rx * 4 + i * 4] = pix0_s;
							pixels[j * rx * 4 + i * 4 + 1] = pix1_s;
							pixels[j * rx * 4 + i * 4 + 2] = pix2_s;

						}
						if (pix00_s == 0 && pix11_s == 0 && pix22_s == 0) {
							//cout << angle << endl;
							pixels[(ry - j - 1) * rx * 4 + i * 4] = textures[down_code][arxf][aryf][0];
							pixels[(ry - j - 1) * rx * 4 + i * 4 + 1] = textures[down_code][arxf][aryf][1];
							pixels[(ry - j - 1) * rx * 4 + i * 4 + 2] = textures[down_code][arxf][aryf][2];
						}
						else {
							pixels[(ry - j - 1) * rx * 4 + i * 4] = pix00_s;
							pixels[(ry - j - 1) * rx * 4 + i * 4 + 1] = pix11_s;
							pixels[(ry - j - 1) * rx * 4 + i * 4 + 2] = pix22_s;
						}
						
					}
						pixels[j * rx * 4 + i * 4 + 3] = 255;
						pixels[(ry - j - 1) * rx * 4 + i * 4 + 3] = 255;
					}
				}
			
		}
	}
	
	sf::Uint8* stage_render(sf::Texture texture_target, sf::Sprite sprite_target) {
		//cout << "stage_render" << endl;
		int i;
		threads_ = new thread[threads_num];
		if (page ==  1) {
			
			for (i = 0; i < threads_num; i++) {
				threads_[i] = thread(&Engine3D::thread_stage_render, this, rx/ threads_num *i, rx/ threads_num *(i+1), threads_num);
			}
			for (i = 0; i < threads_num; i++) {
				threads_[i].join();
			}
			return pixels;
			
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
		bool coll;
		int i;
		x_past = x;
		y_past = y;
		x += cos(radians(angle)) * step;
		y += sin(radians(angle)) * step;
		integer_x = trunc(x);
		integer_y = trunc(y);
		coll = map[integer_x][integer_y] > 0;
		for (i = 0; i < n_of_objects; i++) {
			if (sprites_objects[i].is_collision(x, y)) {
				coll = true;
				break;
			}
		}
		if (coll) {
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