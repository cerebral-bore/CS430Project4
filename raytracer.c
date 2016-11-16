/* 	Jesus Garcia
	Project 4 - Raytracing Illumination - 11/15/16
	CS430 - Prof. Palmer
	
	In the previous project you will wrote code to raycast and shade mathematical primitives based
	on a scene input file into a pixel buffer. In this project you will add recursive raytracing to
	provide reflection and refraction"
	
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>

// Struct and function declarations

typedef struct Object{
  int kind; // 1,2,3,4 = sphere, plane, camera, light respectively
  double diffuse_color[3];
  double specular_color[3];
  double reflectivity;
  double refractivity;
  double ior;
  union {
    struct {
      double center[3];
      double radius;
    } sphere;
    struct {
      double center[3];
      double width;
      double height;
      double normal[3];
    } plane;
    struct {
      double width;
      double height;
    } camera;
    struct light{
      double center[3];
      double color[3];
      double direction[3];
      double theta;
      double radial_a2;
      double radial_a1;
      double radial_a0;
      double angular_a0;
    } light;
  };
} Object;

typedef struct dPixel{
  double r,g,b;
} dPixel;

typedef struct Pixel{
  unsigned char r,g,b;
} Pixel;

typedef Object Scene[128];
int read_scene(char* json, Scene scene);
int raycast (Pixel *buffer, Scene scene, int num_objects, int width, int height);
int ppm_output(Pixel  *buffer, char *output_file_name, int size, int depth, int width, int height);

double* next_vector(FILE* json);
double next_number(FILE* json);
char* next_string(FILE* json);
void skip_ws(FILE* json);
void expect_c(FILE* json, int d);
int next_c(FILE* json);

double sphere_intersection(double* Ro, double* Rd, double* C, double r);
double plane_intersection(double* Ro, double* Rd, double* c, double *n);
static inline double sqr(double v) {
  return v*v;
}

void shoot_ray(Object *scene, int num_objects, double *Ro, double *Rd, double *t, int *o);
void shoot_ray_shadow(Object *scene, int num_objects, int best_obj, double *Ro, double *Rd, double *t, int *o);
void shade_rec(Object *sceneRef, int num_objects, double best_t, int best_obj, double *Ro, double *Rd, double *color, int level);

// Vector math methods
typedef double* V3;

static inline void v3_add(V3 a, V3 b, V3 c){
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
}

static inline void v3_sub(V3 a, V3 b, V3 c){
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

static inline void v3_scale(V3 a, double s, V3 c){
  c[0] = s * a[0];
  c[1] = s * a[1];
  c[2] = s * a[2];  
}

static inline double v3_dot(V3 a, V3 b){
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline void v3_cross(V3 a, V3 b, V3 c){
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

static inline void v3_normalize(double* v) {
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

static inline void v3_reflect(V3 A, V3 N, V3 C){
  double q = 2*v3_dot(A, N);
  double rscale[3];
  v3_scale(N, q, rscale);
  v3_sub(A, rscale, C);
  v3_normalize(C);
}

static inline double v3_mag(double* v){
  double d = (double)sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  return d;
}

static inline int clamp(double d){
  int final;
  if (d > 255.0) final = 255;
  else if (d < 0.0) final = 0;
  else final = floor(d);
  return final;
}

static inline double distance(V3 Ro, V3 p){
  return sqrt(sqr(Ro[0] - p[0]) + sqr(Ro[1] - p[1]) + sqr(Ro[2] - p[2]));
}

static inline double frad(double a0, double a1, double a2, double dl){
  return 1/(a2*sqr(dl) + a1*dl + a0);
}

static inline double fang(V3 Rd, V3 Ld, double theta, double a0){
  v3_normalize(Rd);
  v3_normalize(Ld);
  if(theta == 0)  
    return 1;
  double l_theta = v3_dot(Ld, Rd);
  if (acos(l_theta) > theta)
    return 0;
  else{
    return pow(l_theta, a0);
  }
}

int errCheck(int argc, char *argv[]);

int main(int argc, char* argv[]){
  
  errCheck(argc, argv);

  // Initialize scene
  Object *scene = malloc(sizeof(Object)*128);

  // Populate scene by parsing json
  int num_objects = read_scene(argv[3], scene);

  // Create the pixel buffer and populate it
  // Store the desired picture height and width
  int width = atoi(argv[1]);
  int height = atoi(argv[2]);
  Pixel pixbuffer[width*height];
  raycast(pixbuffer, scene, num_objects, width, height);

  // Write a PPM based on pixel buffer data
  ppm_output(pixbuffer, argv[4], sizeof(Pixel)*(width*height), 255, width, height);
  return 0;
}

int read_scene(char* json_input, Scene scene){
	int c;
	FILE* json = fopen(json_input, "r");
	if (json == NULL) {
		fprintf(stderr, "Error: Could not open file \"%s\"\n", json);
		exit(1);
	}

	skip_ws(json);

	// Find the beginning of the list
	expect_c(json, '[');
	skip_ws(json);

	// Find the objects

	int obj_num = 0;
	while (1) {
		int c = fgetc(json);
		if (c == ']') {
			fprintf(stderr, "Error: This is the worst scene file EVER.\n");
			fclose(json);
			return -1;
		}
		if (c == '{') {
			skip_ws(json);

			// Parse the object
			char* key = next_string(json);
			if (strcmp(key, "type") != 0) {
				fprintf(stderr, "Error: Expected \"type\" key\n");
				exit(1);
			}

			skip_ws(json);
			expect_c(json, ':');
			skip_ws(json);

			char* value = next_string(json);

			//extract the object type
			if (strcmp(value, "sphere") == 0) {
				scene[obj_num].kind = 1;
			} else if (strcmp(value, "plane") == 0) {
				scene[obj_num].kind = 2;
			} else if (strcmp(value, "camera") == 0) {
				scene[obj_num].kind = 3;
			} else if (strcmp(value, "light") == 0) {
				scene[obj_num].kind = 4;
			} else {
				fprintf(stderr, "Error: Unknown type, \"%s\"\n", value);
				exit(1);
			}

			skip_ws(json);

			//extract any other fields
			while (1) {
				// , }
				c = next_c(json);
				if (c == '}') {
					// stop parsing this object
					break;
				}
				else if (c == ',') {
					// read another field
					skip_ws(json);
					char* key = next_string(json);
					skip_ws(json);
					expect_c(json, ':');
					skip_ws(json);

					//check the object type, and insert into the appropriate fields
					if (scene[obj_num].kind == 3){
						if (strcmp(key, "width") == 0){
							double value = next_number(json);              
							memcpy(&scene[obj_num].camera.width, &value, sizeof(double));  
						} else if (strcmp(key, "height") == 0){
							double value = next_number(json);
							scene[obj_num].camera.height = value;
							memcpy(&scene[obj_num].camera.height, &value, sizeof(double));
						} else {
							fprintf(stderr, "Error: Unrecognized field \"%s\" for 'camera'.\n.", key);
							exit(1);
						}
					} else if (scene[obj_num].kind == 1){
						if (strcmp(key, "diffuse_color") == 0){
							double* value = next_vector(json);
							memcpy(&scene[obj_num].diffuse_color, value, sizeof(double)*3);
						} else if (strcmp(key, "specular_color") == 0){
							double* value = next_vector(json);
							memcpy(&scene[obj_num].specular_color, value, sizeof(double)*3);
						} else if (strcmp(key, "radius") == 0){
							scene[obj_num].sphere.radius = next_number(json);
						} else if (strcmp(key, "reflectivity") == 0){
							scene[obj_num].reflectivity = next_number(json);
						} else if (strcmp(key, "refractivity") == 0){
							scene[obj_num].refractivity = next_number(json);
						} else if (strcmp(key, "ior") == 0){
							scene[obj_num].ior = next_number(json);
						} else if (strcmp(key, "position") == 0){
							double* value = next_vector(json);
							memcpy(&scene[obj_num].sphere.center, value, sizeof(double)*3);
						} else{
							fprintf(stderr, "Error: Unrecognized field \"%s\" for 'sphere'.\n.", key);
							exit(1);
						}
					} else if (scene[obj_num].kind == 2){
						if (strcmp(key, "width") == 0){
							scene[obj_num].plane.width = next_number(json);
						} else if (strcmp(key, "height") == 0){
							scene[obj_num].plane.height = next_number(json);
						} else if (strcmp(key, "reflectivity") == 0){
							scene[obj_num].reflectivity = next_number(json);
						} else if (strcmp(key, "refractivity") == 0){
							scene[obj_num].refractivity = next_number(json);
						} else if (strcmp(key, "ior") == 0){
							scene[obj_num].ior = next_number(json);
						} else if (strcmp(key, "diffuse_color") == 0){
							double* value = next_vector(json);
							memcpy(&scene[obj_num].diffuse_color, value, sizeof(double)*3);
						} else if (strcmp(key, "specular_color") == 0){
							double* value = next_vector(json);
							memcpy(&scene[obj_num].specular_color, value, sizeof(double)*3);
						} else if (strcmp(key, "position") == 0){
							double* value = next_vector(json);
							memcpy(&scene[obj_num].plane.center, value, sizeof(double)*3);
						} else if (strcmp(key, "normal") == 0){
							double* value = next_vector(json);
							memcpy(&scene[obj_num].plane.normal, value, sizeof(double)*3);
						} else{
							fprintf(stderr, "Error: Unrecognized field \"%s\" for 'plane'.\n.", key);
							exit(1);
						}
					} else if (scene[obj_num].kind == 4){
						if (strcmp(key, "color") == 0){
							double* value = next_vector(json);
							memcpy(&scene[obj_num].light.color, value, sizeof(double)*3);
						} else if (strcmp(key, "position") == 0){
							double* value = next_vector(json);
							memcpy(&scene[obj_num].light.center, value, sizeof(double)*3);
						} else if (strcmp(key, "direction") == 0){
							double* value = next_vector(json);
							memcpy(&scene[obj_num].light.direction, value, sizeof(double)*3);
						} else if (strcmp(key, "theta") == 0){
							scene[obj_num].light.theta = next_number(json);
						} else if (strcmp(key, "radial-a2") == 0){
							scene[obj_num].light.radial_a2 = next_number(json);
						} else if (strcmp(key, "radial-a1") == 0){
							scene[obj_num].light.radial_a1 = next_number(json);
						} else if (strcmp(key, "radial-a0") == 0){
							scene[obj_num].light.radial_a0 = next_number(json);
						} else if (strcmp(key, "angular-a0") == 0){
							scene[obj_num].light.angular_a0 = next_number(json);    
						} else{
							fprintf(stderr, "Error: Unrecognized field \"%s\" for 'light'.\n.", key);
							exit(1);
						}
					}
				}
				skip_ws(json);
			}
		obj_num++;
		}
		skip_ws(json);
		c = next_c(json);
		
		// Check if end of object or end of json file
		if (c == ',') {
			skip_ws(json);
		} else if (c == ']') {
			fclose(json);
			return obj_num;
		} else {
			fprintf(stderr, "Error: Expecting ',' or ']'.\n");
			exit(1);
		}
	}

	// Return count of parsed objects
	return obj_num;
}

// Main raycast method that will check for collisions and populate pixels into a buffer
int raycast(Pixel *buffer, Scene scene, int num_objects, int width, int height){
	printf("raycast entered\n");
	// Initialize camera and its values
	double cx = 0;
	double cy = 0;
	double h = scene[0].camera.height;
	double w = scene[0].camera.width;

	// These will be the dimensions of the image
	int M = width;
	int N = height;
	
	double pixheight = h / M;
	double pixwidth = w / N;
	printf("\nM:%d N:%d h:%d w:%d pH:%lf pW:%lf\n", M, N, h, w, pixheight, pixwidth);

	int current_pixel = 0;
	int testpixel = 0;
	for (int y = 0; y < N; y += 1) {
		for (int x = 0; x < M; x += 1) {
			double Ro[3] = { 0, 0, 0 };
			double Rd[3] = {
				cx - (w / 2) + pixwidth * (x + 0.5),
				cy - (h / 2) + pixheight * (y + 0.5), 1 };
			double best_t;
			int best_obj;
			shoot_ray(scene, num_objects, Ro, Rd, &best_t, &best_obj);
			
			Pixel pix;
			if (best_t > 0 && best_t != INFINITY) {
				double color[3];

				color[0] = 0;
				color[1] = 0;
				color[2] = 0;

				shade_rec(scene, num_objects, best_t, best_obj, Ro, Rd, color, 1);

				pix.r = (unsigned char)clamp(color[0]*255);
				pix.g = (unsigned char)clamp(color[1]*255);
				pix.b = (unsigned char)clamp(color[2]*255);
			} else {
				pix.r = (unsigned char)clamp(0.1);
				pix.g = (unsigned char)clamp(0.1);
				pix.b = (unsigned char)clamp(0.1);
			}

			buffer[current_pixel++] = pix;
		}
	}
	return 0;
}

int ppm_output(Pixel *buffer, char *output_file_name, int size, int depth, int width, int height){
	printf("\nppm_output entered\n");
	FILE *output_file;
	output_file = fopen(output_file_name, "w");
	if (!output_file){
		fprintf(stderr, "\nError: Could not open file for write.");
		fclose(output_file);
		exit(1);
	} else {
		fprintf(output_file, "P3\n%d %d\n%d\n", width, height, depth);
		for (int i = width-1; i >=0; i--){
			for (int j = 0; j<height; j++){
				fprintf(output_file, "%d ", buffer[(i*width)+j].r);
				fprintf(output_file, "%d ", buffer[(i*width)+j].g);
				fprintf(output_file, "%d ", buffer[(i*width)+j].b);
			}
		fprintf(output_file, "\n");
		}
	}
  
	fclose(output_file);
	return 0;
}

// Error checking function to minimize code in main()
int errCheck(int args, char *argv[]){
	
	// Initial check to see if there are 3 input arguments on launch
	if ((args != 5) || (strlen(argv[3]) <= 5) || (strlen(argv[4]) <= 4)) {
		fprintf(stderr, "Error: Program requires usage: 'raytracer width height input.json output.ppm'");
		exit(1);
	}

	// Check the file extension of input and output files
	char *extIn;
	char *extOut;
	if(strrchr(argv[3],'.') != NULL){
		extIn = strrchr(argv[3],'.');
	}
	if(strrchr(argv[4],'.') != NULL){
		extOut = strrchr(argv[4],'.');
	}
	
	// Check to see if the inputfile is in .ppm format
	if(strcmp(extIn, ".json") != 0){
		printf("Error: Input file not a json");
		exit(1);
	}
	
	// Check to see if the outputfile is in .ppm format
	if(strcmp(extOut, ".ppm") != 0){
		printf("Error: Output file not a PPM");
		exit(1);
	}

	return(0);
}


double* next_vector(FILE* json){
	double* v = malloc(3 * sizeof(double));
	expect_c(json, '[');
	skip_ws(json);
	v[0] = next_number(json);
	skip_ws(json);
	expect_c(json, ',');
	skip_ws(json);
	v[1] = next_number(json);
	skip_ws(json);
	expect_c(json, ',');
	skip_ws(json);
	v[2] = next_number(json);
	skip_ws(json);
	expect_c(json, ']');
	return v;
}

double next_number(FILE* json){
	double value;
	fscanf(json, "%lf", &value);
	return value;
}

char* next_string(FILE* json){
	char buffer[129];
	int c = next_c(json);
	if (c != '"') {
		fprintf(stderr, "Error: Expected string");
		exit(1);
	}
	c = next_c(json);
	int i = 0;
	while (c != '"') {
		if (i >= 128) {
			fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
			exit(1);
		}
		if (c == '\\') {
			fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
			exit(1);
		}
		if (c < 32 || c > 126) {
			fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
			exit(1);
		}
		buffer[i] = c;
		i += 1;
		c = next_c(json);
	}
	buffer[i] = 0;
	return strdup(buffer);
}

void skip_ws(FILE* json){
	int c = next_c(json);
	while (isspace(c)) {
		c = next_c(json);
	}
	ungetc(c, json);
}

void expect_c(FILE* json, int d){
	int c = next_c(json);
	if (c == d){ return; }
	fprintf(stderr, "Error: Expected '%c'");
	exit(1);
}

int next_c(FILE* json){
	int c = fgetc(json);
	if (c == EOF) {
		fprintf(stderr, "Error: Unexpected end of file");
		exit(1);
	}
	return c;
}


double sphere_intersection(double* Ro, double* Rd, double* C, double r){
	double a = (sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]));
	double b = (2 * (Ro[0] * Rd[0] - Rd[0] * C[0] + Ro[1] * Rd[1] - Rd[1] * C[1] + Ro[2] * Rd[2] - Rd[2] * C[2]));
	double c = sqr(Ro[0]) - 2 * Ro[0] * C[0] + sqr(C[0]) + sqr(Ro[1]) - 2 * Ro[1] * C[1] + sqr(C[1]) + sqr(Ro[2]) - 2 * Ro[2] * C[2] + sqr(C[2]) - sqr(r);

	double det = sqr(b) - 4 * a * c;
	if (det <= 0){ return -1; }

	det = sqrt(det);
	
	double t0 = (-b - det) / (2 * a);
	if (t0 > 0){ return t0; }

	double t1 = (-b + det) / (2 * a);
	if (t1 > 0){ return t1; }

	return -1;
}

double plane_intersection(double* Ro, double* Rd, double* c, double *n){
	v3_normalize(n);
	v3_normalize(Rd);
	double d = -(v3_dot(c, n));
	double e = v3_dot(Ro, n);
	double f = v3_dot(Rd, n);

	double t = -(e + d)/f;
  if (t > 0){ return t; }
  return -1;
}

void shoot_ray(Object *scene, int num_objects, double *Ro, double *Rd, double *t, int *o){
	v3_normalize(Rd);

	double best_t = INFINITY;
	int best_obj = 0;
	for (int i = 1; i < num_objects; i += 1) {
		double t = 0;

		switch (scene[i].kind) {
		case 1:
			t = sphere_intersection(Ro, Rd, scene[i].sphere.center, scene[i].sphere.radius);
			break;
		case 2:
			t = plane_intersection(Ro, Rd, scene[i].plane.center, scene[i].plane.normal);
			break;
		case 4:
			break;
		default:
			fprintf(stderr, "Error: Unrecognized object in scene\n");
			exit(1);
		}

		if (t > 0 && t < best_t){
			best_t = t;
			best_obj = i;
		}
	}
	*t = best_t;
	*o = best_obj;
}

void shoot_ray_shadow(Object *scene, int num_objects, int best_obj, double *Ro, double *Rd, double *t, int *o){
	v3_normalize(Rd);

	double best_t_shadow = INFINITY;
	int best_obj_shadow = 0;
	for (int k=1; k < num_objects ; k += 1){
		if (k != best_obj){
			double t_shadow = 0;
			switch(scene[k].kind){
			case 1:
				t_shadow = sphere_intersection(Ro, Rd, scene[k].sphere.center, scene[k].sphere.radius);
				break;
			case 2:
				t_shadow = plane_intersection(Ro, Rd, scene[k].plane.center, scene[k].plane.normal);
				break;
			case 4:
				break;                
			default:
				fprintf(stderr, "Error: Unrecognized object in scene when looking for shadow\n");
				exit(1);
			}
			
			if (t_shadow > 0 && t_shadow < best_t_shadow){
				best_t_shadow = t_shadow;
				best_obj_shadow = k;
			}
		}
	}
	*t = best_t_shadow;
	*o = best_obj_shadow;
}

void shade_rec(Object *scene, int num_objects, double best_t, int best_obj, double *Ro, double *Rd, double *color, int level){

	for (int i = 1; i < num_objects; i += 1) {
		if(scene[i].kind != 4){ continue; }
		double Ron[3];
		double Rdn[3];
		double Ron_temp[3];

		v3_scale(Rd, best_t, Ron_temp);
		v3_add(Ron_temp, Ro, Ron);
		v3_sub(scene[i].light.center, Ron, Rdn);

		double B[3];
		v3_sub(scene[i].light.center, Ron, B);
		double dist = v3_mag(B);

		double best_t_shadow;
		int best_obj_shadow;
		shoot_ray_shadow(scene, num_objects, best_obj, Ron, Rdn, &best_t_shadow, &best_obj_shadow);

		if (best_t_shadow > 0 && best_t_shadow != INFINITY && best_t_shadow < dist ) { /* Do Nothing */ 
		} else {
			double N[3], L[3], R[3], V[3];

			if(scene[best_obj].kind == 2){
				memcpy(N, scene[best_obj].plane.normal, sizeof(double)*3); // if plane
			} else { v3_sub(Ron, scene[best_obj].sphere.center, N); }
			v3_normalize(N);

			v3_sub(scene[i].light.center, Ron, L);
			v3_normalize(L);

			v3_reflect(L, N, R);

			memcpy(V, Rd, sizeof(double)*3);
			v3_scale(V, -1.0, V);
			v3_normalize(V);

			double diff_nl = v3_dot(N, L);
			double diff[3];
			diff[0] = scene[best_obj].diffuse_color[0]*scene[i].light.color[0];
			diff[1] = scene[best_obj].diffuse_color[1]*scene[i].light.color[1];
			diff[2] = scene[best_obj].diffuse_color[2]*scene[i].light.color[2];
			v3_scale(diff, diff_nl, diff);  

			double spec_vr = -1*v3_dot(V, R);
			double spec[3];
			if(spec_vr > 0){
				spec[0] = (scene[best_obj].specular_color[0]*scene[i].light.color[0]*pow(spec_vr, 20));
				spec[1] = (scene[best_obj].specular_color[1]*scene[i].light.color[1]*pow(spec_vr, 20));
				spec[2] = (scene[best_obj].specular_color[2]*scene[i].light.color[2]*pow(spec_vr, 20));
			} else {
				spec[0] = 0;
				spec[1] = 0;
				spec[2] = 0;
			}

			double rad_a = frad(scene[i].light.radial_a0, scene[i].light.radial_a1, scene[i].light.radial_a2, dist);
			double ang_a;
			if(scene[i].light.theta != 0){
				ang_a = fang(Rd, scene[i].light.direction, scene[i].light.theta * M_PI / 180.0, scene[i].light.angular_a0);
			} else { ang_a = 1; }

			color[0] += rad_a * ang_a * (diff[0] + spec[0]);
			color[1] += rad_a * ang_a * (diff[1] + spec[1]);
			color[2] += rad_a * ang_a * (diff[2] + spec[2]);

			if (level < 2) {
			
				if(scene[best_obj].reflectivity > 0){
					double refl_ray[3];
					v3_reflect(Rd, N, refl_ray);
					double rec_t = INFINITY;

					int rec_obj;
					shoot_ray(scene, num_objects, Ron, refl_ray, &rec_t, &rec_obj);

					double refl_color[3] = {0, 0, 0};
			  
					if (rec_t != INFINITY){ shade_rec(scene, num_objects, rec_t, rec_obj, Ron, refl_ray, refl_color, level + 1); }

					v3_scale(refl_color, scene[best_obj].reflectivity, refl_color);

					v3_add(refl_color, color, color);
				}
				if(scene[best_obj].refractivity > 0){/* Do Nothing */}
			}
		}
	}
}
