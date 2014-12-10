/**
 * ray_tracer.cpp
 * CS230
 * -------------------------------
 * Implement ray tracer here.
 */

#define SET_RED(P, C)   (P = (((P) & 0x00ffffff) | ((C) << 24)))
#define SET_GREEN(P, C)  (P = (((P) & 0xff00ffff) | ((C) << 16)))
#define SET_BLUE(P, C) (P = (((P) & 0xffff00ff) | ((C) << 8)))

#include "ray_tracer.h"
#include <iostream>
#include <cmath>
#include <vector>

//WHOOOOOOOOOOOOOOOOOOOOOOOHOOOOOOOOOOOOOOOOOOOOO

using namespace std;

const double Object::small_t=1e-6;
//--------------------------------------------------------------------------------
// utility functions
//--------------------------------------------------------------------------------
double sqr(const double x)
{
    return x*x;
}

Pixel Pixel_Color(const Vector_3D<double>& color)
{
    Pixel pixel=0;
    SET_RED(pixel,(unsigned char)(min(color.x,1.0)*255));
    SET_GREEN(pixel,(unsigned char)(min(color.y,1.0)*255));
    SET_BLUE(pixel,(unsigned char)(min(color.z,1.0)*255));
    return pixel;
}
//--------------------------------------------------------------------------------
// Shader
//--------------------------------------------------------------------------------
Vector_3D<double> Phong_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
	typedef Vector_3D<double> TV;

	TV color; //add diffuse + ambient + specular and return color

	// calculate vector towards the viewer
	TV view_ray = (world.camera.position - intersection_point).Normalized();
	TV normal = same_side_normal.Normalized();

	// Regardless of shadows, add ambiant to color
	// Ambient lighting
	TV ambient = color_ambient*world.lights[0]->Emitted_Light(ray);
	color += ambient;
	
	// get vectors towards the two light sources
	TV l_ray1 = (world.lights[0]->position - intersection_point).Normalized();
	TV l_ray2 = (world.lights[1]->position - intersection_point).Normalized();

	//Create Ray objects to "Cast" shadow rays
	Ray light_ray1(intersection_point, l_ray1);
	Ray light_ray2(intersection_point, l_ray2);
	light_ray1.endpoint += normal*0.05;
	light_ray2.endpoint += normal*0.05;
	vector<Ray> light_vec; // store rays in vector for looping 
	light_vec.push_back(light_ray1);
	light_vec.push_back(light_ray2);

	if (world.enable_shadows == true) { 
		//count the number of obstructed light rays
		int num_shadows = 0;
		for (unsigned int i=0; i<light_vec.size(); i++) {
			// create dummy object to access background color
			TV dummy;
			const Object* obj = new Sphere(dummy, 1.0);
			obj = world.Closest_Intersection(light_vec[i]);
			delete obj;
			if (light_vec[i].semi_infinite == false)
				num_shadows++; 		
		}

		// if one lightray intersects an object, delete from further calcs
		if (num_shadows == 1) {
			if (light_vec[0].semi_infinite == false)
				light_vec.erase(light_vec.begin());
			else
				light_vec.erase(light_vec.begin()+1);
		}	

		// if both intersect, return ambiant light
		if (num_shadows == 2)
			return color;
	}

	// Diffuse
	TV diffuse;
	for (unsigned int i = 0; i < light_vec.size(); i++) {
		double l_n = TV::Dot_Product(light_vec[i].direction, normal);

		TV Ld = world.lights[i]->Emitted_Light(ray)*max((double)0, l_n);
		diffuse += color_diffuse*Ld;
	}

	//Spectral
	TV specular;
	for (unsigned int i = 0; i < light_vec.size(); i++) {
		TV reflected_ray = (normal*2*(TV::Dot_Product(light_vec[i].direction, normal))) - light_vec[i].direction; 
		reflected_ray.Normalize();
		double theta = max((double)0, TV::Dot_Product(reflected_ray, view_ray));
		theta = pow(theta, specular_power);
		TV Ls = world.lights[0]->Emitted_Light(ray)*theta;
		specular += color_specular*Ls; 
	}
	
	color += diffuse + specular;
    return color;
}

Vector_3D<double> Reflective_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
	typedef Vector_3D<double> TV;

	TV color; //add d+a+s and return color

	// calculate vector towards the viewer
	TV view_ray = (world.camera.position - intersection_point).Normalized();
	TV normal = same_side_normal.Normalized();

	// Regardless of shadows, add ambiant to color
	// Ambient lighting
	TV ambient = color_ambient*world.lights[0]->Emitted_Light(ray);
	color += ambient;
	
	// get vectors towards the two light sources
	TV l_ray1 = (world.lights[0]->position - intersection_point).Normalized();
	TV l_ray2 = (world.lights[1]->position - intersection_point).Normalized();

	//Create Ray objects to "Cast" shadow rays
	Ray light_ray1(intersection_point, l_ray1);
	Ray light_ray2(intersection_point, l_ray2);
	light_ray1.endpoint += normal*0.05;
	light_ray2.endpoint += normal*0.05;
	vector<Ray> light_vec;
	light_vec.push_back(light_ray1);
	light_vec.push_back(light_ray2);

	if (world.enable_shadows == true) {
		int num_shadows = 0;
		for (unsigned int i=0; i<light_vec.size(); i++) {	
			TV dummy;
			const Object* obj = new Sphere(dummy, 1.0);
			obj = world.Closest_Intersection(light_vec[i]);
			delete obj;
			if (light_vec[i].semi_infinite == false)
				num_shadows++; 		
		}

		// if one lightray intersects an object, delete from further calcs
		if (num_shadows == 1) {
			if (light_vec[0].semi_infinite == false)
				light_vec.erase(light_vec.begin());
			else
				light_vec.erase(light_vec.begin()+1);
		}	

		// if both intersect, return ambiant light
		if (num_shadows == 2)
			return color;
	}

	// Diffuse
	TV diffuse;
	for (unsigned int i = 0; i < light_vec.size(); i++) {
		double l_n = TV::Dot_Product(light_vec[i].direction, normal);

		TV Ld = world.lights[i]->Emitted_Light(ray)*max((double)0, l_n);
		diffuse += color_diffuse*Ld;
	}

	//Spectral
	TV specular;
	for (unsigned int i = 0; i < light_vec.size(); i++) {
		TV reflected_ray = (normal*2*(TV::Dot_Product(light_vec[i].direction, normal))) - light_vec[i].direction; 
		reflected_ray.Normalize();
		double theta = max((double)0, TV::Dot_Product(reflected_ray, view_ray));
		theta = pow(theta, specular_power);
		TV Ls = world.lights[0]->Emitted_Light(ray)*theta;
		specular += color_specular*Ls; 
	}
	
	// Create Ray object for reflection ray
	TV ray_reflection =  (normal*2*(TV::Dot_Product(view_ray, normal))) - view_ray; 
	ray_reflection.Normalize();
	Ray reflec_ray(intersection_point, ray_reflection);
	reflec_ray.endpoint += normal*0.5;
	reflec_ray.recursion_depth = ray.recursion_depth+1;
	reflec_ray.t_max = 1000000;

	if (ray.recursion_depth <= world.recursion_depth_limit)
		color += world.Cast_Ray(reflec_ray, ray)*reflectivity;

	color += diffuse + specular;
	return color;
}

Vector_3D<double> Flat_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    return color;
}

//--------------------------------------------------------------------------------
// Objects
//--------------------------------------------------------------------------------
// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Sphere::
Intersection(Ray& ray) const
{
	typedef Vector_3D<double> TV;
	
    TV d = ray.direction;
	TV temp = ray.endpoint - center;

	double a, b, c, disc;
	a = TV::Dot_Product(d, d);
	b = 2*(TV::Dot_Product(d, temp));
	c = TV::Dot_Product(temp, temp) - sqr(radius);

	disc = b*b - 4*a*c;

	if (disc < 0) // no intersection
		return false;


	//find closest endpoint and determine if it's the smallest t
	disc = sqrt(disc);
	double t1 = (-b + disc)/(2*a);
	double t2 = (-b - disc)/(2*a);

	// if only one intersection, positive "t1" intersects
	if (disc == 0) {
		if (t1 < ray.t_max) {
			ray.t_max = t1;
			ray.current_object = this; 
		}
		ray.semi_infinite = false;
		return true;
	}
		
	// if both intersect, determine shorter t 
	if (min(t1, t2) > 0) {
		if (min(t1, t2) < ray.t_max) {
			ray.t_max = min(t1, t2);
			ray.current_object = this;
		}
		ray.semi_infinite = false;	
		return true;
	}

	if (max(t1, t2) > 0) {
		if (max(t1, t2) < ray.t_max) {
			ray.t_max = max(t1, t2);
			ray.current_object = this; 
		}
		ray.semi_infinite = false;
		return true;
	}

	return false;
}

Vector_3D<double> Sphere::
Normal(const Vector_3D<double>& location) const
{
    Vector_3D<double> normal;

    normal.x = (location.x - center.x);
    normal.y = (location.y - center.y);
    normal.z = (location.z - center.z);

	normal.Normalize();
    return normal;
}


// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Plane::
Intersection(Ray& ray) const
{
	typedef Vector_3D<double> TV;
	TV d = ray.direction;
	TV P0 = ray.endpoint;
	TV q = x1;
	TV n = normal;

	d.Normalize();	
	n.Normalize();

	double denom = TV::Dot_Product(n, d);

	if (denom != 0)	{
		TV qP0 = q - P0;
		double t = TV::Dot_Product(qP0, n) / denom;
		
		if (t < 0) // t is behind the camera
			return false;

		if (t < ray.t_max) {
			ray.t_max = t;
			ray.current_object = this;
		}
		ray.semi_infinite = false;
		return true;
	}

	return false;
}

Vector_3D<double> Plane::
Normal(const Vector_3D<double>& location) const
{
    return normal;
}



//--------------------------------------------------------------------------------
// Camera
//--------------------------------------------------------------------------------
// Find the world position of the input pixel
Vector_3D<double> Camera::
World_Position(const Vector_2D<int>& pixel_index)
{
	// determine pixel index with -1 <= x,y <= 1 and map to image plane
	Vector_2D<double> indices = film.pixel_grid.X(pixel_index);
	Vector_3D<double> result = focal_point + horizontal_vector*indices.x +
						vertical_vector*indices.y;
	return result;	
}
//--------------------------------------------------------------------------------
// Render_World
//--------------------------------------------------------------------------------
// Find the closest object of intersection and return a pointer to it
//   if the ray intersects with an object, then ray.t_max, ray.current_object, and ray.semi_infinite will be set appropriately
//   if there is no intersection do not modify the ray and return 0
const Object* Render_World::
Closest_Intersection(Ray& ray)
{	
	for (unsigned int i = 0; i < objects.size(); i++)
		objects[i]->Intersection(ray);

	if (ray.semi_infinite == true) // no intersections
		return 0;
	else
		return ray.current_object;

}

// set up the initial view ray and call 
void Render_World::
Render_Pixel(const Vector_2D<int>& pixel_index)
{
	Vector_3D<double> world_pos = camera.World_Position(pixel_index);
	Vector_3D<double> world_pos_direc = (world_pos-camera.position).Normalized();
    Ray ray(camera.position, world_pos_direc); 
	ray.t_max = 1000000;
	ray.recursion_depth = 0;

    Ray dummy_root;
    Vector_3D<double> color=Cast_Ray(ray,dummy_root);
    camera.film.Set_Pixel(pixel_index,Pixel_Color(color));
}

// cast ray and return the color of the closest intersected surface point, 
// or the background color if there is no object intersection
Vector_3D<double> Render_World::
Cast_Ray(Ray& ray,const Ray& parent_ray)
{
 	const Object* obj = Closest_Intersection(ray);
	Vector_3D<double> intersection_point = ray.Point(ray.t_max);

	//if no object intersection, create dummy variables to satisfy Shader function reqs
	//and set color to background
	if (obj == 0) {
		Vector_3D<double> dummy_vec(0,0,0);
		Sphere* obj = new Sphere(dummy_vec, 1.0);
		Vector_3D<double> color = background_shader->Shade_Surface(ray, *obj, dummy_vec, dummy_vec); 
		delete obj;
		return color; 
	}
	
	// call corresponding shader with ray/object information
    Vector_3D<double> color = obj->material_shader->Shade_Surface(ray, *obj, intersection_point, obj->Normal(intersection_point));  

    return color;
}
