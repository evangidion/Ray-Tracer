#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <cfloat>

// #define eps 0.0000001

using namespace std;
using namespace parser;

struct Blaster {
    float t;
    bool shot;
    int material_id;
    Vec3f intersectionP;
    Vec3f normal;
} ;

struct Ray {
    Vec3f start;
    Vec3f direction;
} ;

float vectorLength(const Vec3f &v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vec3f vectorNormalize(const Vec3f &v) {
    float l = vectorLength(v);
    return {v.x / l, v.y / l, v.z / l};
}

Vec3f vectorCrossPro(const Vec3f &v1, const Vec3f &v2) {
    return {v1.y * v2.z - v2.y * v1.z, v2.x * v1.z - v1.x * v2.z, v1.x * v2.y - v2.x * v1.y};
}

Vec3f vectorAdd(const Vec3f &v1, const Vec3f &v2) {
    return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

Vec3f vectorScalarPro(const Vec3f &v, const float i) {
    return {v.x * i, v.y * i, v.z * i};
}

Vec3f vectorSub(const Vec3f &v1, const Vec3f &v2) {
    return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

Ray laserBeams(const Camera &camera, const int i, const int j, const int &width, const int &height, const Vec3f &m) {
    Vec3f v = vectorNormalize(camera.up);
    Vec3f u = vectorCrossPro(vectorNormalize(camera.gaze), v);
    Vec3f q = vectorAdd(vectorAdd(m, vectorScalarPro(u, camera.near_plane.x)),vectorScalarPro(v, camera.near_plane.w));
    float su = (camera.near_plane.y - camera.near_plane.x) * (j + 0.5) / width;
    float sv = (camera.near_plane.w - camera.near_plane.z) * (i + 0.5) / height;
    Vec3f s = vectorSub(vectorAdd(q, vectorScalarPro(u, su)), vectorScalarPro(v, sv));
    Ray ray = {camera.position, vectorNormalize(vectorSub(s, camera.position))};
    return ray;
}

float determinant(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2)
{
	return v0.x * (v1.y * v2.z - v1.z * v2.y) + v0.y * (v1.z * v2.x - v1.x * v2.z) + v0.z * (v1.x * v2.y - v1.y * v2.x);
}

Blaster triangleShooting(const Ray &ray, const Scene &scene, const Triangle &triangle) {
    Vec3f v0 = scene.vertex_data[triangle.indices.v0_id - 1];
    Vec3f v1 = scene.vertex_data[triangle.indices.v1_id - 1];
    Vec3f v2 = scene.vertex_data[triangle.indices.v2_id - 1];
    Vec3f a_b = vectorSub(v0, v1);
    Vec3f a_c = vectorSub(v0, v2);
    Vec3f d = ray.direction;
    float determinantA = determinant(a_b, a_c, d);
    Blaster blaster;
    blaster.shot = false;
    if(determinantA == 0) return blaster;
    Vec3f a_o = vectorSub(v0, ray.start);
    float B = determinant(a_o, a_c, d) / determinantA;
    if(B < 0) return blaster;
    float Y = determinant(a_b, a_o, d) / determinantA;
    if(Y < 0) return blaster;
    if(B + Y > 1) return blaster;
    float t = determinant(a_b, a_c, a_o) / determinantA;
    if(t < 0) return blaster;
    blaster.t = t;
    blaster.shot = true;
    blaster.material_id = triangle.material_id;
    blaster.intersectionP = {ray.start.x + blaster.t * ray.direction.x, ray.start.y + blaster.t * ray.direction.y, ray.start.z + blaster.t * ray.direction.z};
    blaster.normal = vectorNormalize(vectorCrossPro(vectorSub(v1, v0), vectorSub(v2, v0)));
    return blaster;
}

float vectorDotPro(const Vec3f &v1, const Vec3f &v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Blaster sphereShooting(const Ray &ray, const Scene &scene, const Sphere &sphere) {
    Vec3f center = scene.vertex_data[sphere.center_vertex_id - 1];
    float radius = sphere.radius;
    Blaster blaster;
    Vec3f o_c = vectorSub(ray.start, center);
    float A = vectorDotPro(ray.direction, ray.direction);
    float B = 2*vectorDotPro(ray.direction, o_c);
    float C = vectorDotPro(o_c, o_c) - pow(radius, 2);
    float discriminant = B*B - 4*A*C;
    if(discriminant < 0) blaster.shot = false;
    else if (discriminant == 0) {
        blaster.t = -B / (2*A);
        blaster.shot = true;
        blaster.material_id = sphere.material_id;
        blaster.intersectionP = {ray.start.x + blaster.t * ray.direction.x, ray.start.y + blaster.t * ray.direction.y, ray.start.z + blaster.t * ray.direction.z};
        blaster.normal = vectorSub(blaster.intersectionP, center);
        blaster.normal = {blaster.normal.x/radius, blaster.normal.y/radius, blaster.normal.z/radius};
    }
    else {
        float t1 = (-B + sqrt(discriminant)) / (2 * A);
        float t2 = (-B - sqrt(discriminant)) / (2 * A);
        if(t1 <  t2) blaster.t = t1;
        else blaster.t = t2;
        blaster.shot = true;
        blaster.material_id = sphere.material_id;
        blaster.intersectionP = {ray.start.x + blaster.t * ray.direction.x, ray.start.y + blaster.t * ray.direction.y, ray.start.z + blaster.t * ray.direction.z};
        blaster.normal = vectorSub(blaster.intersectionP, center);
        blaster.normal = {blaster.normal.x/radius, blaster.normal.y/radius, blaster.normal.z/radius};
        // blaster.normal = vectorNormalize(blaster.normal);
    }
    return blaster;
}

Blaster meshShooting(const Ray &ray, const Scene &scene, const Mesh &mesh) {
    float tMin = FLT_MAX;
    Blaster blaster;
    blaster.shot = false;
    for(int i = 0; i <  mesh.faces.size(); i++) {
        Vec3f v0 = scene.vertex_data[mesh.faces[i].v0_id - 1];
        Vec3f v1 = scene.vertex_data[mesh.faces[i].v1_id - 1];
        Vec3f v2 = scene.vertex_data[mesh.faces[i].v2_id - 1];
        Vec3f a_b = vectorSub(v0, v1);
        Vec3f a_c = vectorSub(v0, v2);
        Vec3f d = ray.direction;
        float determinantA = determinant(a_b, a_c, d);
        if(determinantA == 0) continue;
        Vec3f a_o = vectorSub(v0, ray.start);
        float B = determinant(a_o, a_c, d) / determinantA;
        if(B < 0) continue; 
        float Y = determinant(a_b, a_o, d) / determinantA;
        if(Y < 0) continue;
        if(B + Y > 1) continue;
        float t = determinant(a_b, a_c, a_o) / determinantA;
        if(t < 0) continue;
        if(t < tMin) {
            tMin = t;
            blaster.t = t;
            blaster.shot = true;
            blaster.material_id = mesh.material_id;
            blaster.intersectionP = {ray.start.x + blaster.t * ray.direction.x, ray.start.y + blaster.t * ray.direction.y, ray.start.z + blaster.t * ray.direction.z};
            blaster.normal = vectorNormalize(vectorCrossPro(vectorSub(v1, v0), vectorSub(v2, v0)));
        }
    }
    return blaster;
}

Vec3f vectorPro(const Vec3f &v1, const Vec3f &v2) {
    return {v1.x * v2.x, v1.y * v2.y, v1.z * v2.z};
}

bool isShadow(const Ray &ray, const Scene &scene, const float &tLight) {
    Triangle triangle;
    Sphere sphere;
    Mesh mesh;
    for(int i = 0; i < scene.triangles.size(); i++) {
        triangle = scene.triangles[i];
        Blaster triangleShot = triangleShooting(ray, scene, triangle);
        // if(triangleShot.shot && triangleShot.t >= scene.shadow_ray_epsilon && triangleShot.t < tLight) return true;
        if(triangleShot.shot && triangleShot.t >= 0 && triangleShot.t < tLight) return true;
    }
    for(int i = 0; i < scene.spheres.size(); i++) {
        sphere = scene.spheres[i];
        Blaster sphereShot = sphereShooting(ray, scene, sphere);
        // if(sphereShot.shot && sphereShot.t >= scene.shadow_ray_epsilon && sphereShot.t < tLight) return true;
        if(sphereShot.shot && sphereShot.t >= 0 && sphereShot.t < tLight) return true;
    }
    for(int i = 0; i < scene.meshes.size(); i++) {
        mesh = scene.meshes[i];
        Blaster meshShot = meshShooting(ray, scene, mesh);
        // if(meshShot.shot && meshShot.t >= scene.shadow_ray_epsilon && meshShot.t < tLight) return true;
        if(meshShot.shot && meshShot.t >= 0 && meshShot.t < tLight) return true;
    }
    return false;
}

float vectorDistanceSquare(const Vec3f &v1, const Vec3f &v2) {
    return pow(v1.x - v2.x, 2) + pow(v1.y - v2.y, 2) + pow(v1.z - v2.z, 2);
}

Vec3f color(const Camera &camera, const Scene &scene, const Ray &ray, const int max_depth, bool main_ray) {
    float tMin = FLT_MAX;
    Vec3f pixel = {0, 0, 0};
    Triangle triangle;
    Sphere sphere;
    Mesh mesh;
    Blaster blaster;
    blaster.shot = false;
    for(int i = 0; i < scene.triangles.size(); i++) {
        triangle = scene.triangles[i];
        Blaster triangleShot = triangleShooting(ray, scene, triangle);
        if(triangleShot.shot && triangleShot.t >= 0 && triangleShot.t < tMin) {
            tMin = triangleShot.t;
            blaster = triangleShot;
        }
    }
    for(int i = 0; i < scene.spheres.size(); i++) {
        sphere = scene.spheres[i];
        Blaster sphereShot = sphereShooting(ray, scene, sphere);
        if(sphereShot.shot && sphereShot.t >= 0 && sphereShot.t < tMin) {
            tMin = sphereShot.t;
            blaster = sphereShot;
        }
    }
    for(int i = 0; i < scene.meshes.size(); i++) {
        mesh = scene.meshes[i];
        Blaster meshShot = meshShooting(ray, scene, mesh);
        if(meshShot.shot && meshShot.t >= 0 && meshShot.t < tMin) {
            tMin = meshShot.t;
            blaster = meshShot;
        }
    }
    if(blaster.shot) {
        pixel = vectorPro(scene.ambient_light, scene.materials[blaster.material_id - 1].ambient);
        bool shadow;
        for(int i = 0; i < scene.point_lights.size(); i++) {
            PointLight pL = scene.point_lights[i];
            Vec3f wi = vectorNormalize(vectorSub(pL.position, blaster.intersectionP));
            // Vec3f wi_s = vectorNormalize(vectorSub(pL.position, vectorAdd(blaster.intersectionP, vectorScalarPro(wi, scene.shadow_ray_epsilon))));
            // Ray lightRay = {vectorAdd(blaster.intersectionP, vectorScalarPro(wi, scene.shadow_ray_epsilon)), wi_s};
            Ray lightRay = {vectorAdd(blaster.intersectionP, vectorScalarPro(wi, scene.shadow_ray_epsilon)), wi};
            // Ray lightRay = {vectorAdd(blaster.intersectionP, vectorScalarPro(blaster.normal, scene.shadow_ray_epsilon)), wi};
            float tLight  = (pL.position.x - lightRay.start.x) / lightRay.direction.x;
            // shadow = isShadow(lightRay, scene, tLight);
            // if(shadow) continue;
            if(isShadow(lightRay, scene, tLight)) continue;
            Vec3f specularComponent = scene.materials[blaster.material_id - 1].specular;
            Vec3f halfVector = vectorNormalize(vectorSub(wi, ray.direction)); 
            float normal_h = vectorDotPro(blaster.normal, halfVector);
            normal_h = 0 <= normal_h ? normal_h : 0;
            float phong_exp = scene.materials[blaster.material_id - 1].phong_exponent;
            Vec3f irradiance = {0, 0, 0};
            float distance_2 = vectorDistanceSquare(pL.position, blaster.intersectionP);
            if(distance_2 != 0) irradiance = {pL.intensity.x / distance_2, pL.intensity.y / distance_2, pL.intensity.z / distance_2};
            Vec3f specular = {specularComponent.x * pow(normal_h, phong_exp) * irradiance.x, specularComponent.y * pow(normal_h, phong_exp) * irradiance.y, specularComponent.z * pow(normal_h, phong_exp) * irradiance.z};
            Vec3f diffuseComponent = scene.materials[blaster.material_id - 1].diffuse;
            float wi_n = vectorDotPro(wi, blaster.normal);
            wi_n = 0 <= wi_n ? wi_n : 0;
            Vec3f diffuse = {diffuseComponent.x * wi_n * irradiance.x, diffuseComponent.y * wi_n * irradiance.y, diffuseComponent.z * wi_n * irradiance.z};
            pixel = vectorAdd(pixel, vectorAdd(specular, diffuse));
        }
        if(scene.materials[blaster.material_id - 1].is_mirror && max_depth > 0) {
            bool isShot = false;
            Vec3f wr = vectorScalarPro(blaster.normal, -2 * vectorDotPro(blaster.normal, ray.direction));
            wr = vectorNormalize(vectorAdd(wr, ray.direction));
            Ray reflection = {vectorAdd(blaster.intersectionP, vectorScalarPro(wr, scene.shadow_ray_epsilon)), wr};
            // Ray reflection = {vectorAdd(blaster.intersectionP, vectorScalarPro(blaster.normal, scene.shadow_ray_epsilon)), wr};
            for(int j = 0; j < scene.triangles.size(); j++) {
                triangle = scene.triangles[j];
                Blaster triangleShot = triangleShooting(reflection, scene, triangle);
                if(triangleShot.shot && triangleShot.t >= 0) {
                    isShot = true;
                    break;
                }
            }
            if(!isShot) {
                for(int j = 0; j < scene.spheres.size(); j++) {
                    sphere = scene.spheres[j];
                    Blaster sphereShot = sphereShooting(reflection, scene, sphere);
                    if(sphereShot.shot && sphereShot.t >= 0) {
                        isShot = true;
                        break;
                    }
                }
            }
            if(!isShot) {
                for(int j = 0; j < scene.meshes.size(); j++) {
                    mesh = scene.meshes[j];
                    Blaster meshShot = meshShooting(reflection, scene, mesh);
                    if(meshShot.shot && meshShot.t >= 0) {
                        isShot = true;
                        break;
                    }
                }
            }
            if(isShot) {
                Vec3f mirror = color(camera, scene, reflection, max_depth - 1, false);
                pixel = vectorAdd(pixel, vectorPro(mirror, scene.materials[blaster.material_id - 1].mirror));
            }
            // Vec3f mirror = color(camera, scene, reflection, max_depth - 1, false);
            // pixel = vectorAdd(pixel, vectorPro(mirror, scene.materials[blaster.material_id - 1].mirror));
        }
    }
    else {
        if(main_ray) pixel = {(float)scene.background_color.x, (float)scene.background_color.y, (float)scene.background_color.z};
    }
    if(pixel.x > 255) pixel.x = 255;
    else pixel.x = floor(pixel.x);
    if(pixel.y > 255) pixel.y = 255;
    else pixel.y = floor(pixel.y);
    if(pixel.z > 255) pixel.z = 255;
    else pixel.z = floor(pixel.z);
    return pixel;
}

void imageGenerate(const Scene &scene, const Camera &camera, unsigned char *&image, const int &width, const int &height, const Vec3f &m) {
    int k = 0;
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            Ray ray = laserBeams(camera, i, j, width, height, m);
            Vec3f pixel = color(camera, scene, ray, scene.max_recursion_depth, true);
            image[k] = pixel.x;
            image[k+1] = pixel.y;
            image[k+2] = pixel.z;
            k += 3;
        }
    }
}

int main(int argc, char* argv[]) {
    Scene scene;
    scene.loadFromXml(argv[1]);
    int numOfCameras = scene.cameras.size();
    for(int i = 0; i < numOfCameras; i++) {
        Camera camera = scene.cameras[i];
        int width = camera.image_width;
        int height = camera.image_height;
        unsigned char *image = new unsigned char[width * height * 3];
        Vec3f m = vectorAdd(camera.position, vectorScalarPro(vectorNormalize(camera.gaze), camera.near_distance));
        imageGenerate(scene, camera, image, width, height, m);
        write_ppm(camera.image_name.c_str(), image, width, height);
        delete[] image;
    }
}
