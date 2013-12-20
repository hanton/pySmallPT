"""Introduction

pysmallpt, a Python path tracer by Hanton Yang, 2013
original smallpt, a path tracer by Kevin Beason, 2009
Style Guide for Python Code, PEP8, http://www.python.org/dev/peps/pep-0008/
Usage: pypy pysmallpt.py or python pysmallpt.py

"""

import sys
import time
from math import sqrt, pi, sin, cos
from random import random


class Vector(object):
    def __init__(self, x, y, z):
        self.x, self.y, self.z = float(x), float(y), float(z)

    def __str__(self):
        return "(%s,%s,%s)" % (self.x, self.y, self.z)

    def __repr__(self):
        return "Vector" + str(self)

    def __add__(self, other):
        return Vector(self.x+other.x, self.y+other.y, self.z+other.z)

    def __sub__(self, other):
        return Vector(self.x-other.x, self.y-other.y, self.z-other.z)

    def __mul__(self, other):
        if type(other) == Vector:
            # Dot Multiple
            return self.x*other.x + self.y*other.y + self.z*other.z
        else:
            # Scale
            return Vector(self.x*other, self.y*other, self.z*other)

    def __rmul__(self, other):
        if type(other) == Vector:
            # dot multiple
            return self.x*other.x + self.y*other.y + self.z*other.z
        else:
            # scale
            return Vector(self.x*other, self.y*other, self.z*other)

    def __div__(self, other):
        return Vector(self.x/other, self.y/other, self.z/other)

    def __xor__(self, other):
        return Vector(self.y*other.z - self.z*other.y,
                      self.z*other.x - self.x*other.z,
                      self.x*other.y - self.y*other.x)

    def length(self):
        return sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

    def normalize(self):
        return self / sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

class Color(object):
    def __init__(self, red, green, blue):
        self.r, self.g, self.b = red, green, blue

    def __str__(self):
        return "(%s,%s,%s)" % (self.r, self.g, self.b)

    def __repr__(self):
        return "Color" + str(self)

    def __add__(self, other):
        return Color(self.r+other.r, self.g+other.g, self.b+other.b)

    def __sub__(self, other):
        return Color(self.r-other.r, self.g-other.g, self.b-other.b)

    def __mul__(self, other):
        if type(other) == Color:
            return Color(self.r*other.r, self.g*other.g, self.b*other.b)
        else:
            return Color(self.r*other, self.g*other, self.b*other)

    def __rmul__(self, other):
        if type(other) == Color:
            return Color(self.r*other.r, self.g*other.g, self.b*other.b)
        else:
            return Color(self.r*other, self.g*other, self.b*other)

    def __div__(self, other):
        return Color(self.r/other, self.g/other, self.b/other)


class Ray(object):
    def __init__(self, origin, direction):
        self.origin = origin
        self.direction = direction


class Sphere(object):
    def __init__(self, radius, position, emission, color, material):
        self.radius = radius
        self.position = position
        self.emission = emission
        self.color = color
        self.material = material

    def intersect(self, ray):
        op = self.position - ray.origin
        eps = 1e-4
        b = op * ray.direction
        det = b*b - op*op + self.radius*self.radius

        if det < 0.0:
            return 0.0
        else:
            det = sqrt(det)

        t = b - det
        if t > eps:
            return t

        t = b + det
        if t > eps:
            return t

        return 0.0


spheres = []
# Left
radius = 1e5
position = Vector(1e5+1, 40.8, 81.6)
emission = Color(0.0, 0.0, 0.0)
color = Color(.75, .25, .25)
material = "diffuse"
sphere = Sphere(radius, position, emission, color, material)
spheres.append(sphere)
# Rght
radius = 1e5
position = Vector(-1e5+99, 40.8, 81.6)
emission = Color(0.0, 0.0, 0.0)
color = Color(.25, .25, .75)
material = "diffuse"
sphere = Sphere(radius, position, emission, color, material)
spheres.append(sphere)
# Back
radius = 1e5
position = Vector(50, 40.8, 1e5)
emission = Color(0.0, 0.0, 0.0)
color = Color(.75, .75, .75)
material = "diffuse"
sphere = Sphere(radius, position, emission, color, material)
spheres.append(sphere)
# Front
radius = 1e5
position = Vector(50, 40.8, -1e5+170)
emission = Color(0.0, 0.0, 0.0)
color = Color(0.0, 0.0, 0.0)
material = "diffuse"
sphere = Sphere(radius, position, emission, color, material)
spheres.append(sphere)
# Bottom
radius = 1e5
position = Vector(50, 1e5, 81.6)
emission = Color(0.0, 0.0, 0.0)
color = Color(.75, .75, .75)
material = "diffuse"
sphere = Sphere(radius, position, emission, color, material)
spheres.append(sphere)
# Top
radius = 1e5
position = Vector(50, -1e5+81.6, 81.6)
emission = Color(0.0, 0.0, 0.0)
color = Color(.75, .75, .75)
material = "diffuse"
sphere = Sphere(radius, position, emission, color, material)
spheres.append(sphere)
# Mirror
radius = 16.5
position = Vector(27, 16.5, 47)
emission = Color(0.0, 0.0, 0.0)
color = Color(1, 1, 1)*.999
material = "specular"
sphere = Sphere(radius, position, emission, color, material)
spheres.append(sphere)
# Glass
radius = 16.5
position = Vector(73, 16.5, 78)
emission = Color(0.0, 0.0, 0.0)
color = Color(1, 1, 1) * .999
material = "refractive"
sphere = Sphere(radius, position, emission, color, material)
spheres.append(sphere)
# Light
radius = 600
position = Vector(50, 681.6-0.27, 81.6)
emission = Color(12, 12, 12)
color = Color(0.0, 0.0, 0.0)
material = "diffuse"
sphere = Sphere(radius, position, emission, color, material)
spheres.append(sphere)


def clamp(x):
    if x < 0.0:
        return 0.0
    else:
        if x > 1.0:
            return 1.0
        else:
            return x


def toInt(x):
    return chr(int(pow(clamp(x), 1.0 / 2.2)*255.0 + 0.5))


def intersect(ray):
    inf = t = 1e20
    d = 0.0
    hit_object = None
    for sphere in spheres:
        d = sphere.intersect(ray)
        if d != 0.0 and d < t:
            t = d
            hit_object = sphere
    hit = (t < inf)
    return hit, t, hit_object


def radiance(ray, depth):
    hit, t, hit_object = intersect(ray)

    if hit is False:
        return Color(0.0, 0.0, 0.0)

    x = ray.origin + t*ray.direction
    n = (x - hit_object.position).normalize()
    if n*ray.direction < 0.0:
        nl = n
    else:
        nl = -1.0 * n
    f = hit_object.color

    p = max(f.r, f.g, f.b)

    depth = depth + 1
    if depth > 5:
        if random() < p:
            f = f * (1.0 / p)
        else:
            return hit_object.emission

    if hit_object.material == 'diffuse':
        r1 = 2.0 * pi * random()
        r2 = random()
        r2s = sqrt(r2)
        w = nl
        if abs(w.x) > 0.1:
            u = Vector(0.0, 1.0, 0.0) ^ w
        else:
            u = Vector(1.0, 1.0, 1.0) ^ w
        u = u.normalize()
        v = w ^ u
        d = u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1.0 - r2)
        d = d.normalize()
        return hit_object.emission + f * radiance(Ray(x, d), depth)
    elif hit_object.material == 'specular':
        ray_direction = ray.direction - n*2.0*(n * ray.direction)
        reflect_ray = Ray(x, ray_direction)
        return hit_object.emission + f*radiance(reflect_ray, depth)

    # Ideal Dielectric Refraction
    reflect_ray_direction = ray.direction - n*2.0*n*ray.direction
    reflect_ray = Ray(x, reflect_ray_direction)
    into = (n*nl > 0.0)
    nc = 1.0
    nt = 1.5
    if into:
        nnt = nc / nt
    else:
        nnt = nt / nc
    ddn = ray.direction * nl
    cos2t = 1.0 - nnt*nnt*(1.0 - ddn*ddn)
    if cos2t < 0.0:
        return hit_object.emission + f*radiance(reflect_ray, depth)
    if into:
        tdir = ray.direction*nnt - n*(ddn*nnt + sqrt(cos2t))
    else:
        tdir = ray.direction*nnt - (-1.0)*n*(ddn*nnt + sqrt(cos2t))
    tdir = tdir.normalize()
    a = nt - nc
    b = nt + nc
    R0 = (a*a) / (b*b)
    if into:
        c = 1.0 - (-ddn)
    else:
        c = 1.0 - (tdir*n)
    Re = R0 + (1.0-R0) * pow(c, 5)
    Tr = 1.0 - Re
    P = 0.25 + 0.5*Re
    RP = Re / P
    TP = Tr / (1.0-P)

    # Russian Roulette
    if depth > 2:
        if random() < P:
            return hit_object.emission + f*radiance(reflect_ray, depth)*RP
        else:
            return hit_object.emission + f*radiance(Ray(x, tdir), depth)*TP
    else:
        f1 = radiance(reflect_ray, depth)
        f2 = radiance(Ray(x, tdir), depth)
        return  f1 * Re + f2 * Tr


def smallpt():
    start_time = time.time()

    width, height = 1024, 768
    samples = 1
    if len(sys.argv) == 2:
        samples = int(sys.argv[1]) / 4

    camera_origin = Vector(50, 52, 295.6)
    camera_direction = Vector(0, -0.042612, -1).normalize()
    cx = width * 0.5135 / height
    cx = Vector(cx, 0.0, 0.0)
    cy = (cx ^ camera_direction).normalize() * 0.5135
    pixels = []
    for y in xrange(height - 1, -1, -1):
        text = "\r" + "Remain: {0:.2f}%".format(100.0 * y / (height-1))
        sys.stdout.write(text)
        sys.stdout.flush()
        for x in xrange(width):
            L = Color(0.0, 0.0, 0.0)
            for sy in xrange(2):
                for sx in xrange(2):
                    r = Color(0.0, 0.0, 0.0)
                    for s in xrange(samples):
                        r1 = 2 * random()
                        if r1 < 1.0:
                            dx = sqrt(r1) - 1
                        else:
                            dx = 1.0 - sqrt(2.0-r1)
                        r2 = 2 * random()
                        if r2 < 1.0:
                            dy = sqrt(r2) - 1
                        else:
                            dy = 1.0 - sqrt(2.0-r2)
                        c_x = cx * (((sx+0.5+dx)/2.0 + x)/width - 0.5)
                        c_y = cy * (((sy+0.5+dy)/2.0 + y)/height - 0.5)
                        d = c_x + c_y + camera_direction
                        ray_origin = camera_origin + d * 135.0
                        ray_direction = d.normalize()
                        r += radiance(Ray(ray_origin, ray_direction), 0)
                    r /= samples
                    L += r
            L /= 4.0
            pixels.append(L)

    with open('pySmallPT.ppm', 'wb') as f:
        f.write('P6 %d %d 255\n' % (width, height))
        for pixel in pixels:
            f.write(toInt(pixel.r) + toInt(pixel.g) + toInt(pixel.b))
        f.close()

    end_time = time.time()
    print "\rRender Time: %i seconds." % (end_time - start_time)


if __name__ == '__main__':
    smallpt()
