// Copyright (C) 2020 Quentin M. Kniep <hello@quentinkniep.com>
// Distributed under terms of the MIT license.

package main

import "math"

type Vec3d struct {
	X, Y, Z float64
}

func Vec() Vec3d {
	return Vec3d{0, 0, 0}
}

func (v Vec3d) mul(s float64) Vec3d {
	return Vec3d{s*v.X, s*v.Y, s*v.Z}
}

func (v1 Vec3d) add(v2 Vec3d) Vec3d {
	return Vec3d{v1.X+v2.X, v1.Y+v2.Y, v1.Z+v2.Z}
}

func (v1 Vec3d) sub(v2 Vec3d) Vec3d {
	return Vec3d{v1.X-v2.X, v1.Y-v2.Y, v1.Z-v2.Z}
}

func (v1 Vec3d) compMul(v2 Vec3d) Vec3d {
	return Vec3d{v1.X*v2.X, v1.Y*v2.Y, v1.Z*v2.Z}
}

func (v1 Vec3d) dot(v2 Vec3d) float64 {
	return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z
}

func (v1 Vec3d) cross(v2 Vec3d) Vec3d {
	return Vec3d{v1.Y*v2.Z-v1.Z*v2.Y, v1.Z*v2.X-v1.X*v2.Z, v1.X*v2.Y-v1.Y*v2.X}
}

func (v Vec3d) norm() Vec3d {
	return v.mul(1/math.Sqrt(v.dot(v)))
}

type Ray struct {
	origin, direction Vec3d
}

type Sphere struct {
	radius  float64
	pos, e, c Vec3d   // position, emission, color
	refl Refl_t     // reflection type (DIFFuse, SPECular, REFRactive)
};

// Solves t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-radius^2 = 0.
// Returns distance from ray's origin, 0 if the sphere is not hit by the ray.
func (s Sphere) intersect(r Ray) float64 {
	var op = s.pos.sub(r.origin);
	var eps, b = 1e-4, op.dot(r.direction)
	var det = b*b - op.dot(op) + s.radius*s.radius

	if det < 0 {
		return 0
	}
	det = math.Sqrt(det)
	if b-det > eps {
		return b-det
	}
	if b+det > eps {
		return b+det
	}
	return 0
}
