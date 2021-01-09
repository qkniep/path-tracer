// Copyright (C) 2020 Quentin M. Kniep <hello@quentinkniep.com>
// Distributed under terms of the MIT license.

package main

import (
	"bufio"
	"fmt"
	"math"
	"math/rand"
	"os"
	"time"
)

func clamp(x float64) float64 {
	return math.Max(0.0, math.Min(1.0, x))
}

func toInt(x float64) int {
	return int(math.Pow(clamp(x), 1.0/2.2) * 255.0 + 0.5)
}

type refl int
const (
	DIFF refl = iota
	SPEC
	REFR
)

var spheres = []Sphere{ //Scene: radius, position, emission, color, material
	{1e5,  Vec3d{ 1e5+1,40.8,81.6},  Vec(),           Vec3d{.75,.25,.25},      DIFF},//Left
	{1e5,  Vec3d{-1e5+99,40.8,81.6}, Vec(),           Vec3d{.25,.25,.75},      DIFF},//Rght
	{1e5,  Vec3d{50,40.8, 1e5},      Vec(),           Vec3d{.75,.75,.75},      DIFF},//Back
	{1e5,  Vec3d{50,40.8,-1e5+170},  Vec(),           Vec3d{.75,.75,.75},      DIFF},//Frnt
	{1e5,  Vec3d{50, 1e5, 81.6},     Vec(),           Vec3d{.75,.75,.75},      DIFF},//Botm
	{1e5,  Vec3d{50,-1e5+81.6,81.6}, Vec(),           Vec3d{.75,.75,.75},      DIFF},//Top
	{16.5, Vec3d{27,16.5,47},        Vec(),           Vec3d{1,1,1}.mul(0.999), SPEC},//Mirr
	{16.5, Vec3d{73,16.5,78},        Vec(),           Vec3d{1,1,1}.mul(0.999), REFR},//Glas
	{600,  Vec3d{50,681.6-.27,81.6}, Vec3d{12,12,12}, Vec(),                   DIFF},//Lite
	//{1.5,  Vec3d{50,81.6-16.5,81.6}, Vec3d{4,4,4}.mul(100),Vec(), DIFF},//Lite
}

// Find the closest object the ray intersects.
// Sets t to the distance, id to the objects index, and returns true iff there was any object.
func intersect(r Ray, t *float64, id *int) bool {
	var d, inf = 0.0, math.Inf(1)
	*t = inf
	for i, s := range spheres {
		d = s.intersect(r)
		if d != 0 && d < *t { // found closer object than before
			*t, *id = d, i
		}
	}
	return !math.IsInf(*t, 1)
}

type rowRender struct {
	y    int
	line []Vec3d
}

func main() {
	var width, height, samples = 1280, 720, 512
	var cam = Ray{Vec3d{50, 52, 295.6}, Vec3d{0, -0.042612, -1}.norm()}  // cam pos, dir
	var cos = math.Cos((62.5/360.0)*(2.0*math.Pi))
	var cx = Vec3d{float64(width)*cos/float64(height), 0, 0}
	var cy = cx.cross(cam.direction).norm().mul(cos)
	var image = make([]Vec3d, width*height);

	rowChannel := make(chan rowRender, 8)
	for offset := 0; offset < 4; offset++ {
		go renderImageRowset(width, height, samples, offset, 4, cam, cx, cy, rowChannel)
	}

	// collect rendered rows
	var progress = 0
	for progress < height {
		rr := <- rowChannel
		progress++
		copy(image[(height-rr.y-1)*width:], rr.line)
		fmt.Printf("\rRendering (%d spp) %.3f%%", samples*4, 100.0*float64(progress)/float64(height));
	}

	f, err := os.Create("/tmp/image.ppm")
	if err != nil {
		panic(err)
	}
	defer f.Close()
	wrt := bufio.NewWriter(f)
	fmt.Fprintf(wrt, "P3\n%d %d\n%d\n", width, height, 255)
	for i := 0; i < len(image); i++ {
		fmt.Fprintf(wrt, "%d %d %d ", toInt(image[i].X), toInt(image[i].Y), toInt(image[i].Z))
	}
	wrt.Flush()
}

func renderImageRowset(w, h, samples, yInit, yStep int, cam Ray, cx, cy Vec3d, rowChannel chan rowRender) {
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	var r = Vec()

	for y := yInit; y < h; y += yStep {
		line := make([]Vec3d, w)
		for x := 0; x < w; x++ {
			// use 2x2 subpixels
			for sy := 0; sy < 2; sy++ {
				for sx := 0; sx < 2; sx++ {
					for s := 0; s < samples; s++ {
						r1, r2 := 2.0 * rng.Float64(), 2.0 * rng.Float64()
						dx, dy := math.Sqrt(r1) - 1.0, math.Sqrt(r2) - 1.0
						if r1 >= 1.0 {
							dx = 1.0 - math.Sqrt(2.0-r1)
						}
						if r2 >= 1.0 {
							dy = 1.0 - math.Sqrt(2.0-r2)
						}
						d := cam.direction.add(
							cx.mul(((float64(sx)+0.5+dx)/2.0 + float64(x))/float64(w) - 0.5)).add(
							cy.mul(((float64(sy)+0.5+dy)/2.0 + float64(y))/float64(h) - 0.5))
						r = r.add(radiance(Ray{cam.origin.add(d.mul(140.0)), d.norm()}, 0, rng, 1.0))
					} // Camera rays are pushed ^^^^^ forward to start in interior
					r = r.mul(1.0/float64(samples))
					line[x] = line[x].add(Vec3d{clamp(r.X), clamp(r.Y), clamp(r.Z)}.mul(0.25))
					r = Vec()
				}
			}
		}
		rowChannel <- rowRender{y, line}
	}
}

func radiance(r Ray, depth int, rng *rand.Rand, E float64) Vec3d {
	var distance float64  // distance to intersection
	var objID = 0     // id of intersected object
	if !intersect(r, &distance, &objID) {
		return Vec() // if miss, return black
	}
	var obj = spheres[objID];  // the hit object

	x := r.origin.add(r.direction.mul(distance))
	n, f := (x.sub(obj.pos)).norm(), obj.c
	nl := n
	if n.dot(r.direction) >= 0 {
		nl = n.mul(-1)
	}
	w := nl
	p := math.Max(f.X, math.Max(f.Y, f.Z))
	depth++
	if depth > 5 {
		if rng.Float64() < p {
			f = f.mul(1.0/p)
		} else {
			return obj.e.mul(E)  //R.R.
		}
	}
	switch (obj.refl) {
	case DIFF:  // ideal DIFFUSE reflection
		r1, r2 := 2.0*math.Pi*rng.Float64(), rng.Float64()
		r2s := math.Sqrt(r2)
		dir := Vec3d{1, 0, 0}
		if math.Abs(w.X) > 0.1 {
			dir = Vec3d{0, 1, 0}
		}
		u := dir.cross(w).norm()
		v := w.cross(u)
		d := u.mul(math.Cos(r1)*r2s).add(v.mul(math.Sin(r1)*r2s)).add(w.mul(math.Sqrt(1-r2))).norm()

		/*// Loop over any lights
		var e Vec3d
		for i, s := range spheres {
			if s.e.X <= 0 && s.e.Y <= 0 && s.e.Z <= 0 {
				continue // skip non-lights
			}

			sw := s.pos.sub(x)
			dir := Vec3d{1, 0, 0}
			if math.Abs(sw.X) > 0.1 {
				dir = Vec3d{0, 1, 0}
			}
			su := dir.norm()
			sv := sw.cross(su)
			cos_a_max := math.Sqrt(1.0-s.radius*s.radius/(x.sub(s.pos)).dot(x.sub(s.pos)))
			eps1, eps2 := rng.Float64(), rng.Float64()
			cos_a := 1.0-eps1+eps1*cos_a_max
			sin_a := math.Sqrt(1.0-cos_a*cos_a)
			phi := 2.0*math.Pi*eps2
			l := su.mul(math.Cos(phi)*sin_a).add(sv.mul(math.Sin(phi)*sin_a)).add(sw.mul(cos_a)).norm()
			if intersect(Ray{x,l}, &t, &id) && id == i {  // shadow ray
				omega := 2.0*math.Pi*(1.0-cos_a_max)
				e = e.add(f.compMul(s.e.mul(l.dot(nl)*omega)).mul(1.0/math.Pi))  // 1/pi for brdf
			}
		}

		return obj.e.mul(E).add(e).add(f.compMul(radiance(Ray{x,d},depth,rng,0.0)))*/
		return obj.e.add(f.compMul(radiance(Ray{x, d}, depth, rng, 0.0)))
	case SPEC:
		reflRay := Ray{x, r.direction.sub(n.mul(2*n.dot(r.direction)))}  // ideal SPECULAR reflection
		return obj.e.add(f.compMul(radiance(reflRay, depth, rng, 1.0)))
	case REFR:
		reflRay := Ray{x, r.direction.sub(n.mul(2*n.dot(r.direction)))}  // ideal dielectric REFRACTION
		into := n.dot(nl) > 0  // ray from outside going in?
		nc, nt := 1.0, 1.5
		nnt, ddn := nt/nc, r.direction.dot(nl)
		if into {
			nnt = nc/nt
		}
		cos2t := 1-nnt*nnt*(1-ddn*ddn)  // total internal reflection
		if cos2t < 0 {
			return obj.e.add(f.compMul(radiance(reflRay, depth, rng, 1.0)))
		}
		scalar := -1.0
		if into {
			scalar = 1.0
		}
		tdir := r.direction.mul(nnt).sub(n.mul(scalar*(ddn*nnt+math.Sqrt(cos2t)))).norm()
		a, b := nt-nc, nt+nc
		R0, c := a*a/(b*b), 1-tdir.dot(n)
		if into {
			c = 1+ddn
		}
		cSquare := c*c
		Re := R0+(1-R0)*cSquare*cSquare
		Tr := 1-Re
		P := 0.25 + 0.5*Re
		RP, TP := Re/P, Tr/(1-P)

		if depth <= 2 {
			return obj.e.add(f.compMul(radiance(reflRay, depth, rng, 1.0).mul(Re).add(radiance(Ray{x, tdir}, depth, rng, 1.0).mul(Tr))))
		} else if rng.Float64()<P {
			return obj.e.add(f.compMul(radiance(reflRay, depth, rng, 1.0).mul(RP)))
		}
		return obj.e.add(f.compMul(radiance(Ray{x, tdir}, depth, rng, 1.0).mul(TP)))
	}
	panic("Unsupported material!")
}
