package s2

import (
	"log"
	"math"
	"math/rand"
	"strconv"
	"strings"
	"testing"

	"github.com/davidreynolds/gos2/r1"
	"github.com/davidreynolds/gos2/r3"
	"github.com/davidreynolds/gos2/s1"
)

func parsePoints(s string) (vertices []Point) {
	if s == "" {
		return []Point{}
	}
	points := strings.Split(s, ",")
	for _, p := range points {
		vertices = append(vertices, makepoint(p))
	}
	return
}

func makepoint(s string) Point {
	s = strings.TrimSpace(s)
	degs := strings.Split(s, ":")
	lat, _ := strconv.ParseFloat(degs[0], 64)
	lng, _ := strconv.ParseFloat(degs[1], 64)
	ll := LatLngFromDegrees(lat, lng)
	return PointFromLatLng(ll)
}

func makeloop(s string) *Loop {
	path := parsePoints(s)
	return NewLoopFromPath(path)
}

var (
	// The northern hemisphere, defined using two pairs of antipodal points.
	north_hemi = makeloop("0:-180, 0:-90, 0:0, 0:90")

	// The northern hemisphere, defined using three points 120 degrees apart.
	north_hemi3 = makeloop("0:-180, 0:-60, 0:60")

	// The southern hemisphere, defined using two pairs of antipodal points.
	south_hemi = makeloop("0:90, 0:0, 0:-90, 0:-180")

	// The western hemisphere, defined using two pairs of antipodal points.
	west_hemi = makeloop("0:-180, -90:0, 0:0, 90:0")

	// The eastern hemisphere, defined using two pairs of antipodal points.
	east_hemi = makeloop("90:0, 0:0, -90:0, 0:-180")

	// The "near" hemisphere, defined using two pairs of antipodal points.
	near_hemi = makeloop("0:-90, -90:0, 0:90, 90:0")

	// The "far" hemisphere, defined using two pairs of antipodal points.
	far_hemi = makeloop("90:0, 0:90, -90:0, 0:-90")

	// A spiral stripe that slightly over-wraps the equator.
	candy_cane = makeloop("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70")

	// A small clockwise loop in the northern & eastern hemispheres.
	small_ne_cw = makeloop("35:20, 45:20, 40:25")

	// Loop around the north pole at 80 degrees.
	arctic_80 = makeloop("80:-150, 80:-30, 80:90")

	// Loop around the south pole at 80 degrees.
	antarctic_80 = makeloop("-80:120, -80:0, -80:-120")

	// A completely degenerate triangle along the equator that RobustCCW()
	// considers to be CCW.
	line_triangle = makeloop("0:1, 0:3, 0:2")

	// A nearly-degenerate CCW chevron near the equator with very long sides
	// (about 80 degrees).  Its area is less than 1e-640, which is too small
	// to represent in double precision.
	skinny_chevron = makeloop("0:0, -1e-320:80, 0:1e-320, 1e-320:80")

	// A diamond-shaped loop around the point 0:180.
	loop_a = makeloop("0:178, -1:180, 0:-179, 1:-180")

	// Another diamond-shaped loop around the point 0:180.
	loop_b = makeloop("0:179, -1:180, 0:-178, 1:-180")

	// The intersection of A and B.
	a_intersect_b = makeloop("0:179, -1:180, 0:-179, 1:-180")

	// The union of A and B.
	a_union_b = makeloop("0:178, -1:180, 0:-178, 1:-180")

	// A minus B (concave).
	a_minus_b = makeloop("0:178, -1:180, 0:179, 1:-180")

	// B minus A (concave).
	b_minus_a = makeloop("0:-179, -1:180, 0:-178, 1:-180")

	// A shape gotten from a by adding one triangle to one edge, and
	// subtracting another triangle on an opposite edge.
	loop_c = makeloop("0:178, 0:180, -1:180, 0:-179, 1:-179, 1:-180")

	// A shape gotten from a by adding one triangle to one edge, and
	// adding another triangle on an opposite edge.
	loop_d = makeloop("0:178, -1:178, -1:180, 0:-179, 1:-179, 1:-180")
)

func TestGetRectBound(t *testing.T) {
	if !candy_cane.bound.Lng.IsFull() {
		t.Errorf("%v.IsFull() == false", candy_cane.bound.Lng)
	}
	deg := candy_cane.bound.Lo().Lat.Degrees()
	if deg >= -20 {
		t.Errorf("%v >= -20", deg)
	}
	deg = candy_cane.bound.Hi().Lat.Degrees()
	if deg <= 10 {
		t.Errorf("%v <= 10", deg)
	}

	if !small_ne_cw.bound.IsFull() {
		t.Errorf("%v.IsFull() == false", small_ne_cw.bound)
	}

	var p1, p2 LatLng
	var rect Rect

	p1 = LatLngFromDegrees(80, -180)
	p2 = LatLngFromDegrees(90, 180)
	rect = Rect{
		Lat: r1.Interval{p1.Lat.Radians(), p2.Lat.Radians()},
		Lng: s1.Interval{p1.Lng.Radians(), p2.Lng.Radians()},
	}

	if !arctic_80.bound.Equal(rect) {
		t.Errorf("%v.Equal(%v) == false", arctic_80.bound, rect)
	}

	p1 = LatLngFromDegrees(-90, -180)
	p2 = LatLngFromDegrees(-80, 180)
	rect = Rect{
		Lat: r1.Interval{p1.Lat.Radians(), p2.Lat.Radians()},
		Lng: s1.Interval{p1.Lng.Radians(), p2.Lng.Radians()},
	}

	if !antarctic_80.bound.Equal(rect) {
		t.Errorf("%v.Equal(%v) == false", antarctic_80.bound, rect)
	}

	// Create a loop that contains the complement of the "arctic_80" loop.
	arctic_80_inv := arctic_80.Clone()
	arctic_80_inv.Invert()
	// The highest altitude of each edge is attained at its midpoint
	mid := arctic_80_inv.vertex(0).Add(arctic_80_inv.vertex(1).Vector).Mul(0.5)
	want := arctic_80_inv.bound.Hi().Lat.Radians()
	got := LatLngFromPoint(Point{mid}).Lat.Radians()
	if math.Abs(got-want) > 1e-14 {
		t.Errorf("%v != %v", want, got)
	}

	if !south_hemi.bound.Lng.IsFull() {
		t.Errorf("%v.IsFull() == false", south_hemi.bound.Lng)
	}

	i := r1.Interval{-math.Pi / 2, 0}
	if !south_hemi.bound.Lat.Equal(i) {
		t.Errorf("%v.Equal(%v) == false", south_hemi.bound.Lat, i)
	}
}

func TestAreaAndCentroid(t *testing.T) {
	if got := north_hemi.Area(); math.Abs(got-2*math.Pi) > 1e-15 {
		t.Errorf("%v.Area() == %v, want %v", north_hemi, got, 2*math.Pi)
	}
	if got := east_hemi.Area(); got > 2*math.Pi+1e-12 {
		t.Errorf("%v.Area() > %v, want <= %v", east_hemi, got, 2*math.Pi+1e-12)
	}
	if got := east_hemi.Area(); got < 2*math.Pi-1e-12 {
		t.Errorf("%v.Area() > %v, want >= %v", east_hemi, got, 2*math.Pi-1e-12)
	}

	// Construct spherical caps of random height, and approximate their
	// boundary with closely spaced vertices. Then check that the area and
	// centroid are correct.
	const kMaxDist = 1e-6
	for i := 0; i < 100; i++ {
		// Choose a coordinate frame for the spherical cap.
		x, y, z := randomFrame()
		// Given two points at latitude phi and whose longitudes differ
		// by dtheta, the geodesic between the two points has a maximum
		// latitude of atan(tan(phi) / cos(dtheta/2)). This can be
		// derived by positioning the two points at (-dtheta/2, phi)
		// and (dtheta/2, phi).
		//
		// We want to position the vertices close enough together so
		// that their maximum distance from the boundary of the
		// spherical cap is kMaxDist. Thus we want
		// fabs(atan(tan(phi) / cos(dtheta/2)) - phi) <= kMaxDist.
		height := 2 * rand.Float64()
		phi := math.Asin(1 - height)
		max_dtheta := 2 * math.Acos(math.Tan(math.Abs(phi))/math.Tan(math.Abs(phi)+kMaxDist))
		max_dtheta = math.Min(math.Pi, max_dtheta)
		var vertices []Point
		for theta := 0.0; theta < 2*math.Pi; theta += rand.Float64() * max_dtheta {
			a := x.Mul(math.Cos(theta) * math.Cos(phi))
			b := y.Mul(math.Sin(theta) * math.Cos(phi))
			c := z.Mul(math.Sin(phi))
			vertices = append(vertices, Point{a.Add(b).Add(c)})
		}
		loop := NewLoopFromPath(vertices)
		area := loop.Area()
		centroid := loop.Centroid()
		expectedArea := 2 * math.Pi * height
		if got := math.Abs(area - expectedArea); got > 2*math.Pi*kMaxDist {
			t.Errorf("%v > %v", got, 2*math.Pi*kMaxDist)
		}
		if got := math.Abs(area - expectedArea); got < 0.01*kMaxDist {
			t.Errorf("%v < %v", got, 0.01*kMaxDist)
		}
		expectedCentroid := z.Mul(expectedArea * (1 - 0.5*height))
		if got := centroid.Sub(expectedCentroid).Norm(); got > 2*kMaxDist {
			t.Errorf("(%v - %v).Norm() == %v, want <= %v", centroid,
				expectedCentroid, got, 2*kMaxDist)
		}
	}
}

func TestTurningAngle(t *testing.T) {
	got := north_hemi3.TurningAngle()
	if math.Abs(got) > 1e-15 {
		t.Errorf("%v.TurningAngle() == %v, want near 0", north_hemi3, got)
	}
	checkTurningAngleInvariants(t, north_hemi3)
	checkTurningAngleInvariants(t, candy_cane)

	want := -2 * math.Pi
	got = line_triangle.TurningAngle()
	if math.Abs(got-want) > 1e-15 {
		t.Errorf("%v.TurningAngle() = %v, want %v", line_triangle, got, want)
	}
	checkTurningAngleInvariants(t, line_triangle)

	want = 2 * math.Pi
	got = skinny_chevron.TurningAngle()
	if math.Abs(got-want) > 1e-15 {
		t.Errorf("%v.TurningAngle() = %v, want %v", skinny_chevron, got, want)
	}
	checkTurningAngleInvariants(t, skinny_chevron)
}

func TestContains(t *testing.T) {
	if got := candy_cane.Contains(PointFromLatLng(LatLngFromDegrees(5, 71))); got != true {
		t.Errorf("%v.Contains(%v) = %v, want %v", candy_cane,
			PointFromLatLng(LatLngFromDegrees(5, 71)), got, true)
	}
	northCopy := north_hemi.Clone()
	southCopy := south_hemi.Clone()
	westCopy := west_hemi.Clone()
	eastCopy := east_hemi.Clone()
	tests := []struct {
		loop    **Loop
		x, y, z float64
		want    bool
	}{
		{&northCopy, 0, 0, 1, true},
		{&northCopy, 0, 0, -1, false},
		{&southCopy, 0, 0, 1, false},
		{&southCopy, 0, 0, -1, true},
		{&westCopy, 0, 1, 0, false},
		{&westCopy, 0, -1, 0, true},
		{&eastCopy, 0, 1, 0, true},
		{&eastCopy, 0, -1, 0, false},
	}
	for i := 0; i < 4; i++ {
		for _, test := range tests {
			p := Point{r3.Vector{test.x, test.y, test.z}}
			if got := (*test.loop).Contains(p); got != test.want {
				t.Errorf("%v.Contains(%v) = %v, want %v", *test.loop, p, got, test.want)
			}
		}
		loopRotate(&northCopy)
		loopRotate(&southCopy)
		loopRotate(&westCopy)
		loopRotate(&eastCopy)
	}

	for level := 0; level < 3; level++ {
		loops := []*Loop{}
		loopVertices := []Point{}
		points := map[Point]bool{}
		for id := CellIDBegin(level); id != CellIDEnd(level); id = id.Next() {
			cell := CellFromCellID(id)
			points[cell.Center()] = true
			for k := 0; k < 4; k++ {
				loopVertices = append(loopVertices, cell.Vertex(k))
				points[cell.Vertex(k)] = true
			}
			loops = append(loops, NewLoopFromPath(loopVertices))
			loopVertices = []Point{}
		}
		for p := range points {
			count := 0
			for j := 0; j < len(loops); j++ {
				if loops[j].Contains(p) {
					count++
				}
			}
			if count != 1 {
				t.Errorf("%v should be 1", count)
			}
		}
		for i, _ := range loops {
			loops[i] = nil
		}
		loops = []*Loop{}
	}
}

func checkTurningAngleInvariants(t *testing.T, loop *Loop) {
	want := loop.TurningAngle()
	loopCopy := loop.Clone()
	for i := 0; i < len(loop.vertices); i++ {
		loopCopy.Invert()
		if got := loopCopy.TurningAngle(); got != -want {
			t.Errorf("%v.TurningAngle() = %v, want %v", loopCopy, got, -want)
		}
		loopCopy.Invert()
		loopRotate(&loopCopy)
		if got := loopCopy.TurningAngle(); got != want {
			t.Errorf("%v.TurningAngle() = %v, want %v", loopCopy, got, want)
		}
	}
}

func loopRotate(loop **Loop) {
	vertices := make([]Point, 0, len((*loop).vertices))
	for i := 1; i <= len((*loop).vertices); i++ {
		vertices = append(vertices, *(*loop).vertex(i))
	}
	*loop = NewLoopFromPath(vertices)
}

func checkRelation(t *testing.T, a, b *Loop, containsOrCrosses int, intersects, nestable bool) {
	if got := a.ContainsLoop(b); got != (containsOrCrosses == 1) {
		t.Errorf("%v.ContainsLoop(%v) = %v, want %v", a, b, got, containsOrCrosses == 1)
	}
	if got := a.Intersects(b); got != intersects {
		t.Errorf("%v.Intersects(%v) = %v, want %v", a, b, got, intersects)
	}
	if nestable {
		if got := a.ContainsNested(b); got != a.ContainsLoop(b) {
			t.Errorf("%v.ContainsNested(%v) = %v, want %v", a, b, got, a.ContainsLoop(b))
		}
	}
	if containsOrCrosses >= -1 {
		if got := a.ContainsOrCrosses(b); got != containsOrCrosses {
			t.Errorf("%v.ContainsOrCrosses(%v) = %v, want %v", a, b, got, containsOrCrosses)
		}
	}
}

func TestLoopRelations(t *testing.T) {
	tests := []struct {
		a                 *Loop
		b                 *Loop
		containsOrCrosses int
		intersects        bool
		nestable          bool
	}{

		{north_hemi, north_hemi, 1, true, false},
		{north_hemi, south_hemi, 0, false, false},
		{north_hemi, east_hemi, -1, true, false},
		{north_hemi, arctic_80, 1, true, true},
		{north_hemi, antarctic_80, 0, false, true},
		{north_hemi, candy_cane, -1, true, false},

		// We can't compare north_hemi3 vs. north_hemi or south_hemi.
		{north_hemi3, north_hemi3, 1, true, false},
		{north_hemi3, east_hemi, -1, true, false},
		{north_hemi3, arctic_80, 1, true, true},
		{north_hemi3, antarctic_80, 0, false, true},
		{north_hemi3, candy_cane, -1, true, false},

		{south_hemi, north_hemi, 0, false, false},
		{south_hemi, south_hemi, 1, true, false},
		{south_hemi, far_hemi, -1, true, false},
		{south_hemi, arctic_80, 0, false, true},
		{south_hemi, antarctic_80, 1, true, true},
		{south_hemi, candy_cane, -1, true, false},

		{candy_cane, north_hemi, -1, true, false},
		{candy_cane, south_hemi, -1, true, false},
		{candy_cane, arctic_80, 0, false, true},
		{candy_cane, antarctic_80, 0, false, true},
		{candy_cane, candy_cane, 1, true, false},

		{near_hemi, west_hemi, -1, true, false},

		{small_ne_cw, south_hemi, 1, true, false},
		{small_ne_cw, west_hemi, 1, true, false},
		{small_ne_cw, north_hemi, -2, true, false},
		{small_ne_cw, east_hemi, -2, true, false},

		{loop_a, loop_a, 1, true, false},
		{loop_a, loop_b, -1, true, false},
		{loop_a, a_intersect_b, 1, true, false},
		{loop_a, a_union_b, 0, true, false},
		{loop_a, a_minus_b, 1, true, false},
		{loop_a, b_minus_a, 0, false, false},

		{loop_b, loop_a, -1, true, false},
		{loop_b, loop_b, 1, true, false},
		{loop_b, a_intersect_b, 1, true, false},
		{loop_b, a_union_b, 0, true, false},
		{loop_b, a_minus_b, 0, false, false},
		{loop_b, b_minus_a, 1, true, false},

		{a_intersect_b, loop_a, 0, true, false},
		{a_intersect_b, loop_b, 0, true, false},
		{a_intersect_b, a_intersect_b, 1, true, false},
		{a_intersect_b, a_union_b, 0, true, true},
		{a_intersect_b, a_minus_b, 0, false, false},
		{a_intersect_b, b_minus_a, 0, false, false},

		{a_union_b, loop_a, 1, true, false},
		{a_union_b, loop_b, 1, true, false},
		{a_union_b, a_intersect_b, 1, true, true},
		{a_union_b, a_union_b, 1, true, false},
		{a_union_b, a_minus_b, 1, true, false},
		{a_union_b, b_minus_a, 1, true, false},

		{a_minus_b, loop_a, 0, true, false},
		{a_minus_b, loop_b, 0, false, false},
		{a_minus_b, a_intersect_b, 0, false, false},
		{a_minus_b, a_union_b, 0, true, false},
		{a_minus_b, a_minus_b, 1, true, false},
		{a_minus_b, b_minus_a, 0, false, true},

		{b_minus_a, loop_a, 0, false, false},
		{b_minus_a, loop_b, 0, true, false},
		{b_minus_a, a_intersect_b, 0, false, false},
		{b_minus_a, a_union_b, 0, true, false},
		{b_minus_a, a_minus_b, 0, false, true},
		{b_minus_a, b_minus_a, 1, true, false},

		// Added in 2010/08: Make sure the relations are correct if the loop
		// crossing happens on two ends of a shared boundary segment.
		{loop_a, loop_c, -1, true, false},
		{loop_c, loop_a, -1, true, false},
		{loop_a, loop_d, 0, true, false},
		{loop_d, loop_a, 1, true, false},
	}
	for _, test := range tests {
		checkRelation(t, test.a, test.b, test.containsOrCrosses, test.intersects, test.nestable)
	}
}

func TestLoopRelations2(t *testing.T) {
	for iter := 0; iter < 1000; iter++ {
		begin := CellID(rand.Int63() | 1)
		if !begin.IsValid() {
			continue
		}
		begin = begin.Parent(rand.Intn(maxLevel))
		a_begin := begin.Advance(int64(skewed(6)))
		a_end := a_begin.Advance(int64(skewed(6)) + 1)
		b_begin := begin.Advance(int64(skewed(6)))
		b_end := b_begin.Advance(int64(skewed(6)) + 1)
		if !a_end.IsValid() || !b_end.IsValid() {
			continue
		}

		a := makeCellLoop(a_begin, a_end)
		b := makeCellLoop(b_begin, b_end)
		if a != nil && b != nil {
			contained := (a_begin <= b_begin && b_end <= a_end)
			intersects := (a_begin < b_end && b_begin < a_end)
			if got := a.ContainsLoop(b); got != contained {
				t.Errorf("%v.ContainsLoop(%v) = %v, want %v", a, b, got, contained)
			}
			if got := a.Intersects(b); got != intersects {
				t.Errorf("%v.Intersects(%v) = %v, want %v", a, b, got, intersects)
			}
		} else {
			log.Println("failed to create loop")
		}
	}
}

func makeCellLoop(begin, end CellID) *Loop {
	edges := map[Point]map[Point]bool{}
	for id := begin; id != end; id = id.Next() {
		cell := CellFromCellID(id)
		for k := 0; k < 4; k++ {
			a := cell.Vertex(k)
			b := cell.Vertex((k + 1) & 3)
			if _, ok := edges[b][a]; !ok {
				if _, ok := edges[a]; !ok {
					edges[a] = map[Point]bool{}
				}
				edges[a][b] = true
			}
			delete(edges[b], a)
			if len(edges[b]) == 0 {
				delete(edges, b)
			}
		}
	}
	vertices := []Point{}
	var p Point
	for p, _ = range edges {
		break
	}

	for len(edges) > 0 {
		if len(edges[p]) != 1 {
			log.Println(edges[p])
		}
		var next Point
		for next, _ = range edges[p] {
			break
		}
		vertices = append(vertices, p)
		delete(edges, p)
		p = next
	}
	return NewLoopFromPath(vertices)
}

func TestBoundaryNear(t *testing.T) {
	degree := s1.Degree.Radians()
	tests := []struct {
		astr     string
		bstr     string
		maxError float64
		want     bool
	}{
		{"0:0, 0:10, 5:5", "0:0.1, -0.1:9.9, 5:5.2", 0.5 * degree, true},
		{"0:0, 0:3, 0:7, 0:10, 3:7, 5:5", "0:0, 0:10, 2:8, 5:5, 4:4, 3:3, 1:1", 1e-3, true},
		{"0:0, 0:2, 2:2, 2:0", "0:0, 1.9999:1, 0:2, 2:2, 2:0", 0.5 * degree, false},
		{"0.1:0, 0.1:1, 0.1:2, 0.1:3, 0.1:4, 1:4, 2:4, 3:4, 2:4.1, 1:4.1, 2:4.2, 3:4.2, 4:4.2, 5:4.2",
			"0:0, 0:1, 0:2, 0:3, 0.1:2, 0.1:1, 0.2:2, 0.2:3, 0.2:4, 1:4.1, 2:4, 3:4, 4:4, 5:4",
			1.5 * degree, true},
		{"0.1:0, 0.1:1, 0.1:2, 0.1:3, 0.1:4, 1:4, 2:4, 3:4, 2:4.1, 1:4.1, 2:4.2, 3:4.2, 4:4.2, 5:4.2",
			"0:0, 0:1, 0:2, 0:3, 0.1:2, 0.1:1, 0.2:2, 0.2:3, 0.2:4, 1:4.1, 2:4, 3:4, 4:4, 5:4",
			0.5 * degree, false},
	}
	for _, test := range tests {
		a := makeloop(test.astr)
		b := makeloop(test.bstr)
		if got := a.BoundaryNear(b, test.maxError); got != test.want {
			t.Errorf("%v.BoundaryNear(%v, %v) = %v, want %v", a, b, test.maxError, got, test.want)
		}
		if got := b.BoundaryNear(a, test.maxError); got != test.want {
			t.Errorf("%v.BoundaryNear(%v, %v) = %v, want %v", b, a, test.maxError, got, test.want)
		}
	}
}

func TestIsValidDetectsInvalidLoops(t *testing.T) {
	tests := []struct {
		str string
	}{
		// Only two vertices.
		{"20:20, 21:21"},
		// There is a duplicate vertex.
		{"20:20, 21:21, 21:20, 20:20, 20:21"},
		// Some edges intersect.
		{"20:20, 21:21, 21:20.5, 21:20, 20:21"},
	}
	for _, test := range tests {
		loop := makeloop(test.str)
		if got := loop.IsValid(); got {
			t.Errorf("%v.IsValid() = %v, want false", loop, got)
		}
	}
}

func TestNormalizedCompatibleWithContains(t *testing.T) {
	checkNormalizeAndContains(t, line_triangle)
	checkNormalizeAndContains(t, skinny_chevron)
}

func checkNormalizeAndContains(t *testing.T, loop *Loop) {
	p := makepoint("40:40")
	flip := loop.Clone()
	flip.Invert()
	if loop.IsNormalized() == loop.Contains(p) {
		t.Errorf("%v.IsNormalized() == %v != %v.Contains(%v) == %v",
			loop, loop.IsNormalized(), loop, p, loop.Contains(p))
	}
	if flip.IsNormalized() == flip.Contains(p) {
		t.Errorf("%v.IsNormalized() == %v != %v.Contains(%v) == %v",
			flip, flip.IsNormalized(), flip, p, flip.Contains(p))
	}
	if loop.IsNormalized() == flip.IsNormalized() {
		t.Errorf("%v.IsNormalized() == %v.IsNormalized()", loop, flip)
	}
	flip.Normalize()
	if flip.Contains(p) {
		t.Errorf("%v.Contains(%v)", flip, p)
	}
}

func TestQuadtreeGetsComputedAutomatically(t *testing.T) {
	const numVerts = 200
	loopCenter := makepoint("42:107")
	loop := makeRegularLoop(loopCenter, numVerts, 7e-3)
	q := makepoint("5:5")
	index := NewLoopIndex(loop)
	it := NewEdgeIndexIterator(index)
	numCandidates := 0
	it.GetCandidates(q, q)
	for ; !it.Done(); it.Next() {
		numCandidates++
	}
	if numVerts != numCandidates {
		t.Errorf("numVerts (%v) != numCandidates (%v)", numVerts, numCandidates)
	}

	computed := false
	for i := 0; i < 500; i++ {
		numCandidates = 0
		for it.GetCandidates(q, q); !it.Done(); it.Next() {
			numCandidates++
		}
		if !computed && numCandidates == 0 {
			computed = true
		}
	}
	numCandidates = 0
	for it.GetCandidates(q, q); !it.Done(); it.Next() {
		numCandidates++
	}
	if 0 != numCandidates {
		t.Errorf("numCandidates (%v) should be 0", numCandidates)
	}
}

func makeRegularLoop(center Point, numVerts int, angleRadius float64) *Loop {
	m := FrameFromPoint(center)
	var vertices []Point
	radianStep := 2 * math.Pi / float64(numVerts)
	planarRadius := math.Tan(angleRadius)
	for vi := 0; vi < numVerts; vi++ {
		angle := float64(vi) * radianStep
		p := PointFromCoords(planarRadius*math.Cos(angle), planarRadius*math.Sin(angle), 1)
		vertices = append(vertices, PointFromFrame(m, p))
	}
	return NewLoopFromPath(vertices)
}
