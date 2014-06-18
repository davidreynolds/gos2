package s2

import (
	"github.com/davidreynolds/gos2/r3"
	"math"
	"testing"
)

func pc(x, y, z float64) Point { return Point{r3.Vector{x, y, z}} }

func TestEquality(t *testing.T) {
	tests := []struct {
		a    Point
		b    Point
		want bool
	}{
		{pc(0, 0, 0), pc(0, 0, 0), true},
		{pc(.5, 0, 0), pc(.5, 0, 0), true},
	}
	for _, test := range tests {
		got := test.a == test.b
		if got != test.want {
			t.Errorf("%v == %v = %v, want %v", test.a, test.b, got, test.want)
		}
	}
}

func TestWedges(t *testing.T) {
	tests := []struct {
		a0            Point
		ab1           Point
		a2            Point
		b0            Point
		b2            Point
		contains      bool
		intersects    bool
		wedgeRelation int
	}{
		// Intersection in one wedge.
		{pc(-1, 0, 10), pc(0, 0, 1), pc(1, 2, 10),
			pc(0, 1, 10), pc(1, -2, 10),
			false, true, WEDGE_PROPERLY_OVERLAPS},

		// Intersection in two wedges.
		{pc(-1, -1, 10), pc(0, 0, 1), pc(1, -1, 10),
			pc(1, 0, 10), pc(-1, 1, 10),
			false, true, WEDGE_PROPERLY_OVERLAPS},

		// Normal containment.
		{pc(-1, -1, 10), pc(0, 0, 1), pc(1, -1, 10),
			pc(-1, 0, 10), pc(1, 0, 10),
			true, true, WEDGE_PROPERLY_CONTAINS},

		// Containment with equality on one side.
		{pc(2, 1, 10), pc(0, 0, 1), pc(-1, -1, 10),
			pc(2, 1, 10), pc(1, -5, 10),
			true, true, WEDGE_PROPERLY_CONTAINS},

		// Containment with equality on the other side.
		{pc(2, 1, 10), pc(0, 0, 1), pc(-1, -1, 10),
			pc(1, -2, 10), pc(-1, -1, 10),
			true, true, WEDGE_PROPERLY_CONTAINS},

		// Containment with equality on both sides.
		{pc(-2, 3, 10), pc(0, 0, 1), pc(4, -5, 10),
			pc(-2, 3, 10), pc(4, -5, 10),
			true, true, WEDGE_EQUALS},

		// Disjoint with equality on one side.
		{pc(-2, 3, 10), pc(0, 0, 1), pc(4, -5, 10),
			pc(4, -5, 10), pc(-2, -3, 10),
			false, false, WEDGE_IS_DISJOINT},

		// Disjoint with equality on the other side.
		{pc(-2, 3, 10), pc(0, 0, 1), pc(0, 5, 10),
			pc(4, -5, 10), pc(-2, 3, 10),
			false, false, WEDGE_IS_DISJOINT},

		// Disjoint with equality on both sides.
		{pc(-2, 3, 10), pc(0, 0, 1), pc(4, -5, 10),
			pc(4, -5, 10), pc(-2, 3, 10),
			false, false, WEDGE_IS_DISJOINT},

		// B contains A with equality on one side.
		{pc(2, 1, 10), pc(0, 0, 1), pc(1, -5, 10),
			pc(2, 1, 10), pc(-1, -1, 10),
			false, true, WEDGE_IS_PROPERLY_CONTAINED},

		// B contains A with equality on the other side.
		{pc(2, 1, 10), pc(0, 0, 1), pc(1, -5, 10),
			pc(-2, 1, 10), pc(1, -5, 10),
			false, true, WEDGE_IS_PROPERLY_CONTAINED},
	}
	for _, test := range tests {
		test.a0 = Point{test.a0.Normalize()}
		test.ab1 = Point{test.ab1.Normalize()}
		test.a2 = Point{test.a2.Normalize()}
		test.b0 = Point{test.b0.Normalize()}
		test.b2 = Point{test.b2.Normalize()}
		got := WedgeContains(test.a0, test.ab1, test.a2, test.b0, test.b2)
		if got != test.contains {
			t.Errorf("WedgeContains(): got %v, want %v", got, test.contains)
		}
		got = WedgeIntersects(test.a0, test.ab1, test.a2, test.b0, test.b2)
		if got != test.intersects {
			t.Errorf("WedgeIntersects(): got %v, want %v", got, test.intersects)
		}
		relation := GetWedgeRelation(test.a0, test.ab1, test.a2, test.b0, test.b2)
		if relation != test.wedgeRelation {
			t.Errorf("GetWedgeRelation(): got %v, want %v", relation, test.wedgeRelation)
		}
	}
}

// Given a point X and an edge AB, check that the distance from X to AB is
// "distance_radians" and the closest point on AB is "expected_closest"
func TestDistanceToEdge(t *testing.T) {
	tests := []struct {
		x                Point
		a                Point
		b                Point
		distance_radians float64
		expected_closest Point
	}{
		{pc(1, 0, 0), pc(1, 0, 0), pc(0, 1, 0), 0, pc(1, 0, 0)},
		{pc(0, 1, 0), pc(1, 0, 0), pc(0, 1, 0), 0, pc(0, 1, 0)},
		{pc(1, 3, 0), pc(1, 0, 0), pc(0, 1, 0), 0, pc(1, 3, 0)},
		{pc(0, 0, 1), pc(1, 0, 0), pc(0, 1, 0), math.Pi / 2, pc(1, 0, 0)},
		{pc(0, 0, -1), pc(1, 0, 0), pc(0, 1, 0), math.Pi / 2, pc(1, 0, 0)},
		{pc(-1, -1, 0), pc(1, 0, 0), pc(0, 1, 0), 0.75 * math.Pi, pc(0, 0, 0)},

		{pc(0, 1, 0), pc(1, 0, 0), pc(1, 1, 0), math.Pi / 4, pc(1, 1, 0)},
		{pc(0, -1, 0), pc(1, 0, 0), pc(1, 1, 0), math.Pi / 2, pc(1, 0, 0)},

		{pc(0, -1, 0), pc(1, 0, 0), pc(-1, 1, 0), math.Pi / 2, pc(1, 0, 0)},
		{pc(-1, -1, 0), pc(1, 0, 0), pc(-1, 1, 0), math.Pi / 2, pc(-1, 1, 0)},

		{pc(1, 1, 1), pc(1, 0, 0), pc(0, 1, 0), math.Asin(math.Sqrt(1. / 3)), pc(1, 1, 0)},
		{pc(1, 1, -1), pc(1, 0, 0), pc(0, 1, 0), math.Asin(math.Sqrt(1. / 3)), pc(1, 1, 0)},

		{pc(-1, 0, 0), pc(1, 1, 0), pc(1, 1, 0), 0.75 * math.Pi, pc(1, 1, 0)},
		{pc(0, 0, -1), pc(1, 1, 0), pc(1, 1, 0), math.Pi / 2, pc(1, 1, 0)},
		{pc(-1, 0, 0), pc(1, 0, 0), pc(1, 0, 0), math.Pi, pc(1, 0, 0)},
	}
	for _, test := range tests {
		test.x = Point{test.x.Normalize()}
		test.a = Point{test.a.Normalize()}
		test.b = Point{test.b.Normalize()}
		test.expected_closest = Point{test.expected_closest.Normalize()}

		got := test.x.DistanceToEdge(test.a, test.b).Radians()
		if math.Abs(got-test.distance_radians) > 1e-14 {
			t.Errorf("%v.DistanceToEdge(%v, %v) = %v, want %v",
				test.x, test.a, test.b, got, test.distance_radians)
		}

		closest := test.x.ClosestPoint(test.a, test.b)
		if test.expected_closest == pc(0, 0, 0) {
			if closest != test.a && closest != test.b {
				t.Errorf("NOT: %v == %v || %v == %v", closest, test.a, closest, test.b)
			}
		} else {
			if !closest.ApproxEqual(test.expected_closest) {
				t.Errorf("%v != %v", closest, test.expected_closest)
			}
		}
	}
}

const kDegen int = -2

func testCrossing(t *testing.T, a, b, c, d Point, robust int, edgeOrVertex, simple bool) {
	a = Point{a.Normalize()}
	b = Point{b.Normalize()}
	c = Point{c.Normalize()}
	d = Point{d.Normalize()}
	CompareResult(t, RobustCrossing(a, b, c, d), robust)
	if simple {
		ok1 := robust > 0
		if ok1 != SimpleCrossing(a, b, c, d) {
			t.Errorf("SimpleCrossing() != %v", ok1)
		}
	}
	crosser := NewEdgeCrosser(&a, &b, &c)
	CompareResult(t, crosser.RobustCrossing(&d), robust)
	CompareResult(t, crosser.RobustCrossing(&c), robust)

	if edgeOrVertex != EdgeOrVertexCrossing(a, b, c, d) {
		t.Errorf("%v != EdgeOrVertexCrossing()", edgeOrVertex)
	}
	if edgeOrVertex != crosser.EdgeOrVertexCrossing(&d) {
		t.Errorf("%v != crosser.EdgeOrVertexCrossing(&d)", edgeOrVertex)
	}
	if edgeOrVertex != crosser.EdgeOrVertexCrossing(&c) {
		t.Errorf("%v != crosser.EdgeOrVertexCrossing(&c)", edgeOrVertex)
	}
}

func TestCrossings(t *testing.T) {
	tests := []struct {
		a, b, c, d   Point
		robust       int
		edgeOrVertex bool
		simple       bool
	}{
		// Two regular edges that cross.
		{pc(1, 2, 1), pc(1, -3, 0.5),
			pc(1, -0.5, -3), pc(0.1, 0.5, 3), 1, true, true},

		// Two regular edges that cross antipodal points.
		{pc(1, 2, 1), pc(1, -3, 0.5),
			pc(-1, 0.5, 3), pc(-0.1, -0.5, -3), -1, false, true},

		// Two edges on the same great circle.
		{pc(0, 0, -1), pc(0, 1, 0),
			pc(0, 1, 1), pc(0, 0, 1), -1, false, true},

		// Two edges that cross where one vertex is OriginPoint().
		{pc(1, 0, 0), OriginPoint(),
			pc(1, -0.1, 1), pc(1, 1, -0.1), 1, true, true},

		// Two edges that cross antipodal points where one vertex is S2::Origin().
		{pc(1, 0, 0), pc(0, 1, 0),
			pc(0, 0, -1), pc(-1, -1, 1), -1, false, true},

		// Two edges that share an endpoint.  The Ortho() direction is (-4,0,2),
		// and edge CD is further CCW around (2,3,4) than AB.
		{pc(2, 3, 4), pc(-1, 2, 5),
			pc(7, -2, 3), pc(2, 3, 4), 0, false, true},

		// Two edges that barely cross each other near the middle of one edge.  The
		// edge AB is approximately in the x=y plane, while CD is approximately
		// perpendicular to it and ends exactly at the x=y plane.
		{pc(1, 1, 1), pc(1, math.Nextafter(1, 0), -1),
			pc(11, -12, -1), pc(10, 10, 1), 1, true, false},

		// In this version, the edges are separated by a distance of about 1e-15.
		{pc(1, 1, 1), pc(1, math.Nextafter(1, 2), -1),
			pc(1, -1, 0), pc(1, 1, 0), -1, false, false},

		// Two edges that barely cross each other near the end of both edges.  This
		// example cannot be handled using regular double-precision arithmetic due
		// to floating-point underflow.
		{pc(0, 0, 1), pc(2, -1e-323, 1),
			pc(1, -1, 1), pc(1e-323, 0, 1), 1, true, false},

		// In this version, the edges are separated by a distance of about 1e-640.
		{pc(0, 0, 1), pc(2, 1e-323, 1),
			pc(1, -1, 1), pc(1e-323, 0, 1), -1, false, false},

		// Two edges that barely cross each other near the middle of one edge.
		// Computing the exact determinant of some of the triangles in this test
		// requires more than 2000 bits of precision.
		{pc(1, -1e-323, -1e-323), pc(1e-323, 1, 1e-323),
			pc(1, -1, 1e-323), pc(1, 1, 0),
			1, true, false},

		// In this version, the edges are separated by a distance of about 1e-640.
		{pc(1, 1e-323, -1e-323), pc(-1e-323, 1, 1e-323),
			pc(1, -1, 1e-323), pc(1, 1, 0),
			-1, false, false},
	}
	for _, test := range tests {
		a, b, c, d := test.a, test.b, test.c, test.d
		robust := test.robust
		edgeOrVertex := test.edgeOrVertex
		simple := test.simple

		r := false
		if robust == 0 {
			r = true
		}
		testCrossing(t, a, b, c, d, robust, edgeOrVertex, simple)
		testCrossing(t, b, a, c, d, robust, edgeOrVertex, simple)
		testCrossing(t, a, b, d, c, robust, edgeOrVertex, simple)
		testCrossing(t, b, a, d, c, robust, edgeOrVertex, simple)
		testCrossing(t, a, a, c, d, kDegen, false, false)
		testCrossing(t, a, b, c, c, kDegen, false, false)
		testCrossing(t, a, b, a, b, 0, true, false)
		testCrossing(t, c, d, a, b, robust, edgeOrVertex != r, simple)
	}
}

func CompareResult(t *testing.T, actual, expected int) {
	if expected == kDegen {
		if actual > 0 {
			t.Errorf("expected == kDegen: %d > 0", actual)
		}
	} else {
		if expected != actual {
			t.Errorf("expected = %v, actual = %v", expected, actual)
		}
	}
}

func TestCollinearEdgesThatDontTouch(t *testing.T) {
	for iter := 0; iter < 1000; iter++ {
		a := randomPoint()
		d := randomPoint()
		b := EdgeInterpolate(0.05, a, d)
		c := EdgeInterpolate(0.95, a, d)
		if got := RobustCrossing(a, b, c, d); got >= 0 {
			t.Errorf("RobustCrossing() == %v >= 0", got)
		}
		crosser := NewEdgeCrosser(&a, &b, &c)
		if got := crosser.RobustCrossing(&d); got >= 0 {
			t.Errorf("crosser.RobustCrossing(%v) == %v >= 0", d, got)
		}
		if got := crosser.RobustCrossing(&c); got >= 0 {
			t.Errorf("crosser.RobustCrossing(%v) == %v >= 0", c, got)
		}
	}
}
