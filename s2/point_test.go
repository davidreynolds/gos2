package s2

import (
	"github.com/davidreynolds/gos2/r3"
	"log"
	"math"
	"testing"
)

// float64Near reports if the two values are within the given epsilon.
func float64Near(x, y, ε float64) bool {
	return math.Abs(x-y) <= ε
}

func TestOriginPoint(t *testing.T) {
	if math.Abs(OriginPoint().Norm()-1) > 1e-16 {
		t.Errorf("Origin point norm = %v, want 1", OriginPoint().Norm())
	}
}

func TestPointCross(t *testing.T) {
	tests := []struct {
		p1x, p1y, p1z, p2x, p2y, p2z float64
	}{
		{1, 0, 0, 1, 0, 0},
		{1, 0, 0, 0, 1, 0},
		{0, 1, 0, 1, 0, 0},
		{1, 2, 3, -4, 5, -6},
	}
	for _, test := range tests {
		p1 := PointFromCoords(test.p1x, test.p1y, test.p1z)
		p2 := PointFromCoords(test.p2x, test.p2y, test.p2z)
		result := p1.PointCross(p2)
		if !float64Eq(result.Norm(), 1) {
			t.Errorf("|%v ⨯ %v| = %v, want 1", p1, p2, result.Norm())
		}
		if x := result.Dot(p1.Vector); !float64Eq(x, 0) {
			t.Errorf("|(%v ⨯ %v) · %v| = %v, want 0", p1, p2, p1, x)
		}
		if x := result.Dot(p2.Vector); !float64Eq(x, 0) {
			t.Errorf("|(%v ⨯ %v) · %v| = %v, want 0", p1, p2, p2, x)
		}
	}
}

func TestCCW(t *testing.T) {
	tests := []struct {
		p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z float64
		want                                        bool
	}{
		{1, 0, 0, 0, 1, 0, 0, 0, 1, true},
		{0, 1, 0, 0, 0, 1, 1, 0, 0, true},
		{0, 0, 1, 1, 0, 0, 0, 1, 0, true},
		{1, 1, 0, 0, 1, 1, 1, 0, 1, true},
		{-3, -1, 4, 2, -1, -3, 1, -2, 0, true},

		// All degenerate cases of CCW(). Let M_1, M_2, ... be the sequence of
		// submatrices whose determinant sign is tested by that function. Then the
		// i-th test below is a 3x3 matrix M (with rows A, B, C) such that:
		//
		//    det(M) = 0
		//    det(M_j) = 0 for j < i
		//    det(M_i) != 0
		//    A < B < C in lexicographic order.
		// det(M_1) = b0*c1 - b1*c0
		{-3, -1, 0, -2, 1, 0, 1, -2, 0, false},
		// det(M_2) = b2*c0 - b0*c2
		{-6, 3, 3, -4, 2, -1, -2, 1, 4, false},
		// det(M_3) = b1*c2 - b2*c1
		{0, -1, -1, 0, 1, -2, 0, 2, 1, false},
		// From this point onward, B or C must be zero, or B is proportional to C.
		// det(M_4) = c0*a1 - c1*a0
		{-1, 2, 7, 2, 1, -4, 4, 2, -8, false},
		// det(M_5) = c0
		{-4, -2, 7, 2, 1, -4, 4, 2, -8, false},
		// det(M_6) = -c1
		{0, -5, 7, 0, -4, 8, 0, -2, 4, false},
		// det(M_7) = c2*a0 - c0*a2
		{-5, -2, 7, 0, 0, -2, 0, 0, -1, false},
		// det(M_8) = c2
		{0, -2, 7, 0, 0, 1, 0, 0, 2, false},
	}

	for _, test := range tests {
		p1 := PointFromCoords(test.p1x, test.p1y, test.p1z)
		p2 := PointFromCoords(test.p2x, test.p2y, test.p2z)
		p3 := PointFromCoords(test.p3x, test.p3y, test.p3z)
		result := CCW(p1, p2, p3)
		if result != test.want {
			t.Errorf("CCW(%v, %v, %v) = %v, want %v", p1, p2, p3, result, test.want)
		}
		if test.want {
			// For these cases we can test the reversibility condition
			result = CCW(p3, p2, p1)
			if result == test.want {
				t.Errorf("CCW(%v, %v, %v) = %v, want %v", p3, p2, p1, result, !test.want)
			}
		}
	}
}

func TestPointDistance(t *testing.T) {
	tests := []struct {
		x1, y1, z1 float64
		x2, y2, z2 float64
		want       float64 // radians
	}{
		{1, 0, 0, 1, 0, 0, 0},
		{1, 0, 0, 0, 1, 0, math.Pi / 2},
		{1, 0, 0, 0, 1, 1, math.Pi / 2},
		{1, 0, 0, -1, 0, 0, math.Pi},
		{1, 2, 3, 2, 3, -1, 1.2055891055045298},
	}
	for _, test := range tests {
		p1 := PointFromCoords(test.x1, test.y1, test.z1)
		p2 := PointFromCoords(test.x2, test.y2, test.z2)
		if a := p1.Distance(p2).Radians(); !float64Eq(a, test.want) {
			t.Errorf("%v.Distance(%v) = %v, want %v", p1, p2, a, test.want)
		}
		if a := p2.Distance(p1).Radians(); !float64Eq(a, test.want) {
			t.Errorf("%v.Distance(%v) = %v, want %v", p2, p1, a, test.want)
		}
	}
}

func TestApproxEqual(t *testing.T) {
	epsilon := 1e-14
	tests := []struct {
		x1, y1, z1 float64
		x2, y2, z2 float64
		want       bool
	}{
		{1, 0, 0, 1, 0, 0, true},
		{1, 0, 0, 0, 1, 0, false},
		{1, 0, 0, 0, 1, 1, false},
		{1, 0, 0, -1, 0, 0, false},
		{1, 2, 3, 2, 3, -1, false},
		{1, 0, 0, 1 * (1 + epsilon), 0, 0, true},
		{1, 0, 0, 1 * (1 - epsilon), 0, 0, true},
		{1, 0, 0, 1 + epsilon, 0, 0, true},
		{1, 0, 0, 1 - epsilon, 0, 0, true},
		{1, 0, 0, 1, epsilon, 0, true},
		{1, 0, 0, 1, epsilon, epsilon, false},
		{1, epsilon, 0, 1, -epsilon, epsilon, false},
	}
	for _, test := range tests {
		p1 := PointFromCoords(test.x1, test.y1, test.z1)
		p2 := PointFromCoords(test.x2, test.y2, test.z2)
		if got := p1.ApproxEqual(p2); got != test.want {
			t.Errorf("%v.ApproxEqual(%v), got %v want %v", p1, p2, got, test.want)
		}
	}
}

func TestColinearPoints(t *testing.T) {
	// The following points happen to be *exactly colinear* along a line
	// that is approximately tangent to the surface of the unit sphere.
	// In fact, "c" is the exact midpoint of the line segment "ab". All of
	// these points are close enough to unit length to satisfy IsUnitLength.
	a := Point{r3.Vector{0.72571927877036835, 0.46058825605889098, 0.51106749730504852}}
	b := Point{r3.Vector{0.7257192746638208, 0.46058826573818168, 0.51106749441312738}}
	c := Point{r3.Vector{0.72571927671709457, 0.46058826089853633, 0.51106749585908795}}
	c_sub_a := c.Sub(a.Vector)
	b_sub_c := b.Sub(c.Vector)
	if c_sub_a != b_sub_c {
		t.Errorf("%v != %v", c_sub_a, b_sub_c)
	}

	if RobustCCW(a, b, c) == 0 {
		t.Errorf("%v == %v", RobustCCW(a, b, c), 0)
	}

	if RobustCCW(a, b, c) != RobustCCW(b, c, a) {
		t.Errorf("%v != %v", RobustCCW(a, b, c), RobustCCW(b, c, a))
	}

	if RobustCCW(a, b, c) != -RobustCCW(c, b, a) {
		t.Errorf("%v != %v", RobustCCW(a, b, c), -RobustCCW(c, b, a))
	}

	// The points "x1" and "x2" are exactly proportional, i.e. they both
	// lie on a common line through the origin. Both points are considered
	// to be normalized, and in fact they both satisfy (x == x.Normalize()).
	// Therefore the triangle (x1, x2, -x1) consists of three distinct
	// points that all lie on a common line through the origin.
	x1 := pc(0.99999999999999989, 1.4901161193847655e-08, 0)
	x2 := pc(1, 1.4901161193847656e-08, 0)
	neg_x1 := Point{x1.Neg()}
	if x1.Vector != x1.Normalize() {
		t.Errorf("%v != %v", x1, x1.Normalize())
	}
	if x2.Vector != x2.Normalize() {
		t.Errorf("%v != %v", x2, x2.Normalize())
	}

	if RobustCCW(x1, x2, neg_x1) == 0 {
		t.Errorf("%v == %v", RobustCCW(x1, x2, neg_x1), 0)
	}

	if RobustCCW(x1, x2, neg_x1) != RobustCCW(x2, neg_x1, x1) {
		t.Errorf("%v != %v", RobustCCW(x1, x2, neg_x1), RobustCCW(x2, neg_x1, x1))
	}

	if RobustCCW(x1, x2, neg_x1) != -RobustCCW(neg_x1, x2, x1) {
		t.Errorf("%v != %v", RobustCCW(x1, x2, neg_x1), -RobustCCW(neg_x1, x2, x1))
	}

	// Here are two more points that are distinct, exactly proportional,
	// and that satisfy (x == x.Normalize()).
	x3 := PointFromCoords(1, 1, 1)
	x4 := Point{x3.Mul(0.99999999999999989)}
	neg_x3 := Point{x3.Neg()}
	if x3.Vector != x3.Normalize() {
		t.Errorf("%v != %v", x3, x3.Normalize())
	}
	if x4.Vector != x4.Normalize() {
		t.Errorf("%v != %v", x4, x4.Normalize())
	}
	if RobustCCW(x3, x4, neg_x3) == 0 {
		t.Errorf("%v == %v", RobustCCW(x3, x4, neg_x3), 0)
	}
}

func TestSymbolicPerturbation(t *testing.T) {
	// The purpose of this test is simply to get code coverage of
	// SymbolicallyPerturbedCCW(). Let M_1, M_2, ... be the sequence of
	// submatrices whose determinant sign is tested by that function. Then
	// the i-th test below is a 3x3 matrix M (with rows A, B, C) such that:
	//
	//  det(M)    = 0
	//  det(M_j)  = 0 for j < i
	//  det(M_i) != 0
	//  A < B < C in lexicographic order.
	tests := []struct {
		p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z float64
		want                                        int
	}{
		// All degenerate cases of CCW().  Let M_1, M_2, ... be the sequence of
		// submatrices whose determinant sign is tested by that function.  Then the
		// i-th test below is a 3x3 matrix M (with rows A, B, C) such that:
		//
		//    det(M) = 0
		//    det(M_j) = 0 for j < i
		//    det(M_i) != 0
		//    A < B < C in lexicographic order.
		// det(M_1) = b0*c1 - b1*c0
		{-3, -1, 0, -2, 1, 0, 1, -2, 0, 1},
		// det(M_2) = b2*c0 - b0*c2
		{-6, 3, 3, -4, 2, -1, -2, 1, 4, 1},
		// det(M_3) = b1*c2 - b2*c1
		{0, -1, -1, 0, 1, -2, 0, 2, 1, 1},
		// From this point onward, B or C must be zero, or B is proportional to C.
		// det(M_4) = c0*a1 - c1*a0
		{-1, 2, 7, 2, 1, -4, 4, 2, -8, 1},
		// det(M_5) = c0
		{-4, -2, 7, 2, 1, -4, 4, 2, -8, 1},
		// det(M_6) = -c1
		{0, -5, 7, 0, -4, 8, 0, -2, 4, 1},
		// det(M_7) = c2*a0 - c0*a2
		{-5, -2, 7, 0, 0, -2, 0, 0, -1, 1},
		// det(M_8) = c2
		{0, -2, 7, 0, 0, 1, 0, 0, 2, 1},

		// From this point onward, C must be zero
		// det(M_9) = a0*b1 - a1*b0
		{-3, 1, 7, -1, -4, 1, 0, 0, 0, 1},
		// det(M_10) = -b0
		{-6, -4, 7, -3, -2, 1, 0, 0, 0, 1},
		// det(M_11) = b1
		{0, -4, 7, 0, -2, 1, 0, 0, 0, -1},
		// det(M_12) = a0
		{-1, -4, 5, 0, 0, -3, 0, 0, 0, -1},
		// det(M_13) = 1
		{0, -4, 5, 0, 0, -5, 0, 0, 0, 1},
	}
	for _, test := range tests {
		a := Point{r3.Vector{test.p1x, test.p1y, test.p1z}}
		b := Point{r3.Vector{test.p2x, test.p2y, test.p2z}}
		c := Point{r3.Vector{test.p3x, test.p3y, test.p3z}}
		checkSymbolicCCW(t, test.want, a, b, c)
	}
}

func checkSymbolicCCW(t *testing.T, want int, a, b, c Point) {
	if !a.LessThan(b.Vector) {
		t.Errorf("%v >= %v", a, b)
	}
	if !b.LessThan(c.Vector) {
		t.Errorf("%v >= %v", b, c)
	}
	if got := ExpensiveCCW(a, b, c); got != want {
		log.Printf("ExpensiveCCW(%v, %v, %v) = %v, want %v", a, b, c, got, want)
		t.Fatalf("ExpensiveCCW(%v, %v, %v) = %v, want %v", a, b, c, got, want)
	}
	if got := ExpensiveCCW(b, c, a); got != want {
		log.Printf("ExpensiveCCW(%v, %v, %v) = %v, want %v", b, c, a, got, want)
		t.Fatalf("ExpensiveCCW(%v, %v, %v) = %v, want %v", b, c, a, got, want)
	}
	if got := ExpensiveCCW(c, a, b); got != want {
		log.Printf("ExpensiveCCW(%v, %v, %v) = %v, want %v", c, a, b, got, want)
		t.Fatalf("ExpensiveCCW(%v, %v, %v) = %v, want %v", c, a, b, got, want)
	}
	if got := ExpensiveCCW(c, b, a); got != -want {
		log.Printf("ExpensiveCCW(%v, %v, %v) = %v, want %v", c, b, a, got, -want)
		t.Fatalf("ExpensiveCCW(%v, %v, %v) = %v, want %v", c, b, a, got, -want)
	}
	if got := ExpensiveCCW(b, a, c); got != -want {
		log.Printf("ExpensiveCCW(%v, %v, %v) = %v, want %v", b, a, c, got, -want)
		t.Fatalf("ExpensiveCCW(%v, %v, %v) = %v, want %v", b, a, c, got, -want)
	}
	if got := ExpensiveCCW(a, c, b); got != -want {
		log.Printf("ExpensiveCCW(%v, %v, %v) = %v, want %v", a, c, b, got, -want)
		t.Fatalf("ExpensiveCCW(%v, %v, %v) = %v, want %v", a, c, b, got, -want)
	}
}

var (
	pz   = PointFromCoords(0, 0, 1)
	p000 = PointFromCoords(1, 0, 0)
	p045 = PointFromCoords(1, 1, 0)
	p090 = PointFromCoords(0, 1, 0)
	p180 = PointFromCoords(-1, 0, 0)
	// Degenerate triangles.
	pr = PointFromCoords(0.257, -0.5723, 0.112)
	pq = PointFromCoords(-0.747, 0.401, 0.2235)

	// For testing the Girard area fall through case.
	g1 = PointFromCoords(1, 1, 1)
	g2 = Point{g1.Add(pr.Mul(1e-15)).Normalize()}
	g3 = Point{g1.Add(pq.Mul(1e-15)).Normalize()}
)

func TestPointArea(t *testing.T) {
	epsilon := 1e-10
	tests := []struct {
		a, b, c  Point
		want     float64
		nearness float64
	}{
		{p000, p090, pz, math.Pi / 2.0, 0},
		// This test case should give 0 as the epsilon, but either Go or C++'s value for Pi,
		// or the accuracy of the multiplications along the way, cause a difference ~15 decimal
		// places into the result, so it is not quite a difference of 0.
		{p045, pz, p180, 3.0 * math.Pi / 4.0, 1e-14},
		// Make sure that Area has good *relative* accuracy even for very small areas.
		{PointFromCoords(epsilon, 0, 1), PointFromCoords(0, epsilon, 1), pz, 0.5 * epsilon * epsilon, 1e-14},
		// Make sure that it can handle degenerate triangles.
		{pr, pr, pr, 0.0, 0},
		{pr, pq, pr, 0.0, 1e-15},
		{p000, p045, p090, 0.0, 0},
		// Try a very long and skinny triangle.
		{p000, PointFromCoords(1, 1, epsilon), p090, 5.8578643762690495119753e-11, 1e-9},
		// TODO(roberts):
		// C++ includes a 10,000 loop of perterbations to test out the Girard area
		// computation is less than some noise threshold.
		// Do we need that many? Will one or two suffice?
		{g1, g2, g3, 0.0, 1e-15},
	}
	for _, test := range tests {
		if got := PointArea(test.a, test.b, test.c); !float64Near(got, test.want, test.nearness) {
			t.Errorf("PointArea(%v, %v, %v), got %v want %v", test.a, test.b, test.c, got, test.want)
		}
	}
}

func TestPointAreaQuarterHemisphere(t *testing.T) {
	epsilon := 1e-14
	tests := []struct {
		a, b, c, d, e Point
		want          float64
	}{
		// Triangles with near-180 degree edges that sum to a quarter-sphere.
		{PointFromCoords(1, 0.1*epsilon, epsilon), p000, p045, p180, pz, math.Pi},
		// Four other triangles that sum to a quarter-sphere.
		{PointFromCoords(1, 1, epsilon), p000, p045, p180, pz, math.Pi},
		// TODO(roberts):
		// C++ Includes a loop of 100 perturbations on a hemisphere for more tests.
	}
	for _, test := range tests {
		area := PointArea(test.a, test.b, test.c) +
			PointArea(test.a, test.c, test.d) +
			PointArea(test.a, test.d, test.e) +
			PointArea(test.a, test.e, test.b)

		if !float64Eq(area, test.want) {
			t.Errorf("Adding up 4 quarter hemispheres with PointArea(), got %v want %v", area, test.want)
		}
	}
}

func BenchmarkPointArea(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PointArea(p000, p090, pz)
	}
}

func BenchmarkPointAreaGirardCase(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PointArea(g1, g2, g3)
	}
}
