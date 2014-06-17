package s2

import (
	"log"
	"math/rand"
	"strings"
	"testing"
)

func makepolygon(s string) *Polygon {
	loopStrs := strings.Split(s, ";")
	loops := []*Loop{}
	for _, str := range loopStrs {
		if str != "" {
			loop := makeloop(str)
			loop.Normalize()
			loops = append(loops, loop)
		}
	}
	return NewPolygonFromLoops(&loops)
}

var (
	// A set of nested loops around the point 0:0 (lat:lng).
	// Every vertex of kNear0 is a vertex of kNear1.
	kNearPoint = "0:0"
	kNear0     = "-1:0, 0:1, 1:0, 0:-1;"
	kNear1     = "-1:-1, -1:0, -1:1, 0:1, 1:1, 1:0, 1:-1, 0:-1;"
	kNear2     = "5:-2, -2:5, -1:-2;"
	kNear3     = "6:-3, -3:6, -2:-2;"
	kNearHemi  = "0:-90, -90:0, 0:90, 90:0;"

	// A set of nested loops around the point 0:180 (lat:lng).
	// Every vertex of kFar0 and kFar2 belongs to kFar1, and all
	// the loops except kFar2 are non-convex.
	kFar0    = "0:179, 1:180, 0:-179, 2:-180;"
	kFar1    = "0:179, -1:179, 1:180, -1:-179, 0:-179, 3:-178, 2:-180, 3:178;"
	kFar2    = "-1:-179, -1:179, 3:178, 3:-178;" // opposite direction
	kFar3    = "-3:-178, -2:179, -3:178, 4:177, 4:-177;"
	kFarHemi = "0:-90, 60:90, -60:90;"

	// A set of nested loops around the point -90:0 (lat:lng).
	kSouthPoint = "-89.9999:0.001"
	kSouth0a    = "-90:0, -89.99:0, -89.99:0.01;"
	kSouth0b    = "-90:0, -89.99:0.02, -89.99:0.03;"
	kSouth0c    = "-90:0, -89.99:0.04, -89.99:0.05;"
	kSouth1     = "-90:0, -89.9:-0.1, -89.9:0.1;"
	kSouth2     = "-90:0, -89.8:-0.2, -89.8:0.2;"
	kSouthHemi  = "0:-180, 0:60, 0:-60;"

	// Two different loops that surround all the Near and Far loops except
	// for the hemispheres.
	kNearFar1 = "-1:-9, -9:-9, -9:9, 9:9, 9:-9, 1:-9, 1:-175, 9:-175, 9:175, -9:175, -9:-175, -1:-175;"
	kNearFar2 = "-8:-4, 8:-4, 2:15, 2:170, 8:-175, -8:-175, -2:170, -2:15;"
	// Loops that result from intersection of other loops.
	kFarHSouthH = "0:-180, 0:90, -60:90, 0:-90;"

	// Rectangles that form a cross, with only shared vertices, no crossing edges.
	// Optional holes outside the intersecting region.
	kCross1          = "-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1;"
	kCross1SideHole  = "-1.5:0.5, -1.2:0.5, -1.2:-0.5, -1.5:-0.5;"
	kCross2          = "1:-2, 1:-1, 1:1, 1:2, -1:2, -1:1, -1:-1, -1:-2;"
	kCross2SideHole  = "0.5:-1.5, 0.5:-1.2, -0.5:-1.2, -0.5:-1.5;"
	kCrossCenterHole = "-0.5:0.5, 0.5:0.5, 0.5:-0.5, -0.5:-0.5;"

	// Two rectangles that intersect, but no edges cross and there's always
	// local containment (rather than crossing) at each shared vertex.
	// In this ugly ASCII art, 1 is A+B, 2 is B+C:
	//      +---+---+---+
	//      | A | B | C |
	//      +---+---+---+
	kOverlap1          = "0:1, 1:1, 2:1, 2:0, 1:0, 0:0;"
	kOverlap1SideHole  = "0.2:0.8, 0.8:0.8, 0.8:0.2, 0.2:0.2;"
	kOverlap2          = "1:1, 2:1, 3:1, 3:0, 2:0, 1:0;"
	kOverlap2SideHole  = "2.2:0.8, 2.8:0.8, 2.8:0.2, 2.2:0.2;"
	kOverlapCenterHole = "1.2:0.8, 1.8:0.8, 1.8:0.2, 1.2:0.2;"

	// Actual polygons
	near_0     = makepolygon(kNear0)
	near_10    = makepolygon(kNear0 + kNear1)
	near_30    = makepolygon(kNear3 + kNear0)
	near_32    = makepolygon(kNear2 + kNear3)
	near_3210  = makepolygon(kNear0 + kNear2 + kNear3 + kNear1)
	near_H3210 = makepolygon(kNear0 + kNear2 + kNear3 + kNearHemi + kNear1)

	far_10    = makepolygon(kFar0 + kFar1)
	far_21    = makepolygon(kFar2 + kFar1)
	far_321   = makepolygon(kFar2 + kFar3 + kFar1)
	far_H20   = makepolygon(kFar2 + kFarHemi + kFar0)
	far_H3210 = makepolygon(kFar2 + kFarHemi + kFar0 + kFar1 + kFar3)

	south_0ab          = makepolygon(kSouth0a + kSouth0b)
	south_2            = makepolygon(kSouth2)
	south_210b         = makepolygon(kSouth2 + kSouth0b + kSouth1)
	south_H21          = makepolygon(kSouth2 + kSouthHemi + kSouth1)
	south_H20abc       = makepolygon(kSouth2 + kSouth0b + kSouthHemi + kSouth0a + kSouth0c)
	nf1_n10_f2_s10abc  = makepolygon(kSouth0c + kFar2 + kNear1 + kNearFar1 + kNear0 + kSouth1 + kSouth0b + kSouth0a)
	nf2_n2_f210_s210ab = makepolygon(kFar2 + kSouth0a + kFar1 + kSouth1 + kFar0 + kSouth0b + kNearFar2 + kSouth2 + kNear2)
	f32_n0             = makepolygon(kFar2 + kNear0 + kFar3)
	n32_s0b            = makepolygon(kNear3 + kSouth0b + kNear2)

	cross1             = makepolygon(kCross1)
	cross1_side_hole   = makepolygon(kCross1 + kCross1SideHole)
	cross1_center_hole = makepolygon(kCross1 + kCrossCenterHole)
	cross2             = makepolygon(kCross2)
	cross2_side_hole   = makepolygon(kCross2 + kCross2SideHole)
	cross2_center_hole = makepolygon(kCross2 + kCrossCenterHole)

	overlap1             = makepolygon(kOverlap1)
	overlap1_side_hole   = makepolygon(kOverlap1 + kOverlap1SideHole)
	overlap1_center_hole = makepolygon(kOverlap1 + kOverlapCenterHole)
	overlap2             = makepolygon(kOverlap2)
	overlap2_side_hole   = makepolygon(kOverlap2 + kOverlap2SideHole)
	overlap2_center_hole = makepolygon(kOverlap2 + kOverlapCenterHole)

	far_H         = makepolygon(kFarHemi)
	south_H       = makepolygon(kSouthHemi)
	far_H_south_H = makepolygon(kFarHSouthH)
)

func TestRelations(t *testing.T) {
	tests := []struct {
		a          *Polygon
		b          *Polygon
		contains   int
		intersects bool
	}{
		{near_10, near_30, -1, true},
		{near_10, near_32, 0, false},
		{near_10, near_3210, -1, true},
		{near_10, near_H3210, 0, false},
		{near_30, near_32, 1, true},
		{near_30, near_3210, 1, true},
		{near_30, near_H3210, 0, true},
		{near_32, near_3210, -1, true},
		{near_32, near_H3210, 0, false},
		{near_3210, near_H3210, 0, false},

		{far_10, far_21, 0, false},
		{far_10, far_321, -1, true},
		{far_10, far_H20, 0, false},
		{far_10, far_H3210, 0, false},
		{far_21, far_321, 0, false},
		{far_21, far_H20, 0, false},
		{far_21, far_H3210, -1, true},
		{far_321, far_H20, 0, true},
		{far_321, far_H3210, 0, true},
		{far_H20, far_H3210, 0, true},

		{south_0ab, south_2, -1, true},
		{south_0ab, south_210b, 0, true},
		{south_0ab, south_H21, -1, true},
		{south_0ab, south_H20abc, -1, true},
		{south_2, south_210b, 1, true},
		{south_2, south_H21, 0, true},
		{south_2, south_H20abc, 0, true},
		{south_210b, south_H21, 0, true},
		{south_210b, south_H20abc, 0, true},
		{south_H21, south_H20abc, 1, true},

		{nf1_n10_f2_s10abc, nf2_n2_f210_s210ab, 0, true},
		{nf1_n10_f2_s10abc, near_32, 1, true},
		{nf1_n10_f2_s10abc, far_21, 0, false},
		{nf1_n10_f2_s10abc, south_0ab, 0, false},
		{nf1_n10_f2_s10abc, f32_n0, 1, true},

		{nf2_n2_f210_s210ab, near_10, 0, false},
		{nf2_n2_f210_s210ab, far_10, 1, true},
		{nf2_n2_f210_s210ab, south_210b, 1, true},
		{nf2_n2_f210_s210ab, south_0ab, 1, true},
		{nf2_n2_f210_s210ab, n32_s0b, 1, true},

		{cross1, cross2, 0, true},
		{cross1_side_hole, cross2, 0, true},
		{cross1_center_hole, cross2, 0, true},
		{cross1, cross2_side_hole, 0, true},
		{cross1, cross2_center_hole, 0, true},
		{cross1_side_hole, cross2_side_hole, 0, true},
		{cross1_center_hole, cross2_side_hole, 0, true},
		{cross1_side_hole, cross2_center_hole, 0, true},
		{cross1_center_hole, cross2_center_hole, 0, true},

		// These cases, when either polygon has a hole, test a different code path
		// from the other cases.
		{overlap1, overlap2, 0, true},
		{overlap1_side_hole, overlap2, 0, true},
		{overlap1_center_hole, overlap2, 0, true},
		{overlap1, overlap2_side_hole, 0, true},
		{overlap1, overlap2_center_hole, 0, true},
		{overlap1_side_hole, overlap2_side_hole, 0, true},
		{overlap1_center_hole, overlap2_side_hole, 0, true},
		{overlap1_side_hole, overlap2_center_hole, 0, true},
		{overlap1_center_hole, overlap2_center_hole, 0, true},
	}
	for _, test := range tests {
		want := test.contains > 0
		got := test.a.ContainsPolygon(test.b)
		if want != got {
			t.Errorf("%v != %v.ContainsPolygon(%v) == %v", want, test.a, test.b, got)
		}
		want = test.contains < 0
		got = test.b.ContainsPolygon(test.a)
		if want != got {
			t.Errorf("%v != %v.ContainsPolygon(%v) == %v", want, test.b, test.a, got)
		}
		if got := test.a.Intersects(test.b); got != test.intersects {
			log.Printf("%v.Intersect(%v) == %v, want %v", test.a, test.b, got, test.intersects)
			t.Errorf("%v.Intersects(%v) == %v, want %v", test.a, test.b, got, test.intersects)
		}
		// TODO: check contains and disjoint relations
	}
}

func TestSplitting(t *testing.T) {
	tests := []struct {
		poly *Polygon
	}{
		{near_H3210},
		{far_H3210},
		{south_0ab},
		{south_210b},
		{south_H20abc},
		{nf1_n10_f2_s10abc},
		{nf2_n2_f210_s210ab},
		{far_H},
		{south_H},
		{far_H_south_H},
	}
	for _, test := range tests {
		SplitAndAssemble(t, test.poly)
	}
}

func SplitAndAssemble(t *testing.T, polygon *Polygon) {
	builder := NewPolygonBuilder(DIRECTED_XOR())
	builder.AddPolygon(polygon)
	var expected Polygon
	if !builder.AssemblePolygon(&expected, nil) {
		// TODO: use Fatalf()
		t.Errorf("%v.AssemblePolygon() failed", builder)
	}
	for iter := 0; iter < 10; iter++ {
		coverer := NewRegionCoverer()
		diameter := 2 * polygon.CapBound().Radius().Radians()
		minLevel := MaxWidth.MinLevel(diameter)
		level := minLevel + rand.Intn(4)
		coverer.SetMinLevel(minLevel)
		coverer.SetMaxLevel(level)
		coverer.SetMaxCells(500)

		cells := coverer.Covering(*polygon)
		var covering CellUnion
		covering.Init(cells)
		CheckCompleteCovering(t, *polygon, covering, false, CellID(0))
	}
}

func CheckContains(t *testing.T, astr, bstr string) {
	a := makepolygon(astr)
	b := makepolygon(bstr)
	if !a.ContainsPolygon(b) {
		t.Errorf("!%v.Contains(%v)", a, b)
	}
}

func CheckContainsPoint(t *testing.T, astr, bstr string) {
	a := makepolygon(astr)
	if !a.ContainsPoint(makepoint(bstr)) {
		t.Errorf("%v.ContainsPoint(%v) == false", a, makepoint(bstr))
	}
}

func TestPolygonInit(t *testing.T) {
	CheckContains(t, kNear1, kNear0)
	CheckContains(t, kNear2, kNear1)
	CheckContains(t, kNear3, kNear2)
	CheckContains(t, kNearHemi, kNear3)
	CheckContains(t, kFar1, kFar0)
	CheckContains(t, kFar2, kFar1)
	CheckContains(t, kFar3, kFar2)
	CheckContains(t, kFarHemi, kFar3)
	CheckContains(t, kSouth1, kSouth0a)
	CheckContains(t, kSouth1, kSouth0b)
	CheckContains(t, kSouth1, kSouth0c)
	CheckContains(t, kSouthHemi, kSouth2)
	CheckContains(t, kNearFar1, kNear3)
	CheckContains(t, kNearFar1, kFar3)
	CheckContains(t, kNearFar2, kNear3)
	CheckContains(t, kNearFar2, kFar3)

	CheckContainsPoint(t, kNear0, kNearPoint)
	CheckContainsPoint(t, kNear1, kNearPoint)
	CheckContainsPoint(t, kNear2, kNearPoint)
	CheckContainsPoint(t, kNear3, kNearPoint)
	CheckContainsPoint(t, kNearHemi, kNearPoint)
	CheckContainsPoint(t, kSouth0a, kSouthPoint)
	CheckContainsPoint(t, kSouth1, kSouthPoint)
	CheckContainsPoint(t, kSouth2, kSouthPoint)
	CheckContainsPoint(t, kSouthHemi, kSouthPoint)
}

func TestOperations(t *testing.T) {
	//	var farSouth Polygon
	//	farSouth.InitToIntersection(far_H, south_H)
	//	CheckEqual(t, &farSouth, far_H_south_H, 1e-31)
}

func CheckEqual(t *testing.T, a, b *Polygon, maxError float64) {
	if a.IsNormalized() && b.IsNormalized() {
		//		t.Fatalf(a.BoundaryApproxEquals(b, maxError))
	}
}
