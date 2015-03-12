package s2

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"strings"
	"testing"

	"github.com/davidreynolds/gos2/r3"
	"github.com/davidreynolds/gos2/s1"
)

func makepolygon(s string) *Polygon {
	loops := makeloops(s)
	return NewPolygonFromLoops(loops)
}

func makeloops(s string) []*Loop {
	loopStrs := strings.Split(s, ";")
	loops := []*Loop{}
	for _, str := range loopStrs {
		if str != "" {
			loop := makeloop(str)
			loop.Normalize()
			loops = append(loops, loop)
		}
	}
	return loops
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

		if test.contains > 0 {
			checkContains2(t, test.a, test.b)
		} else if test.contains < 0 {
			checkContains2(t, test.b, test.a)
		}
		if !test.intersects {
			checkDisjoint(t, test.a, test.b)
		}
	}
}

var benchLoops = makeloops(kSouth0c + kFar2 + kNear1 + kNearFar1 + kNear0 + kSouth1 + kSouth0b + kSouth0a)

func BenchmarkInit(b *testing.B) {
	for i := 0; i < b.N; i++ {
		NewPolygonFromLoops(benchLoops)
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
		t.Fatalf("%v.AssemblePolygon() failed", builder)
	}
	for iter := 0; iter < 10; iter++ {
		coverer := NewRegionCoverer()
		diameter := 2 * polygon.CapBound().Radius().Radians()
		minLevel := MaxWidth.MinLevel(diameter)
		level := minLevel + rand.Intn(4)
		coverer.SetMinLevel(minLevel)
		coverer.SetMaxLevel(level)
		coverer.SetMaxCells(500)

		cells := coverer.Covering(polygon)
		var covering CellUnion
		covering.Init(cells)
		CheckCompleteCovering(t, polygon, covering, false, CellID(0))
		var pieces []*Polygon
		for i := 0; i < len(cells); i++ {
			cell := CellFromCellID(cells[i])
			window := NewPolygonFromCell(cell)
			var piece Polygon
			piece.InitToIntersection(polygon, window)
			pieces = append(pieces, &piece)
			log.Printf("\nPiece %d:\n  Window: %s\n  Piece: %s\n",
				i, polygonToString(window), polygonToString(&piece))
		}

		for len(pieces) > 1 {
			a := choosePiece(&pieces)
			b := choosePiece(&pieces)
			var c Polygon
			c.InitToUnion(a, b)
			pieces = append(pieces, &c)
			log.Printf("\nJoining piece a: %s\n  With piece b: %s\n  To get piece c: %s\n",
				polygonToString(a), polygonToString(b), polygonToString(&c))
		}
		result := pieces[0]
		pieces = pieces[:len(pieces)-1]
		if got := expected.BoundaryNear(result, 1e-15); !got {
			t.Errorf("\nActual:\n%s\nExpected:\n%s\n",
				polygonToString(result), polygonToString(&expected))
		}
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

type testCase struct {
	a         string
	b         string
	a_and_b   string
	a_or_b    string
	a_minus_b string
}

var testCases = []testCase{
	// Two triangles that share an edge.
	{
		"4:2, 3:1, 3:3;",
		"3:1, 2:2, 3:3;",
		"",
		"4:2, 3:1, 2:2, 3:3;",
		"4:2, 3:1, 3:3;",
	},

	// Two vertical bars and a horizontal bar connecting them.
	{
		"0:0, 0:2, 3:2, 3:0;   0:3, 0:5, 3:5, 3:3;",
		"1:1, 1:4, 2:4, 2:1;",
		"1:1, 1:2, 2:2, 2:1;   1:3, 1:4, 2:4, 2:3;",
		"0:0, 0:2, 1:2, 1:3, 0:3, 0:5, 3:5, 3:3, 2:3, 2:2, 3:2, 3:0;",
		"0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0; 0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3;",
	},

	// Two vertical bars and two horizontal bars centered around S2::Origin().
	{
		"1:88, 1:93, 2:93, 2:88;   -1:88, -1:93, 0:93, 0:88;",
		"-2:89, -2:90, 3:90, 3:89;   -2:91, -2:92, 3:92, 3:91;",
		"1:89, 1:90, 2:90, 2:89;   1:91, 1:92, 2:92, 2:91; -1:89, -1:90, 0:90, 0:89;   -1:91, -1:92, 0:92, 0:91;",
		"-1:88, -1:89, -2:89, -2:90, -1:90, -1:91, -2:91, -2:92, -1:92, -1:93, 0:93, 0:92, 1:92, 1:93, 2:93, 2:92, 3:92, 3:91, 2:91, 2:90, 3:90, 3:89, 2:89, 2:88, 1:88, 1:89, 0:89, 0:88; 0:90, 0:91, 1:91, 1:90;",

		"1:88, 1:89, 2:89, 2:88; 1:90, 1:91, 2:91, 2:90; 1:92, 1:93, 2:93, 2:92; -1:88, -1:89, 0:89, 0:88; -1:90, -1:91, 0:91, 0:90; -1:92, -1:93, 0:93, 0:92;",
	},

	// Two interlocking square doughnuts centered around -S2::Origin().
	{
		"-1:-93, -1:-89, 3:-89, 3:-93;   0:-92, 0:-90, 2:-90, 2:-92;",
		"-3:-91, -3:-87, 1:-87, 1:-91;   -2:-90, -2:-88, 0:-88, 0:-90;",
		"-1:-91, -1:-90, 0:-90, 0:-91;   0:-90, 0:-89, 1:-89, 1:-90;",
		"-1:-93, -1:-91, -3:-91, -3:-87, 1:-87, 1:-89, 3:-89, 3:-93; 0:-92, 0:-91, 1:-91, 1:-90, 2:-90, 2:-92; -2:-90, -2:-88, 0:-88, 0:-89, -1:-89, -1:-90;",

		"-1:-93, -1:-91, 0:-91, 0:-92, 2:-92, 2:-90, 1:-90, 1:-89, 3:-89, 3:-93; -1:-90, -1:-89, 0:-89, 0:-90;",
	},

	// An incredibly thin triangle intersecting a square, such that the two
	// intersection points of the triangle with the square are identical.
	// This results in a degenerate loop that needs to be handled correctly.
	{
		"10:44, 10:46, 12:46, 12:44;",
		"11:45, 89:45.00000000000001, 90:45;",
		"", // Empty intersection!
		// Original square with extra vertex, and triangle disappears (due to
		// default vertex_merge_radius of S2EdgeUtil::kIntersectionTolerance).
		"10:44, 10:46, 12:46, 12:45, 12:44;",
		"10:44, 10:46, 12:46, 12:45, 12:44;",
	},
}

func TestOperations(t *testing.T) {
	var farSouth Polygon
	farSouth.InitToIntersection(far_H, south_H)
	CheckEqual(t, &farSouth, far_H_south_H, 1e-31)

	for _, test := range testCases {
		a := makepolygon(test.a)
		b := makepolygon(test.b)
		expected_a_and_b := makepolygon(test.a_and_b)
		expected_a_or_b := makepolygon(test.a_or_b)
		expected_a_minus_b := makepolygon(test.a_minus_b)
		maxErr := 1e-4
		var a_and_b, a_or_b, a_minus_b Polygon
		a_and_b.InitToIntersection(a, b)
		CheckEqual(t, &a_and_b, expected_a_and_b, maxErr)
		a_or_b.InitToUnion(a, b)
		CheckEqual(t, &a_or_b, expected_a_or_b, maxErr)
		a_minus_b.InitToDifference(a, b)
		CheckEqual(t, &a_minus_b, expected_a_minus_b, maxErr)
	}
}

func TestPolylineIntersection(t *testing.T) {
	for v := 0; v < 3; v++ {
		polylineIntersectionSharedEdgeTest(t, cross1, v, 1)
		polylineIntersectionSharedEdgeTest(t, cross1, v+1, -1)
		polylineIntersectionSharedEdgeTest(t, cross1_side_hole, v, 1)
		polylineIntersectionSharedEdgeTest(t, cross1_side_hole, v+1, -1)
	}

	const maxErr = 1e-4
	for _, test := range testCases {
		a := makepolygon(test.a)
		b := makepolygon(test.b)
		expected_a_and_b := makepolygon(test.a_and_b)

		var points []Point
		var polylines []*Polyline
		for ab := 0; ab < 2; ab++ {
			var tmp0, tmp1 *Polygon
			if ab != 0 {
				tmp0 = a
				tmp1 = b
			} else {
				tmp0 = b
				tmp1 = a
			}
			for l := 0; l < tmp0.NumLoops(); l++ {
				points = []Point{}
				if tmp0.Loop(l).IsHole() {
					for v := tmp0.Loop(l).NumVertices(); v >= 0; v-- {
						points = append(points, *tmp0.Loop(l).Vertex(v))
					}
				} else {
					for v := 0; v <= tmp0.Loop(l).NumVertices(); v++ {
						points = append(points, *tmp0.Loop(l).Vertex(v))
					}
				}
				polyline := PolylineFromPoints(points)
				var lines []*Polyline
				tmp1.IntersectWithPolyline(polyline, &lines)
				polylines = append(polylines, lines...)
			}
		}

		builder := NewPolygonBuilder(DIRECTED_XOR())
		for i := 0; i < len(polylines); i++ {
			for j := 0; j < polylines[i].NumVertices()-1; j++ {
				builder.AddEdge(polylines[i].Vertex(j), polylines[i].Vertex(j+1))
			}
		}
		var a_and_b Polygon
		if !builder.AssemblePolygon(&a_and_b, nil) {
			t.Fatalf("AssemblePolygon() failed")
		}
		CheckEqual(t, &a_and_b, expected_a_and_b, maxErr)
	}
}

func polylineIntersectionSharedEdgeTest(t *testing.T, p *Polygon, startVertex, dir int) {
	var points []Point
	points = append(points, *p.Loop(0).Vertex(startVertex))
	points = append(points, *p.Loop(0).Vertex(startVertex + dir))
	polyline := PolylineFromPoints(points)
	var polylines []*Polyline
	if dir < 0 {
		p.IntersectWithPolyline(polyline, &polylines)
		if len(polylines) != 0 {
			t.Errorf("len(%v) != 0", polylines)
		}
		polylines = []*Polyline{}
		p.SubtractFromPolyline(polyline, &polylines)
		if len(polylines) != 1 {
			t.Fatalf("len(%v) != 1", polylines)
		}
		if got := polylines[0].NumVertices(); got != 2 {
			t.Fatalf("%v.NumVertices() == %v, want 2", polylines[0], got)
		}
		if points[0] != polylines[0].Vertex(0) {
			t.Errorf("%v != %v", points[0], polylines[0].Vertex(0))
		}
		if points[1] != polylines[0].Vertex(1) {
			t.Errorf("%v != %v", points[1], polylines[0].Vertex(1))
		}
	} else {
		p.IntersectWithPolyline(polyline, &polylines)
		if len(polylines) != 1 {
			t.Fatalf("len(%v) != 1", polylines)
		}
		if got := polylines[0].NumVertices(); got != 2 {
			t.Fatalf("%v.NumVertices() == %v, want 2", polylines[0], got)
		}
		if points[0] != polylines[0].Vertex(0) {
			t.Errorf("%v != %v", points[0], polylines[0].Vertex(0))
		}
		if points[1] != polylines[0].Vertex(1) {
			t.Errorf("%v != %v", points[1], polylines[0].Vertex(1))
		}
		polylines = []*Polyline{}
		p.SubtractFromPolyline(polyline, &polylines)
		if len(polylines) != 0 {
			t.Errorf("len(%v) != 0", polylines)
		}
	}
}

func TestCellConstructorAndContains(t *testing.T) {
	ll := LatLng{(40565459 * s1.E6), (-74645276 * s1.E6)}
	cell := CellFromLatLng(ll)
	cell_as_poly := NewPolygonFromCell(cell)
	var empty Polygon
	var poly_copy Polygon
	poly_copy.InitToUnion(cell_as_poly, &empty)
	if !poly_copy.ContainsPolygon(cell_as_poly) {
		t.Errorf("%v.ContainsPolygon(%v) failed", poly_copy, cell_as_poly)
	}
	if !poly_copy.ContainsCell(cell) {
		t.Errorf("%v.ContainsCell(%v) failed", poly_copy, cell)
	}
}

func CheckEqual(t *testing.T, a, b *Polygon, maxError float64) {
	if a.IsNormalized() && b.IsNormalized() {
		if !a.BoundaryApproxEquals(b, maxError) {
			t.Fatalf("\nBoundaryApproxEquals() failed for:\n  A: %s\n\n  B: %s\n\n", polygonToString(a), polygonToString(b))
		}
	} else {
		var a2, b2 Polygon
		builder := NewPolygonBuilder(DIRECTED_XOR())
		builder.AddPolygon(a)
		if !builder.AssemblePolygon(&a2, nil) {
			t.Fatalf("Expected true")
		}
		builder.AddPolygon(b)
		if !builder.AssemblePolygon(&b2, nil) {
			t.Fatalf("Expected true")
		}
		if !a2.BoundaryApproxEquals(&b2, maxError) {
			t.Fatalf("Expected approx equals")
		}
	}
}

func checkDisjoint(t *testing.T, a, b *Polygon) {
	builder := NewPolygonBuilder(DIRECTED_XOR())
	builder.AddPolygon(a)
	builder.AddPolygon(b)
	var ab Polygon
	if !builder.AssemblePolygon(&ab, nil) {
		log.Fatalf("%v.AssemblePolygon(%v, nil) failed", builder, ab)
	}

	var c, d, e, f Polygon
	c.InitToUnion(a, b)
	CheckEqual(t, &c, &ab, 1e-31)
	checkUnion(t, a, b)

	d.InitToIntersection(a, b)
	if got := d.NumLoops(); got != 0 {
		t.Errorf("%v.NumLoop() == %v, expected 0", d, got)
	}

	e.InitToDifference(a, b)
	CheckEqual(t, &e, a, 1e-31)

	f.InitToDifference(b, a)
	CheckEqual(t, &f, b, 1e-31)
}

func checkUnion(t *testing.T, a, b *Polygon) {
	var cUnion Polygon
	cUnion.InitToUnion(a, b)
	cDestructiveUnion := DestructiveUnion(&[]*Polygon{a, b})
	CheckEqual(t, &cUnion, cDestructiveUnion, 1e-31)
}

func checkContains2(t *testing.T, a, b *Polygon) {
	var c, d, e Polygon
	c.InitToUnion(a, b)
	CheckEqual(t, &c, a, 1e-31)
	checkUnion(t, &c, a)

	d.InitToIntersection(a, b)
	CheckEqual(t, &d, b, 1e-31)

	e.InitToDifference(b, a)
	if got := e.NumLoops(); got != 0 {
		t.Errorf("%v.NumLoops() == %v, want 0", e, got)
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

func appendVertex(v Point, s *string) {
	ll := LatLngFromPoint(v)
	*s += fmt.Sprintf("%.17f:%.17f", ll.Lat.Degrees(), ll.Lng.Degrees())
}

func appendLoopVertices(loop *Loop, s *string) {
	for i := 0; i < loop.NumVertices(); i++ {
		if i > 0 {
			*s += ", "
		}
		appendVertex(*loop.Vertex(i), s)
	}
}

func polygonToString(poly *Polygon) string {
	var s string
	for i := 0; i < poly.NumLoops(); i++ {
		loop := poly.Loop(i)
		if i > 0 {
			s += ";\n"
		}
		appendLoopVertices(loop, &s)
	}
	return s
}

func choosePiece(pieces *[]*Polygon) *Polygon {
	i := rand.Intn(len(*pieces))
	res := (*pieces)[i]
	*pieces = append((*pieces)[:i], (*pieces)[i+1:]...)
	return res
}

func concentricLoops(center Point, numLoops, numVertsPerLoop int, poly *Polygon) {
	m := FrameFromPoint(center)
	var loops []*Loop
	for i := 0; i < numLoops; i++ {
		var vertices []Point
		radius := 0.005 * float64((i+1)/numLoops)
		radianStep := 2 * math.Pi / float64(numVertsPerLoop)
		for j := 0; j < numVertsPerLoop; j++ {
			angle := float64(j) * radianStep
			p := PointFromCoords(radius*math.Cos(angle), radius*math.Sin(angle), 1)
			vertices = append(vertices, PointFromFrame(m, p))
		}
		loops = append(loops, NewLoopFromPath(vertices))
	}
	poly.Init(loops)
}

func unionOfPolygons(b *testing.B, numVertsPerLoop int, offset float64) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		var p1, p2 Polygon
		center := randomPoint()
		concentricLoops(center, 10, numVertsPerLoop, &p1)
		center2 := PointFromCoords(offset, offset, offset)
		concentricLoops(Point{center.Add(r3.Vector{center2.X, center2.Y, center2.Z}).Normalize()}, 10, numVertsPerLoop, &p2)
		b.StartTimer()
		var p3 Polygon
		p3.InitToUnion(&p1, &p2)
	}
}

func benchmarkDeepPolygonUnion(b *testing.B, numVertsPerLoop int) {
	unionOfPolygons(b, numVertsPerLoop, 0.000001)
}

func benchmarkShallowPolygonUnion(b *testing.B, numVertsPerLoop int) {
	unionOfPolygons(b, numVertsPerLoop, 0.004)
}

func benchmarkDisjointPolygonUnion(b *testing.B, numVertsPerLoop int) {
	unionOfPolygons(b, numVertsPerLoop, 0.3)
}

func BenchmarkDeepPolygonUnion8(b *testing.B)    { benchmarkDeepPolygonUnion(b, 8) }
func BenchmarkDeepPolygonUnion64(b *testing.B)   { benchmarkDeepPolygonUnion(b, 64) }
func BenchmarkDeepPolygonUnion128(b *testing.B)  { benchmarkDeepPolygonUnion(b, 128) }
func BenchmarkDeepPolygonUnion256(b *testing.B)  { benchmarkDeepPolygonUnion(b, 256) }
func BenchmarkDeepPolygonUnion512(b *testing.B)  { benchmarkDeepPolygonUnion(b, 512) }
func BenchmarkDeepPolygonUnion1024(b *testing.B) { benchmarkDeepPolygonUnion(b, 1024) }
func BenchmarkDeepPolygonUnion4096(b *testing.B) { benchmarkDeepPolygonUnion(b, 4096) }
func BenchmarkDeepPolygonUnion8192(b *testing.B) { benchmarkDeepPolygonUnion(b, 8192) }

func BenchmarkShallowPolygonUnion8(b *testing.B)    { benchmarkShallowPolygonUnion(b, 8) }
func BenchmarkShallowPolygonUnion64(b *testing.B)   { benchmarkShallowPolygonUnion(b, 64) }
func BenchmarkShallowPolygonUnion128(b *testing.B)  { benchmarkShallowPolygonUnion(b, 128) }
func BenchmarkShallowPolygonUnion256(b *testing.B)  { benchmarkShallowPolygonUnion(b, 256) }
func BenchmarkShallowPolygonUnion512(b *testing.B)  { benchmarkShallowPolygonUnion(b, 512) }
func BenchmarkShallowPolygonUnion1024(b *testing.B) { benchmarkShallowPolygonUnion(b, 1024) }
func BenchmarkShallowPolygonUnion4096(b *testing.B) { benchmarkShallowPolygonUnion(b, 4096) }
func BenchmarkShallowPolygonUnion8192(b *testing.B) { benchmarkShallowPolygonUnion(b, 8192) }

func BenchmarkDisjointPolygonUnion8(b *testing.B)    { benchmarkDisjointPolygonUnion(b, 8) }
func BenchmarkDisjointPolygonUnion64(b *testing.B)   { benchmarkDisjointPolygonUnion(b, 64) }
func BenchmarkDisjointPolygonUnion128(b *testing.B)  { benchmarkDisjointPolygonUnion(b, 128) }
func BenchmarkDisjointPolygonUnion256(b *testing.B)  { benchmarkDisjointPolygonUnion(b, 256) }
func BenchmarkDisjointPolygonUnion512(b *testing.B)  { benchmarkDisjointPolygonUnion(b, 512) }
func BenchmarkDisjointPolygonUnion1024(b *testing.B) { benchmarkDisjointPolygonUnion(b, 1024) }
func BenchmarkDisjointPolygonUnion4096(b *testing.B) { benchmarkDisjointPolygonUnion(b, 4096) }
func BenchmarkDisjointPolygonUnion8192(b *testing.B) { benchmarkDisjointPolygonUnion(b, 8192) }
