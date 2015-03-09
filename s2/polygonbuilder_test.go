package s2

import (
	"fmt"
	"math"
	"math/rand"
	"testing"

	"github.com/davidreynolds/gos2/r3"
	"github.com/davidreynolds/gos2/s1"
)

func TestAddEdge(t *testing.T) {
	pb := NewPolygonBuilder(DIRECTED_XOR())
	p1 := PointFromLatLng(LatLngFromDegrees(47, -122))
	p2 := PointFromLatLng(LatLngFromDegrees(48, -122))
	if ok := pb.AddEdge(p1, p2); !ok {
		t.Errorf("Failed to add edge")
	}
}

func TestExistingDirectedEdge(t *testing.T) {
	pb := NewPolygonBuilder(DIRECTED_XOR())
	p1 := PointFromLatLng(LatLngFromDegrees(47, -122))
	p2 := PointFromLatLng(LatLngFromDegrees(48, -122))
	pb.AddEdge(p1, p2)
	pb.AddEdge(p2, p1) // Old edge is removed and new one isn't added
	if pb.HasEdge(p1, p2) || pb.HasEdge(p2, p1) {
		t.Errorf("Undirected edge found in directed graph")
	}
}

func TestHasEdge(t *testing.T) {
	pb := NewPolygonBuilder(DIRECTED_XOR())
	p1 := PointFromLatLng(LatLngFromDegrees(47, -122))
	p2 := PointFromLatLng(LatLngFromDegrees(48, -122))
	pb.AddEdge(p1, p2)
	if !pb.HasEdge(p1, p2) {
		t.Errorf("Edge not found")
	}
}

func TestSingleLoop(t *testing.T) {
	rect := []float64{
		-172.08984375,
		73.1758971742261,
		-21.4453125,
		-25.641526373065755,
	}
	vertices := []Point{
		PointFromLatLng(LatLngFromDegrees(rect[1], rect[0])),
		PointFromLatLng(LatLngFromDegrees(rect[3], rect[0])),
		PointFromLatLng(LatLngFromDegrees(rect[3], rect[2])),
		PointFromLatLng(LatLngFromDegrees(rect[1], rect[2])),
	}
	pb := NewPolygonBuilder(DIRECTED_XOR())
	for i := 0; i < 4; i++ {
		pb.AddEdge(vertices[i], vertices[(i+1)%len(vertices)])
	}

	loops := []*Loop{}
	edges := []Edge{}
	if !pb.AssembleLoops(&loops, &edges) {
		t.Errorf("unused edges: %v (loops: %v)", edges, loops)
	}

	if len(loops) != 1 {
		t.Errorf("len(%v) == %v, want 1", loops, len(loops))
	}
}

func BenchmarkBuilder(b *testing.B) {
	b.StopTimer()
	test := TestCase{-1, 0, true, 0.0, 0.9, 30.0, // Directed edges required for unique result.
		[]Chain{
			Chain{"0:0, 0:2, 0:4, 0:6, 1:5, 2:4, 3:3, 2:2, 1:1", true},
			Chain{"0:2, 1:1, 1:3", true},
			Chain{"0:4, 1:3, 1:5", true},
			Chain{"1:3, 2:2, 2:4", true},
			Chain{"0:0, -1:1", false},
			Chain{"3:3, 5:5", false},
		}, []string{
			"0:0, 0:2, 1:1",
			"0:2, 0:4, 1:3",
			"0:4, 0:6, 1:5",
			"1:1, 1:3, 2:2",
			"1:3, 1:5, 2:4",
			"2:2, 2:4, 3:3",
		}, 2}
	builder := NewPolygonBuilder(DIRECTED_XOR())
	for _, chain := range test.chainsIn {
		var vertices []Point
		line := makePolyline(chain.str)
		for _, v := range line.vertices {
			vertices = append(vertices, Point{v.Vector})
		}
		for i := 1; i < len(vertices); i++ {
			builder.AddEdge(vertices[i-1], vertices[i])
		}
	}
	b.StartTimer()
	var poly Polygon
	for i := 0; i < b.N; i++ {
		builder.AssemblePolygon(&poly, nil)
	}
}

type Chain struct {
	str    string
	closed bool
}

type TestCase struct {
	undirectedEdges    int      // +1 = undirected, -1 = directed, 0 = either
	xorEdges           int      // +1 = XOR, -1 = don't XOR, 0 = either
	canSplit           bool     // Can edges be split for this test case?
	minMerge, maxMerge float64  // Min and max vertex merge radius for this test case in degrees.
	minVertexAngle     float64  // Min angle in degrees between any two edges *after* vertex merging.
	chainsIn           []Chain  // Each test case consists of a set of input loops and polylines.
	loopsOut           []string // The expected set of output loops, directed appropriately.
	numUnusedEdges     int      // The expected number of unused edges.
}

func TestAssembleLoops(t *testing.T) {
	rand.Seed(4)
	tests := []TestCase{
		// 0: No loops.
		{0, 0, true, 0.0, 10.0, 90.0, []Chain{Chain{"", false}}, []string{""}, 0},

		// 1: One loop with some extra edges.
		{0, 0, true, 0.0, 4.0, 15.0,
			[]Chain{
				Chain{"0:0, 0:10, 10:5", true},
				Chain{"0:0, 5:5", false},
				Chain{"10:5, 20:7, 30:10, 40:15, 50:3, 60:-20", false},
			},
			[]string{
				"0:0, 0:10, 10:5",
			}, 6},

		// 2: One loop that has an edge removed by XORing, plus lots of extra edges.
		{0, 1, true, 0.0, 1.0, 45.0, // XOR
			[]Chain{
				Chain{"0:0, 0:10, 5:15, 10:10, 10:0", true},
				Chain{"10:10, 12:12, 14:14, 16:16, 18:18", false},
				Chain{"14:14, 14:16, 14:18, 14:20", false},
				Chain{"14:18, 16:20, 18:22", false},
				Chain{"18:12, 16:12, 14:12, 12:12", false},
				Chain{"20:18, 18:16, 16:14, 14:12", false},
				Chain{"20:14, 18:14, 16:14", false},
				Chain{"5:15, 0:10", false},
			}, []string{
				"",
			}, 21},

		// 3: Three loops (two shells and one hole) that combine into one.
		{0, 1, true, 0.0, 4.0, 90.0, // XOR
			[]Chain{
				Chain{"0:0, 0:10, 5:10, 10:10, 10:5, 10:0", true},
				Chain{"0:10, 0:15, 5:15, 5:10", true},
				Chain{"10:10, 5:10, 5:5, 10:5", true},
			}, []string{
				"0:0, 0:10, 0:15, 5:15, 5:10, 5:5, 10:5, 10:0",
			}, 0},

		// 4: A big CCW triangle contained 3 CW triangular holes. The
		// whole thing looks like a pyramid of nine small triangles
		// (with two extra edges).
		{-1, 0, true, 0.0, 0.9, 30.0, // Directed edges required for unique result.
			[]Chain{
				Chain{"0:0, 0:2, 0:4, 0:6, 1:5, 2:4, 3:3, 2:2, 1:1", true},
				Chain{"0:2, 1:1, 1:3", true},
				Chain{"0:4, 1:3, 1:5", true},
				Chain{"1:3, 2:2, 2:4", true},
				Chain{"0:0, -1:1", false},
				Chain{"3:3, 5:5", false},
			}, []string{
				"0:0, 0:2, 1:1",
				"0:2, 0:4, 1:3",
				"0:4, 0:6, 1:5",
				"1:1, 1:3, 2:2",
				"1:3, 1:5, 2:4",
				"2:2, 2:4, 3:3",
			}, 2},

		// 5: A square divided into four subsquares. In this case we
		// want to extract the four loops rather than taking their
		// union.
		{0, -1, true, 0.0, 4.0, 90.0,
			[]Chain{
				Chain{"0:0, 0:5, 5:5, 5:0", true},
				Chain{"0:5, 0:10, 5:10, 5:5", true},
				Chain{"5:0, 5:5, 10:5, 10:0", true},
				Chain{"5:5, 5:10, 10:10, 10:5", true},
				Chain{"0:10, 0:15, 0:20", false},
				Chain{"20:0, 15:0, 10:0", false},
			}, []string{
				"0:0, 0:5, 5:5, 5:0",
				"0:5, 0:10, 5:10, 5:5",
				"5:0, 5:5, 10:5, 10:0",
				"5:5, 5:10, 10:10, 10:5",
			}, 4},

		// 6: Five nested loops that touch at a point.
		{0, 0, true, 0.0, 0.8, 5.0,
			[]Chain{
				Chain{"0:0, 0:10, 10:10, 10:0", true},
				Chain{"0:0, 1:9, 9:9, 9:1", true},
				Chain{"0:0, 2:8, 8:8, 8:2", true},
				Chain{"0:0, 3:7, 7:7, 7:3", true},
				Chain{"0:0, 4:6, 6:6, 6:4", true},
			}, []string{
				"0:0, 0:10, 10:10, 10:0",
				"0:0, 1:9, 9:9, 9:1",
				"0:0, 2:8, 8:8, 8:2",
				"0:0, 3:7, 7:7, 7:3",
				"0:0, 4:6, 6:6, 6:4",
			}, 0},

		// 7: Four diamonds nested within each other touching at two points.
		{-1, 0, true, 0.0, 4.0, 15.0, // Directed edges required for unique result.
			[]Chain{
				Chain{"0:-20, -10:0, 0:20, 10:0", true},
				Chain{"0:10, -10:0, 0:-10, 10:0", true},
				Chain{"0:-10, -5:0, 0:10, 5:0", true},
				Chain{"0:5, -5:0, 0:-5, 5:0", true},
			}, []string{
				"0:-20, -10:0, 0:-10, 10:0",
				"0:-10, -5:0, 0:-5, 5:0",
				"0:5, -5:0, 0:10, 5:0",
				"0:10, -10:0, 0:20, 10:0",
			}, 0},

		// 8: Seven diamonds nested within each other touching at one
		// point between each nested pair.
		{0, 0, true, 0.0, 9.0, 4.0,
			[]Chain{
				Chain{"0:-70, -70:0, 0:70, 70:0", true},
				Chain{"0:-70, -60:0, 0:60, 60:0", true},
				Chain{"0:-50, -60:0, 0:50, 50:0", true},
				Chain{"0:-40, -40:0, 0:50, 40:0", true},
				Chain{"0:-30, -30:0, 0:30, 40:0", true},
				Chain{"0:-20, -20:0, 0:30, 20:0", true},
				Chain{"0:-10, -20:0, 0:10, 10:0", true},
			}, []string{
				"0:-70, -70:0, 0:70, 70:0",
				"0:-70, -60:0, 0:60, 60:0",
				"0:-50, -60:0, 0:50, 50:0",
				"0:-40, -40:0, 0:50, 40:0",
				"0:-30, -30:0, 0:30, 40:0",
				"0:-20, -20:0, 0:30, 20:0",
				"0:-10, -20:0, 0:10, 10:0",
			}, 0},

		// 9: A triangle and a self-intersecting bowtie.
		{0, 0, false, 0.0, 4.0, 45.0,
			[]Chain{
				Chain{"0:0, 0:10, 5:5", true},
				Chain{"0:20, 0:30, 10:20", false},
				Chain{"10:20, 10:30, 0:20", false},
			}, []string{
				"0:0, 0:10, 5:5",
			}, 4},

		// 10: Two triangles that intersect each other.
		{0, 0, false, 0.0, 2.0, 45.0,
			[]Chain{
				Chain{"0:0, 0:12, 6:6", true},
				Chain{"3:6, 3:18, 9:12", true},
			}, []string{
				"",
			}, 6},

		// 11: Four squares that combine to make a big square.  The nominal edges of
		// the square are at +/-8.5 degrees in latitude and longitude.  All vertices
		// except the center vertex are perturbed by up to 0.5 degrees in latitude
		// and/or longitude.  The various copies of the center vertex are misaligned
		// by more than this (i.e. they are structured as a tree where adjacent
		// vertices are separated by at most 1 degree in latitude and/or longitude)
		// so that the clustering algorithm needs more than one iteration to find
		// them all.  Note that the merged position of this vertex doesn't matter
		// because it is XORed away in the output.  However, it's important that
		// all edge pairs that need to be XORed are separated by no more than
		// 'min_merge' below.

		{0, 1, true, 1.7, 5.8, 70.0, // XOR, min_merge > sqrt(2), max_merge < 6.
			[]Chain{
				Chain{"-8:-8, -8:0", false},
				Chain{"-8:1, -8:8", false},
				Chain{"0:-9, 1:-1", false},
				Chain{"1:2, 1:9", false},
				Chain{"0:8, 2:2", false},
				Chain{"0:-2, 1:-8", false},
				Chain{"8:9, 9:1", false},
				Chain{"9:0, 8:-9", false},
				Chain{"9:-9, 0:-8", false},
				Chain{"1:-9, -9:-9", false},
				Chain{"8:0, 1:0", false},
				Chain{"-1:1, -8:0", false},
				Chain{"-8:1, -2:0", false},
				Chain{"0:1, 8:1", false},
				Chain{"-9:8, 1:8", false},
				Chain{"0:9, 8:8", false},
			}, []string{
				"8.5:8.5, 8.5:0.5, 8.5:-8.5, 0.5:-8.5, -8.5:-8.5, -8.5:0.5, -8.5:8.5, 0.5:8.5",
			}, 0},
	}
	for i, test := range tests {
		if !runTestCase(t, test) {
			t.Errorf("#%d: runTestCase(t, %v) == false", i, test)
		}
	}
}

func runTestCase(t *testing.T, test TestCase) bool {
	for iter := 0; iter < 250; iter++ {
		options := PolygonBuilderOptions{
			edge_splice_fraction: .866,
			validate:             false,
			vertex_merge_radius:  s1.Angle(0),
		}
		options.undirected_edges = evalTristate(test.undirectedEdges)
		options.xor_edges = evalTristate(test.xorEdges)
		minMerge := (s1.Angle(test.minMerge) * s1.Degree).Radians()
		maxMerge := (s1.Angle(test.maxMerge) * s1.Degree).Radians()
		minSin := math.Sin((s1.Angle(test.minVertexAngle) * s1.Degree).Radians())

		// Half of the time we allow edges to be split into smaller
		// pieces (up to 5 levels, i.e. up to 32 pieces).
		maxSplits := max(0, rand.Intn(10)-4)
		if !test.canSplit {
			maxSplits = 0
		}

		// We choose randomly among two different values for the edge
		// fraction, just to exercise that code.
		edgeFraction := options.edge_splice_fraction
		var vertexMerge, maxPerturb float64
		if minSin < edgeFraction && oneIn(2) {
			edgeFraction = minSin
		}
		if maxSplits == 0 && oneIn(2) {
			// Turn off edge splicing completely.
			edgeFraction = 0
			vertexMerge = minMerge + smallFraction()*(maxMerge-minMerge)
			maxPerturb = 0.5 * math.Min(vertexMerge-minMerge, maxMerge-vertexMerge)
		} else {
			// Splice edges. These bounds also assume that edges
			// may be split.
			//
			// If edges are actually split, need to bump up the
			// minimum merge radius to ensure that split edges
			// in opposite directions are unified. Otherwise
			// there will be tiny degenerate loops created.
			if maxSplits > 0 {
				minMerge += 1e-15
			}
			minMerge /= edgeFraction
			maxMerge *= minSin
			if maxMerge < minMerge {
				t.Errorf("%v < %v", maxMerge, minMerge)
			}
			vertexMerge = minMerge + smallFraction()*(maxMerge-minMerge)
			maxPerturb = 0.5 * math.Min(edgeFraction*(vertexMerge-minMerge), maxMerge-vertexMerge)
		}

		// We can perturb by any amount up to the maximum, but choosing
		// a lower maximum decreases the error bounds when checking the
		// output.
		maxPerturb *= smallFraction()

		// This is the minimum length of a split edge to prevent
		// unexpected merging and/or splicing.
		minEdge := minMerge + (vertexMerge+2*maxPerturb)/minSin
		options.vertex_merge_radius = s1.Angle(vertexMerge)
		options.edge_splice_fraction = edgeFraction
		options.validate = true
		builder := NewPolygonBuilder(options)

		// On each iteration we randomly rotate the test case around
		// the sphere. This causes the PolygonBuilder to choose
		// different first edges when trying to build loops.
		x, y, z := randomFrame()
		m := r3.MatrixFromCols(x.Vector, y.Vector, z.Vector)
		for _, chain := range test.chainsIn {
			addChain(chain, m, maxSplits, maxPerturb, minEdge, &builder)
		}

		var loops []*Loop
		var unusedEdges []Edge
		if test.xorEdges < 0 {
			builder.AssembleLoops(&loops, &unusedEdges)
		} else {
			var polygon Polygon
			builder.AssemblePolygon(&polygon, &unusedEdges)
			polygon.Release(&loops)
		}

		expected := []*Loop{}
		for _, str := range test.loopsOut {
			if str != "" {
				vertices := getVertices(str, m)
				expected = append(expected, NewLoopFromPath(vertices))
			}
		}

		// We assume that the vertex locations in the expected output
		// polygon are separated from the corresponding vertex
		// locations in the input edges by at most half of the
		// minimum merge radius. Essentially this means that the
		// expected output vertices should be near the centroid of the
		// various input vertices.
		//
		// If any edges were split, we need to allow a bit more error
		// due to inaccuracies in the interpolated positions.
		// Similarly, if any vertices were perturbed, we need to bump
		// up the error to allow for numerical errors in the actual
		// perturbation.
		maxError := 0.5*minMerge + maxPerturb
		if maxSplits > 0 || maxPerturb > 0 {
			maxError += 1e-15
		}

		ok0 := findMissingLoops(loops, expected, m, maxSplits, maxError, "Actual")
		ok1 := findMissingLoops(expected, loops, m, maxSplits, maxError, "Expected")
		ok2 := unexpectedUnusedEdgeCount(len(unusedEdges), test.numUnusedEdges, maxSplits)
		if ok0 || ok1 || ok2 {
			// We found a problem.
			dumpUnusedEdges(unusedEdges, m, test.numUnusedEdges)
			fmt.Printf(`During iteration %d:
  undirected: %v
  xor: %v
  maxSplits: %d
  maxPerturb: %.6g
  vertexMergeRadius: %.6g
  edgeSpliceFraction: %.6g
  minEdge: %.6g
  maxError: %.6g

`, iter, options.undirected_edges, options.xor_edges, maxSplits,
				s1.Angle(maxPerturb).Degrees(),
				options.vertex_merge_radius.Degrees(),
				options.edge_splice_fraction,
				s1.Angle(minEdge).Degrees(),
				s1.Angle(maxError).Degrees())
			return false
		}
	}
	return true
}

func evalTristate(state int) bool {
	if state > 0 {
		return true
	}
	if state < 0 {
		return false
	}
	return oneIn(2)
}

func smallFraction() float64 {
	r := rand.Float64()
	u := rand.Float64()
	if r < 0.3 {
		return 0.0
	}
	if r < 0.6 {
		return u
	}
	return math.Pow(1e-10, u)
}

func randomFrame() (x, y, z Point) {
	x = randomPoint()
	y = Point{x.Cross(randomPoint().Vector).Normalize()}
	z = Point{x.Cross(y.Vector).Normalize()}
	return
}

func samplePoint(s2cap Cap) Point {
	// We consider the cap axis to be the "z" axis. We choose other
	// axes to complete the coordinate frame.
	m := FrameFromPoint(s2cap.center)
	// The surface area of a spherical cap is directly proportional to its
	// height. First we choose a random height, and then we choose a random
	// point along the circle at that height.
	h := rand.Float64() * s2cap.height
	theta := 2 * math.Pi * rand.Float64()
	r := math.Sqrt(h * (2 - h)) // radius of circle
	// The result should already be very close to unit-length, but we might
	// as well make it as accurate as possible.
	return PointFromFrame(m, PointFromCoords(math.Cos(theta)*r, math.Sin(theta)*r, 1-h))
}

func Perturb(x Point, maxPerturb float64) Point {
	// Perturb "x" randomly within the radius of maxPerturb.
	if maxPerturb == 0 {
		return x
	}
	s2cap := CapFromCenterAngle(Point{x.Normalize()}, s1.Angle(maxPerturb))
	return samplePoint(s2cap)
}

func addChain(chain Chain, m r3.Matrix, maxSplits int,
	maxPerturb, minEdge float64, builder *PolygonBuilder) {
	// Transform the given edge chain to the frame (x, y, z), optionally
	// split each edge into pieces and/or perturb the vertices up to the
	// given radius, and add them to the builder.
	vertices := getVertices(chain.str, m)
	if chain.closed {
		vertices = append(vertices, vertices[0])
	}
	for i := 1; i < len(vertices); i++ {
		addEdge(vertices[i-1], vertices[i], maxSplits, maxPerturb, minEdge, builder)
	}
}

func addEdge(v0, v1 Point, maxSplits int, maxPerturb, minEdge float64, builder *PolygonBuilder) {
	length := v0.Distance(v1).Radians()
	if maxSplits > 0 && oneIn(2) && length >= 2*minEdge {
		// Choose an interpolation parameter such that the length of
		// each piece is at least minEdge.
		f := minEdge / length
		t := f + (1-2*f)*rand.Float64()

		// Now add the two sub-edges recursively.
		vmid := EdgeInterpolate(t, v0, v1)
		addEdge(v0, vmid, maxSplits-1, maxPerturb, minEdge, builder)
		addEdge(vmid, v1, maxSplits-1, maxPerturb, minEdge, builder)
	} else {
		builder.AddEdge(Perturb(v0, maxPerturb), Perturb(v1, maxPerturb))
	}
}

func getVertices(s string, m r3.Matrix) (vertices []Point) {
	line := makePolyline(s)
	for _, v := range line.vertices {
		vertices = append(vertices, Point{m.MulVector(v.Vector).Normalize()})
	}
	return
}

func makePolyline(s string) *Polyline {
	vertices := parsePoints(s)
	return PolylineFromPoints(vertices)
}

func findMissingLoops(actual, expected []*Loop,
	m r3.Matrix, maxSplits int, maxError float64, label string) bool {
	// Dump any loops from "actual" that are not present in "expected".
	found := false
	for i, loop := range actual {
		if findLoop(loop, expected, maxSplits, maxError) {
			continue
		}
		fmt.Printf("%s loop %d:\n", label, i)
		for _, v := range loop.vertices {
			ll := LatLngFromPoint(Point{m.Transpose().MulVector(v.Vector)})
			fmt.Printf("  [%.6f, %.6f]\n", ll.Lat.Degrees(), ll.Lng.Degrees())
		}
		found = true
	}
	return found
}

func unexpectedUnusedEdgeCount(numActual, numExpected, maxSplits int) bool {
	if maxSplits == 0 {
		return numActual != numExpected
	} else {
		return (numActual > 0) != (numExpected > 0)
	}
}

func dumpUnusedEdges(unusedEdges []Edge, m r3.Matrix, numExpected int) {
	if len(unusedEdges) == numExpected {
		return
	}
	fmt.Printf("Wrong number of unused edges (%d expected, %d actual):\n",
		numExpected, len(unusedEdges))
	for _, edge := range unusedEdges {
		p0 := LatLngFromPoint(Point{m.Transpose().MulVector(edge.v0.Vector)})
		p1 := LatLngFromPoint(Point{m.Transpose().MulVector(edge.v1.Vector)})
		fmt.Printf("  [%.6f, %.6f] -> [%.6f, %.6f]\n",
			p0.Lat.Degrees(), p0.Lng.Degrees(),
			p1.Lat.Degrees(), p1.Lng.Degrees())
	}
}

func findLoop(loop *Loop, candidates []*Loop, maxSplits int, maxError float64) bool {
	for _, c := range candidates {
		if maxSplits == 0 {
			if loop.BoundaryApproxEquals(c, maxError) {
				return true
			}
		} else {
			if loop.BoundaryNear(c, maxError) {
				return true
			}
		}
	}
	return false
}
