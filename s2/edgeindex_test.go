package s2

import (
	"fmt"
	"math/rand"
	"testing"

	"github.com/davidreynolds/gos2/s1"
)

type EdgeVectorIndex struct {
	*EdgeIndex
	edges *[]Edge
}

func (idx *EdgeVectorIndex) NumEdges() int        { return len(*idx.edges) }
func (idx *EdgeVectorIndex) IncrementQueryCount() { idx.EdgeIndex.IncrementQueryCount() }
func (idx *EdgeVectorIndex) EdgeFromTo(i int) (*Point, *Point) {
	return &(*idx.edges)[i].v1, &(*idx.edges)[i].v0
}

func EdgesToString(a, b Edge) string {
	_, u1, v1 := xyzToFaceUV(a.v0.Vector)
	_, u2, v2 := xyzToFaceUV(a.v1.Vector)
	_, u3, v3 := xyzToFaceUV(b.v0.Vector)
	_, u4, v4 := xyzToFaceUV(b.v1.Vector)
	u2 = u2 - u1
	u3 = u3 - u1
	u4 = u4 - u1
	v2 = v2 - v1
	v3 = v3 - v1
	v4 = v4 - v1
	return fmt.Sprintf("[%6.2f,%6.2f] 0,0 -> %3.3e,%3.3e) -- (%3.3e,%3.3e -> %3.3e %3.3e)",
		u1, v1, u2, v2, u3, v3, u4, v4)
}

func CheckAllCrossings(t *testing.T, allEdges []Edge, minCrossings, maxChecksCrossingsRatio int) {
	index := EdgeVectorIndex{NewEdgeIndex(), &allEdges}
	ComputeIndex(&index)
	it := NewEdgeIndexIterator(&index)
	totalCrossings := 0
	totalChecks := 0

	for in := 0; in < len(allEdges); in++ {
		e := allEdges[in]
		candidates := make(map[int]struct{})
		for it.GetCandidates(e.v0, e.v1); !it.Done(); it.Next() {
			candidates[it.Index()] = struct{}{}
			totalChecks++
		}
		for i := 0; i < len(allEdges); i++ {
			crossing := RobustCrossing(e.v0, e.v1, allEdges[i].v0, allEdges[i].v1)
			if crossing >= 0 {
				if _, ok := candidates[i]; !ok {
					t.Errorf("Edge %d is not a candidate of edge %d (%s)",
						i, in, EdgesToString(allEdges[i], e))
				}
				totalCrossings++
			}
		}
	}
	if minCrossings > totalCrossings {
		t.Errorf("minCrossings (%v) > totalCrossings (%v)", minCrossings, totalCrossings)
	}
	if got := totalCrossings * maxChecksCrossingsRatio; got < totalChecks {
		t.Errorf("%v * %v == %v, want >= %v", totalCrossings, maxChecksCrossingsRatio, got, totalChecks)
	}
}

func TestLoopCandidateOfItself(t *testing.T) {
	ps := []Point{}
	ps = append(ps, makepoint("0:178"))
	ps = append(ps, makepoint("-1:180"))
	ps = append(ps, makepoint("0:-179"))
	ps = append(ps, makepoint("1:-180"))
	allEdges := []Edge{}
	for i := 0; i < 4; i++ {
		allEdges = append(allEdges, Edge{ps[i], ps[(i+1)%4]})
	}
	CheckAllCrossings(t, allEdges, 0, 16)
}

const kEarthRadiusMeters = 6371000.0

func randomEdgeCrossingCap(maxLengthMeters float64, s2cap Cap) Edge {
	center := samplePoint(s2cap)
	edgeCap := CapFromCenterAngle(center, s1.Angle(maxLengthMeters/kEarthRadiusMeters/2))
	p0 := samplePoint(edgeCap)
	p1 := samplePoint(edgeCap)
	return Edge{p0, p1}
}

func GenerateRandomEarthEdges(edgeLengthMetersMax, capSpanMeters float64, numEdges int, edges *[]Edge) {
	s2cap := CapFromCenterAngle(randomPoint(), s1.Angle(capSpanMeters/kEarthRadiusMeters))
	for i := 0; i < numEdges; i++ {
		*edges = append(*edges, randomEdgeCrossingCap(edgeLengthMetersMax, s2cap))
	}
}

func CheckCrossingsRandomInCap(t *testing.T, numEdges int, edgeLengthMax, capSpanMeters float64,
	minCrossings, maxChecksCrossingsRatio int) {
	var allEdges []Edge
	GenerateRandomEarthEdges(edgeLengthMax, capSpanMeters, numEdges, &allEdges)
	CheckAllCrossings(t, allEdges, minCrossings, maxChecksCrossingsRatio)
}

func TestRandomEdgeCrossings(t *testing.T) {
	CheckCrossingsRandomInCap(t, 2000, 30, 5000, 500, 2)
	CheckCrossingsRandomInCap(t, 1000, 100, 5000, 500, 3)
	CheckCrossingsRandomInCap(t, 1000, 1000, 5000, 1000, 40)
	CheckCrossingsRandomInCap(t, 500, 5000, 5000, 5000, 20)
}

func TestRandomEdgeCrossingsSparse(t *testing.T) {
	for i := 0; i < 5; i++ {
		CheckCrossingsRandomInCap(t, 2000, 100, 5000, 500, 8)
		CheckCrossingsRandomInCap(t, 2000, 300, 50000, 1000, 10)
	}
}

func ComputeCrossings(b *testing.B, numEdges int, maxEdgeLen, capSpanMeters float64, bruteForce bool) {
	b.StopTimer()
	var allEdges []Edge
	GenerateRandomEarthEdges(maxEdgeLen, capSpanMeters, numEdges, &allEdges)
	b.StartTimer()
	if bruteForce {
		for _, e0 := range allEdges {
			for _, e1 := range allEdges {
				RobustCrossing(e0.v0, e0.v1, e1.v0, e1.v1)
			}
		}
	} else {
		index := EdgeVectorIndex{NewEdgeIndex(), &allEdges}
		ComputeIndex(&index)
		it := NewEdgeIndexIterator(&index)
		for _, e := range allEdges {
			for it.GetCandidates(e.v0, e.v1); !it.Done(); it.Next() {
				in := it.Index()
				RobustCrossing(allEdges[in].v0, allEdges[in].v1, e.v0, e.v1)
			}
		}
	}
}

func benchmarkCrossingsSparse(b *testing.B, numEdges int, bruteForce bool) {
	for i := 0; i < b.N; i++ {
		ComputeCrossings(b, numEdges, 30, 5000, bruteForce)
		ComputeCrossings(b, numEdges, 100, 5000, bruteForce)
	}
}

func BenchmarkCrossingsSparse10t(b *testing.B)    { benchmarkCrossingsSparse(b, 10, true) }
func BenchmarkCrossingsSparse10f(b *testing.B)    { benchmarkCrossingsSparse(b, 10, false) }
func BenchmarkCrossingsSparse50t(b *testing.B)    { benchmarkCrossingsSparse(b, 50, true) }
func BenchmarkCrossingsSparse50f(b *testing.B)    { benchmarkCrossingsSparse(b, 50, false) }
func BenchmarkCrossingsSparse100t(b *testing.B)   { benchmarkCrossingsSparse(b, 100, true) }
func BenchmarkCrossingsSparse100f(b *testing.B)   { benchmarkCrossingsSparse(b, 100, false) }
func BenchmarkCrossingsSparse300t(b *testing.B)   { benchmarkCrossingsSparse(b, 300, true) }
func BenchmarkCrossingsSparse300f(b *testing.B)   { benchmarkCrossingsSparse(b, 300, false) }
func BenchmarkCrossingsSparse600t(b *testing.B)   { benchmarkCrossingsSparse(b, 600, true) }
func BenchmarkCrossingsSparse600f(b *testing.B)   { benchmarkCrossingsSparse(b, 600, false) }
func BenchmarkCrossingsSparse1000t(b *testing.B)  { benchmarkCrossingsSparse(b, 1000, true) }
func BenchmarkCrossingsSparse1000f(b *testing.B)  { benchmarkCrossingsSparse(b, 1000, false) }
func BenchmarkCrossingsSparse2000f(b *testing.B)  { benchmarkCrossingsSparse(b, 2000, false) }
func BenchmarkCrossingsSparse5000f(b *testing.B)  { benchmarkCrossingsSparse(b, 5000, false) }
func BenchmarkCrossingsSparse10000f(b *testing.B) { benchmarkCrossingsSparse(b, 10000, false) }
func BenchmarkCrossingsSparse20000f(b *testing.B) { benchmarkCrossingsSparse(b, 20000, false) }

func TestDegenerateEdgeOnCellVertexIsItsOwnCandidate(t *testing.T) {
	for i := 0; i < 100; i++ {
		cell := CellFromCellID(randomCellID())
		e := Edge{cell.Vertex(0), cell.Vertex(0)}
		allEdges := []Edge{e}
		index := EdgeVectorIndex{NewEdgeIndex(), &allEdges}
		ComputeIndex(&index)
		it := NewEdgeIndexIterator(&index)
		it.GetCandidates(e.v0, e.v1)
		if it.Done() {
			t.Fatalf("it.Done()")
		}
		it.Next()
		if !it.Done() {
			t.Fatalf("!it.Done()")
		}
	}
}

func TestLongEdgeCrossesLoop(t *testing.T) {
	var allEdges []Edge
	center := makepoint("42:107")
	testEdge := Edge{OriginPoint(), center}
	loop := makeRegularLoop(center, 4, 7e-3)
	for i := 0; i < len(loop.vertices); i++ {
		allEdges = append(allEdges, Edge{*loop.vertex(i), *loop.vertex(i + 1)})
	}
	index := EdgeVectorIndex{NewEdgeIndex(), &allEdges}
	ComputeIndex(&index)
	it := NewEdgeIndexIterator(&index)
	it.GetCandidates(testEdge.v0, testEdge.v1)
	if it.Done() {
		t.Fatalf("Should have at least one candidate")
	}
}

func TestCollinearEdgesOnCellBoundaries(t *testing.T) {
	alwaysRecurseOnChildren = true
	const numPointsOnEdge = 8
	for level := 0; level <= maxLevel; level++ {
		cell := CellFromCellID(randomCellIDForLevel(level))
		v1 := rand.Intn(4)
		v2 := (v1 + 1) & 3
		p1 := cell.Vertex(v1)
		p2 := cell.Vertex(v2)
		p2_p1 := Point{p2.Sub(p1.Vector).Mul(1. / numPointsOnEdge)}
		var allEdges []Edge
		points := make([]Point, numPointsOnEdge+1)
		for i := 0; i <= numPointsOnEdge; i++ {
			points[i] = Point{(p1.Add(p2_p1.Mul(float64(i)))).Normalize()}
			for j := 0; j < i; j++ {
				allEdges = append(allEdges, Edge{points[i], points[j]})
			}
		}
		CheckAllCrossings(t, allEdges, numPointsOnEdge*numPointsOnEdge,
			numPointsOnEdge*numPointsOnEdge)
	}
}

func BenchmarkEdgeCovering(b *testing.B) {
	b.StopTimer()
	const numVerts int = 1000
	loopCenter := makepoint("42:107")
	loop := makeRegularLoop(loopCenter, numVerts, 7e-3)
	var cover []CellID
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		idx := NewLoopIndex(loop)
		for j := 0; j < idx.NumEdges(); j++ {
			from, to := idx.EdgeFromTo(j)
			EdgeCovering(*from, *to, true, &cover)
		}
	}
}

func BenchmarkQuadTreeInsertionCost(b *testing.B) {
	b.StopTimer()
	const numVerts int = 1000
	loopCenter := makepoint("42:107")
	loop := makeRegularLoop(loopCenter, numVerts, 7e-3)
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		index := NewLoopIndex(loop)
		ComputeIndex(index)
	}
}

func benchmarkQuadTreeFindCost(b *testing.B, numVerts int) {
	b.StopTimer()
	loopCenter := makepoint("42:107")
	loop := makeRegularLoop(loopCenter, numVerts, 7e-3)
	p := makepoint("42:106.99")
	q := makepoint("42:107.01")
	index := NewLoopIndex(loop)
	iter := NewEdgeIndexIterator(index)
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		iter.GetCandidates(p, q)
	}
}

func BenchmarkQuadTreeFindCost10(b *testing.B)     { benchmarkQuadTreeFindCost(b, 10) }
func BenchmarkQuadTreeFindCost100(b *testing.B)    { benchmarkQuadTreeFindCost(b, 100) }
func BenchmarkQuadTreeFindCost1000(b *testing.B)   { benchmarkQuadTreeFindCost(b, 1000) }
func BenchmarkQuadTreeFindCost10000(b *testing.B)  { benchmarkQuadTreeFindCost(b, 10000) }
func BenchmarkQuadTreeFindCost100000(b *testing.B) { benchmarkQuadTreeFindCost(b, 100000) }
