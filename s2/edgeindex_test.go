package s2

import (
	"fmt"
	//	"log"
	"math/rand"
	"testing"

	"github.com/davidreynolds/gos2/s1"
)

type EdgeVectorIndex struct {
	EdgeIndex
	edges *[]Edge
}

func (idx EdgeVectorIndex) NumEdges() int          { return len(*idx.edges) }
func (idx EdgeVectorIndex) edge_from(i int) *Point { return &(*idx.edges)[i].v1 }
func (idx EdgeVectorIndex) edge_to(i int) *Point   { return &(*idx.edges)[i].v0 }
func (idx *EdgeVectorIndex) IncrementQueryCount()  { idx.EdgeIndex.IncrementQueryCount() }

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
		candidates := map[int]bool{}
		for it.GetCandidates(e.v0, e.v1); !it.Done(); it.Next() {
			candidates[it.Index()] = true
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
	allEdges := []Edge{}
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

func TestDegenerateEdgeOnCellVertexIsItsOwnCandidate(t *testing.T) {
	for i := 0; i < 100; i++ {
		cell := CellFromCellID(randomCellID())
		e := Edge{cell.Vertex(0), cell.Vertex(0)}
		allEdges := []Edge{}
		allEdges = append(allEdges, e)
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
	allEdges := []Edge{}
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
		p2_p1 := Point{p2.Sub(p1.Vector).Mul(1 / numPointsOnEdge)}
		allEdges := []Edge{}
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
