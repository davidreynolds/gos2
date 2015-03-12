package s2

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	//	"unsafe"
)

func oneIn(x int) bool { return rand.Intn(x) == 0 }

type LevelStats struct {
	count                                    float64
	minArea, maxArea, avgArea                float64
	minWidth, maxWidth, avgWidth             float64
	minEdge, maxEdge, avgEdge, maxEdgeAspect float64
	minDiag, maxDiag, avgDiag, maxDiagAspect float64
	minAngleSpan, maxAngleSpan, avgAngleSpan float64
	minApproxRatio, maxApproxRatio           float64
}

func NewLevelStats() LevelStats {
	return LevelStats{
		minArea:        100,
		minWidth:       100,
		minEdge:        100,
		minDiag:        100,
		minAngleSpan:   100,
		minApproxRatio: 100,
	}
}

var (
	level_stats []LevelStats
)

func init() {
	level_stats = make([]LevelStats, maxLevel+1)
	for i := 0; i <= maxLevel; i++ {
		level_stats[i] = NewLevelStats()
	}
}

func GatherStats(cell Cell) {
	s := &level_stats[cell.level]
	exactArea := cell.ExactArea()
	approxArea := cell.ApproxArea()
	var minEdge, maxEdge, avgEdge float64
	var minDiag, maxDiag float64
	var minWidth, maxWidth float64
	var minAngleSpan, maxAngleSpan float64

	minEdge = 100
	minDiag = 100
	minWidth = 100
	minAngleSpan = 100

	for i := 0; i < 4; i++ {
		edge := float64(cell.VertexRaw(i).Angle(cell.VertexRaw((i + 1) & 3).Vector))
		minEdge = math.Min(edge, minEdge)
		maxEdge = math.Max(edge, maxEdge)
		avgEdge += 0.25 * edge
		mid := cell.VertexRaw(i).Add(cell.VertexRaw((i + 1) & 3).Vector)
		width := M_PI_2 - float64(mid.Angle(cell.EdgeRaw(i^2).Vector))
		minWidth = math.Min(width, minWidth)
		maxWidth = math.Max(width, maxWidth)
		if i < 2 {
			diag := float64(cell.VertexRaw(i).Angle(cell.VertexRaw(i ^ 2).Vector))
			minDiag = math.Min(diag, minDiag)
			maxDiag = math.Max(diag, maxDiag)
			angleSpan := float64(cell.EdgeRaw(i).Angle(cell.EdgeRaw(i ^ 2).Neg()))
			minAngleSpan = math.Min(angleSpan, minAngleSpan)
			maxAngleSpan = math.Max(angleSpan, maxAngleSpan)
		}
	}

	s.count++
	s.minArea = math.Min(exactArea, s.minArea)
	s.maxArea = math.Max(exactArea, s.maxArea)
	s.avgArea += exactArea

	s.minWidth = math.Min(minWidth, s.minWidth)
	s.maxWidth = math.Max(maxWidth, s.maxWidth)
	s.avgWidth += 0.5 * (minWidth + maxWidth)

	s.minEdge = math.Min(minEdge, s.minEdge)
	s.maxEdge = math.Max(maxEdge, s.maxEdge)
	s.avgEdge += avgEdge
	s.maxEdgeAspect = math.Max(maxEdge/minEdge, s.maxEdgeAspect)

	s.minDiag = math.Min(minDiag, s.minDiag)
	s.maxDiag = math.Max(maxDiag, s.maxDiag)
	s.avgDiag += 0.5 * (minDiag + maxDiag)
	s.maxDiagAspect = math.Max(maxDiag/minDiag, s.maxDiagAspect)

	s.minAngleSpan = math.Min(minAngleSpan, s.minAngleSpan)
	s.maxAngleSpan = math.Max(maxAngleSpan, s.maxAngleSpan)
	s.avgAngleSpan += 0.5 * (minAngleSpan + maxAngleSpan)
	approxRatio := approxArea / exactArea
	s.minApproxRatio = math.Min(approxRatio, s.minApproxRatio)
	s.maxApproxRatio = math.Max(approxRatio, s.maxApproxRatio)
}

func DoSubdivide(t *testing.T, cell Cell) {
	GatherStats(cell)
	if cell.IsLeaf() {
		return
	}
	children := make([]Cell, 0, 4)
	if !cell.Subdivide(&children) {
		t.Errorf("cell.Subdivide failed")
	}
	childId := cell.id.ChildBegin()
	var exactArea, approxArea, avgArea float64
	for i := 0; i < 4; i, childId = i+1, childId.Next() {
		child := children[i]
		exactArea += child.ExactArea()
		approxArea += child.ApproxArea()
		avgArea += child.AverageArea()

		// Check that the child geometry is consistent with its cellid.
		if childId != child.id {
			t.Errorf("%v != %v", childId, child.id)
		}
		if !child.Center().ApproxEqual(childId.Point()) {
			t.Errorf("%v.ApproxEqual(%v) == false", child.Center(), childId.Point())
		}
		direct := CellFromCellID(childId)
		if direct.face != child.face {
			t.Errorf("%v != %v", direct.face, child.face)
		}
		if direct.level != child.level {
			t.Errorf("%v != %v", direct.level, child.level)
		}
		if direct.orientation != child.orientation {
			t.Errorf("%v != %v", direct.orientation, child.orientation)
		}
		if direct.CenterRaw() != child.CenterRaw() {
			t.Errorf("%v != %v", direct.CenterRaw(), child.CenterRaw())
		}

		for k := 0; k < 4; k++ {
			dv := direct.VertexRaw(k)
			cv := child.VertexRaw(k)
			if dv != cv {
				t.Errorf("%v != %v", dv, cv)
			}

			de := direct.EdgeRaw(k)
			ce := child.EdgeRaw(k)
			if de != ce {
				t.Errorf("%v != %v", de, ce)
			}
		}

		// Test Contains() and MayIntersect()
		if !cell.ContainsCell(child) {
			t.Errorf("!%v.ContainsCell(%v)", cell, child)
		}
		if !cell.MayIntersect(child) {
			t.Errorf("!%v.MayIntersect(%v)", cell, child)
		}
		if child.ContainsCell(cell) {
			t.Errorf("%v.ContainsCell(%v)", child, cell)
		}
		if !cell.ContainsPoint(child.CenterRaw()) {
			t.Errorf("!%v.ContainsPoint(%v)", cell, child.CenterRaw())
		}
		for j := 0; j < 4; j++ {
			if !cell.ContainsPoint(child.VertexRaw(j)) {
				t.Errorf("!%v.ContainsPoint(%v)", cell, child.VertexRaw(j))
			}
			if j != i {
				if child.ContainsPoint(children[j].CenterRaw()) {
					t.Errorf("%v.ContainsPoint(%v)", child, children[j].CenterRaw())
				}
				if child.MayIntersect(children[j]) {
					t.Errorf("%v.MayIntersect(%v)", child, children[j])
				}
			}
		}

		// Test CapBound and RectBound
		parentRect := cell.RectBound()
		if cell.ContainsPoint(PointFromCoords(0, 0, 1)) ||
			cell.ContainsPoint(PointFromCoords(0, 0, -1)) {
			if !parentRect.Lng.IsFull() {
				t.Errorf("!%v.Lng.IsFull()", parentRect)
			}
		}

		parentCap := cell.CapBound()
		childCap := child.CapBound()
		childRect := child.RectBound()
		tests := []struct {
			a    Region
			b    Point
			want bool
		}{
			{childCap, child.Center(), true},
			{childRect, child.CenterRaw(), true},
			{parentCap, child.Center(), true},
			{parentRect, child.CenterRaw(), true},
		}
		for i, test := range tests {
			got := test.a.ContainsPoint(test.b)
			if got != test.want {
				t.Errorf("(%d) %v.ContainsPoint(%v), got = %v, want = %v",
					i, test.a, test.b, got, test.want)
			}
		}

		for j := 0; j < 4; j++ {
			tests := []struct {
				a    Region
				b    Point
				want bool
			}{
				{childCap, child.Vertex(j), true},
				{childRect, child.Vertex(j), true},
				{childRect, child.VertexRaw(j), true},
				{parentCap, child.Vertex(j), true},
				{parentRect, child.Vertex(j), true},
				{parentRect, child.VertexRaw(j), true},
			}
			for _, test := range tests {
				got := test.a.ContainsPoint(test.b)
				if got != test.want {
					t.Errorf("%v.ContainsPoint(%v), got = %v, want = %v",
						test.a, test.b, got, test.want)
				}
			}

			if j != i {
				// The bounding caps and rectangles should be tight enough so that
				// they exclude at least two vertices of each adjacent cell.
				var capCount, rectCount int
				for k := 0; k < 4; k++ {
					if childCap.ContainsPoint(children[j].Vertex(k)) {
						capCount++
					}
					if childRect.ContainsPoint(children[j].VertexRaw(k)) {
						rectCount++
					}
				}
				if capCount > 2 {
					t.Errorf("%v > 2", capCount)
				}
				if childRect.Lat.Lo > -M_PI_2 &&
					childRect.Lat.Hi < M_PI_2 {
					// Bounding rectangles may be too large at the poles because the
					// pole itself has an arbitrary fixed longitude.
					if rectCount > 2 {
						t.Errorf("%v > 2", rectCount)
					}
				}
			}
		}

		// Check all children for the first few levels, and then sample randomly.
		// Also subdivide one corner cell, one edge cell, and one center cell.
		forceSubdivide := false
		face := int(child.face)
		center := Point{faceNorm(face)}
		edge := Point{center.Add(uAxis(face))}
		corner := Point{edge.Add(vAxis(face))}
		for j := 0; j < 4; j++ {
			p := child.VertexRaw(j)
			if p == center || p == edge || p == corner {
				forceSubdivide = true
			}
		}
		// This run takes awhile at non-debug levels...
		var level int8
		var chance int
		level = 5
		chance = 4
		if forceSubdivide || cell.level < level || oneIn(chance) {
			DoSubdivide(t, child)
		}
	}

	// Check sum of child areas equals parent area.
	//
	// For ExactArea(), the best relative error we can expect is about
	// 1e-6 because the precision of the unit vector coordinates is only
	// about 1e-15 and the edge length of a leaf cell is about 1e-9.
	//
	// For ApproxArea(), the areas are accurate within a few percent.
	//
	// For AverageArea(), the areas themselves are not very accurate, but
	// the average area of a parent is exactly 4 times the area of a child.
	tests := []struct {
		name string
		a    float64
		b    float64
		want float64
	}{
		{"ExactArea", exactArea, cell.ExactArea(), math.Abs(math.Log(1 + 1e-6))},
		{"ApproxArea", approxArea, cell.ApproxArea(), math.Abs(math.Log(1.03))},
		{"AverageArea", avgArea, cell.AverageArea(), math.Abs(math.Log(1 + 1e-15))},
	}
	for _, test := range tests {
		got := math.Abs(math.Log(test.a / test.b))
		if got > test.want {
			t.Errorf("%s: %v > %v", test.name, got, test.want)
		}
	}
}

func TestSubdivide(t *testing.T) {
	for face := 0; face < 6; face++ {
		id := CellIDFromFacePosLevel(face, 0, 0)
		cell := CellFromCellID(id)
		DoSubdivide(t, cell)
	}

	// The max edge *ratio* is the ratio of the longest edge of any cell to
	// the shortest edge of any cell at the same level (and similarly for
	// the max diagonal ratio).
	//
	// The max edge *aspect* is the max ratio of the longest edge of a cell
	// to the shortest edge of that same cell (and similarly for the max
	// diagonal aspect).
	fmt.Println("Level     Area      Edge          Diag          Approx       Average")
	fmt.Println("         Ratio  Ratio Aspect  Ratio Aspect    Min    Max   Min     Max")
	for i := 0; i <= maxLevel; i++ {
		s := &level_stats[i]
		if s.count > 0 {
			s.avgArea /= s.count
			s.avgWidth /= s.count
			s.avgEdge /= s.count
			s.avgDiag /= s.count
			s.avgAngleSpan /= s.count
		}
		fmt.Printf("%5d  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
			i, s.maxArea/s.minArea,
			s.maxEdge/s.minEdge, s.maxEdgeAspect,
			s.maxDiag/s.minDiag, s.maxDiagAspect,
			s.minApproxRatio, s.maxApproxRatio,
			AverageArea(i)/s.maxArea,
			AverageArea(i)/s.minArea)
	}

	// Now check the validity of the S2 length and area metrics.
	for i := 0; i <= maxLevel; i++ {
		s := &level_stats[i]
		if s.count == 0 {
			continue
		}
		fmt.Printf("Level %2d - metric (error/actual : error/tolerance)\n", i)

		// The various length calculations are only accurate to 1e-15
		// or so, so we need to allow for this amount of discrepancy
		// with the theoretical minimums and maximums. The area
		// calculation is accurate to about 1e-15 times the cell width.
		CheckMinMaxAvg(t, "area", i, s.count, 1e-15*s.minWidth,
			s.minArea, s.maxArea, s.avgArea,
			MinArea.Metric, MaxArea.Metric, AvgArea.Metric)
		CheckMinMaxAvg(t, "width", i, s.count, 1e-15,
			s.minWidth, s.maxWidth, s.avgWidth,
			MinWidth.Metric, MaxWidth.Metric, AvgWidth.Metric)
		CheckMinMaxAvg(t, "edge", i, s.count, 1e-15,
			s.minEdge, s.maxEdge, s.avgEdge,
			MinEdge.Metric, MaxEdge.Metric, AvgEdge.Metric)
		CheckMinMaxAvg(t, "diagonal", i, s.count, 1e-15,
			s.minDiag, s.maxDiag, s.avgDiag,
			MinDiag.Metric, MaxDiag.Metric, AvgDiag.Metric)
		CheckMinMaxAvg(t, "angle span", i, s.count, 1e-15,
			s.minAngleSpan, s.maxAngleSpan, s.avgAngleSpan,
			MinAngleSpan.Metric, MaxAngleSpan.Metric, AvgAngleSpan.Metric)

		var want float64
		shift := 1 << uint(i)
		want = MaxEdgeAspect + 1e-15*float64(shift)
		if s.maxEdgeAspect > want {
			t.Errorf("maxEdgeAspect: %v > %v", s.maxEdgeAspect, want)
		}

		want = MaxDiagAspect + 1e-15*float64(shift)
		if s.maxDiagAspect > want {
			t.Errorf("maxDiagAspect: %v > %v", s.maxDiagAspect, want)
		}
	}
}

func CheckMinMaxAvg(t *testing.T, label string, level int,
	count, absError, minValue, maxValue, avgValue float64,
	minMetric, maxMetric, avgMetric Metric) {
	// All metrics are minimums, maximums, or averages of differential
	// quantities, and therefore will not be exact for cells at any finite
	// level. The differential minimum is always a lower bound, and the
	// maximum is always an upper bound, but these minimums and maximums
	// may not be achieved for two different reasons. First, the cells at
	// each level are sampled and we may miss the most extreme examples.
	// Second, the actual metric for a cell is obtained by integrating the
	// differential quantity, which is no constant across the cell.
	// Therefore cells at low levels (bigger cells) have smaller variations.
	//
	// The "tolerance" below is an attempt to model both of these effects.
	// At low levels, error is dominated by the variation of differential
	// quantities across the cells, while at high levels error is dominated
	// by the effects of random sampling.
	shift := 1 << uint(level)
	tolerance := maxMetric.Value(level) - minMetric.Value(level)
	tolerance /= math.Sqrt(math.Min(count, 0.5*float64(shift)))
	if tolerance == 0 {
		tolerance = absError
	}
	minError := minValue - minMetric.Value(level)
	maxError := maxMetric.Value(level) - maxValue
	avgError := math.Abs(avgMetric.Value(level) - avgValue)
	fmt.Printf("%-10s (%6.0f samples, tolerance %8.3g) - min (%9.3g : %9.3g) max (%9.3g : %9.3g), avg (%9.3g : %9.3g)\n",
		label, count, tolerance,
		minError/minValue, minError/tolerance,
		maxError/maxValue, maxError/tolerance,
		avgError/avgValue, avgError/tolerance)

	if minMetric.Value(level) > minValue+absError {
		t.Errorf("%v > %v", minMetric.Value(level), minValue+absError)
	}
	if minMetric.Value(level) < minValue-tolerance {
		t.Errorf("%v < %v", minMetric.Value(level), minValue-tolerance)
	}
	if maxMetric.Value(level) > maxValue+tolerance {
		fmt.Printf("%v > %v\n", maxMetric.Value(level), maxValue+tolerance)
		t.Errorf("%v > %v", maxMetric.Value(level), maxValue+tolerance)
	}
	if maxMetric.Value(level) < maxValue-absError {
		t.Errorf("%v < %v", maxMetric.Value(level), maxValue-absError)
	}
	if math.Abs(avgMetric.Value(level)-avgValue) > 10*tolerance {
		t.Errorf("%v - %v > %v", avgMetric.Value(level), avgValue, 10*tolerance)
	}
}

// maxCellSize is the upper bounds on the number of bytes we want the Cell object to ever be.
/*
const maxCellSize = 48

func TestCellObjectSize(t *testing.T) {
	if sz := unsafe.Sizeof(Cell{}); sz > maxCellSize {
		t.Errorf("Cell struct too big: %d bytes > %d bytes", sz, maxCellSize)
	}
}
*/

func TestCellFaces(t *testing.T) {
	edgeCounts := make(map[Point]int)
	vertexCounts := make(map[Point]int)

	for face := 0; face < 6; face++ {
		id := CellIDFromFace(face)
		cell := CellFromCellID(id)

		if cell.id != id {
			t.Errorf("cell.id != id; %v != %v", cell.id, id)
		}

		if cell.face != int8(face) {
			t.Errorf("cell.face != face: %v != %v", cell.face, face)
		}

		if cell.level != 0 {
			t.Errorf("cell.level != 0: %v != 0", cell.level)
		}

		// Top-level faces have alternating orientations to get RHS coordinates.
		if cell.orientation != int8(face&swapMask) {
			t.Errorf("cell.orientation != orientation: %v != %v", cell.orientation, face&swapMask)
		}

		if cell.IsLeaf() {
			t.Errorf("cell should not be a leaf: IsLeaf = %v", cell.IsLeaf())
		}
		for k := 0; k < 4; k++ {
			edgeCounts[cell.Edge(k)]++
			vertexCounts[cell.Vertex(k)]++
			if d := cell.Vertex(k).Dot(cell.Edge(k).Vector); !float64Eq(0.0, d) {
				t.Errorf("dot product of vertex and edge failed, got %v, want 0", d)
			}
			if d := cell.Vertex((k + 1) & 3).Dot(cell.Edge(k).Vector); !float64Eq(0.0, d) {
				t.Errorf("dot product for edge and next vertex failed, got %v, want 0", d)
			}
			if d := cell.Vertex(k).Vector.Cross(cell.Vertex((k + 1) & 3).Vector).Normalize().Dot(cell.Edge(k).Vector); !float64Eq(1.0, d) {
				t.Errorf("dot product of cross product for vertices failed, got %v, want 1.0", d)
			}
		}
	}

	// Check that edges have multiplicity 2 and vertices have multiplicity 3.
	for k, v := range edgeCounts {
		if v != 2 {
			t.Errorf("edge %v counts wrong, got %d, want 2", k, v)
		}
	}
	for k, v := range vertexCounts {
		if v != 3 {
			t.Errorf("vertex %v counts wrong, got %d, want 3", k, v)
		}
	}
}
