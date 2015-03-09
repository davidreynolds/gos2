package s2

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
)

/*
func r64() uint64 {
	bits := uint(math.Log2(1 + (1 << 15)))
	res := uint64(0)
	for i := uint(0); i < 64; i += bits {
		res <<= bits
		res += uint64(rand.Int63())
	}
	return res
}

func r32() uint32 {
	return uint32(r64() & ((1 << 32) - 1))
}

func skewed(maxLog int) int32 {
	base := r32() % uint32(maxLog+1)
	return int32(r32() & ((1 << base) - 1))
}
*/

func uniform(upperBound int) int {
	return int(rand.Float64() * float64(upperBound))
}

func skewed(maxLog int32) int32 {
	base := rand.Int31() % (maxLog + 1)
	return rand.Int31() & ((1 << uint(base)) - 1)
}

func randomCellIDForLevel(level int) CellID {
	face := rand.Intn(numFaces)
	pos := uint64(rand.Int63() & ((1 << (2 * maxLevel)) - 1))
	return CellIDFromFacePosLevel(face, pos, level)
}

func randomCellID() CellID {
	return randomCellIDForLevel(uniform(maxLevel + 1))
}

func randomPoint() Point {
	x := 2*rand.Float64() - 1
	y := 2*rand.Float64() - 1
	z := 2*rand.Float64() - 1
	return PointFromCoords(x, y, z)
}

func randomCap(min_area, max_area float64) Cap {
	cap_area := max_area * math.Pow(min_area/max_area, rand.Float64())
	if cap_area < min_area {
		fmt.Printf("%v < %v\n", cap_area, min_area)
	}
	if cap_area > max_area {
		fmt.Printf("%v > %v\n", cap_area, max_area)
	}
	return CapFromCenterArea(randomPoint(), cap_area)
}

func CheckCompleteCovering(t *testing.T, region Region, covering CellUnion, check_tight bool, id CellID) {
	if !id.IsValid() {
		for face := 0; face < 6; face++ {
			newId := CellIDFromFacePosLevel(face, 0, 0)
			CheckCompleteCovering(t, region, covering, check_tight, newId)
		}
		return
	}
	if !region.MayIntersect(CellFromCellID(id)) {
		// If region does not intersect id, then neither should the
		// covering.
		if check_tight {
			if covering.IntersectsCellID(id) {
				fmt.Printf("%v.IntersectsCellID(%v)\n", covering, id)
				t.Errorf("%v.IntersectsCellID(%v)", covering, id)
			}
		}
	} else if !covering.ContainsCellID(id) {
		// The region may intersect id, but we can't assert that the
		// covering intersects id because we may discover that the
		// region does not actually intersect upon further subdivision.
		// (MayIntersect is not exact.)
		if region.ContainsCell(CellFromCellID(id)) {
			t.Errorf("%v.ContainsCell(%v)", region, CellFromCellID(id))
		}
		if id.IsLeaf() {
			t.Errorf("%d.IsLeaf()", id)
		}
		for ci := id.ChildBegin(); ci != id.ChildEnd(); ci = ci.Next() {
			CheckCompleteCovering(t, region, covering, check_tight, ci)
		}
	}
}

func CheckCovering(t *testing.T, coverer RegionCoverer, region Region, covering []CellID, interior bool) {
	// Keep track of how many cells have the same coverer.minLevel ancestor.
	min_level_cells := map[CellID]int{}
	for _, id := range covering {
		level := id.Level()
		if level < coverer.minLevel {
			t.Errorf("%v < %v", level, coverer.minLevel)
		}
		if level > coverer.maxLevel {
			t.Errorf("%v > %v", level, coverer.maxLevel)
		}
		if (level-coverer.minLevel)%coverer.levelMod != 0 {
			t.Errorf("%v % %v != 0", level-coverer.minLevel, coverer.levelMod)
		}
		min_level_cells[id.Parent(coverer.minLevel)]++
	}
	if len(covering) > coverer.maxCells {
		// If the covering has more than the requested number of cells,
		// then check that the cell count cannot be reduced by using
		// the parent of some cell.
		for k, v := range min_level_cells {
			if v != 1 {
				t.Errorf("(levels = %2d - %2d, cells = %d/%d) %v => %v",
					coverer.minLevel, coverer.maxLevel, len(covering), coverer.maxCells, k, v)
			}
		}
	}

	if interior {
		for _, id := range covering {
			if !region.ContainsCell(CellFromCellID(id)) {
				t.Errorf("!%v.ContainsCell(%v)", region, CellFromCellID(id))
			}
		}
	} else {
		var cell_union CellUnion
		cell_union.Init(covering)
		CheckCompleteCovering(t, region, cell_union, true, CellID(0))
	}
}

func TestRandomCaps(t *testing.T) {
	rand.Seed(4)
	coverer := NewRegionCoverer()
	for i := 0; i < 1000; i++ {
		for {
			coverer.SetMinLevel(rand.Intn(maxLevel + 1))
			coverer.SetMaxLevel(rand.Intn(maxLevel + 1))
			if coverer.minLevel <= coverer.maxLevel {
				break
			}
		}

		coverer.SetMaxCells(int(skewed(10)))
		coverer.SetLevelMod(1 + rand.Intn(3))
		max_area := math.Min(4*math.Pi, float64(3*coverer.maxCells+1)*AverageArea(coverer.minLevel))
		s2cap := randomCap(0.1*AverageArea(maxLevel), max_area)
		covering := coverer.Covering(s2cap)
		CheckCovering(t, coverer, s2cap, covering, false)
		// Check that Covering is deterministic
		covering2 := coverer.Covering(s2cap)
		if len(covering) == len(covering2) {
			for i, v := range covering {
				if v != covering2[i] {
					t.Errorf("%v != %v", covering, covering2)
					break
				}
			}
		} else {
			t.Errorf("len(%v) != len(%v)", covering, covering2)
		}
		// Also check Denormalize(). The denormalized covering may
		// still be different and smaller than "covering" because
		// RegionCoverer does not guarantee that it will not output
		// all four children of the same parent.
		var cells CellUnion
		cells.Init(covering)
		denormalized := cells.Denormalize(coverer.minLevel, coverer.levelMod)
		CheckCovering(t, coverer, s2cap, denormalized, false)
	}
}

func TestSimpleCoverings(t *testing.T) {
	rand.Seed(4)
	coverer := NewRegionCoverer()
	coverer.SetMaxCells(math.MaxInt32)
	for i := 0; i < 1000; i++ {
		level := uniform(maxLevel + 1)
		coverer.SetMinLevel(level)
		coverer.SetMaxLevel(level)
		max_area := math.Min(4*math.Pi, 1000*AverageArea(level))
		s2cap := randomCap(0.1*AverageArea(maxLevel), max_area)
		covering := SimpleCovering(s2cap, s2cap.center, level)
		CheckCovering(t, coverer, s2cap, covering, false)
	}
}

func TestRandomCells(t *testing.T) {
	coverer := NewRegionCoverer()
	coverer.SetMaxCells(1)
	// Test random cell ids at all levels.
	for i := 0; i < 10000; i++ {
		id := randomCellID()
		covering := coverer.Covering(CellFromCellID(id))
		if len(covering) != 1 {
			t.Errorf("%v != 1", len(covering))
		}
		if covering[0] != id {
			t.Errorf("%v != %v", covering[0], id)
		}
	}
}

func benchmarkCovering(b *testing.B, maxCells, numVerts int) {
	b.StopTimer()
	coverer := NewRegionCoverer()
	coverer.SetMaxCells(maxCells)
	for i := 0; i < b.N; i++ {
		center := randomPoint()
		loop := makeRegularLoop(center, numVerts, 0.005)
		b.StartTimer()
		coverer.Covering(loop)
		b.StopTimer()
	}
}

func BenchmarkCovering0(b *testing.B) { benchmarkCovering(b, 8, 1<<10) }
func BenchmarkCovering1(b *testing.B) { benchmarkCovering(b, 8, 1<<11) }
func BenchmarkCovering2(b *testing.B) { benchmarkCovering(b, 8, 1<<12) }
func BenchmarkCovering3(b *testing.B) { benchmarkCovering(b, 8, 1<<13) }
func BenchmarkCovering4(b *testing.B) { benchmarkCovering(b, 8, 1<<14) }
func BenchmarkCovering5(b *testing.B) { benchmarkCovering(b, 8, 1<<15) }
func BenchmarkCovering6(b *testing.B) { benchmarkCovering(b, 8, 1<<16) }
func BenchmarkCovering7(b *testing.B) { benchmarkCovering(b, 8, 1<<17) }
