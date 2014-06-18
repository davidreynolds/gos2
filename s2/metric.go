package s2

import (
	"math"
)

type Metric struct {
	Deriv float64
	dim   int
}

type LengthMetric struct {
	Metric
}

type AreaMetric struct {
	Metric
}

func NewMetric(deriv float64, dim int) Metric    { return Metric{deriv, dim} }
func NewLengthMetric(deriv float64) LengthMetric { return LengthMetric{Metric{deriv, 1}} }
func NewAreaMetric(deriv float64) AreaMetric     { return AreaMetric{Metric{deriv, 2}} }

// Return the length on the unit sphere for cells at the given level.
func (m Metric) Value(level int) float64 {
	return math.Ldexp(m.Deriv, -m.dim*level)
}

// Return the level at which the length is approximately the given value.
// The return value is always a valid level.
func (m Metric) ClosestLevel(value float64) int {
	var val float64
	if m.dim == 1 {
		val = math.Sqrt2 * value
	} else {
		val = 2 * value
	}
	return m.MinLevel(val)
}

// Return the minimum level such that the length is at most the given value,
// or maxLevel if there is no such level. The return value is always a valid
// level.
func (m Metric) MinLevel(value float64) int {
	if value <= 0 {
		return maxLevel
	}
	// This code is equivalent to computing a floating-point "level"
	// value and rounding up. Frexp() returns a fraction in the range
	// [0.5, 1) and the corresponding exponent.
	_, level := math.Frexp(value / m.Deriv)
	level = max(0, min(maxLevel, -((level-1)>>uint(m.dim-1))))
	return level
}

// Return the maximum level such that the length is at least the given value,
// or zero if there is no such level. The return value is always a valid level.
func (m Metric) MaxLevel(value float64) int {
	if value <= 0 {
		return maxLevel
	}
	_, level := math.Frexp(m.Deriv / value)
	level = max(0, min(maxLevel, (level-1)>>uint(m.dim-1)))
	return level
}

var (
	MinAngleSpan  LengthMetric
	MaxAngleSpan  LengthMetric
	AvgAngleSpan  LengthMetric
	MinWidth      LengthMetric
	MaxWidth      LengthMetric
	AvgWidth      LengthMetric
	MinEdge       LengthMetric
	MaxEdge       LengthMetric
	AvgEdge       LengthMetric
	MinDiag       LengthMetric
	MaxDiag       LengthMetric
	AvgDiag       LengthMetric
	MinArea       AreaMetric
	MaxArea       AreaMetric
	AvgArea       AreaMetric
	MaxEdgeAspect float64
	MaxDiagAspect float64
)

// All of the values below were obtained by a combination of hand analysis and
// Mathematica.
// Note that I've only implemented the quadratic projections since stToUV uses
// the quadratic transform.
func init() {
	MinAngleSpan = NewLengthMetric(4. / 3)               // 1.333
	MaxAngleSpan = NewLengthMetric(1.704897179199218452) // 1.705
	AvgAngleSpan = NewLengthMetric(math.Pi / 2)          // 1.571 (true for all projections)
	MinWidth = NewLengthMetric(2 * math.Sqrt(2) / 3)     // 0.943
	MaxWidth = NewLengthMetric(MaxAngleSpan.Deriv)       // (true for all projections)
	AvgWidth = NewLengthMetric(1.434523672886099389)     // 1.435
	MinEdge = NewLengthMetric(2 * math.Sqrt(2) / 3)      // 0.943
	MaxEdge = NewLengthMetric(MaxAngleSpan.Deriv)        // (true for all projections)
	AvgEdge = NewLengthMetric(1.459213746386106062)      // 1.459
	MinDiag = NewLengthMetric(8 * math.Sqrt(2) / 9)      // 1.257
	MaxDiag = NewLengthMetric(2.438654594434021032)      // 2.439
	AvgDiag = NewLengthMetric(2.060422738998471683)      // 2.060
	MinArea = NewAreaMetric(8 * math.Sqrt(2) / 9)        // 1.257
	MaxArea = NewAreaMetric(2.635799256963161491)        // 2.636
	AvgArea = NewAreaMetric(4 * math.Pi / 6)             // 2.094 (true for all projections)
	MaxEdgeAspect = 1.442615274452682920                 // 1.443
	MaxDiagAspect = math.Sqrt(3)                         // 1.732 (true for all projections)
}
