package s2

import (
	"math"
	"testing"
)

func TestDefaults(t *testing.T) {
	tests := []struct {
		metric LengthMetric
		want   float64
	}{
		{MinAngleSpan, 4. / 3},
		{MaxAngleSpan, 1.704897179199218452},
		{AvgAngleSpan, math.Pi / 2},
		{MinWidth, 2 * math.Sqrt(2) / 3},
		{MaxWidth, MaxAngleSpan.Deriv},
		{AvgWidth, 1.434523672886099389},
		{MinEdge, 2 * math.Sqrt(2) / 3},
		{MaxEdge, MaxAngleSpan.Deriv},
		{AvgEdge, 1.459213746386106062},
		{MinDiag, 8 * math.Sqrt(2) / 9},
		{MaxDiag, 2.438654594434021032},
		{AvgDiag, 2.060422738998471683},
	}
	for _, test := range tests {
		got := test.metric.Deriv
		if math.Abs(got-test.want) > 1e-14 {
			t.Errorf("%v.Deriv = %v, want %v", test.metric, got, test.want)
		}
	}
}

func TestLengthMetricValue(t *testing.T) {
	tests := []struct {
		metric LengthMetric
		level  int
		want   float64
	}{
		{NewLengthMetric(1.0), 1, 0.5},
		{NewLengthMetric(1.0), 2, 0.25},
		{NewLengthMetric(2.0), 1, 1.0},
		{NewLengthMetric(2.0), 2, 0.5},
	}
	for _, test := range tests {
		got := test.metric.Value(test.level)
		if math.Abs(got-test.want) > 1e-14 {
			t.Errorf("%v.Value(%d) = %v, want %v", test.metric, test.level, got, test.want)
		}
	}
}

func TestLengthMetricClosestLevel(t *testing.T) {
	tests := []struct {
		metric LengthMetric
		value  float64
		want   int
	}{
		{NewLengthMetric(1.0), 1.0, 0},
		{NewLengthMetric(1.0), .25, 2},
		{NewLengthMetric(2.0), 1.0, 1},
		{NewLengthMetric(2.0), 0.5, 2},
	}
	for _, test := range tests {
		got := test.metric.ClosestLevel(test.value)
		if got != test.want {
			t.Errorf("%v.ClosestLevel(%f) = %v, want %v", test.metric, test.value, got, test.want)
		}
	}
}
