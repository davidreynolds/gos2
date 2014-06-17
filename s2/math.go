package s2

import "math"

const (
	M_1_PI = 1. / math.Pi
	M_PI_2 = math.Pi / 2
	M_PI_4 = math.Pi / 4
)

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
