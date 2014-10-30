package exactfloat

import (
	"math"
	"testing"
)

func TestSign(t *testing.T) {
	tests := []struct {
		v    float64
		want int
	}{
		{0, 1},
		{1, 1},
		{2, 1},
		{-1, -1},
		{1.2345, 1},
		{-1.2345, -1},
	}
	for _, test := range tests {
		f := NewExactFloat(test.v)
		if f.sign != test.want {
			t.Errorf("got %v, want %v", f.sign, test.want)
		}
	}
}

func TestSignedZeroAndInfinity(t *testing.T) {
	tests := []struct {
		f    ExactFloat
		want float64
	}{
		{SignedZero(1), math.Copysign(0, 1)},
		{SignedZero(-1), math.Copysign(0, -1)},
		{Infinity(1), math.Inf(1)},
		{Infinity(-1), math.Inf(-1)},
	}
	for _, test := range tests {
		var wantsign int
		if math.Signbit(test.want) {
			wantsign = -1
		} else {
			wantsign = 1
		}
		if test.f.sign != wantsign {
			t.Errorf("sign %v: got %d, want %d", test.f, test.f.sign, wantsign)
		}
	}
}

func TestToDouble(t *testing.T) {
	tests := []struct {
		f    ExactFloat
		want float64
	}{
		{NewExactFloat(0.0), 0.0},
		{NewExactFloat(math.Copysign(0, -1)), math.Copysign(0, -1)},
		{NewExactFloat(1.0), 1.0},
		{NewExactFloat(-1.0), -1.0},
		{NewExactFloat(2.5), 2.5},
		{NewExactFloat(-2.5), -2.5},
		{NewExactFloat(math.SmallestNonzeroFloat64), math.SmallestNonzeroFloat64},
		{NewExactFloat(math.MaxFloat64), math.MaxFloat64},
		{NewExactFloat(12345.6789), 12345.6789},
		{NewExactFloat(-12345.6789), -12345.6789},
	}
	for _, test := range tests {
		got := test.f.ToDouble()
		if got != test.want {
			t.Errorf("%v.ToDouble() = %v, want %v", test.f, got, test.want)
		}
	}
}

func TestAdd(t *testing.T) {
	tests := []struct {
		a    ExactFloat
		b    ExactFloat
		want float64
	}{
		{NewExactFloat(0), NewExactFloat(0), 0},
		{NewExactFloat(1), NewExactFloat(-1), 0},
		{NewExactFloat(5), NewExactFloat(5), 10},
		{NewExactFloat(1.25), NewExactFloat(1.25), 2.5},
		{NewExactFloat(1.25), NewExactFloat(-0.25), 1.0},
		{NewExactFloat(math.MaxFloat64), NewExactFloat(-math.MaxFloat64), 0},
	}
	for _, test := range tests {
		got := test.a.Add(test.b).ToDouble()
		if math.Abs(got-test.want) > 1e-14 {
			t.Errorf("%v.Add(%v) = %v, want %v", test.a, test.b, got, test.want)
		}
	}
}

func TestSub(t *testing.T) {
	tests := []struct {
		a    ExactFloat
		b    ExactFloat
		want float64
	}{
		{NewExactFloat(0), NewExactFloat(0), 0},
		{NewExactFloat(1), NewExactFloat(2), -1},
		{NewExactFloat(1), NewExactFloat(-1), 2},
		{NewExactFloat(5), NewExactFloat(3), 2},
		{NewExactFloat(1.25), NewExactFloat(.25), 1},
		{NewExactFloat(1.25), NewExactFloat(1), .25},
		{NewExactFloat(1), NewExactFloat(0), 1},
		{NewExactFloat(math.MaxFloat64), NewExactFloat(math.MaxFloat64), 0},
	}
	for _, test := range tests {
		got := test.a.Sub(test.b).ToDouble()
		if math.Abs(got-test.want) > 1e-14 {
			t.Errorf("%v.Sub(%v) = %v, want %v", test.a, test.b, got, test.want)
		}
	}
}

func TestMul(t *testing.T) {
	tests := []struct {
		a    ExactFloat
		b    ExactFloat
		want ExactFloat
	}{
		{NewExactFloat(0), NewExactFloat(0), NewExactFloat(0)},
		{NewExactFloat(1.25), NewExactFloat(1.25), NewExactFloat(1.5625)},
		{NewExactFloat(10), NewExactFloat(10), NewExactFloat(100)},
		{NewExactFloat(-2), NewExactFloat(2), NewExactFloat(-4)},
		{NewExactFloat(.5), NewExactFloat(.5), NewExactFloat(.25)},
		{NewExactFloat(-1), NewExactFloat(-1), NewExactFloat(1)},
	}
	for _, test := range tests {
		got := test.a.Mul(test.b)
		if got.ToString() != test.want.ToString() {
			t.Errorf("%s * %s = %s, want %s", test.a.ToString(), test.b.ToString(),
				got.ToString(), test.want.ToString())
		}
	}
}

func TestLargeAddEqual(t *testing.T) {
	tests := []struct {
		a, b, want ExactFloat
	}{
		{NewExactFloat(math.MaxFloat64), NewExactFloat(math.MaxFloat64), NewExactFloat(math.MaxFloat64).Add(NewExactFloat(math.MaxFloat64))},
	}
	for _, test := range tests {
		c := test.a.Add(test.b)
		if !c.Eq(test.want) {
			t.Errorf("%v + %v = %v, want %v", test.a, test.b, c, test.want)
		}
	}
}

func TestLargeMul(t *testing.T) {
	tests := []struct {
		a, b ExactFloat
		want string
	}{
		{NewExactFloat(math.MaxFloat64), NewExactFloat(math.MaxFloat64), "3.23170060713110001248980312245796e+616"},
	}
	for _, test := range tests {
		got := test.a.Mul(test.b)
		if got.ToString() != test.want {
			t.Errorf("%s * %s = %s, want %s", test.a.ToString(), test.b.ToString(), got.ToString(), test.want)
		}
	}
}

func TestCmp(t *testing.T) {
	tests := []struct {
		a    ExactFloat
		b    ExactFloat
		want bool
	}{
		{NewExactFloat(math.NaN()), NewExactFloat(math.NaN()), false},
		{NewExactFloat(0), NewExactFloat(-0), true},
		{NewExactFloat(-1), NewExactFloat(1), false},
		{NewExactFloat(math.SmallestNonzeroFloat64), NewExactFloat(math.SmallestNonzeroFloat64), true},
		{NewExactFloat(math.MaxFloat64), NewExactFloat(math.MaxFloat64), true},
	}
	for _, test := range tests {
		got := test.a.Eq(test.b)
		if got != test.want {
			t.Errorf("%v.Eq(%v) = %v, want %v", test.a, test.b, got, test.want)
		}
	}
}

func TestToString(t *testing.T) {
	tests := []struct {
		a    ExactFloat
		want string
	}{
		{NewExactFloat(math.NaN()), "nan"},
		{NewExactFloat(math.Inf(1)), "inf"},
		{NewExactFloat(math.Inf(-1)), "-inf"},
		{NewExactFloat(0), "0"},
		{NewExactFloat(math.Copysign(0, -1)), "-0"},
		{NewExactFloat(1.0), "1"},
		{NewExactFloat(1.5), "1.5"},
		{NewExactFloat(1. / 512), "0.001953125"},
		{NewExactFloat(1.23456789), "1.2345678899999999"},
		{NewExactFloat(math.SmallestNonzeroFloat64), "4.9406564584124654e-324"},
		{NewExactFloat(math.MaxFloat64), "1.7976931348623157e+308"},
		{NewExactFloat(math.MaxFloat64).Add(NewExactFloat(math.MaxFloat64)), "3.59538626972463142e+308"},
	}
	for _, test := range tests {
		got := test.a.ToString()
		if got != test.want {
			t.Errorf("%v.ToString() = %s, want %s", test.a, got, test.want)
		}
	}
}
