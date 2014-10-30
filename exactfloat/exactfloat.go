package exactfloat

import (
	"fmt"
	"math"
	"math/big"
	"strings"
)

const (
	maxExp             = 200 * 1000 * 1000
	minExp             = -maxExp
	maxPrec            = 64 << 20
	expNaN             = math.MaxInt32
	expInfinity        = math.MaxInt32 - 1
	expZero            = math.MaxInt32 - 2
	doubleMantissaBits = 53
)

const (
	roundTiesToEven = iota
	roundTiesAwayFromZero
	roundTowardZero
	roundAwayFromZero
	roundTowardPositive
	roundTowardNegative
)

type ExactFloat struct {
	sign   int
	bn_exp int
	bn     *big.Int
}

func NewExactFloat(v float64) ExactFloat {
	f := ExactFloat{bn: big.NewInt(0)}
	sb := math.Signbit(v)
	if sb {
		f.sign = -1
	} else {
		f.sign = 1
	}
	if math.IsNaN(v) {
		f.set_nan()
	} else if math.IsInf(v, int(f.sign)) {
		f.set_inf(f.sign)
	} else {
		frac, exp := math.Frexp(math.Abs(v))
		m := uint64(math.Ldexp(frac, doubleMantissaBits))
		f.bn = f.bn.SetUint64(m)
		f.bn_exp = exp - doubleMantissaBits
		f.Canonicalize()
	}
	return f
}

func Abs(a ExactFloat) ExactFloat {
	return a.CopyWithSign(+1)
}

func SignedZero(sign int) ExactFloat {
	f := NewExactFloat(math.Copysign(0, float64(sign)))
	f.set_zero(sign)
	return f
}

func Infinity(sign int) ExactFloat {
	f := NewExactFloat(math.Inf(sign))
	f.set_inf(sign)
	return f
}

func NaN() ExactFloat {
	f := NewExactFloat(0)
	f.set_nan()
	return f
}

func (a ExactFloat) LessThan(b ExactFloat) bool {
	// NaN is unordered compared to everything, including itself.
	if a.IsNaN() || b.IsNaN() {
		return false
	}
	// Positive and negative zero are equal.
	if a.is_zero() && b.is_zero() {
		return false
	}
	// Otherwise, anything negative is less than anything positive.
	if a.sign != b.sign {
		return a.sign < b.sign
	}
	// Now we just compare absolute values.
	if a.sign > 0 {
		return a.UnsignedLess(b)
	}
	return b.UnsignedLess(a)
}

func (a ExactFloat) UnsignedLess(b ExactFloat) bool {
	// Handle the zero/infinity cases (NaN has already been done).
	if a.is_inf() || b.is_zero() {
		return false
	}
	if a.is_zero() || b.is_inf() {
		return true
	}
	// If the high-order bit positions differ, we are done.
	cmp := a.exp() - b.exp()
	if cmp != 0 {
		return cmp < 0
	}
	// Otherwise shift one of the two values so that they both have
	// the same bn_exp and then compare the mantissas.
	if a.bn_exp >= b.bn_exp {
		return a.ScaleAndCompare(b) < 0
	}
	return b.ScaleAndCompare(a) > 0
}

func (a ExactFloat) ScaleAndCompare(b ExactFloat) int {
	tmp := a
	tmp.bn = tmp.bn.Lsh(tmp.bn, uint(a.bn_exp-b.bn_exp))
	return tmp.bn.Cmp(b.bn)
}

func (f *ExactFloat) Canonicalize() {
	if !f.is_normal() {
		return
	}

	// Underflow/overflow occurs if exp() is not in [MinExp, MaxExp].
	// We also convert a zero mantissa to signed zero.
	my_exp := f.exp()
	if my_exp < minExp || f.bn.BitLen() == 0 {
		f.set_zero(f.sign)
	} else if my_exp > maxExp {
		f.set_inf(f.sign)
	} else if f.bn.BitLen() > 0 && f.bn.Bit(0) != 0 {
		shift := count_low_zero_bits(f.bn)
		if shift > 0 {
			f.bn_exp += shift
		}
	}
	if f.prec() > maxPrec {
		f.set_nan()
	}
}

func (f ExactFloat) Eq(b ExactFloat) bool {
	// NaN is not equal to anything, not even itself.
	if f.is_nan() || b.is_nan() {
		return false
	}
	// Since Canonicalize() strips low-order zero bits, all other cases
	// (including non-normal values) require bn_exp to be equal.
	if f.bn_exp != b.bn_exp {
		return false
	}
	// Positive and negative zero are equal.
	if f.is_zero() && b.is_zero() {
		return true
	}
	// Otherwise, the signs and mantissas must match. Note that non-normal
	// values such as infinity have a mantissa of zero.
	if f.sign != b.sign {
		return false
	}
	return f.bn.Cmp(b.bn) == 0
}

func (f ExactFloat) Add(b ExactFloat) ExactFloat {
	return SignedSum(f.sign, &f, b.sign, &b)
}

func (f ExactFloat) Sub(b ExactFloat) ExactFloat {
	return SignedSum(f.sign, &f, -b.sign, &b)
}

func (f ExactFloat) Mul(b ExactFloat) ExactFloat {
	result_sign := f.sign * b.sign
	if !f.is_normal() || !b.is_normal() {
		// Handle zero, inf, and NaN according to IEEE 754-2008.
		if f.is_nan() {
			return f
		}
		if b.is_nan() {
			return b
		}
		if f.is_inf() {
			// Infinity times zero yields NaN.
			if b.is_zero() {
				return NaN()
			}
			return Infinity(result_sign)
		}
		if b.is_inf() {
			if f.is_zero() {
				return NaN()
			}
			return Infinity(result_sign)
		}
		return SignedZero(result_sign)
	}
	r := NewExactFloat(0)
	r.sign = result_sign
	r.bn_exp = f.bn_exp + b.bn_exp
	r.bn = r.bn.Mul(f.bn, b.bn)
	r.Canonicalize()
	return r
}

func SignedSum(a_sign int, a *ExactFloat, b_sign int, b *ExactFloat) ExactFloat {
	if !a.is_normal() || !b.is_normal() {
		// Handle zero, inf, and NaN according to IEEE 754-2008.
		if a.is_nan() {
			return *a
		}
		if b.is_nan() {
			return *b
		}
		if a.is_inf() {
			// Adding two infinities with opposite signs yields NaN.
			if b.is_inf() && a_sign != b_sign {
				return NaN()
			}
			return Infinity(a_sign)
		}
		if b.is_inf() {
			return Infinity(b_sign)
		}
		if a.is_zero() {
			if !b.is_zero() {
				return b.CopyWithSign(b_sign)
			}
			// Adding two zeros with the same sign preserves the sign
			if a_sign == b_sign {
				return SignedZero(a_sign)
			}
			return SignedZero(+1)
		}
		return a.CopyWithSign(a_sign)
	}
	// Swap the numbers if necessary so that "a" has the larger bn_exp.
	if a.bn_exp < b.bn_exp {
		a_sign, b_sign = b_sign, a_sign
		a, b = b, a
	}
	// Shift "a" if necessary so that both values have the same bn_exp.
	r := NewExactFloat(0)
	if a.bn_exp > b.bn_exp {
		r.bn = r.bn.Lsh(a.bn, uint(a.bn_exp-b.bn_exp))
		a = &r // The only field of "a" used below is bn.
	}
	r.bn_exp = b.bn_exp
	if a_sign == b_sign {
		r.bn = r.bn.Add(a.bn, b.bn)
		r.sign = a_sign
	} else {
		r.bn = r.bn.Sub(a.bn, b.bn)
		if r.bn.BitLen() == 0 {
			r.sign = +1
		} else if r.bn.Sign() == -1 {
			// The magnitude of "b" was larger.
			r.sign = b_sign
			r.bn = r.bn.Mul(r.bn, big.NewInt(-1))
		} else {
			// The were equal, or the magnitude of "a" was larger.
			r.sign = a_sign
		}
	}
	r.Canonicalize()
	return r
}

func (f ExactFloat) ToDouble() float64 {
	if f.prec() <= doubleMantissaBits {
		return f.ToDoubleHelper()
	} else {
		r := f.RoundToMaxPrec(doubleMantissaBits, roundTiesToEven)
		return r.ToDoubleHelper()
	}
}

func (f ExactFloat) ToDoubleHelper() float64 {
	sign := float64(f.sign)
	if !f.is_normal() {
		if f.is_zero() {
			return math.Copysign(0, sign)
		}
		if f.is_inf() {
			return math.Inf(f.sign)
		}
		return math.Copysign(math.NaN(), sign)
	}
	mantissa := f.bn.Uint64()
	return sign * math.Ldexp(float64(mantissa), f.bn_exp)
}

func (f ExactFloat) NumSignificantDigitsForPrec(prec int) int {
	// The simplest bound is:
	//
	//    d <= 1 + ceil(prec * log10(2))
	//
	// The following bound is tighter by 0.5 digits on average, but requires
	// the exponent to be known as well:
	//
	//    d <= ceil(exp * log10(2)) - floor((exp - prec) * log10(2))
	//
	// Since either of these bounds can be too large by 0, 1, or 2 digits,
	// we stick with the simpler first bound.
	return 1 + int(math.Ceil(float64(prec)*(math.Ln2/math.Ln10)))
}

// Numbers are always formatted with at least this many significant digits.
// This prevents small integers from being formatted in exponential notation
// (e.g. 1024 formatted as 1e+03), and also avoids the confusion of having
// supposedly "high precision" numbers formatted with just 1 or 2 digits
// (e.g. 1/512 == 0.001953125 formatted as 0.002).
const minSignificantDigits = 10

func (f ExactFloat) ToString() string {
	max_digits := int(math.Max(minSignificantDigits, float64(f.NumSignificantDigitsForPrec(f.prec()))))
	return f.ToStringWithMaxDigits(max_digits)
}

func (f ExactFloat) ToStringWithMaxDigits(max_digits int) string {
	var str string
	if !f.is_normal() {
		if f.is_nan() {
			return "nan"
		}
		if f.is_zero() {
			if f.sign < 0 {
				return "-0"
			} else {
				return "0"
			}
		}
		if f.sign < 0 {
			return "-inf"
		} else {
			return "inf"
		}
	}
	digits, exp10 := f.GetDecimalDigits(max_digits)
	if f.sign < 0 {
		str += "-"
	}

	// We use the standard '%g' formatting rules. If the exponent is less
	// than -4 or greater than or equal to the requested precision
	// (i.e., max_digits) then we use exponential notation.
	//
	// But since "exp10" is the base-10 exponent corresponding to a
	// mantissa in the range [0.1, 1), whereas the '%g' rules assum a
	// mantissa in the range [1.0, 10), we need to adjust these parameters
	// by 1.
	if exp10 <= -4 || exp10 > max_digits {
		// Use exponential format.
		str += digits[0:1]
		if len(digits) > 1 {
			str += "."
			str += digits[1:]
		}
		str += fmt.Sprintf("e%+02d", exp10-1)
	} else {
		// Use fixed format. We split this into two cases depending on
		// whether the integer portion is non-zero or not.
		if exp10 > 0 {
			if exp10 >= len(digits) {
				str += digits
				for i := exp10 - len(digits); i > 0; i-- {
					str += "0"
				}
			} else {
				str += digits[:exp10]
				str += "."
				str += digits[exp10:]
			}
		} else {
			str += "0."
			for i := exp10; i < 0; i++ {
				str += "0"
			}
			str += digits
		}
	}
	return str
}

func (f ExactFloat) GetDecimalDigits(max_digits int) (string, int) {
	// Convert the value to the form (bn * (10 ** bn_exp10)) where "bn"
	// is a positive integer (big.Int).
	var digits string
	var bn_exp10 int
	bn := big.NewInt(0)
	if f.bn_exp >= 0 {
		// The easy case: bn = f.bn * (2 ** f.bn_exp), bn_exp10 = 0.
		bn = bn.Lsh(f.bn, uint(f.bn_exp))
		bn_exp10 = 0
	} else {
		// Set bn = f.bn * (5 ** -f.bn_exp) and bn_exp10 = f.bn_exp.
		// This is equivalent to the original value of
		// (f.bn * (2 ** f.bn_exp)).
		power := big.NewInt(int64(-f.bn_exp))
		mul := big.NewInt(5)
		bn = bn.Mul(f.bn, mul.Exp(mul, power, nil))
		bn_exp10 = f.bn_exp
	}
	// Now convert "bn" to a decimal string.
	all_digits := bn.String()
	num_digits := len(all_digits)
	if num_digits <= max_digits {
		digits = all_digits
	} else {
		digits = all_digits[:max_digits]
		// Standard "printf" formatting rounds ties to an even number.
		// This means that we round up (away from zero) if highest
		// discarded digit is '5' or more, unless all other discarded
		// digits are zero, in which case we round up only if the lowest
		// kept digit is odd.
		odd := all_digits[max_digits-1] & 1
		idx := strings.IndexAny(all_digits[max_digits+1:], "123456789")
		if all_digits[max_digits] >= '5' && (odd == 1 || idx != -1) {
			digits = incrementDecimalDigits(digits)
		}
		bn_exp10 += num_digits - max_digits
	}

	// Strip any trailing zeros.
	end := len(digits)
	digits = strings.TrimRight(digits, "0")
	if len(digits) < end {
		bn_exp10 += end - len(digits)
	}
	return digits, bn_exp10 + len(digits)
}

// Increment an unsigned integer represented as a string of ASCII digits.
func incrementDecimalDigits(digits string) string {
	dbytes := []byte(digits)
	for pos := len(digits) - 1; pos >= 0; pos-- {
		if dbytes[pos] < '9' {
			dbytes[pos]++
			return string(dbytes)
		}
		dbytes[pos] = '0'
	}
	digits = "1" + string(dbytes)
	return digits
}

func (f ExactFloat) CopyWithSign(sign int) ExactFloat {
	r := f
	r.sign = sign
	return r
}

func (f ExactFloat) RoundToMaxPrec(max_prec, mode int) ExactFloat {
	// The "roundTiesToEven" mode requires at least 2 bits of precision
	// (otherwise both adjacent representable values may be odd).

	// The following test also catches zero, inf, and NaN.
	shift := f.prec() - max_prec
	if shift <= 0 {
		return f
	}
	// Round by removing the appropriate number of bits from the mantissa.
	// Note that if the value is rounded up to a power of 2, the high-order
	// bit position may increase, but in that case Canonicalize() will
	// remove at least one zero bit and so the output will still have
	// prec() <= max_prec.
	return f.RoundToPowerOf2(f.bn_exp+shift, mode)
}

func (f ExactFloat) RoundToPowerOf2(bit_exp, mode int) ExactFloat {
	// If the exponent is already large enough, or the value is zero, inf,
	// or NaN, then there is nothing to do.
	shift := bit_exp - f.bn_exp
	if shift <= 0 {
		return f
	}

	// Convert rounding up/down to toward/away from zero, so that we
	// don't need to consider the sign of the number from this point onward.
	if mode == roundTowardPositive {
		if f.sign > 0 {
			mode = roundAwayFromZero
		} else {
			mode = roundTowardZero
		}
	} else if mode == roundTowardNegative {
		if f.sign > 0 {
			mode = roundTowardZero
		} else {
			mode = roundAwayFromZero
		}
	}

	// Rounding consists of right-shifting the mantissa by "shift", and then
	// possibly incrementing the result (depending on the rounding mode, the
	// bits that were discarded, and sometimes the lowest kept bit). The
	// following code figures out whether we need to increment.
	r := NewExactFloat(0)
	increment := false
	if mode == roundTowardZero {
		// Never increment.
	} else if mode == roundTiesAwayFromZero {
		// Increment if the highest discarded bit is 1.
		if f.bn.Bit(shift-1) != 0 {
			increment = true
		}
	} else if mode == roundAwayFromZero {
		// Increment unless all discarded bits are zero.
		if count_low_zero_bits(f.bn) < shift {
			increment = true
		}
	} else {
		// Let "w/xyz" denote a mantissa where "w" is the lowest kept
		// bit and "xyz" are the discarded bits. Then using regexp
		// notation:
		//   ./0.*    -> Don't increment (fraction < 1/2)
		//   0/10*    -> Don't increment (fraction = 1/2, kept part even)
		//   1/10*    -> Increment (fraction = 1/2, kept part odd)
		//   ./1.*1.* -> Increment (fraction > 1/2)
		if (f.bn.Bit(shift-1) != 0 && f.bn.Bit(shift) != 0) || count_low_zero_bits(f.bn) < shift-1 {
			increment = true
		}
	}
	r.bn_exp = f.bn_exp + shift
	r.bn = r.bn.Rsh(f.bn, uint(shift))
	if increment {
		r.bn.Add(r.bn, big.NewInt(1))
	}
	r.sign = f.sign
	r.Canonicalize()
	return r
}

func (f *ExactFloat) set_nan() {
	f.sign = 1
	f.bn_exp = expNaN
	f.bn = f.bn.SetUint64(0)
}

func (f *ExactFloat) set_zero(sign int) {
	f.sign = sign
	f.bn_exp = expZero
	f.bn = f.bn.SetUint64(0)
}

func (f *ExactFloat) set_inf(sign int) {
	f.sign = sign
	f.bn_exp = expInfinity
	f.bn = f.bn.SetUint64(0)
}

func (f ExactFloat) prec() int {
	return f.bn.BitLen()
}

func (f ExactFloat) Prec() int    { return f.prec() }
func (f ExactFloat) IsNaN() bool  { return f.is_nan() }
func (f ExactFloat) MaxPrec() int { return maxPrec }

func (f ExactFloat) exp() int {
	return int(f.bn_exp) + f.bn.BitLen()
}

func (f ExactFloat) is_zero() bool {
	return f.bn_exp == expZero
}

func (f ExactFloat) is_inf() bool {
	return f.bn_exp == expInfinity
}

func (f ExactFloat) is_nan() bool {
	return f.bn_exp == expNaN
}

func (f ExactFloat) is_normal() bool {
	return f.bn_exp < expZero
}

func (f ExactFloat) is_finite() bool {
	return f.bn_exp <= expZero
}

func (f ExactFloat) sign_bit() bool {
	return f.sign < 0
}

// Return +1 if this ExactFloat is positive, -1 if it is negative, and 0
// if it is zero or NaN. Note that unlike sign_bit(), Sgn() returns 0 for
// both positive and negative zero.
func (f ExactFloat) Sgn() int {
	if f.is_nan() || f.is_zero() {
		return 0
	}
	return f.sign
}

// XXX: I don't like this code. I _think_ it matches BN_ext_count_low_zero_bits
// in the C++ exactfloat.cc version. Needs more testing.
func count_low_zero_bits(bn *big.Int) int {
	count := 0
	words := bn.Bits()
	for i := 0; i < len(words); i++ {
		if words[i] == 0 {
			count += 64 //8 * int(unsafe.Sizeof(&words[i]))
		} else {
			for j := 0; j < bn.BitLen(); j++ {
				if bn.Bit(j) == 0 {
					count++
				} else {
					break
				}
			}
			break
		}
	}
	return count
}
