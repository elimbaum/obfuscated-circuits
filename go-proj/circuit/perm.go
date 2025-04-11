package circuit

import (
	"fmt"
	"sync"

	"golang.org/x/exp/slices"
	"gonum.org/v1/gonum/stat/combin"
)

type PermStore struct {
	Perm           Permutation
	Ckts           sync.Map
	Count          int
	HaveAnyCircuit bool
	HaveNonCanon   bool
}

func NewPermStore(s Permutation) *PermStore {
	ps := new(PermStore)
	(*ps).Perm = s
	// (*ps).Ckts = make(map[string]bool)
	(*ps).Count = 0
	return ps
}

func (p *PermStore) AddCircuit(repr string) {
	(*p).Count += 1
	(*p).Ckts.Store(repr, true)
}

func (p *PermStore) Replace(repr string) {
	(*p).Count += 1
	(*p).Ckts.Clear()
	(*p).Ckts.Store(repr, true)
}

func (p *PermStore) Increment() {
	(*p).Count += 1
}

type Permutation []int

func (p Permutation) Invert() Permutation {
	inv := make([]int, len(p))

	for i, q := range p {
		inv[q] = i
	}
	return inv
}

func (p Permutation) Repr() string {
	// This will have perfect correctness up to 8 wires (perm will be 0-255)
	if len(p) > 256 {
		panic("perm too large!")
	}
	b := make([]byte, len(p))
	for i, v := range p {
		b[i] = byte(v)
	}
	return string(b)
}

func (p Permutation) String() string {
	return fmt.Sprintf("%v", []int(p))
}

var bit_shuf [][]int

// All possible permutations of bits over n wires
// Skip the first one, because that's the identity, and we implicitly consider
// it below by just setting `m = p` by default.
func Init(n int) {
	bit_shuf = combin.Permutations(n, n)[1:]
}

// save time by memoizing canonical perms
// note: uses a lot of memory
var memoized_canonical sync.Map

var MemoizedStored int
var Memoized int

// Find the canonical bit-shuffling of p
// TODO: only consider active wires
// TODO: consider if control pins prevent this from working
func (p Permutation) Canonical() Permutation {
	n := len(p)

	if bit_shuf == nil {
		panic("Call Init first!")
	}

	ps := p.String()
	if c, ok := memoized_canonical.Load(ps); ok {
		if c == nil {
			return p
		}
		cp := c.(Permutation)
		return cp
	}

	// store minimal bit permutation in here
	// to start, load with the current perm (unshuffled)
	m := make([]int, n)
	copy(m, p)
	// temporary to reconstruct shuffled bits
	t := make([]int, n)
	// temporary to reconstruct shuffled indices
	idx := make([]int, n)
	// temporary to shuffle t into, according to idx
	s := make([]int, n)

	// Note: it may still make sense to check the powers of two first.
	// this cuts the work in half, on average, but then (n/2)! is exponentially
	// smaller. basically same thing i implemented in python, but just running
	// at indices 2^i.

	// skip the first one, because that is the identity
	for _, r := range bit_shuf {
		for src, dst := range r {
			// map bit `src` to bit `dst`
			for i := range p {
				t[i] |= ((p[i] >> src) & 1) << dst
				idx[i] |= ((i >> src) & 1) << dst
			}
		}

		for i := range t {
			s[idx[i]] = t[i]
		}

		// lexicographical sort
		for i := range m {
			// tie - check next slot
			if s[i] == m[i] {
				continue
			}

			// new winner
			if s[i] < m[i] {
				copy(m, s)
			}

			// if s[i] != m[i], sort is over, either way
			break
		}

		// clear slices out for the next round
		clear(t)
		clear(idx)
	}

	pm := Permutation(m)
	if slices.Equal(m, p) {
		// already canonical. don't store perm
		memoized_canonical.Store(ps, nil)
	} else {
		memoized_canonical.Store(ps, pm)
		MemoizedStored += 1
	}

	Memoized += 1

	return pm
}
