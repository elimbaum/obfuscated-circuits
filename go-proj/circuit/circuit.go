package circuit

import (
	"fmt"
	"math/bits"
	"strings"

	"golang.org/x/exp/slices"
	"gonum.org/v1/gonum/stat/combin"
)

// A gate with three pins
type Gate struct {
	Active int
	Ctrl1  int
	Ctrl2  int
}

func MakeGate(a, c1, c2 int) Gate {
	return Gate{Active: a, Ctrl1: c1, Ctrl2: c2}
}

func (g Gate) top() int {
	return min(g.Active, g.Ctrl1, g.Ctrl2)
}

func (g Gate) bottom() int {
	return max(g.Active, g.Ctrl1, g.Ctrl2)
}

// NIMPLY universal control function (rule 57 / 99)
func (g Gate) eval(input int) int {
	// hack XOR
	c1 := (input >> g.Ctrl1) & 1
	c2 := (input >> g.Ctrl2) & 1
	input ^= (c1 | ((^c2) & 1)) << g.Active
	return input
}

func (g Gate) Repr() string {
	// return fmt.Sprintf("%d %d %d", g.Active, g.Ctrl1, g.Ctrl2)
	return fmt.Sprintf("%d%d%d", g.Active, g.Ctrl1, g.Ctrl2)
}

func (g Gate) Collides(h Gate) bool {
	return (g.Active == h.Ctrl1 || g.Active == h.Ctrl2 ||
		h.Active == g.Ctrl1 || h.Active == g.Ctrl2)
}

func (g Gate) Ordered(h Gate) bool {
	if g.Active > h.Active {
		return false
	} else if g.Active == h.Active {
		if g.Ctrl1 > h.Ctrl1 {
			return false
		} else if g.Ctrl1 == h.Ctrl1 {
			return g.Ctrl2 < h.Ctrl2
		}
	}

	return true
}

// A collection of gates on `n` wires
type Circuit struct {
	Wires int
	Gates []Gate
}

func MakeCircuit(gates []Gate) Circuit {
	w := 0
	for _, g := range gates {
		w = max(w, g.bottom())
	}
	return Circuit{Wires: w + 1, Gates: gates}
}

func (c Circuit) String() string {
	CKT_WIRE_CHAR := "─"
	CKT_CONTROL_CHAR := "│"
	CKT_DOT_CHAR := "●"
	CKT_OPEN_DOT := "○"

	var sb strings.Builder
	for w := 0; w < c.Wires; w++ {
		sb.WriteString(fmt.Sprintf("%2d ", w) + CKT_WIRE_CHAR)
		for _, g := range c.Gates {
			switch {
			case w == g.Active:
				sb.WriteString("( )")
			case w == g.Ctrl1:
				sb.WriteString(CKT_WIRE_CHAR + CKT_DOT_CHAR + CKT_WIRE_CHAR)
			case w == g.Ctrl2:
				sb.WriteString(CKT_WIRE_CHAR + CKT_OPEN_DOT + CKT_WIRE_CHAR)
			case w > g.top() && w < g.bottom():
				sb.WriteString(CKT_WIRE_CHAR + CKT_CONTROL_CHAR + CKT_WIRE_CHAR)
			default:
				sb.WriteString(strings.Repeat(CKT_WIRE_CHAR, 3))
			}

			sb.WriteString(CKT_WIRE_CHAR)
		}
		sb.WriteString("\n    ")

		for _, g := range c.Gates {
			if w >= g.top() && w < g.bottom() {
				sb.WriteString(" " + CKT_CONTROL_CHAR + "  ")
			} else {
				sb.WriteString(strings.Repeat(" ", 4))
			}
		}

		if w < c.Wires-1 {
			sb.WriteString("\n")
		}
	}

	return sb.String()
}

func (c Circuit) eval(state int) int {
	for _, g := range c.Gates {
		state = g.eval(state)
	}

	return state
}

type Permutation []int

func (c Circuit) Perm() Permutation {
	output := make(Permutation, 1<<c.Wires)

	for i := range output {
		output[i] = c.eval(i)
	}

	return output
}

func (c Circuit) Repr() string {
	var sb strings.Builder
	for _, g := range c.Gates {
		// sb.WriteString(g.Repr() + ";")
		sb.WriteString(g.Repr())
	}

	return sb.String()
}

func FromString(s string) Circuit {
	var gates []Gate
	for _, gs := range strings.Split(s, ";") {
		if gs == "" {
			continue
		}

		var g Gate
		fmt.Sscanf(gs, "%d %d %d", &g.Active, &g.Ctrl1, &g.Ctrl2)
		gates = append(gates, g)
	}

	return MakeCircuit(gates)
}

func FromStringCompressed(s string) Circuit {
	var gates []Gate
	for i := 0; i+2 < len(s); i += 3 {

		a := int(s[i] - '0')
		c1 := int(s[i+1] - '0')
		c2 := int(s[i+2] - '0')

		gates = append(gates, Gate{Active: a, Ctrl1: c1, Ctrl2: c2})
	}
	return MakeCircuit(gates)
}

func (c Circuit) Canonicalize() {
	// insertion sort-based
	for i := 1; i < len(c.Gates); i++ {
		gi := c.Gates[i]
		j := i - 1
		to_swap := -1
		for ; j >= 0; j-- {
			if c.Gates[j].Collides(gi) {
				break
			} else if !c.Gates[j].Ordered(gi) {
				to_swap = j
			}
		}

		if to_swap >= 0 {
			// fmt.Printf("Move %d to %d\n", i, j)
			g := c.Gates[i]
			remove := append(c.Gates[:i], c.Gates[i+1:]...)
			c.Gates = slices.Insert(remove, to_swap, g)
		}
	}
}

// Does this circuit contain two adjacent gates? (and thus a trivial identity)
func (c Circuit) AdjacentId() bool {
	for i := 0; i < len(c.Gates)-1; i++ {
		if c.Gates[i] == c.Gates[i+1] {
			return true
		}
	}
	return false
}

func adjacentRepeats(x []int) bool {
	for i := 0; i < len(x)-1; i++ {
		if x[i] == x[i+1] {
			return true
		}
	}
	return false
}

func BaseGates(n int) [][]int {
	return combin.Permutations(n, 3)
}

func AllCircuits(wires, gates int) chan []int {
	// Represent circuits as indices into an array of gates (represented as an
	// array of pins) -- minimize object construction for memory efficiency.
	B := BaseGates(wires)

	z := make([]int, gates)
	for i := 0; i < gates; i++ {
		z[i] = len(B)
	}

	ch := make(chan []int)

	go func() {
		gen := combin.NewCartesianGenerator(z)
		for gen.Next() {
			p := gen.Product(nil)

			if adjacentRepeats(p) {
				continue
			}
			ch <- p
		}
		close(ch)
	}()

	return ch
}

// Find the canonical bit-shuffling of p
func (p Permutation) Canonical() []int {
	n := len(p)
	bw := bits.Len(uint(n - 1))

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

	// skip the first one, because that is the identity
	for _, r := range combin.Permutations(bw, bw)[1:] {
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
			if s[i] == m[i] {
				continue
			}

			if s[i] < m[i] {
				copy(m, s)
			}

			// if s[i] != m[i], sort is over, either way
			break
		}

		// clear slices out for the next round
		clear(t)
		clear(idx)
		// clear(s)
	}

	return m
}
