package circuit

import (
	"fmt"
	"strings"

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
	return fmt.Sprintf("%d %d %d", g.Active, g.Ctrl1, g.Ctrl2)
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

func (c Circuit) Perm() []int {
	output := make([]int, 1<<c.Wires)

	for i := range output {
		output[i] = c.eval(i)
	}

	return output
}

func (c Circuit) Repr() string {
	var sb strings.Builder
	for _, g := range c.Gates {
		sb.WriteString(g.Repr() + ";")
	}

	return sb.String()
}

func adjacentRepeats(x []int) bool {
	for i := 0; i < len(x)-1; i++ {
		if x[i] == x[i+1] {
			return true
		}
	}
	return false
}

func BasePerms(n int) [][]int {
	return combin.Permutations(n, 3)
}

func AllCircuits(wires, gates int) chan []int {
	// Represent circuits as indices into an array of gates (represented as an
	// array of pins) -- minimize object construction for memory efficiency.
	B := BasePerms(wires)

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
