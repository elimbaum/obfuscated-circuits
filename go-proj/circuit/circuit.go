package circuit

import (
	"fmt"
	"strings"
	"sync"

	"golang.org/x/exp/slices"
	"gonum.org/v1/gonum/stat/combin"
)

// A gate with three pins
type Gate struct {
	Active int
	Ctrl1  int
	Ctrl2  int
	id     int // the gate's identity in the gate library
}

func (g Gate) Id() int {
	return g.id
}

func MakeGate(a, c1, c2, id int) Gate {
	if id > 255 {
		panic(fmt.Sprintf("id %d too large!", id))
	}
	return Gate{Active: a, Ctrl1: c1, Ctrl2: c2, id: id}
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
	return string(byte(g.id))
	// return fmt.Sprintf("%d%d%d", g.Active, g.Ctrl1, g.Ctrl2)
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
// Canonicalized circuits will necessarily sort identical gates at the same topo
// level next to eachother
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

func BuildFrom(wires, gates int, store map[string]PersistPermStore) chan []int {
	ch := make(chan []int)

	B := BaseGates(wires)
	rev := make(map[string]int)

	for i, g := range B {
		rev[fmt.Sprintf("%d%d%d", g[0], g[1], g[2])] = i
	}

	go func() {
		for _, perm := range store {
			s := make([]int, gates-1)

			for _, c := range perm.Ckts {
				// reconstruct prefix circuit
				for i, g := range c {
					// fmt.Println(i, c[j:j+3])
					s[i] = int(g)
				}
				// fmt.Println("Load", c, s)

				// now create new ones
				for g := 0; g < len(B); g++ {
					// no duplicates
					if g == s[len(s)-1] {
						continue
					}

					q1 := append(s, g)
					ch <- q1

					q2 := append([]int{g}, s...)
					ch <- q2

					// TODO: also send reverse circuit here?
				}
			}
		}

		close(ch)
	}()

	return ch
}

// TODO catalan numbers/dyck trees?
func ParAllCircuits(wires, gates int) chan []int {
	ch := make(chan []int)

	B := BaseGates(wires)
	z := int64(len(B))
	total := int64(1)
	for i := 0; i < gates; i++ {
		total *= z
	}

	go func() {
		WORKERS := 1
		work_ch := make(chan int64)

		defer close(ch)

		var wg sync.WaitGroup

		fmt.Println(WORKERS, "build workers launching")

		// Represent gates as `gates`-digit, base-`z` numbers
		for w := 0; w < WORKERS; w++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				send := true
				for i := range work_ch {
					send = true
					s := make([]int, gates)
					for j := 0; j < gates; j++ {
						s[j] = int(i % z)
						i /= z

						if j > 0 && s[j] == s[j-1] {
							send = false
							break
						}
					}
					if send {
						ch <- s
					}
				}
			}()
		}

		for i := int64(0); i < total; i++ {
			work_ch <- i
		}

		close(work_ch)

		wg.Wait()
	}()

	return ch
}
