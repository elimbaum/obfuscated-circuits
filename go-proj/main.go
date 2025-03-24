package main

import (
	"fmt"
	ckt "local-mixing/circuit"
)

func main() {
	n := 3
	m := 10
	ch := ckt.AllCircuits(n, m)

	bp := ckt.BasePerms(n)
	B := make([]ckt.Gate, len(bp))
	for i, b := range bp {
		B[i] = ckt.MakeGate(b[0], b[1], b[2])
	}

	C := make([]ckt.Gate, m)

	for cc := range ch {
		for i, g := range cc {
			C[i] = B[g]
		}
		circ := ckt.MakeCircuit(C)
		// fmt.Println(circ)
		fmt.Println(circ.Perm())
	}
}
