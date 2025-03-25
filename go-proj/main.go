package main

import (
	"fmt"
	ckt "local-mixing/circuit"
	"math"
	"time"
)

func PermString(slice []int) string {
	// This will have perfect correctness up to 8 wires
	b := make([]byte, len(slice))
	for i, v := range slice {
		b[i] = byte(v)
	}
	return string(b)
}

type PermStore struct {
	Perm []int
	Ckts map[string]bool
}

func NewPermStore(s []int) *PermStore {
	ps := new(PermStore)
	(*ps).Perm = s
	(*ps).Ckts = make(map[string]bool)
	return ps
}

func (p *PermStore) add(repr string) {
	(*p).Ckts[repr] = true
}

func invertPerm(p []int) []int {
	inv := make([]int, len(p))

	for i, q := range p {
		inv[q] = i
	}
	return inv
}

// func checkCircuit(s string) {
// 	c := ckt.FromString(s)
// 	fmt.Println(c)
// 	p := c.Perm()
// 	ph := PermString(p)
// 	fmt.Println(p, ph)
// 	ip := invertPerm(p)
// 	iph := PermString(ip)
// 	fmt.Println(ip, iph)
// }

func main2() {
	c := ckt.FromString("0 2 1; 1 2 0; 0 2 1; 1 2 0; 2 0 1; 0 1 2; 2 1 0; 2 0 1; 1 2 0; 2 0 1")
	fmt.Println(c)
	c.Canonicalize()
	fmt.Println(c)
	fmt.Println(c.Perm())

	// fmt.Println("--\n")

	// c = ckt.FromString("0 1 2;1 3 2;1 2 0;1 0 2;")
	// // fmt.Println(c)
	// c.Canonicalize()
	// fmt.Println(c)

	// 0 2 1;1 3 2;1 3 0;1 0 2;
	// 0 2 1;1 3 0;1 3 2;1 0 2;
	// 0 2 1;1 2 0;1 2 3;1 0 3;
}

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

	Counter := make(map[string]*PermStore)
	c2p := make(map[string]string)

	total_ckt := float64(len(B)) * math.Pow(float64(len(B)-1), float64(m-1))
	fmt.Println("Circuits:", total_ckt)

	i := 0

	go func() {
		start := time.Now()

		for {
			time.Sleep(1 * time.Second)
			elapsed := time.Since(start)
			if elapsed.Seconds() < 1 {
				continue
			}
			kper_second := float64(i/int(elapsed.Seconds())) / 1000
			eta := int((total_ckt - float64(i)) / kper_second / 1000)
			fmt.Printf("@ %dM, %.1fk/s ETA: %d sec \n", i/1000000, kper_second, eta)
		}
	}()

	for cc := range ch {
		for i, g := range cc {
			C[i] = B[g]
		}

		c := ckt.MakeCircuit(C)
		c.Canonicalize()
		r := c.Repr()

		// Already seen this repr?
		if ph, ok := c2p[r]; ok {
			Counter[ph].add(r)
		} else {
			p := c.Perm()
			ph := PermString(p)

			// Already seen this perm's inverse?
			if _, ok := Counter[ph]; !ok {
				// Add the inverse to the map
				ip := invertPerm(p)
				iph := PermString(ip)
				if _, ok := Counter[ph]; !ok {
					Counter[iph] = NewPermStore(p)
				}
				Counter[iph].add(r)
				c2p[r] = iph
			}
		}
		i += 1
	}

	// CktCountCount := make(map[int]int)

	for _, pp := range Counter {
		fmt.Println("====")
		fmt.Println(pp.Perm)
		for c, _ := range pp.Ckts {
			fmt.Println(c)
		}
	}

	fmt.Println(Counter[PermString([]int{7, 6, 5, 4, 3, 2, 1, 0})])

	fmt.Println(len(c2p))
	fmt.Println(len(Counter))

	// fmt.Println("NumCktPerPerm Occurences")
	// for i, p := range CktCountCount {
	// 	fmt.Printf("%d\t%d\n", i, p)
	// }
}
