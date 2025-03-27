package main

import (
	"fmt"
	ckt "local-mixing/circuit"
	"math"
	"runtime"
	"sync"
	"time"

	"golang.org/x/exp/slices"
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
	Perm  []int
	Ckts  map[string]bool
	Count int
}

func NewPermStore(s []int) *PermStore {
	ps := new(PermStore)
	(*ps).Perm = s
	(*ps).Ckts = make(map[string]bool)
	(*ps).Count = 0
	return ps
}

func (p *PermStore) addCircuit(repr string) {
	(*p).Count += 1
	(*p).Ckts[repr] = true
}

func (p *PermStore) increment() {
	(*p).Count += 1
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
	// c := ckt.FromString("0 2 1; 1 2 0; 0 2 1; 1 2 0; 2 0 1; 0 1 2; 2 1 0; 2 0 1; 1 2 0; 2 0 1")
	// fmt.Println(c)
	// c.Canonicalize()
	// fmt.Println(c)
	// fmt.Println(c.Perm())
	// fmt.Println(c.Perm().Canonical())
	// fmt.Println(ckt.Permutation([]int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}).Canonical())
	// ckt.CanonicalPerm([]int{0, 1, 3, 2, 4, 7, 6, 5})

	c := ckt.FromStringCompressed("210230023130")
	// c := ckt.FromStringCompressed("210230")
	fmt.Println(c)
	fmt.Println(c.Perm())
	c.Canonicalize()
	fmt.Println(c)
	fmt.Println(c.Perm())

	fmt.Println("=====")

	c = ckt.FromStringCompressed("230210023130")
	fmt.Println(c)
	fmt.Println(c.Perm())
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

type PR struct {
	P         ckt.Permutation
	R         string
	Canonical bool
}

func main() {
	n := 4
	m := 5
	ch := ckt.AllCircuits(n, m)

	bp := ckt.BaseGates(n)
	B := make([]ckt.Gate, len(bp))
	for i, b := range bp {
		B[i] = ckt.MakeGate(b[0], b[1], b[2])
	}

	C := make([]ckt.Gate, m)

	// Counter := make(map[string]*PermStore)
	var CircuitStore sync.Map
	// c2p := make(map[string]string)

	total_ckt := float64(len(B)) * math.Pow(float64(len(B)-1), float64(m-1))
	fmt.Println("Circuits:", total_ckt)

	ckt_i := 0
	skip_id := 0
	skip_inv := 0

	go func() {
		start := time.Now()

		var m runtime.MemStats

		for {
			time.Sleep(1 * time.Second)
			elapsed := time.Since(start)
			if elapsed.Seconds() < 1 {
				continue
			}

			runtime.ReadMemStats(&m)
			mib_used := m.Alloc / 1024 / 1024
			kper_second := float64(ckt_i/int(elapsed.Seconds())) / 1000
			eta := int((total_ckt - float64(ckt_i)) / kper_second / 1000)
			fmt.Printf("@ %.1fM, %.1fk/s ETA: %d sec (mem: %d MiB)\n", float64(ckt_i)/1000000, kper_second, eta, mib_used)
		}
	}()

	const BATCH_SIZE = 16384
	pch := make(chan []PR)

	go func() {
		index := 0
		batch := make([]PR, BATCH_SIZE)

		for cc := range ch {

			ckt_i += 1
			for i, g := range cc {
				C[i] = B[g]
			}

			// Build the circuit...
			c := ckt.MakeCircuit(C)
			// Enforce `n` wires in case terminal wires have no gates on them
			// (in that case auto-sizing would shrink the circuit)
			c.Wires = n

			// ...and compute its canonical representation
			c.Canonicalize()
			if c.AdjacentId() {
				// Skip circuits with a trivial identity = pair of adjacent
				// identical gates
				skip_id += 1
				continue
			}

			// Compute the permutation, and its bit-shuffled canonicalization
			isCanonicalPerm := false
			var p []int
			{
				p_raw := c.Perm()
				p = p_raw.Canonical()
				isCanonicalPerm = slices.Equal(p_raw, p)
			}

			batch[index] = PR{P: p, R: c.Repr(), Canonical: isCanonicalPerm}

			index++
			if index == BATCH_SIZE {
				pch <- batch
				batch = make([]PR, BATCH_SIZE)
				index = 0
			}

			// pch <- PR{P: p, R: c.Repr(), Canonical: isCanonicalPerm}
		}
		close(pch)
	}()

	WORKERS := runtime.NumCPU() - 2
	fmt.Println(WORKERS, "workers launching")

	n_perms := 0

	var wg sync.WaitGroup

	for w := 0; w < WORKERS; w++ {
		wg.Add(1)

		go func() {
			defer wg.Done()
			for batch := range pch {
				for _, pr := range batch {
					// Check if this permutation is its own inverse
					p := pr.P

					ip := invertPerm(p)
					ownInv := slices.Equal(p, ip)

					// Generate the string repr
					ph := PermString(p)

					// Check if we've already seen this perm's (unique) inverse.
					// If so, skip
					if _, ok := CircuitStore.Load(ph); ok && !ownInv {
						skip_inv += 1
						continue
					}

					// At this point: either we haven't yet seen this perm's inverse, or it
					// is its own inverse. Add it to the record.
					iph := PermString(ip)
					_st, ok := CircuitStore.Load(iph)
					var store *PermStore
					if !ok {
						store = NewPermStore(p)
						CircuitStore.Store(iph, store)
						n_perms += 1
					} else {
						store = _st.(*PermStore)
					}

					if pr.Canonical {
						// Only store circuits for which P = Canon(P). All other circuits
						// can be easily generated from that set, and we save O(n!) space
						store.addCircuit(pr.R)
					} else {
						// Not a canonical circuit, so just bump the counter
						store.increment()
					}
				}
			}
		}()
	}

	wg.Wait()

	totalStore := 0
	CircuitStore.Range(func(k any, v any) bool {
		pp, _ := v.(*PermStore)
		totalStore += len(pp.Ckts)

		if len(pp.Ckts) <= 1 {
			return true
		}

		fmt.Println("====")
		fmt.Printf("%v %v total; %v canon\n", pp.Perm, pp.Count, len(pp.Ckts))

		for c, _ := range pp.Ckts {
			fmt.Println(c)
			// fmt.Println(ckt.FromStringCompressed(c))
		}
		return true
	})

	fmt.Println("===========================")
	fmt.Printf("Total n=%d,m=%d circuits: %d\n", n, m, int(total_ckt))

	fmt.Println("  Circuits stored:", totalStore)

	// fmt.Println(Counter[PermString([]int{7, 6, 5, 4, 3, 2, 1, 0})])

	// fmt.Println(len(c2p))
	fmt.Println("  Canonical perms:", n_perms)

	fmt.Println("  Skipped b/c trivial Id:", skip_id)
	fmt.Println("  Skipped b/c saw inverse:", skip_inv)

	// fmt.Println("NumCktPerPerm Occurences")
	// for i, p := range CktCountCount {
	// 	fmt.Printf("%d\t%d\n", i, p)
	// }
}
