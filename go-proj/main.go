package main

import (
	"fmt"
	ckt "local-mixing/circuit"
	"math"
	"runtime"
	"sync"
	"sync/atomic"
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
	Ckts  sync.Map
	Count int
}

func NewPermStore(s []int) *PermStore {
	ps := new(PermStore)
	(*ps).Perm = s
	// (*ps).Ckts = make(map[string]bool)
	(*ps).Count = 0
	return ps
}

func (p *PermStore) addCircuit(repr string) {
	(*p).Count += 1
	(*p).Ckts.Store(repr, true)
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

var BaseGates []ckt.Gate

var n_perms atomic.Uint64
var ckt_check atomic.Uint64
var skip_inv atomic.Uint64
var ckt_i atomic.Uint64
var skip_id atomic.Uint64

// Take as input a bunch of circuits and output a batch of permutations
func buildCircuit(wires, gates, workers int, ckt_ch <-chan []int) <-chan []PR {
	var wg sync.WaitGroup

	const BATCH_SIZE = 1024

	perm_ch := make(chan []PR)

	go func() {
		defer close(perm_ch)
		for w := 0; w < workers; w++ {
			wg.Add(1)

			// Each worker fills up a batch
			batch := make([]PR, BATCH_SIZE)
			bidx := 0

			go func() {
				defer wg.Done()

				// Temporary to hold gates
				C := make([]ckt.Gate, gates)

				for cc := range ckt_ch {
					ckt_i.Add(1)

					for i, g := range cc {
						C[i] = BaseGates[g]
					}

					// Build the circuit, and manually specify number of wires
					// This is necessary in case terminal wires are not touched by
					// any gate (they they will be auto-removed)
					c := ckt.MakeCircuit(C)
					c.Wires = wires

					// Compute circuit's canonical repr; skip if any adjacent
					// duplicates
					c.Canonicalize()
					if c.AdjacentId() {
						skip_id.Add(1)
						continue
					}

					p := c.Perm()
					cp := p.Canonical()
					isCanonicalPerm := slices.Equal(p, cp)

					// if slices.Equal(cp, []int{0, 1, 3, 2, 7, 5, 6, 4}) {
					// 	fmt.Println("====")
					// 	fmt.Println(c)
					// 	fmt.Println(p)
					// 	fmt.Println(cp, isCanonicalPerm)
					// }

					batch[bidx] = PR{P: cp, R: c.Repr(), Canonical: isCanonicalPerm}
					bidx++

					if bidx >= BATCH_SIZE {
						perm_ch <- batch
						batch = make([]PR, BATCH_SIZE)
						bidx = 0
					}
				}

				// might have some left
				perm_ch <- batch[:bidx]
			}()
		}
		wg.Wait()
	}()

	return perm_ch
}

func main() {
	n := 5
	m := 5
	ch := ckt.ParAllCircuits(n, m)

	ckt.Init(n)

	bp := ckt.BaseGates(n)
	BaseGates = make([]ckt.Gate, len(bp))
	for i, b := range bp {
		BaseGates[i] = ckt.MakeGate(b[0], b[1], b[2])
	}
	var CircuitStore sync.Map
	// c2p := make(map[string]string)

	total_ckt := float64(len(BaseGates)) * math.Pow(float64(len(BaseGates)-1), float64(m-1))
	fmt.Printf("n=%d, m=%d. circuits: %d\n", n, m, int64(total_ckt))

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
			ci := ckt_i.Load()
			kper_second := float64(ci/(uint64(elapsed.Seconds()))) / 1000
			eta := int((total_ckt - float64(ci)) / kper_second / 1000)
			fmt.Printf("@ %.1fM, %.1fk/s ETA: %d sec (mem: %d MiB)\n", float64(ci)/1000000, kper_second, eta, mib_used)
		}
	}()

	perm_ch := buildCircuit(n, m, runtime.NumCPU(), ch)

	WORKERS := 1
	fmt.Println(WORKERS, "workers launching")

	var wg sync.WaitGroup

	for w := 0; w < WORKERS; w++ {
		wg.Add(1)

		go func() {
			defer wg.Done()
			for batch := range perm_ch {
				// _ = batch
				for _, pr := range batch {
					ckt_check.Add(1)
					if pr.P == nil {
						fmt.Println("nil perm")
						continue
					}
					// Check if this permutation is its own inverse
					p := pr.P

					ip := invertPerm(p)
					ownInv := slices.Equal(p, ip)

					if !ownInv {
						// Generate the string repr
						ph := PermString(p)

						// Check if we've already seen this perm's (unique) inverse.
						// If so, skip
						if _, ok := CircuitStore.Load(ph); ok {
							skip_inv.Add(1)
							continue
						}
					}

					// At this point: either we haven't yet seen this perm's inverse, or it
					// is its own inverse. Add it to the record.
					iph := PermString(ip)
					_st, ok := CircuitStore.Load(iph)
					var store *PermStore
					if !ok {
						store = NewPermStore(p)
						CircuitStore.Store(iph, store)
						n_perms.Add(1)
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
	all := 0
	CircuitStore.Range(func(k any, v any) bool {
		pp, _ := v.(*PermStore)
		canonical := 0
		pp.Ckts.Range(func(k any, v any) bool { canonical += 1; return true })

		totalStore += canonical
		all += pp.Count

		// if len(pp.Ckts) <= 1 {
		// 	return true
		// }

		if canonical == 0 {
			fmt.Println("ERROR! No canonical")
			fmt.Println(" ", pp.Perm, pp.Count)
		}

		// fmt.Println("====")
		// fmt.Printf("%v %v total; %v canon\n", pp.Perm, pp.Count, canonical)

		// for c, _ := range pp.Ckts {
		// fmt.Println(c)
		// fmt.Println(ckt.FromStringCompressed(c))
		// }
		return true
	})

	fmt.Println("===========================")
	fmt.Printf("Total n=%d,m=%d circuits: %d\n", n, m, int(total_ckt))

	fmt.Println("  Received from circuit gen:", ckt_i.Load())
	fmt.Println("  Skipped b/c trivial Id:", skip_id.Load())
	fmt.Println("  Circuits checked:", ckt_check.Load())
	fmt.Println("  Skipped b/c saw inverse:", skip_inv.Load())
	fmt.Println("  Circuits stored:", totalStore)
	fmt.Println("  Total count:", all)

	// fmt.Println(Counter[PermString([]int{7, 6, 5, 4, 3, 2, 1, 0})])

	// fmt.Println(len(c2p))
	fmt.Println("  Canonical perms:", n_perms.Load())

	// fmt.Println("NumCktPerPerm Occurences")
	// for i, p := range CktCountCount {
	// 	fmt.Printf("%d\t%d\n", i, p)
	// }
}
