package main

import (
	"flag"
	"fmt"
	ckt "local-mixing/circuit"
	"math"
	"os"
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
	Perm           []int
	Ckts           sync.Map
	Count          int
	HaveAnyCircuit bool
	HaveNonCanon   bool
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

func (p *PermStore) Replace(repr string) {
	(*p).Count += 1
	(*p).Ckts.Clear()
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

type PR struct {
	P         ckt.Permutation
	R         string
	Canonical bool
}

var BaseGates []ckt.Gate

var n_perms atomic.Int64
var ckt_check atomic.Int64
var skip_inv atomic.Int64
var ckt_i atomic.Int64
var skip_id atomic.Int64
var own_inv_count atomic.Int64

// Take as input a bunch of circuits and output a batch of permutations
func buildCircuit(wires, gates, workers int, ckt_ch <-chan []int) <-chan []PR {
	var wg sync.WaitGroup

	const BATCH_SIZE = 1024

	perm_ch := make(chan []PR)

	go func() {
		defer close(perm_ch)
		fmt.Println(workers, "perm workers launching")
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
	var n, m int
	var load string

	flag.IntVar(&n, "n", 0, "number of wires")
	flag.IntVar(&m, "m", 0, "number of gates")
	flag.StringVar(&load, "load", "", "load from pregenerated circuit database")
	fresh := flag.Bool("new", false, "generate fresh circuit database")

	flag.Parse()

	if n < 3 || m < 1 {
		fmt.Println("Invalid circuit size.")
		os.Exit(1)
	}

	if (load != "") == (*fresh) {
		fmt.Println("Specify one of --load or --new but not both")
		os.Exit(1)
	}

	var total_ckt int64

	bp := ckt.BaseGates(n)
	BaseGates = make([]ckt.Gate, len(bp))
	for i, b := range bp {
		BaseGates[i] = ckt.MakeGate(b[0], b[1], b[2])
	}

	ckt.Init(n)

	var ch chan []int

	generateNew := load == ""

	n_scratch := int64(len(BaseGates)) * int64(math.Pow(float64(len(BaseGates)-1), float64(m-1)))

	if generateNew {
		total_ckt = n_scratch
		ch = ckt.ParAllCircuits(n, m)
	} else {
		fmt.Println("Loading existing database:", load)
		store := ckt.Load(n, m, load)
		// storage should be a struct and include n, m
		ch = ckt.BuildFrom(n, m, store)

		prev_count := 0
		for _, c := range store {
			prev_count += len(c.Ckts)
		}

		total_ckt = int64(prev_count) * int64(len(BaseGates))
	}

	var CircuitStore sync.Map
	// c2p := make(map[string]string)

	fmt.Printf("n=%d, m=%d. circuits: %d\n", n, m, total_ckt)

	go func() {
		start := time.Now()

		var m runtime.MemStats
		var last int64

		for {
			time.Sleep(1 * time.Second)
			elapsed := time.Since(start)
			ci := ckt_i.Load()
			if elapsed.Seconds() < 1 || ci < 10 || last == 0 {
				last = ci
				continue
			}

			runtime.ReadMemStats(&m)
			mib_used := m.Alloc / 1024 / 1024

			kper_second := float64(ci) / elapsed.Seconds() / 1000
			eta := int(float64(total_ckt-ci) / kper_second / 1000)
			now_kper := float64(ci-last) / 1000
			fmt.Printf("@ %.1fM, now %.1fk/s, avg %.1fk/s ETA: %d sec (mem: %d MiB)\n", float64(ci)/1000000, now_kper, kper_second, eta, mib_used)
			last = ci
		}
	}()

	// -2 for MBP efficiency cores
	// -1 for breathing room
	perm_ch := buildCircuit(n, m, runtime.NumCPU()-3, ch)

	WORKERS := 1
	// fmt.Println(WORKERS, "workers launching")

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
					// ownInv := slices.Equal(p, ip)

					// if !ownInv {
					// 	if generateNew {
					// 		// Generate the string repr
					// 		ph := PermString(p)

					// 		// Check if we've already seen this perm's (unique) inverse.
					// 		// If so, skip
					// 		if _, ok := CircuitStore.Load(ph); ok {
					// 			skip_inv.Add(1)
					// 			continue
					// 		}
					// 	}
					// } else {
					// 	own_inv_count.Add(1)
					// }

					// At this point: either we haven't yet seen this perm's inverse, or it
					// is its own inverse. Add it to the record.
					iph := PermString(ip)
					_st, ok := CircuitStore.Load(iph)
					var store *PermStore
					if !ok {
						store = NewPermStore(ip)
						CircuitStore.Store(iph, store)
						n_perms.Add(1)
					} else {
						store = _st.(*PermStore)
					}

					// Don't store all circuits. Only store canonical circuits,
					// or the first circuit seen.
					//
					// Canonical circuits are those for which P = Canon(P). All
					// other circuits can be easily generated from this set.
					//
					// However, if we are loading from an existing DB, we will
					// only check a small fraction of the total circuits, and
					// thus will not see a canonical circuit. So add the _first_
					// circuit seen, along with any later canonical circuits.
					//
					// This saves O(n!) space.

					if pr.Canonical || !store.HaveAnyCircuit {
						if store.HaveNonCanon {
							store.Replace(pr.R)
							store.HaveNonCanon = false
						} else {
							store.addCircuit(pr.R)
						}

						store.HaveAnyCircuit = true
						store.HaveNonCanon = !pr.Canonical
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

	saveMap := make(map[string]ckt.PersistPermStore)

	CircuitStore.Range(func(k any, v any) bool {
		pp, _ := v.(*PermStore)
		s, _ := k.(string)
		canonical := 0

		var ckts []string

		pp.Ckts.Range(func(k any, v any) bool {
			canonical += 1

			repr, _ := k.(string)
			ckts = append(ckts, repr)

			return true
		})

		totalStore += canonical
		all += pp.Count

		ip := invertPerm(pp.Perm)

		// not own inverse?
		if !slices.Equal(ip, pp.Perm) {
			// if inv_, ok := CircuitStore.Load(PermString(ip)); ok {
			// 	inv, _ := inv_.(*PermStore)
			// 	fmt.Println("WARNING: Perm and Inverse!", pp.Perm, inv.Perm)
			// }
		}

		if canonical == 0 {
			// fmt.Println("ERROR! No canonical")
			// fmt.Println(" ", pp.Perm, pp.Count)
		}

		saveMap[s] = ckt.PersistPermStore{pp.Perm, ckts, pp.Count}

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
	if !generateNew {
		fmt.Printf("  If gen from scratch: %d\n", n_scratch)
	}

	fmt.Println("  Received from circuit gen:", ckt_i.Load())
	fmt.Println("    Skipped b/c trivial Id:", skip_id.Load())
	fmt.Println("  Circuits checked:", ckt_check.Load())
	fmt.Println("    Skipped b/c saw inverse:", skip_inv.Load())
	fmt.Println("    Circuits self inverse:", own_inv_count.Load())
	fmt.Println("  Circuits stored:", totalStore)
	fmt.Println("  Total count:", all)
	fmt.Println("  Canonical perms:", n_perms.Load())

	fmt.Println("===========================")

	// fmt.Println(Counter[PermString([]int{7, 6, 5, 4, 3, 2, 1, 0})])

	// fmt.Println(len(c2p))

	ckt.Save(n, m, saveMap)

	// fmt.Println("NumCktPerPerm Occurences")
	// for i, p := range CktCountCount {
	// 	fmt.Printf("%d\t%d\n", i, p)
	// }
}
