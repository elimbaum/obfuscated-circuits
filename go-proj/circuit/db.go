package circuit

import (
	"encoding/gob"
	"fmt"
	"os"
)

type PersistPermStore struct {
	Perm  []int
	Ckts  []string
	Count int
}

// A PermStore that can be persisted. (no sync.Map)
type Persist struct {
	Wires int
	Gates int
	Store map[string]PersistPermStore
}

// TODO: compression?
func Save(n int, m int, store map[string]PersistPermStore) {
	f, _ := os.Create(fmt.Sprintf("db/n%dm%d.gob", n, m))
	defer f.Close()

	enc := gob.NewEncoder(f)
	persist := Persist{Wires: n, Gates: m, Store: store}
	err := enc.Encode(persist)
	if err != nil {
		panic(err)
	}
}

func Load(n int, m int, fn string) map[string]PersistPermStore {
	f, err := os.Open(fn)
	if err != nil {
		panic(err)
	}

	dec := gob.NewDecoder(f)
	var p Persist
	err = dec.Decode(&p)
	if err != nil {
		panic(err)
	}

	if p.Gates != m-1 || p.Wires != n {
		fmt.Printf("Database size does not match: load has n=%d, m=%d; requested n=%d, m=%d\n", p.Wires, p.Gates, n, m)
		os.Exit(1)
	}

	return p.Store
}
