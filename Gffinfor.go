package main

import (
	"bufio"
	"errors"
	"fmt"
	gzip "github.com/klauspost/pgzip" //"compress/gzip"
	"log"
	"net/url"
	"os"
	"reflect"
	"sort"
	"strconv"
	"strings"
	"text/tabwriter"
)

const HELP = `
Summary gff/gtf (.gz) and extract attributions, usage:
    summary sequences, sources, types
    $ Gffinfor  <gff>
    note: stdin ("-") will be treated as gff format

    summary types' attributions
    $ Gffinfor  <gff>  <type1,type2...>
    note: "" for any type

    extract attributions and Dbxref (tsv format)
    $ Gffinfor  <gff>  <type1,type2...>  <attr1,attr2,dbxref1,dbxref2...>
	note: "position" for "chrom:start:end:strand", "type" for the third column


author: d2jvkpn
version: 1.0
release: 2019-04-02
project: https://github.com/d2jvkpn/Gffinfor
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

var parseAttr func(string, map[string]string) error

func main() {
	if len(os.Args) == 1 || os.Args[1] == "-h" || os.Args[1] == "--help" {
		fmt.Println(HELP)
		os.Exit(2)
	}

	if strings.HasSuffix(os.Args[1], ".gtf") ||
		strings.HasSuffix(os.Args[1], ".gtf.gz") {
		parseAttr = gtfattr
	} else {
		parseAttr = gffattr
	}

	ci, err := NewCmdInput(os.Args[1])
	if err != nil {
		log.Fatal(err)
	}
	defer ci.Close()

	switch len(os.Args) {
	case 2:
		P1(ci.Scanner)
	case 3:
		P2(ci.Scanner, strings.SplitN(os.Args[2], ",", -1))
	case 4:
		P3(ci.Scanner, strings.SplitN(os.Args[2], ",", -1),
			strings.SplitN(os.Args[3], ",", -1))
	default:
		fmt.Println(HELP)
	}
}

//
func gtfattr(s string, kv map[string]string) (err error) {
	tmp := make(map[string]string)

	s = strings.TrimRight(strings.Trim(s, " "), ";")

	for _, i := range strings.Split(s, "\"; ") {
		if i == "" {
			continue
		}
		ii := strings.SplitN(i, " ", 2)
		ii[1], _ = url.QueryUnescape(ii[1])
		if len(ii) != 2 {
			err = errors.New(fmt.Sprintf("failed to split %s", i))
			return
		}

		tmp[ii[0]] = strings.Trim(ii[1], "\"")
	}

	for k, v := range tmp {
		kv[k] = v
	}

	return
}

func gffattr(s string, kv map[string]string) (err error) {
	tmp := make(map[string]string)

	for _, i := range strings.Split(s, ";") {
		if i == "" {
			continue
		}
		ii := strings.SplitN(i, "=", 2)
		ii[1], _ = url.QueryUnescape(ii[1])
		if len(ii) != 2 {
			err = errors.New(fmt.Sprintf("failed to split %s", i))
			return
		}

		tmp[ii[0]] = ii[1]
	}

	for k, v := range tmp {
		kv[k] = v
	}

	return
}

//
func P1(scanner *bufio.Scanner) {
	Sequences := make(map[string]int)
	Sources := make(map[string]int)
	Types := make(map[string]int)
	var fds []string

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fds = strings.SplitN(line, "\t", 9)
		if fds[0] == "" {
			continue
		}
		Sequences[fds[0]]++
		Sources[fds[1]]++
		Types[fds[2]]++
	}

	var array [][]string
	array = append(array, []string{"NAME", "COUNT"})

	array = append(array,
		[]string{"Sequences", strconv.Itoa(len(Sequences))})

	var sKeys []string
	for k, _ := range Sources {
		sKeys = append(sKeys, k)
	}

	SortStrSlice(sKeys)

	for _, k := range sKeys {
		x := []string{"source: " + k, strconv.Itoa(Sources[k])}
		array = append(array, x)
	}

	var tKeys []string
	for k, _ := range Types {
		tKeys = append(tKeys, k)
	}
	SortStrSlice(tKeys)

	for _, k := range tKeys {
		x := []string{"type: " + k, strconv.Itoa(Types[k])}
		array = append(array, x)
	}

	PrintStrSlice(array)
}

//
func P2(scanner *bufio.Scanner, types []string) {
	TypeAttrs := make(map[string]map[string]int)
	var fds []string
	var err error
	var i int

	for scanner.Scan() {
		i++
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fds = strings.SplitN(line, "\t", 9)
		if fds[0] == "" {
			continue
		}
		if types[0] != "" && !HasElem(types, fds[2]) {
			continue
		}

		kv := make(map[string]string)
		err = parseAttr(fds[8], kv)

		if err != nil {
			log.Printf("failed to parse attributions at line %d:"+
				"\n    %s\n    %s\n\n", i, err, fds[8])

			continue
		}

		for k, v := range kv {
			nk := fds[2] + "\t" + k

			if _, ok := TypeAttrs[nk]; !ok {
				TypeAttrs[nk] = make(map[string]int)
			} else {
				TypeAttrs[nk][v]++
			}
		}
	}

	var array [][]string
	array = append(array, []string{"TYPE\tATTRIBUTION", "TOTAL", "UNIQUE"})

	var keys []string
	for k, _ := range TypeAttrs {
		keys = append(keys, k)
	}
	SortStrSlice(keys)

	for _, v := range keys {
		u := 0
		for _, c := range TypeAttrs[v] {
			u += c
		}

		array = append(array,
			[]string{v, strconv.Itoa(u), strconv.Itoa(len(TypeAttrs[v]))})
	}

	PrintStrSlice(array)
}

//
func P3(scanner *bufio.Scanner, types, attrs []string) {
	var fds []string
	fmt.Println(strings.Join(attrs, "\t"))
	var err error
	var i int

	for scanner.Scan() {
		i++
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fds = strings.SplitN(line, "\t", 9)
		if fds[0] == "" {
			continue
		}
		if types[0] != "" && !HasElem(types, fds[2]) {
			continue
		}

		kv := make(map[string]string)
		err = parseAttr(fds[8], kv)

		if err != nil {
			log.Printf("failed to parse attributions at line %d:"+
				"\n    %s\n    %s\n\n", i, err, fds[8])

			continue
		}

		for _, d := range strings.Split(kv["Dbxref"], ",") {
			if d == "" {
				continue
			}
			x := strings.SplitN(d, ":", 2)
			kv[x[0]] = x[1]
		}

		kv["position"] = strings.Join(
			[]string{fds[0], fds[3], fds[4], fds[6]}, ":")

		kv["type"] = fds[2]

		values := []string{}
		for _, k := range attrs {
			values = append(values, kv[k])
		}

		fmt.Println(strings.Join(values, "\t"))
	}
}

func HasElem(s interface{}, elem interface{}) bool {
	arrV := reflect.ValueOf(s)

	if arrV.Kind() == reflect.Slice {
		for i := 0; i < arrV.Len(); i++ {
			// XXX - panics if slice element points to an unexported struct field
			// see https://golang.org/pkg/reflect/#Value.Interface
			if arrV.Index(i).Interface() == elem {
				return true
			}
		}
	}

	return false
}

func SortStrSlice(s []string) {
	sort.Slice(s, func(i, j int) bool {
		return strings.ToLower(s[i]) < strings.ToLower(s[j])
	})
}

func PrintStrSlice(array [][]string) {
	w := tabwriter.NewWriter(os.Stdout, 4, 0, 4, ' ', tabwriter.StripEscape)
	for _, r := range array {
		fmt.Fprintln(w, "    "+strings.Join(r, "\t"))
	}
	w.Flush()

}

type CmdInput struct {
	Name    string
	File    *os.File
	Reader  *gzip.Reader
	Scanner *bufio.Scanner
}

func (ci *CmdInput) Close() {
	if ci.Reader != nil {
		ci.Reader.Close()
	}

	if ci.File != nil {
		ci.File.Close()
	}
}

func NewCmdInput(name string) (ci *CmdInput, err error) {
	ci = new(CmdInput)
	ci.Name = name

	if ci.Name == "-" {
		ci.Scanner = bufio.NewScanner(os.Stdin)
		return
	}

	ci.File, err = os.Open(ci.Name)

	if err != nil {
		return
	}

	if strings.HasSuffix(ci.Name, ".gz") {
		if ci.Reader, err = gzip.NewReader(ci.File); err != nil {
			return
		}
		ci.Scanner = bufio.NewScanner(ci.Reader)
	} else {
		ci.Scanner = bufio.NewScanner(ci.File)
	}

	return
}
