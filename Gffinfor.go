package main

import (
	"bufio"
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
GFF/GTF (.gz) summary and attributions extraction, usage:
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
version: 1.2
release: 2019-06-03
project: https://github.com/d2jvkpn/Gffinfor
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
`

var parseAttr func(string, map[string]string, int) bool

func main() {
	if len(os.Args) == 1 || os.Args[1] == "-h" || os.Args[1] == "--help" {
		fmt.Println(HELP)
		os.Exit(2)
	}

	if strings.HasSuffix(os.Args[1], ".gtf") ||
		strings.HasSuffix(os.Args[1], ".gtf.gz") {
		log.Printf("parsing GTF file: %s\n", os.Args[1])
		parseAttr = gtfattr
	} else {
		log.Println("parsing GFF file: %s\n", os.Args[1])
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
func gtfattr(s string, kv map[string]string, line int) bool {
	tmp := make(map[string]string)
	var msg string

	s = strings.TrimRight(strings.Trim(s, " "), ";")
	msg = "failed to parse attribution of gtf at line %d\n"

	for _, i := range strings.Split(s, "\"; ") {
		if i == "" {
			continue
		}
		ii := strings.SplitN(i, " ", 2)

		if len(ii) != 2 {
			fmt.Fprintf(os.Stderr, msg, line)
			return false
		}

		ii[1], _ = url.QueryUnescape(ii[1])
		tmp[ii[0]] = strings.Trim(ii[1], "\"")
	}

	for k, v := range tmp {
		kv[k] = v
	}

	return true
}

func gffattr(s string, kv map[string]string, line int) bool {
	tmp := make(map[string]string)
	var msg string
	msg = "failed to parse attribution of gff at line %d\n"

	for _, i := range strings.Split(s, ";") {
		if i == "" {
			continue
		}
		ii := strings.SplitN(i, "=", 2)

		if len(ii) != 2 {
			fmt.Fprintf(os.Stderr, msg, line)
			return false
		}

		ii[1], _ = url.QueryUnescape(ii[1])

		tmp[ii[0]] = ii[1]
	}

	for k, v := range tmp {
		kv[k] = v
	}

	return true
}

//
func P1(scanner *bufio.Scanner) {
	Sequences := make(map[string]int)
	Sources := make(map[string]int)
	Types := make(map[string]int)
	var fds []string
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

		if len(fds) != 9 {
			fmt.Fprintf(os.Stderr, "invalid record at line %d\n", i)
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

		if len(fds) != 9 {
			fmt.Fprintf(os.Stderr, "invalid record at line %d\n", i)
			continue
		}

		if types[0] != "" && !HasElem(types, fds[2]) {
			continue
		}

		kv := make(map[string]string)
		if !parseAttr(fds[8], kv, i) {
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

		if len(fds) != 9 {
			fmt.Fprintf(os.Stderr, "invalid record at line %d\n", i)
			continue
		}

		if types[0] != "" && !HasElem(types, fds[2]) {
			continue
		}

		kv := make(map[string]string)
		if !parseAttr(fds[8], kv, i) {
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
