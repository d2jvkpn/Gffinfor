package main

import (
	"bufio"
	"fmt"
	gzip "github.com/klauspost/pgzip" //"compress/gzip"
	"log"
	"net/url"
	"os"
	"sort"
	"strconv"
	"strings"
	"text/tabwriter"
)

const HELP = `GFF/GTF(.gz) summary and attributions extraction, usage:
  sequences, sources, types summary
    $ Gffinfor  <gff|gtf>
    # stdin ("-") will be treated as gff format

  selected types' attributions statistics
    $ Gffinfor  <gff|gtf>  <type1,type2...>
    # "" for any type

  extract attributions and Dbxref (tsv format)
    $ Gffinfor  <gff|gtf>  <type1,type2...>  <attr1,attr2,dbxref1...> \
        [h1:H1,h2:H2]
    # attr "4" for position in "sequence_id:start:end:strand" format, 
    #   colname: position
    # attr "3" for the type(3rd column), colname: type
    # attr "0" for all columns(whole line), colnames: 1 2 3 4 5 6 7 8 9
    # [h1:H1,h2:H2] for rename tsv header

author: d2jvkpn
version: 1.4
release: 2019-06-10
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
		log.Printf("parsing GFF file: %s\n", os.Args[1])
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
		P2(ci.Scanner, strings.Split(os.Args[2], ","))
	case 4, 5:
		var renames []string
		if len(os.Args) == 5 {
			renames = strings.Split(os.Args[4], ",")
		}
		P3(ci.Scanner, strings.Split(os.Args[2], ","),
			strings.SplitN(os.Args[3], ",", -1), renames)
	default:
		fmt.Println(HELP)
	}
}

//
func gtfattr(s string, kv map[string]string, line int) bool {
	tmp := make(map[string]string)
	var msg string
	var ok bool

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
		ii[1] = strings.Trim(ii[1], "\"")

		if _, ok = tmp[ii[0]]; !ok {
			tmp[ii[0]] = ii[1] // make sure using the first k-v paire
		}
	}

	for k, v := range tmp {
		kv[k] = v
	}

	return true
}

func gffattr(s string, kv map[string]string, line int) bool {
	tmp := make(map[string]string)
	var msg string
	var ok bool
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

		if _, ok = tmp[ii[0]]; !ok {
			tmp[ii[0]] = ii[1] // make sure using the first k-v paire
		}
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

	array = append(array, []string{"sequences", strconv.Itoa(len(Sequences))})

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

		if types[0] != "" && StrSliceIndex(types, fds[2]) == -1 {
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
func P3(scanner *bufio.Scanner, types, attrs, renames []string) {
	var fds, x, rn []string
	var i, j int

	for i = 0; i < len(attrs)-1; i++ {
		if j = StrSliceIndex(attrs[i+1:], attrs[i]); j != -1 {
			log.Printf("duplicated attr field: %s\n", attrs[j])
			attrs = append(attrs[:j], attrs[j+1:]...)
		}
	}

	x = make([]string, len(attrs), len(attrs))
	copy(x, attrs)

	if j = StrSliceIndex(x, "0"); j > -1 {
		x[j] = "1\t2\t3\t4\t5\t6\t7\t8\t9"
	}

	if j = StrSliceIndex(x, "3"); j > -1 {
		x[j] = "type"
	}

	if j = StrSliceIndex(x, "4"); j > -1 {
		x[j] = "position"
	}

	if len(renames) > 0 {
		for j = range renames {
			rn = strings.SplitN(renames[j], ":", 2)

			if len(rn) != 2 {
				log.Fatalf("invalid rename field: %s\n", renames[j])
			}

			if StrSliceIndex(x, rn[0]) == -1 {
				log.Fatalf("head \"%s\" isn't in selected attrs: %s\n", rn[0])
			}

			x[StrSliceIndex(x, rn[0])] = rn[1]
		}
	}

	fmt.Println(strings.Join(x, "\t"))

	i = 0

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

		if types[0] != "" && StrSliceIndex(types, fds[2]) == -1 {
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

		kv["0"], kv["3"] = line, fds[2]
		kv["4"] = strings.Join([]string{fds[0], fds[3], fds[4], fds[6]}, ":")

		values := []string{}
		for _, k := range attrs {
			values = append(values, kv[k])
		}

		fmt.Println(strings.Join(values, "\t"))
	}
}

func SortStrSlice(s []string) {
	sort.Slice(s, func(i, j int) bool {
		return strings.ToLower(s[i]) < strings.ToLower(s[j])
	})
}

func StrSliceIndex(slice []string, value string) (p int) {
	for p = range slice {
		if slice[p] == value {
			return
		}
	}

	p = -1
	return
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
