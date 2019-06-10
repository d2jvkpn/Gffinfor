<pre>
GFF/GTF(.gz) summary and attributions extraction, usage:
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
lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html
</pre>
