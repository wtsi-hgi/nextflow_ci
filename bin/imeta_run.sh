#!/usr/bin/env bash
id_run=$1
guides=$2

rm -f samples.tsv
rm -f guides_library.tsv

printf 'sample\tobject\tsample_supplier_name\tid_run\tis_paired_read\tstudy_id\tstudy\n' > samples.tsv
printf 'samplename,library,includeG\n' > guides_library.tsv

jq --arg id_run $id_run -n '{avus: [
       {attribute: "target", value: "1", o: "="},
       {attribute: "manual_qc", value: "1", o: "="}, 
       {attribute: "id_run", value: $id_run, o: "="}]}' |\
baton-metaquery \
		--zone seq --obj --avu |\
jq '.[] as $a| 
"\($a.avus | .[] | select(.attribute == "sample") | .value)____\($a.collection)/\($a.data_object)____\($a.avus | .[] | select(.attribute == "sample_supplier_name") | .value)____\($a.avus | .[] | select(.attribute == "id_run") | .value)____\($a.avus | .[] | select(.attribute == "is_paired_read") | .value)____\($a.avus | .[] | select(.attribute == "study_id") | .value)____\($a.avus | .[] | select(.attribute == "study") | .value)"' |\
    sed s"/$(printf '\t')//"g |\
    sed s"/\"//"g |\
    sed s"/____/$(printf '\t')/"g |\
    sort | uniq >> samples.tsv       

tail -n +2 samples.tsv | awk -v guides=$guides 'BEGIN{FS="\t";OFS=","}{print$1,guides,"False"}' | sort -u >> guides_library.tsv

echo jq search study id done
