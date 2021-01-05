process test_lustre_access {
    tag "${lustre_location}"
    publishDir "${params.outdir}/test_lustre_access/", mode: 'copy', pattern: "test.txt", overwrite: true

    input: 
    val(lustre_location)

    output: 
    tuple val(lustre_location), path('test.txt'), emit: test_txt

    script:
    """
ls -ltra ${lustre_location} > test.txt
    """
}
