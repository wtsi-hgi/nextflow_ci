process test_lustre_access {
    tag "check ${lustre_location}"
    publishDir "${params.outdir}/test_lustre_access/", mode: "${params.copy_mode}", pattern: "test.txt", overwrite: true

    when: 
    params.test_lustre_access.run_test

    input: 
    val(lustre_location)

    output: 
    tuple val(lustre_location), path('test.txt'), emit: test_txt

    script:
    """
umask 2 # make files group_writable

ls -ltra ${lustre_location} > test.txt
    """
}
