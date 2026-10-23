// Test harness: print the task pod spec conf/site_nrp.config builds (tests/test_k8s_avoid_nodes.sh).
workflow { println "POD=${params.nrp_task_pod}" }
