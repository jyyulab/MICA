psub -K -P test_no_install_LSF -J test_no_install_LSF_MICA_Transform.0 -q compbio -M 2000 -i ./test_data/outputs/test_no_install_LSF/.tmp/02_Transform_test_no_install_LSF_mds.sh -oo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Transform.0.%J.%I.out -eo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Transform.0.%J.%I.err 
psub -K -P test_no_install_LSF -J test_no_install_LSF_MICA_Transform.1 -q compbio -M 2000 -i ./test_data/outputs/test_no_install_LSF/.tmp/03_Share_test_no_install_LSF_mds.sh -oo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Transform.1.%J.%I.out -eo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Transform.0.%J.%I.err 
psub -K -P test_no_install_LSF -J test_no_install_LSF_MICA_Kmeans.0 -q compbio -M 2000 -i ./test_data/outputs/test_no_install_LSF/.tmp/04_Prep_test_no_install_LSF_mds.sh -oo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Kmeans.0.%J.%I.out -eo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Kmeans.0.%J.%I.err 
psub -K -P test_no_install_LSF -J test_no_install_LSF_MICA_Kmeans.1 -q compbio -M 2000 -i ./test_data/outputs/test_no_install_LSF/.tmp/05_Kmeans_test_no_install_LSF_mds.sh -oo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Kmeans.1.%J.%I.out -eo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Kmeans.1.%J.%I.err 
psub -K -P test_no_install_LSF -J test_no_install_LSF_MICA_CClust -q compbio -M 2000 -i ./test_data/outputs/test_no_install_LSF/.tmp/06_CClust_test_no_install_LSF_mds.sh -oo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Cclust.%J.%I.out -eo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Cclust.%J.%I.err 
psub -K -P test_no_install_LSF -J test_no_install_LSF_MICA_GGplot -q compbio -M 2000 -i ./test_data/outputs/test_no_install_LSF/.tmp/07_GGplot_test_no_install_LSF_mds.sh -oo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Ggplot.%J.%I.out -eo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Ggplot.%J.%I.err 
psub -K -P test_no_install_LSF -J test_no_install_LSF_MICA_Clean -q compbio -M 2000 -i ./test_data/outputs/test_no_install_LSF/.tmp/08_Clean_test_no_install_LSF_mds.sh -oo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Clean.%J.%I.out -eo ./test_data/outputs/test_no_install_LSF/.tmp/.log/test_no_install_LSF_MICA_Clean.%J.%I.err 
