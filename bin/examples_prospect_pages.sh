#!/bin/bash

#- Set of examples/tests to run `prospect_pages`

#- To be run with DESI data (eg at NERSC computing center)
# Requires to have DESI_SPECTRO_REDUX, DESICONDA and DESI_ROOT environment variables set.

#- Usage:
#  ./examples_prospect_pages     => run all examples/tests
#  ./examples_prospect_pages N   => run example/test N only

#- Base directory for output html pages:
OUTPUT_ROOT=${DESI_ROOT}/spectro/prospect/tests
#- With the default NERSC DESI environment, redrock templates used by prospect
# are as given from the module redrock-templates/main.
# For redrock fits produced before Aug 2022, one needs to use an earlier
# version of templates, matching those from module redrock-templates/0.7.2:
OLD_TEMPLATE_DIR=${DESICONDA}/../code/redrock-templates/0.7.2

# --------------------------
# Mode: Explicit input files
# --------------------------

#- 1) Simplest use: run prospect from a single file with spectra
#     Does not display any redrock fit result
#     424 spectra
if [[ $1 == 1 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 1 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/denali/tiles/cumulative/81067/20210327
    OUTPUTDIR=${OUTPUT_ROOT}/1
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   -o ${OUTPUTDIR}
fi

#- 2) Now make pages with redrock fit results in addition to spectra
if [[ $1 == 2 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 2 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/denali/tiles/cumulative/81067/20210327
    OUTPUTDIR=${OUTPUT_ROOT}/2
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/zbest-5-81067-thru20210327.fits \
                   --redrock_details_files ${DATAPATH}/redrock-5-81067-thru20210327.h5 \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   --outputdir ${OUTPUTDIR}
fi

#- 3) Options to modify html pages
#     Use --mask_type SV2_DESI_TARGET as this is SV2 data (needed to display targeting info)
#     Use --top_metadata to highlight some metadata (from fibermap)
#     Use --no-clean_fiberstatus => displays all 500 spectra
if [[ $1 == 3 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 3 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/denali/tiles/cumulative/81067/20210327
    OUTPUTDIR=${OUTPUT_ROOT}/3
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/zbest-5-81067-thru20210327.fits \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   -o ${OUTPUTDIR} \
                   --mask_type SV2_DESI_TARGET \
                   --top_metadata TARGETID TILEID mag_G mag_R \
                   --titlepage_prefix prospect_example3 \
                   --nspecperfile 20 \
                   --vi_countdown 5 \
                   --with_thumbnail_only_pages \
                   --no-clean_fiberstatus
fi

#- 4) Filter spectra: from metadata
#     Select BGS_ANY targets, redrock's deltachi2 and photometric r_mag
if [[ $1 == 4 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 4 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/denali/tiles/cumulative/81067/20210327
    OUTPUTDIR=${OUTPUT_ROOT}/4
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/zbest-5-81067-thru20210327.fits \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   -o ${OUTPUTDIR} \
                   --mask_type SV2_DESI_TARGET \
                   --targeting_mask BGS_ANY \
                   --rmag_min 19.8 \
                   --chi2_min 100
fi

#- 5) Process several spectra/zcat/redrock files (matched one-to-one)
#     Metadata filter: select some bright, high-SNR MWS objects
if [[ $1 == 5 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 5 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/denali/tiles/cumulative/81067/20210327
    OUTPUTDIR=${OUTPUT_ROOT}/5
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                                   ${DATAPATH}/coadd-7-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/zbest-5-81067-thru20210327.fits \
                                ${DATAPATH}/zbest-7-81067-thru20210327.fits \
                   --redrock_details_files ${DATAPATH}/redrock-5-81067-thru20210327.h5 \
                                           ${DATAPATH}/redrock-7-81067-thru20210327.h5 \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   -o ${OUTPUTDIR} \
                   --titlepage prospect_example_5 \
                   --mask_type SV2_DESI_TARGET \
                   --targeting_mask MWS_ANY \
                   --gmag_max 18 \
                   --snr_min 5
fi

#- 5a) Similar to (5), with guadalupe data release
if [[ $1 == 5a ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 5a ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/guadalupe/tiles/cumulative/1776/20210516
    OUTPUTDIR=${OUTPUT_ROOT}/5a
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    SP_LIST_FILE=tmp_sp_list.txt
    ZCAT_LIST_FILE=tmp_zcat_list.txt
    RR_LIST_FILE=tmp_rr_list.txt
    for i in {6..7}; do
        echo ${DATAPATH}/coadd-${i}-1776-thru20210516.fits >> ${SP_LIST_FILE}
        echo ${DATAPATH}/redrock-${i}-1776-thru20210516.fits >> ${ZCAT_LIST_FILE}
        echo ${DATAPATH}/rrdetails-${i}-1776-thru20210516.h5 >> ${RR_LIST_FILE}
    done
    prospect_pages --spectra_file_list ${SP_LIST_FILE} \
                   --zcat_file_list ${ZCAT_LIST_FILE} \
                   --redrock_details_file_list ${RR_LIST_FILE} \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   -o ${OUTPUTDIR} \
                   --titlepage prospect_example_5a \
                   --mask_type DESI_TARGET \
                   --targeting_mask QSO \
                   --snr_min 5
    rm -f ${SP_LIST_FILE} ${ZCAT_LIST_FILE} ${RR_LIST_FILE}
fi

#- 6) Filter spectra from a list of targets
#     Input target(s) by hand
if [[ $1 == 6 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 6 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/denali/tiles/cumulative/81067/20210327
    OUTPUTDIR=${OUTPUT_ROOT}/6
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/zbest-5-81067-thru20210327.fits \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   -o ${OUTPUTDIR} \
                   --targets 39633289673182563 39633289673182038 39633289673182074
fi

#- 6a) Same as (6), with a more recent data release (guadalupe)
if [[ $1 == 6a ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 6a ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/guadalupe/healpix/main/dark/309/30973
    OUTPUTDIR=${OUTPUT_ROOT}/6a
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-main-dark-30973.fits \
                   --zcat_files ${DATAPATH}/redrock-main-dark-30973.fits \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   -o ${OUTPUTDIR} \
                   --targets 39627770770234272 39627764738818698 39627764734628731
fi



# -------------------------
# Mode: Scan directory tree
# -------------------------

#- 7) Parse tiles/nights in andes release
#     Mandatory to specify dirtree_type, 'pernight' here
#     Filter over tiles, petals (here we want to process a small nb of files)
#     (in andes, tile 70004 has a single night, CMX targeting, petal 5 is missing)
#     No redrock info, no filter
if [[ $1 == 7 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 7 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/andes/tiles
    OUTPUTDIR=${OUTPUT_ROOT}/7
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type pernight \
                   -o ${OUTPUTDIR} \
                   --tiles 70004 \
                   --petals 5 6 7 \
                   --mask_type CMX_TARGET
fi

#- 8) Inspect deep coadds tiles/nights in blanc release
#     => Set dirtree_type = pernight and night = deep
#     Select QSO in COSMOS + LYNX tiles
if [[ $1 == 8 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 8 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/blanc/tiles
    OUTPUTDIR=${OUTPUT_ROOT}/8
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type pernight \
                   -o ${OUTPUTDIR} \
                   --with_zcatalog \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   --tiles 80607 80609 --nights deep \
                   --petals 6 7 \
                   --mask_type SV1_DESI_TARGET --targeting_mask QSO
fi

#- 9) Inspect "perexp" data in denali release
#     Input tile list from file
#     Stop creating pages after 300 spectra (last exposure is not processed)
if [[ $1 == 9 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 9 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/denali/tiles/perexp
    OUTPUTDIR=${OUTPUT_ROOT}/9
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    TILE_LIST_FILE=tmp_tile_list.txt
    echo -e ' 81058 \n 81059 \n 81060' > ${TILE_LIST_FILE}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type perexp \
                   -o ${OUTPUTDIR} \
                   --with_zcatalog \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   --tile_list ${TILE_LIST_FILE} \
                   --mask_type SV2_DESI_TARGET --targeting_mask LRG \
                   --chi2_min 200 --chi2_max 500 \
                   --nmax_spectra 300
    rm -f ${TILE_LIST_FILE}
fi

#- 10) Inspect "cumulative" data in everest release
#      Input target list from file (target order is kept in html pages)
#      Display multiple models (=> read readrock .h5 files), customize displayed metadata
if [[ $1 == 10 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 10 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/everest/tiles/cumulative
    OUTPUTDIR=${OUTPUT_ROOT}/10
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    TARGET_LIST_FILE=tmp_target_list.txt
    echo -e ' 616094114412233066 \n 39632930179383639 \n 39632930179384420 \n 39632930179384518 \n 616094111199396534 \n 616094114420622028 \n 616094111195201658 ' > ${TARGET_LIST_FILE}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type cumulative \
                   -o ${OUTPUTDIR} \
                   --with_zcatalog --with_multiple_models \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   --tiles 81062 80654 \
                   --target_list ${TARGET_LIST_FILE} \
                   --mask_type DESI_TARGET \
                   --top_metadata TARGETID TILEID COADD_EXPTIME
    rm -f ${TARGET_LIST_FILE}
fi

#- 11) Inspect cframes in "exposure" directory trees
#      No zcatalog in that case
if [[ $1 == 11 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 11 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/daily/exposures
    OUTPUTDIR=${OUTPUT_ROOT}/11
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type exposures --spectra_type cframe \
                   -o ${OUTPUTDIR} \
                   --nights 20210628 --expids 00096449 00096450 --petals 2 \
                   --titlepage prospect-frames
fi

#- 12) Inspect data in healpix directories, everest
if [[ $1 == 12 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 12 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/everest/healpix
    OUTPUTDIR=${OUTPUT_ROOT}/12
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type healpix \
                   --with_zcatalog --with_multiple_models \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   -o ${OUTPUTDIR} \
                   --pixels 9557 --survey_program main dark \
                   --targeting_mask QSO --chi2_min 200
fi
