#!/bin/bash

#- Set of examples/tests to run `prospect_pages`
#- To be run with DESI data (eg at NERSC computing center)

#- Usage:
#  ./examples_prospect_pages     => run all examples/tests
#  ./examples_prospect_pages N   => run example/test N only

# All of the 15 test cases can serve as examples.
# Suggestions to consider first:
#  2) single input coadd + redrock files
#  14) parse a healpix directory

#- Base directory for output html pages:
OUTPUT_ROOT=${DESI_ROOT}/spectro/prospect/tests

#- With the default NERSC DESI environment, redrock templates used by prospect
# are as given from the module redrock-templates/main.
# For redrock fits produced before Aug 2022, one needs to use an earlier
# version of templates, matching those from module redrock-templates/0.7.2:
OLD_TEMPLATE_DIR=${DESICONDA}/../code/redrock-templates/0.7.2
if [ ! -d $OLD_TEMPLATE_DIR ] ; then
    OLD_TEMPLATE_DIR=${DESICONDA}/../../20220119-2.0.1/code/redrock-templates/0.7.2
fi

#- Check env vars / paths:
# The script requires to have TMPDIR, DESI_SPECTRO_REDUX, DESICONDA and DESI_ROOT environment variables set.
CheckPaths=(TMPDIR DESI_SPECTRO_REDUX DESICONDA DESI_ROOT OUTPUT_ROOT OLD_TEMPLATE_DIR)
for var in "${CheckPaths[@]}"; do
    if [ ! -d "${!var}" ] ; then echo "Path does not exist: "${!var}" ("$var")" && exit ; fi
done

# --------------------------
# Mode: Explicit input files
# --------------------------

#- 1) Simplest use: run prospect from a single file with spectra
#     Does not display any redrock fit result
if [[ $1 == 1 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 1 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/iron/tiles/cumulative/81067/20210327
    OUTPUTDIR=${OUTPUT_ROOT}/1
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   -o ${OUTPUTDIR}
fi

#- 2) Now make pages with redrock fit results in addition to spectra
if [[ $1 == 2 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 2 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/iron/tiles/cumulative/81067/20210327
    OUTPUTDIR=${OUTPUT_ROOT}/2
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/redrock-5-81067-thru20210327.fits \
                   --redrock_details_files ${DATAPATH}/rrdetails-5-81067-thru20210327.h5 \
                   --outputdir ${OUTPUTDIR}
fi

#- 3) Options to modify html pages, eg:
#     Use --mask_type SV2_DESI_TARGET as this is SV2 data (needed to display targeting info)
#     Use --top_metadata to highlight some metadata (from fibermap)
#     Use --no_clean_fiberstatus => displays all 500 spectra
if [[ $1 == 3 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 3 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/iron/tiles/cumulative/81067/20210327
    OUTPUTDIR=${OUTPUT_ROOT}/3
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/redrock-5-81067-thru20210327.fits \
                   -o ${OUTPUTDIR} \
                   --mask_type SV2_DESI_TARGET \
                   --top_metadata TARGETID TILEID mag_G mag_R \
                   --titlepage_prefix prospect_example3 \
                   --nspec_per_page 20 \
                   --vi_countdown 5 \
                   --with_thumbnail_only_pages \
                   --no_clean_fiberstatus
fi

#- 4) Filter spectra: from metadata
#     Select BGS_ANY targets, redrock's deltachi2 and photometric r_mag
#     Using fuji (~ EDR)
if [[ $1 == 4 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 4 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/fuji/tiles/cumulative/81067/20210327
    OUTPUTDIR=${OUTPUT_ROOT}/4
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/redrock-5-81067-thru20210327.fits \
                   -o ${OUTPUTDIR} \
                   --mask_type SV2_DESI_TARGET \
                   --targeting_mask BGS_ANY \
                   --rmag_min 19.8 \
                   --chi2_min 100 \
                   --titlepage_prefix example_fuji_bgs_badchi2
fi

#- 5) Process several spectra/zcat/redrock files (matched one-to-one)
#     Metadata filter: select some bright, high-SNR MWS objects
if [[ $1 == 5 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 5 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/iron/healpix/main/bright/280
    OUTPUTDIR=${OUTPUT_ROOT}/5
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/28012/coadd-main-bright-28012.fits \
                                   ${DATAPATH}/28013/coadd-main-bright-28013.fits \
                   --zcat_files ${DATAPATH}/28012/redrock-main-bright-28012.fits \
                                ${DATAPATH}/28013/redrock-main-bright-28013.fits \
                   --redrock_details_files ${DATAPATH}/28012/rrdetails-main-bright-28012.h5 \
                                           ${DATAPATH}/28013/rrdetails-main-bright-28013.h5 \
                   -o ${OUTPUTDIR} \
                   --titlepage prospect_example_5 \
                   --targeting_mask MWS_ANY \
                   --gmag_max 19 \
                   --snr_min 4
fi

#- 6) Similar to (5), with guadalupe data release
#      The lists of spectra/redrock files is provided in the form of ASCII files
if [[ $1 == 6 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 6 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/guadalupe/tiles/cumulative/1776/20210516
    OUTPUTDIR=${OUTPUT_ROOT}/6
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    SP_LIST_FILE=${TMPDIR}/tmp_sp_list.txt
    ZCAT_LIST_FILE=${TMPDIR}/tmp_zcat_list.txt
    RR_LIST_FILE=${TMPDIR}/tmp_rr_list.txt
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
                   --titlepage prospect_example_6 \
                   --mask_type DESI_TARGET \
                   --targeting_mask QSO \
                   --snr_min 5
    rm -f ${SP_LIST_FILE} ${ZCAT_LIST_FILE} ${RR_LIST_FILE}
fi

#- 7) Filter spectra from a list of targets
#     Input targets with an ASCII file
if [[ $1 == 7 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 7 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/iron/tiles/cumulative/81067/20210327
    OUTPUTDIR=${OUTPUT_ROOT}/7
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    TARGET_LIST_FILE=${TMPDIR}/tmp_target_list.txt
    echo -e "39633289673182563\n39633289673182038\n39633289673182074" > ${TARGET_LIST_FILE}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/redrock-5-81067-thru20210327.fits \
                   -o ${OUTPUTDIR} \
                   --targets 39633289673182563 39633289673182038 39633289673182074
    rm -f ${TARGET_LIST_FILE}
fi

#- 8) Similar to (7), checking himalayas data release
#     Input targets by hand
if [[ $1 == 8 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 8 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/himalayas/healpix/main/dark/309/30973
    OUTPUTDIR=${OUTPUT_ROOT}/8
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-main-dark-30973.fits \
                   --zcat_files ${DATAPATH}/redrock-main-dark-30973.fits \
                   --redrock_details_files ${DATAPATH}/rrdetails-main-dark-30973.h5 \
                   -o ${OUTPUTDIR} \
                   --targets 39627770770234272 39627764738818698 39627764734628731
fi

#- 9) A sample with high-PM stars (A. Dey),
#     A second cross-hair is then visible in imaging thumbs
if [[ $1 == 9 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 9 ------"
    OUTPUTDIR=${OUTPUT_ROOT}/9
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    DATAPATH=${DESI_SPECTRO_REDUX}/iron/healpix/main/bright
    TARGET_LIST_FILE=${TMPDIR}/tmp_target_list.txt
    SPECTRA_LIST_FILE=${TMPDIR}/tmp_spectra_list.txt
    echo -e "2376933706825728 \n 2376982100705280 \n 2377019157381120 \n 2376644622811136 \n 2376978715901952" > ${TARGET_LIST_FILE}
    echo -e "${DATAPATH}/7/733/coadd-main-bright-733.fits \n ${DATAPATH}/21/2110/coadd-main-bright-2110.fits \n ${DATAPATH}/22/2202/coadd-main-bright-2202.fits \n ${DATAPATH}/42/4285/coadd-main-bright-4285.fits \n ${DATAPATH}/52/5264/coadd-main-bright-5264.fits" > ${SPECTRA_LIST_FILE}
    prospect_pages --spectra_file_list ${SPECTRA_LIST_FILE} \
                   --titlepage_prefix high-pm-stars \
                   --target_list ${TARGET_LIST_FILE} \
                   -o ${OUTPUTDIR}
    rm -f ${TARGET_LIST_FILE} ${SPECTRA_LIST_FILE}
fi

# -------------------------
# Mode: Scan directory tree
# -------------------------

#- 10) Parse tiles/nights in andes release
#     Mandatory to specify dirtree_type, 'pernight' here
#     Filter over tiles, petals
#     (in andes, tile 70004 has a single night, CMX targeting, petal 5 is missing)
#     No redrock info, no filter
if [[ $1 == 10 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 10 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/andes/tiles
    OUTPUTDIR=${OUTPUT_ROOT}/10
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type pernight \
                   -o ${OUTPUTDIR} \
                   --tiles 70004 \
                   --petals 5 6 7 \
                   --mask_type CMX_TARGET
fi

#- 11) Inspect deep coadds tiles/nights in blanc release
#     => Set dirtree_type = pernight and night = deep
#     Select QSO in COSMOS + LYNX tiles
if [[ $1 == 11 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 11 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/blanc/tiles
    OUTPUTDIR=${OUTPUT_ROOT}/11
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

#- 12) Inspect "cumulative" data in iron release
#      Input target list from file (target order is kept in html pages)
#      Customize displayed metadata
if [[ $1 == 12 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 12 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/iron/tiles/cumulative
    OUTPUTDIR=${OUTPUT_ROOT}/12
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    TARGET_LIST_FILE=${TMPDIR}/tmp_target_list.txt
    echo -e ' 616094114412233066 \n 39632930179383639 \n 39632930179384420 \n 39632930179384518 \n 616094111199396534 \n 616094114420622028 \n 616094111195201658 ' > ${TARGET_LIST_FILE}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type cumulative \
                   -o ${OUTPUTDIR} \
                   --with_zcatalog --with_multiple_models \
                   --tiles 81062 80654 \
                   --target_list ${TARGET_LIST_FILE} \
                   --mask_type DESI_TARGET \
                   --top_metadata TARGETID TILEID COADD_EXPTIME
    rm -f ${TARGET_LIST_FILE}
fi

#- 13) Inspect cframes in exposure-based directory trees (no zcatalog in that case)
#      From daily production
if [[ $1 == 13 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 13 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/daily/exposures
    OUTPUTDIR=${OUTPUT_ROOT}/13
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type exposures --spectra_type cframe \
                   -o ${OUTPUTDIR} \
                   --nights 20230405 --expids 00175085 00175084 --petals 2 \
                   --titlepage prospect-frames
fi

#- 14) Inspect data in healpix directories, iron
if [[ $1 == 14 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 14 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/iron/healpix
    OUTPUTDIR=${OUTPUT_ROOT}/14
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type healpix \
                   --with_zcatalog --with_multiple_models \
                   -o ${OUTPUTDIR} \
                   --pixels 9557 --survey_program main dark \
                   --targeting_mask QSO --chi2_min 200
fi

#- 15) Inspect data in cascades release, a subsample of the SV0 LRG VI set
#     Input tile list from file
#     Stop creating pages after 300 spectra
if [[ $1 == 15 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 15 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/cascades/tiles
    OUTPUTDIR=${OUTPUT_ROOT}/15
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    TILE_LIST_FILE=${TMPDIR}/tmp_tile_list.txt
    echo -e '80609 \n 80607' > ${TILE_LIST_FILE}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type pernight \
                   -o ${OUTPUTDIR} \
                   --with_zcatalog \
                   --template_dir ${OLD_TEMPLATE_DIR} \
                   --tile_list ${TILE_LIST_FILE} \
                   --mask_type SV1_DESI_TARGET --targeting_mask LRG \
                   --chi2_min 200 --chi2_max 300 \
                   --petals 2 4 6 \
                   --nmax_spectra 300
    rm -f ${TILE_LIST_FILE}
fi

#- 16) Inspect bright stars in tiles/cumulative directories, iron
#       Customize html page displays
if [[ $1 == 16 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 16 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/iron/tiles/cumulative
    OUTPUTDIR=${OUTPUT_ROOT}/16
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type cumulative \
                   --with_zcatalog \
                   -o ${OUTPUTDIR} \
                   --tiles 23937 24780 \
                   --targeting_mask STD_BRIGHT \
                   --no_imaging --no_noise --no_vi_widgets
fi

# -------------------------
# Some other tests
#  (those are more tests than useful examples)
# -------------------------

#- 17) Inspect a single spectrum from daily spectra
#      Which turns out to have a missing 'r' band (spectra.bands = ['b', 'z'])
if [[ $1 == 17 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 17 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/daily/tiles/cumulative/1823/20210616
    OUTPUTDIR=${OUTPUT_ROOT}/17
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-9-1823-thru20210616.fits \
                   --zcat_files ${DATAPATH}/zbest-9-1823-thru20210616.fits \
                   --redrock_details_files ${DATAPATH}/redrock-9-1823-thru20210616.h5 \
                   --targets 39627747886112324 \
                   --outputdir ${OUTPUTDIR}
fi

#- 18) A single spectrum from iron spectra
#      With a completely masked, though existing, 'r' band
if [[ $1 == 18 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 18 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/iron/healpix/main/dark/134/13455
    OUTPUTDIR=${OUTPUT_ROOT}/18
    [ ! -d ${OUTPUTDIR} ] && mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-main-dark-13455.fits \
                   --zcat_files ${DATAPATH}/redrock-main-dark-13455.fits \
                   --redrock_details_files ${DATAPATH}/rrdetails-main-dark-13455.h5 \
                   --targets 39628472435347978 \
                   --no_clean_fiberstatus \
                   --outputdir ${OUTPUTDIR}
fi
