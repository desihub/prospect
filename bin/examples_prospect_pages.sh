#!/bin/bash

#- Set of examples/tests to run `prospect_pages`

#- To be run with DESI data (eg at NERSC computing center)
# Requires to have DESI_SPECTRO_REDUX and DESI_ROOT environment variables set.

#- Usage:
#  ./examples_prospect_pages     => run all examples/tests
#  ./examples_prospect_pages N   => run example/test N only


# --------------------------
# Mode: Explicit input files
# --------------------------

#- 1) Simplest use: run prospect from a single file with spectra
#     Does not display any redrock fit result
#     424 spectra
if [[ $1 == 1 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 1 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/denali/tiles/cumulative/81067/20210327
    OUTPUTDIR=${DESI_ROOT}/spectro/prospect/tests/1
    mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   -o ${OUTPUTDIR}
fi

#- 2) Now make pages with redrock fit results in addition to spectra
if [[ $1 == 2 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 2 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/denali/tiles/cumulative/81067/20210327
    OUTPUTDIR=${DESI_ROOT}/spectro/prospect/tests/2
    mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/zbest-5-81067-thru20210327.fits \
                   --redrock_files ${DATAPATH}/redrock-5-81067-thru20210327.h5 \
                   --outputdir ${OUTPUTDIR}
fi

#- 3) Options to modify html pages
#     Use --mask_type SV2_DESI_TARGET as this is SV2 data (needed to display targeting info)
#     Use --top_metadata to highlight some metadata (from fibermap)
#     Use --no-clean_fiberstatus => displays all 500 spectra
if [[ $1 == 3 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 3 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/denali/tiles/cumulative/81067/20210327
    OUTPUTDIR=${DESI_ROOT}/spectro/prospect/tests/3
    mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/zbest-5-81067-thru20210327.fits \
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
    OUTPUTDIR=${DESI_ROOT}/spectro/prospect/tests/4
    mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/zbest-5-81067-thru20210327.fits \
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
    OUTPUTDIR=${DESI_ROOT}/spectro/prospect/tests/5
    mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                                   ${DATAPATH}/coadd-7-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/zbest-5-81067-thru20210327.fits \
                                ${DATAPATH}/zbest-7-81067-thru20210327.fits \
                   --redrock_files ${DATAPATH}/redrock-5-81067-thru20210327.h5 \
                                   ${DATAPATH}/redrock-7-81067-thru20210327.h5 \
                   -o ${OUTPUTDIR} \
                   --titlepage prospect_example_5 \
                   --mask_type SV2_DESI_TARGET \
                   --targeting_mask MWS_ANY \
                   --gmag_max 18 \
                   --snr_min 5
fi

#- 6) Filter spectra from a list of targets
#     Input target(s) by hand
if [[ $1 == 6 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 6 ------"
    DATAPATH=${DESI_SPECTRO_REDUX}/denali/tiles/cumulative/81067/20210327
    OUTPUTDIR=${DESI_ROOT}/spectro/prospect/tests/6
    mkdir ${OUTPUTDIR}
    prospect_pages --spectra_files ${DATAPATH}/coadd-5-81067-thru20210327.fits \
                   --zcat_files ${DATAPATH}/zbest-5-81067-thru20210327.fits \
                   -o ${OUTPUTDIR} \
                   --targets 39633289673182563 39633289673182038 39633289673182074
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
    OUTPUTDIR=${DESI_ROOT}/spectro/prospect/tests/7
    mkdir ${OUTPUTDIR}
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
    OUTPUTDIR=${DESI_ROOT}/spectro/prospect/tests/8
    mkdir ${OUTPUTDIR}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type pernight \
                   -o ${OUTPUTDIR} \
                   --with_zcatalog \
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
    OUTPUTDIR=${DESI_ROOT}/spectro/prospect/tests/9
    mkdir ${OUTPUTDIR}
    TILE_LIST_FILE=tmp_tile_list.txt
    echo -e ' 81058 \n 81059 \n 81060' > ${TILE_LIST_FILE}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type perexp \
                   -o ${OUTPUTDIR} \
                   --with_zcatalog \
                   --tile_list ${TILE_LIST_FILE} \
                   --mask_type SV2_DESI_TARGET --targeting_mask LRG \
                   --chi2_min 200 --chi2_max 500 \
                   --nmax_spectra 300
    rm -f ${TILE_LIST_FILE}
fi

#- 10) Inspect "cumulative" data in denali release
#      Input target list from file (target order is kept in html pages)
#      Display multiple models (=> read readrock .h5 files), customize displayed metadata
if [[ $1 == 10 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 10 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/denali/tiles/cumulative
    OUTPUTDIR=${DESI_ROOT}/spectro/prospect/tests/10
    mkdir ${OUTPUTDIR}
    TARGET_LIST_FILE=tmp_target_list.txt
    echo -e ' 616094114412233066 \n 39632930179383639 \n 39632930179384420 \n 39632930179384518 \n 616094111199396534 \n 616094114420622028 \n 616094111195201658 ' > ${TARGET_LIST_FILE}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type cumulative \
                   -o ${OUTPUTDIR} \
                   --with_zcatalog --with_multiple_models \
                   --tiles 81062 80654 \
                   --target_list ${TARGET_LIST_FILE} \
                   --mask_type DESI_TARGET \
                   --top_metadata TARGETID TILEID EXPID
    rm -f ${TARGET_LIST_FILE}
fi

#- 11) Inspect cframes in "exposure" directory trees
#      No zcatalog in that case
if [[ $1 == 11 ]] || [[ $1 == '' ]]; then
    echo "------ Example/Test 11 ------"
    DATADIR=${DESI_SPECTRO_REDUX}/daily/exposures
    OUTPUTDIR=${DESI_ROOT}/spectro/prospect/tests/11
    mkdir ${OUTPUTDIR}
    prospect_pages --datadir ${DATADIR} \
                   --dirtree_type exposures --spectra_type cframe \
                   -o ${OUTPUTDIR} \
                   --nights 20210628 --expids 00096449 00096450 --petals 2 \
                   --titlepage prospect-frames
fi
