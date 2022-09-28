#!/bin/bash

# Source the env.sh file for specific paths, etc
. "$(dirname $0)/env.sh"

# Create a temp space if it does not exist
if [[ ! -d $TMPDIR ]]; then
  TMPDIR=$(mktemp -d /tmp/histreg_XXXXXX)
fi

# This function defines variables for the whole experiment
function set_experiment_vars()
{
  # User must pass parameter json file as input
  JSON=${1?}

  # Read important fields from the JSON file
  ID=$(jq -r '.subject_id' < $JSON)
  MRI_ID=$(jq -r '.mri_id' < $JSON)
  SLAB=$(jq -r '.slab' < $JSON)
  RSTAIN=$(jq -r '.reference_stain' < $JSON)

  # Read the list of stains provided for registration
  STAIN_LIST=$(echo $(jq -r '.histology_set | keys[]' < $JSON))

  # Work directory for this experiment
  WDIR=$ROOT/work/$ID/slab$(printf %02d $SLAB)

  # Manual input directory for this experiment
  MANUAL_DIR=$ROOT/manual/$ID

  # MRI of the hemisphere, resliced into correct orientation
  MRI=$ROOT/input/$ID/${MRI_ID}_reslice.nii.gz

  # Slab mask with dot segmentations
  SLAB_ID="slab$(printf %02d $SLAB)"
  SLAB_MASK=$ROOT/input/$ID/${MRI_ID}_${SLAB_ID}_mask_with_dots.nii.gz

  # MRI of the slab
  MRI_SLAB=$WDIR/${MRI_ID}_slab$(printf %02d $SLAB)_slab.nii.gz

  # Input workspace for the manual registration between MRI and reference slide
  MRI_TO_REFSLIDE_INPUT_WORKSPACE=$MANUAL_DIR/mri_to_refslide_manual_input.itksnap

  # Result ITK-SNAP workspace (after completing registration)
  MRI_TO_REFSLIDE_MANUAL_WORKSPACE=$MANUAL_DIR/mri_to_refslide_manual_result.itksnap

  # Initial manual registration between reference slide (as fixed)
  # and MRI slab (as moving)
  MRI_TO_REFSLIDE_INIT=$WDIR/mri_to_refslide_init.mat
  MRI_TO_REFSLIDE_MANUAL=$WDIR/mri_to_refslide_manual_rigid.mat
  MRI_TO_REFSLIDE_RIGID=$WDIR/mri_to_refslide_rigid.mat
  MRI_TO_REFSLIDE_AFFINE=$WDIR/mri_to_refslide_affine.mat
  MRI_TO_REFSLIDE_WARP=$WDIR/mri_to_refslide_warp.mhd
  MRI_TO_REFSLIDE_RESLICE_AFFINE=$WDIR/mri_to_refslide_reslice_affine.nii.gz
  MRI_TO_REFSLIDE_RESLICE_DEFORM=$WDIR/mri_to_refslide_reslice_deform.nii.gz
  MRI_TO_REFSLIDE_RESULT_WORKSPACE=$WDIR/mri_to_refslide_greedy_result.itksnap

}

# This function defines variables for a particular stain
function set_stain_vars()
{
  # User must pass a stain and a parameter json file as input
  local JSON=${1?}
  local STAIN=${2?}

  # Read the id of the slide with this stain
  SLIDE_ID=$(jq -r ".histology_set.${STAIN}" < $JSON)

  # Filename of the raw slide
  SLIDE_RAW=$ROOT/input/${ID}/${SLIDE_ID}.nii.gz

  # Work directory for this stain
  SLIDE_DIR=$WDIR/slides/${SLIDE_ID}

  # Scalar version of the slide
  SLIDE_SCALAR=$SLIDE_DIR/${SLIDE_ID}_gray.nii.gz

  # Global rigid transform between reference and this stain
  SLIDE_TO_REF_GLOBAL_RIGID=$SLIDE_DIR/${STAIN}_to_${RSTAIN}_global_rigid.mat

  # Piecewise rigid transform between reference and this stain. This is a
  # pattern, each %02d is replaced by the chunk number
  SLIDE_TO_REF_CHUNK_RIGID=$SLIDE_DIR/${STAIN}_to_${RSTAIN}_chunk_rigid_%02d.mat

  # Piecewise deformable transformation
  SLIDE_TO_REF_CHUNK_WARP=$SLIDE_DIR/${STAIN}_to_${RSTAIN}_chunk_warp_%02d.nii.gz
  SLIDE_TO_REF_CHUNK_INVWARP=$SLIDE_DIR/${STAIN}_to_${RSTAIN}_chunk_invwarp_%02d.nii.gz
  SLIDE_TO_REF_CHUNK_METRIC=$SLIDE_DIR/${STAIN}_to_${RSTAIN}_metric_deform_%02d.nii.gz

  # Slide resampled to reference slide
  SLIDE_TO_REF_RESLICE_DEFORM_RAW=$SLIDE_DIR/${STAIN}_to_${RSTAIN}_reslice_deform_raw.nii.gz

  # For the reference stain we define additional variables
  if [[ $STAIN == $RSTAIN ]]; then
    REFSLIDE_RAW=$SLIDE_RAW
    REFSLIDE_SCALAR=$SLIDE_SCALAR
    REFSLIDE_MASK=$ROOT/input/${ID}/${SLIDE_ID}_mask.nii.gz
    REFSLIDE_CHUNK_MASK=$SLIDE_DIR/${SLIDE_ID}_chunk_mask.nii.gz

    # Lower-resolution version of the reference slide for matching to MRI
    REFSLIDE_MRILIKE=$SLIDE_DIR/${SLIDE_ID}_gray_lores.nii.gz
    REFSLIDE_MRILIKE_MASK=$SLIDE_DIR/${SLIDE_ID}_mask_lores.nii.gz

  fi
}

# Helper function - check single input file for existence
function check_file()
{
  local RED='\033[0;31m'
  local GREEN='\033[0;32m'
  local NC='\033[0m' # No Color
  if [[ -f ${1?} ]]; then
    echo -e "${GREEN}$(echo $1 | sed -e "s|$PWD/|./|")${NC}"
  else
    echo -e "${RED}Not Found${NC}"
  fi
}

# Check if all the inputs are located
function check_inputs()
{
  # User must pass a stain and a parameter json file as input
  JSON=${1?}

  # Set the variables
  set_experiment_vars $JSON
  set_stain_vars $JSON $RSTAIN

  # Check the MRI inputs
  echo "=== Checking inputs for $ID-$SLAB_ID ==="
  printf "%-30s%s\n" "MRI:" $(check_file $MRI)
  printf "%-30s%s\n" "MRI slab mask:" "$(check_file $SLAB_MASK)"

  # Check the reference stain
  printf "%-30s%s\n" "Stain $RSTAIN slide:" $(check_file $REFSLIDE_RAW)
  printf "%-30s%s\n" "Stain $RSTAIN mask:" $(check_file $REFSLIDE_MASK)

  # Check each of the other stains
  for stain in $STAIN_LIST; do
    if [[ $stain != $RSTAIN ]]; then
      printf "%-30s%s\n" "Stain $stain slide:" $(check_file $SLIDE_RAW)
    fi
  done
}

# Preprocess the MRI
function mri_preprocess()
{
  # User must pass a stain and a parameter json file as input
  JSON=${1?}

  # Set the variables
  set_experiment_vars $JSON
  mkdir -p $WDIR

  # Generate the slab MRI (useful for manual initialization of registration)
  c3d $SLAB_MASK $MRI -reslice-identity -o $MRI_SLAB
}

# This function applies preprocessing to a stain
function stain_preprocess()
{
  # User must pass a stain and a parameter json file as input
  local JSON=${1?}
  local STAIN=${2?}

  # Set the variables
  set_experiment_vars $JSON
  set_stain_vars $JSON $STAIN
  mkdir -p $SLIDE_DIR

  # Extract the grayscale component from the stain. For now, we just pick one
  # color channel that makes it look most like the MRI. This is really basic
  # and in the future we should replace this with a translation module that maps
  # histology intensity to MRI-like intensity
  local COMPONENT=0
  if [[ $STAIN=="HE" || $STAIN=="LBF" ]]; then
    # For HE and LBF, green works best
    COMPONENT=1
  elif $STAIN=="IBA"; then
    # For IBA, red works best
    COMPONENT=0
  fi

  # Extract the scalar component
  c2d -mcs $SLIDE_RAW -pick $COMPONENT -o $SLIDE_SCALAR

  # For the reference stain, perform additional processing
  if [[ $STAIN == $RSTAIN ]]; then

    # Determine the number of "chunks" into which to break the binary mask. These
    # chunks are used for piecewise registration between stains. It may be good to
    # let the user specify the number of chunks optionally in the JSON file
    local VOL=$(c2d $REFSLIDE_MASK  -dup -lstat | awk '$1==1 {print $7}')
    local N_CHUNKS=$(echo "import math; print(int(math.ceil($VOL/120.)))" | python)

    # Cut the mask into pieces
    if [[ $N_CHUNKS -gt 1 ]]; then
      image_graph_cut -u 1.2 -n 100 -c 4 0.1 $REFSLIDE_MASK $REFSLIDE_CHUNK_MASK $N_CHUNKS
    else
      cp -av $REFSLIDE_MASK $REFSLIDE_CHUNK_MASK
    fi 

    # Generate MRI-like appearance. For this we multiply the scalar image by the mask, but
    # first we refine the mask, making it a scalar image from 0 to 1 (probability of a voxel
    # being a foreground voxel). To do this, we train a random forest classifier with the 
    # input mask as training data.
    c2d -mcs $REFSLIDE_RAW $REFSLIDE_MASK \
      -shift 1 -rf-param-patch 2x2 -rf-train $TMPDIR/mask.rf \
      $REFSLIDE_RAW -rf-apply $TMPDIR/mask.rf -popas PFG \
      $REFSLIDE_SCALAR -push PFG -times \
      -smooth-fast 0.15mm -resample-mm 0.3mm -o $REFSLIDE_MRILIKE \
      -push PFG -smooth-fast 0.15mm -reslice-identity \
      -thresh 0.5 inf 1 0 -type uchar -o $REFSLIDE_MRILIKE_MASK

    # Generate an initial transformation between the slide and the MRI that simply
    # accounts for image centers and differences in orientation - the rest has to be
    # done by the user
    c3d_affine_tool \
      -tran $(c3d $MRI_SLAB -probe 50% | awk '{print $5,$6,$7}') \
      -rot 90 1 0 0 \
      -tran $(c3d $REFSLIDE_MRILIKE -probe 50% | awk '{print $5,$6,$7}') \
      -inv -mult -mult -info -o $MRI_TO_REFSLIDE_INIT

    # Create the initial workspace for the registration between the reference stain and
    # the MRI - user opens workspace, performs initial alignment, and saves as the
    # result workspace 
    mkdir -p $MANUAL_DIR
    itksnap-wt \
      -laa $REFSLIDE_MRILIKE -psn "Reference slide (MRI-like)" \
      -laa $MRI_SLAB -psn "Slab MRI" -props-set-transform $MRI_TO_REFSLIDE_INIT \
      -las $REFSLIDE_MRILIKE_MASK -psn "Mask" \
      -o $MRI_TO_REFSLIDE_INPUT_WORKSPACE

  fi
}

# Apply reslicing for chunked transformations, where the outside of the mask
# has to be set to a background value
function chunked_warp_reslice()
{
  local fixed_mask moving mov_coord_ref warp result background dilation greedy_opts
  read -r fixed_mask moving mov_coord_ref warp result background dilation greedy_opts <<< "$@"

  # Fix the header of the moving image to match a reference image
  local mov_header_fix=$TMPDIR/moving_header_fix.nii.gz
  c2d $mov_coord_ref -popas R -mcs $moving \
    -foreach -insert R 1 -mbb -endfor -omc $mov_header_fix

  # Apply the transformation
  local reslice_tmp=$TMPDIR/tmp_reslice.nii.gz
  greedy -d 2 $greedy_opts -rf $fixed_mask -rm $mov_header_fix $reslice_tmp -r $warp

  # Background defaults to zero
  if [[ ! $background ]]; then
    background="0"
  fi

  # Erosion defaults to zero too
  if [[ ! $dilation ]]; then
    local dilation_cmd=""
  else
    local dilation_cmd="-dilate 1 ${dilation}x${dilation}vox"
  fi

  # Clean up outside of the mask
  c2d \
    $fixed_mask -thresh 1 inf 1 0 $dilation_cmd -as M \
    -thresh 0 0 $background 0 -popas B \
    -mcs $reslice_tmp -foreach -push M -times -push B -add -endfor \
    -omc $result
}

# This function performs registration between reference stain (as fixed image)
# and another stain (as moving image)
function stain_register_to_reference()
{
  # User must pass a stain and a parameter json file as input
  JSON=${1?}
  STAIN=${2?}

  # Set the variables
  set_experiment_vars $JSON
  set_stain_vars $JSON $RSTAIN
  set_stain_vars $JSON $STAIN

  # Nothing to do for the reference stain
  if [[ $STAIN == $RSTAIN ]]; then return; fi

  # Get the number of chunks in the chunk mask
  local N_ROWS=$(c2d $REFSLIDE_CHUNK_MASK -dup -lstat | wc -l)
  local N_CHUNKS=$((N_ROWS-2))

  # Perform global rigid registration
  greedy -d 2 -a -dof 6 -i $REFSLIDE_SCALAR $SLIDE_SCALAR -gm $REFSLIDE_MASK \
    -ia-image-centers -m WNCC 4x4 -bg NaN -wncc-mask-dilate \
    -n 200x200x40x0x0 -search 20000 any 10 -o $SLIDE_TO_REF_GLOBAL_RIGID

  # Perform piecewise rigid registration
  if [[ $N_CHUNKS -gt 1 ]]; then
    multi_chunk_greedy -d 2 -cm $REFSLIDE_CHUNK_MASK -o $SLIDE_TO_REF_CHUNK_RIGID \
      -i $REFSLIDE_SCALAR $SLIDE_SCALAR -ia $SLIDE_TO_REF_GLOBAL_RIGID \
      -m WNCC 4x4 -bg NaN -wncc-mask-dilate -n 600x600x200x0 \
      -a -dof 6 -search 10000 10 5 -wreg 0.05
  else
    cp $SLIDE_TO_REF_GLOBAL_RIGID $(printf $SLIDE_TO_REF_CHUNK_RIGID 1)
  fi

  # Perform piecewise deformable registration
  multi_chunk_greedy -d 2 -cm $REFSLIDE_CHUNK_MASK \
    -sv -o $SLIDE_TO_REF_CHUNK_WARP -oinv $SLIDE_TO_REF_CHUNK_INVWARP \
    -it $SLIDE_TO_REF_CHUNK_RIGID \
    -i $REFSLIDE_SCALAR $SLIDE_SCALAR \
    -m WNCC 4x4 -bg NaN -wncc-mask-dilate -n 400x200x100x20 -s 0.6mm 0.1mm -e 0.25

  # Generate metric images
  multi_chunk_greedy -d 2 -cm $REFSLIDE_CHUNK_MASK \
    -metric -o $SLIDE_TO_REF_CHUNK_METRIC \
    -it $SLIDE_TO_REF_CHUNK_WARP $SLIDE_TO_REF_CHUNK_RIGID \
    -i $REFSLIDE_SCALAR $SLIDE_SCALAR \
    -m WNCC 4x4 -bg NaN -wncc-mask-dilate 

  # Apply the piecewise deformation to the input
  multi_chunk_greedy -d 2 -cm $REFSLIDE_CHUNK_MASK \
    -rf $REFSLIDE_SCALAR \
    -rb 255 -rm $SLIDE_RAW $SLIDE_TO_REF_RESLICE_DEFORM_RAW \
    -r $SLIDE_TO_REF_CHUNK_WARP $SLIDE_TO_REF_CHUNK_RIGID
}

# A single function that preprocesses the MRI and all slides
function preprocess_main()
{
  # User must pass a stain and a parameter json file as input
  JSON=${1?}

  # Set the variables
  set_experiment_vars $JSON

  # Process the MRI
  mri_preprocess $JSON

  # Process each of the stains
  for stain in $STAIN_LIST; do
    stain_preprocess $JSON $stain
  done
}

# This function performs registration between reference stain (as fixed image)
# and the MRI scan
function register_to_mri()
{
  # User must pass a stain and a parameter json file as input
  JSON=${1?}

  # Set the variables
  set_experiment_vars $JSON
  set_stain_vars $JSON $RSTAIN

  # Check the manual initialization
  if [[ ! -f $MRI_TO_REFSLIDE_MANUAL_WORKSPACE ]]; then
    echo "Missing manual registration workspace: $(check_file $MRI_TO_REFSLIDE_MANUAL_WORKSPACE)"
    exit 255
  fi

  # Extract the registration as a matrix
  itksnap-wt -i $MRI_TO_REFSLIDE_MANUAL_WORKSPACE \
    -lp 1 -props-get-transform \
    | grep '3>' | sed -e "s/3> //g" \
    > $MRI_TO_REFSLIDE_MANUAL

  # Perform affine registration based on the initial manual matrix
  greedy -d 3 -z -a -dof 7 \
    -i $REFSLIDE_MRILIKE $MRI -gm $REFSLIDE_MRILIKE_MASK \
    -ia $MRI_TO_REFSLIDE_MANUAL -n 100x100x40 \
    -m WNCC 2x2x0 -bg NaN -o $MRI_TO_REFSLIDE_RIGID

  greedy -d 3 -z -a -dof 12 \
    -i $REFSLIDE_MRILIKE $MRI -gm $REFSLIDE_MRILIKE_MASK \
    -ia $MRI_TO_REFSLIDE_RIGID -n 100x100x40 \
    -m WNCC 2x2x0 -bg NaN -o $MRI_TO_REFSLIDE_AFFINE

  # Perform deformable registration. The key flags here are -z (2D/3D)
  # registration, all smoothing operations are only performed in x,y 
  # and -ref-pad, which adds padding to the 2D slide before the 3D image
  # is resampled to its space via the -it command. Also the -WNCC flag
  # has 2x2x0 radius, so that the metric is computed correctly, in 2D.
  greedy -d 3 -z -ref-pad 0x0x2 \
    -i $REFSLIDE_MRILIKE $MRI -gm $REFSLIDE_MRILIKE_MASK \
    -it $MRI_TO_REFSLIDE_AFFINE -n 200x200x100 \
    -m WNCC 2x2x0 -bg NaN \
    -s 1.5mm 0.5mm -sv -e 0.5 \
    -o $MRI_TO_REFSLIDE_WARP 

  # Warp reference slice to the MRI
  greedy -d 3 \
    -rf $REFSLIDE_MRILIKE \
    -rm $MRI $MRI_TO_REFSLIDE_RESLICE_DEFORM \
    -r $MRI_TO_REFSLIDE_WARP $MRI_TO_REFSLIDE_AFFINE
  
  # Also reslice using affine tranform only, for comparison
  greedy -d 3 \
    -rf $REFSLIDE_MRILIKE \
    -rm $MRI $MRI_TO_REFSLIDE_RESLICE_AFFINE \
    -r $MRI_TO_REFSLIDE_AFFINE

  # Create a workspace to inspect the registration result
  itksnap-wt \
    -laa $REFSLIDE_MRILIKE -psn "Reference slide (MRI-like)" \
    -laa $MRI_TO_REFSLIDE_RESLICE_DEFORM -psn "MRI (deformed)" \
    -laa $MRI_TO_REFSLIDE_RESLICE_AFFINE -psn "MRI (affine)" \
    -o $MRI_TO_REFSLIDE_RESULT_WORKSPACE


}


function usage()
{
  echo "runme.sh : Dots project histology/MRI registration script"
  echo "Usage:"
  echo "  runme.sh [options] <function> [args]"
  echo "Options:"
  echo "  -d            : Turn on command echoing (for debugging)"
  echo "  -s <level>    : Skip expensive operations if results exist"
  echo "Primary functions:"
  echo "  check_inputs <json>                        : Check inputs for a given case"
  echo "  preprocess_main <json>                     : Run all preprocessing steps for a case,"
  echo "                                               calls mri_preprocess and stain_preprocess"
  echo "  mri_preprocess <json>                      : Run all preprocessing steps for a case"
  echo "  stain_preprocess <json> <stain>            : Preprocess a histology for given stain"
  echo "  register_to_mri <json>                     : Register reference stain to MRI"
  echo "  stain_register_to_reference <json> <stain> : Register stain to reference stain"
}

# Read the command-line options
while getopts "dhs:" opt; do
  case $opt in
    d) set -x -e;;
    h) usage; exit 0;;
    s) SKIPLEVEL=$OPTARG;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;
  esac
done

# Get remaining args
shift $((OPTIND - 1))

# No parameters? Show usage
if [[ "$#" -lt 1 ]]; then
  usage
  exit 255
fi

# Main entrypoint into script
COMMAND=$1
shift
$COMMAND "$@"





