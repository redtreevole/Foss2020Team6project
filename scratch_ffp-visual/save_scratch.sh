## icommands
# https://cyverse-data-store-guide.readthedocs-hosted.com/en/latest/step2.html#icommands-installation-for-linux
# https://learning.cyverse.org/projects/cyverse-cyverse-reproducbility-tutorial/en/latest/step2.html#make-a-local-clone-of-your-github-repository

$SCRATCH_DIR="/scratch/reproducibility-tutorial" #target path
$USER_NAME=""

#copy whole /scratch folder

if [ ! -d "$SCRATCH_DIR" ]; then
  # Control will enter here if $DIRECTORY doesn't exist.
  #imkdir "$SCRATCH_DIR"
  mkdir "$SCRATCH_DIR"
  iget -rPV /iplant/home/jaejinchoi/reproducibility-tutorial/reproducibility-tutorial $SCRATCH_DIR
  #exact iplant branch can be seen from Discovery environment

fi


## save / put to data storage
#iput -rPV /scratch

iput -rPV $SCRATCH_DIR ./ #to iplant homefolder
