#!/bin/bash
# Check if the script is being run from the "slamv2" folder
if [[ "$(basename "$(pwd)")" != "slamv2" ]]; then
  echo "Please run the script from the 'slamv2' folder."
  exit 1
fi
# Access the first command-line argument
folder="$1"
# Check if the argument is provided and if the folder exists
if [ -z "$folder" ]; then
  echo "No folder location provided."
  exit 1
elif [ ! -d "$folder" ]; then
  echo "Folder does not exist: $folder"
  exit 1
elif [ "$(ls -A "$folder")" ]; then
  echo "Folder is not empty: $folder"
  exit 1
fi
# Use the folder location in your script
echo "Copying to: $folder"
cp ./build*/*png $folder
cp ./data/shared_data/cur_env.ttt $folder
cp -r ./data/shared_data/scans/ $folder
