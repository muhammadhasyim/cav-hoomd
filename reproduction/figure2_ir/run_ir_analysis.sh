#!/bin/bash

# Base directory where all numbered directories are located
BASE_DIR=$(pwd)

# Loop through directories (1-20)
for idx in {0..499}; do
    for dir in {0..20}; do
        # Check if directory exists
        if [ ! -d "$dir" ]; then
            echo "Directory $dir does not exist, skipping..."
            continue
        fi

        echo "Processing directory $dir..."
        #rm $dir/spectrum_hoomd_${idx}.txt
        #rm $dir/dacf_hoomd_${idx}.txt
        
        # Loop through indices (1-500)
    
        # Check if input file exists
        if [ ! -f "$dir/ir-${idx}.h5" ]; then
            continue
        fi
        
        # Check if output files already exist
        if [ -f "$dir/spectrum_hoomd_${idx}.txt" ] && [ -f "$dir/dacf_hoomd_${idx}.txt" ]; then
            echo "Output files for dir $dir, idx $idx already exist, skipping..."
            continue
        fi
        
        echo "Processing dir $dir, idx $idx..."
        python compute_irspectrum.py --dir "$BASE_DIR/$dir" --idx "$idx"
        
        # Check if the processing was successful
        if [ $? -ne 0 ]; then
            echo "Error processing dir $dir, idx $idx"
        fi
    done
done

echo "Processing complete!" 
