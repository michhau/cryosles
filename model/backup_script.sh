#!/bin/bash

# Get the current directory name
current_dir=$(basename "$PWD")

# Define backup directory name (one level up)
backup_dir="../bckp_${current_dir}"

# Check if backup directory already exists
if [ -d "$backup_dir" ]; then
    echo "Error: Backup directory '$backup_dir' already exists!"
    echo "Please remove or rename it before running this script."
    exit 1
fi

# Create backup directory
echo "Creating backup directory: $backup_dir"
mkdir -p "$backup_dir"

# Copy files (continue even if some don't exist)
echo "Copying files..."

# Copy text files
cp *.txt "$backup_dir/" 2>/dev/null && echo "  ✓ Copied *.txt" || echo "  ⚠ No *.txt files found"

# Copy .in files
cp *.in "$backup_dir/" 2>/dev/null && echo "  ✓ Copied *.in" || echo "  ⚠ No *.in files found"

# Copy launch script
cp launch_eiger.sh "$backup_dir/" 2>/dev/null && echo "  ✓ Copied launch_eiger.sh" || echo "  ⚠ launch_eiger.sh not found"

# Copy total_time.dat
cp total_time.dat "$backup_dir/" 2>/dev/null && echo "  ✓ Copied total_time.dat" || echo "  ⚠ total_time.dat not found"

# Copy src/ directory
if [ -d "src" ]; then
    cp -r src "$backup_dir/" && echo "  ✓ Copied src/"
else
    echo "  ⚠ src/ directory not found"
fi

# Copy PARAMETERS.py_0
cp PARAMETERS.py_0 "$backup_dir/" 2>/dev/null && echo "  ✓ Copied PARAMETERS.py_0" || echo "  ⚠ PARAMETERS.py_0 not found"

# Copy program.exe
cp program.exe "$backup_dir/" 2>/dev/null && echo "  ✓ Copied program.exe" || echo "  ⚠ program.exe not found"

# Copy log files
cp *.log "$backup_dir/" 2>/dev/null && echo "  ✓ Copied *.log" || echo "  ⚠ No *.log files found"

# Handle restart_files_0/ with filtering
if [ -d "restart_files_0" ]; then
    echo "Processing restart_files_0/..."
    
    # Find the largest xxx number in restart_xxx.cyyy files
    max_xxx=$(ls restart_files_0/restart_*.c* 2>/dev/null | \
              sed 's/.*restart_\([0-9]*\)\.c.*/\1/' | \
              sort -n | tail -1)
    
    if [ -n "$max_xxx" ]; then
        echo "  Found maximum timestep: $max_xxx"
        
        # Create restart_files_0 in backup
        mkdir -p "$backup_dir/restart_files_0"
        
        # Copy only files with the largest xxx
        cp restart_files_0/restart_${max_xxx}.c* "$backup_dir/restart_files_0/" 2>/dev/null
        
        num_files=$(ls "$backup_dir/restart_files_0/" 2>/dev/null | wc -l)
        echo "  ✓ Copied $num_files restart files with timestep $max_xxx"
    else
        echo "  ⚠ No restart files found in restart_files_0/"
    fi
else
    echo "  ⚠ restart_files_0/ directory not found"
fi

echo ""
echo "Backup completed successfully!"
echo "Backup location: $backup_dir"
