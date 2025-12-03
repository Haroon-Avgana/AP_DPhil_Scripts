# Set folders
$data_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Data'
$archive_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Archive\to_archive'

# Create move log file
$datestring = (Get-Date).ToString("yyyy-MM-dd")
$move_log_path = "\\qnap-ap001.dpag.ox.ac.uk\APlab\Archive\logs\move_to_archive_log\move_log_$datestring.csv"

# Set widefield file formats
# ('widefield_(rec time)_data.bin is from plab.rig.widefield)
# (old/unused: tif/dcimg is from HCImage)
$raw_wf_format = 'widefield_*_data.bin'

# Find all raw widefield files
Write-Output("Looking for raw widefield files on server to archive...")
$raw_wf_files = Get-ChildItem -Path $data_path -Filter $raw_wf_format -Recurse | Where-Object { $_.FullName -like "*\widefield\*" }

# Select files by relative date (2 weeks before today)
$end_date = (Get-Date).AddDays(-14)
$archive_file_idx = $raw_wf_files | Where-Object { $_.LastWriteTime -lt $end_date }

# Print size of all/selected raw files
$raw_size_total = ($raw_wf_files | Measure-Object -Property Length -Sum).Sum / 1e12
Write-Output ("Total raw size: {0:N2} TB" -f $raw_size_total)

$raw_size_archive = ($archive_file_idx | Measure-Object -Property Length -Sum).Sum / 1e12
Write-Output ("To-archive raw size: {0:N2} TB" -f $raw_size_archive)

# Initialize an array to hold the source and destination paths
$file_movements = @()

# Move files to archive folder
foreach ($i in 0..($archive_file_idx.Count - 1)) {
    $move_file = $archive_file_idx[$i]

    # Get data folder
    $source_folder = $move_file.DirectoryName

    # Pull out subfolders (e.g. animal, day, 'widefield')
    $source_subfolder = $source_folder.Replace("$($data_path)\", "").Split('\')

    # Construct destination folder from subfolders
    $destination_folder = Join-Path -Path $archive_path -ChildPath ($source_subfolder -join '_')
    if (-not (Test-Path -Path $destination_folder -PathType Container)) {
        New-Item -ItemType Directory -Path $destination_folder | Out-Null
    }

    $curr_source = $move_file.FullName
    $curr_destination = Join-Path -Path $destination_folder -ChildPath $move_file.Name

    # Move file (if it exists), print status
    if (Test-Path -Path $curr_source -PathType Leaf) {
        Write-Output ("Moving: {0} --> {1}" -f $curr_source, $curr_destination)
        Move-Item -Path $curr_source -Destination $curr_destination  
        
        # Add source and destination to the array
        $file_movements += [PSCustomObject]@{
            Source = $curr_source
            Destination = $curr_destination   
        }
    }
}

# Write log of what was moved (source and destination)
$file_movements | Export-Csv -Path $move_log_path -NoTypeInformation
