$todayDate = Get-Date -Format "yyyy-MM-dd"

$tapePath = 'T:\'
$outputFile = "\\qnap-ap001.dpag.ox.ac.uk\APlab\Archive\tape_contents\tape_contents_$($todayDate).csv"
                

Write-Output "Writing tape contents log: $outputFile"

# Initialize an array to hold the custom objects
$items = @()

# Get all items recursively
Get-ChildItem -Recurse $tapePath | ForEach-Object {
    # Create a custom object for each item
    $item = [PSCustomObject]@{
        Path = $_.FullName
        Name = $_.Name
        ParentFolder = $_.DirectoryName
        ItemType = if ($_.PSIsContainer) { 'Directory' } else { 'File' }
        SizeGB = if ($_.PSIsContainer) { 0 } else { [math]::Round($_.Length / 1GB, 4) }
        CreationTime = $_.CreationTime
    }
    # Add the custom object to the array
    $items += $item
}

# Export the array to CSV
$items | Export-Csv -Path $outputFile -NoTypeInformation