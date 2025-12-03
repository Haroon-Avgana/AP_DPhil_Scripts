# Set script path as current script path
$tapeScriptPath = $PSScriptRoot

# Display function
Write-Output 'Running daily data archive...' 

# Move raw data from data folder to to 'to_archive' folder
$serverToArchiveScript = Join-Path -Path $tapeScriptPath -ChildPath "server_to_archive_widefield.ps1"
Write-Output 'Moving raw data to archive folder...' 
& $serverToArchiveScript


# Move 'to_archive' contents to tape
$archiveToTapeScript = Join-Path -Path $tapeScriptPath -ChildPath "archive_to_tape.ps1"
Write-Output 'Moving archive folder contents to tape...'
& $archiveToTapeScript

