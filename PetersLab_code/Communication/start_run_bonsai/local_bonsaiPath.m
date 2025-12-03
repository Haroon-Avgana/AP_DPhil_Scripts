function pathout = local_bonsaiPath(x64x86)

%%version of bonsai
if nargin <1; x64x86 = 64; end
switch x64x86
    case {'64', 64, 'x64', 1, '64bit'}
        bonsaiEXE = 'Bonsai.exe';
    case {32, '32', 'x86', 0, '32bit'}
        bonsaiEXE = 'Bonsai32.exe';
end

dirs = {'appdata', 'localappdata', 'programfiles', 'programfiles(x86)',...
    'C:\Software\Bonsai.Packages\Externals\Bonsai\Bonsai.Editor\bin\x86\Release',...
    'C:\Software\Bonsai.Packages\Externals\Bonsai\Bonsai.Editor\bin\x64\Release',...
    };
foundBonsai = 0;
dirIDX = 1;
while ~foundBonsai && (dirIDX <= length(dirs))
    if dirIDX<=4
        dir = getenv(dirs{dirIDX});
    else
        dir = dirs{dirIDX};
    end
    pathout = fullfile(dir,'Bonsai', bonsaiEXE);    
    foundBonsai = exist(pathout, 'file');
    dirIDX = dirIDX +1;
end

if ~foundBonsai
    warning('could not find bonsai executable, please insert it manually');
    [fname fpath] = uigetfile( '*.exe', 'Provide the path to Bonsai executable');
    pathout = fullfile(fpath, fname);
end

end


