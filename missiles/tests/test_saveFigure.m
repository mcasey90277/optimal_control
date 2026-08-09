function test_saveFigure()
%% Purpose:
%
%  Verify coorbital.viz.saveFigure: that it writes a real image and not an
%  empty stub, that the RESOLUTION it claims is the resolution the file has,
%  that the BACKGROUND is the figure's own and not a white card laid behind
%  it, that a missing directory is created, that a re-run overwrites rather
%  than accumulates, that the figure it was given survives untouched, and that
%  every refusal path raises with its own identifier.
%
%  THE RESOLUTION IS MEASURED, NOT TRUSTED. exportgraphics crops a figure to
%  its content, so the pixel count is not the figure's nominal size times the
%  dots per inch and cannot be predicted from the figure alone. What is fixed
%  is the CONTENT SIZE IN INCHES, which does not depend on the resolution, so
%  the same figure is exported twice -- once at an explicit 100 dpi, which
%  measures the content, and once at the default -- and the default must come
%  out at twice the pixels each way. A default changed from 200, or a
%  Resolution never forwarded and left at the exportgraphics default of 150,
%  both miss that by ten times the crop rounding.
%
%  THE BACKGROUND IS MEASURED THE SAME WAY, by reading the corner pixel of a
%  BLACK figure back off the disk. That is the check that bites the mistake
%  worth guarding against: exportgraphics will happily white-card a dark
%  figure, and the resulting picture is a white rectangle with invisible white
%  text on it. Existence and file size cannot see that; a pixel can.
%
%  THE WRITE IS PROVEN BY CONTENT, not by dir(). A function that returned a
%  plausible path and wrote nothing would satisfy every existence check ever
%  written; imread of a file that is not there, or is empty, throws, and the
%  image dimensions asserted afterwards are the evidence that the export ran.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none                                         Throws on any failed assertion
%
%% Revision History:
%  Michael Casey                                                08/09/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

             nFig0 = numel(findall(groot,'Type','figure'));

%% A scratch directory of this run's own, removed at the end, so the test
%% neither reads a file an earlier run left nor leaves one for the next:
            outDir = [tempname '_saveFigure'];
    assert(~isfolder(outDir),'scratch directory %s already exists.',outDir);
            mkdir(outDir);
             clean = onCleanup(@() rmdir(outDir,'s'));

%% The subject. Sized in INCHES and coloured BLACK on purpose: the size makes
%% the resolution check an exact one and the colour makes the background check
%% a real one:
            wInch  = 4;
            hInch  = 3;
             hFig  = figure('Visible','off','Color',[0 0 0], ...
                            'Units','inches','Position',[1 1 wInch hInch]);
              hAx  = axes('Parent',hFig,'Color',[0 0 0], ...
                          'XColor',[1 1 1],'YColor',[1 1 1]);
            line(hAx,linspace(0,1,50),linspace(0,1,50).^2, ...
                 'Color',[0 1 0],'LineWidth',2);

%% 1. A plain save. The extension is not given, so .png must be appended, and
%% the returned path must be absolute -- a caller prints it and a reader has
%% to be able to find the file from it:
              stem = fullfile(outDir,'plain');
           outPlan = coorbital.viz.saveFigure(hFig,stem);
    assert(strcmp(outPlan,[stem '.png']), ...
        'expected "%s.png" back, got "%s".',stem,outPlan);
    assert(isfile(outPlan),'saveFigure returned "%s" but wrote no file.',outPlan);
              dPln = dir(outPlan);
    assert(dPln.bytes > 0,'saveFigure left a zero-byte file at %s.',outPlan);

%% 2. The image is a real RGB picture and not a stub. imread throws on an
%% absent or empty file, so reaching the next line is already evidence:
           imgPlan = imread(outPlan);
           [nRow,nCol,nChn] = size(imgPlan);
    assert(nChn == 3,'expected an RGB PNG, got %d channel(s).',nChn);
    assert(nRow > 100 && nCol > 100, ...
        'the export is %d x %d pixels, which is not a rendered figure.',nCol,nRow);

%% 3. The background, read off the corner pixel. The figure is black, so an
%% export that kept the figure's own colour has a black corner and one that
%% laid a white card behind it has a white one. Nothing else in this test can
%% tell those two files apart:
            corner = double(squeeze(imgPlan(1,1,:)))';
    assert(all(corner <= 12), ...
        ['background not preserved: the corner of a BLACK figure came back ' ...
         'as RGB [%g %g %g]. exportgraphics must be given ' ...
         '''BackgroundColor'',''current''.'],corner(1),corner(2),corner(3));

%% 4. THE RESOLUTION, measured against a second export of the SAME figure at a
%% KNOWN 100 dpi. exportgraphics crops a figure to its content, so the written
%% pixel count is not the figure's nominal 4 x 3 inches times the resolution;
%% what IS exact is that the CONTENT has one size in inches whatever the dots
%% per inch, so the 100 dpi image measures it and the shipped default must
%% then come out at 200/100 times as many pixels each way. A default silently
%% changed, or a Resolution never forwarded and left at the exportgraphics
%% default of 150 dpi, both miss this by far more than the crop rounding:
           outHalf = coorbital.viz.saveFigure(hFig,fullfile(outDir,'half.png'), ...
                                              struct('Resolution',100));
           imgHalf = imread(outHalf);
            resDef = 200;
            wCntIn = size(imgHalf,2)./100;
            hCntIn = size(imgHalf,1)./100;
    assert(abs(wCntIn - wInch) < 1 && abs(hCntIn - hInch) < 1, ...
        ['the 100 dpi export measures the content at %.3f x %.3f inches ' ...
         'against a %g x %g inch figure; the crop is not what this test ' ...
         'assumes.'],wCntIn,hCntIn,wInch,hInch);
            wWantP = wCntIn.*resDef;
            hWantP = hCntIn.*resDef;
    assert(abs(nCol - wWantP) <= 0.02.*wWantP && ...
           abs(nRow - hWantP) <= 0.02.*hWantP, ...
        ['default resolution wrong: the default export is %d x %d pixels of ' ...
         'content that the 100 dpi export measures at %.3f x %.3f inches, ' ...
         'which is %.1f dpi across and %.1f down, not %g.'], ...
        nCol,nRow,wCntIn,hCntIn,nCol./wCntIn,nRow./hCntIn,resDef);

%% 5. The background is an option too. Asking for white on the same black
%% figure must actually produce white, which proves the option is forwarded
%% rather than that 'current' happens to be exportgraphics' own default:
           outWhit = coorbital.viz.saveFigure(hFig,fullfile(outDir,'white.png'), ...
                                              struct('BackgroundColor','white'));
           imgWhit = imread(outWhit);
            cornWh = double(squeeze(imgWhit(1,1,:)))';
    assert(all(cornWh >= 243), ...
        ['BackgroundColor is not forwarded: asking for white gave a corner of ' ...
         'RGB [%g %g %g].'],cornWh(1),cornWh(2),cornWh(3));

%% 6. A missing directory is CREATED, two levels deep, because a user who
%% points the stem at a folder of their own should not have to make it first:
            deepDr = fullfile(outDir,'made','up');
    assert(~isfolder(deepDr),'the deep directory exists before the call.');
           outDeep = coorbital.viz.saveFigure(hFig,fullfile(deepDr,'deep'));
    assert(isfolder(deepDr),'saveFigure did not create the missing directory.');
    assert(isfile(outDeep),'nothing written into the created directory.');

%% 7. A re-run OVERWRITES. Same path in, same path out, one file on disk --
%% the property that lets an entry script be run twice without piling up a
%% numbered heap of pictures:
           outAgan = coorbital.viz.saveFigure(hFig,stem);
    assert(strcmp(outAgan,outPlan), ...
        'the second save went to "%s", not back to "%s".',outAgan,outPlan);
              nPng = numel(dir(fullfile(outDir,'plain*.png')));
    assert(nPng == 1,'re-running left %d files matching plain*.png.',nPng);

%% 8. A relative path is made absolute against the current directory. Checked
%% by going there, saving with a bare name, and requiring the answer to name
%% the file that appeared:
              cwd0 = pwd;
           backCwd = onCleanup(@() cd(cwd0));
            cd(outDir);
            outRel = coorbital.viz.saveFigure(hFig,'relative');
            clear backCwd;
    assert(~isempty(regexp(outRel,'^([/\\]|[A-Za-z]:)','once')), ...
        'a relative name gave a relative answer: "%s".',outRel);
    assert(isfile(fullfile(outDir,'relative.png')), ...
        'the relative save did not land in the working directory.');

%% 9. THE FIGURE IS NOT TOUCHED. Ownership belongs to the caller, and an
%% export that closed or recoloured its subject would make saving a figure and
%% then looking at it mutually exclusive:
    assert(isgraphics(hFig,'figure'),'saveFigure closed the figure it was given.');
    assert(isequal(get(hFig,'Color'),[0 0 0]), ...
        'saveFigure left the figure Color at %s.',mat2str(get(hFig,'Color')));
    assert(strcmp(get(hFig,'Visible'),'off'), ...
        'saveFigure made an invisible figure visible.');

%% 10. The refusals, each by its own identifier so a test cannot be satisfied
%% by the wrong error. An axes is not a figure; a deleted handle is not a
%% figure; [] is not a figure:
    assert(strcmp(idOf(@() coorbital.viz.saveFigure(hAx,fullfile(outDir,'ax'))), ...
                  'coorbital:saveFigure:badFigure'), ...
        'an axes handle was accepted as a figure.');
    assert(strcmp(idOf(@() coorbital.viz.saveFigure([],fullfile(outDir,'nul'))), ...
                  'coorbital:saveFigure:badFigure'), ...
        '[] was accepted as a figure.');
            hDead  = figure('Visible','off');
            close(hDead);
    assert(strcmp(idOf(@() coorbital.viz.saveFigure(hDead,fullfile(outDir,'dead'))), ...
                  'coorbital:saveFigure:badFigure'), ...
        'a deleted figure handle was accepted.');

%% 11. Bad names and bad options:
    assert(strcmp(idOf(@() coorbital.viz.saveFigure(hFig,'')), ...
                  'coorbital:saveFigure:badFileName'), ...
        'an empty file name was accepted.');
    assert(strcmp(idOf(@() coorbital.viz.saveFigure(hFig,42)), ...
                  'coorbital:saveFigure:badFileName'), ...
        'a numeric file name was accepted.');
    assert(strcmp(idOf(@() coorbital.viz.saveFigure(hFig,[outDir filesep])), ...
                  'coorbital:saveFigure:badFileName'), ...
        'a directory with no file name was accepted.');
    assert(strcmp(idOf(@() coorbital.viz.saveFigure(hFig,fullfile(outDir,'r'), ...
                                                    struct('Resolution',-5))), ...
                  'coorbital:saveFigure:badResolution'), ...
        'a negative resolution was accepted.');

%% 12. A directory that CANNOT be created. A regular file is put where the
%% directory component would have to go, so mkdir has no way to succeed:
            blockF = fullfile(outDir,'blocker');
             fidBk = fopen(blockF,'w');
    assert(fidBk > 0,'could not create the blocking file.');
            fclose(fidBk);
    assert(strcmp(idOf(@() coorbital.viz.saveFigure(hFig, ...
                          fullfile(blockF,'sub','x.png'))), ...
                  'coorbital:saveFigure:dirFailed'), ...
        'a directory under a regular file was not refused.');

%% 13. A path that cannot be WRITTEN, the directory being fine. An existing
%% DIRECTORY is put at the exact output path, which exportgraphics cannot
%% open, and the failure must be reported as a write failure naming the file:
            dirAsF = fullfile(outDir,'occupied.png');
            mkdir(dirAsF);
    assert(strcmp(idOf(@() coorbital.viz.saveFigure(hFig,dirAsF)), ...
                  'coorbital:saveFigure:writeFailed'), ...
        'writing over an existing directory was not refused.');

%% 14. The self-demo runs, writes, and leaves no figure behind. A header claim
%% that "nargin == 0 demonstrates it" is worth exactly as much as a test that
%% calls it:
            nFigDm = numel(findall(groot,'Type','figure'));
           demFile = '';
            outDem = evalc('demFile = coorbital.viz.saveFigure();');
    assert(ischar(demFile) && isfile(demFile), ...
        'the self-demo returned "%s", which is not a file.',demFile);
    assert(contains(outDem,demFile), ...
        'the self-demo did not print the path it wrote.');
    assert(numel(findall(groot,'Type','figure')) == nFigDm, ...
        'the self-demo left %d figure(s) open.', ...
        numel(findall(groot,'Type','figure')) - nFigDm);
            delete(demFile);

%% Leave the figure count as it was found:
            close(hFig);
    assert(numel(findall(groot,'Type','figure')) == nFig0, ...
        'test_saveFigure leaked %d figure(s).', ...
        numel(findall(groot,'Type','figure')) - nFig0);
            clear clean;
end

function idStr = idOf(fcn)
%% Purpose:
%
%  Run a function handle that is expected to raise and return the identifier
%  of the error it raised, or '' when it did not raise at all. Factored out
%  because every refusal check above is the same four lines otherwise.
%
%% Inputs:
%
%  fcn              Function handle             Call expected to raise
%
%% Outputs:
%
%  idStr            Char [1 x n]                Error identifier, or '' if the
%                                               call returned normally
%
%% Revision History:
%  Michael Casey                                                08/09/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

             idStr = '';
    try
        fcn();
    catch err
             idStr = err.identifier;
    end
end
