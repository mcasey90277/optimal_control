function outFile = saveFigure(hFig,fileName,opts)
%% Purpose:
%
%  Write ONE figure to an image file at a fixed resolution and with the
%  figure's own background preserved. Every plotting function in
%  coorbital.viz returns the figure it drew into, and every entry script that
%  wants those figures on disk needs the same three decisions made the same
%  way -- what resolution, what background, and what happens when the target
%  directory does not exist. Making them here means the shipped pictures from
%  BM/run_ballistic_target, HGV/run_target and anything added later are
%  comparable rather than each the product of one inline exportgraphics call
%  with whatever arguments were to hand.
%
%  THE BACKGROUND DEFAULTS TO 'current', NOT TO WHITE. exportgraphics is
%  perfectly willing to lay a white card behind a dark figure, and the globe
%  movie's frames and any dark-styled still are ruined by exactly that: the
%  starfield, the white legend text and the white annotation boxes all end up
%  invisible. 'current' takes the figure's own Color property, so a light
%  figure exports light and a dark figure exports dark without the caller
%  having to know which it has.
%
%  THE EXTENSION IS OPTIONAL and defaults to .png, because the callers hold a
%  path STEM -- one user-block entry per script from which several figure
%  names are derived by suffix -- and appending '.png' at each call site would
%  be the same string written once per figure.
%
%  A MISSING DIRECTORY IS CREATED rather than refused. The alternative was
%  considered and rejected: the entry scripts default their stem into tempdir,
%  where the directory always exists, so the only time this branch is reached
%  is when a user has deliberately pointed the stem somewhere of their own
%  choosing, and refusing that with "make the folder first" is a chore and not
%  a safety property. A directory that cannot be created still raises.
%
%  THE WRITE IS CHECKED. exportgraphics returns nothing, so a caller that only
%  looks for the file existing can be satisfied by a zero-byte stub left
%  behind by a failed write. The size on disk is read back and a non-positive
%  one raises, which is the cheapest available evidence that pixels actually
%  reached the file.
%
%% Inputs:
%
%  hFig             [1 x 1] figure              Figure to export. Must be a
%                                               valid, live figure handle;
%                                               'Visible','off' is fine and is
%                                               what the tests use
%
%  fileName         Char [1 x n] or string      Output path. An extension of
%                                               any kind is honoured as given;
%                                               with none, '.png' is appended.
%                                               A relative path is resolved
%                                               against the current directory
%                                               so the returned path is always
%                                               unambiguous
%
%  opts             Struct, optional            All fields optional:
%                                               Resolution      dpi of the
%                                                               written image
%                                                               (dots/inch);
%                                                               default 200
%                                               BackgroundColor exportgraphics
%                                                               background:
%                                                               'current'
%                                                               (default),
%                                                               'none', a
%                                                               colour name or
%                                                               a [1 x 3] RGB
%
%% Outputs:
%
%  outFile          Char [1 x m]                Absolute path actually
%                                               written, extension included.
%                                               Print this rather than the
%                                               caller's stem: it is what a
%                                               reader has to type to find the
%                                               picture
%
%% Notes:
%
%  THE FIGURE IS NOT CLOSED and nothing about it is modified. Ownership of a
%  figure belongs to whoever created it, which is the contract the rest of
%  coorbital.viz already keeps, and an export that quietly disposed of its
%  subject would make saving and looking at a figure mutually exclusive.
%
%  200 dpi is chosen as the resolution a reader can zoom into and a document
%  can embed without resampling, while staying small enough that three
%  figures per run are not a nuisance on disk. A 6.5 x 4.5 inch figure lands
%  near 1300 x 900 pixels. The number is an option, not a constant, precisely
%  because a figure destined for a slide and one destined for a paper want
%  different answers.
%
%  RE-RUNNING OVERWRITES. The file name is fully determined by the caller's
%  argument, so a second run of the same script replaces its own pictures
%  instead of accumulating a numbered pile of them.
%
%% Revision History:
%  Michael Casey                                                08/09/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo. A dark figure on purpose, because the background rule is the one
%% thing about this function a reader is most likely to get wrong, and a demo
%% that wrote a white-on-white picture would be evidence for the wrong claim:
if nargin == 0
            hDemo  = figure('Color',[0.05 0.05 0.08],'Visible','off', ...
                            'Units','inches','Position',[1 1 6 4], ...
                            'Name','coorbital.viz.saveFigure self-demo', ...
                            'NumberTitle','off');
             hAxD  = axes('Parent',hDemo,'Color',[0.05 0.05 0.08], ...
                          'XColor',[1 1 1],'YColor',[1 1 1]);
              tDem = linspace(0,4.*pi,400);
            line(hAxD,tDem,sin(tDem).*exp(-tDem./10), ...
                 'Color',[0.2 1 0.4],'LineWidth',2);
            grid(hAxD,'on');
            xlabel(hAxD,'t (s)');
            ylabel(hAxD,'amplitude (-)');
            title(hAxD,'coorbital.viz.saveFigure self-demo','Color',[1 1 1]);
           outFile = coorbital.viz.saveFigure(hDemo, ...
                         fullfile(tempdir,'saveFigure_selfdemo'));
            close(hDemo);
            fprintf('  saveFigure self-demo wrote\n    %s\n',outFile);
    return;
end

%% Options, defaulted before anything touches the disk:
    if nargin < 3 || isempty(opts)
              opts = struct();
    end
    assert(isstruct(opts),'opts must be a struct of options.');
               res = vizOption(opts,'Resolution',200);
            bgSpec = vizOption(opts,'BackgroundColor','current');
    if ~(isscalar(res) && isnumeric(res) && isfinite(res) && res > 0)
        error('coorbital:saveFigure:badResolution', ...
              ['Resolution must be a positive finite scalar in dots per inch; ' ...
               '%s given.'],mat2str(res));
    end

%% The figure. Checked before the path is built, so a caller who passed an
%% axes or a deleted handle is told that and not told about a directory:
    if ~(isscalar(hFig) && isgraphics(hFig,'figure'))
        error('coorbital:saveFigure:badFigure', ...
              ['hFig must be a single valid figure handle. Every plotting ' ...
               'function in coorbital.viz returns one as its first output; ' ...
               'an axes handle, a deleted figure or [] will not do.']);
    end

%% The path. A string is taken as a char, an absent extension becomes .png,
%% and a relative path is made absolute so that the value handed back names
%% the file no matter where the caller's directory goes next:
    if isstring(fileName) && isscalar(fileName)
          fileName = char(fileName);
    end
    if ~(ischar(fileName) && ~isempty(fileName) && size(fileName,1) == 1)
        error('coorbital:saveFigure:badFileName', ...
              ['fileName must be a non-empty character row vector or a ' ...
               'scalar string naming the output file.']);
    end
    [fDir,fBase,fExt] = fileparts(fileName);
    if isempty(fExt)
              fExt = '.png';
    end
    if isempty(fBase)
        error('coorbital:saveFigure:badFileName', ...
              '"%s" has no file name, only a directory.',fileName);
    end
    if isempty(fDir)
              fDir = pwd;
    elseif ~isAbsPath(fDir)
              fDir = fullfile(pwd,fDir);
    end
           outFile = fullfile(fDir,[fBase fExt]);

%% The directory, created when it is not there. mkdir's own message is kept in
%% the raised error, since "permission denied" and "name too long" are
%% different problems and only it knows which one happened:
    if ~isfolder(fDir)
        [okDir,msgDir] = mkdir(fDir);
        if ~okDir
            error('coorbital:saveFigure:dirFailed', ...
                  'Cannot create the output directory "%s": %s',fDir,msgDir);
        end
    end

%% The export itself. Wrapped so the failure names the path: exportgraphics
%% reports the underlying cause but not always the file it was trying, and the
%% caller supplied a stem rather than this name:
    try
        exportgraphics(hFig,outFile, ...
                       'Resolution',res, ...
                       'BackgroundColor',bgSpec);
    catch errExp
        error('coorbital:saveFigure:writeFailed', ...
              'Cannot write the figure to "%s": %s',outFile,errExp.message);
    end

%% Proof that pixels arrived. A file can exist and be empty, which is what a
%% write interrupted part way through leaves behind, and that state has fooled
%% this project before:
              dInf = dir(outFile);
    if isempty(dInf) || dInf(1).bytes <= 0
        error('coorbital:saveFigure:emptyFile', ...
              ['The export left no bytes at "%s". The figure was reached and ' ...
               'exportgraphics did not raise, so treat this as a failed write ' ...
               'rather than a bad figure.'],outFile);
    end
end

function tf = isAbsPath(pathStr)
%% Purpose:
%
%  Report whether a directory string is already absolute. Factored out because
%  the answer differs by platform and the test for it reads badly inline.
%
%% Inputs:
%
%  pathStr          Char [1 x n]                Directory part of a path
%
%% Outputs:
%
%  tf               [1 x 1] logical             True when the path is absolute
%
%% Revision History:
%  Michael Casey                                                08/09/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    if ispc
                tf = ~isempty(regexp(pathStr,'^([A-Za-z]:[\\/]|\\\\)','once'));
    else
                tf = strncmp(pathStr,'/',1);
    end
end
