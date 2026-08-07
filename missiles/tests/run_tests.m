function run_tests()
%% Purpose:
%
%  Run every test_*.m in this directory and report pass/fail. Each test file
%  is a script or function that throws on failure and returns silently on
%  success.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none                                         Prints a summary; sets a
%                                               non-zero exit code on failure
%                                               when run under -batch
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Locate this directory and the library root:
              here = fileparts(mfilename('fullpath'));
              root = fullfile(here,'..');
    addpath(root);
    addpath(fullfile(root,'HGV'));

%% Discover the test files:
             files = dir(fullfile(here,'test_*.m'));
                nT = numel(files);
            passed = 0;
            failed = 0;
          failures = cell(nT,1);

%% Run each test in turn, recording the message of any that throws:
    for k = 1:nT
             name = erase(files(k).name,'.m');
        try
            feval(name);
            passed = passed + 1;
            fprintf('  PASS  %s\n',name);
        catch err
            failed = failed + 1;
            failures{k} = sprintf('%s: %s',name,err.message);
            fprintf('  FAIL  %s\n',name);
        end
    end
          failures = failures(~cellfun(@isempty,failures));

%% Report:
    fprintf('\n%d passed, %d failed\n',passed,failed);
    for k = 1:numel(failures)
        fprintf('\n----- %s\n',failures{k});
    end
    if failed > 0
        error('missiles:testsFailed','%d test(s) failed',failed);
    end
end
