function run_smooth(adapt_sigma, sigma_frac, cad_order, out_order, rlam_override)
%RUN_SMOOTH Experimental fixed-CAD two-stage smoothing of dish_fd023_m.
%
%   run_smooth                       % defaults: as=1, sigma=wall/4, CAD 12, p4
%   run_smooth(1, 0.50, 12, 4)       % sigma = wall/2
%   run_smooth(0, 0.25, 12, 4)       % uniform mollifier
%   run_smooth(1, 0.25, 16, 4)       % finer CAD quadrature for phi
%   run_smooth(1, 0.25, 12, 4, 50)   % force rlam = 50, ignoring sigma_frac
%
% Arguments
%   adapt_sigma   0 = uniform mollifier  (sigma = mean(D)/rlam, grad_sigma = 0)
%                 1 = longest-side adaptive (sigma_j = D_j/rlam)   [default 1]
%   sigma_frac    target sigma as a fraction of the 0.01 m wall    [default 0.25]
%   cad_order     order the fixed CAD skeleton is oversampled to   [default 12]
%   out_order     order of the output smoothed surface             [default 4]
%   rlam_override optional explicit rlam; overrides sigma_frac
%
% WHY THESE DEFAULTS (established on dish_fd021_m_rimfree / _rimcap015):
%   1. The dish is a CLOSED SHELL 0.01 m thick. sigma must stay <= wall/2 or
%      the two paraboloid sheets merge and phi never reaches 1/2. This is
%      physics, converged to 1e-13 in the CAD quadrature - not a resolution
%      artifact.
%   2. The CAD skeleton must be oversampled to resolve that sigma. Measured
%      phi error at a blocking node, max|phi - phi(CAD 20)|:
%          CAD  8 -> 3.8e-02     CAD 12 -> 9.8e-03     CAD 16 -> 7.4e-04
%      Order 4 (the raw export) cannot see the shell at all.
%   Oversampling adds quadrature points, not geometric fidelity: the patches
%   remain the order-4 polynomial, so the surface itself is unchanged.
%
% Results are written under <mesh>/smooth/<tag>/. Repeating a configuration
% replaces its result files. Each run uses a separate scratch directory.

if nargin < 1 || isempty(adapt_sigma), adapt_sigma = 1;    end
if nargin < 2 || isempty(sigma_frac),  sigma_frac  = 0.25; end
if nargin < 3 || isempty(cad_order),   cad_order   = 12;   end
if nargin < 4 || isempty(out_order),   out_order   = 4;    end
if nargin < 5, rlam_override = []; end
if ischar(adapt_sigma)||isstring(adapt_sigma), adapt_sigma = str2double(adapt_sigma); end
if ischar(sigma_frac) ||isstring(sigma_frac),  sigma_frac  = str2double(sigma_frac);  end
if ischar(cad_order)  ||isstring(cad_order),   cad_order   = str2double(cad_order);   end
if ischar(out_order)  ||isstring(out_order),   out_order   = str2double(out_order);   end
if ischar(rlam_override)||isstring(rlam_override), rlam_override = str2double(rlam_override); end

validateattributes(adapt_sigma,{'numeric'},{'real','scalar','finite','integer','>=',0,'<=',1});
validateattributes(sigma_frac,{'numeric'},{'real','scalar','finite','positive'});
validateattributes(cad_order,{'numeric'},{'real','scalar','finite','integer','positive'});
validateattributes(out_order,{'numeric'},{'real','scalar','finite','integer','positive'});
if ~isempty(rlam_override)
    validateattributes(rlam_override,{'numeric'},{'real','scalar','finite','positive'});
end

stem    = 'dish_fd023_m';
[meshDir,repoRoot] = dish_smoothing_paths(stem);
go3In   = fullfile(meshDir,[stem '_order4.go3']);
gidIn   = fullfile(meshDir,[stem '_scaffold.gidmsh']);
assert(isfile(go3In),'missing %s',go3In);
assert(isfile(gidIn),'missing %s',gidIn);
assert(isfile(fullfile(repoRoot,'matlab','setup_surfsmooth3d.m')), ...
    'Missing SurfSmooth3D; set SURFSMOOTH3D_ROOT to its repository root.');
addpath(fullfile(repoRoot,'matlab')); setup_surfsmooth3d();

% From ../cad/dish_fd023_m.json (units m, f/D 0.23)
areaExact = 71.25967667149959;
wall      = 0.01;

tag = sprintf('as%d_sig%03.0fw_cad%02d_p%02d', ...
    adapt_sigma, sigma_frac*100, cad_order, out_order);
if ~isempty(rlam_override)
    tag = [tag sprintf('_rlam%.17g',rlam_override)];
end
outDir = fullfile(meshDir,'smooth',tag);
if ~isfolder(outDir), mkdir(outDir); end
workRoot = fullfile(meshDir,'smooth','work');
if ~isfolder(workRoot), mkdir(workRoot); end
workDir = tempname(workRoot);
mkdir(workDir);
scaffold = fullfile(workDir,'scaffold.gidmsh');
copyfile(gidIn,scaffold);

diary(fullfile(outDir,'run.log')); diary on;
cleanupDiary = onCleanup(@() diary('off'));
fprintf('\n===== %s : %s =====\n',stem,tag);
fprintf('started %s on %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),getenv('HOSTNAME'));
fprintf('MATLAB %s, maxNumCompThreads %d\n\n',version,maxNumCompThreads);

%% Scaffold statistics
[nodes,tris] = read_gid(scaffold);
P = reshape(nodes(:,tris),3,3,[]);
e = [vecnorm(squeeze(P(:,2,:)-P(:,1,:))); ...
     vecnorm(squeeze(P(:,3,:)-P(:,2,:))); ...
     vecnorm(squeeze(P(:,1,:)-P(:,3,:)))];
D = max(e,[],1);
Dmax = max(D); Dmin = min(D); Dmean = mean(D);
fprintf('scaffold: %d triangles\n',size(tris,2));
fprintf('  longest side  min %.5g  mean %.5g  max %.5g m   (spread %.1fx)\n', ...
    Dmin,Dmean,Dmax,Dmax/Dmin);
fprintf('  wall %.4g m -> h_max/wall = %.2f\n\n',wall,Dmax/wall);

sigmaTarget = wall*sigma_frac;
if ~isempty(rlam_override) && isfinite(rlam_override)
    rlam = rlam_override;
    if adapt_sigma == 0, sigmaTarget = Dmean/rlam; else, sigmaTarget = Dmax/rlam; end
    fprintf('rlam forced to %.6g -> sigma %s = %.5g m (%.3f x wall)\n\n', ...
        rlam,modeword(adapt_sigma),sigmaTarget,sigmaTarget/wall);
else
    if adapt_sigma == 0, rlam = Dmean/sigmaTarget; else, rlam = Dmax/sigmaTarget; end
    fprintf('adapt_sigma %d -> rlam = %s/sigma = %.6g  (sigma %s = %.5g m = %.3f x wall)\n\n', ...
        adapt_sigma,ternary(adapt_sigma==0,'mean(D)','max(D)'),rlam, ...
        modeword(adapt_sigma),sigmaTarget,sigmaTarget/wall);
end
if sigmaTarget > wall/2
    warning('run_smooth:sigmaTooLarge', ...
        ['sigma = %.5g m exceeds wall/2 = %.5g m. The two shell sheets will ', ...
         'merge and phi will never reach 1/2; expect failure.'],sigmaTarget,wall/2);
end

%% Fixed CAD source, oversampled for phi
S4 = surfsmooth3d.surfer.load_from_file(go3In);
fprintf('CAD source: p=%d, npatches=%d, npts=%d\n', ...
    double(S4.norders(1)),S4.npatches,S4.npts);
areaCad = sum(S4.wts);
fprintf('  quadrature area %.10f m^2 (rel err vs analytic %.3e)\n',areaCad, ...
    abs(areaCad-areaExact)/areaExact);
rep = surfsmooth3d.edgepreserve.validity_report(S4,'fixed CAD source');
surfsmooth3d.edgepreserve.print_validity_report(rep);

t0 = tic;
if cad_order == double(S4.norders(1)), Sp = S4; else, Sp = S4.oversample(cad_order); end
cadFile = fullfile(workDir,sprintf('cad_p%02d.txt',cad_order));
surfsmooth3d.edgepreserve.write_point_source(Sp,cadFile);
dd = dir(cadFile);
fprintf('\nCAD oversampled to order %d: %d pts, area %.10f, %.0f MB, %.1f s\n\n', ...
    cad_order,Sp.npts,sum(Sp.wts),dd.bytes/1e6,toc(t0));
nCadPts = Sp.npts; clear Sp;

%% Solve
opts = struct('filetype',3,'fcad',cadFile,'nquad',min(cad_order,20),'nrefine',0, ...
    'rlam',rlam,'adapt_sigma',adapt_sigma,'two_stage_smoother',true, ...
    'newton_recovery',true);
fprintf('opts: '); disp(opts);

r = struct('tag',string(tag),'adapt_sigma',adapt_sigma,'sigma_frac',sigma_frac, ...
    'rlam',rlam,'sigma',sigmaTarget,'sigma_over_wall',sigmaTarget/wall, ...
    'cad_order',cad_order,'cad_pts',nCadPts,'out_order',out_order, ...
    'status',"",'identifier',"",'fail_stage',"",'fail_reason',"", ...
    'deferred',nan,'recovered',nan,'unresolved',nan, ...
    'validity',"",'area',nan,'rel_area_err',nan,'max_disp',nan, ...
    'bad_nodes',nan,'ncd_min',nan,'jac_min',nan,'seconds',nan,'file',"");

tsolve = tic;
try
    [S,info] = surfsmooth3d.multiscale_mesher(scaffold,out_order,opts);
    r.seconds = toc(tsolve); Sout = S{1};
    r.status = "solved";
    if isstruct(info) && isfield(info,'recovery')
        rc = info.recovery;
        r.deferred = rc.deferred; r.recovered = rc.recovered; r.unresolved = rc.unresolved;
        fprintf('\nnewton recovery: deferred %d, recovered %d, unresolved %d\n', ...
            rc.deferred,rc.recovered,rc.unresolved);
        if isfield(rc,'report_file') && ~isempty(rc.report_file) && isfile(rc.report_file)
            copyfile(rc.report_file,fullfile(outDir,'newton_recovery.txt'));
            fprintf('recovery report -> %s\n',fullfile(outDir,'newton_recovery.txt'));
        end
    end
    r.area = sum(Sout.wts);
    r.rel_area_err = abs(r.area-areaExact)/areaExact;
    ref = [];
    if Sout.npts == S4.npts && Sout.npatches == S4.npatches, ref = S4; end
    rp = surfsmooth3d.edgepreserve.validity_report(Sout,sprintf('smoothed %s',tag),ref);
    surfsmooth3d.edgepreserve.print_validity_report(rp);
    if rp.is_valid && ~rp.has_warnings, r.validity = "PASS";
    elseif rp.is_valid, r.validity = "WARN"; else, r.validity = "FAIL"; end
    r.bad_nodes = rp.stats.bad_node_count;
    r.ncd_min   = rp.stats.normal_cross_dot_min;
    r.jac_min   = rp.stats.jac_min;
    if rp.has_reference
        r.max_disp = rp.stats.displacement_max;
        fprintf('  max displacement from CAD %.6e m (%.3f x wall)\n',r.max_disp,r.max_disp/wall);
    end
    fprintf('\n*** SOLVED *** %.1f s, area %.10f m^2 (rel err %.3e)\n', ...
        r.seconds,r.area,r.rel_area_err);
    fn = fullfile(outDir,sprintf('%s_smooth_%s.go3',stem,tag));
    surfsmooth3d.edgepreserve.write_go3(fn,Sout);
    r.file = string(fn);
    fprintf('wrote %s\n',fn);
    clear S Sout;
catch ME
    r.seconds = toc(tsolve); r.status = "failed";
    r.identifier = string(ME.identifier);
    fprintf(2,'\n*** FAILED *** after %.1f s: %s\n',r.seconds,ME.identifier);
    fprintf(2,'%s\n',strtrim(regexprep(ME.message,'\s+',' ')));
    % Guarded-projection diagnostic sits next to the scaffold.
    dg = fullfile(workDir,'tmp_newton_failure.txt');
    if isfile(dg)
        copyfile(dg,fullfile(outDir,'newton_failure.txt'));
        L = readlines(dg);
        for k = 1:min(10,numel(L))
            t = strtrim(L(k));
            if startsWith(t,'stage '),          r.fail_stage  = extractAfter(t,'stage ');  end
            if startsWith(t,'failure_reason '), r.fail_reason = extractAfter(t,'failure_reason '); end
        end
        fprintf(2,'stage: %s | reason: %s\n',r.fail_stage,r.fail_reason);
        fprintf(2,'newton_failure.txt saved (%d lines)\n',numel(L));
    end
    % The recovery report is preserved to a tempname and attached as a cause.
    c = ME.cause;
    for k = 1:numel(c)
        if strcmp(c{k}.identifier,'MULTISCALE_MESHER:NewtonRecoveryFile') && isfile(c{k}.message)
            copyfile(c{k}.message,fullfile(outDir,'newton_recovery.txt'));
            rr = fullfile(outDir,'newton_recovery.txt');
            fprintf(2,'recovery report -> %s\n',rr);
            txt = readlines(rr);
            r.deferred   = sum(startsWith(strtrim(txt),'TARGET'));
            r.recovered  = sum(strtrim(txt) == "OUTCOME recovered");
            r.unresolved = r.deferred - r.recovered;
            fprintf(2,'recovery: deferred %d, recovered %d, unresolved %d\n', ...
                r.deferred,r.recovered,r.unresolved);
            outs = txt(startsWith(strtrim(txt),'OUTCOME'));
            [u,~,ic] = unique(outs); cnt = accumarray(ic,1);
            for q = 1:numel(u), fprintf(2,'   %-40s %d\n',u(q),cnt(q)); end
        end
    end
end

writetable(struct2table(r),fullfile(outDir,'manifest.csv'));
fprintf('\nmanifest: %s\n',fullfile(outDir,'manifest.csv'));
fprintf('finished %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
% Keep the oversampled skeleton out of the way; it is large and reproducible.
try, delete(cadFile); catch, end
diary off;
end

function s = modeword(as)
if as == 0, s = '(uniform)'; else, s = '(max)'; end
end
function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end

function [nodes,tris] = read_gid(f)
L = readlines(f);
iC = find(strcmpi(strtrim(L),'Coordinates'));
eC = find(strcmpi(strtrim(L),'End Coordinates'));
iE = find(strcmpi(strtrim(L),'Elements'));
eE = find(strcmpi(strtrim(L),'End Elements'));
C = cell2mat(cellfun(@(s) sscanf(s,'%f').',cellstr(L(iC+1:eC-1)),'UniformOutput',false));
E = cell2mat(cellfun(@(s) sscanf(s,'%f').',cellstr(L(iE+1:eE-1)),'UniformOutput',false));
nodes = C(:,2:4).'; tris = E(:,2:4).';
end
