function smooth_dish_rimfree()
%SMOOTH_DISH_RIMFREE Historical experimental smoothing of the rimfree dish.
%
% Target : <DISH_MESH_ROOT>/dish_fd021_m_rimfree
% Source : dish_fd021_m_rimfree_order4.go3     (fixed CAD quadrature, p=4)
%          dish_fd021_m_rimfree_scaffold.gidmsh (launch scaffold, 6010 tris)
% Settings: orders 4 and 8, nrefine 0, rlam 10, adapt_sigma 1, two-stage.

stem    = 'dish_fd021_m_rimfree';
[dishDir,repoRoot] = dish_smoothing_paths(stem);
go3In   = fullfile(dishDir,[stem '_order4.go3']);
gidIn   = fullfile(dishDir,[stem '_scaffold.gidmsh']);
assert(isfile(go3In),'missing %s',go3In);
assert(isfile(gidIn),'missing %s',gidIn);
assert(isfile(fullfile(repoRoot,'matlab','setup_surfsmooth3d.m')), ...
    'Missing SurfSmooth3D; set SURFSMOOTH3D_ROOT to its repository root.');
addpath(fullfile(repoRoot,'matlab'));
setup_surfsmooth3d();

% Exact analytic area for the same CAD, from the mesher's own report.
areaExact = 73.74006397;   % m^2 (cad.area_total_exact)
wall      = 0.01;          % m   (cad.wall_thickness)

% Stage inputs in a scratch work dir: multiscale_mesher writes its tmp_*.go3
% and tmp_newton_failure.txt next to the input mesh.
outDir = fullfile(dishDir,'smooth');
if ~isfolder(outDir), mkdir(outDir); end
workRoot = fullfile(outDir,'work');
if ~isfolder(workRoot), mkdir(workRoot); end
workDir = tempname(workRoot);
mkdir(workDir);
scaffold = fullfile(workDir,'scaffold.gidmsh');
copyfile(gidIn,scaffold);

diary(fullfile(outDir,'smooth_log.txt')); diary on;
cleanupDiary = onCleanup(@() diary('off'));

fprintf('\n=== %s : fixed-CAD two-stage smoothing ===\n',stem);

%% Fixed CAD source -> 7-column point skeleton
S_cad = surfsmooth3d.surfer.load_from_file(go3In);
sourceOrder = double(S_cad.norders(1));
fprintf('CAD source: p=%d, npatches=%d, npts=%d\n', ...
    sourceOrder,S_cad.npatches,S_cad.npts);

rep = surfsmooth3d.edgepreserve.validity_report(S_cad,'fixed CAD source');
surfsmooth3d.edgepreserve.print_validity_report(rep);
areaCad = sum(S_cad.wts);
fprintf('  CAD quadrature area  %.10e m^2  (rel err vs exact %.3e)\n', ...
    areaCad,abs(areaCad-areaExact)/areaExact);

cadFile = fullfile(workDir,'cad_skeleton.txt');
surfsmooth3d.edgepreserve.write_point_source(S_cad,cadFile);
d = dir(cadFile);
fprintf('  wrote CAD skeleton: %d points, %.1f MB\n\n', ...
    S_cad.npts,d.bytes/1e6);

%% Sigma budget against the shell thickness
rlam = 10; adapt_sigma = 1;
[~,tris,nodes] = local_read_gidmsh(scaffold);
P = reshape(nodes(:,tris),3,3,[]);
e = [vecnorm(squeeze(P(:,2,:)-P(:,1,:))); ...
     vecnorm(squeeze(P(:,3,:)-P(:,2,:))); ...
     vecnorm(squeeze(P(:,1,:)-P(:,3,:)))];
dmax = max(e,[],1);
sig  = dmax/rlam;   % Mod_Fast_Sigma.f90: sgmas = dmax/rlam
fprintf('scaffold: %d triangles, longest side %.4g - %.4g m\n', ...
    size(tris,2),min(dmax),max(dmax));
fprintf('sigma = D/rlam with rlam=%g: median %.4g, max %.4g m\n',rlam, ...
    median(sig),max(sig));
fprintf('shell wall thickness %.4g m -> %.1f%% of triangles have sigma > wall/2\n\n', ...
    wall,100*mean(sig > wall/2));

%% Solve
orders  = [4 8];
nrefine = 0;
results = struct('order',{},'status',{},'message',{},'npatches',{}, ...
    'npts',{},'area',{},'rel_area_err',{},'validity',{}, ...
    'normal_cross_dot_min',{},'jac_ratio_min',{},'bad_nodes',{}, ...
    'max_disp',{},'seconds',{},'file',{});

for order = orders
    fprintf('--- order %d, nrefine %d, rlam %g, adapt_sigma %d, two-stage ---\n', ...
        order,nrefine,rlam,adapt_sigma);
    opts = struct('filetype',3,'fcad',cadFile,'nquad',sourceOrder, ...
        'nrefine',nrefine,'rlam',rlam,'adapt_sigma',adapt_sigma, ...
        'two_stage_smoother',true);

    r = struct('order',order,'status',"",'message',"",'npatches',nan, ...
        'npts',nan,'area',nan,'rel_area_err',nan,'validity',"", ...
        'normal_cross_dot_min',nan,'jac_ratio_min',nan,'bad_nodes',nan, ...
        'max_disp',nan,'seconds',nan,'file',"");
    t0 = tic;
    try
        S = surfsmooth3d.multiscale_mesher(scaffold,order,opts);
        r.seconds = toc(t0);
        Sout = S{1};
        r.status = "solved";
        r.npatches = Sout.npatches; r.npts = Sout.npts;
        r.area = sum(Sout.wts);
        r.rel_area_err = abs(r.area-areaExact)/areaExact;
        fprintf('solved in %.1f s: npatches=%d npts=%d\n', ...
            r.seconds,Sout.npatches,Sout.npts);

        % Compare against the CAD source on its own layout when orders match.
        ref = [];
        if order == sourceOrder && Sout.npts == S_cad.npts && ...
                Sout.npatches == S_cad.npatches
            ref = S_cad;
        end
        rp = surfsmooth3d.edgepreserve.validity_report(Sout, ...
            sprintf('smoothed p%d r%d',order,nrefine),ref);
        surfsmooth3d.edgepreserve.print_validity_report(rp);
        if rp.is_valid && ~rp.has_warnings
            r.validity = "PASS";
        elseif rp.is_valid
            r.validity = "WARN";
        else
            r.validity = "FAIL";
        end
        r.normal_cross_dot_min = rp.stats.normal_cross_dot_min;
        r.bad_nodes = rp.stats.bad_node_count;
        if rp.has_reference
            r.jac_ratio_min = rp.stats.jac_ratio_min;
            r.max_disp = rp.stats.displacement_max;
            fprintf('  max displacement from CAD  %.6e m  (= %.2f x wall)\n', ...
                r.max_disp,r.max_disp/wall);
        end
        fprintf('  smoothed area        %.10e m^2  (rel err vs exact %.3e)\n', ...
            r.area,r.rel_area_err);

        outFile = fullfile(outDir,sprintf('%s_smooth_p%02d_r%02d_rlam%g_as%d.go3', ...
            stem,order,nrefine,rlam,adapt_sigma));
        surfsmooth3d.edgepreserve.write_go3(outFile,Sout);
        r.file = string(outFile);
        fprintf('  wrote %s\n\n',outFile);
        clear S Sout;
    catch ME
        r.seconds = toc(t0);
        r.status = "failed";
        r.message = string(ME.message);
        fprintf(2,'FAILED after %.1f s\n  %s: %s\n',r.seconds,ME.identifier,ME.message);
        % Preserve the Newton diagnostic the mesher leaves next to the input.
        diagFile = fullfile(workDir,'tmp_newton_failure.txt');
        if isfile(diagFile)
            keep = fullfile(outDir,sprintf('newton_failure_p%02d_r%02d.txt',order,nrefine));
            copyfile(diagFile,keep);
            fprintf(2,'  Newton diagnostic preserved: %s\n',keep);
            fprintf(2,'--- diagnostic head ---\n');
            txt = readlines(keep);
            fprintf(2,'%s\n',txt(1:min(20,numel(txt))));
        end
        fprintf(2,'\n');
    end
    results(end+1) = r; %#ok<AGROW>
end

%% Summary
fprintf('=== summary ===\n');
T = struct2table(results);
disp(T(:,{'order','status','validity','npatches','area','rel_area_err','max_disp','seconds'}));
writetable(T,fullfile(outDir,'smooth_manifest.csv'));
fprintf('manifest: %s\n',fullfile(outDir,'smooth_manifest.csv'));
diary off;
end

function [npts,tris,nodes] = local_read_gidmsh(f)
L = readlines(f);
iC = find(strcmpi(strtrim(L),'Coordinates'));
eC = find(strcmpi(strtrim(L),'End Coordinates'));
iE = find(strcmpi(strtrim(L),'Elements'));
eE = find(strcmpi(strtrim(L),'End Elements'));
C = cell2mat(cellfun(@(s) sscanf(s,'%f').',cellstr(L(iC+1:eC-1)),'UniformOutput',false));
E = cell2mat(cellfun(@(s) sscanf(s,'%f').',cellstr(L(iE+1:eE-1)),'UniformOutput',false));
nodes = C(:,2:4).';
tris  = E(:,2:4).';
npts  = size(nodes,2);
end
