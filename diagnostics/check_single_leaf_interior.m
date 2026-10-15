function check_single_leaf_interior
% Check vanishing interior total fields for the closed PEC representation.
root = fileparts(fileparts(mfilename('fullpath')));
run('../fmm3dbie-hirax-dev/matlab/startup.m')
orders = [6 8 10 12];
targets = struct('r',[0 40 60;0 0 0;0 0 0]);
results = cell(size(orders));
for k = 1:numel(orders)
    d = load(fullfile(root,'data',sprintf('run_hirax_single_leaf_order%d.mat',orders(k))));
    c = d.case_settings;
    source = struct('r',c.x0,'edips',-c.p0);
    [Ei,Hi] = em3d.incoming_sources(c.zk,source,targets,'electric dipole');
    timer = tic;
    [Es,Hs] = em3d.pec.eval(d.S,d.densities,targets,1e-8,c.zk,c.alpha, ...
        struct('rep','nrccie'));
    a = struct('order',orders(k),'targets',targets.r,'Etotal',Ei+Es, ...
        'Htotal',Hi+Hs,'Einc',Ei,'Hinc',Hi);
    a.relative_E = vecnorm(Ei+Es)./vecnorm(Ei);
    a.relative_H = vecnorm(Hi+Hs)./vecnorm(Hi);
    a.elapsed = toc(timer);
    results{k} = a;
    disp(a)
    save(fullfile(root,'diagnostics','single_leaf_interior.mat'),'results');
end
fprintf('INTERIOR_CHECK_DONE\n');
end
