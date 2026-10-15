function fig = plot_toy_array_scattering(result,visible)
% PLOT_TOY_ARRAY_SCATTERING Summarize the toy symmetry/coupling rehearsal.
%   fig = plot_toy_array_scattering(result,'off')
%   fig = plot_toy_array_scattering('/path/to/toy_array_result.mat','on')
%
%  result is the struct returned by TEST_TOY_ARRAY_SCATTERING or its MAT
%  filename. Only numeric result data are loaded; no surfer, native library,
%  geometry generation, or solver is needed. Station-only checkpoints and
%  failed validation runs can also be plotted without hiding their status.
%  PNG and FIG files are saved beside result.output_file, with the suffix
%  _summary. The default figure visibility is 'off'.

    if nargin<2 || isempty(visible), visible = 'off'; end
    visible = validatestring(visible,{'on','off'});
    source_file = '';
    if ischar(result) || (isstring(result) && isscalar(result))
      source_file = char(result);
      requested = {'model','stats','errors','symmetry','array','output_file', ...
          'accuracy_limit','norder','nproxy','tol','passed','run_array'};
      contents = whos('-file',source_file);
      selected = intersect(requested,{contents.name},'stable');
      result = load(source_file,selected{:});
    end
    assert(isstruct(result) && isfield(result,'stats') && isfield(result,'errors'), ...
        'Expected a saved toy-array result or station checkpoint.');
    if isfield(result,'output_file') && ~isempty(result.output_file)
      output_file = result.output_file;
    else
      assert(~isempty(source_file),'The result struct needs an output_file field.');
      output_file = source_file;
    end
    [folder,stem] = fileparts(output_file);
    if isempty(folder), folder = pwd; end
    assert(isfolder(folder),'The result output directory does not exist.');
    s = result.stats;
    e = result.errors;
    has_array = isfield(result,'array') && isstruct(result.array) ...
        && ~isempty(result.array) && isfield(e,'array');
    limit = value_or(result,'accuracy_limit',5e-3);
    status = 'Station checkpoint';
    if isfield(result,'passed')
      if result.passed
        status = 'Validation passed';
      else
        status = 'Validation threshold not met';
      end
    end
    fig = figure('Visible',visible,'Color','w', ...
        'Name','Toy array scattering rehearsal','NumberTitle','off', ...
        'Position',[60,100,1500,600]);
    caption(fig,[0.02,0.93,0.96,0.055], ...
        'Symmetry and nested scattering on the toy geometry',18,'bold');
    subtitle = sprintf('Physical order %d | proxy order %d | ID tolerance %.0e | %s', ...
        value_or(result,'norder',NaN),value_or(result,'nproxy',NaN), ...
        value_or(result,'tol',NaN),status);
    caption(fig,[0.02,0.88,0.96,0.04],subtitle,11,'normal');
    caption(fig,[0.025,0.81,0.29,0.045],'Construction and reuse',13,'bold');
    caption(fig,[0.36,0.81,0.28,0.045],'One station: scattered E/H',13,'bold');
    caption(fig,[0.70,0.81,0.28,0.045],'Four stations: scattered E/H',13,'bold');

    ax = axes(fig,'Position',[0.035,0.14,0.27,0.64]);
    hold(ax,'on'); axis(ax,[0,1,0,1]); axis(ax,'off');
    nquad = 4;
    nquad_full = 16;
    if isfield(result,'symmetry')
      nquad = value_or(result.symmetry,'quadrature_blocks_computed',4);
      nquad_full = value_or(result.symmetry,'quadrature_blocks_without_reuse',16);
    end
    feed_size = s.feed_scattering_size(1);
    outer_size = s.outer_scattering_size(1);
    factor_size = s.feed_physical_dofs/4;
    box(ax,[0,0.71,1,0.27],[0.92,0.96,1.00], ...
        'Four rotated leaves: C4 symmetry', ...
        {sprintf('Quadrature blocks: %d computed / %d',nquad,nquad_full), ...
         sprintf('Four %d x %d Fourier-mode factors',factor_size,factor_size), ...
         sprintf('Feed scattering matrix: %d x %d',feed_size,feed_size)});
    text(ax,0.5,0.665,'\downarrow','HorizontalAlignment','center','FontSize',21);
    box(ax,[0,0.37,1,0.25],[0.93,0.97,0.94], ...
        'Feed interface + physical toy dish', ...
        {sprintf('Full station: %d physical DOFs',s.station_physical_dofs), ...
         sprintf('Hybrid system: %d DOFs',s.hybrid_dofs), ...
         sprintf('Station scattering matrix: %d x %d',outer_size,outer_size)});
    text(ax,0.5,0.325,'\downarrow','HorizontalAlignment','center','FontSize',21);
    if has_array
      as = result.array.stats;
      reuse_lines = {sprintf('Station interface builds: %d',as.station_scattering_constructions), ...
          sprintf('Translated reuses: %d',as.translated_reuses), ...
          sprintf('Array DOFs: %d compressed / %d physical', ...
              as.compressed_array_dofs,as.physical_array_dofs)};
      reuse_title = 'The same station matrix, four locations';
    else
      reuse_lines = {'One station interface constructed', ...
          'Ready for three translated copies', ...
          'Array solve not included in this result'};
      reuse_title = 'Whole-station interface ready for reuse';
    end
    box(ax,[0,0.02,1,0.26],[1.00,0.96,0.89],reuse_title,reuse_lines);

    ax = axes(fig,'Position',[0.365,0.29,0.275,0.48]);
    hold(ax,'on');
    names = {'feed_field','hybrid_field','outer_field','outer_green_reference'};
    labels = {'Feed interface alone','Coupled feed + body','Outer station interface', ...
        'Outer Green formula only'};
    colors = [0.18,0.43,0.73;0.00,0.52,0.38;0.57,0.31,0.69;0.78,0.46,0.12];
    markers = {'o','s','d','^'};
    handles = gobjects(0);
    used_labels = {};
    all_errors = [];
    nrhs = 0;
    for j = 1:numel(names)
      if ~isfield(e,names{j}), continue; end
      values = reshape(e.(names{j}),1,[]);
      nrhs = max(nrhs,numel(values));
      handles(end+1) = semilogy(ax,1:numel(values),max(values,1e-16), ...
          ['-',markers{j}],'Color',colors(j,:),'LineWidth',1.4,'MarkerSize',5); %#ok<AGROW>
      used_labels{end+1} = labels{j}; %#ok<AGROW>
      all_errors = [all_errors,values]; %#ok<AGROW>
    end
    handles(end+1) = yline(ax,limit,'--','Color',[0.4,0.4,0.4],'LineWidth',1);
    used_labels{end+1} = 'Validation threshold';
    style_error_axes(ax,nrhs,[all_errors,limit]);
    xlabel(ax,'Independent incident field'); ylabel(ax,'Relative field error');
    leg = legend(ax,handles,used_labels,'Box','off','FontSize',10);
    leg.Units = 'normalized'; leg.Position = [0.36,0.075,0.29,0.15];

    ax = axes(fig,'Position',[0.705,0.29,0.275,0.48]);
    if has_array
      hold(ax,'on');
      ae = e.array;
      coupled = reshape(ae.target_vs_physical_reference,1,[]);
      control = reshape(ae.no_array_coupling_vs_physical_reference,1,[]);
      h1 = semilogy(ax,1:numel(coupled),max(coupled,1e-16),'-o', ...
          'Color',colors(3,:),'LineWidth',1.5,'MarkerSize',6);
      h2 = semilogy(ax,1:numel(control),max(control,1e-16),'-s', ...
          'Color',colors(4,:),'LineWidth',1.5,'MarkerSize',6);
      array_limit = value_or(result.array,'accuracy_limit',limit);
      h3 = yline(ax,array_limit,'--','Color',[0.4,0.4,0.4],'LineWidth',1);
      style_error_axes(ax,numel(coupled),[coupled,control,array_limit]);
      xlabel(ax,'Independent incident field'); ylabel(ax,'Relative field error');
      leg = legend(ax,[h1,h2,h3],{'Coupled station interfaces', ...
          'Inter-station coupling omitted','Validation threshold'}, ...
          'Box','off','FontSize',10);
      leg.Units = 'normalized'; leg.Position = [0.70,0.105,0.29,0.12];
      caption(fig,[0.70,0.06,0.29,0.035],sprintf('Full-array NRCCIE residual: %.1e', ...
          max(ae.physical_reference_residual)),10,'normal');
    else
      axis(ax,'off');
      text(ax,0.5,0.55,{'Array coupling was not run.', ...
          'The station interface can be translated', ...
          'and reused at all four locations.'}, ...
          'Units','normalized','HorizontalAlignment','center', ...
          'FontSize',12,'Color',[0.35,0.35,0.35]);
    end
    caption(fig,[0.02,0.005,0.96,0.04], ...
        ['Reference: corrected NRCCIE on the same physical mesh. ', ...
         'Free space; ground reflection and mesh-convergence claims are not included.'], ...
        10,'normal');
    drawnow;
    png_file = fullfile(folder,[stem,'_summary.png']);
    fig_file = fullfile(folder,[stem,'_summary.fig']);
    exportgraphics(fig,png_file,'Resolution',150);
    savefig(fig,fig_file);
    fprintf('Toy scattering summary saved: %s\n',png_file);
end

function box(ax,position,color,title_text,lines)
    rectangle(ax,'Position',position,'Curvature',0.06, ...
        'FaceColor',color,'EdgeColor',[0.72,0.76,0.79],'LineWidth',0.8);
    top = position(2)+position(4);
    text(ax,0.045,top-0.035,title_text,'FontSize',10.8, ...
        'FontWeight','bold','VerticalAlignment','top','Interpreter','none');
    for j = 1:numel(lines)
      text(ax,0.045,top-0.095-(j-1)*0.054,lines{j},'FontSize',10.3, ...
          'VerticalAlignment','top','Interpreter','none');
    end
end

function style_error_axes(ax,nrhs,values)
    set(ax,'YScale','log','FontSize',10,'Box','on','XGrid','on','YGrid','on');
    xlim(ax,[0.65,max(1,nrhs)+0.35]);
    xticks(ax,1:max(1,nrhs));
    values = values(isfinite(values) & values>0);
    if isempty(values), values = 1e-12; end
    low = max(1e-16,10^floor(log10(min(values))));
    high = 10^ceil(log10(max(values)));
    if high<=low, high = 10*low; end
    ylim(ax,[low,high]);
end

function caption(fig,position,label,size,weight)
    annotation(fig,'textbox',position,'String',label, ...
        'EdgeColor','none','HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontSize',size,'FontWeight',weight, ...
        'Interpreter','none');
end

function value = value_or(data,name,fallback)
    if isfield(data,name) && ~isempty(data.(name))
      value = data.(name);
    else
      value = fallback;
    end
end
