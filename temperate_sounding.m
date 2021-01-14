% TEMPERATE_SOUNDING Script to perform analysis and figure generation for MacGregor et al. (in review, The Cryosphere Discussions).
%
% Joe MacGregor (NASA/GSFC)
% Last updated: 14 January 2020

clear

% switches
plotting                    = false; % do plots

% directories
dir_GTD                     = '/Users/jamacgre/Documents/data/glathida-3.1.0/data/'; % GlaThiDa v3.1.0 directory
dir_F19                     = '/Users/jamacgre/Documents/data/consensus_thickness/'; % Farinotti et al. (2019, Nature Geoscience) GeoTIFF directory, with sub-directories by RGI Level-1 region
dir_RGI                     = '/Users/jamacgre/Documents/data/00_rgi60/'; % Randolph Glacier Inventory version 6 directory, with sub-directories by RGI Level-1 region

%% GlaThiDA (GTD) loading and analysis

disp('Load and analyze GlaTHiDa v3.1.0...')

% load GlaThiDa meta-database
GTD_T                       = readtable([dir_GTD 'T.csv']); % read in GlaThiDa .csv
num_GTD_T                   = size(GTD_T, 1);
GTD_T.SURVEY_YEAR           = NaN(num_GTD_T, 1);
for ii = find(~isnan(GTD_T.SURVEY_DATE))'
    tmp_str                 = num2str(GTD_T.SURVEY_DATE(ii));
    GTD_T.SURVEY_YEAR(ii)   = str2double(tmp_str(1:4)); % extract survey year
end

% read in CSV file supplementing GlaThiDa with survey frequency in MHz, updated survey method (adjusted to include GPRh, helicopter-based sounding)
GTD_T_supp                  = readtable([dir_GTD 'T_supplement.csv']);

% add in missing mean/max thickness values from GlaThiDa raw database (TTT), mostly relevant to WISE data
GTD_T_adj                   = struct('MAXIMUM_THICKNESS', GTD_T.MAXIMUM_THICKNESS);
GTD_TTT                     = readtable([dir_GTD 'TTT']);
for ii = find(isnan(GTD_T.MAXIMUM_THICKNESS))'
    idx_GTD_curr            = find(GTD_TTT.GlaThiDa_ID == GTD_T.GlaThiDa_ID(ii)); % current GlaThiDa_ID of interest
    if isempty(idx_GTD_curr) % skip if ID not found
        continue
    end
    if isnan(GTD_T.MAXIMUM_THICKNESS(ii)) % add a maximum value that's missing
        GTD_T_adj.MAXIMUM_THICKNESS(ii) ...
                            = max(GTD_TTT.THICKNESS(idx_GTD_curr));
    end
end
clear GTD_TTT % save some memory
GTD_T_adj.MAXIMUM_THICKNESS(GTD_T_adj.MAXIMUM_THICKNESS == 0) ...
                            = NaN; % NaN our zeros

% GTD indices for likely temperate glacier surveys
idx_GTD_temperate           = setdiff(1:num_GTD_T,       find(abs(GTD_T.LAT) >= 67)); % exclude polar regions based on latitude
idx_GTD_temperate           = setdiff(idx_GTD_temperate, find(ismember(GTD_T.POLITICAL_UNIT, {'AQ' 'GL' 'SJ'}))); % exclude Antarctica, Greenland, Svalbard, non-western Canada based on political unit/country
idx_GTD_temperate           = setdiff(idx_GTD_temperate, find((GTD_T.LON >= -112) & ismember(GTD_T.POLITICAL_UNIT, 'CA'))); % exclude Arctic Canada based on longitude
idx_GTD_temperate           = setdiff(idx_GTD_temperate, find(strcmp(GTD_T.GLACIER_NAME, 'HAZARD'))); % exclude Hazard Glacier, based on Narod and Clarke (1980) using UHF over identified cold/polar ice
idx_GTD_temperate           = setdiff(idx_GTD_temperate, find(strcmp(GTD_T.GLACIER_NAME, 'GORNER'))); % exclude Gornergletscher, which has been previously identified as polythermal
idx_GTD_temperate           = setdiff(idx_GTD_temperate, find(strcmp(GTD_T.GLACIER_NAME, 'KHUK NURU UUL'))); % exclude Khukh Nuru Uul, which has been previously identified as cold (Herren et al., 2013, QSR)

% temperate region glacier surveys done by radar for which a maximum thickness is available
idx_GTD_temperate_radar_good= intersect(idx_GTD_temperate, find(~isnan(GTD_T_supp.SURVEY_FREQUENCY) & contains(GTD_T_supp.SURVEY_METHOD, {'GPRt' 'GPRh' 'GPRa'}) & ~isnan(GTD_T_adj.MAXIMUM_THICKNESS)));
num_GTD_temperate_radar_good= length(idx_GTD_temperate_radar_good);

% number by survey type
num_GTD_temperate_radar_ground ...
                            = length(find(contains(GTD_T_supp.SURVEY_METHOD(idx_GTD_temperate_radar_good), 'GPRt')));
num_GTD_temperate_radar_helo ...
                            = length(find(contains(GTD_T_supp.SURVEY_METHOD(idx_GTD_temperate_radar_good), 'GPRh')));
num_GTD_temperate_radar_fixed_wing ...
                            = length(find(contains(GTD_T_supp.SURVEY_METHOD(idx_GTD_temperate_radar_good), 'GPRa')));

% indices of non-temperate surveys
idx_GTD_other_ground        = setdiff(find(~isnan(GTD_T_supp.SURVEY_FREQUENCY) & contains(GTD_T_supp.SURVEY_METHOD, 'GPRt') & ~isnan(GTD_T_adj.MAXIMUM_THICKNESS)), idx_GTD_temperate)';
idx_GTD_other_helo          = setdiff(find(~isnan(GTD_T_supp.SURVEY_FREQUENCY) & contains(GTD_T_supp.SURVEY_METHOD, 'GPRh') & ~isnan(GTD_T_adj.MAXIMUM_THICKNESS)), idx_GTD_temperate)';
idx_GTD_other_fixed_wing    = setdiff(find(~isnan(GTD_T_supp.SURVEY_FREQUENCY) & contains(GTD_T_supp.SURVEY_METHOD, 'GPRa') & ~isnan(GTD_T_adj.MAXIMUM_THICKNESS)), idx_GTD_temperate)';

%% ALL RGI PREDICTED THICKNESS ANALYSIS

disp('Load RGI shapefiles for assumed temperate regions...')

% RGI region names
RGI_str                     = {'01 Alaska';
                               '02 Western Canada / U.S.A.';
                               '03 Arctic Canada (North)';
                               '04 Arctic Canada (South)';
                               '05 Greenland periphery';
                               '06 Iceland';
                               '07 Svalbard';
                               '08 Scandinavia';
                               '09 Russian Arctic';
                               '10 North Asia';
                               '11 Central Europe';
                               '12 Caucasus / Middle East';
                               '13 Central Asia';
                               '14 South Asia (West)';
                               '15 South Asia (East)';
                               '16 Low latitudes';
                               '17 Southern Andes';
                               '18 New Zealand';
                               '19 Antarctic periphery / Subantarctic'};
num_RGI_reg                 = length(RGI_str);

idx_RGI_temperate           = setdiff(1:num_RGI_reg, [3:5 7 9 19]); % exclude regions assumed to contain mostly polar or polythermal glaciers

% RGI region directory names
RGI_set                     = dir(dir_RGI);
RGI_set                     = RGI_set([RGI_set(:).isdir]);
RGI_set                     = {RGI_set.name};
RGI_set                     = RGI_set(~contains(RGI_set, {'.' '..'}));

% load RGI shapefiles for temperate regions only
[RGI, RGI_name_cell]        = deal(cell(1, num_RGI_reg));
for ii = idx_RGI_temperate
    disp([RGI_str{ii} '...'])
    RGI_str_curr            = num2str(ii);
    if (ii < 10) % add leading zero if less than 10
        RGI_str_curr        = ['0' RGI_str_curr]; %#ok<AGROW>
    end
    RGI{ii}                 = shaperead([dir_RGI RGI_set{contains(RGI_set, RGI_str_curr)} '/' RGI_set{contains(RGI_set, RGI_str_curr)} '.shp'], 'Attributes', {'Area' 'CenLat', 'CenLon', 'Name' 'RGIId', 'Zmed'}); % load from sub-directory for each region
    RGI{ii}                 = rmfield(RGI{ii}, {'BoundingBox' 'Geometry' 'X' 'Y'}); % don't need glacier outlines or most geographic information, so remove to save memory
    for jj = 1:length(RGI{ii})
        RGI{ii}(jj).RGIIdNum= str2double(RGI{ii}(jj).RGIId((end - 4):end)); % generate integer value RGI ID number
    end
    RGI_name_cell{ii}       = {RGI{ii}(:).Name}; % cells that are just glacier names (easier for name matching later on)
end

%% F19 MODELED THICKNESS ANALYSIS

disp('Load F19 modeled thicknesses and calculate maximum value for each glacier...')

% loop through all F19 modeled glaciers and extract maximum modeled thickness
[F19_RGIId, thick_F19_max]  = deal(cell(1, num_RGI_reg));
for ii = idx_RGI_temperate
    disp(ii)
    RGI_str_curr            = num2str(ii);
    if (ii < 10) % add leading zero if less than 10
        RGI_str_curr        = ['0' RGI_str_curr]; %#ok<AGROW>
    end
    file_F19_curr           = dir([dir_F19 'RGI60-' RGI_str_curr '/*.tif']); % all GeoTIFFs in current directory
    file_F19_curr           = {file_F19_curr.name};
    num_F19_curr            = length(file_F19_curr);
    [F19_RGIId{ii}, thick_F19_max{ii}] ...
                            = deal(NaN(num_F19_curr, 1));
    for jj = 1:num_F19_curr % loop through all glaciers for RGI region modeled by F19 (for some regions not all glaciers modeled)
        if ~mod(jj, 500)
            disp([num2str(1e2 * (jj / num_F19_curr)) '%']) % report percentage complete every 500 glaciers
        end
        F19_RGIId{ii}(jj)   = str2double(file_F19_curr{jj}(10:14)); % RGI ID for current glacier
        thick_F19_max{ii}(jj) ...
                            = max(readgeoraster([dir_F19 'RGI60-' RGI_str_curr '/' file_F19_curr{jj}]), [], 'all', 'omitnan'); % maximum modeled consensus ice thickness for current glacier
    end
end

% find large glaciers in F19 set
area_RGI_large              = 5; % threshold for "large" glaciers, km^2
idx_F19_large               = cell(1, num_RGI_reg);
for ii = idx_RGI_temperate
    [~, ~, idx_F19_large{ii}] ...
                            = intersect([RGI{ii}([RGI{ii}(:).Area] >= area_RGI_large).RGIIdNum], F19_RGIId{ii}); % indices in F19 set whose RGI IDs match with large glaciers
end

%% MATCH GLATHIDA GLACIERS WITH RGI/F19 GLACIER IDS

disp('Match each survey of interest in GlaThiDa to an RGI ID to then match to F19 modeled thickness...')

dist_match_max              = 5e3; % maximum match distance for glaciers, km

RGI_reg                     = shaperead([dir_RGI '00_rgi60_regions/00_rgi60_O1Regions.shp']); % RGI Level-1 region outlines

[dist_GTD_RGI, GTD_RGIId, GTD_match_type, idx_GTD_RGI, reg_GTD_RGI] ...
                            = deal(NaN(1, num_GTD_temperate_radar_good)); 

% loop through every temperate glacier suveyed by radar and find matching glacier in RGI
for ii = 1:num_GTD_temperate_radar_good
    jj                      = idx_GTD_temperate_radar_good(ii); % index within GTD list (NOT GTD ID)
    for kk = 1:num_RGI_reg
        if inpolygon(GTD_T.LON(jj), GTD_T.LAT(jj), RGI_reg(kk).X, RGI_reg(kk).Y)
            reg_GTD_RGI(ii) = RGI_reg(kk).RGI_CODE; % identify RGI region glacier belongs to 
            break
        end
    end
    if isnan(reg_GTD_RGI(ii)) % couldn't find RGI region, stop search for RGIv6 ID
        continue
    end
    if strncmp(GTD_T.GLACIER_ID{jj}, 'RGI60', 5) % easy: GTD Glacier ID is an RGIv6 ID
        GTD_RGIId(ii)       = str2double(GTD_T.GLACIER_ID{jj}((end - 4):end));
        idx_GTD_RGI(ii)     = find(GTD_RGIId(ii) == [RGI{reg_GTD_RGI(ii)}(:).RGIIdNum]);
        GTD_match_type(ii)  = 1;
    elseif isscalar(find(contains(RGI_name_cell{reg_GTD_RGI(ii)}, GTD_T.GLACIER_NAME{jj}, 'ignorecase', true))) % not too hard: only one glacier name match within region
        idx_GTD_RGI(ii)     = find(contains(RGI_name_cell{reg_GTD_RGI(ii)}, GTD_T.GLACIER_NAME{jj}, 'ignorecase', true));
        GTD_RGIId(ii)       = RGI{reg_GTD_RGI(ii)}(idx_GTD_RGI(ii)).RGIIdNum;
        GTD_match_type(ii)  = 2;
    elseif ~isempty(find(contains(RGI_name_cell{reg_GTD_RGI(ii)}, GTD_T.GLACIER_NAME{jj}, 'ignorecase', true), 1)) % harder: more than one glacier name match within region, so pick closest within threshold
        idx_RGI_test        = find(contains(RGI_name_cell{reg_GTD_RGI(ii)}, GTD_T.GLACIER_NAME{jj}, 'ignorecase', true));
        [dist_GTD_RGI(ii), idx_GTD_RGI_best] ...
                            = min(distance([RGI{reg_GTD_RGI(ii)}(idx_RGI_test).CenLat], [RGI{reg_GTD_RGI(ii)}(idx_RGI_test).CenLon], GTD_T.LAT(jj .* ones(1, length(idx_RGI_test)))', GTD_T.LON(jj .* ones(1, length(idx_RGI_test)))', wgs84Ellipsoid));
        if (dist_GTD_RGI(ii) > dist_match_max)
            [dist_GTD_RGI(ii), idx_GTD_RGI_best] ...
                            = deal(NaN);
            continue
        end
        idx_GTD_RGI(ii)     = idx_RGI_test(idx_GTD_RGI_best);
        GTD_RGIId(ii)       = RGI{reg_GTD_RGI(ii)}(idx_RGI_test(idx_GTD_RGI_best)).RGIIdNum;
        GTD_match_type(ii)  = 3;
    else % give up on matching via id or name and pick closest glacier
        [dist_GTD_RGI(ii), idx_GTD_RGI(ii)] ...
                            = min(distance([RGI{reg_GTD_RGI(ii)}(:).CenLat], [RGI{reg_GTD_RGI(ii)}(:).CenLon], GTD_T.LAT(jj .* ones(1, length(RGI{reg_GTD_RGI(ii)})))', GTD_T.LON(jj .* ones(1, length(RGI{reg_GTD_RGI(ii)})))', wgs84Ellipsoid));
        if (dist_GTD_RGI(ii) > dist_match_max)
            [dist_GTD_RGI(ii), idx_GTD_RGI(ii)] ...
                            = deal(NaN); % distance too great so NaN out
            continue
        else
            GTD_RGIId(ii)   = RGI{reg_GTD_RGI(ii)}(idx_GTD_RGI(ii)).RGIIdNum; % distance close enough so assign RGIId
        end
        GTD_match_type(ii)  = 4;
    end
end

% concatenate measured/modeled thicknesses
GTD_F19_thick_max_cat       = NaN(num_GTD_temperate_radar_good, 2);
GTD_RGI_area                = NaN(1, num_GTD_temperate_radar_good);
for ii = find(~isnan(GTD_RGIId))
    GTD_F19_thick_max_cat(ii, :) ...
                            = [GTD_T_adj.MAXIMUM_THICKNESS(idx_GTD_temperate_radar_good(ii)) thick_F19_max{reg_GTD_RGI(ii)}(F19_RGIId{reg_GTD_RGI(ii)} == GTD_RGIId(ii))]; % GTD then F19 maximum thickness for surveyed/modeled glaciers
    GTD_RGI_area(ii)        = RGI{reg_GTD_RGI(ii)}(GTD_RGIId(ii)).Area; % RGI area for identified glaciers
end

% mean and standard deviation for well-matched measured/modeled thicknesses
GTD_RGI_diff_thick_mean     = mean(diff(GTD_F19_thick_max_cat((GTD_match_type <= 3), :), [], 2), 'omitnan');
GTD_RGI_diff_thick_std      = std(diff(GTD_F19_thick_max_cat((GTD_match_type <= 3), :), [], 2), 'omitnan');

%%
if plotting

%% GLATHIDA FREQUENCY VS. MAXIMUM THICKNESS MEASURED
    
    year_range              = 1970:10:2020;
    colors                  = [1 1 1; flipud(hot(length(year_range)))];
    colors                  = [colors(2, :); 1 1 0.8; colors(4:end, :)];
    idx_GTD_year_color      = 1 + discretize(GTD_T.SURVEY_YEAR, year_range);
    idx_GTD_year_color(isnan(idx_GTD_year_color)) ...
                            = 1;
    [~, idx_ord]            = sort(GTD_T.SURVEY_YEAR(idx_GTD_temperate_radar_good), 'ascend');
    figure('position', [100 100 960 540], 'color', 'w')
    colormap(colors)
    hold on
    axis([1 1e3 0 1.5e3])
    ax                      = gca;
    fill(ax.XLim([1 2 2]), ax.YLim([2 1 2]), [0.9 0.9 1], 'linestyle', 'none', 'facealpha', 0.5)
    fill(ax.XLim([1 2 1]), ax.YLim([2 1 1]), [1 0.9 0.9], 'linestyle', 'none', 'facealpha', 0.5)
    for ii = idx_GTD_other_ground
        line(GTD_T_supp.SURVEY_FREQUENCY(ii), GTD_T_adj.MAXIMUM_THICKNESS(ii), 'color', 'k', 'linewidth', 0.5, 'marker', 'v', 'markersize', 8, 'markerfacecolor', [0.9 0.9 1])
    end
    for ii = idx_GTD_other_helo
        line(GTD_T_supp.SURVEY_FREQUENCY(ii), GTD_T_adj.MAXIMUM_THICKNESS(ii), 'color', 'k', 'linewidth', 0.5, 'marker', 'o', 'markersize', 8, 'markerfacecolor', [0.9 0.9 1])
    end
    for ii = idx_GTD_other_fixed_wing
        line(GTD_T_supp.SURVEY_FREQUENCY(ii), GTD_T_adj.MAXIMUM_THICKNESS(ii), 'color', 'k', 'linewidth', 0.5, 'marker', '^', 'markersize', 8, 'markerfacecolor', [0.9 0.9 1])
    end
    line([1 1e3], [1.5e3 0], 'color', 'k', 'linestyle', '--', 'linewidth', 2)
    for ii = idx_GTD_temperate_radar_good(idx_ord)'
        switch GTD_T_supp.SURVEY_METHOD{ii}
            case 'GPRt' % "towed", i.e., ground-based
                line(GTD_T_supp.SURVEY_FREQUENCY(ii), GTD_T_adj.MAXIMUM_THICKNESS(ii), 'color', 'k', 'linewidth', 0.5, 'marker', 'v', 'markersize', 12, 'markerfacecolor', colors(idx_GTD_year_color(ii), :), ...
                                                                                            'tag', [num2str(ii) ': GlaThiDa ID ' num2str(GTD_T.GlaThiDa_ID(ii)) ' ' GTD_T.GLACIER_NAME{ii}])
            case 'GPRh' % helicopter-based, new distinction for this study
                line(GTD_T_supp.SURVEY_FREQUENCY(ii), GTD_T_adj.MAXIMUM_THICKNESS(ii), 'color', 'k', 'linewidth', 0.5, 'marker', 'o', 'markersize', 12, 'markerfacecolor', colors(idx_GTD_year_color(ii), :), ...
                                                                                            'tag', [num2str(ii) ': GlaThiDa ID ' num2str(GTD_T.GlaThiDa_ID(ii)) ' ' GTD_T.GLACIER_NAME{ii}])
            case 'GPRa' % "airborne", specified by this study to mean fixed-wing
                line(GTD_T_supp.SURVEY_FREQUENCY(ii), GTD_T_adj.MAXIMUM_THICKNESS(ii), 'color', 'k', 'linewidth', 0.5, 'marker', '^', 'markersize', 12, 'markerfacecolor', colors(idx_GTD_year_color(ii), :), ...
                                                                                            'tag', [num2str(ii) ': GlaThiDa ID ' num2str(GTD_T.GlaThiDa_ID(ii)) ' ' GTD_T.GLACIER_NAME{ii}])
        end
    end
    line(3.5, 440, 'color', 'k', 'linewidth', 0.5, 'marker', '^', 'markersize', 12, 'markerfacecolor', colors(end, :), 'tag', 'Pritchard et al. (2020, Annals of Glaciology)') 
    line(10, 318, 'color', 'k', 'linewidth', 0.5, 'marker', '^', 'markersize', 12, 'markerfacecolor', colors(end, :), 'tag', 'Pelto et al. (2020, Journal of Glaciology)')
    line(1.875, 1460, 'color', 'k', 'linewidth', 0.5, 'marker', 'v', 'markersize', 12, 'markerfacecolor', colors(end, :), 'tag', 'M. Truffer and J.W. Holt (2020, pers. comm.) for UAF HF radar sounder on Bagley Icefield')
    pg                      = plot(NaN, NaN, 'kv', 'linewidth', 0.5, 'markersize', 12);
    ph                      = plot(NaN, NaN, 'ko', 'linewidth', 0.5, 'markersize', 12);
    pa                      = plot(NaN, NaN, 'k^', 'linewidth', 0.5, 'markersize', 12);
    pht                     = plot(NaN, NaN, 'ko', 'linewidth', 0.5, 'markersize', 12, 'markerfacecolor', colors((end - 1), :));
    phc                     = plot(NaN, NaN, 'ko', 'linewidth', 0.5, 'markersize', 8, 'markerfacecolor', [0.9 0.9 1]);
    caxis([1 (length(year_range) + 1)])
    set(gca, 'fontsize', 20, 'fontweight', 'bold', 'xscale', 'log', 'xticklabel', {'1' '10' '100' '1000'}, 'linewidth', 1)
    xlabel('Center frequency (MHz)')
    ylabel({'Maximum ice thickness measured (m)'})
    lg1                      = legend([pa ph pg pht phc], 'fixed-wing', 'helicopter', 'ground', 'temperate', 'cold');
    set(lg1, 'location', 'northeast', 'fontsize', 20)
    colorbar('yticklabel', {'unknown' '1970' '1980' '1990' '2000' '2010' '2018'}, 'ticklength', 0.035, 'fontsize', 20, 'color', 'k')
    text(3.5112e3, 0.9409e3, 'Survey year', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'rotation', 270)
    text(2.2, 1445, {'empirical envelope of'; 'maximum measurable thickness'}, 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'rotation', -33)    
    grid on
    box on
    
%%  MAXIMUM MODELED THICKNESS DISTRIBUTIONS WITH FREQUENCY FOR TEMPERATE RGIs

    colors                  = [0   0   0;
                               0.5 0.5 0.5;
                               1   0   0;
                               0.8 0   0;
                               0   1   0;
                               0   0   1;
                               0   0   0.5;
                               1   0   1;
                               0.8 0   0.8;
                               0.5 0   0.5;
                               0   1   1;
                               0   0.8 0.8;
                               0.8 0.8 0];
    thick_int               = 25; % thickness bin interval
    thick_max               = 1500; % maximum thickness to consider
    thick_bin               = 0:thick_int:thick_max; % thickness bin vector
    letters                 = 'a':'b';
    ylims                   = [4e2 2e3];
    ylabels                 = {'Number of glaciers ≥ 5 km^2' 'Cumulative number of glaciers ≥ 5 km^2'};
    typ                     = [650 4000];
    figure('position', [100 100 1920 540], 'color', 'w', 'renderer', 'painters')
    for ii = 1:2
        subplot('position', [(0.08 + (0.50 * (ii - 1))) 0.12 0.39 0.75])
        hold on
        pm                      = length(idx_RGI_temperate);
        for jj = fliplr(1:length(idx_RGI_temperate))
            switch ii
                case 1
                    histogram(thick_F19_max{idx_RGI_temperate(jj)}(idx_F19_large{idx_RGI_temperate(jj)}), thick_bin, 'normalization', 'count', 'displaystyle', 'stairs', 'edgecolor', colors(jj, :), 'linewidth', 3)
                case 2
                    h = histogram(thick_F19_max{idx_RGI_temperate(jj)}(idx_F19_large{idx_RGI_temperate(jj)}), thick_bin, 'normalization', 'cumcount', 'displaystyle', 'stairs', 'edgecolor', colors(jj, :), 'linewidth', 3);
            end
            pm(jj)              = line(NaN, NaN, 'color', colors(jj, :), 'linewidth', 3);
        end
        axis([0 thick_max 0 ylims(ii)])
        set(gca, 'fontsize', 20, 'fontweight', 'bold', 'xtick', 0:100:1500, 'yscale', 'log');
        switch ii
            case 1
                set(gca, 'ytick', [1 10 100 400], 'yticklabel', [1 10 100 400])
            case 2
                set(gca, 'ytick', [1 10 100 1000 2000], 'yticklabel', [1 10 100 1000 2000])
        end     
        xlabel('Maximum modeled ice thickness (m)')
        ylabel(ylabels{ii})
        text(-200, typ(ii), ['(' letters(ii) ')'], 'fontsize', 20, 'fontweight', 'bold')
        grid on
        axes('position', get(gca, 'position'), 'color', 'none', 'xaxislocation', 'top', 'yaxislocation', 'right', 'fontsize', 20, 'fontweight', 'bold', 'xscale', 'log', 'xdir', 'reverse', 'yticklabel', {})
        xlabel('Maximum plausible center frequency (MHz)')
        axis([1e0 1e3 0 ylims(ii)])
        set(gca, 'xticklabel', {'1' '10' '100' '1000'}, 'linewidth', 1)
        if (ii == 2)
            lg_str                  = RGI_str(idx_RGI_temperate);
            for jj = 1:length(idx_RGI_temperate)
                if (jj == 1)
                    lg_str{jj}     = [lg_str{jj} ' [n = ' num2str(length(idx_F19_large{idx_RGI_temperate(jj)})) ']'];
                else
                    lg_str{jj}     = [lg_str{jj} ' [' num2str(length(idx_F19_large{idx_RGI_temperate(jj)})) ']'];
                end
            end
            lg                      = legend(pm, lg_str);
            set(lg, 'location', 'southeast', 'fontsize', 18)
        end
    end

%%
end