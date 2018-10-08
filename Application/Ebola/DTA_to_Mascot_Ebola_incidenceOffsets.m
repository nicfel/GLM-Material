% Extract covariate information from dta xml
clear

rng(1)


start_predictors = false;
start_locations = false;
all_locations = cell(0,0);
% build dictionary to convert matrix indices
index = 1;
for a = 1 : 56
    for b = a + 1 : 56
        indices(a,b) = index;
        index = index+1;
    end
end

for a = 1 : 56
    for b = a + 1 : 56
        indices(b,a) = index;
        index = index+1;
    end
end



lc = 1;
st = false;


% data of the taxa
taxa_name = cell(0,0);
taxa_date = zeros(0,0);
taxa_count = cell(0,0);
taxa_loc = cell(0,0);
taxa_seq = cell(0,0);
f = fopen('data/Makona_1610_cds_ig.joint_GLM.xml');

while ~feof(f)
    line = fgets(f);
    if ~isempty(strfind(line,'<taxa')) && ~st
        while isempty(strfind(line,'</taxa>'))
            line = fgets(f);       
            % get the id
            if ~isempty(strfind(line,'<taxon'))
                tmp = strsplit(line,'id="');
                tmp = strsplit(tmp{2},'"');
                tmp = tmp{1};
                taxa_name{end+1,1} = tmp;
            end
            % get the sampling time    
            if ~isempty(strfind(line,'<date'))
                tmp = strsplit(line,'value="');
                tmp = strsplit(tmp{2},'"');
                tmp = tmp{1};
                taxa_date(end+1,1) = str2double(tmp);
            end
            % get the country    
            if ~isempty(strfind(line,'<attr name="country">'))
                line = strtrim(fgets(f));
                taxa_count{end+1,1} = line;
            end
            % get the country    
            if ~isempty(strfind(line,'<attr name="location">'))
                line = strtrim(fgets(f));
                taxa_loc{end+1,1} = line;
            end
        end
        st = true;
    end
    
    % read in the sequences
    if ~isempty(strfind(line,'<sequence'))
        line = fgets(f);
        name = strsplit(line,'idref="');
        name = strsplit(name{2},'"');
        
        ind = find(ismember(taxa_name,name{1}));
        
        if ~strcmp(name{1}, taxa_name{ind})
            error('wrong sequence assingment')
        end
        line = fgets(f);
        tmp = '';
        while isempty(strfind(line,'</sequence>'))
            tmp = [tmp strtrim(line)]; 
            line = fgets(f);
        end
        taxa_seq{ind,1} = tmp;
    end
    
    
    % read in the possible locations
    if ~isempty(strfind(line,'<generalDataType id="location.dataType">'))
        start_locations = true;
    end
    if ~isempty(strfind(line,'</generalDataType>'))
        start_locations = false;
    end
    if start_locations
        if  ~isempty(strfind(line,'<state code'))
            tmp = strsplit(line, '"');
            all_locations{end+1,1} = tmp{2};
        end
    end

    
    % read in the predictors
    if ~isempty(strfind(line,'<designMatrix'))
        start_predictors = true;
    end
    if ~isempty(strfind(line,'</designMatrix>'))
        start_predictors = false;
    end
    if start_predictors
        if ~isempty(strfind(line,'parameter'))
            name = regexp(line,'id="(\w*)"','match');
            name = regexprep(name, 'id="(\w*)"','$1');
            name = regexprep(name, '(\d*)(\w*)','$2$1');
            % make a vector of the values
            tmp = regexprep(line, '<parameter id=\"(\w*)\" value=\"','');
            tmp = regexprep(strtrim(tmp), '" />','');
            tmp = regexprep(strtrim(tmp), '"/>','');
            split_tmp = strsplit(strtrim(tmp));
            for i = 1 : length(split_tmp)
                [a,b] = find(indices==i);
                m_covariate.(name{1})(a,b) = str2double(split_tmp{i});
            end
        end
    end
end
fclose(f);



% read in the phylogenetic trees split into parts
f = fopen('data/Ebola.tree');
line = fgets(f);
tmp = strsplit(line, ';');
tree = strrep(tmp{5}, '''','');
fclose(f);
% get all the leafnames
ptree = phytreeread(tree);
all_leafs = get(ptree, 'leafnames');
leafs = cell(0,0);
leaf_locations = cell(0,0);

dont_use = {'WesternArea GIN LBR SLE'};

for i = 1 : length(all_leafs)
    if isempty(strfind(all_leafs{i},'?'));
        ind = find(ismember(taxa_name, all_leafs{i}));
        if taxa_date(ind) < 2015
            tmp = strsplit(all_leafs{i}, '|');
            if length(tmp) == 6
                if isempty(strfind(dont_use{1}, tmp{5}))
                    leafs{end+1,1} = all_leafs{i};
                    leaf_locations{end+1,1} = tmp{5};
                end
            else
                if isempty(strfind(dont_use{1}, tmp{4}))
                    leafs{end+1,1} = all_leafs{i};
                    leaf_locations{end+1,1} = tmp{4};
                end
            end
        end
    end 
end
locations = unique(leaf_locations);
loc_indices = zeros(0,0);

% count the number of samples per location
for i = 1 : length(locations)
    nr_samples(i) = sum(ismember(leaf_locations,locations{i}));
end

% truncate the migration covariates to only include "locations"
for i = 1 : length(locations)
    loc_indices(end+1,1) = find(ismember(all_locations, locations{i}));
end
loc_indices = sort(loc_indices);
don_use_indices = true(56,1);
don_use_indices(loc_indices) = false;
cov_names = fieldnames(m_covariate);
for i = 1 : length(cov_names)    
    truncated_m_covariates.(cov_names{i}) = m_covariate.(cov_names{i});
    truncated_m_covariates.(cov_names{i})(don_use_indices,:) = [];
    truncated_m_covariates.(cov_names{i})(:,don_use_indices) = [];
end

cov_names = fieldnames(truncated_m_covariates);
for i = 1 : length(cov_names) 
    tmp = truncated_m_covariates.(cov_names{i})(:);
    tmp(1:(length(locations)+1):end) = [];
    if length(unique(tmp))==1
        truncated_m_covariates = rmfield(truncated_m_covariates, cov_names{i});
    end
end

% add samples from and to
mat_samples = zeros(length(nr_samples));
for a = 1 : length(nr_samples);
    for b = 1 : length(nr_samples);
        if a~=b
            mat_samples(a,b) = log(nr_samples(a));
        end
    end
end

truncated_m_covariates.Samples_from = mat_samples;
truncated_m_covariates.Samples_to = mat_samples';


% add an unsampled state without covariates
cov_names = fieldnames(truncated_m_covariates);

% make a one dimension vector from the covariates
for i = 1 : length(cov_names);
    c = 1;
    for a = 1 : size(truncated_m_covariates.(cov_names{i}),1)
        for b = 1 : size(truncated_m_covariates.(cov_names{i}),2)
            if a~=b
                m_cov.(cov_names{i})(c) = truncated_m_covariates.(cov_names{i})(a,b);
                c = c+1;
            end
        end
    end
end

% renormalize covariates
for i = 1 : length(cov_names);
    if length(unique(m_cov.(cov_names{i})))>2
        m_cov.(cov_names{i}) = (m_cov.(cov_names{i}) - mean(m_cov.(cov_names{i})))/std(m_cov.(cov_names{i}));
    end
end

% get the covariates for the effective population sizes
cov_names = fieldnames(m_cov);
for i = 1 : length(cov_names);
    if ~isempty(strfind(cov_names{i},'origin'))
        Ne_cov.(strrep(cov_names{i},'origin','')) = max(truncated_m_covariates.(cov_names{i})')';
    end
end
ne_cov_names = fieldnames(Ne_cov);

% renormalize covariates
for i = 1 : length(ne_cov_names);
    Ne_cov.(ne_cov_names{i}) = (Ne_cov.(ne_cov_names{i}) - mean(Ne_cov.(ne_cov_names{i})))/std(Ne_cov.(ne_cov_names{i}));
end


% get the subtree with the used leafs

tree_leafnames = get(ptree, 'leafnames');

sel = ismember(tree_leafnames,leafs);
subtree =  prune(ptree,~sel);

% import the number of cases over time
csv = importdata('data/EBOV_maxCases.csv');
start_day = csv.textdata(1,4:end);
csv.data(csv.data==0) = 0.1;

ind = find(ismember(taxa_name,leafs));
mr_sample = max(taxa_date(ind));

% convert the epi weeks into decimals
for i = 1 : length(start_day)
    tmp = strsplit(start_day{i},'-');
    decimal_day(i) = -mr_sample + str2double(tmp{1}) + (datenum(start_day{i},'yyyy-mmm-dd') - datenum(tmp{1}, 'yyyy'))/(datenum(num2str(str2double(tmp{1})+1), 'yyyy') -datenum(tmp{1}, 'yyyy'));
end
decimal_day = -decimal_day;
first_interval = find(decimal_day<0,1)-1;
decimal_day(first_interval+1:end) = [];

%

% get the indices
indices = zeros(0,0);
for a = 1 : length(locations)
    indices(a,1) = find(ismember(csv.textdata(:,3), locations{a})) - 1;
end

% get the incidences over time and the incidence ratios over time
inc_Ne = zeros(0,0);
inc_m = zeros(0,0);

for t = first_interval : -1 : 1
    for a = 1 : length(indices)
        inc_Ne(end+1,1) = csv.data(indices(a),t);
        for b = 1 : length(indices)
            if a~=b
                inc_m(end+1,1) = csv.data(indices(a),t)/csv.data(indices(b),t);
            end
        end
    end
end

% add the covariates and multiply the constant ones
cov_names = fieldnames(m_cov);
for i = 1 : length(cov_names);
    m = repmat(m_cov.(cov_names{i})', 1, first_interval);
    m_t_cov.(cov_names{i}) = m(:);
end

lme = log(inc_m);


m_t_cov.incidenceRatio = lme-mean(lme(~isinf(lme)));

% add the covariates and multiply the constant ones
cov_names = fieldnames(Ne_cov);
for i = 1 : length(cov_names);
    m = repmat(Ne_cov.(cov_names{i})', 1, first_interval);
    Ne_t_cov.(cov_names{i}) = m(:);
end


lne = log(inc_Ne);
Ne_t_cov.incidence = lne-mean(lne(~isinf(lne)));



%%
% add incidence with offsets of weeks, but do not renormalize!
offset = [3 6 9];
for i = 1 : length(offset)
    tmp_incidence_earlier = [ones(length(locations)*offset(i),1)*min(Ne_t_cov.incidence);
        Ne_t_cov.incidence(1:end-length(locations)*offset(i))];
    
    tmp_incidence_later = [Ne_t_cov.incidence(length(locations)*offset(i)+1:end);
         ones(length(locations)*offset(i),1)*min(Ne_t_cov.incidence)];       
%     
    Ne_t_cov.(sprintf('incidence%dweeklater',offset(i))) = tmp_incidence_later;
    Ne_t_cov.(sprintf('incidence%dweekearlier',offset(i))) = tmp_incidence_earlier;
    
    plot(Ne_t_cov.(sprintf('incidence%dweeklater',offset(i)))(1:14:end)); hold on

end


%%
system('rm -r xmls/mascot');
system('mkdir xmls/mascot');


for utn = 1 : 10
    i = utn;
    print_tree = getnewickstr(subtree);

    % make tripletts of all runs with different random initial values
    for tr = 1 : 1    
        % make the xmls for the structcoal
        flog = ['Ebola_glm_mascot_rep' num2str(i-1)];
        fname = sprintf('xmls/mascot/%s.xml',flog);
        g = fopen(fname,'w');


        fprintf(g,'<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate=''Standard'' beautistatus='''' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.mascot.dynamics:beast.mascot.distribution:beast.mascot.logger" version="2.0">\n');

        fprintf(g,'\t<data id="sequences" name="alignment">\n');
        for j = 1 : length(leafs)
            seq_ind = find(ismember(taxa_name, leafs{j}));
            fprintf(g,'\t\t<sequence id="%s" taxon="%s" totalcount="4" value="%s"/>\n',leafs{j},leafs{j}, taxa_seq{seq_ind});
        end
        fprintf(g,'\t</data>\n');
        fprintf(g,'<map name="prior">beast.math.distributions.Prior</map>\n');
        fprintf(g,'<map name="maxRate">beast.mascot.glmmodel.MaxRate</map>\n');

        fprintf(g,'\t<run id="mcmc" spec="MCMC" chainLength="5000000">\n');
        fprintf(g,'\t\t<state id="state" storeEvery="200000">\n');
        fprintf(g,'\t\t\t<stateNode id="tree" spec="beast.app.mascot.beauti.TreeWithTrait">\n');
        fprintf(g,'\t\t\t\t<typeTrait id="typeTraitSet.t" spec="beast.evolution.tree.TraitSet" traitname="type" value="');
        for j = 1 : length(leafs)-1
            seq_ind = find(ismember(taxa_name, leafs{j}));
            fprintf(g,'%s=%s,',leafs{j},taxa_loc{seq_ind});
        end
        seq_ind = find(ismember(taxa_name, leafs{length(leafs)}));
        fprintf(g,'%s=%s">\n',leafs{length(leafs)},taxa_loc{seq_ind});

        fprintf(g,'\t\t\t\t\t<taxa id="TaxonSet.0" spec="TaxonSet">\n');
        fprintf(g,'\t\t\t\t\t\t<alignment idref="sequences"/>\n');
        fprintf(g,'\t\t\t\t\t</taxa>\n');
        fprintf(g,'\t\t\t\t</typeTrait>\n');
        fprintf(g,'\t\t\t</stateNode>\n');        
        fprintf(g,'\t\t\t<parameter id="migrationScaler" name="stateNode" dimension="0">0</parameter>\n');
        fprintf(g,'\t\t\t<stateNode id="migrationIndicator" spec="parameter.BooleanParameter" dimension="0">true</stateNode>\n');
        start_val = lognrnd(0,1,1);
        while start_val > 5
            start_val = lognrnd(0,1,1);
        end
        fprintf(g,'\t\t\t<parameter id="migrationClock" name="stateNode" dimension="0">%.12f</parameter>\n', start_val);
%         fprintf(g,'\t\t\t<parameter id="migrationError" name="stateNode" dimension="0">0</parameter>\n');

        
        
        fprintf(g,'\t\t\t<parameter id="NeScaler" name="stateNode" dimension="0">0</parameter>\n');
        fprintf(g,'\t\t\t<stateNode id="NeIndicator" spec="parameter.BooleanParameter" dimension="0">true</stateNode>\n');
        start_val = lognrnd(0,1,1);
        while start_val > 5
            start_val = lognrnd(0,1,1);
        end
        fprintf(g,'\t\t\t<parameter id="NeClock" name="stateNode" dimension="0">%.12f</parameter>\n', start_val);
%         fprintf(g,'\t\t\t<parameter id="NeError" name="stateNode" dimension="0">0</parameter>\n');

        fprintf(g,'\t\t\t<parameter id="cds.CP1.mu" name="stateNode">0.536</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="cds.CP2.mu" name="stateNode">0.503</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="cds.CP3.mu" name="stateNode">1.412</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="ig.mu" name="stateNode">1.597</parameter>\n');
               
        fprintf(g,'\t\t\t<parameter id="cds.CP1.kappa" name="stateNode">9.016</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="cds.CP2.kappa" name="stateNode">12.929</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="cds.CP3.kappa" name="stateNode">13.646</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="ig.kappa" name="stateNode">8.735</parameter>\n');

        fprintf(g,'\t\t\t<parameter id="cds.CP1.alpha" name="stateNode">0.169</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="cds.CP2.alpha" name="stateNode">5.424E-2</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="cds.CP3.alpha" name="stateNode">0.605</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="ig.alpha" name="stateNode">0.354</parameter>\n');
        
        fprintf(g,'\t\t\t<parameter id="clockRate.c" name="stateNode">1.21E-3</parameter>\n');

        fprintf(g,'\t\t\t<parameter id="freqParameter.s" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n');
        fprintf(g,'\t\t</state>\n');



        fprintf(g,'\t\t<init spec="beast.util.TreeParser" id="NewickTree.t:Species" adjustTipHeights="false"\n');
        fprintf(g,'\t\t\tinitial="@tree" taxa="@sequences"\n');
        fprintf(g,'\t\t\tIsLabelledNewick="true"\n');
        fprintf(g,'\t\t\tnewick="%s"/>\n',print_tree);
        fprintf(g,' \t\t<distribution id="posterior" spec="util.CompoundDistribution">\n');
        fprintf(g,'\t\t\t<distribution id="prior" spec="util.CompoundDistribution">\n');

        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@NeClock">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.OneOnX"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@migrationClock">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.OneOnX"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        
        
%         fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@NeError">\n');
%         fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal" mean="0" sigma="0.5"/>\n');
%         fprintf(g,'\t\t\t\t</distribution>\n');
%         fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@migrationError">\n');
%         fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal" mean="0" sigma="0.5"/>\n');
%         fprintf(g,'\t\t\t\t</distribution>\n');


        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@migrationScaler">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal" mean="0" sigma="1"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@NeScaler">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal" mean="0" sigma="1"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
%         fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@kappa.s">\n');
%         fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.LogNormalDistributionModel" M="1" S="1.25"/>\n');
%         fprintf(g,'\t\t\t\t</distribution>\n');

%         fprintf(g,'\t\t\t\t<maxRate id="maxRate" name="distribution">\n');
%         fprintf(g,'\t\t\t\t\t<GLMmodel idref="constant"/>\n');
%         fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Uniform" upper=''100''/>\n');
%         fprintf(g,'\t\t\t\t</maxRate>\n');
        
%         fprintf(g,'\t\t\t\t<maxRate id="rateExponential" name="distribution" migrationOnly="true">\n');
%         fprintf(g,'\t\t\t\t\t<GLMmodel idref="constant"/>\n');
%         fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Exponential" mean=''0.1''/>\n');
%         fprintf(g,'\t\t\t\t</maxRate>\n');



        fprintf(g,'\t\t\t\t<prior id="nonZeroRatePrior.Migration" name="distribution">\n');
        fprintf(g,'\t\t\t\t\t<x id="nonZeroCovarates.Migration" spec="util.Sum">\n');
        fprintf(g,'\t\t\t\t\t\t<arg idref="migrationIndicator"/>\n');
        fprintf(g,'\t\t\t\t\t</x>\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Poisson" lambda=''3''/>\n');
        fprintf(g,'\t\t\t\t</prior>\n');

        fprintf(g,'\t\t\t\t<prior id="nonZeroRatePrior.Ne" name="distribution">\n');
        fprintf(g,'\t\t\t\t\t<x id="nonZeroCovarates.Ne" spec="util.Sum">\n');
        fprintf(g,'\t\t\t\t\t\t<arg idref="NeIndicator"/>\n');
        fprintf(g,'\t\t\t\t\t</x>\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Poisson" lambda=''3''/>\n');
        fprintf(g,'\t\t\t\t</prior>\n');

        
        fprintf(g,'\t\t\t\t<distribution id="coalescent" spec="Mascot">\n');
        fprintf(g,'\t\t\t\t\t<structuredTreeIntervals spec="StructuredTreeIntervals" id="TreeIntervals" tree="@tree"/>\n');
        fprintf(g,'\t\t\t\t\t<dynamics spec="GLM" id="constant" typeTrait="@typeTraitSet.t" maxRate="100" types="Bo Bombali Bonthe Kailahun Kambia Kenema Koinadugu Kono Moyamba PortLoko Pujehun Tonkolili WesternRural WesternUrban">\n');
        fprintf(g,'\t\t\t\t\t\t<rateShifts spec="beast.core.parameter.RealParameter" id="relativeRateShifts">%s</rateShifts>\n', sprintf('%f ',sort(decimal_day)));
        fprintf(g,'\t\t\t\t\t\t<migrationGLM spec="beast.mascot.glmmodel.LogLinear" id="migrationGLM">\n');

        cov_names = fieldnames(m_t_cov);
        for k = 1 : length(cov_names)
            fprintf(g,'\t\t\t\t\t\t\t<covariates spec="beast.core.parameter.RealParameter" id="%s">%s</covariates>\n', cov_names{k}, sprintf('%f ',m_t_cov.(cov_names{k})));
        end
        fprintf(g,'\t\t\t\t\t\t\t<scaler idref="migrationScaler"/>\n');
        fprintf(g,'\t\t\t\t\t\t\t<indicator idref="migrationIndicator"/>\n');
        fprintf(g,'\t\t\t\t\t\t\t<clock idref="migrationClock"/>\n');
%         fprintf(g,'\t\t\t\t\t\t\t<error idref="migrationError"/>\n');
        fprintf(g,'\t\t\t\t\t\t</migrationGLM>\n');
        fprintf(g,'\t\t\t\t\t\t<NeGLM spec="beast.mascot.glmmodel.LogLinear" id="NeGLM">\n');
        cov_names = fieldnames(Ne_t_cov);
        for k = 1 : length(cov_names)
            fprintf(g,'\t\t\t\t\t\t\t<covariates spec="beast.core.parameter.RealParameter" id="%s">%s</covariates>\n', cov_names{k}, sprintf('%f ',Ne_t_cov.(cov_names{k})));
        end
        fprintf(g,'\t\t\t\t\t\t\t<scaler idref="NeScaler"/>\n');
        fprintf(g,'\t\t\t\t\t\t\t<indicator idref="NeIndicator"/>\n');
        fprintf(g,'\t\t\t\t\t\t\t<clock idref="NeClock"/>\n');
%         fprintf(g,'\t\t\t\t\t\t\t<error idref="NeError"/>\n');
        fprintf(g,'\t\t\t\t\t\t</NeGLM>\n');
        fprintf(g,'\t\t\t\t\t</dynamics>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');

        fprintf(g,'\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t<distribution id="likelihood" spec="util.CompoundDistribution">\n');
        
        % CP1
        fprintf(g,'\t\t\t\t<distribution id="treeLikelihood.CP1" spec="TreeLikelihood" tree="@tree">\n');
        fprintf(g,'\t\t\t\t\t<data id="sequences.CP1" spec="FilteredAlignment" filter="1:14517:3">\n');
		fprintf(g,'\t\t\t\t\t\t<data idref="sequences"/>\n');
		fprintf(g,'\t\t\t\t\t</data>\n');
        fprintf(g,'\t\t\t\t\t<siteModel id="SiteModel.CP1" spec="SiteModel" gammaCategoryCount="4" shape="@cds.CP1.alpha">\n');
        fprintf(g,'\t\t\t\t\t\t<mutationRate idref="cds.CP1.mu"/>\n');
        fprintf(g,'\t\t\t\t\t\t<parameter id="proportionInvariant.CP1" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n');
        fprintf(g,'\t\t\t\t\t\t<substModel id="hky.CP1" spec="HKY" kappa="@cds.CP1.kappa">\n');
        fprintf(g,'\t\t\t\t\t\t\t<frequencies id="empiricalFreqs.CP1" spec="Frequencies" data="@sequences.CP1"/>\n');
        fprintf(g,'\t\t\t\t\t\t</substModel>\n');
        fprintf(g,'\t\t\t\t\t</siteModel>\n');
        fprintf(g,'\t\t\t\t\t<branchRateModel id="StrictClock.CP1" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        
        % CP2
        fprintf(g,'\t\t\t\t<distribution id="treeLikelihood.CP2" spec="TreeLikelihood" tree="@tree">\n');
        fprintf(g,'\t\t\t\t\t<data id="sequences.CP2" spec="FilteredAlignment" filter="2:14517:3">\n');
		fprintf(g,'\t\t\t\t\t\t<data idref="sequences"/>\n');
		fprintf(g,'\t\t\t\t\t</data>\n');
        fprintf(g,'\t\t\t\t\t<siteModel id="SiteModel.CP2" spec="SiteModel" gammaCategoryCount="4" shape="@cds.CP2.alpha">\n');
        fprintf(g,'\t\t\t\t\t\t<mutationRate idref="cds.CP2.mu"/>\n');
        fprintf(g,'\t\t\t\t\t\t<parameter id="proportionInvariant.CP2" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n');
        fprintf(g,'\t\t\t\t\t\t<substModel id="hky.CP2" spec="HKY" kappa="@cds.CP2.kappa">\n');
        fprintf(g,'\t\t\t\t\t\t\t<frequencies id="empiricalFreqs.CP2" spec="Frequencies" data="@sequences.CP2"/>\n');
        fprintf(g,'\t\t\t\t\t\t</substModel>\n');
        fprintf(g,'\t\t\t\t\t</siteModel>\n');
        fprintf(g,'\t\t\t\t\t<branchRateModel id="StrictClock.CP2" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        
        % CP3
        fprintf(g,'\t\t\t\t<distribution id="treeLikelihood.CP3" spec="TreeLikelihood" tree="@tree">\n');
        fprintf(g,'\t\t\t\t\t<data id="sequences.CP3" spec="FilteredAlignment" filter="3:14517:3">\n');
		fprintf(g,'\t\t\t\t\t\t<data idref="sequences"/>\n');
		fprintf(g,'\t\t\t\t\t</data>\n');
        fprintf(g,'\t\t\t\t\t<siteModel id="SiteModel.CP3" spec="SiteModel" gammaCategoryCount="4" shape="@cds.CP3.alpha">\n');
        fprintf(g,'\t\t\t\t\t\t<mutationRate idref="cds.CP3.mu"/>\n');
        fprintf(g,'\t\t\t\t\t\t<parameter id="proportionInvariant.CP3" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n');
        fprintf(g,'\t\t\t\t\t\t<substModel id="hky.CP3" spec="HKY" kappa="@cds.CP3.kappa">\n');
        fprintf(g,'\t\t\t\t\t\t\t<frequencies id="empiricalFreqs.CP3" spec="Frequencies" data="@sequences.CP3"/>\n');
        fprintf(g,'\t\t\t\t\t\t</substModel>\n');
        fprintf(g,'\t\t\t\t\t</siteModel>\n');
        fprintf(g,'\t\t\t\t\t<branchRateModel id="StrictClock.CP3" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        
        % Ig
        fprintf(g,'\t\t\t\t<distribution id="treeLikelihood.ig" spec="TreeLikelihood" tree="@tree">\n');
        fprintf(g,'\t\t\t\t\t<data id="sequences.ig" spec="FilteredAlignment" filter="14518::1">\n');
		fprintf(g,'\t\t\t\t\t\t<data idref="sequences"/>\n');
		fprintf(g,'\t\t\t\t\t</data>\n');
        fprintf(g,'\t\t\t\t\t<siteModel id="SiteModel.ig" spec="SiteModel" gammaCategoryCount="4" shape="@ig.alpha">\n');
        fprintf(g,'\t\t\t\t\t\t<mutationRate idref="ig.mu"/>\n');
        fprintf(g,'\t\t\t\t\t\t<parameter id="proportionInvariant.ig" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n');
        fprintf(g,'\t\t\t\t\t\t<substModel id="hky.ig" spec="HKY" kappa="@ig.kappa">\n');
        fprintf(g,'\t\t\t\t\t\t\t<frequencies id="empiricalFreqs.ig" spec="Frequencies" data="@sequences.ig"/>\n');
        fprintf(g,'\t\t\t\t\t\t</substModel>\n');
        fprintf(g,'\t\t\t\t\t</siteModel>\n');
        fprintf(g,'\t\t\t\t\t<branchRateModel id="StrictClock.ig" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');

       
        fprintf(g,'\t\t\t</distribution>\n');


        
        
        fprintf(g,'\t\t</distribution>\n');
        
%         fprintf(g,'\t\t<operator id="FrequenciesExchanger.s" spec="DeltaExchangeOperator" delta="0.01" weight="1">\n');
%         fprintf(g,'\t\t\t<parameter idref="freqParameter.s"/>\n');
%         fprintf(g,'\t\t</operator>\n');
%         fprintf(g,'\t\t<operator id="KappaScaler.s" spec="ScaleOperator" parameter="@kappa.s" scaleFactor="0.5" weight="1"/>\n');
        fprintf(g,'\t\t<operator id="TreeScaler.t" spec="ScaleOperator" tree="@tree" scaleFactor="0.9" optimise="false" weight="1.0"/>\n');
        fprintf(g,'\t\t<operator id="RootScaler.t" spec="ScaleOperator" tree="@tree" scaleFactor="0.9" rootOnly="true" optimise="false" weight="1.0"/>\n');
        fprintf(g,'\t\t<operator id="Uniform.t" spec="Uniform" tree="@tree" weight="1.0"/>\n');
        fprintf(g,'\t\t<operator id="SubtreeSlide.t" spec="SubtreeSlide" tree="@tree" size="0.05" weight="20.0"/>\n');
        fprintf(g,'\t\t<operator id="Narrow.t" spec="Exchange" tree="@tree" weight="5.0"/>\n');
        fprintf(g,'\t\t<operator id="Wide.t" spec="Exchange" tree="@tree" isNarrow="false" weight="30.0"/>\n');
        fprintf(g,'\t\t<operator id="WilsonBalding.t" spec="WilsonBalding" tree="@tree" weight="30.0"/>\n');
%         fprintf(g,'\t\t<operator id="NeMigUpDownOperator.t" spec="UpDownOperator" scaleFactor="0.98" weight="1">\n');
%         fprintf(g,'\t\t\t<up idref="tree"/>\n');
%         fprintf(g,'\t\t\t<up idref="NeClock"/>\n');
%         fprintf(g,'\t\t\t<down idref="migrationClock"/>\n');
%         fprintf(g,'\t\t</operator>\n');

        
        
        fprintf(g,'\t\t<operator id="NeClockScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@NeClock" weight="1.0"/>\n');
		fprintf(g,'\t\t<operator id="MigrationClockScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@migrationClock" weight="1.0"/>\n');

        
		fprintf(g,'\t\t<operator id="indicatorFlip.s:Migration" spec="BitFlipOperator" parameter="@migrationIndicator" weight="10"/>\n');
% 		fprintf(g,'\t\t<operator id="BSSVSoperator.c:Migration" spec="BitFlipBSSVSOperator" indicator="@migrationIndicator" mu="@migrationClock" weight="30"/>\n');
		fprintf(g,'\t\t<operator id="indicatorFlip.s:Ne" spec="BitFlipOperator" parameter="@NeIndicator" weight="10"/>\n');
% 		fprintf(g,'\t\t<operator id="BSSVSoperator.c:Ne" spec="BitFlipBSSVSOperator" indicator="@NeIndicator" mu="@NeClock" weight="30"/>\n');

		fprintf(g,'\t\t<operator id="NeScalerScaler" spec="beast.mascot.operators.MultiRealRandomWalkOperator" windowSize="0.02" parameter="@NeScaler" useGaussian="true" weight="20.0"/>\n');
        fprintf(g,'\t\t<operator id="MigrationScalerScaler" spec="beast.mascot.operators.MultiRealRandomWalkOperator" windowSize="0.02" parameter="@migrationScaler" useGaussian="true" weight="20.0"/>\n');

 		fprintf(g,'\t\t<operator id="indicatorSwap.s:Ne" spec="beast.mascot.operators.BooleanSwapOperator" indicator="@NeIndicator" weight="10"/>\n');
 		fprintf(g,'\t\t<operator id="indicatorSwap.s:Migration" spec="beast.mascot.operators.BooleanSwapOperator" indicator="@migrationIndicator" weight="10"/>\n');

        fprintf(g,'\t\t<operator id="JointSwap.s:Ne" spec="beast.mascot.operators.BooleanSwapOperator" indicator="@NeIndicator" parameter="@NeScaler" weight="10"/>\n');
 		fprintf(g,'\t\t<operator id="JointSwap.s:Migration" spec="beast.mascot.operators.BooleanSwapOperator" indicator="@migrationIndicator" parameter="@migrationScaler"  weight="10"/>\n');

        
        fprintf(g,'\t\t<logger id="probtreelog" fileName="%s.trees" logEvery="200000" mode="tree">\n',flog);
        fprintf(g,'\t\t\t<log id="logTrees" spec="StructuredTreeLogger" mascot="@coalescent"/>\n');
        fprintf(g,'\t\t</logger>\n');


        fprintf(g,'\t\t<logger id="tracelog" fileName="%s.log" logEvery="200000" model="@posterior" sanitiseHeaders="true" sort="smart">\n',flog);
        fprintf(g,'\t\t\t<log idref="posterior"/>\n');
        fprintf(g,'\t\t\t<log idref="prior"/>\n');
        fprintf(g,'\t\t\t<log idref="treeLikelihood.CP1"/>\n');
        fprintf(g,'\t\t\t<log idref="treeLikelihood.CP2"/>\n');
        fprintf(g,'\t\t\t<log idref="treeLikelihood.CP3"/>\n');
        fprintf(g,'\t\t\t<log idref="treeLikelihood.ig"/>\n');
        fprintf(g,'\t\t\t<log idref="constant"/>\n');
        fprintf(g,'\t\t\t<log idref="nonZeroCovarates.Migration"/>\n');
        fprintf(g,'\t\t\t<log idref="nonZeroCovarates.Ne"/>\n');
        fprintf(g,'\t\t\t<log idref="migrationClock"/>\n');
        fprintf(g,'\t\t\t<log idref="NeClock"/>\n');
        fprintf(g,'\t\t\t<log idref="migrationGLM"/>\n');
        fprintf(g,'\t\t\t<log idref="NeGLM"/>\n');
%         fprintf(g,'\t\t\t<log idref="migrationScaler"/>\n');
%         fprintf(g,'\t\t\t<log idref="migrationIndicator"/>\n');
%         fprintf(g,'\t\t\t<log idref="NeScaler"/>\n');
%         fprintf(g,'\t\t\t<log idref="NeIndicator"/>\n');
        fprintf(g,'\t\t</logger>\n');
        fprintf(g,'\t\t<logger id="screenlog" logEvery="20000">\n');
        fprintf(g,'\t\t\t<log idref="posterior"/>\n');
        fprintf(g,'\t\t</logger>\n');
        fprintf(g,'\t</run>\n');
        fprintf(g,'</beast>\n');
        fclose(f);
    end
end