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

% build the mapping matrix to get correct order of predictors
index=1;
for a = 1 : length(locations)
    for b = a + 1 : length(locations)
        mapper(a,b) = index;
        index = index+1;
    end
end

for a = 1 : length(locations)
    for b = a + 1 : length(locations)
        mapper(b,a) = index;
        index = index+1;
    end
end



% make a one dimension vector from the covariates
for i = 1 : length(cov_names);
    c = 1;
    for a = 1 : size(truncated_m_covariates.(cov_names{i}),1)
        for b = 1 : size(truncated_m_covariates.(cov_names{i}),2)
            if a~=b
                m_cov.(cov_names{i})(mapper(a,b)) = truncated_m_covariates.(cov_names{i})(a,b);
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


% get the subtree with the used leafs
tree_leafnames = get(ptree, 'leafnames');

sel = ismember(tree_leafnames,leafs);
subtree =  prune(ptree,~sel);



%%
system('rm -r xmls/dta');
system('mkdir xmls/dta');


for utn = 1 : 5
    i = utn;
    print_tree = getnewickstr(subtree);

    % make tripletts of all runs with different random initial values
    for tr = 1 : 1    
        % make the xmls for the structcoal
        flog = ['Ebola_glm_dta_rep' num2str(i)];
        fname = sprintf('xmls/dta/%s.xml',flog);
        g = fopen(fname,'w');


        
        
        
        fprintf(g,'<?xml version="1.0" standalone="yes"?>\n');
        
        fprintf(g,'<beast>\n');
        fprintf(g,'\t<taxa id="taxa">\n');
        for j = 1 : length(leafs)
            seq_ind = find(ismember(taxa_name, leafs{j}));
            fprintf(g,'\t\t<taxon id="%s">\n',leafs{j});
            fprintf(g,'\t\t\t<date value="%f" direction="forwards" units="years"/>\n', taxa_date(seq_ind));
            fprintf(g,'\t\t\t<attr name="location">\n');
            fprintf(g,'\t\t\t\t%s\n', taxa_loc{seq_ind});
            fprintf(g,'\t\t\t</attr>\n');
            fprintf(g,'\t\t</taxon>\n');
        end
        fprintf(g,'\t</taxa>\n');
        fprintf(g,'\t<alignment id="alignment" dataType="nucleotide">\n');

        for j = 1 : length(leafs)
            seq_ind = find(ismember(taxa_name, leafs{j}));
            fprintf(g,'\t\t<sequence>\n');
            fprintf(g,'\t\t\t<taxon idref="%s"/>\n', leafs{j});
            fprintf(g,'\t\t\t\t%s\n', taxa_seq{seq_ind});
            fprintf(g,'\t\t</sequence>\n');
        end
        
        fprintf(g,'\t</alignment>\n');

        fprintf(g,'\t<patterns id="cds.CP1.patterns" from="1" to="14517" every="3" strip="false">\n');
        fprintf(g,'\t\t<alignment idref="alignment"/>\n');
        fprintf(g,'\t</patterns>\n');
        fprintf(g,'\t<patterns id="cds.CP2.patterns" from="2" to="14517" every="3" strip="false">\n');
        fprintf(g,'\t\t<alignment idref="alignment"/>\n');
        fprintf(g,'\t</patterns>\n');
        fprintf(g,'\t<patterns id="cds.CP3.patterns" from="3" to="14517" every="3" strip="false">\n');
        fprintf(g,'\t\t<alignment idref="alignment"/>\n');
        fprintf(g,'\t</patterns>\n');
        fprintf(g,'\t<patterns id="ig.patterns" from="14518" strip="false">\n');
        fprintf(g,'\t\t<alignment idref="alignment"/>\n');
        fprintf(g,'\t</patterns>\n');
        fprintf(g,'\t<constantSize id="initialDemo" units="years">\n');
        fprintf(g,'\t\t<populationSize>\n');
        fprintf(g,'\t\t\t<parameter id="initialDemo.popSize" value="0.5"/>\n');
        fprintf(g,'\t\t</populationSize>\n');
        fprintf(g,'\t</constantSize>\n');
        fprintf(g,'\t<newick id="startingTree">\n');
        fprintf(g,'\t%s\n',print_tree);
        fprintf(g,'\t</newick>\n');
        
        fprintf(g,'\t<treeModel id="treeModel">\n');
        fprintf(g,'\t\t<coalescentTree idref="startingTree"/>\n');
        fprintf(g,'\t\t<rootHeight>\n');
        fprintf(g,'\t\t\t<parameter id="treeModel.rootHeight"/>\n');
        fprintf(g,'\t\t</rootHeight>\n');
        fprintf(g,'\t\t<nodeHeights internalNodes="true">\n');
        fprintf(g,'\t\t\t<parameter id="treeModel.internalNodeHeights"/>\n');
        fprintf(g,'\t\t</nodeHeights>\n');
        fprintf(g,'\t\t<nodeHeights internalNodes="true" rootNode="true">\n');
        fprintf(g,'\t\t\t<parameter id="treeModel.allInternalNodeHeights"/>\n');
        fprintf(g,'\t\t</nodeHeights>\n');
        fprintf(g,'\t</treeModel>\n');
        fprintf(g,'\t<gmrfSkyGridLikelihood id="skygrid">\n');
        fprintf(g,'\t\t<populationSizes>\n');
        fprintf(g,'\t\t\t<parameter id="skygrid.logPopSize" dimension="100" value="0"/>\n');
        fprintf(g,'\t\t</populationSizes>\n');
        fprintf(g,'\t\t<precisionParameter>\n');
        fprintf(g,'\t\t\t<parameter id="skygrid.precision" value="0.1" lower="0.0"/>\n');
        fprintf(g,'\t\t</precisionParameter>\n');
        fprintf(g,'\t\t<numGridPoints>\n');
        fprintf(g,'\t\t\t<parameter id="skygrid.numGridPoints" value="99.0"/>\n');
        fprintf(g,'\t\t</numGridPoints>\n');
        fprintf(g,'\t\t<cutOff>\n');
        fprintf(g,'\t\t\t<parameter id="skygrid.cutOff" value="1.5"/>\n');
        fprintf(g,'\t\t</cutOff>\n');
        fprintf(g,'\t\t<populationTree>\n');
        fprintf(g,'\t\t\t<treeModel idref="treeModel"/>\n');
        fprintf(g,'\t\t</populationTree>\n');
        fprintf(g,'\t</gmrfSkyGridLikelihood>\n');
        
        
        fprintf(f, '\t<strictClockBranchRates id="branchRates">\n');
        fprintf(f, '\t\t<rate>\n');
        fprintf(f, '\t\t\t<parameter id="Ebola.clock.rate" value="0.0036" lower="0.0"/>\n');
        fprintf(f, '\t\t</rate>\n');
        fprintf(f, '\t</strictClockBranchRates>\n');

        fprintf(f, '\n');
        fprintf(f, '\t<HKYModel id="cds.CP1.hky">\n');
        fprintf(f, '\t\t<frequencies>\n');
        fprintf(f, '\t\t\t<frequencyModel dataType="nucleotide">\n');
        fprintf(f, '\t\t\t\t<patterns idref="cds.CP1.patterns"/>\n');
        fprintf(f, '\t\t\t\t<frequencies>\n');
        fprintf(f, '\t\t\t\t\t<parameter id="cds.CP1.frequencies" dimension="4"/>\n');
        fprintf(f, '\t\t\t\t</frequencies>\n');
        fprintf(f, '\t\t\t</frequencyModel>\n');
        fprintf(f, '\t\t</frequencies>\n');
        fprintf(f, '\t\t<kappa>\n');
        fprintf(f, '\t\t\t<parameter id="cds.CP1.kappa" lower="0.0" value="9.016"/>\n');
        fprintf(f, '\t\t</kappa>\n');
        fprintf(f, '\t</HKYModel>\n');
        fprintf(f, '\n');
        fprintf(f, '\t<HKYModel id="cds.CP2.hky">\n');
        fprintf(f, '\t\t<frequencies>\n');
        fprintf(f, '\t\t\t<frequencyModel dataType="nucleotide">\n');
        fprintf(f, '\t\t\t\t<patterns idref="cds.CP2.patterns"/>\n');
        fprintf(f, '\t\t\t\t<frequencies>\n');
        fprintf(f, '\t\t\t\t\t<parameter id="cds.CP2.frequencies" dimension="4"/>\n');
        fprintf(f, '\t\t\t\t</frequencies>\n');
        fprintf(f, '\t\t\t</frequencyModel>\n');
        fprintf(f, '\t\t</frequencies>\n');
        fprintf(f, '\t\t<kappa>\n');
        fprintf(f, '\t\t\t<parameter id="cds.CP2.kappa" lower="0.0" value="12.929"/>\n');
        fprintf(f, '\t\t</kappa>\n');
        fprintf(f, '\t</HKYModel>\n');
        fprintf(f, '\n');
        fprintf(f, '\t<HKYModel id="cds.CP3.hky">\n');
        fprintf(f, '\t\t<frequencies>\n');
        fprintf(f, '\t\t\t<frequencyModel dataType="nucleotide">\n');
        fprintf(f, '\t\t\t\t<patterns idref="cds.CP3.patterns"/>\n');
        fprintf(f, '\t\t\t\t<frequencies>\n');
        fprintf(f, '\t\t\t\t\t<parameter id="cds.CP3.frequencies" dimension="4"/>\n');
        fprintf(f, '\t\t\t\t</frequencies>\n');
        fprintf(f, '\t\t\t</frequencyModel>\n');
        fprintf(f, '\t\t</frequencies>\n');
        fprintf(f, '\t\t<kappa>\n');
        fprintf(f, '\t\t\t<parameter id="cds.CP3.kappa" lower="0.0" value="13.646"/>\n');
        fprintf(f, '\t\t</kappa>\n');
        fprintf(f, '\t</HKYModel>\n');
        fprintf(f, '\n');    
        fprintf(f, '\t<HKYModel id="ig.hky">\n');
        fprintf(f, '\t\t<frequencies>\n');
        fprintf(f, '\t\t\t<frequencyModel dataType="nucleotide">\n');
        fprintf(f, '\t\t\t\t<patterns idref="ig.patterns"/>\n');
        fprintf(f, '\t\t\t\t<frequencies>\n');
        fprintf(f, '\t\t\t\t\t<parameter id="ig.frequencies" dimension="4"/>\n');
        fprintf(f, '\t\t\t\t</frequencies>\n');
        fprintf(f, '\t\t\t</frequencyModel>\n');
        fprintf(f, '\t\t</frequencies>\n');
        fprintf(f, '\t\t<kappa>\n');
        fprintf(f, '\t\t\t<parameter id="ig.kappa" lower="0.0" value="8.735"/>\n');
        fprintf(f, '\t\t</kappa>\n');
        fprintf(f, '\t</HKYModel>\n');
        fprintf(f, '\n');
        fprintf(f, '\t<!-- site model                                                              -->\n');
        fprintf(f, '\t<siteModel id="cds.CP1.siteModel">\n');
        fprintf(f, '\t\t<substitutionModel>\n');
        fprintf(f, '\t\t\t<HKYModel idref="cds.CP1.hky"/>\n');
        fprintf(f, '\t\t</substitutionModel>\n');
        fprintf(f, '\t\t<relativeRate>\n');
        fprintf(f, '\t\t\t<parameter id="cds.CP1.mu" lower="0.0" value="0.536"/>\n');
        fprintf(f, '\t\t</relativeRate>\n');
        fprintf(f, '\t\t<gammaShape gammaCategories="4">\n');
        fprintf(f, '\t\t\t<parameter id="cds.CP1.alpha" lower="0.0" value="0.169"/>\n');
        fprintf(f, '\t\t</gammaShape>\n');
        fprintf(f, '\t</siteModel>\n');
        fprintf(f, '\n');
        fprintf(f, '\t<!-- site model                                                              -->\n');
        fprintf(f, '\t<siteModel id="cds.CP2.siteModel">\n');
        fprintf(f, '\t\t<substitutionModel>\n');
        fprintf(f, '\t\t\t<HKYModel idref="cds.CP2.hky"/>\n');
        fprintf(f, '\t\t</substitutionModel>\n');
        fprintf(f, '\t\t<relativeRate>\n');
        fprintf(f, '\t\t\t<parameter id="cds.CP2.mu" lower="0.0" value="0.503"/>\n');
        fprintf(f, '\t\t</relativeRate>\n');
        fprintf(f, '\t\t<gammaShape gammaCategories="4">\n');
        fprintf(f, '\t\t\t<parameter id="cds.CP2.alpha" lower="0.0" value="5.424E-2"/>\n');
        fprintf(f, '\t\t</gammaShape>\n');
        fprintf(f, '\t</siteModel>\n');
        fprintf(f, '\n');
        fprintf(f, '\t<!-- site model                                                              -->\n');
        fprintf(f, '\t<siteModel id="cds.CP3.siteModel">\n');
        fprintf(f, '\t\t<substitutionModel>\n');
        fprintf(f, '\t\t\t<HKYModel idref="cds.CP3.hky"/>\n');
        fprintf(f, '\t\t</substitutionModel>\n');
        fprintf(f, '\t\t<relativeRate>\n');
        fprintf(f, '\t\t\t<parameter id="cds.CP3.mu" lower="0.0" value="1.412"/>\n');
        fprintf(f, '\t\t</relativeRate>\n');
        fprintf(f, '\t\t<gammaShape gammaCategories="4">\n');
        fprintf(f, '\t\t\t<parameter id="cds.CP3.alpha" lower="0.0" value="0.605"/>\n');
        fprintf(f, '\t\t</gammaShape>\n');
        fprintf(f, '\t</siteModel>\n');
        fprintf(f, '\n');    
        fprintf(f, '\t<siteModel id="ig.siteModel">\n');
        fprintf(f, '\t\t<substitutionModel>\n');
        fprintf(f, '\t\t\t<HKYModel idref="ig.hky"/>\n');
        fprintf(f, '\t\t</substitutionModel>\n');
        fprintf(f, '\t\t<relativeRate>\n');
        fprintf(f, '\t\t\t<parameter id="ig.mu" lower="0.0" value="1.597"/>\n');
        fprintf(f, '\t\t</relativeRate>\n');
        fprintf(f, '\t\t<gammaShape gammaCategories="4">\n');
        fprintf(f, '\t\t\t<parameter id="ig.alpha" lower="0.0" value="0.354"/>\n');
        fprintf(f, '\t\t</gammaShape>\n');
        fprintf(f, '\t</siteModel>\n');
        fprintf(f, '\n');               
        fprintf(f, '\t<compoundParameter id="allMus">\n');
        fprintf(f, '\t\t<parameter idref="cds.CP1.mu"/>\n');
        fprintf(f, '\t\t<parameter idref="cds.CP2.mu"/>\n');
        fprintf(f, '\t\t<parameter idref="cds.CP3.mu"/>\n');
        fprintf(f, '\t\t<parameter idref="ig.mu"/>\n');
        fprintf(f, '\t</compoundParameter>\n');
        fprintf(f, '\n');    
        fprintf(f, '\t<!-- Likelihood for tree given sequence data                                 -->\n');
        fprintf(f, '\t<treeLikelihood id="cds.CP1.treeLikelihood" useAmbiguities="false">\n');
        fprintf(f, '\t\t<patterns idref="cds.CP1.patterns"/>\n');
        fprintf(f, '\t\t<treeModel idref="treeModel"/>\n');
        fprintf(f, '\t\t<siteModel idref="cds.CP1.siteModel"/>\n');
        fprintf(f, '\t\t<strictClockBranchRates idref="branchRates"/>\n');
        fprintf(f, '\t</treeLikelihood>\n');
        fprintf(f, '\n');
        fprintf(f, '\t<!-- Likelihood for tree given sequence data                                 -->\n');
        fprintf(f, '\t<treeLikelihood id="cds.CP2.treeLikelihood" useAmbiguities="false">\n');
        fprintf(f, '\t\t<patterns idref="cds.CP2.patterns"/>\n');
        fprintf(f, '\t\t<treeModel idref="treeModel"/>\n');
        fprintf(f, '\t\t<siteModel idref="cds.CP2.siteModel"/>\n');
        fprintf(f, '\t\t<strictClockBranchRates idref="branchRates"/>\n');
        fprintf(f, '\t</treeLikelihood>\n');
        fprintf(f, '\n');
        fprintf(f, '\t<!-- Likelihood for tree given sequence data                                 -->\n');
        fprintf(f, '\t<treeLikelihood id="cds.CP3.treeLikelihood" useAmbiguities="false">\n');
        fprintf(f, '\t\t<patterns idref="cds.CP3.patterns"/>\n');
        fprintf(f, '\t\t<treeModel idref="treeModel"/>\n');
        fprintf(f, '\t\t<siteModel idref="cds.CP3.siteModel"/>\n');
        fprintf(f, '\t\t<strictClockBranchRates idref="branchRates"/>\n');
        fprintf(f, '\t</treeLikelihood>\n');
        fprintf(f, '\n');
        fprintf(f, '\t<!-- Likelihood for tree given sequence data                                 -->\n');
        fprintf(f, '\t<treeLikelihood id="ig.treeLikelihood" useAmbiguities="false">\n');
        fprintf(f, '\t\t<patterns idref="ig.patterns"/>\n');
        fprintf(f, '\t\t<treeModel idref="treeModel"/>\n');
        fprintf(f, '\t\t<siteModel idref="ig.siteModel"/>\n');
        fprintf(f, '\t\t<strictClockBranchRates idref="branchRates"/>\n');
        fprintf(f, '\t</treeLikelihood>\n');
        fprintf(f, '\n');
        fprintf(f, '\n');
        fprintf(f, '\t<generalDataType id="location.dataType">\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t<!-- SLE -->\n');
        fprintf(f, '\t\t<state code="Kailahun"/>\n');fprintf(f, '\t\t\t\t<!-- 43 -->\n');
        fprintf(f, '\t\t<state code="Kenema"/>\n');
        fprintf(f, '\t\t<state code="Kono"/>\n');
        fprintf(f, '\t\t<state code="Bombali"/>\n');
        fprintf(f, '\t\t<state code="Kambia"/>\n');
        fprintf(f, '\t\t<state code="Koinadugu"/>\n');
        fprintf(f, '\t\t<state code="PortLoko"/>\n');
        fprintf(f, '\t\t<state code="Tonkolili"/>\n');
        fprintf(f, '\t\t<state code="Bo"/>\n');
        fprintf(f, '\t\t<state code="Bonthe"/>\n');
        fprintf(f, '\t\t<state code="Moyamba"/>\n');
        fprintf(f, '\t\t<state code="Pujehun"/>\n');
        fprintf(f, '\t\t<state code="WesternRural"/>\n');
        fprintf(f, '\t\t<state code="WesternUrban"/>\n');
        fprintf(f, '\t</generalDataType>\n');
        fprintf(f, '\n');    
        fprintf(f, '\t<attributePatterns id="location.pattern" attribute="location">\n');
        fprintf(f, '\t\t<taxa idref="taxa"/>\n');
        fprintf(f, '\t\t<generalDataType idref="location.dataType"/>\n');
        fprintf(f, '\t</attributePatterns>\n');
        fprintf(f, '\n');
        fprintf(f, '\t<glmSubstitutionModel id="location.model">\n');
        fprintf(f, '\t\t<dataType idref="location.dataType"/>\n');
        fprintf(f, '\t\t<rootFrequencies>\n');
        fprintf(f, '\t\t\t<frequencyModel id="location.frequencyModel" normalize="true">\n');
        fprintf(f, '\t\t\t<dataType idref="location.dataType"/>\n');
        fprintf(f, '\t\t\t\t<frequencies>\n');
        fprintf(f, '\t\t\t\t\t<parameter id="location.frequencies" dimension="%d"/>\n', length(locations));
        fprintf(f, '\t\t\t\t</frequencies>\n');
        fprintf(f, '\t\t\t</frequencyModel>\n');
        fprintf(f, '\t\t</rootFrequencies>\n');
        cov_names = fieldnames(m_cov);
        fprintf(f, '\t\t<glmModel id="location.glmModel" checkFullRank="false" family="logLinear" checkIdentifiability="false">\n');
        fprintf(f, '\t\t\t<independentVariables>\n');
        fprintf(f, '\t\t\t\t<parameter id="location.glmCoefficients" dimension="%d" value="0.1"/>\n', length(cov_names));
        fprintf(f, '\t\t\t\t<indicator>\n');
        fprintf(f, '\t\t\t\t\t<parameter id="location.coefIndicator" dimension="%d" value="1"/>\n', length(cov_names)); 
        fprintf(f, '\t\t\t\t</indicator>\n');
        fprintf(f, '\t\t\t\t<designMatrix id="location.designMatrix">\n');        
        for k = 1 : length(cov_names)
             fprintf(f, '\t\t\t\t\t<parameter id="%s" value="%s"/>\n', cov_names{k}, sprintf('%f ',m_cov.(cov_names{k})));
        end
        fprintf(f, '\t\t\t\t</designMatrix>\n');
        fprintf(f, '\t\t\t</independentVariables>\n');
        fprintf(f, '\t\t</glmModel>\n');
        fprintf(f, '\t</glmSubstitutionModel>\n');
        fprintf(f, '\n');  
        fprintf(f, '\t<sumStatistic id="location.nonZeroIndicators" name="nonZeroIndicatorCount" elementwise="true">\n');
        fprintf(f, '\t\t<parameter idref="location.coefIndicator"/>\n');
        fprintf(f, '\t</sumStatistic>\n');
        fprintf(f, '\t<productStatistic id="location.coefficientsTimesIndicators" elementwise="false">\n');
        fprintf(f, '\t\t<parameter idref="location.glmCoefficients"/>\n');
        fprintf(f, '\t\t<parameter idref="location.coefIndicator"/>\n');
        fprintf(f, '\t</productStatistic>\n');
        fprintf(f, '\n');
        
        fprintf(f, '\t<siteModel id="location.siteModel">\n');
        fprintf(f, '\t\t<substitutionModel>\n');
        fprintf(f, '\t\t\t<generalSubstitutionModel idref="location.model"/>\n');
        fprintf(f, '\t\t</substitutionModel>\n');
        fprintf(f, '\t</siteModel>\n');

        
        fprintf(f, '\t<strictClockBranchRates id="location.branchRates">\n');
        fprintf(f, '\t\t<rate>\n');
        fprintf(f, '\t\t\t<parameter id="location.rate" value="1.0" lower="0.0"/>\n');
        fprintf(f, '\t\t</rate>\n');
        fprintf(f, '\t</strictClockBranchRates>\n');
        
        fprintf(f, '\t<ancestralTreeLikelihood id="location.treeLikelihood" stateTagName="location.states" useUniformization="true" saveCompleteHistory="false" logCompleteHistory="false">\n');
        fprintf(f, '\t\t<attributePatterns idref="location.pattern"/>\n');
        fprintf(f, '\t\t<treeModel idref="treeModel"/>\n');
        fprintf(f, '\t\t<siteModel idref="location.siteModel"/>\n');
        fprintf(f, '\t\t<glmSubstitutionModel idref="location.model"/>\n');
        fprintf(f, '\t\t<frequencyModel id="location.root.frequencyModel" normalize="true">\n');
        fprintf(f, '\t\t\t<generalDataType idref="location.dataType"/>\n');
        fprintf(f, '\t\t\t<frequencies>\n');
        fprintf(f, '\t\t\t\t<parameter id="location.root.frequencies" dimension="%d"/>\n', 14);
        fprintf(f, '\t\t\t</frequencies>\n');
        fprintf(f, '\t\t</frequencyModel>\n');
        fprintf(f, '\t\t<strictClockBranchRates idref="location.branchRates"/>\n');
        fprintf(f, '\t</ancestralTreeLikelihood>\n');

        
        
%         fprintf(f, '\t<markovJumpsTreeLikelihood id="location.treeLikelihood" useAmbiguities="true" stateTagName="location.states" useUniformization="true" numberOfSimulants="1" saveCompleteHistory="true" logCompleteHistory="true">\n');
%         fprintf(f, '\t\t<attributePatterns idref="location.pattern"/>\n');
%         fprintf(f, '\t\t<treeModel idref="treeModel"/>\n');
%         fprintf(f, '\t\t<siteModel idref="location.siteModel"/>\n');
%         fprintf(f, '\t\t<glmSubstitutionModel idref="location.model"/>\n');
%         fprintf(f, '\t\t<strictClockBranchRates idref="location.branchRates"/>\n');
%         fprintf(f, '\t\t<frequencyModel id="location.root.frequencyModel" normalize="true">\n');
%         fprintf(f, '\t\t\t<generalDataType idref="location.dataType"/>\n');
%         fprintf(f, '\t\t\t<frequencies>\n');
%         fprintf(f, '\t\t\t\t<parameter id="location.root.frequencies" dimension="14"/>\n');
%         fprintf(f, '\t\t\t</frequencies>\n');
%         fprintf(f, '\t\t</frequencyModel>\n');
% 
%         fprintf(f, '\t\t<rewards>\n');
%         fprintf(f, '\t\t\t<parameter id="Kailahun_R"			value="1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="Kenema_R"			value="0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="Kono_R"				value="0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="Bombali_R"			value="0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="Kambia_R"			value="0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="Koinadugu_R"			value="0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="PortLoko_R"			value="0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="Tonkolili_R"			value="0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="Bo_R"				value="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="Bonthe_R"			value="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="Moyamba_R"			value="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="Pujehun_R"			value="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="WesternRural_R"		value="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="WesternUrban_R"		value="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0" />\n');
%         fprintf(f, '\t\t</rewards>\n');
%         fprintf(f, '\t</markovJumpsTreeLikelihood>\n');

        
        fprintf(f, '\t<sumStatistic id="location.nonZeroRates" elementwise="true">\n');
		fprintf(f, '\t\t<parameter idref="location.coefIndicator"/>\n');
        fprintf(f, '\t</sumStatistic>\n');
        
        
        fprintf(f, '\n');    
        fprintf(f, '\t<operators id="operators" optimizationSchedule="log">\n');
        fprintf(f, '\t\t<subtreeLeap size="1.0" weight="100">\n');
        fprintf(f, '\t\t\t<treeModel idref="treeModel"/>\n');
        fprintf(f, '\t\t</subtreeLeap>\n');
        fprintf(f, '\t\t<upDownOperator scaleFactor="0.75" weight="5">\n');
        fprintf(f, '\t\t\t<up>\n');
        fprintf(f, '\t\t\t\t<parameter idref="Ebola.clock.rate"/>\n');
        fprintf(f, '\t\t\t</up>\n');
        fprintf(f, '\t\t\t<down>\n');
        fprintf(f, '\t\t\t\t<parameter idref="treeModel.allInternalNodeHeights"/>\n');
        fprintf(f, '\t\t\t</down>\n');
        fprintf(f, '\t\t</upDownOperator>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t<gmrfGridBlockUpdateOperator scaleFactor="1.0" weight="100">\n');
        fprintf(f, '\t\t\t<gmrfSkyrideLikelihood idref="skygrid"/>\n');
        fprintf(f, '\t\t</gmrfGridBlockUpdateOperator>\n');
        fprintf(f, '\t\t<scaleOperator scaleFactor="0.75" weight="5">\n');
        fprintf(f, '\t\t\t<parameter idref="skygrid.precision"/>\n');
        fprintf(f, '\t\t</scaleOperator>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t<scaleOperator scaleFactor="0.75" weight="3">\n');
        fprintf(f, '\t\t\t<parameter idref="location.rate"/>\n');
        fprintf(f, '\t\t</scaleOperator>\n');
        fprintf(f, '\t\t<deltaExchange delta="0.75" weight="1">\n');
        fprintf(f, '\t\t\t<parameter idref="location.root.frequencies"/>\n');
        fprintf(f, '\t\t</deltaExchange>\n');
        fprintf(f, '\t\t<bitFlipOperator weight="6" usesPriorOnSum="true">\n');
        fprintf(f, '\t\t\t<parameter idref="location.coefIndicator"/>\n');
        fprintf(f, '\t\t</bitFlipOperator>\n');
        fprintf(f, '\t\t<rateBitExchangeOperator weight="3" usesPriorOnSum="true">\n');
        fprintf(f, '\t\t\t<bits>\n');
        fprintf(f, '\t\t\t\t<parameter idref="location.coefIndicator"/>\n');
        fprintf(f, '\t\t\t</bits>\n');
        fprintf(f, '\t\t\t<rates>\n');
        fprintf(f, '\t\t\t\t<parameter idref="location.glmCoefficients"/>\n');
        fprintf(f, '\t\t\t</rates>\n');
        fprintf(f, '\t\t</rateBitExchangeOperator>\n');
        fprintf(f, '\t\t<bitMoveOperator weight="3.0" numBitsToMove="1" usesPriorOnSum="true">\n');
        fprintf(f, '\t\t\t<bits>\n');
        fprintf(f, '\t\t\t\t<parameter idref="location.coefIndicator"/>\n');
        fprintf(f, '\t\t\t</bits>\n');
        fprintf(f, '\t\t</bitMoveOperator>\n');
        fprintf(f, '\t\t<randomWalkOperator windowSize="0.5" weight="20">\n');
        fprintf(f, '\t\t\t<parameter idref="location.glmCoefficients"/>\n');
        fprintf(f, '\t\t</randomWalkOperator>\n');
%         fprintf(f, '\t\t<mvnOperator scaleFactor="1.0" weight="10" formXtXInverse="true">\n');
%         fprintf(f, '\t\t\t<parameter idref="location.glmCoefficients"/>\n');
%         fprintf(f, '\t\t\t<varMatrix>\n');
%         fprintf(f, '\t\t\t\t<parameter idref="location.designMatrix"/>\n');
%         fprintf(f, '\t\t\t</varMatrix>\n');
%         fprintf(f, '\t\t</mvnOperator>\n');
        fprintf(f, '\n'); 
%         fprintf(f, '\t\t<randomWalkOperator windowSize="0.5" weight="20">\n');
%         fprintf(f, '\t\t\t<parameter idref="glmRandCoefficients"/>\n');
%         fprintf(f, '\t\t</randomWalkOperator>\n');
        fprintf(f, '\t\t</operators>\n');
        fprintf(f, '\n');        
        fprintf(f, '\t<mcmc id="mcmc" chainLength="100000000" autoOptimize="true" operatorAnalysis="Makona_1610_cds_ig.GLM.ops">\n');
        fprintf(f, '\t\t<posterior id="posterior">\n');
        fprintf(f, '\t\t\t<prior id="prior">\n');
        fprintf(f, '\t\t\t\t<logNormalPrior mean="1.0" stdev="1.25" offset="0.0" meanInRealSpace="false">\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="cds.CP1.kappa"/>\n');
        fprintf(f, '\t\t\t\t</logNormalPrior>\n');
        fprintf(f, '\t\t\t\t<logNormalPrior mean="1.0" stdev="1.25" offset="0.0" meanInRealSpace="false">\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="cds.CP2.kappa"/>\n');
        fprintf(f, '\t\t\t\t</logNormalPrior>\n');
        fprintf(f, '\t\t\t\t<logNormalPrior mean="1.0" stdev="1.25" offset="0.0" meanInRealSpace="false">\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="cds.CP3.kappa"/>\n');
        fprintf(f, '\t\t\t\t</logNormalPrior>\n');
        fprintf(f, '\t\t\t\t<exponentialPrior mean="0.5" offset="0.0">\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="cds.CP1.alpha"/>\n');
        fprintf(f, '\t\t\t\t</exponentialPrior>\n');
        fprintf(f, '\t\t\t\t<exponentialPrior mean="0.5" offset="0.0">\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="cds.CP2.alpha"/>\n');
        fprintf(f, '\t\t\t\t</exponentialPrior>\n');
        fprintf(f, '\t\t\t\t<exponentialPrior mean="0.5" offset="0.0">\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="cds.CP3.alpha"/>\n');
        fprintf(f, '\t\t\t\t</exponentialPrior>\n');
        fprintf(f, '\t\t\t\t<logNormalPrior mean="1.0" stdev="1.25" offset="0.0" meanInRealSpace="false">\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="ig.kappa"/>\n');
        fprintf(f, '\t\t\t\t</logNormalPrior>\n');
        fprintf(f, '\t\t\t\t<exponentialPrior mean="0.5" offset="0.0">\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="ig.alpha"/>\n');
        fprintf(f, '\t\t\t\t</exponentialPrior>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t\t<gmrfSkyGridLikelihood idref="skygrid"/>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t\t<oneOnXPrior>\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="location.rate"/>\n');
        fprintf(f, '\t\t\t\t</oneOnXPrior>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t\t<normalPrior mean="0" stdev="2">\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="location.glmCoefficients"/>\n');
        fprintf(f, '\t\t\t\t</normalPrior>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t\t<poissonPrior mean="5.0" offset="0.0">\n');
        fprintf(f, '\t\t\t\t\t<statistic idref="location.nonZeroRates"/>\n');
		fprintf(f, '\t\t\t\t</poissonPrior>\n');

%         fprintf(f, '\t\t\t\t<binomialLikelihood>\n');
%         fprintf(f, '\t\t\t\t\t<proportion>\n');
%         fprintf(f, '\t\t\t\t\t\t<parameter value="0.0273"/>\n');
%         fprintf(f, '\t\t\t\t\t</proportion>\n');
%         fprintf(f, '\t\t\t\t\t<trials>\n');
%         fprintf(f, '\t\t\t\t\t\t<parameter dimension="14" value="1"/>\n'); 
%         fprintf(f, '\t\t\t\t\t</trials>\n');
%         fprintf(f, '\t\t\t\t\t<counts>\n');
%         fprintf(f, '\t\t\t\t\t\t<parameter idref="location.coefIndicator"/>\n');
%         fprintf(f, '\t\t\t\t\t</counts>\n');
%         fprintf(f, '\t\t\t\t</binomialLikelihood>\n');
        fprintf(f, '\n');
%         fprintf(f, '\t\t\t\t<distributionLikelihood idref="randomEffects.prior"/>\n');
%         fprintf(f, '\t\t\t\t<gammaPrior idref="gammaPrecisionRandEffectPrior"/>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t\t<!-- START Discrete Traits Model                                             -->\n');
        fprintf(f, '\t\t\t\t<glmSubstitutionModel idref="location.model"/>\n');
        fprintf(f, '\t\t\t\t<!-- END Discrete Traits Model                                               -->\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t</prior>\n');
        fprintf(f, '\t\t\t<likelihood id="likelihood">\n');
        fprintf(f, '\t\t\t\t<treeLikelihood idref="cds.CP1.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t\t<treeLikelihood idref="cds.CP2.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t\t<treeLikelihood idref="cds.CP3.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t\t<treeLikelihood idref="ig.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t\t<strictClockBranchRate idref="branchRates"/>\n');
        fprintf(f, '\t\t\t\t<!-- START Discrete Traits Model                                             -->\n');
        fprintf(f, '\t\t\t\t<markovJumpsTreeLikelihood idref="location.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t\t<!-- END Discrete Traits Model                                               -->\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t</likelihood>\n');
        fprintf(f, '\t\t</posterior>\n');
        fprintf(f, '\t\t<operators idref="operators"/>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t<!-- write log to screen                                                     -->\n');
        fprintf(f, '\t\t<log id="screenLog" logEvery="10000">\n');
        fprintf(f, '\t\t\t<column label="Posterior" dp="4" width="12">\n');
        fprintf(f, '\t\t\t\t<posterior idref="posterior"/>\n');
        fprintf(f, '\t\t\t</column>\n');
        fprintf(f, '\t\t\t<column label="Prior" dp="4" width="12">\n');
        fprintf(f, '\t\t\t\t<prior idref="prior"/>\n');
        fprintf(f, '\t\t\t</column>\n');
        fprintf(f, '\t\t\t<column label="Likelihood" dp="4" width="12">\n');
        fprintf(f, '\t\t\t\t<likelihood idref="likelihood"/>\n');
        fprintf(f, '\t\t\t</column>\n');
        fprintf(f, '\t\t\t<column label="rootHeight" sf="6" width="12">\n');
        fprintf(f, '\t\t\t\t<parameter idref="treeModel.rootHeight"/>\n');
        fprintf(f, '\t\t\t</column>\n');
        fprintf(f, '\t\t\t<column label="location.nonZeroPredictors" sf="6" width="12">\n');
        fprintf(f, '\t\t\t\t<sumStatistic idref="location.nonZeroIndicators"/>\n');
        fprintf(f, '\t\t\t</column>\n');
        fprintf(f, '\t\t</log>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t<!-- write log to file                                                       -->\n');
        fprintf(f, '\t\t<log id="fileLog1" logEvery="500000" fileName="%s.GLM.log" overwrite="false">\n',flog);
        fprintf(f, '\t\t\t<posterior idref="posterior"/>\n');
        fprintf(f, '\t\t\t<prior idref="prior"/>\n');
        fprintf(f, '\t\t\t<likelihood idref="likelihood"/>\n');
        fprintf(f, '\t\t\t<parameter idref="treeModel.rootHeight"/>\n');
        fprintf(f, '\t\t\t<parameter idref="cds.CP1.kappa"/>\n');
        fprintf(f, '\t\t\t<parameter idref="cds.CP2.kappa"/>\n');
        fprintf(f, '\t\t\t<parameter idref="cds.CP3.kappa"/>\n');
        fprintf(f, '\t\t\t<parameter idref="ig.kappa"/>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t<parameter idref="cds.CP1.alpha"/>\n');
        fprintf(f, '\t\t\t<parameter idref="cds.CP2.alpha"/>\n');
        fprintf(f, '\t\t\t<parameter idref="cds.CP3.alpha"/>\n');
        fprintf(f, '\t\t\t<parameter idref="ig.alpha"/>\n');
        fprintf(f, '\t\t\t<compoundParameter idref="allMus"/>\n');
        fprintf(f, '\t\t\t<parameter idref="skygrid.precision"/>\n');
        fprintf(f, '\t\t\t<parameter idref="skygrid.logPopSize"/>\n');
        fprintf(f, '\t\t\t<parameter idref="skygrid.cutOff"/>\n');
        fprintf(f, '\t\t\t<treeLikelihood idref="cds.CP1.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t<treeLikelihood idref="cds.CP2.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t<treeLikelihood idref="cds.CP3.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t<treeLikelihood idref="ig.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t<gmrfSkyGridLikelihood idref="skygrid"/>\n');
        fprintf(f, '\n');
%         fprintf(f, '\t\t\t<distributionLikelihood idref="randomEffects.prior"/>\n');
%         fprintf(f, '\t\t\t<gammaPrior idref="gammaPrecisionRandEffectPrior"/>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t<parameter idref="location.rate"/>\n');
        fprintf(f, '\t\t\t<parameter idref="location.root.frequencies"/>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t<!-- START Discrete Traits Model                                             -->\n');
        fprintf(f, '\t\t\t<sumStatistic idref="location.nonZeroIndicators"/>\n');
        fprintf(f, '\t\t\t<parameter idref="location.coefIndicator"/>\n');
        fprintf(f, '\t\t\t<parameter idref="location.glmCoefficients"/>\n');
        fprintf(f, '\t\t\t<productStatistic idref="location.coefficientsTimesIndicators"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="LNprecRandEffect"/>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t<markovJumpsTreeLikelihood idref="location.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t<glmSubstitutionModel idref="location.model"/>\n');
        fprintf(f, '\t\t\t<!-- END Discrete Traits Model                                               -->\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t</log>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t<log id="Makona_1610_cds_ig.GLM.LocationLog" logEvery="500000" fileName="%s.GLM.location.log">\n',flog);
        fprintf(f, '\t\t\t<parameter idref="location.coefIndicator"/>\n');
        fprintf(f, '\t\t\t<parameter idref="location.glmCoefficients"/>\n');
        fprintf(f, '\t\t\t<productStatistic idref="location.coefficientsTimesIndicators"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="coefRandIndicator"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="glmRandCoefficients"/>\n');
        fprintf(f, '\t\t\t<glmModel idref="location.glmModel"/>\n');
        fprintf(f, '\t\t\t<sumStatistic idref="location.nonZeroIndicators"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="LNprecRandEffect"/>\n');
        fprintf(f, '\t\t</log>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t<!-- END Discrete Traits Model                                               -->\n');
        fprintf(f, '\n');
%         fprintf(f, '\t\t<log id="completeJumpHistory" logEvery="500000" fileName="%s.GLM.history.log">\n',flog);
%         fprintf(f, '\t\t\t<completeHistoryLogger>\n');
%         fprintf(f, '\t\t\t\t<markovJumpsTreeLikelihood idref="location.treeLikelihood"/>\n');
%         fprintf(f, '\t\t\t</completeHistoryLogger>\n');
%         fprintf(f, '\t\t</log>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t<!-- write tree log to file                                                  -->\n');
        fprintf(f, '\t\t<logTree id="treeFileLog" logEvery="500000" nexusFormat="true" fileName="%s.trees" sortTranslationTable="true">\n',flog);
        fprintf(f, '\t\t\t<treeModel idref="treeModel"/>\n');
        fprintf(f, '\t\t\t<trait name="rate" tag="rate">\n');
        fprintf(f, '\t\t\t\t<strictClockBranchRate idref="branchRates"/>\n');
        fprintf(f, '\t\t\t</trait>\n');
        fprintf(f, '\t\t\t<markovJumpsTreeLikelihood idref="location.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t<posterior idref="posterior"/>\n');
        fprintf(f, '\t\t</logTree>\n');
        fprintf(f, '\t</mcmc>\n');
        fprintf(f, '\n');
        fprintf(f, '\t<report>\n');
        fprintf(f, '\t\t<property name="timer">\n');
        fprintf(f, '\t\t\t<mcmc idref="mcmc"/>\n');
        fprintf(f, '\t\t</property>\n');
        fprintf(f, '\t</report>\n');
        fprintf(f, '</beast>\n');
        
        fclose(f);
    end
end