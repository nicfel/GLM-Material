%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates structured coalescent xmls from the master trees. Always creates
% 3 xmls per tree with different initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% create lisco files
tree_files = dir('master/*.tree');

system('rm -r xmls');
system('mkdir xmls');

states = 5;
% read in covariates
f = fopen('mCovariates.txt');
while ~feof(f)
    line = fgets(f);
    tmp = strsplit(line, '\t');
    m_cov{str2num(tmp{1})} = str2num(tmp{2});
    % use the forward covariates
    for i = 1 : size(m_cov{str2num(tmp{1})},2)
        count_a = 1;
        count_b = 1;
        for j = 1 : size(m_cov{str2num(tmp{1})},1)/(states*(states-1))
            m = zeros(states, states);
            for a = 1 : states
                for b = 1 : states
                    if a~=b
                        m(a,b) = m_cov{str2num(tmp{1})}(count_a, i);
                        count_a = count_a+1;
                    end
                end
            end
            m = m';
            %
            for a = 1 : states
                for b = 1 : states
                    if a~=b
                        m_cov{str2num(tmp{1})}(count_b, i) = m(a,b);
                        count_b = count_b+1;
                    end
                end
            end
        end

    end
end
fclose(f);

f = fopen('NeCovariates.txt');
while ~feof(f)
    line = fgets(f);
    tmp = strsplit(line, '\t');
    Ne_cov{str2num(tmp{1})} = str2num(tmp{2});
end
fclose(f);

%%

system('rm -r xmls');
system('mkdir xmls');

for i = 1 : length(tree_files)
    disp(tree_files(i).name)
    % read tree files
    g = fopen(['master/' tree_files(i).name],'r');
    t = textscan(g,'%s'); fclose(g);
    % read in the sequences
    fasta = fastaread(['Sequences/' strrep(tree_files(i).name, '.tree', '.fasta')]);
    seq_header = cell(0,0);
    for j = 1 : length(fasta)
        seq_header{j} = strrep(fasta(j).Header, 'loc', 'loc_');
    end
    
    % coalescing
    tree_tmp2 = regexprep(t{1}(end-1),'&type="L",location="(\d*)",reaction="Coalescence",time=(\d*).(\d*)','');

    %migrating
    tree_tmp1 = regexprep(tree_tmp2,'&type="L",location="(\d*)",reaction="Migration",time=(\d*).(\d*)','');

    
    % make the MASTER tree compatible with BEAST2
    % sampling
    tree_tmp1 = regexprep(tree_tmp1,'E[-](\d)]',']');
    tip_locs = regexp(tree_tmp1,'[&type="L",location="(\d*)",time=(\d*)\.(\d*)\]','match');
     
    for j = 1 : length(tip_locs{1})
        loc = regexprep(tip_locs{1}{j},'[&type="L",location="(\d*)",time=(\d*)\.(\d*)\]','$1');
        tree_tmp1 = strrep(tree_tmp1,tip_locs{1}{j},['loc_' loc 'kickout']);
        tree_tmp1 = strrep(tree_tmp1,'kickout','');
    end

    tree_tmp = regexprep(tree_tmp1,'(\d*)loc_','inv$1loc_');
    
    tree = strrep(tree_tmp{1},'[]','');
    if ~isempty(strfind(tree,']'))
        b = strfind(tree,']');
        c = tree((b-50):(b+50));
        disp(tree_files(i).name)
    end

    % get the leafnames
    ptree = phytreeread(tree);
    leafnames = get(ptree,'leafnames');
    
    print_tree = tree;
    
    % get the covariates
    
    % make tripletts of all runs with different random initial values
    for tr = 1 : 1    
        % make the xmls for the structcoal
        flog = strrep(tree_files(i).name,'master.tree',sprintf('%dmascot',tr));
        fname = sprintf('xmls/%s.xml',flog);
        f = fopen(fname,'w');


        fprintf(g,'<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate=''Standard'' beautistatus='''' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.mascot.dynamics:beast.mascot.distribution:beast.mascot.logger" version="2.0">\n');

        fprintf(g,'\t<data id="sequences" name="alignment">\n');
        for j = 1 : length(leafnames)
            seq_ind = find(ismember(seq_header, leafnames{j}));
            fprintf(g,'\t\t<sequence id="%s" taxon="%s" totalcount="4" value="%s"/>\n',leafnames{j},leafnames{j},fasta(seq_ind).Sequence);
        end
        fprintf(g,'\t</data>\n');
        fprintf(g,'<map name="prior">beast.math.distributions.Prior</map>\n');
        fprintf(g,'<map name="maxRate">beast.mascot.glmmodel.MaxRate</map>\n\n');

        fprintf(g,'\t<run id="mcmc" spec="MCMC" chainLength="1000000">\n');
        fprintf(g,'\t\t<state id="state" storeEvery="5000">\n');
        fprintf(g,'\t\t\t<stateNode id="tree" spec="beast.app.mascot.beauti.TreeWithTrait">\n');
        fprintf(g,'\t\t\t\t<typeTrait id="typeTraitSet.t" spec="beast.evolution.tree.TraitSet" traitname="type" value="');
        for j = 1 : length(leafnames)-1
            tmp1 = strsplit(leafnames{j},'_');
            fprintf(g,'%s=state%s,',leafnames{j},tmp1{end});
        end
        tmp1 = strsplit(leafnames{end},'_');
        fprintf(g,'%s=state%s">\n',leafnames{end},tmp1{end});

        fprintf(g,'\t\t\t\t\t<taxa id="TaxonSet.0" spec="TaxonSet">\n');
        fprintf(g,'\t\t\t\t\t\t<alignment idref="sequences"/>\n');
        fprintf(g,'\t\t\t\t\t</taxa>\n');
        fprintf(g,'\t\t\t\t</typeTrait>\n');
        fprintf(g,'\t\t\t</stateNode>\n');        
        fprintf(g,'\t\t\t<parameter id="migrationScaler" name="stateNode" dimension="0">0</parameter>\n');
        fprintf(g,'\t\t\t<stateNode id="migrationIndicator" spec="parameter.BooleanParameter" dimension="0">true</stateNode>\n');
        fprintf(g,'\t\t\t<parameter id="migrationClock" name="stateNode" dimension="0">0.01</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="NeScaler" name="stateNode" dimension="0">0</parameter>\n');
        fprintf(g,'\t\t\t<stateNode id="NeIndicator" spec="parameter.BooleanParameter" dimension="0">true</stateNode>\n');
        fprintf(g,'\t\t\t<parameter id="NeClock" name="stateNode" dimension="0">10</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="clockrate" name="stateNode" dimension="1">0.001</parameter>\n');
        fprintf(g,'\t\t</state>\n');
        

                
        fprintf(g,'\t\t<init spec="beast.util.TreeParser" id="NewickTree.t:Species" adjustTipHeights="false"\n');
        fprintf(g,'\t\t\tinitial="@tree" taxa="@sequences"\n');
        fprintf(g,'\t\t\tIsLabelledNewick="true"\n');
        fprintf(g,'\t\t\tnewick="%s"/>\n',print_tree);
        fprintf(g,' \t\t<distribution id="posterior" spec="util.CompoundDistribution">\n');
        fprintf(g,'\t\t\t<distribution id="prior" spec="util.CompoundDistribution">\n');
        
        % define params of the lognormal distribution of the Ne
        mean_ne = 10;
        sigma_ne = 0.5;
        mu_ne = log(mean_ne) - sigma_ne^2/2;

        % define params of the lognormal distribution of the migration rate
        mean_rea = 0.1;
        sigma_rea = 0.5;
        mu_rea = log(mean_rea) - sigma_rea^2/2;

        
        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@NeClock">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.LogNormalDistributionModel" M="%.5f" S="0.5"/>\n', mu_ne);
        fprintf(g,'\t\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@migrationClock">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.LogNormalDistributionModel" M="%.5f" S="0.5"/>\n', mu_rea);
        fprintf(g,'\t\t\t\t</distribution>\n');
                
        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@migrationScaler">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal" mean="0" sigma="1"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@NeScaler">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal" mean="0" sigma="1"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        
        fprintf(g,'\t\t\t\t<prior id="nonZeroRatePrior.Migration" name="distribution">\n');
        fprintf(g,'\t\t\t\t\t<x id="nonZeroCovarates.Migration" spec="util.Sum">\n');
        fprintf(g,'\t\t\t\t\t\t<arg idref="migrationIndicator"/>\n');
        fprintf(g,'\t\t\t\t\t</x>\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Poisson" lambda=''0.693''/>\n');
        fprintf(g,'\t\t\t\t</prior>\n');
        
        fprintf(g,'\t\t\t\t<prior id="nonZeroRatePrior.Ne" name="distribution">\n');
        fprintf(g,'\t\t\t\t\t<x id="nonZeroCovarates.Ne" spec="util.Sum">\n');
        fprintf(g,'\t\t\t\t\t\t<arg idref="NeIndicator"/>\n');
        fprintf(g,'\t\t\t\t\t</x>\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Poisson" lambda=''0.693''/>\n');
        fprintf(g,'\t\t\t\t</prior>\n');


        fprintf(g,'\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t<distribution id="likelihood" spec="util.CompoundDistribution">\n');
        fprintf(g,'\t\t\t\t<distribution id="coalescent" spec="Mascot">\n');
        fprintf(g,'\t\t\t\t\t<structuredTreeIntervals spec="StructuredTreeIntervals" id="TreeIntervals" tree="@tree"/>\n');
        fprintf(g,'\t\t\t\t\t<dynamics spec="GLM" id="constant" maxRate="100" typeTrait="@typeTraitSet.t" types="%s">\n', sprintf('state%d ',[0:(states-1)]));
        fprintf(g,'\t\t\t\t\t\t<parameter name="rateShifts">%s</parameter>\n',sprintf('%.2f ', [20:20:200]));
        fprintf(g,'\t\t\t\t\t\t<migrationGLM spec="beast.mascot.glmmodel.LogLinear" id="migrationGLM">\n');
        
        splitn = strsplit(tree_files(i).name,'_');
        splitn = strsplit(splitn{2},'S');
        runnumber = str2double(splitn{2});
       
        
        for k = 1 : length(m_cov{runnumber}(1,:))
            fprintf(g,'\t\t\t\t\t\t\t<covariates spec="beast.core.parameter.RealParameter" id="m_co%d">%s</covariates>\n', k, sprintf('%f ',m_cov{runnumber}(:,k)));
        end
        fprintf(g,'\t\t\t\t\t\t\t<scaler idref="migrationScaler"/>\n');
        fprintf(g,'\t\t\t\t\t\t\t<indicator idref="migrationIndicator"/>\n');
        fprintf(g,'\t\t\t\t\t\t\t<clock idref="migrationClock"/>\n');
        fprintf(g,'\t\t\t\t\t\t</migrationGLM>\n');
        fprintf(g,'\t\t\t\t\t\t<NeGLM spec="beast.mascot.glmmodel.LogLinear" id="NeGLM">\n');
        for k = 1 : length(Ne_cov{runnumber}(1,:))
            fprintf(g,'\t\t\t\t\t\t\t<covariates spec="beast.core.parameter.RealParameter" id="Ne_co%d">%s</covariates>\n', k, sprintf('%f ',Ne_cov{runnumber}(:,k)));
        end
        fprintf(g,'\t\t\t\t\t\t\t<scaler idref="NeScaler"/>\n');
        fprintf(g,'\t\t\t\t\t\t\t<indicator idref="NeIndicator"/>\n');
        fprintf(g,'\t\t\t\t\t\t\t<clock idref="NeClock"/>\n');
        fprintf(g,'\t\t\t\t\t\t</NeGLM>\n');
        fprintf(g,'\t\t\t\t\t</dynamics>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t<distribution id="TreeLikelihood" spec="TreeLikelihood">\n');
        fprintf(g,'\t\t\t\t<data idref="sequences"/>\n');
        fprintf(g,'\t\t\t\t<tree idref="tree"/>\n');
        fprintf(g,'\t\t\t\t<siteModel spec="SiteModel" mutationRate="@clockrate">\n');
        fprintf(g,'\t\t\t\t\t<substModel spec="JukesCantor"/>\n');
        fprintf(g,'\t\t\t\t</siteModel>\n');
        fprintf(g,'\t\t\t</distribution>\n');        
        fprintf(g,'\t\t</distribution>\n');
        
        fprintf(g,'\t\t<operator id="NeClockScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@NeClock" scaleAll="true" scaleAllIndependently="true" weight="2.0"/>\n');
		fprintf(g,'\t\t<operator id="MigrationClockScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@migrationClock" scaleAll="true" scaleAllIndependently="true" weight="2.0"/>\n');

		fprintf(g,'\t\t<operator id="NeScalerScaler" spec="beast.mascot.operators.RealRandomWalkOperator" windowSizeOn="0.1" windowSizeOff="0.5" indicator="@NeIndicator" parameter="@NeScaler" useGaussian="true" weight="10.0"/>\n');
        fprintf(g,'\t\t<operator id="MigrationScalerScaler" spec="beast.mascot.operators.RealRandomWalkOperator" windowSizeOn="0.1" windowSizeOff="0.5" indicator="@migrationIndicator" parameter="@migrationScaler" useGaussian="true" weight="10.0"/>\n');
        
		fprintf(g,'\t\t<operator id="indicatorFlip.s:Migration" spec="BitFlipOperator" parameter="@migrationIndicator" weight="30"/>\n');
		fprintf(g,'\t\t<operator id="BSSVSoperator.c:Migration" spec="BitFlipBSSVSOperator" indicator="@migrationIndicator" mu="@migrationClock" weight="30"/>\n');
		fprintf(g,'\t\t<operator id="indicatorFlip.s:Ne" spec="BitFlipOperator" parameter="@NeIndicator" weight="30"/>\n');
		fprintf(g,'\t\t<operator id="BSSVSoperator.c:Ne" spec="BitFlipBSSVSOperator" indicator="@NeIndicator" mu="@NeClock" weight="30"/>\n');
       
%         fprintf(g,'\t\t<logger id="probtreelognud" fileName="%s.nud.trees" logEvery="5000" mode="tree">\n',flog);
%         fprintf(g,'\t\t\t<log id="logTreesnud" spec="StructuredTreeLogger" upDown="false" mascot="@coalescent"/>\n');
%         fprintf(g,'\t\t</logger>\n');
        fprintf(g,'\t\t<logger id="probtreelog" fileName="%s.trees" logEvery="5000" mode="tree">\n',flog);
        fprintf(g,'\t\t\t<log id="logTrees" spec="StructuredTreeLogger" mascot="@coalescent"/>\n');
        fprintf(g,'\t\t</logger>\n');

        fprintf(g,'\t\t<logger id="tracelog" fileName="%s.log" logEvery="200" model="@posterior" sanitiseHeaders="true" sort="smart">\n',flog);
        fprintf(g,'\t\t\t<log idref="posterior"/>\n');
        fprintf(g,'\t\t\t<log idref="prior"/>\n');
        fprintf(g,'\t\t\t<log idref="TreeLikelihood"/>\n');
        fprintf(g,'\t\t\t<log spec="RootStateLogger" id="RootStateLogger" mascot="@coalescent"/>\n');
        fprintf(g,'\t\t\t<log idref="nonZeroCovarates.Migration"/>\n');
        fprintf(g,'\t\t\t<log idref="nonZeroCovarates.Ne"/>\n');
        fprintf(g,'\t\t\t<log idref="migrationClock"/>\n');
        fprintf(g,'\t\t\t<log idref="NeClock"/>\n');
        fprintf(g,'\t\t\t<log idref="migrationGLM"/>\n');
        fprintf(g,'\t\t\t<log idref="NeGLM"/>\n');
%         fprintf(g,'\t\t\t<log idref="migrationScaler"/>\n');
%         fprintf(g,'\t\t\t<log idref="constant"/>\n');
%         fprintf(g,'\t\t\t<log idref="NeScaler"/>\n');
%         fprintf(g,'\t\t\t<log idref="NeIndicator"/>\n');
        fprintf(g,'\t\t</logger>\n');
        fprintf(g,'\t\t<logger id="screenlog" logEvery="1000">\n');
        fprintf(g,'\t\t\t<log idref="posterior"/>\n');
        fprintf(g,'\t\t</logger>\n');
        fprintf(g,'\t</run>\n');
        fprintf(g,'</beast>\n');
        fclose(f);
    end
end
