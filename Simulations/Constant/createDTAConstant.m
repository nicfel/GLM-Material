%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates structured coalescent xmls from the master trees. Always creates
% 3 xmls per tree with different initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% create lisco files
tree_files = dir('master/*.tree');

system('rm -r dta');
system('mkdir dta');

states = 10;
% build the mapping matrix to get correct order of predictors
index=1;
for a = 1 : states
    for b = a + 1 : states
        mapper(a,b) = index;
        index = index+1;
    end
end

for a = 1 : states
    for b = a + 1 : states
        mapper(b,a) = index; 
        index = index+1;
    end
end



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
            %TODO make this the correct order for DTA
            for a = 1 : states
                for b = 1 : states
                    if a~=b
                        m_cov{str2num(tmp{1})}(mapper(a,b), i) = m(a,b);
                        count_b = count_b+1;
                    end
                end
            end
        end

    end
end
fclose(f);

% f = fopen('NeCovariates.txt');
% while ~feof(f)
%     line = fgets(f);
%     tmp = strsplit(line, '\t');
%     Ne_cov{str2num(tmp{1})} = str2num(tmp{2});
% end
% fclose(f);



%%

for i = 1 : length(tree_files)
    disp(tree_files(i).name)
    % read tree files
    g = fopen(['master/' tree_files(i).name],'r');
    t = textscan(g,'%s'); fclose(g);
    
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
    
    % get the tip heights
    distances = pdist(ptree, 'SquareForm',true, 'Nodes','all');
    for j = 1 : size(leafnames,1)
        heights(j) = distances(end,j);
    end
    
    heights = abs(heights-max(heights));
    
    
    print_tree = tree;
    
    % get the covariates
    
    % make tripletts of all runs with different random initial values
    for tr = 1 : 1    
        % make the xmls for dta
        flog = strrep(tree_files(i).name,'master.tree','dta');
        
        fname = sprintf('dta/%s.xml',flog);
        g = fopen(fname,'w');

        
        fprintf(g,'<?xml version="1.0" standalone="yes"?>\n');
        
        fprintf(g,'<beast>\n');
        fprintf(g,'\t<taxa id="taxa">\n');
        for j = 1 : length(leafnames)
%             seq_ind = find(ismember(taxa_name, leafs{j}));
            fprintf(g,'\t\t<taxon id="%s">\n',leafnames{j});
            fprintf(g,'\t\t\t<date value="%.8f" direction="backwards" units="years"/>\n', heights(j));
            fprintf(g,'\t\t\t<attr name="location">\n');
            tmp1 = strsplit(leafnames{j},'_');
            fprintf(g,'\t\t\t\t%s\n', tmp1{end});
            fprintf(g,'\t\t\t</attr>\n');
            fprintf(g,'\t\t</taxon>\n');
        end
        fprintf(g,'\t</taxa>\n');
        fprintf(g,'\t<alignment id="alignment" dataType="nucleotide">\n');

        for j = 1 : length(leafnames)
%             seq_ind = find(ismember(taxa_name, leafs{j}));
            fprintf(g,'\t\t<sequence>\n');
            fprintf(g,'\t\t\t<taxon idref="%s"/>\n', leafnames{j});
            fprintf(g,'\t\t\t\t???\n');
            fprintf(g,'\t\t</sequence>\n');
        end

        fprintf(g,'\t</alignment>\n');
        
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
        

        fprintf(f, '\n');
        fprintf(f, '\n');
        fprintf(f, '\n');
        fprintf(f, '\t<generalDataType id="location.dataType">\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t<state code="0"/>\n');
        fprintf(f, '\t\t<state code="1"/>\n');
        fprintf(f, '\t\t<state code="2"/>\n');
        fprintf(f, '\t\t<state code="3"/>\n');
        fprintf(f, '\t\t<state code="4"/>\n');
        fprintf(f, '\t\t<state code="5"/>\n');
        fprintf(f, '\t\t<state code="6"/>\n');
        fprintf(f, '\t\t<state code="7"/>\n');
        fprintf(f, '\t\t<state code="8"/>\n');
        fprintf(f, '\t\t<state code="9"/>\n');
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
        fprintf(f, '\t\t\t\t\t<parameter id="location.frequencies" dimension="%d"/>\n', states);
        fprintf(f, '\t\t\t\t</frequencies>\n');
        fprintf(f, '\t\t\t</frequencyModel>\n');
        fprintf(f, '\t\t</rootFrequencies>\n');
        fprintf(f, '\t\t<glmModel id="location.glmModel" checkFullRank="false" family="logLinear" checkIdentifiability="false">\n');
        fprintf(f, '\t\t\t<independentVariables>\n');
        fprintf(f, '\t\t\t\t<parameter id="location.glmCoefficients" dimension="%d" value="0.1"/>\n', size(m_cov{1},2));
        fprintf(f, '\t\t\t\t<indicator>\n');
        fprintf(f, '\t\t\t\t\t<parameter id="location.coefIndicator" dimension="%d" value="1"/>\n', size(m_cov{1},2)); 
        fprintf(f, '\t\t\t\t</indicator>\n');
        fprintf(f, '\t\t\t\t<designMatrix id="location.designMatrix">\n');        
        
        splitn = strsplit(tree_files(i).name,'_');
        splitn = strsplit(splitn{2},'S');
        runnumber = str2double(splitn{2});

        for k = 1 : size(m_cov{runnumber},2)
             fprintf(f, '\t\t\t\t\t<parameter id="mcov%d" value="%s"/>\n', k, sprintf('%f ',m_cov{runnumber}(:,k)));
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
        fprintf(f, '\t\t\t\t<parameter id="location.root.frequencies" dimension="%d"/>\n', states);
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
%         fprintf(f, '\t\t\t\t<parameter id="location.root.frequencies" dimension="%d"/>\n', states);
%         fprintf(f, '\t\t\t</frequencies>\n');
%         fprintf(f, '\t\t</frequencyModel>\n');
% 
%         fprintf(f, '\t\t<rewards>\n');
%         fprintf(f, '\t\t\t<parameter id="0_R" value="1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="1_R" value="0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="2_R" value="0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="3_R" value="0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="4_R" value="0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="5_R" value="0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="6_R" value="0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="7_R" value="0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="8_R" value="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0" />\n');
%         fprintf(f, '\t\t\t<parameter id="9_R" value="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0" />\n');
%         fprintf(f, '\t\t</rewards>\n');
%         fprintf(f, '\t</markovJumpsTreeLikelihood>\n');

        
        fprintf(f, '\t<sumStatistic id="location.nonZeroRates" elementwise="true">\n');
		fprintf(f, '\t\t<parameter idref="location.coefIndicator"/>\n');
        fprintf(f, '\t</sumStatistic>\n');
        
        
        
        fprintf(f, '\n');    
        fprintf(f, '\t<operators id="operators" optimizationSchedule="log">\n');
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
        fprintf(f, '\t<mcmc id="mcmc" chainLength="10000000" autoOptimize="true">\n');
        fprintf(f, '\t\t<posterior id="posterior">\n');
        fprintf(f, '\t\t\t<prior id="prior">\n');
        fprintf(f, '\t\t\t\t<exponentialPrior mean="1.0">\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="location.rate"/>\n');
        fprintf(f, '\t\t\t\t</exponentialPrior>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t\t<normalPrior mean="0" stdev="2">\n');
        fprintf(f, '\t\t\t\t\t<parameter idref="location.glmCoefficients"/>\n');
        fprintf(f, '\t\t\t\t</normalPrior>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t\t\t<poissonPrior mean="3.0" offset="0.0">\n');
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
        fprintf(f, '\t\t\t\t<markovJumpsTreeLikelihood idref="location.treeLikelihood"/>\n');
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
        fprintf(f, '\t\t<log id="fileLog1" logEvery="5000" fileName="%s.GLM.log" overwrite="false">\n',flog);
        fprintf(f, '\t\t\t<posterior idref="posterior"/>\n');
        fprintf(f, '\t\t\t<prior idref="prior"/>\n');
        fprintf(f, '\t\t\t<likelihood idref="likelihood"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="treeModel.rootHeight"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="cds.CP1.kappa"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="cds.CP2.kappa"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="cds.CP3.kappa"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="ig.kappa"/>\n');
%         fprintf(f, '\n');
%         fprintf(f, '\t\t\t<parameter idref="cds.CP1.alpha"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="cds.CP2.alpha"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="cds.CP3.alpha"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="ig.alpha"/>\n');
%         fprintf(f, '\t\t\t<compoundParameter idref="allMus"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="skygrid.precision"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="skygrid.logPopSize"/>\n');
%         fprintf(f, '\t\t\t<parameter idref="skygrid.cutOff"/>\n');
%         fprintf(f, '\t\t\t<treeLikelihood idref="cds.CP1.treeLikelihood"/>\n');
%         fprintf(f, '\t\t\t<treeLikelihood idref="cds.CP2.treeLikelihood"/>\n');
%         fprintf(f, '\t\t\t<treeLikelihood idref="cds.CP3.treeLikelihood"/>\n');
%         fprintf(f, '\t\t\t<treeLikelihood idref="ig.treeLikelihood"/>\n');
%         fprintf(f, '\t\t\t<gmrfSkyGridLikelihood idref="skygrid"/>\n');
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
        fprintf(f, '\t\t<log id="Makona_1610_cds_ig.GLM.LocationLog" logEvery="5000" fileName="%s.GLM.location.log">\n',flog);
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
%         fprintf(f, '\t\t<log id="completeJumpHistory" logEvery="500000" fileName="%s.GLM.history.log">\n',flog);
%         fprintf(f, '\t\t\t<completeHistoryLogger>\n');
%         fprintf(f, '\t\t\t\t<markovJumpsTreeLikelihood idref="location.treeLikelihood"/>\n');
%         fprintf(f, '\t\t\t</completeHistoryLogger>\n');
%         fprintf(f, '\t\t</log>\n');
        fprintf(f, '\n');
        fprintf(f, '\t\t<!-- write tree log to file                                                  -->\n');
        fprintf(f, '\t\t<logTree id="treeFileLog" logEvery="50000" nexusFormat="true" fileName="%s.trees" sortTranslationTable="true">\n',flog);
        fprintf(f, '\t\t\t<treeModel idref="treeModel"/>\n');
%         fprintf(f, '\t\t\t<trait name="rate" tag="rate">\n');
%         fprintf(f, '\t\t\t\t<strictClockBranchRate idref="branchRates"/>\n');
%         fprintf(f, '\t\t\t</trait>\n');
        fprintf(f, '\t\t\t<markovJumpsTreeLikelihood idref="location.treeLikelihood"/>\n');
        fprintf(f, '\t\t\t<posterior idref="posterior"/>\n');
        fprintf(f, '\t\t</logTree>\n');
        fprintf(f, '\t</mcmc>\n');
        fprintf(f, '\n');
        fprintf(f, '</beast>\n');
        
        fclose(f);          
        
        
        
    end
end
