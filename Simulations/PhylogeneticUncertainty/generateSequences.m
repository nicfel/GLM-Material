%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script genereates Sequences for the simulated trees. The
% sequences are generated using an HKY model and seq-gen v1.3.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% add software folder to path
addpath('../../../Software')

% delete Sequences folder and create it newly
system('rm -r Sequences');
system('mkdir Sequences/')


for i = 201 : 500
    % read in the tree file
    disp(i)
    
    
    g = fopen(sprintf('master/Rates_S%d_master.tree', i),'r');    
    t = textscan(g,'%s'); fclose(g);
    
    % coalescing
    tree_tmp1 = regexprep(t{1}(end-1),'&type="L",location="(\d*)",reaction="Coalescence",time=(\d*).(\d*)','');
    
    % make the MASTER tree compatible with BEAST2
    % sampling
    tree_tmp1 = regexprep(tree_tmp1,'E[-](\d)]',']');
    tip_locs = regexp(tree_tmp1,'[&type="L",location="(\d*)",time=(\d*)\.(\d*)\]','match');
     
    for j = 1 : length(tip_locs{1})
        loc = regexprep(tip_locs{1}{j},'[&type="L",location="(\d*)",time=(\d*)\.(\d*)\]','$1');
        tree_tmp1 = strrep(tree_tmp1,tip_locs{1}{j},['loc_' loc 'kickout']);
        tree_tmp1 = strrep(tree_tmp1,'kickout','');
    end

    tree_tmp = regexprep(tree_tmp1,'(\d*)loc_','inv$1loc');
    
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
    
    g = fopen('tree.tree','w');
    fprintf(g, '#nexus\n');

    fprintf(g, 'Begin trees;\n');
    fprintf(g, 'tree TREE = ');
    fprintf(g, print_tree);
    fprintf(g, '\nEND;');

    fclose(g);
    
    disp(i)
    command = sprintf('%s -mHKY %s -l 1000 -s 0.001 < tree.tree > %s',...
        '../../Software/seq-gen1.3.3',...
            '-f0.25,0.25,0.25,0.25','seq.fasta');
    system(command);

    % read the sequence file
    s = fopen('seq.fasta', 'r');
    % read file
    seq = textscan(s, '%s');
    fclose(s);
    
    s = fopen(sprintf('Sequences/Rates_S%d_master.fasta',i), 'w');
    
    for j = 3 : length(seq{1})
        if length(seq{1}{j})>250
            line = regexprep(seq{1}{j},'inv(\d*)loc(\d*)','inv$1loc$2\n');            
        else
            line = seq{1}{j};
        end
        if ~isempty(strfind(line, 'inv'))
            fprintf(s, '>%s\n',line);
        else
            fprintf(s, '%s\n',line);
        end
    end
    fclose(s);
    
    system('rm tree.tree');system('rm seq.fasta');
end
% fclose(f_rates);
