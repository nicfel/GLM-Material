% takes the xml without inferring the tree and adds operators to infer the
% tree as well
clear
xmls = dir('xmls/*.xml');
system('rm -r treeinference');
system('mkdir treeinference');

for i = 1 : length(xmls)
    f = fopen(['xmls/' xmls(i).name]);
    g = fopen(['treeinference/' strrep(xmls(i).name, '.xml', '_treeinf.xml')], 'w');
    while ~feof(f)
        line = fgets(f);
        if ~isempty(strfind(line, '	<operator id="BSSVSoperator.c:Ne"'))
            fprintf(g, line);
            % add all tre operators
            fprintf(g, '<operator id="TreeScaler.t" spec="ScaleOperator" tree="@tree" scaleFactor="0.9" optimise="false" weight="3.0"/>\n');
            fprintf(g, '<operator id="RootScaler.t" spec="ScaleOperator" tree="@tree" scaleFactor="0.9" rootOnly="true" optimise="false" weight="3.0"/>\n');
            fprintf(g, '<operator id="Uniform.t" spec="Uniform" tree="@tree" weight="10.0"/>\n');
            fprintf(g, '<operator id="SubtreeSlide.t" spec="SubtreeSlide" tree="@tree" size="0.05" weight="20.0"/>\n');
            fprintf(g, '<operator id="Narrow.t" spec="Exchange" tree="@tree" weight="20.0"/>\n');
            fprintf(g, '<operator id="Wide.t" spec="Exchange" tree="@tree" isNarrow="false" weight="3.0"/>\n');
            fprintf(g, '<operator id="WilsonBalding.t" spec="WilsonBalding" tree="@tree" weight="3.0"/>\n');                      
            
        elseif ~isempty(strfind(line, 'chainLength="1000000"'))
            fprintf(g, strrep(line, 'chainLength="1000000"','chainLength="5000000"'));
        elseif ~isempty(strfind(line, 'logEvery="200"'))
            fprintf(g, strrep(line, 'logEvery="200"','logEvery="10000"'));
        elseif ~isempty(strfind(line, 'logEvery="5000"'))
            fprintf(g, strrep(line, 'logEvery="5000"','logEvery="10000"'));
        else
            fprintf(g, line);
        end
    end
end