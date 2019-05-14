%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates master files from the template with the rates being linear
% combinations of covariates plus and error term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% open the template
f = fopen('Rates_master.xml','r');   

% cell to save the lines
temp_lines = cell(0,0);  

system('rm -r master');
system('mkdir master');


% while there are still lines
while ~feof(f) 
    % read line by line
    temp_lines{end+1,1} = fgets(f);   
end

% close the template xml    
fclose(f); 

states = 5;
nr_covariates = 10;
use_covariates = 2;

time_intervals = 0:20:180;

Ne_cov = fopen('NeCovariates.txt','w');
m_cov = fopen('mCovariates.txt','w');

Ne_params = fopen('NeParameters.txt','w');
m_params = fopen('mParameters.txt','w');

fprintf(Ne_params, 'Nr');
fprintf(m_params, 'Nr');
for i = 1 : nr_covariates
   fprintf(Ne_params, '\tNeGLMscaler.Ne_co%d',i);
   fprintf(m_params, '\tmigrationGLMscaler.m_co%d',i);
end
fprintf(Ne_params, '\n');
fprintf(m_params, '\n');



S = 201;


scalers = zeros(0,0);
while S <= 500
    filename = sprintf('Rates_S%d_master',S);
    
    % randomly draw nr_covariates covariates for migration and Ne's
    cov_migration = normrnd(0,1,states*(states-1)*length(time_intervals),nr_covariates);
    cov_Ne = normrnd(0,1,states*length(time_intervals),nr_covariates);
        
    ori_m_cov = cov_migration;
	ori_Ne_cov = cov_Ne;
    
    
    % standardise covariates
    for i = 1 : nr_covariates
        log_migration(:,i) = cov_migration(:,i) - mean(cov_migration(:,i));
        log_Ne(:,i) = cov_Ne(:,i) - mean(cov_Ne(:,i));
    end
    
    % standardise covariates
    for i = 1 : nr_covariates
        log_migration(:,i) = log_migration(:,i)./std(log_migration(:,i));
        log_Ne(:,i) = log_Ne(:,i)./std(log_Ne(:,i));
    end
   
    % randomly draw number of covariates to use
    active_migration = poissrnd(0.693);
    active_Ne =poissrnd(0.693);
        
    % randomly draw which indicators to use
    ind_migration = zeros(1,nr_covariates);
    ind_migration(randsample(nr_covariates,active_migration,false)) = 1;
    ind_Ne = zeros(1,nr_covariates);
    ind_Ne(randsample(nr_covariates,active_Ne,false)) = 1;
    
    % randomly draw the scalers
    scaler_migration = normrnd(0,1,1,nr_covariates);
    scaler_Ne = normrnd(0,1,1,nr_covariates);
    
    scaler_migration(~ind_migration) = 0;
    scaler_Ne(~ind_Ne) = 0;
    
    
    migration = exp(log_migration*scaler_migration');
    Ne = exp(log_Ne*scaler_Ne');
    
    % define params of the lognormal distribution of the Ne
    mean_ne = 10;
    sigma_ne = 0.5;
    mu_ne = log(mean_ne) - sigma_ne^2/2;

    % define params of the lognormal distribution of the migration rate
    mean_rea = 0.1;
    sigma_rea = 0.5;
    mu_rea = log(mean_rea) - sigma_rea^2/2;

    
    % sample the Ne and the migration scalers
    m_tot_scaler = lognrnd(mu_rea,sigma_rea);
    Ne_tot_scaler = lognrnd(mu_ne,sigma_ne);
    
    migration = migration * m_tot_scaler;
    Ne = Ne * Ne_tot_scaler;
    
    counter = 1;
    for a = 1 : states
        for b = 1 : states
            if a ~= b
                m_in = 'abcdef';
                for j = 1 : length(time_intervals)
                    migration((states*(states-1))*(j-1) + counter) =...
                        Ne(states*(j-1) + b)*migration((states*(states-1))*(j-1) + counter)/Ne(states*(j-1) + a);
                end
                counter = counter + 1;
           end
        end
    end

 
    ori_ne_scaler = scaler_Ne;
    ori_m_scaler = scaler_migration;
    
    % check if the maximal rates are below a 5
    if max(migration) < 1000 && max(1./(Ne)) < 1000
        fprintf(Ne_cov,'%d\t%s\n', S, mat2str(log_Ne));
        fprintf(m_cov,'%d\t%s\n', S, mat2str(log_migration));
        fprintf(Ne_params,'%d\t%s\n', S, sprintf('%.12f\t',ori_ne_scaler));
        fprintf(m_params,'%d\t%s\n', S, sprintf('%.12f\t',ori_m_scaler));
        S = S + 1;
        
        fname = sprintf('./master/%s.xml',filename);   % set the file name
        p = fopen(fname,'w');



        % print all the glm data to file

        ri = randi(states,20,1);
        for i = 1 : states
            sample_nr(i) =  sum(ri==i);
        end

        while min(sample_nr)==0
            ri = randi(states,20,1);
            for i = 1 : states
                sample_nr(i) =  sum(ri==i);
            end
        end


        sample_nr = sample_nr*20;

        Ne_nr = 1;
        m_nr = 1;

        for l = 1 : length(temp_lines)
            if ~isempty(strfind(temp_lines{l},'insert_coalescent'));
                for i = 1 : states
                    ne_in = 'abcdef';
                    for j = 1 : length(time_intervals)
                        ne_in = [ne_in ',' num2str(1/(2*Ne(states*(j-1) + i))) ':' num2str(time_intervals(j))];
                    end
                    ne_in = strrep(ne_in, 'abcdef,','');
                    fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%s">\n',ne_in);
                    fprintf(p, '\t\t\t\t\t2L[%d]:1 -> L[%d]:1\n',i-1,i-1);
                    fprintf(p, '\t\t\t\t</reaction>\n');
                end
            elseif  ~isempty(strfind(temp_lines{l},'insert_migration'));
                counter = 1;
                for a = 1 : states
                    for b = 1 : states
                        if a ~= b
                            m_in = 'abcdef';
                            for j = 1 : length(time_intervals)
                                m_in = [m_in ',' num2str(...
                                    migration((states*(states-1))*(j-1) + counter)...
                                    ) ':' num2str(time_intervals(j))];
                            end
                            m_in = strrep(m_in, 'abcdef,','');
                            fprintf(p, '\t\t\t\t<reaction spec=''Reaction'' rate="%s">\n',m_in);
                            fprintf(p, '\t\t\t\t\tL[%d]:1 -> L[%d]:1\n',a-1,b-1);
                            fprintf(p, '\t\t\t\t</reaction>\n');
                            counter = counter + 1;
                        end
                    end
                end
            elseif  ~isempty(strfind(temp_lines{l},'insert_samples')); 
                fprintf(p,'\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="1" time="0">\n');
                fprintf(p,'\t\t\t\t<population spec="Population" type="@L" location="0"/>\n');
                fprintf(p,'\t\t\t</lineageSeedMultiple>\n');

                for i = 1 : states
                    rest_samples = sample_nr(i);
                    next_samples = poissrnd(10);
                    while rest_samples > 0
                        time = 50*rand;
                        fprintf(p,'\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="%d" time="%.4f">\n',next_samples, time);
                        fprintf(p,'\t\t\t\t<population spec="Population" type="@L" location="%d"/>\n',i-1);
                        fprintf(p,'\t\t\t</lineageSeedMultiple>\n');
                        rest_samples = rest_samples - next_samples;
                        next_samples = poissrnd(10);
                    end
                end
            elseif ~isempty(strfind(temp_lines{l},'insert_dimension'));
                fprintf(p,'%s',strrep(temp_lines{l},'insert_dimension',num2str(states)));
            elseif ~isempty(strfind(temp_lines{l},'insert_filename'));
                fprintf(p,'%s',strrep(temp_lines{l},'insert_filename',filename));
            else
                fprintf(p,'%s',temp_lines{l});  % print line unchanged
            end
        end
        fclose(p); %close file again
    end
    
end
fclose('all')
