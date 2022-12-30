clear

f = fopen('nextstrain_analysis/t1_nt_muts.json');
t1muts = zeros(260000,1);
while ~feof(f)
    line = fgets(f);    
    if contains(line, '"muts": [')
        if contains(line, ']')
            continue
        end
        line = fgets(f);
        while ~contains(line, '],')
            tmp = strrep(strtrim(line),'"','');
            tmp = strrep(tmp, ',','');
            val = str2double(tmp(2:end-1));
            t1muts(val) = t1muts(val)+1;
            line = fgets(f);
        end
    else
    end
end

f = fopen('nextstrain_analysis/t2_nt_muts.json');
t2muts = zeros(260000,1);
while ~feof(f)
    line = fgets(f);    
    if contains(line, '"muts": [')
        if contains(line, ']')
            continue
        end
        line = fgets(f);
        while ~contains(line, '],')
            tmp = strrep(strtrim(line),'"','');
            tmp = strrep(tmp, ',','');
            val = str2double(tmp(2:end-1));
            t2muts(val) = t2muts(val)+1;
            line = fgets(f);
        end
    else
    end
end



ts = {'T1', 'T2', 'T1red'};

newstates = {'t1f2=loc1,t1f4=loc0,t1f9=loc0,t1f11=loc0,t1f14=loc0,t1f23=loc0,t1f24=loc1,t1l1=loc1,t1l3=loc0,t1l6=loc0,t1l8=loc0,t1l10=loc0,t1l13=loc1,t1z1=loc1,t1z3=loc1,t1z5=loc1',...
               't2f2=loc0,t2f5=loc0,t2f9=loc1,t2f13=loc0,t2z1=loc1,t2z6=loc1,t2z9=loc1,t2z11=loc1,t2z13=loc1'};
rng(1111); % ensures all subsets have the same sequences


clockmodel = {'state','strict'};

for dataset=1:3
    for tt = 1:length(ts)
        first = true;
        for clock = 1:2
            for sym = 1:4
                for rep = 0:2
                    if tt<3
                        f = fopen(sprintf('templates/%s_wgs_state_clocks.xml', ts{tt}));
                    elseif sym==3                         
                        f = fopen(sprintf('templates/%s_wgs_state_clocks.xml', ts{1}));
                    else
                        continue
                    end

                    if sym==1 && tt<3
                        g = fopen(sprintf('xmls/%s_wgs_oristates_bidir_%d_%s_rep%d.xml', ts{tt}, dataset, clockmodel{clock}, rep), 'w');
                    elseif sym==2 && tt<3
                        g = fopen(sprintf('xmls/%s_wgs_newstates_bidir_%d_%s_rep%d.xml', ts{tt}, dataset, clockmodel{clock}, rep), 'w');
                    elseif sym==3 
                        g = fopen(sprintf('xmls/%s_wgs_oristates_unidir_%d_%s_rep%d.xml', ts{tt}, dataset, clockmodel{clock}, rep), 'w');
                    elseif sym==4 && tt<3
                        g = fopen(sprintf('xmls/%s_wgs_newstates_unidir_%d_%s_rep%d.xml', ts{tt}, dataset, clockmodel{clock}, rep), 'w');
                    end

                    while ~feof(f)
                        line = fgets(f);
                        if contains(line, '<sequence')
                            tmp = strsplit(line, '"');
                            if first
                                if contains(ts{tt}, 'T1')
                                    use_muts = find(t1muts==1);
                                else
                                    use_muts = find(t2muts==1);
                                end
                                ind = sort(randsample(use_muts, 25000));
                                disp(length(tmp{10}))
                                first = false;
                            end
                            if tt < 3 || ~contains(line, 't1l13')
                                fprintf(g, strrep(line, tmp{10}, tmp{10}(ind)));
                            end
                        elseif contains(line, 'id="typeTraitSet.t:HCCtumor') && (sym==2 || sym==4)
                            tmp = strsplit(line, '"');
                            fprintf(g, strrep(line, tmp{8}, newstates{tt}));
                        elseif contains(line, 'id="typeTraitSet.t:HCCtumor') && tt==3
%                             tmp = strsplit(line, '"');
                            
                            fprintf(g, strrep(line, 't1l13=loc1,', ''));

                        elseif contains(line, '<parameter id="migrationRateCanonical.t:HCCtumor"') && sym>2                   
                            fprintf(g, strrep(line, '1 1','0 1'));
                        elseif contains(line, '<branchRateModel id="StrictClock.c:HCCtumor"') && clock==2
                            fprintf(g, strrep(line, 'sdevo.bdmm.BdmmStateClock" normalize="true"', 'beast.evolution.branchratemodel.StrictClockModel"'));
                            line = fgets(f);
                            fprintf(g, line);
                            line = fgets(f);
                            line = fgets(f);
                            line = fgets(f);
                        elseif contains(line, 'beast.util.Transform$LogTransform" f="@migrationRateCanonical.t:HCCtumor') && sym>2       
                            
                        else
                            fprintf(g, line);
                        end
                    end
                    fclose('all');
                end
            end
        end
    end
end