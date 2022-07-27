clear


ts = {'T1', 'T2'};

newstates = {'t1f2=loc1,t1f4=loc0,t1f9=loc0,t1f11=loc0,t1f14=loc0,t1f23=loc0,t1f24=loc1,t1l1=loc1,t1l3=loc0,t1l6=loc0,t1l8=loc0,t1l10=loc0,t1l13=loc1,t1z1=loc1,t1z3=loc1,t1z5=loc1',...
               't2f2=loc0,t2f5=loc0,t2f9=loc1,t2f13=loc0,t2z1=loc1,t2z6=loc1,t2z9=loc1,t2z11=loc1,t2z13=loc1'};
rng(1); % ensures all subsets have the same sequences


clockmodel = {'state','strict'};

for tt = 1:length(ts)
    first = true;
    for clock = 1:2
        for sym = 1:4
            for rep = 0:2

                f = fopen(sprintf('templates/%s_wgs_state_clocks.xml', ts{tt}));

                if sym==1
                    g = fopen(sprintf('xmls/%s_wgs_oristates_bidir_%s_rep%d.xml', ts{tt}, clockmodel{clock}, rep), 'w');
                elseif sym==2
                    g = fopen(sprintf('xmls/%s_wgs_newstates_bidir_%s_rep%d.xml', ts{tt}, clockmodel{clock}, rep), 'w');
                elseif sym==3
                    g = fopen(sprintf('xmls/%s_wgs_oristates_unidir_%s_rep%d.xml', ts{tt}, clockmodel{clock}, rep), 'w');
                elseif sym==4
                    g = fopen(sprintf('xmls/%s_wgs_newstates_unidir_%s_rep%d.xml', ts{tt}, clockmodel{clock}, rep), 'w');
                end

                while ~feof(f)
                    line = fgets(f);
                    if contains(line, '<sequence')
                        tmp = strsplit(line, '"');
                        if first
                            ind = sort(randsample(length(tmp{10}), 10000,1));
                            first = false;
                        end
                        fprintf(g, strrep(line, tmp{10}, tmp{10}(ind)));
                    elseif contains(line, 'id="typeTraitSet.t:HCCtumor') && (sym==2 || sym==4)
                        tmp = strsplit(line, '"');
                        fprintf(g, strrep(line, tmp{8}, newstates{tt}));
                    elseif contains(line, '<parameter id="migrationRateCanonical.t:HCCtumor"') && sym>2                   
                        fprintf(g, strrep(line, '1 1','0 1'));
                    elseif contains(line, '<branchRateModel id="StrictClock.c:HCCtumor"') && clock==2
                        fprintf(g, strrep(line, 'bdmm.stateclock.BdmmStateClock" normalize="true"', 'beast.evolution.branchratemodel.StrictClockModel"'));
                        line = fgets(f);
                        fprintf(g, line);
                        line = fgets(f);
                        line = fgets(f);
                        line = fgets(f);
                    else
                        fprintf(g, line);
                    end
                end
                fclose('all');
            end
        end
    end
end