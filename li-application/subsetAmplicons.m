clear
ts = {'T1', 'T2'};

nr_amplicons = 500

for tt = 1:length(ts)
    for sym = 1:2
        for rep = 0:2
            rng(1)
            f = fopen(sprintf('templates/%s_state_clocks_amplicon_only_weighted_fixed_topology.xml', ts{tt}));
            if sym==1
                g = fopen(sprintf('xmls/%s_amplicon_bidir_rep%d.xml', ts{tt},rep), 'w');
            else
                g = fopen(sprintf('xmls/%s_amplicon_unidir_rep%d.xml', ts{tt},rep), 'w');
            end

            seq = cell(0,0);
            rem_seq = cell(0,0);
            wgs = {'t1f4','f1l6','t1f2','t1l3','t1l10','t1f9','t1f14','t1l8','t1f11','t1l1', 't1z3','t1f23','t1f24','t1l13','t1z1','t1l6','t1z5',...
                't2z1','t2z6','t2z9', 't2z11','t2z13','t2f2','t2f5','t2f9','t2f13'};
            % wgs = {'t1f4','f1l6','t1f2','t1l3','t1l10','t1f9','t1f14','t1l8','t1f11','t1l1','t1l6', 't1f24','t1z1'};

            while ~feof(f)
                line = fgets(f);

                if contains(line, '<data id="HCCtumor"')
                    tmp = strsplit(line, '"');
                    weights = str2double(strsplit(strrep(tmp{8}, '1,','0,'), ','));
                    uni_weights = unique(weights);
                    
                    wgs_trait = 'rem';
                    
                    while ~contains(line, '</data>')
                        if contains(line, '<sequence')
                            tmp = strsplit(line, '"');
                            if rand<0.15 || ismember(tmp{6}, wgs)
                                seq{end+1} = line;
                                if ismember(tmp{6}, wgs)
                                    wgs_trait = [wgs_trait ',' tmp{6} '=1'];
                                else
                                    wgs_trait = [wgs_trait ',' tmp{6} '=0'];
                                end
                            else
                                rem_seq{end+1} = tmp{6};
                            end
                        end  
                        line = fgets(f);
                    end    
%                     weights(weights<10)=1;
                    uni_weights=uni_weights(1);
                    
                    
                    


                    for i = [1:length(uni_weights)]
                        variable_sites = find(weights<1000 & weights>0);
                        variable_weights = weights(weights<1000 & weights>0);

                        ratio = sum(weights(weights>1000))/sum(weights<1000);
                          
%                         ind_constant=sort(randsample(find(weights>1000), floor(nr_amplicons*ratio), true, weights(weights>1000)));
                        ind_amplicons=sort(randsample(variable_sites, nr_amplicons, true, variable_weights));
%                         ind = [ind_amplicons unique(ind_constant)];
                        ind = ind_amplicons;
                        newweights = ones(1,length(ind_amplicons));
%                         newweights(end+1) = sum(ind_constant==ind(end-3));
%                         newweights(end+1) = sum(ind_constant==ind(end-1));
%                         newweights(end+1) = sum(ind_constant==ind(end-2));
%                         newweights(end+1) = sum(ind_constant==ind(end));
                        
                        weightsstring = [sprintf('%d,', newweights) 'rem'];
                        
                        fprintf(g, '<data id="HCCtumor.%d" spec="Alignment" name="alignment" weights="%s">\n', i, strrep(weightsstring, ',rem',''));

                        
                        for j = 1:length(seq)
                            lprint = seq{j};
                            lprint = strrep(lprint, 'seq_',['seq' num2str(i) '_']);
                            tmp = strsplit(lprint, '"');
                            if uni_weights(i)<0
                                lprint = strrep(lprint, tmp{10}, tmp{10}(weights==uni_weights(i)) );
                            else
                                lprint = strrep(lprint, tmp{10}, tmp{10}(ind) );
                            end
                            fprintf(g, lprint);
                        end
                        fprintf(g, '</data>\n');
                    end

                  fprintf(g, '<data id="HCCtumor.forinit" spec="Alignment" name="alignment">\n');
                  for j = 1:length(seq)
                      lprint = seq{j};
                      lprint = strrep(lprint, 'seq_',['seqforinit_']);
                      tmp = strsplit(lprint, '"');
                      fprintf(g, lprint);
                  end
                  fprintf(g, '</data>\n');
                elseif contains(line, '<parameter id="migrationRateCanonical.t:HCCtumor"') && sym==2
                    
                    fprintf(g, strrep(line, '1 1','0 1') );

                elseif contains(line, '<typeTrait')
                    for i = 1 :length(rem_seq)
                        line = strrep(line, [rem_seq{i} '=loc0'], '');
                        line = strrep(line, [rem_seq{i} '=loc1'], '');
                        line = strrep(line, ',,', ',');
                        line = strrep(line, ',"', '"');
                        line = strrep(line, '",', '"');
                    end
                    fprintf(g, line);
                elseif contains(line, '<distribution id="treeLikelihood.HCCtumor.1"')
                    while ~contains(line, '</distribution>')
                        line = fgets(f);
                    end
                    inv_weights = 1./uni_weights;
                    for i = [1:length(uni_weights)]
                        disp(i)
%                         rate = inv_weights(i)/mean(inv_weights)/ sum(weights==uni_weights(i))*100;
%                         rate=1/ratio;
                        rate = 0.1;
            %             disp(inv_weights(i)/mean(inv_weights) / sum(weights==uni_weights(i)));
                        fprintf(g, '\t\t\t<distribution id="treeLikelihood.HCCtumor.%d" spec="ThreadedTreeLikelihood" data="@HCCtumor.%d" tree="@Tree.t:HCCtumor">\n', i, i);
                        fprintf(g, '\t\t\t\t<siteModel id="SiteModel.s:HCCtumor.%d" spec="SiteModel">\n', i);
            %             fprintf(g, '\t\t\t\t<siteModel id="SiteModel.s:HCCtumor.%d" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape.s:t1_wgs">\n', i);
                        fprintf(g, '\t\t\t\t\t<parameter id="mutationRate.s:t1_wgs.%d" spec="parameter.RealParameter" estimate="false" name="mutationRate">%f</parameter>\n', i, rate);
                        if i==1
                            fprintf(g, '\t\t\t\t\t<parameter id="proportionInvariant.s:t1_wgs" spec="parameter.RealParameter" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n');
                            fprintf(g, '\t\t\t\t\t<substModel id="gtr.s:t1_wgs" spec="GTR" rateAC="@rateAC.s:t1_wgs" rateAG="@rateAG.s:t1_wgs" rateAT="@rateAT.s:t1_wgs" rateCG="@rateCG.s:t1_wgs" rateGT="@rateGT.s:t1_wgs">\n');
                            fprintf(g, '\t\t\t\t\t\t<parameter id="rateCT.s:t1_wgs" spec="parameter.RealParameter" estimate="false" lower="0.0" name="rateCT">1.0</parameter>\n');
                            fprintf(g, '\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:t1_wgs" spec="Frequencies" frequencies="@freqParameter.s:t1_wgs"/>\n');
                            fprintf(g, '\t\t\t\t\t</substModel>\n');
                        else
                            fprintf(g, '\t\t\t\t\t<proportionInvariant idref="proportionInvariant.s:t1_wgs"/>\n');
                            fprintf(g, '\t\t\t\t\t<substModel idref="gtr.s:t1_wgs"/>\n');
                        end
                        fprintf(g, '\t\t\t\t</siteModel>\n');
                        if i==1
            %                 fprintf(g, '\t\t\t\t<branchRateModel id="StrictClock.c:HCCtumor" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c:HCCtumor"/>\n');

                            fprintf(g, '\t\t\t\t<branchRateModel id="StrictClock.c:HCCtumor" spec="bdmm.stateclock.BdmmStateClock" normalize="true">\n');
                            
                            fprintf(g, '\t\t\t\t\t\t<parameter idref="clockRate.c:HCCtumor" name = "clock.rate"/>\n');
                            fprintf(g, '\t\t\t\t\t\t<relativeClockRates idref="birthRateCanonical.t:HCCtumor"/>\n');
                            fprintf(g, '\t\t\t\t\t\t<relativeClockScalerRates idref="relativeScaler"/>\n');
                            fprintf(g, '\t\t\t\t\t\t<typedTree spec="bdmm.stateclock.FlatTypeMappedTree" parameterization="@CanonicalBDMMPrimeParameterization.t:HCCtumor" typeLabel="loc" frequencies="@typeFrequencies.t:HCCtumor" untypedTree="@Tree.t:HCCtumor" typeTraitSet="@typeTraitSet.t:HCCtumor">\n');
                            
                            
                            fprintf(g, '\t\t\t\t\t\t\t<isWGS id="iswgs.t:HCCtumor" spec="bdmmprime.util.InitializedTraitSet" traitname="wgs" value="%s">\n', strrep(wgs_trait, 'rem,',''));
                            fprintf(g, '\t\t\t\t\t\t\t\t<taxa id="TaxonSet.8" spec="TaxonSet" alignment="@HCCtumor.1"/>\n');
                            fprintf(g, '\t\t\t\t\t\t\t</isWGS>\n');
                            fprintf(g, '\t\t\t\t\t\t</typedTree>\n');
                            
                            fprintf(g, '\t\t\t\t</branchRateModel>\n');
                        else
                            fprintf(g, '\t\t\t\t<branchRateModel idref="StrictClock.c:HCCtumor"/>\n');
                        end
                        fprintf(g, '\t\t\t</distribution>\n');
                    end
                elseif contains(line, '@HCCtumor"') || contains(line, 'idref="HCCtumor"')
                    line = strrep(line, '@HCCtumor"','@HCCtumor.1"');
                    line = strrep(line, '"HCCtumor"','"HCCtumor.1"');
                    fprintf(g, line);
                elseif contains(line, '<log idref="treeLikelihood.HCCtumor"/>')
                else
                    fprintf(g, line);
                end
            end
            fclose('all')
        end
    end
end
        